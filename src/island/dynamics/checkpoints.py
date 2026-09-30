"""Versioned JSON boundaries and explicit segment continuation. No live objects persisted."""

import os
import tempfile
from copy import deepcopy
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from island.chemistry.coordinate_stereo import (
    assigned_cip_labels,
    validate_coordinate_stereochemistry,
)
from island.core import Coordinates
from island.core.coordinate_provenance import (
    coordinate_hash,
    previous_coordinate_source,
    updated_coordinate_metadata,
)
from island.evaluation import OpenMMBoundPotential, PotentialEvaluator
from island.exceptions import (
    DynamicsCheckpointCompatibilityError as Incompatible,
)
from island.exceptions import (
    DynamicsCheckpointIOError,
    DynamicsInputError,
    DynamicsUnavailableError,
    EvaluationInputError,
    EvaluationUnavailableError,
    StereochemistryError,
)
from island.exceptions import (
    InvalidDynamicsCheckpointError as Invalid,
)
from island.minimization.integrity import number

from . import _checkpoint_data as data
from ._steps import baoab_drift, kick, verlet_drift
from .integrity import checked_evaluation, retained_count, same_model, verify_final
from .langevin_models import BAOAB, LangevinFrame, LangevinOptions
from .models import (
    DynamicsFrame,
    DynamicsOptions,
    kinetic_energy,
    owned_vectors,
    velocity_hash,
)
from .thermal import (
    GAS_CONSTANT,
    RNG_ALGORITHM,
    _normal_increments,
    _rng,
    mass_inventory,
)


@dataclass(frozen=True)
class DynamicsSegmentOptions:
    """Segment-local positive budgets only. Resume physical settings cannot be changed."""

    steps: int
    max_evaluations: int
    max_frames: int
    recording_interval: int = 1

    def __post_init__(self):
        if not all(data.integer(v, 1) for v in asdict(self).values()):
            raise DynamicsInputError("Segment budgets must be positive integers")


@dataclass(frozen=True)
class DynamicsCheckpoint:
    """Immutable canonical JSON plus SHA-256; payload access returns an owned copy."""

    payload_json: str
    content_checksum: str

    def __post_init__(self):
        self.validate_integrity()

    def validate_integrity(self):
        _validate_record(self, False)

    @property
    def payload(self):
        self.validate_integrity()
        return data.strict_load(self.payload_json)

    @property
    def state(self):
        p = self.payload
        return data.frame_from(p["state"], p["physical"]["integrator"] == BAOAB)

    @property
    def rng_state(self):
        p = self.payload
        return None if p["rng"] is None else p["rng"]["state"]

    @property
    def absolute_step(self):
        return self.payload["state"]["step"]

    @property
    def time_ps(self):
        return self.payload["state"]["time_ps"]

    def __deepcopy__(self, memo):
        self.validate_integrity()
        return self


@dataclass(frozen=True)
class DynamicsSegment:
    """Owned segment diagnostic, absolute frames, cumulative counters and RNG snapshot."""

    payload_json: str
    content_checksum: str

    def __post_init__(self):
        self.validate_integrity()

    def validate_integrity(self):
        _validate_record(self, True)

    @property
    def payload(self):
        self.validate_integrity()
        return data.strict_load(self.payload_json)

    @property
    def frames(self):
        p = self.payload
        return tuple(
            data.frame_from(f, p["physical"]["integrator"] == BAOAB)
            for f in p["frames"]
        )

    @property
    def final_state(self):
        p = self.payload
        return data.frame_from(p["state"], p["physical"]["integrator"] == BAOAB)

    @property
    def completed(self):
        return self.payload["diagnostic"]["completed"]

    @property
    def evaluations(self):
        return self.payload["lineage"][-1]["evaluations"]

    @property
    def completed_steps(self):
        row = self.payload["lineage"][-1]
        return row["end_step"] - row["start_step"]

    @property
    def termination_reason(self):
        return self.payload["lineage"][-1]["termination_reason"]

    def __deepcopy__(self, memo):
        self.validate_integrity()
        return self

    def to_system(self, system, *, allow_incomplete=False):
        p = self.payload
        if type(allow_incomplete) is not bool:
            raise DynamicsInputError("Explicit boolean allow_incomplete required")
        if not self.completed and not allow_incomplete:
            raise DynamicsInputError(
                "Incomplete segment requires allow_incomplete=True"
            )
        _validate_system(system, p)
        if not p["diagnostic"]["startup_verified"]:
            raise DynamicsInputError("Startup was not independently verified")
        frame = self.final_state
        try:
            expected = assigned_cip_labels(system)
            for f in self.frames:
                validate_coordinate_stereochemistry(
                    system, f.coordinates, expected, stage="apply dynamics segment"
                )
        except ImportError as error:
            raise DynamicsUnavailableError(
                "RDKit required for assigned stereochemistry"
            ) from error
        copied = deepcopy(system)
        copied.coordinates = Coordinates(frame.coordinates)
        source = (
            ("langevin" if p["physical"]["integrator"] == BAOAB else "nve")
            + "_dynamics"
            + ("" if self.completed else "_diagnostic")
        )
        provenance = {
            "coordinate_source": source,
            "original_coordinate_source": previous_coordinate_source(system.metadata),
            "initial_coordinate_fingerprint": p["frames"][0]["coordinate_fingerprint"],
            "coordinate_fingerprint": frame.coordinate_fingerprint,
            "velocity_fingerprint": frame.velocity_fingerprint,
            "trajectory_fingerprint": p["origin"]["trajectory_fingerprint"],
            "segment_checksum": self.content_checksum,
            "absolute_step": frame.step,
            "time_ps": frame.time_ps,
            "completed": self.completed,
            "molecular_system_is_complete_restart": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        copied.metadata = updated_coordinate_metadata(
            system,
            provenance,
            previous_fingerprint=coordinate_hash(
                {s: tuple(system.coordinates.get(s)) for s in frame.coordinates}
            ),
            dynamics=provenance,
        )
        copied.metadata.update(
            production_validated=False, simulation_readiness="not_established"
        )
        copied.validate()
        return copied


def _validate_record(record, segment):
    try:
        data.require(
            type(record.payload_json) is str and type(record.content_checksum) is str,
            "Invalid checkpoint envelope",
        )
        payload = data.strict_load(record.payload_json)
        data.require(
            data.canonical(payload) == record.payload_json,
            "Payload is not canonical strict JSON",
        )
        data.require(
            data.checksum(payload) == record.content_checksum,
            "Checkpoint content checksum mismatch",
        )
        data.validate_payload(payload, segment=segment)
    except Invalid:
        raise
    except Exception as error:
        raise Invalid(f"Malformed dynamics checkpoint/segment: {error}") from error


def _record(payload, cls):
    return cls(data.canonical(payload), data.checksum(payload))


def create_dynamics_checkpoint(segment):
    """Only verified completed/budget boundaries from the segment API are resumable."""
    if type(segment) is not DynamicsSegment:
        raise Invalid(
            "Checkpoint requires DynamicsSegment; old results lack captured RNG state"
        )
    p = segment.payload
    if (
        p["lineage"][-1]["termination_reason"]
        not in ("completed", "maximum_evaluations", "maximum_frames")
        or not p["diagnostic"]["final_verified"]
    ):
        raise Invalid(
            "Checkpoint requires independently verified completion or pre-trial budget stop"
        )
    p.pop("frames")
    p.pop("diagnostic")
    p["schema"] = data.SCHEMA
    return _record(p, DynamicsCheckpoint)


def save_dynamics_checkpoint(checkpoint, path, *, replace=False):
    """Strict JSON, exclusive create by default; explicit atomic replacement."""
    if type(checkpoint) is not DynamicsCheckpoint:
        raise Invalid("Expected DynamicsCheckpoint")
    checkpoint.validate_integrity()
    if type(replace) is not bool:
        raise DynamicsInputError("replace must be boolean")
    temporary = None
    try:
        path = Path(path)
        envelope = (
            data.canonical(
                {"payload": checkpoint.payload, "sha256": checkpoint.content_checksum}
            )
            + "\n"
        )
        fd, temporary = tempfile.mkstemp(
            prefix="." + path.name + ".", suffix=".tmp", dir=path.parent
        )
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(envelope)
            stream.flush()
            os.fsync(stream.fileno())
        if replace:
            os.replace(temporary, path)
        else:
            os.link(
                temporary, path
            )  # exclusive publication, never truncates existing file
    except (OSError, TypeError, ValueError) as error:
        raise DynamicsCheckpointIOError(f"Cannot save checkpoint: {error}") from error
    finally:
        if temporary is not None:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass


def load_dynamics_checkpoint(path):
    try:
        text = Path(path).read_text(encoding="utf-8")
    except (OSError, TypeError, ValueError) as error:
        raise DynamicsCheckpointIOError(f"Cannot read checkpoint: {error}") from error
    try:
        envelope = data.strict_load(text)
        data.keys(envelope, {"payload", "sha256"})
        return DynamicsCheckpoint(
            data.canonical(envelope["payload"]), envelope["sha256"]
        )
    except Invalid:
        raise
    except Exception as error:
        raise Invalid(f"Malformed checkpoint JSON: {error}") from error


def _validate_system(system, p):
    try:
        masses = mass_inventory(system)
        if (
            data.rows(masses) != p["masses"]
            or data.compatibility(system) != p["system_fingerprint"]
        ):
            raise Incompatible(
                "Checkpoint chemical system, masses, or non-coordinate provenance differ"
            )
    except Incompatible:
        raise
    except Exception as error:
        raise Incompatible(str(error)) from error


def run_dynamics_segment(
    system, evaluator, velocities, options, *, velocity_unit="angstrom/ps"
):
    """Start a checkpoint-capable trajectory with existing NVE or Langevin options."""
    if type(options) not in (DynamicsOptions, LangevinOptions):
        raise DynamicsInputError("Supply DynamicsOptions or LangevinOptions")
    if velocity_unit != "angstrom/ps":
        raise DynamicsInputError("Initial velocities must be explicitly in angstrom/ps")
    options.__post_init__()
    budgets = DynamicsSegmentOptions(
        options.steps,
        options.max_evaluations,
        options.max_frames,
        options.recording_interval,
    )
    return _run(system, evaluator, velocities, data.physical(options), budgets, None)


def resume_dynamics(checkpoint, system, evaluator, options):
    """Continue using external compatible inputs and new execution budgets only."""
    if type(checkpoint) is not DynamicsCheckpoint:
        raise Invalid("Expected DynamicsCheckpoint")
    p = checkpoint.payload
    if p["environment"] != data.environment():
        raise Incompatible(
            "Checkpoint execution environment differs (strict v1 policy)"
        )
    if type(options) is not DynamicsSegmentOptions:
        raise DynamicsInputError("Resume accepts execution-only DynamicsSegmentOptions")
    _validate_system(system, p)
    return _run(
        system,
        evaluator,
        data.inventory(p["state"]["velocities"]),
        p["physical"],
        options,
        checkpoint,
    )


def _run(system, evaluator, velocities, physical, options, checkpoint):
    options.__post_init__()
    prior = None if checkpoint is None else checkpoint.payload
    thermal = physical["integrator"] == BAOAB
    dt = physical["timestep_fs"] * 0.001
    try:
        if not isinstance(evaluator, PotentialEvaluator):
            raise DynamicsInputError("Expected PotentialEvaluator")
        owned = deepcopy(system)
        masses = mass_inventory(owned)
        ids = tuple(masses)
        if isinstance(evaluator, OpenMMBoundPotential):
            evaluator.validate_system(owned)
        identity = data.compatibility(owned)
        start = 0 if prior is None else prior["state"]["step"]
        if (
            not number((start + options.steps) * dt)
            or (start + options.steps) * dt <= start * dt
        ):
            raise DynamicsInputError("Requested absolute duration is not representable")
        coords = owned_vectors(
            {s: owned.coordinates.get(s) for s in ids}
            if prior is None
            else data.inventory(prior["state"]["coordinates"])
        )
        speeds = owned_vectors(velocities)
        if set(speeds) != set(ids):
            raise DynamicsInputError("Exact velocity coverage required")
        expected = assigned_cip_labels(owned)
        validate_coordinate_stereochemistry(
            owned, coords, expected, stage="segment start"
        )
        kinetic = kinetic_energy(masses, speeds)
        if not number(kinetic) or not number(
            2 * kinetic / (3 * len(ids) * GAS_CONSTANT)
        ):
            raise DynamicsInputError("Unsafe initial kinetic energy")
        mass_array = np.array(list(masses.values()))[:, None]
        rng = None
        if thermal:
            if prior is None:
                rng = _rng(physical["thermostat_seed"])
            else:
                # The actual captured state is installed; counters are never used to advance it.
                rng = np.random.Generator(np.random.PCG64(0))
                rng.bit_generator.state = deepcopy(prior["rng"]["state"])
            product = physical["friction_per_ps"] * dt
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                c = np.exp(-product)
                sigma = np.sqrt(
                    -np.expm1(-2 * product)
                    * (100 * GAS_CONSTANT * physical["temperature_kelvin"] / mass_array)
                )
            if not np.isfinite(sigma).all():
                raise DynamicsInputError("Unsafe thermostat scale")
    except ImportError as error:
        raise DynamicsUnavailableError(
            "RDKit required for assigned stereochemistry"
        ) from error
    except DynamicsInputError:
        raise
    except Exception as error:
        raise DynamicsInputError(str(error)) from error
    calls = 0
    failed_calls = 0
    random_steps = 0
    startup = False
    final_attempted = False
    final = None
    final_error = None
    reason = "completed"
    message = "Requested segment completed"
    stage = None
    attempt = None
    initial_energy = None if prior is None else prior["origin"]["initial_total_energy"]
    maximum = None if prior is None else prior["max_abs_energy_deviation"]
    saved = (
        None if prior is None else data.evaluation_from(prior["state"]["evaluation"])
    )
    model = saved

    def evaluate(x, *, fresh=False):
        nonlocal calls
        calls += 1
        method = (
            evaluator.evaluate_fresh
            if fresh and isinstance(evaluator, OpenMMBoundPotential)
            else evaluator.evaluate
        )
        record = method(dict(x), coordinate_unit="angstrom")
        checked_evaluation(record, x)
        if model is not None:
            try:
                same_model(model, record)
            except Exception as error:
                if prior is not None and calls == 1:
                    incompatible = Incompatible(
                        "Startup model/parameter/backend/settings identity differs"
                    )
                    incompatible.evaluations = calls
                    raise incompatible from error
                raise
        return record

    def frame(step, x, v, record):
        k = kinetic_energy(masses, v)
        total = None if record is None else record.potential_energy + k
        deviation = (
            None
            if total is None
            else total - (total if initial_energy is None else initial_energy)
        )
        if (
            not number(k)
            or not number(2 * k / (3 * len(ids) * GAS_CONSTANT))
            or (total is not None and (not number(total) or not number(deviation)))
        ):
            raise EvaluationInputError("Nonfinite full-step energy")
        args = (
            step,
            step * dt,
            x,
            v,
            record,
            k,
            total,
            deviation,
            coordinate_hash(x),
            velocity_hash(v),
        )
        return (
            LangevinFrame(*args, 3 * len(ids), 2 * k / (3 * len(ids) * GAS_CONSTANT))
            if thermal
            else DynamicsFrame(*args)
        )

    current = frame(start, coords, speeds, saved)
    try:
        checked = evaluate(coords, fresh=True)
        if saved is not None:
            verify_final(saved, checked)
        current = frame(start, coords, speeds, checked)
        model = checked
        startup = True
        if initial_energy is None:
            initial_energy = current.total_energy
        maximum = max(maximum or 0.0, abs(current.energy_deviation))
        if not thermal and maximum > physical["max_energy_deviation"]:
            raise EvaluationInputError("Startup violates original NVE guard")
    except EvaluationUnavailableError as error:
        raise DynamicsUnavailableError(str(error)) from error
    except Incompatible:
        raise
    except Exception as error:  # noqa: BLE001 -- retained runtime diagnostic
        # Preserve the saved boundary on failed resume verification; never propagate it.
        current = frame(start, coords, speeds, saved)
        startup = False
        reason = "startup_verification_failed"
        stage = "startup"
        attempt = start
        message = str(error) or type(error).__name__
    stored = [current]
    if startup:
        while current.step - start < options.steps:
            next_step = current.step + 1
            if calls >= options.max_evaluations - 1:
                reason = "maximum_evaluations"
                stage = "budget"
                attempt = next_step
                message = "Final evaluator check reserved"
                break
            if (
                retained_count(next_step - start, options.recording_interval)
                > options.max_frames
            ):
                reason = "maximum_frames"
                stage = "budget"
                attempt = next_step
                message = "Frame capacity reached before trial"
                break
            before = calls
            if thermal:
                random_steps += 1
            try:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    noise = _normal_increments(rng, (len(ids), 3)) if thermal else None
                    x = np.array([current.coordinates[s] for s in ids])
                    v = np.array([current.velocities[s] for s in ids])
                    force = np.array([current.evaluation.forces[s] for s in ids])
                    if thermal:
                        x, half = baoab_drift(
                            x, v, force, mass_array, dt, c, sigma, noise
                        )
                    else:
                        x, half = verlet_drift(x, v, force, mass_array, dt)
                    trial_x = owned_vectors(dict(zip(ids, x, strict=True)))
                    validate_coordinate_stereochemistry(
                        owned, trial_x, expected, stage=f"segment trial {next_step}"
                    )
                    record = evaluate(trial_x)
                    next_force = np.array([record.forces[s] for s in ids])
                    trial_v = owned_vectors(
                        dict(
                            zip(
                                ids,
                                kick(half, next_force, mass_array, dt),
                                strict=True,
                            )
                        )
                    )
                    trial = frame(next_step, trial_x, trial_v, record)
                if (
                    not thermal
                    and abs(trial.energy_deviation) > physical["max_energy_deviation"]
                ):
                    reason = "energy_guard_exceeded"
                    raise EvaluationInputError(
                        "Original trajectory NVE energy guard exceeded"
                    )
                current = trial
                maximum = max(maximum, abs(current.energy_deviation))
                if (current.step - start) % options.recording_interval == 0:
                    stored.append(current)
            except (ImportError, EvaluationUnavailableError) as error:
                raise DynamicsUnavailableError(str(error)) from error
            except Exception as error:  # noqa: BLE001 -- retained runtime diagnostic
                if reason != "energy_guard_exceeded":
                    reason = (
                        "stereochemistry_changed"
                        if isinstance(error, StereochemistryError)
                        else "invalid_geometry"
                        if isinstance(
                            error,
                            (EvaluationInputError, FloatingPointError, OverflowError),
                        )
                        else "evaluation_failed"
                    )
                stage = "trial"
                attempt = next_step
                message = str(error) or type(error).__name__
                failed_calls = calls - before
                break
    if stored[-1].step != current.step:
        stored.append(current)
    if startup:
        if calls < options.max_evaluations:
            final_attempted = True
            try:
                checked = evaluate(current.coordinates, fresh=True)
                verify_final(current.evaluation, checked)
                final = checked
            except EvaluationUnavailableError as error:
                raise DynamicsUnavailableError(str(error)) from error
            except Exception as error:  # noqa: BLE001 -- retained runtime diagnostic
                final_error = str(error) or type(error).__name__
        else:
            final_error = "No capacity for independent final check"
        if final is None and reason == "completed":
            reason = "final_evaluation_failed"
            stage = "final"
            attempt = current.step
            message = final_error
    origin = (
        prior["origin"]
        if prior is not None
        else {
            "trajectory_fingerprint": "pending",
            "initial_total_energy": initial_energy,
            "initial_state": data.frame_data(stored[0]),
        }
    )
    lineage = [] if prior is None else deepcopy(prior["lineage"])
    lineage.append(
        {
            "index": len(lineage),
            "start_step": start,
            "end_step": current.step,
            "evaluations": calls,
            "parent_checksum": None
            if checkpoint is None
            else checkpoint.content_checksum,
            "options": asdict(options),
            "termination_reason": reason,
        }
    )
    counters = {
        "evaluations": calls
        + (0 if prior is None else prior["counters"]["evaluations"]),
        "random_steps": random_steps
        + (0 if prior is None else prior["counters"]["random_steps"]),
        "normal_draws": random_steps * 3 * len(ids)
        + (0 if prior is None else prior["counters"]["normal_draws"]),
    }
    payload = {
        "schema": data.SEGMENT_SCHEMA,
        "units": data.UNITS,
        "physical": physical,
        "environment": data.environment(),
        "system_fingerprint": identity,
        "origin": origin,
        "lineage": lineage,
        "state": data.frame_data(current),
        "frames": [data.frame_data(f) for f in stored],
        "masses": data.rows(masses),
        "rng": None
        if not thermal
        else {
            "algorithm": RNG_ALGORITHM,
            "numpy_version": np.__version__,
            "state": deepcopy(rng.bit_generator.state),
        },
        "counters": counters,
        "max_abs_energy_deviation": maximum,
        "final_evaluation": data.evaluation_data(final),
        "diagnostic": {
            "completed": reason == "completed",
            "startup_verified": startup,
            "final_verified": final is not None,
            "final_attempted": final_attempted,
            "failed_trial_evaluations": failed_calls,
            "attempted_step": attempt,
            "failure_stage": stage,
            "message": message,
            "final_error": final_error,
        },
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    if prior is None:
        payload["origin"]["trajectory_fingerprint"] = data.trajectory_identity(payload)
    return _record(payload, DynamicsSegment)
