"""Versioned prepared-bundle workflows; no construction or scientific preparation."""

from copy import deepcopy
from dataclasses import asdict, dataclass, field, fields
from functools import wraps
from pathlib import Path
from uuid import uuid4

from island.dynamics import (
    DynamicsSegmentOptions,
    LangevinOptions,
    create_dynamics_checkpoint,
    initialize_velocities,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    save_dynamics_checkpoint,
)
from island.dynamics._checkpoint_data import checksum as payload_checksum
from island.exceptions import WorkflowError
from island.forcefields import create_evaluator, load_prepared_forcefield
from island.minimization import MinimizationOptions, minimize_geometry

from . import bundle, storage
from .chain import _load_segment, _manifest, _put, _segment_count
from .consistency import (
    compatible_boundary,
    validate_prepared_attempt,
    validate_prepared_trajectory,
)

SCHEMA = "island_prepared_workflow_v1"
SETUP_SCHEMA = "island_prepared_workflow_setup_v1"


def _boundary(fn):
    @wraps(fn)
    def call(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except WorkflowError:
            raise
        except Exception as error:
            raise WorkflowError(f"{fn.__name__}: {error}") from error

    return call


@dataclass(frozen=True)
class PreparedWorkflowConfig:
    """Finite execution budgets only; source/charge/pair policies come from I2."""

    minimization: MinimizationOptions = field(default_factory=MinimizationOptions)
    temperature_kelvin: float = 300.0
    friction_per_ps: float = 5.0
    timestep_fs: float = 0.1
    velocity_seed: int = 78123
    thermostat_seed: int = 99181
    total_steps: int = 200
    segment_steps: int = 100
    recording_interval: int = 20
    max_evaluations_per_segment: int = 102
    max_frames_per_segment: int = 6
    max_segments: int = 100
    max_artifact_bytes: int = 100_000_000
    schema: str = "island_prepared_workflow_config_v1"

    @_boundary
    def __post_init__(self):
        if self.schema != "island_prepared_workflow_config_v1":
            raise WorkflowError("Unsupported prepared workflow configuration")
        for name in (
            "total_steps",
            "segment_steps",
            "max_segments",
            "max_artifact_bytes",
        ):
            if type(getattr(self, name)) is not int or getattr(self, name) <= 0:
                raise WorkflowError(f"{name} must be a positive integer")
        if type(self.velocity_seed) is not int or not 0 <= self.velocity_seed < 2**128:
            raise WorkflowError("Invalid velocity seed")
        if type(self.minimization) is not MinimizationOptions:
            raise WorkflowError("MinimizationOptions required")
        self.minimization.__post_init__()
        self.langevin(self.segment_steps)

    def langevin(self, steps):
        return LangevinOptions(
            self.timestep_fs,
            self.temperature_kelvin,
            self.friction_per_ps,
            self.thermostat_seed,
            steps,
            self.max_evaluations_per_segment,
            self.max_frames_per_segment,
            self.recording_interval,
        )

    def to_dict(self):
        return asdict(self)

    @classmethod
    @_boundary
    def from_dict(cls, value):
        if type(value) is not dict or set(value) != {f.name for f in fields(cls)}:
            raise WorkflowError("Malformed prepared workflow configuration")
        value = deepcopy(value)
        value["minimization"] = MinimizationOptions(**value["minimization"])
        return cls(**value)


def _identities(system, prepared):
    # The exact existing native fingerprint definitions, without Context creation.
    from island.evaluation.models import fingerprint

    native = prepared._native
    family = prepared._family
    if family in ("gaff", "gaff2"):
        from island.evaluation.openmm import SETTINGS, _graph

        imported = native.imported_result
        return imported.content_signature(), fingerprint(
            {
                "schema": "island_openmm_singlepoint_v1",
                "graph": _graph(system),
                "parameters": imported.model_content_signature(),
                "settings": SETTINGS,
            }
        )
    if family == "oplsaa":
        from island.evaluation.oplsaa import SETTINGS

        return native.identity, fingerprint(
            {"parameters": native.identity, "settings": SETTINGS}
        )
    from island.evaluation.pcff_identity import pcff_evaluation_identity

    _, parameter, model = pcff_evaluation_identity(native, system=system)
    return parameter, model


def _load_inputs(root, m, sources):
    if (
        storage.decode(storage.read_json(storage.child(root, m["configuration"])))
        != m["config"]
    ):
        raise WorkflowError("Configuration contradicts immutable configuration record")
    loaded = load_prepared_forcefield(Path(root) / "prepared", sources=sources)
    expected = {"prepared/manifest.json"} | {
        "prepared/" + row["path"]
        for row in loaded.manifest["payload"]["files"].values()
    }
    if set(m["input_files"]) != expected or not expected <= m["files"].keys():
        raise WorkflowError("Prepared input artifact inventory mismatch")
    if loaded.prepared.identity != m["prepared_identity"]:
        raise WorkflowError("Prepared facade identity mismatch")
    return loaded


def _minimum_record(payload):
    """Reconstruct even unsuccessful diagnostics without granting application eligibility."""
    from island.evaluation.models import EvaluationResult
    from island.minimization.models import MinimizationResult, MinimizationStep

    p = deepcopy(payload)
    for name in ("initial_evaluation", "final_evaluation"):
        if p[name] is not None:
            p[name] = EvaluationResult(**p[name])
    p["options"] = MinimizationOptions(**p["options"])
    p["history"] = tuple(MinimizationStep(**r) for r in p["history"])
    result = MinimizationResult(**p)
    result.validate_integrity()
    return result


def _load_setup(root, m, sources, *, loaded=None, with_minimum=False):
    if m["bundle"] is None:
        raise WorkflowError(
            "No verified setup published; failed setup requires a new run"
        )
    p = storage.decode(storage.read_json(storage.child(root, m["bundle"])))
    if (
        set(p)
        != {
            "schema",
            "config",
            "prepared_identity",
            "minimum",
            "starting_system",
            "initialization",
        }
        or p["schema"] != SETUP_SCHEMA
    ):
        raise WorkflowError("Malformed workflow setup")
    if p["config"] != m["config"] or p["prepared_identity"] != m["prepared_identity"]:
        raise WorkflowError("Setup/config/prepared identity mismatch")
    loaded = loaded or _load_inputs(root, m, sources)
    original, prepared = loaded.system, loaded.prepared
    minimum = bundle.minimum_from(p["minimum"])
    if dict(minimum.initial_coordinates) != {
        s: tuple(original.coordinates.get(s)) for s in original.topology.sites
    }:
        raise WorkflowError("Minimum original coordinates differ from preparation")
    starting = minimum.to_system(original)
    if bundle.system_data(starting) != p["starting_system"]:
        raise WorkflowError("Starting system differs from applied minimum")
    prepared.validate_integrity(starting)
    initialization = bundle.initialization_from(p["initialization"])
    from island.dynamics.thermal import mass_inventory
    from island.minimization.models import system_identity

    config = PreparedWorkflowConfig.from_dict(m["config"])
    if (
        initialization.system_fingerprint != system_identity(starting)
        or dict(initialization.masses) != mass_inventory(starting)
        or initialization.temperature_kelvin != config.temperature_kelvin
        or initialization.velocity_seed != config.velocity_seed
        or minimum.options != config.minimization
    ):
        raise WorkflowError("Minimum/initialization/configuration mismatch")
    parameter, model = _identities(starting, prepared)
    if (
        minimum.final_evaluation.parameter_fingerprint != parameter
        or minimum.final_evaluation.model_fingerprint != model
    ):
        raise WorkflowError("Minimum native parameter/model mismatch")
    result = starting, prepared, initialization
    return (*result, minimum) if with_minimum else result


def _publish_manifest(root, m, sources):
    _validate_manifest(root, m, sources)
    _manifest(root, m)


def _failure(root, m, stage, error, sources):
    candidate = deepcopy(m)
    candidate["status"] = "stage_failed"
    candidate["stages"].append(
        {
            "stage": stage,
            "status": "failed",
            "error_type": type(error).__name__,
            "message": str(error),
        }
    )
    _publish_manifest(root, candidate, sources)
    return candidate


@_boundary
def start_prepared_bundle_workflow(
    prepared_directory, directory, config, *, sources=None, segments=1
):
    """Copy checked I2 artifacts, minimize once, initialize once, then propagate."""
    if type(config) is not PreparedWorkflowConfig:
        raise WorkflowError("PreparedWorkflowConfig required")
    config.__post_init__()
    _segment_count(segments)
    loaded = load_prepared_forcefield(prepared_directory, sources=sources)
    original, prepared = loaded.system, loaded.prepared
    # Snapshot only manifest-listed artifacts. Reload the copied bytes before setup.
    input_root = Path(prepared_directory)
    names = [
        "manifest.json",
        *(r["path"] for r in loaded.manifest["payload"]["files"].values()),
    ]
    raw = {n: storage.child(input_root, n).read_bytes() for n in names}
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=False)
    with storage.writer(root):
        (root / "prepared").mkdir()
        m = {
            "schema": SCHEMA,
            "run_id": uuid4().hex,
            "config": config.to_dict(),
            "prepared_identity": prepared.identity,
            "input_files": [],
            "configuration": None,
            "status": "starting",
            "stages": [],
            "files": {},
            "bundle": None,
            "segments": [],
            "checkpoint": None,
            "accepted_step": 0,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        for name, value in raw.items():
            name = "prepared/" + name
            storage.publish(storage.child(root, name), value)
            m["files"][name] = {"bytes": len(value), "sha256": storage.checksum(value)}
            m["input_files"].append(name)
        m["configuration"] = _put(
            root,
            m,
            "configuration",
            storage.json_bytes(storage.encode(config.to_dict())),
        )
        _publish_manifest(root, m, sources)
        stage = "minimization"
        try:
            evaluator = create_evaluator(original, prepared)
            with evaluator.open_session() as session:
                minimum = minimize_geometry(original, session, config.minimization)
            minimum.validate_integrity()
            record = _put(
                root,
                m,
                "minimum",
                storage.json_bytes(storage.encode(bundle.record(minimum))),
            )
            passed = minimum.converged and minimum.final_evaluation_verified
            m["stages"].append(
                {
                    "stage": stage,
                    "status": "passed" if passed else "failed",
                    "record": record,
                    "evaluations": minimum.evaluations,
                    "reason": minimum.termination_reason,
                }
            )
            if not passed:
                m["status"] = "stage_failed"
                _publish_manifest(root, m, sources)
                return deepcopy(m)
            starting = minimum.to_system(original)
            stage = "initialization"
            initialization = initialize_velocities(
                starting,
                temperature_kelvin=config.temperature_kelvin,
                seed=config.velocity_seed,
            )
            p = {
                "schema": SETUP_SCHEMA,
                "config": config.to_dict(),
                "prepared_identity": prepared.identity,
                "minimum": bundle.record(minimum),
                "starting_system": bundle.system_data(starting),
                "initialization": bundle.record(initialization),
            }
            m["bundle"] = _put(root, m, "setup", storage.json_bytes(storage.encode(p)))
            m["status"] = "ready"
            m["stages"].append({"stage": stage, "status": "passed"})
            _publish_manifest(root, m, sources)
        except Exception as error:  # noqa: BLE001 -- preserve diagnostic stage outcome
            return _failure(root, m, stage, error, sources)
        return _advance(root, m, segments, sources)


@_boundary
def prepared_workflow_status(directory, *, sources=None):
    """Deep offline validation. No Context, preparation or energy evaluation."""
    root = Path(directory)
    envelope = storage.read_json(root / "manifest.json")
    if (
        set(envelope) != {"payload", "sha256"}
        or payload_checksum(envelope["payload"]) != envelope["sha256"]
    ):
        raise WorkflowError("Manifest envelope/checksum mismatch")
    return deepcopy(_validate_manifest(root, envelope["payload"], sources))


@_boundary
def read_prepared_workflow_frames(directory, *, sources=None):
    m = prepared_workflow_status(directory, sources=sources)
    frames = []
    for name in m["segments"]:
        for frame in _load_segment(directory, name).frames:
            if frames and frame.step == frames[-1].step:
                compatible_boundary(frames[-1], frame)
                continue
            if frames and frame.step < frames[-1].step:
                raise WorkflowError("Out-of-order frame")
            frames.append(frame)
    return tuple(frames)


@_boundary
def resume_prepared_workflow(directory, *, sources=None, segments=1):
    _segment_count(segments)
    root = Path(directory)
    with storage.writer(root):
        m = prepared_workflow_status(root, sources=sources)
        if m["status"] == "completed":
            return m
        if m["status"] == "stage_failed":
            raise WorkflowError(
                "Failed stage requires inspection; automatic retry is disabled"
            )
        return _advance(root, m, segments, sources)


def _validate_manifest(root, m, sources):
    try:
        expected = {
            "schema",
            "run_id",
            "config",
            "prepared_identity",
            "input_files",
            "configuration",
            "status",
            "stages",
            "files",
            "bundle",
            "segments",
            "checkpoint",
            "accepted_step",
            "production_validated",
            "simulation_readiness",
        }
        if set(m) != expected or m["schema"] != SCHEMA:
            raise ValueError("Manifest schema/checksum mismatch")
        config = PreparedWorkflowConfig.from_dict(m["config"])
        if (
            type(m["input_files"]) is not list
            or any(type(n) is not str for n in m["input_files"])
            or len(set(m["input_files"])) != len(m["input_files"])
            or type(m["run_id"]) is not str
            or not m["run_id"]
        ):
            raise ValueError("Malformed input inventory/run identity")
        if (
            m["status"] in ("ready", "paused", "completed", "budget_exhausted")
            and m["bundle"] is None
        ):
            raise ValueError("Executable state without verified setup")
        if (
            m["production_validated"] is not False
            or m["simulation_readiness"] != "not_established"
            or type(m["accepted_step"]) is not int
            or not 0 <= m["accepted_step"] <= config.total_steps
            or m["status"]
            not in (
                "starting",
                "ready",
                "paused",
                "completed",
                "budget_exhausted",
                "stage_failed",
            )
        ):
            raise ValueError("Invalid workflow status")
        if (
            type(m["stages"]) is not list
            or type(m["segments"]) is not list
            or len(m["segments"]) > config.max_segments
        ):
            raise ValueError("Invalid stage/segment records")
        if len(set(m["segments"])) != len(m["segments"]):
            raise ValueError("Duplicate segment reference")
        for stage in m["stages"]:
            if (
                type(stage) is not dict
                or stage.get("stage")
                not in (
                    "preflight",
                    "load",
                    "minimization",
                    "initialization",
                    "dynamics",
                )
                or stage.get("status")
                not in (
                    "passed",
                    "supplied",
                    "failed",
                    "unavailable",
                    "paused",
                    "completed",
                    "budget_exhausted",
                    "stage_failed",
                )
            ):
                raise ValueError("Malformed stage outcome")
            if "record" in stage and stage["record"] not in m["files"]:
                raise ValueError("Unlisted stage record")
            if "evaluations" in stage and (
                type(stage["evaluations"]) is not int or stage["evaluations"] < 0
            ):
                raise ValueError("Malformed stage call count")
        if sum(s["stage"] == "dynamics" for s in m["stages"]) > config.max_segments:
            raise ValueError("Segment attempt budget exceeded")
        size = 0
        for name, row in m["files"].items():
            if (
                set(row) != {"sha256", "bytes"}
                or type(row["bytes"]) is not int
                or row["bytes"] < 0
            ):
                raise ValueError("Malformed artifact descriptor")
            path = storage.child(root, name)
            if path.is_symlink() or path.stat().st_size != row["bytes"]:
                raise ValueError(f"Artifact size mismatch: {name}")
            size += row["bytes"]
            if (
                size > config.max_artifact_bytes
                or storage.checksum(path.read_bytes()) != row["sha256"]
            ):
                raise ValueError(f"Artifact checksum/storage limit mismatch: {name}")
        for name in (m["configuration"], m["bundle"], m["checkpoint"], *m["segments"]):
            if name is not None and name not in m["files"]:
                raise ValueError("Unlisted operational artifact")
        dynamics_stages = [
            s for s in m["stages"] if s["stage"] == "dynamics" and "record" in s
        ]
        accepted_stage_records = []
        attempts = []
        for stage in dynamics_stages:
            diagnostic = _load_segment(root, stage["record"])
            if (
                stage.get("evaluations") != diagnostic.evaluations
                or stage.get("reason") != diagnostic.termination_reason
            ):
                raise ValueError("Stage accounting differs from signed diagnostic")
            attempts.append((stage["record"], diagnostic))
            is_accepted = stage["record"] in m["segments"]
            expected_status = (
                "completed"
                if is_accepted and diagnostic.final_state.step == config.total_steps
                else "budget_exhausted"
                if diagnostic.termination_reason
                in ("maximum_frames", "maximum_evaluations")
                else "stage_failed"
                if not is_accepted
                else "paused"
            )
            if stage["status"] != expected_status:
                raise ValueError(
                    "Dynamics stage status contradicts diagnostic/progress"
                )
            if is_accepted:
                accepted_stage_records.append(stage["record"])
        if accepted_stage_records != m["segments"]:
            raise ValueError("Accepted history differs from dynamics stage records")
        segments = []
        previous = None
        accepted = None
        for name in m["segments"]:
            segment = _load_segment(root, name)
            segments.append(segment)
            if previous is not None:
                compatible_boundary(previous, segment.frames[0])
            parent = segment.payload["lineage"][-1]["parent_checksum"]
            if parent != (None if accepted is None else accepted.content_checksum):
                raise ValueError("Segment lineage does not extend published boundary")
            if accepted is None and segment.frames[0].step != 0:
                raise ValueError("Workflow must start at step zero")
            previous = segment.final_state
            if segment.payload["diagnostic"][
                "final_verified"
            ] and segment.termination_reason in (
                "completed",
                "maximum_frames",
                "maximum_evaluations",
            ):
                accepted = create_dynamics_checkpoint(segment)
            else:
                raise ValueError("Ineligible diagnostic in accepted segment sequence")
        loaded = _load_inputs(root, m, sources)
        minima = [
            s for s in m["stages"] if s["stage"] == "minimization" and "record" in s
        ]
        if len(minima) > 1:
            raise ValueError("Multiple minimization records")
        for stage in minima:
            diagnostic = _minimum_record(
                storage.decode(storage.read_json(storage.child(root, stage["record"])))
            )
            if (
                diagnostic.options != config.minimization
                or diagnostic.evaluations != stage.get("evaluations")
                or diagnostic.termination_reason != stage.get("reason")
                or (stage["status"] == "passed")
                != (diagnostic.converged and diagnostic.final_evaluation_verified)
            ):
                raise ValueError("Minimum stage/config/diagnostic mismatch")
            from island.minimization.models import system_identity

            original = loaded.system
            if diagnostic.system_fingerprint != system_identity(original) or dict(
                diagnostic.initial_coordinates
            ) != {
                site: tuple(original.coordinates.get(site))
                for site in original.topology.sites
            }:
                raise ValueError("Minimum diagnostic original system mismatch")
            parameter, model = _identities(original, loaded.prepared)
            for evaluation in (
                diagnostic.initial_evaluation,
                diagnostic.final_evaluation,
            ):
                if evaluation is not None and (
                    evaluation.parameter_fingerprint != parameter
                    or evaluation.model_fingerprint != model
                ):
                    raise ValueError("Minimum diagnostic model mismatch")
        if m["bundle"] is not None:
            setup = _load_setup(root, m, sources, loaded=loaded, with_minimum=True)
            system, prepared, initialization, minimum = setup
            if len(minima) != 1 or bundle.record(minimum) != bundle.record(diagnostic):
                raise ValueError("Setup minimum differs from stage record")
            parameter, model = _identities(system, prepared)
            if minimum.final_evaluation.model_fingerprint != model:
                raise ValueError("Minimum model identity differs from prepared model")
            boundary = None
            for name, attempt in attempts:
                validate_prepared_attempt(
                    config, attempt, system, initialization, minimum, boundary
                )
                if name in m["segments"]:
                    boundary = create_dynamics_checkpoint(attempt)
            validate_prepared_trajectory(
                config, segments, system, initialization, minimum, parameter
            )
        elif attempts or segments or m["checkpoint"] is not None:
            raise ValueError("Trajectory without durable setup bundle")
        if m["checkpoint"] is not None:
            saved = load_dynamics_checkpoint(storage.child(root, m["checkpoint"]))
            if saved != accepted or saved.absolute_step != m["accepted_step"]:
                raise ValueError("Manifest checkpoint/progress contradicts segments")
        elif m["accepted_step"] != 0 or accepted is not None:
            raise ValueError("Progress without published checkpoint")
        last = m["stages"][-1] if m["stages"] else None
        expected_status = (
            "starting"
            if last is None
            else (
                "stage_failed"
                if last["status"] in ("failed", "unavailable", "stage_failed")
                else "ready"
                if last["stage"] == "initialization" and last["status"] == "passed"
                else last["status"]
            )
        )
        allowed = {expected_status}
        if (
            expected_status == "paused"
            and sum(s["stage"] == "dynamics" for s in m["stages"])
            >= config.max_segments
        ):
            allowed.add("budget_exhausted")
        if m["status"] not in allowed:
            raise ValueError("Manifest status contradicts recorded stage outcomes")
        if (m["status"] == "completed") != (m["accepted_step"] == config.total_steps):
            raise ValueError("False workflow completion")
        return m
    except WorkflowError:
        raise
    except Exception as error:
        raise WorkflowError(f"Invalid workflow manifest: {error}") from error


def _advance(root, manifest, count, sources):
    config = PreparedWorkflowConfig.from_dict(manifest["config"])
    system, prepared, initialization = _load_setup(root, manifest, sources)
    try:
        evaluator = create_evaluator(system, prepared)
    except Exception as error:  # noqa: BLE001 -- preserve dependency/backend diagnostic
        return _failure(root, manifest, "dynamics", error, sources)
    checkpoint = (
        None
        if manifest["checkpoint"] is None
        else load_dynamics_checkpoint(storage.child(root, manifest["checkpoint"]))
    )
    if checkpoint is not None:
        physical = checkpoint.payload["physical"]
        for name in (
            "temperature_kelvin",
            "friction_per_ps",
            "timestep_fs",
            "thermostat_seed",
        ):
            if physical[name] != getattr(config, name):
                raise WorkflowError(
                    f"Checkpoint/config physical setting differs: {name}"
                )
    for _ in range(count):
        if manifest["accepted_step"] == config.total_steps:
            break
        if (
            sum(s["stage"] == "dynamics" for s in manifest["stages"])
            >= config.max_segments
        ):
            manifest["status"] = "budget_exhausted"
            _publish_manifest(root, manifest, sources)
            break
        steps = min(
            config.segment_steps, config.total_steps - manifest["accepted_step"]
        )
        try:
            with evaluator.open_session() as session:
                if checkpoint is None:
                    result = run_dynamics_segment(
                        system,
                        session,
                        initialization.velocities,
                        config.langevin(steps),
                    )
                else:
                    result = resume_dynamics(
                        checkpoint,
                        system,
                        session,
                        DynamicsSegmentOptions(
                            steps,
                            config.max_evaluations_per_segment,
                            config.max_frames_per_segment,
                            config.recording_interval,
                        ),
                    )
        except Exception as error:  # noqa: BLE001 -- preserve stage outcome, no retry
            return _failure(root, manifest, "dynamics", error, sources)
        # Publish all immutable data before exposing progress. Any exception
        # leaves the previous manifest/boundary authoritative; orphan files ignored.
        candidate = deepcopy(manifest)
        name = _put(
            root,
            candidate,
            "segment",
            storage.json_bytes(
                {"payload": result.payload, "sha256": result.content_checksum}
            ),
        )
        segment_name = name
        eligible = result.payload["diagnostic"][
            "final_verified"
        ] and result.termination_reason in (
            "completed",
            "maximum_frames",
            "maximum_evaluations",
        )
        if eligible:
            candidate["segments"].append(segment_name)
            next_checkpoint = create_dynamics_checkpoint(result)
            name = f"checkpoint-{uuid4().hex}.json"
            raw = storage.json_bytes(
                {
                    "payload": next_checkpoint.payload,
                    "sha256": next_checkpoint.content_checksum,
                }
            )
            if (
                sum(row["bytes"] for row in candidate["files"].values()) + len(raw)
                > config.max_artifact_bytes
            ):
                raise WorkflowError(
                    "Published artifact byte budget exhausted; prior boundary retained"
                )
            save_dynamics_checkpoint(next_checkpoint, storage.child(root, name))
            storage.sync_directory(root)
            candidate["files"][name] = {
                "sha256": storage.checksum(storage.child(root, name).read_bytes()),
                "bytes": storage.child(root, name).stat().st_size,
            }
            candidate["checkpoint"] = name
            candidate["accepted_step"] = next_checkpoint.absolute_step
        candidate["status"] = (
            "completed"
            if candidate["accepted_step"] == config.total_steps
            else "budget_exhausted"
            if result.termination_reason in ("maximum_frames", "maximum_evaluations")
            else "stage_failed"
            if not eligible
            else "paused"
            if result.completed
            else "budget_exhausted"
        )
        candidate["stages"].append(
            {
                "stage": "dynamics",
                "status": candidate["status"],
                "reason": result.termination_reason,
                "evaluations": result.evaluations,
                "record": segment_name,
            }
        )
        _publish_manifest(root, candidate, sources)
        manifest = candidate
        if not eligible or not result.completed:
            break
        checkpoint = next_checkpoint
    return deepcopy(manifest)
