"""Bounded BAOAB propagation with synchronized endpoints and private random streams."""

from copy import deepcopy

import numpy as np

from island.chemistry.coordinate_stereo import (
    assigned_cip_labels,
    validate_coordinate_stereochemistry,
)
from island.core import AtomSite, MolecularSystem
from island.core.coordinate_provenance import coordinate_hash
from island.evaluation import OpenMMBoundPotential, PotentialEvaluator
from island.exceptions import (
    DynamicsInputError,
    DynamicsUnavailableError,
    EvaluationInputError,
    EvaluationUnavailableError,
    IslandError,
    StereochemistryError,
)
from island.minimization.integrity import number
from island.minimization.models import system_identity

from .integrity import checked_evaluation, retained_count, same_model, verify_final
from .langevin_models import (
    LangevinFrame,
    LangevinOptions,
    LangevinResult,
    langevin_identity,
)
from .models import (
    kinetic_energy,
    owned_vectors,
    velocity_hash,
)
from .thermal import GAS_CONSTANT, _normal_increments, _rng


class _Stop(Exception):
    def __init__(self, reason, message, details=None):
        super().__init__(message)
        self.reason, self.details = reason, details or {}


def run_langevin(
    system, evaluator, velocities, options, *, velocity_unit="angstrom/ps"
):
    """Advance full-step states without altering the bound potential or inputs.

    Budgets include initial and final checks. Any failed trial leaves both
    coordinates and synchronized velocities at the last accepted state.
    """
    if (
        type(options) is not LangevinOptions
        or not isinstance(system, MolecularSystem)
        or not isinstance(evaluator, PotentialEvaluator)
    ):
        raise DynamicsInputError(
            "Supply MolecularSystem, PotentialEvaluator and LangevinOptions"
        )
    if velocity_unit != "angstrom/ps":
        raise DynamicsInputError("Initial velocities must be explicitly in angstrom/ps")
    try:
        options.__post_init__()
        owned = deepcopy(system)
        owned.validate()
        if (
            owned.representation != "atomistic"
            or not owned.topology.sites
            or owned.box is not None
        ):
            raise DynamicsInputError(
                "Require nonempty atomistic system without a simulation cell"
            )
        if any(type(site) is not AtomSite for site in owned.topology.sites.values()):
            raise DynamicsInputError("Dynamics supports atomistic sites only")
        ids = tuple(sorted(owned.topology.sites))
        masses = {site: owned.topology.sites[site].mass for site in ids}
        if any(
            type(site) is not int or not number(m) or m <= 0
            for site, m in masses.items()
        ):
            raise DynamicsInputError("Require positive finite atomic masses in dalton")
        masses = {site: float(mass) for site, mass in masses.items()}
        coords = owned_vectors({site: owned.coordinates.get(site) for site in ids})
        initial_velocities = owned_vectors(velocities)
        if set(initial_velocities) != set(ids):
            raise DynamicsInputError("Velocity coverage must exactly match stable IDs")
        initial_kinetic = kinetic_energy(masses, initial_velocities)
        if not number(initial_kinetic) or not number(
            2 * initial_kinetic / (3 * len(ids) * GAS_CONSTANT)
        ):
            raise DynamicsInputError("Initial kinetic energy must be finite")
        if isinstance(evaluator, OpenMMBoundPotential):
            evaluator.validate_system(owned)
        identity = system_identity(owned)
        expected = assigned_cip_labels(owned)
        validate_coordinate_stereochemistry(
            owned, coords, expected, stage="initial dynamics coordinates"
        )
    except ImportError as error:
        raise DynamicsUnavailableError(
            "RDKit is required for assigned tetrahedral stereochemistry"
        ) from error
    except EvaluationUnavailableError as error:
        raise DynamicsUnavailableError(str(error)) from error
    except (
        IslandError,
        TypeError,
        ValueError,
        AttributeError,
        KeyError,
        OverflowError,
        FloatingPointError,
    ) as error:
        raise DynamicsInputError(str(error)) from error
    dt = (
        options.timestep_fs * 0.001
    )  # fs -> ps, once; public forces already per angstrom
    mass_array = np.array([masses[site] for site in ids])[:, None]
    evaluations = 0
    initial_record = None
    rng = _rng(options.thermostat_seed)
    random_steps = 0
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            gamma_dt = options.friction_per_ps * dt
            decay = np.exp(-gamma_dt)
            variance_fraction = -np.expm1(-2 * gamma_dt)
            noise_scale = np.sqrt(
                variance_fraction
                * (100 * GAS_CONSTANT * options.temperature_kelvin / mass_array)
            )
        if not np.isfinite(noise_scale).all():
            raise ValueError("Nonfinite thermal velocity scale")
    except (FloatingPointError, OverflowError, ValueError) as error:
        raise DynamicsInputError(f"Unsafe thermostat coefficients: {error}") from error

    def evaluate(coordinates, *, final=False):
        nonlocal evaluations
        evaluations += 1  # includes failures; capacity is checked by caller
        try:
            method = (
                evaluator.evaluate_fresh
                if final and isinstance(evaluator, OpenMMBoundPotential)
                else evaluator.evaluate
            )
            record = method(dict(coordinates), coordinate_unit="angstrom")
        except EvaluationUnavailableError as error:
            raise DynamicsUnavailableError(str(error)) from error
        except EvaluationInputError as error:
            raise _Stop(
                "invalid_geometry", str(error) or type(error).__name__
            ) from error
        except Exception as error:
            raise _Stop(
                "evaluation_failed", str(error) or type(error).__name__
            ) from error
        try:
            checked_evaluation(record, coordinates)
            if initial_record is not None:
                same_model(initial_record, record)
            return record
        except (IslandError, TypeError, ValueError, AttributeError) as error:
            raise _Stop(
                "evaluation_failed", str(error) or type(error).__name__
            ) from error

    def frame(step, coordinates, speeds, record, initial_total):
        kinetic = kinetic_energy(masses, speeds)
        total = None if record is None else record.potential_energy + kinetic
        deviation = (
            None
            if total is None
            else (0.0 if initial_total is None else total - initial_total)
        )
        if (
            not number(kinetic)
            or not number(2 * kinetic / (3 * len(ids) * GAS_CONSTANT))
            or (total is not None and (not number(total) or not number(deviation)))
        ):
            raise _Stop(
                "invalid_geometry", "Nonfinite synchronized kinetic/total energy"
            )
        return LangevinFrame(
            step,
            step * dt,
            coordinates,
            speeds,
            record,
            kinetic,
            total,
            deviation,
            coordinate_hash(coordinates),
            velocity_hash(speeds),
            3 * len(ids),
            2 * kinetic / (3 * len(ids) * GAS_CONSTANT),
        )

    current = frame(0, coords, initial_velocities, None, None)
    stored = [current]
    failed_calls = 0
    reason, message, failure_stage, attempted = (
        "completed",
        "Requested BAOAB steps completed; no equilibration claim",
        None,
        None,
    )
    details = {}
    maximum_deviation = None
    try:
        initial_record = evaluate(coords)
        current = frame(0, coords, initial_velocities, initial_record, None)
        stored = [current]
        maximum_deviation = 0.0
    except _Stop as error:
        reason, message, details = error.reason, str(error), error.details
        failure_stage, attempted = "initial_evaluation", 0
    if current.evaluation is not None:
        initial_total = current.total_energy
        while current.step < options.steps:
            next_step = current.step + 1
            if evaluations >= options.max_evaluations - 1:
                reason, message, failure_stage, attempted = (
                    "maximum_evaluations",
                    "Evaluator-call budget reached; final check reserved",
                    "budget",
                    next_step,
                )
                break
            if (
                retained_count(next_step, options.recording_interval)
                > options.max_frames
            ):
                reason, message, failure_stage, attempted = (
                    "maximum_frames",
                    "Stored-frame budget prevents another step",
                    "budget",
                    next_step,
                )
                break
            before = evaluations
            random_steps += 1
            try:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    noise = _normal_increments(rng, (len(ids), 3))
                    x = np.array([current.coordinates[site] for site in ids])
                    v = np.array([current.velocities[site] for site in ids])
                    forces = np.array([current.evaluation.forces[site] for site in ids])
                    half = v + 0.5 * dt * (100.0 * forces / mass_array)
                    midpoint = x + 0.5 * dt * half
                    thermal_velocity = decay * half + noise_scale * noise
                    trial_coordinates = owned_vectors(
                        dict(
                            zip(
                                ids, midpoint + 0.5 * dt * thermal_velocity, strict=True
                            )
                        )
                    )
                    validate_coordinate_stereochemistry(
                        owned,
                        trial_coordinates,
                        expected,
                        stage=f"dynamics trial step {next_step}",
                    )
                    trial_record = evaluate(trial_coordinates)
                    next_forces = np.array([trial_record.forces[site] for site in ids])
                    next_velocities = owned_vectors(
                        dict(
                            zip(
                                ids,
                                thermal_velocity
                                + 0.5 * dt * (100.0 * next_forces / mass_array),
                                strict=True,
                            )
                        )
                    )
                    trial = frame(
                        next_step,
                        trial_coordinates,
                        next_velocities,
                        trial_record,
                        initial_total,
                    )
                current = trial  # commit only after synchronized velocities and all checks pass
                maximum_deviation = max(
                    maximum_deviation, abs(current.energy_deviation)
                )
                if current.step % options.recording_interval == 0:
                    stored.append(current)
            except DynamicsUnavailableError:
                raise
            except ImportError as error:
                raise DynamicsUnavailableError(
                    "RDKit is required for assigned tetrahedral stereochemistry"
                ) from error
            except (
                StereochemistryError,
                _Stop,
                IslandError,
                FloatingPointError,
                ValueError,
                OverflowError,
            ) as error:
                reason = (
                    error.reason
                    if isinstance(error, _Stop)
                    else "stereochemistry_changed"
                    if isinstance(error, StereochemistryError)
                    else "invalid_geometry"
                )
                message, details = (
                    str(error),
                    error.details if isinstance(error, _Stop) else {},
                )
                failure_stage, attempted = "trial", next_step
                failed_calls = evaluations - before
                break
    if stored[-1].step != current.step:
        stored.append(current)
    final_attempted, verified, final_record, final_error = False, False, None, None
    if current.evaluation is not None:
        if evaluations < options.max_evaluations:
            final_attempted = True
            try:
                checked = evaluate(current.coordinates, final=True)
                verify_final(current.evaluation, checked)
                final_record, verified = checked, True
            except DynamicsUnavailableError:
                raise
            except (IslandError, _Stop) as error:
                final_error = str(error) or type(error).__name__
        else:
            final_error = "No evaluator capacity for independent final check"
        if not verified and reason == "completed":
            reason, message, failure_stage, attempted = (
                "final_evaluation_failed",
                final_error,
                "final_evaluation",
                current.step,
            )
    result = LangevinResult(
        options,
        masses,
        tuple(stored),
        current.step,
        evaluations,
        failed_calls,
        final_attempted,
        verified,
        final_record,
        reason == "completed",
        reason,
        attempted,
        failure_stage,
        message,
        details,
        final_error,
        maximum_deviation,
        identity,
        langevin_identity(options, masses, stored[0], identity, np.__version__),
        "passed" if expected else "not_assigned",
        normal_draws=random_steps * 3 * len(ids),
        random_steps=random_steps,
        numpy_version=np.__version__,
    )
    result.validate_integrity()
    return result
