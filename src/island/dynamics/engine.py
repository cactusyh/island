"""Bounded fixed-step velocity Verlet; trial states are committed atomically."""

from copy import deepcopy

import numpy as np

from island.chemistry.coordinate_stereo import (
    assigned_cip_labels,
    validate_coordinate_stereochemistry,
)
from island.core import AtomSite, MolecularSystem
from island.core.coordinate_provenance import coordinate_hash
from island.evaluation import OpenMMSinglePointEvaluator, PotentialEvaluator
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
from .models import (
    DynamicsFrame,
    DynamicsOptions,
    DynamicsResult,
    dynamics_identity,
    kinetic_energy,
    owned_vectors,
    velocity_hash,
)


class _Stop(Exception):
    def __init__(self, reason, message, details=None):
        super().__init__(message)
        self.reason, self.details = reason, details or {}


def run_nve(system, evaluator, velocities, options, *, velocity_unit="angstrom/ps"):
    """Advance full-step states without altering the bound potential or inputs.

    Budgets include initial and final checks. Any failed trial leaves both
    coordinates and synchronized velocities at the last accepted state.
    """
    if (
        type(options) is not DynamicsOptions
        or not isinstance(system, MolecularSystem)
        or not isinstance(evaluator, PotentialEvaluator)
    ):
        raise DynamicsInputError(
            "Supply MolecularSystem, PotentialEvaluator and DynamicsOptions"
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
        if not number(initial_kinetic):
            raise DynamicsInputError("Initial kinetic energy must be finite")
        if isinstance(evaluator, OpenMMSinglePointEvaluator):
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

    def evaluate(coordinates):
        nonlocal evaluations
        evaluations += 1  # includes failures; capacity is checked by caller
        try:
            record = evaluator.evaluate(dict(coordinates), coordinate_unit="angstrom")
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
        if not number(kinetic) or (
            total is not None and (not number(total) or not number(deviation))
        ):
            raise _Stop(
                "invalid_geometry", "Nonfinite synchronized kinetic/total energy"
            )
        return DynamicsFrame(
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
        )

    current = frame(0, coords, initial_velocities, None, None)
    stored = [current]
    failed_calls = 0
    reason, message, failure_stage, attempted = (
        "completed",
        "Requested fixed steps completed",
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
            try:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    x = np.array([current.coordinates[site] for site in ids])
                    v = np.array([current.velocities[site] for site in ids])
                    forces = np.array([current.evaluation.forces[site] for site in ids])
                    half = v + 0.5 * dt * (100.0 * forces / mass_array)
                    trial_coordinates = owned_vectors(
                        dict(zip(ids, x + dt * half, strict=True))
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
                                half + 0.5 * dt * (100.0 * next_forces / mass_array),
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
                if abs(trial.energy_deviation) > options.max_energy_deviation:
                    raise _Stop(
                        "energy_guard_exceeded",
                        "Absolute total-energy deviation exceeds declared guard",
                        {
                            "trial_total_energy": trial.total_energy,
                            "trial_energy_deviation": trial.energy_deviation,
                            "guard_kj_mol": options.max_energy_deviation,
                        },
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
                checked = evaluate(current.coordinates)
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
    result = DynamicsResult(
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
        dynamics_identity(options, masses, stored[0], identity),
        "passed" if expected else "not_assigned",
    )
    result.validate_integrity()
    return result
