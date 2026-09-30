"""Validate state content and bookkeeping; this is not computational authenticity."""

from math import isclose
from types import MappingProxyType

import numpy as np

from island.core.coordinate_provenance import coordinate_hash
from island.exceptions import InvalidDynamicsResultError, IslandError
from island.minimization.integrity import evaluation, number, vectors


def require(condition, message):
    if not condition:
        raise InvalidDynamicsResultError(message)


def checked_evaluation(record, coordinates):
    try:
        evaluation(record, coordinates, "Dynamics evaluation")
    except (
        IslandError,
        AttributeError,
        KeyError,
        TypeError,
        ValueError,
        OverflowError,
        FloatingPointError,
    ) as error:
        raise InvalidDynamicsResultError(str(error)) from error


def same_model(left, right):
    require(
        left.model_fingerprint == right.model_fingerprint
        and left.parameter_fingerprint == right.parameter_fingerprint,
        "Dynamics model/parameter identity changed",
    )
    require(
        (left.backend_name, left.backend_version, left.platform, dict(left.settings))
        == (
            right.backend_name,
            right.backend_version,
            right.platform,
            dict(right.settings),
        ),
        "Dynamics backend/settings changed",
    )


def verify_final(accepted, checked):
    same_model(accepted, checked)
    require(
        isclose(
            accepted.potential_energy,
            checked.potential_energy,
            rel_tol=1e-10,
            abs_tol=1e-8,
        ),
        "Independent final energy differs from accepted state",
    )
    for site in accepted.forces:
        require(
            np.allclose(
                accepted.forces[site], checked.forces[site], rtol=1e-10, atol=1e-8
            ),
            "Independent final forces differ from accepted state",
        )
    require(
        set(accepted.energy_components) == set(checked.energy_components),
        "Final energy component names changed",
    )
    for key in accepted.energy_components:
        require(
            isclose(
                accepted.energy_components[key],
                checked.energy_components[key],
                rel_tol=1e-10,
                abs_tol=1e-8,
            ),
            "Independent final component differs from accepted state",
        )


def validate_frame_content(frame):
    from .langevin_models import LangevinFrame
    from .models import DynamicsFrame, velocity_hash
    from .thermal import GAS_CONSTANT

    try:
        require(
            type(frame) in (DynamicsFrame, LangevinFrame), "Expected dynamics frame"
        )
        if type(frame) is LangevinFrame:
            require(
                type(frame.degrees_of_freedom) is int
                and frame.degrees_of_freedom > 0
                and frame.degrees_of_freedom == 3 * len(frame.coordinates),
                "DOF must be 3N",
            )
            require(
                frame.temperature_unit == "kelvin"
                and number(frame.instantaneous_temperature_kelvin, nonnegative=True)
                and isclose(
                    frame.instantaneous_temperature_kelvin,
                    2
                    * frame.kinetic_energy
                    / (frame.degrees_of_freedom * GAS_CONSTANT),
                    rel_tol=1e-12,
                    abs_tol=1e-10,
                ),
                "Incorrect kinetic temperature",
            )
        require(type(frame.step) is int and frame.step >= 0, "Invalid frame step")
        require(number(frame.time_ps, nonnegative=True), "Invalid frame time")
        require(
            (
                frame.coordinate_unit,
                frame.velocity_unit,
                frame.energy_unit,
                frame.time_unit,
            )
            == ("angstrom", "angstrom/ps", "kJ/mol", "ps"),
            "Unsupported frame units",
        )
        require(
            isinstance(frame.coordinates, MappingProxyType)
            and isinstance(frame.velocities, MappingProxyType),
            "Frame vectors must be owned immutable mappings",
        )
        require(
            vectors(frame.coordinates, "Coordinates")
            == vectors(frame.velocities, "Velocities"),
            "Coordinate/velocity coverage differs",
        )
        require(
            type(frame.coordinate_fingerprint) is str
            and frame.coordinate_fingerprint == coordinate_hash(frame.coordinates),
            "Coordinate fingerprint mismatch",
        )
        require(
            type(frame.velocity_fingerprint) is str
            and frame.velocity_fingerprint == velocity_hash(frame.velocities),
            "Velocity fingerprint mismatch",
        )
        require(
            number(frame.kinetic_energy, nonnegative=True), "Invalid kinetic energy"
        )
        if frame.evaluation is None:
            require(
                frame.step == 0
                and frame.time_ps == 0
                and frame.total_energy is None
                and frame.energy_deviation is None,
                "Unevaluated frame must be initial diagnostic",
            )
        else:
            checked_evaluation(frame.evaluation, frame.coordinates)
            require(
                number(frame.total_energy)
                and isclose(
                    frame.total_energy,
                    frame.evaluation.potential_energy + frame.kinetic_energy,
                    rel_tol=1e-12,
                    abs_tol=1e-10,
                ),
                "Incorrect total energy arithmetic",
            )
            require(number(frame.energy_deviation), "Invalid total-energy deviation")
    except InvalidDynamicsResultError:
        raise
    except (
        IslandError,
        AttributeError,
        KeyError,
        TypeError,
        ValueError,
        OverflowError,
        FloatingPointError,
    ) as error:
        raise InvalidDynamicsResultError(f"Malformed frame: {error}") from error


def retained_count(completed, interval):
    return completed // interval + 1 + int(completed % interval != 0)


def retained_steps(completed, interval):
    return sorted(set(range(0, completed + 1, interval)) | {completed})


def validate_result(result):
    from .langevin_models import (
        BAOAB,
        LangevinFrame,
        LangevinOptions,
        LangevinResult,
        langevin_identity,
    )
    from .models import (
        INTEGRATOR,
        DynamicsFrame,
        DynamicsOptions,
        DynamicsResult,
        dynamics_identity,
        kinetic_energy,
    )
    from .thermal import RNG_ALGORITHM

    thermal = type(result) is LangevinResult
    require(type(result) in (DynamicsResult, LangevinResult), "Invalid result type")
    require(
        type(result.options) is (LangevinOptions if thermal else DynamicsOptions),
        "Invalid dynamics options",
    )
    result.options.__post_init__()
    require(
        result.integrator == (BAOAB if thermal else INTEGRATOR)
        and result.mass_unit == "dalton",
        "Unsupported integrator/mass units",
    )
    require(
        result.production_validated is False
        and result.simulation_readiness == "not_established",
        "Unsupported readiness declarations",
    )
    require(
        isinstance(result.masses, MappingProxyType) and bool(result.masses),
        "Masses must be an owned nonempty mapping",
    )
    require(
        all(
            type(site) is int and number(mass) and mass > 0
            for site, mass in result.masses.items()
        ),
        "Masses require integer IDs and positive finite values",
    )
    require(
        type(result.completed_steps) is int
        and 0 <= result.completed_steps <= result.options.steps,
        "Invalid completed steps",
    )
    require(
        type(result.frames) is tuple
        and 1 <= len(result.frames) <= result.options.max_frames,
        "Invalid frame storage count",
    )
    for frame in result.frames:
        require(
            type(frame) is (LangevinFrame if thermal else DynamicsFrame),
            "Frame/integrator mismatch",
        )
        validate_frame_content(frame)
        require(
            set(frame.coordinates) == set(result.masses), "Mass/state coverage differs"
        )
        require(
            isclose(
                frame.time_ps,
                frame.step * (result.options.timestep_fs * 0.001),
                rel_tol=0,
                abs_tol=1e-14,
            ),
            "Frame time differs from fixed timestep",
        )
        require(
            isclose(
                frame.kinetic_energy,
                kinetic_energy(result.masses, frame.velocities),
                rel_tol=1e-12,
                abs_tol=1e-10,
            ),
            "Kinetic energy differs from masses and full-step velocities",
        )
    require(
        len(result.frames)
        == retained_count(result.completed_steps, result.options.recording_interval),
        "Stored frame count contradicts schedule",
    )
    require(
        [f.step for f in result.frames]
        == retained_steps(result.completed_steps, result.options.recording_interval),
        "Frames must include initial, interval samples and final without duplicates",
    )
    initial, final = result.frames[0], result.frames[-1]
    for field in (
        "completed",
        "final_evaluation_attempted",
        "final_evaluation_verified",
    ):
        require(type(getattr(result, field)) is bool, f"{field} must be boolean")
    require(
        type(result.failed_trial_evaluations) is int
        and result.failed_trial_evaluations in (0, 1),
        "Invalid failed trial count",
    )
    require(
        type(result.evaluations) is int
        and 1 <= result.evaluations <= result.options.max_evaluations,
        "Invalid evaluator-call count",
    )
    require(
        result.evaluations
        == 1
        + result.completed_steps
        + result.failed_trial_evaluations
        + int(result.final_evaluation_attempted),
        "Evaluator-call accounting inconsistent",
    )
    require(
        type(result.message) is str and bool(result.message),
        "Missing termination message",
    )
    require(
        type(result.system_fingerprint) is str and bool(result.system_fingerprint),
        "Missing system identity",
    )
    require(
        type(result.stereochemistry) is str
        and result.stereochemistry in ("passed", "not_assigned"),
        "Invalid stereo declaration",
    )
    require(
        isinstance(result.failure_details, MappingProxyType)
        and all(
            type(k) is str and (type(v) in (str, bool) or number(v))
            for k, v in result.failure_details.items()
        ),
        "Diagnostics must be immutable finite scalar data",
    )
    require(
        result.final_check_error is None
        or (type(result.final_check_error) is str and bool(result.final_check_error)),
        "Invalid final-check diagnostic",
    )
    require(
        type(result.termination_reason) is str
        and result.termination_reason
        in {
            "completed",
            "maximum_evaluations",
            "maximum_frames",
            "evaluation_failed",
            "invalid_geometry",
            "stereochemistry_changed",
            "energy_guard_exceeded",
            "final_evaluation_failed",
        },
        "Invalid termination reason",
    )
    require(
        result.completed == (result.termination_reason == "completed"),
        "Completion flag contradicts termination",
    )
    require(
        result.failure_stage is None
        or (
            type(result.failure_stage) is str
            and result.failure_stage
            in {"initial_evaluation", "trial", "budget", "final_evaluation"}
        ),
        "Invalid failure stage",
    )
    if result.completed:
        require(
            result.completed_steps == result.options.steps
            and result.final_evaluation_verified
            and result.attempted_step is None
            and result.failure_stage is None
            and not result.failure_details
            and result.failed_trial_evaluations == 0,
            "Completion requires requested steps and independent final verification",
        )
    else:
        require(
            type(result.attempted_step) is int and result.failure_stage is not None,
            "Partial result requires attempted step and failure stage",
        )
        expected = (
            0
            if result.failure_stage == "initial_evaluation"
            else result.completed_steps
            if result.failure_stage == "final_evaluation"
            else result.completed_steps + 1
        )
        require(
            result.attempted_step == expected
            and 0 <= result.attempted_step <= result.options.steps,
            "Attempted step is inconsistent",
        )
    stages = {
        "completed": {None},
        "maximum_evaluations": {"budget"},
        "maximum_frames": {"budget"},
        "evaluation_failed": {"initial_evaluation", "trial"},
        "invalid_geometry": {"initial_evaluation", "trial"},
        "stereochemistry_changed": {"trial"},
        "energy_guard_exceeded": {"trial"},
        "final_evaluation_failed": {"final_evaluation"},
    }
    require(
        result.failure_stage in stages[result.termination_reason],
        "Termination contradicts failure stage",
    )
    require(
        not result.failed_trial_evaluations or result.failure_stage == "trial",
        "Rejected evaluator call must belong to a failed trial",
    )
    if result.termination_reason == "final_evaluation_failed":
        require(
            not result.final_evaluation_verified
            and result.completed_steps == result.options.steps,
            "Final-check failure must follow all requested steps",
        )
    if initial.evaluation is not None and not result.final_evaluation_attempted:
        require(
            result.evaluations == result.options.max_evaluations == 1
            and result.completed_steps == 0,
            "Reserved final check was omitted despite available capacity",
        )
    if result.termination_reason in {"maximum_evaluations", "maximum_frames"}:
        require(
            result.failure_stage == "budget" and result.failed_trial_evaluations == 0,
            "Budget stop cannot commit an attempted trial",
        )
    if result.termination_reason == "maximum_evaluations":
        require(
            result.evaluations == result.options.max_evaluations,
            "Evaluator budget not exhausted",
        )
    if result.termination_reason == "maximum_frames":
        require(
            retained_count(
                result.completed_steps + 1, result.options.recording_interval
            )
            > result.options.max_frames,
            "Storage limit did not prevent next step",
        )
    if initial.evaluation is None:
        require(
            result.completed_steps == 0
            and not result.final_evaluation_attempted
            and not result.final_evaluation_verified
            and result.final_evaluation is None
            and result.evaluations == 1
            and result.failure_stage == "initial_evaluation"
            and result.max_abs_energy_deviation is None
            and result.termination_reason in {"evaluation_failed", "invalid_geometry"},
            "Invalid initial-failure diagnostic",
        )
    else:
        require(
            number(result.max_abs_energy_deviation, nonnegative=True)
            and (
                thermal
                or result.max_abs_energy_deviation
                <= result.options.max_energy_deviation
            ),
            "Invalid maximum accepted energy deviation",
        )
        observed = []
        for frame in result.frames:
            require(frame.evaluation is not None, "Accepted frames require evaluations")
            same_model(initial.evaluation, frame.evaluation)
            require(
                isclose(
                    frame.energy_deviation,
                    frame.total_energy - initial.total_energy,
                    rel_tol=1e-12,
                    abs_tol=1e-10,
                ),
                "Incorrect energy deviation",
            )
            require(
                thermal
                or abs(frame.energy_deviation) <= result.options.max_energy_deviation,
                "Accepted frame violates energy guard",
            )
            observed.append(abs(frame.energy_deviation))
        require(
            result.max_abs_energy_deviation >= max(observed),
            "Maximum deviation omits a stored frame",
        )
        if len(result.frames) == result.completed_steps + 1:
            require(
                result.max_abs_energy_deviation == max(observed),
                "Maximum deviation contradicts complete trajectory",
            )
    if result.final_evaluation_verified:
        require(
            result.final_evaluation_attempted
            and result.final_evaluation is not None
            and result.final_check_error is None
            and final.evaluation is not None,
            "Missing verified final record",
        )
        checked_evaluation(result.final_evaluation, final.coordinates)
        verify_final(final.evaluation, result.final_evaluation)
    else:
        require(
            not result.completed and result.final_evaluation is None,
            "Unverified final state cannot claim completion",
        )
        if initial.evaluation is not None:
            require(
                result.final_check_error is not None,
                "Missing final-check failure diagnostic",
            )
    if thermal:
        require(
            result.termination_reason != "energy_guard_exceeded",
            "Langevin has no NVE conservation guard",
        )
        require(
            result.rng_algorithm == RNG_ALGORITHM
            and type(result.numpy_version) is str
            and bool(result.numpy_version),
            "Invalid RNG identity",
        )
        require(
            type(result.random_steps) is int
            and result.random_steps
            == result.completed_steps + int(result.failure_stage == "trial"),
            "Random step accounting differs",
        )
        require(
            type(result.normal_draws) is int
            and result.normal_draws == result.random_steps * 3 * len(result.masses),
            "Normal draw accounting differs",
        )
    require(
        type(result.dynamics_fingerprint) is str
        and result.dynamics_fingerprint
        == (
            langevin_identity(
                result.options,
                result.masses,
                initial,
                result.system_fingerprint,
                result.numpy_version,
            )
            if thermal
            else dynamics_identity(
                result.options, result.masses, initial, result.system_fingerprint
            )
        ),
        "Dynamics fingerprint mismatch",
    )
