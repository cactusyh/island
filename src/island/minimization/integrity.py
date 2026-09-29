"""Content consistency checks for public records, not computational authenticity."""

from collections.abc import Mapping
from math import isclose, isfinite
from types import MappingProxyType

import numpy as np

from island.evaluation.models import EvaluationResult
from island.exceptions import InvalidMinimizationResultError


def require(condition, message):
    if not condition:
        raise InvalidMinimizationResultError(message)


def number(value, *, nonnegative=False):
    return (
        type(value)
        in (int, float, np.float16, np.float32, np.float64, np.int32, np.int64)
        and isfinite(value)
        and (not nonnegative or value >= 0)
    )


def vectors(mapping, label):
    require(
        isinstance(mapping, Mapping) and bool(mapping),
        f"{label} must be a nonempty mapping",
    )
    for site, xyz in mapping.items():
        require(type(site) is int, f"{label} requires integer stable IDs")
        require(
            type(xyz) is tuple and len(xyz) == 3 and all(number(v) for v in xyz),
            f"{label} requires finite three-component numeric vectors",
        )
    return set(mapping)


def evaluation(record, coords, label):
    from .models import coordinate_hash, force_metrics

    require(type(record) is EvaluationResult, f"{label} must be EvaluationResult")
    require(
        record.coordinate_unit == "angstrom"
        and record.energy_unit == "kJ/mol"
        and record.force_unit == "kJ/(mol*angstrom)",
        f"{label} has unsupported units",
    )
    require(
        record.production_validated is False
        and record.simulation_readiness == "not_established",
        f"{label} has unsupported readiness declarations",
    )
    require(number(record.potential_energy), f"{label} energy must be finite")
    require(
        isinstance(record.forces, MappingProxyType),
        f"{label} forces must be owned immutable data",
    )
    require(
        vectors(record.forces, f"{label} forces") == set(coords),
        f"{label} force-site coverage mismatch",
    )
    require(
        isinstance(record.energy_components, MappingProxyType)
        and bool(record.energy_components)
        and all(
            type(k) is str and k and number(v)
            for k, v in record.energy_components.items()
        ),
        f"{label} components must be finite named values",
    )
    require(
        isclose(
            sum(record.energy_components.values()),
            record.potential_energy,
            rel_tol=1e-10,
            abs_tol=1e-8,
        ),
        f"{label} components contradict total energy",
    )
    for name in (
        "coordinate_fingerprint",
        "parameter_fingerprint",
        "model_fingerprint",
        "evaluation_fingerprint",
        "backend_name",
        "backend_version",
        "platform",
    ):
        require(
            type(getattr(record, name)) is str and bool(getattr(record, name)),
            f"{label} lacks {name}",
        )
    require(
        record.coordinate_fingerprint == coordinate_hash(coords),
        f"{label} coordinate fingerprint mismatch",
    )
    require(
        isinstance(record.settings, MappingProxyType)
        and all(
            type(k) is str
            and type(v) in (str, bool, int, float)
            and (not isinstance(v, (int, float)) or isfinite(v))
            for k, v in record.settings.items()
        ),
        f"{label} settings must be finite immutable scalars",
    )
    require(
        all(number(v, nonnegative=True) for v in force_metrics(record)),
        f"{label} has invalid force metrics",
    )


def validate(result):
    from .models import (
        MinimizationOptions,
        MinimizationStep,
        coordinate_hash,
        force_metrics,
    )

    ids = vectors(result.initial_coordinates, "Initial coordinates")
    require(
        vectors(result.coordinates, "Returned coordinates") == ids,
        "Initial/returned site coverage differs",
    )
    require(
        isinstance(result.initial_coordinates, MappingProxyType)
        and isinstance(result.coordinates, MappingProxyType),
        "Coordinates must be owned immutable mappings",
    )
    require(
        type(result.options) is MinimizationOptions, "Invalid minimization options type"
    )
    result.options.__post_init__()
    for name in (
        "force_tolerance",
        "energy_change_tolerance",
        "energy_increase_tolerance",
    ):
        require(
            number(getattr(result.options, name)),
            f"Options {name} must be an immutable numeric scalar",
        )
    require(
        type(result.iterations) is int
        and 0 <= result.iterations <= result.options.max_iterations,
        "Invalid iteration count",
    )
    require(
        type(result.evaluations) is int
        and 1 <= result.evaluations <= result.options.max_evaluations,
        "Invalid evaluation count",
    )
    for name in ("converged", "final_evaluation_verified"):
        require(type(getattr(result, name)) is bool, f"{name} must be boolean")
    require(
        result.optimizer_success is None or type(result.optimizer_success) is bool,
        "Invalid optimizer success flag",
    )
    for name in (
        "optimizer_name",
        "optimizer_version",
        "optimizer_message",
        "system_fingerprint",
    ):
        require(
            type(getattr(result, name)) is str and bool(getattr(result, name)),
            f"Invalid {name}",
        )
    require(result.optimizer_name == "SciPy L-BFGS-B", "Unsupported optimizer")
    require(
        result.coordinate_unit == "angstrom"
        and result.energy_unit == "kJ/mol"
        and result.force_unit == "kJ/(mol*angstrom)",
        "Unsupported result units",
    )
    require(
        result.production_validated is False
        and result.simulation_readiness == "not_established",
        "Unsupported readiness declarations",
    )
    require(
        type(result.termination_reason) is str
        and result.termination_reason
        in {
            "force_converged",
            "maximum_iterations",
            "maximum_evaluations",
            "energy_stagnation",
            "line_search_failed",
            "invalid_geometry",
            "evaluation_failed",
            "stereochemistry_changed",
        },
        "Invalid termination reason",
    )
    require(
        type(result.stereochemistry) is str
        and result.stereochemistry in {"passed", "not_assigned"},
        "Invalid stereochemistry outcome",
    )
    require(
        type(result.returned_state) is str
        and result.returned_state
        in {
            "accepted_iterate",
            "best_admissible_evaluated_trial",
            "unevaluated_initial_coordinates",
        },
        "Invalid returned state",
    )
    require(type(result.history) is tuple, "History must be immutable")
    initial, final = result.initial_evaluation, result.final_evaluation
    for record, coords, label in (
        (initial, result.initial_coordinates, "Initial evaluation"),
        (final, result.coordinates, "Final evaluation"),
    ):
        if record is not None:
            evaluation(record, coords, label)
    if initial is None:
        require(
            final is None
            and not result.history
            and result.iterations == 0
            and result.evaluations == 1
            and not result.converged
            and not result.final_evaluation_verified
            and result.optimizer_success is None
            and result.returned_state == "unevaluated_initial_coordinates"
            and result.termination_reason in {"evaluation_failed", "invalid_geometry"}
            and result.coordinates == result.initial_coordinates,
            "Inconsistent initial-failure diagnostic",
        )
        return
    require(
        final is not None,
        "An evaluated initial state requires an available diagnostic final state",
    )
    require(
        initial.model_fingerprint == final.model_fingerprint
        and initial.parameter_fingerprint == final.parameter_fingerprint,
        "Evaluation identities differ",
    )
    require(
        result.returned_state != "unevaluated_initial_coordinates",
        "Evaluated result cannot claim unevaluated coordinates",
    )
    require(
        len(result.history) == result.iterations + 1,
        "History does not cover accepted iterations",
    )
    previous_calls = 0
    for index, step in enumerate(result.history):
        require(type(step) is MinimizationStep, "Invalid history record type")
        require(
            type(step.iteration) is int and step.iteration == index,
            "Invalid history iteration",
        )
        require(
            type(step.evaluations) is int
            and previous_calls < step.evaluations <= result.evaluations,
            "Invalid history evaluation count",
        )
        previous_calls = step.evaluations
        require(
            number(step.energy)
            and number(step.fmax, nonnegative=True)
            and number(step.rms_force, nonnegative=True)
            and step.rms_force <= step.fmax + 1e-12,
            "Invalid history energy/force metrics",
        )
        require(
            type(step.coordinate_fingerprint) is str
            and bool(step.coordinate_fingerprint),
            "Invalid history fingerprint",
        )
    first = result.history[0]
    require(
        first.evaluations == 1
        and first.energy == initial.potential_energy
        and (first.fmax, first.rms_force) == force_metrics(initial)
        and first.coordinate_fingerprint == coordinate_hash(result.initial_coordinates),
        "Initial history contradicts initial evaluation",
    )
    if result.final_evaluation_verified:
        require(
            result.evaluations >= 2
            and result.history[-1].evaluations < result.evaluations,
            "Final verification requires an additional evaluation",
        )
    else:
        require(
            not result.converged
            and result.termination_reason
            in {
                "maximum_evaluations",
                "evaluation_failed",
                "invalid_geometry",
                "stereochemistry_changed",
            },
            "Unverified final state must remain a failure diagnostic",
        )
    if not result.final_evaluation_verified:
        require(
            result.returned_state == "best_admissible_evaluated_trial",
            "Unverified result must identify diagnostic coordinates",
        )
    if result.returned_state == "accepted_iterate":
        require(
            result.termination_reason in {"force_converged", "energy_stagnation"},
            "Accepted-iterate status contradicts termination",
        )
        require(
            result.history[-1].coordinate_fingerprint
            == coordinate_hash(result.coordinates),
            "Returned accepted iterate does not match history",
        )
    require(
        result.converged == (result.termination_reason == "force_converged"),
        "Convergence and termination disagree",
    )
    if result.termination_reason == "maximum_iterations":
        require(
            result.iterations == result.options.max_iterations,
            "Iteration limit was not reached",
        )
    if result.termination_reason == "maximum_evaluations":
        require(
            result.evaluations == result.options.max_evaluations,
            "Evaluation limit was not reached",
        )
    if result.converged:
        require(
            result.final_evaluation_verified
            and force_metrics(final)[0] <= result.options.force_tolerance
            and final.potential_energy
            <= initial.potential_energy + result.options.energy_increase_tolerance,
            "Claimed convergence fails final force/energy/verification criteria",
        )
