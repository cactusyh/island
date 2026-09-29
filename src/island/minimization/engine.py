"""Strictly budgeted Cartesian L-BFGS-B using a bound PotentialEvaluator."""

from copy import deepcopy

import numpy as np

from island.chemistry.coordinate_stereo import (
    assigned_cip_labels,
    validate_coordinate_stereochemistry,
)
from island.core import AtomSite, MolecularSystem
from island.evaluation import (
    EvaluationResult,
    OpenMMSinglePointEvaluator,
    PotentialEvaluator,
)
from island.exceptions import (
    EvaluationInputError,
    IslandError,
    MinimizationInputError,
    MinimizationUnavailableError,
    StereochemistryError,
)

from .models import (
    MinimizationOptions,
    MinimizationResult,
    MinimizationStep,
    coordinate_hash,
    force_metrics,
    system_identity,
)


def _scipy():
    try:
        import scipy
        from scipy.optimize import minimize
    except ImportError as error:
        raise MinimizationUnavailableError(
            "Local minimization requires pip install 'island[minimization]'"
        ) from error
    return scipy.__version__, minimize


class _Stop(Exception):
    def __init__(self, reason, message):
        super().__init__(message)
        self.reason = reason


def minimize_geometry(
    system: MolecularSystem,
    evaluator: PotentialEvaluator,
    options: MinimizationOptions | None = None,
) -> MinimizationResult:
    """Find a bounded local minimum; fmax is the maximum atomic force norm.

    Every real evaluator invocation is counted, including failed calls. One
    evaluation is reserved for final verification. Failed runs return the best
    admissible evaluated trial, which need not be a SciPy-accepted iterate.
    """
    options = options if options is not None else MinimizationOptions()
    if not isinstance(options, MinimizationOptions):
        raise MinimizationInputError("options must be MinimizationOptions")
    if not isinstance(system, MolecularSystem) or not isinstance(
        evaluator, PotentialEvaluator
    ):
        raise MinimizationInputError(
            "Supply a MolecularSystem and bound PotentialEvaluator"
        )
    owned = deepcopy(system)
    try:
        owned.validate()
        if owned.representation != "atomistic" or not owned.topology.sites:
            raise ValueError("Require a nonempty atomistic molecule")
        if owned.box is not None and any(owned.box.periodic):
            raise ValueError("Periodic minimization is unsupported")
        if any(type(site) is not AtomSite for site in owned.topology.sites.values()):
            raise ValueError("Only atomistic sites are supported")
        ids = tuple(sorted(owned.topology.sites))
        if any(type(site) is not int for site in ids):
            raise ValueError("Stable site IDs must be integers")
        xyz = np.array([owned.coordinates.get(site) for site in ids], dtype=float)
        if xyz.shape != (len(ids), 3) or not np.isfinite(xyz).all():
            raise ValueError("Initial coordinates must be finite N x 3 angstrom values")
        if isinstance(evaluator, OpenMMSinglePointEvaluator):
            evaluator.validate_system(owned)  # no energy call, graph cannot be bypassed
        identity = system_identity(owned)
    except (IslandError, TypeError, ValueError, AttributeError, KeyError) as error:
        raise MinimizationInputError(f"Invalid minimization input: {error}") from error
    initial_coordinates = {site: tuple(xyz[i]) for i, site in enumerate(ids)}
    try:
        expected_cip = assigned_cip_labels(owned)
        validate_coordinate_stereochemistry(
            owned,
            initial_coordinates,
            expected_cip,
            stage="initial minimization coordinates",
        )
    except ImportError as error:
        raise MinimizationUnavailableError(
            "RDKit is required to validate assigned tetrahedral stereochemistry"
        ) from error
    except IslandError as error:
        raise MinimizationInputError(str(error)) from error
    version, minimize = _scipy()
    evaluations = 0
    iterations = 0
    initial = None
    best = None
    latest = None
    history = []
    optimizer_success = None
    reason = "evaluation_failed"
    message = "Optimizer not started"
    verified = False

    def coordinates(x):
        array = np.asarray(x, dtype=float)
        if array.shape != (3 * len(ids),) or not np.isfinite(array).all():
            raise _Stop(
                "invalid_geometry",
                "Optimizer proposed nonfinite/malformed Cartesian coordinates",
            )
        return {site: tuple(array.reshape(-1, 3)[i]) for i, site in enumerate(ids)}

    def evaluate(x, *, final=False):
        nonlocal evaluations, best, latest
        coords = coordinates(x)
        limit = options.max_evaluations if final else options.max_evaluations - 1
        if evaluations >= limit:
            raise _Stop(
                "maximum_evaluations",
                "ISLAND evaluator budget reached (final check reserved)",
            )
        try:
            validate_coordinate_stereochemistry(
                owned, coords, expected_cip, stage="minimization trial"
            )
        except StereochemistryError as error:
            raise _Stop("stereochemistry_changed", str(error)) from error
        evaluations += 1
        try:
            result = evaluator.evaluate(dict(coords), coordinate_unit="angstrom")
        except EvaluationInputError as error:
            raise _Stop("invalid_geometry", str(error)) from error
        except Exception as error:
            raise _Stop("evaluation_failed", f"Evaluator failed: {error}") from error
        try:
            if not isinstance(result, EvaluationResult):
                raise ValueError("Evaluator must return EvaluationResult")
            if (
                result.energy_unit != "kJ/mol"
                or result.force_unit != "kJ/(mol*angstrom)"
                or result.coordinate_unit != "angstrom"
            ):
                raise ValueError("Evaluator changed or supplied unsupported units")
            if any(type(site) is not int for site in result.forces) or set(
                result.forces
            ) != set(ids):
                raise ValueError(
                    "Evaluator force-site coverage differs from stable IDs"
                )
            forces = np.array([result.forces[site] for site in ids], dtype=float)
            if (
                forces.shape != (len(ids), 3)
                or not np.isfinite(forces).all()
                or not np.isfinite(result.potential_energy)
                or not all(
                    np.isfinite(value) for value in result.energy_components.values()
                )
            ):
                raise ValueError(
                    "Evaluator supplied nonfinite/malformed energy or forces"
                )
            if result.coordinate_fingerprint != coordinate_hash(coords):
                raise ValueError(
                    "Evaluation fingerprint does not match the evaluated coordinates"
                )
            if not result.model_fingerprint or not result.parameter_fingerprint:
                raise ValueError("Evaluator must identify its model and parameters")
            if initial is not None and (
                result.model_fingerprint != initial.model_fingerprint
                or result.parameter_fingerprint != initial.parameter_fingerprint
            ):
                raise ValueError(
                    "Evaluator model/parameter fingerprint changed during minimization"
                )
            if not all(np.isfinite(value) for value in force_metrics(result)):
                raise ValueError("Force metrics are nonfinite")
        except (ValueError, TypeError, AttributeError) as error:
            raise _Stop("evaluation_failed", str(error)) from error
        latest = (np.array(x, copy=True), result)
        if best is None or result.potential_energy < best[1].potential_energy:
            best = latest
        return result

    def objective(x):
        if latest is not None and np.array_equal(x, latest[0]):
            result = latest[1]
        else:
            result = evaluate(x)
        gradient = -np.array(
            [result.forces[site] for site in ids], dtype=float
        ).reshape(-1)
        return result.potential_energy, gradient

    def accepted(x):
        nonlocal iterations
        iterations += 1
        if latest is None or not np.array_equal(x, latest[0]):
            raise _Stop(
                "evaluation_failed", "Accepted iterate lacks an exact evaluated state"
            )
        result = latest[1]
        fmax, rms = force_metrics(result)
        history.append(
            MinimizationStep(
                iterations,
                evaluations,
                result.potential_energy,
                fmax,
                rms,
                result.coordinate_fingerprint,
            )
        )
        if fmax <= options.force_tolerance:
            raise _Stop(
                "force_converged",
                "ISLAND force criterion reached at accepted iterate; SciPy interrupted by callback",
            )

    selected = None
    try:
        # Initial invocation may use the entire budget if the caller specified 1.
        initial = evaluate(xyz.reshape(-1), final=True)
        fmax, rms = force_metrics(initial)
        history.append(
            MinimizationStep(
                0,
                evaluations,
                initial.potential_energy,
                fmax,
                rms,
                initial.coordinate_fingerprint,
            )
        )
        if fmax <= options.force_tolerance:
            reason, message = (
                "force_converged",
                "Optimizer not started: initial geometry meets force tolerance",
            )
            selected = latest
        else:
            if evaluations >= options.max_evaluations - 1:
                raise _Stop(
                    "maximum_evaluations",
                    "Insufficient budget for an optimization trial and final check",
                )
            optimized = minimize(
                objective,
                xyz.reshape(-1),
                method="L-BFGS-B",
                jac=True,
                callback=accepted,
                options={
                    "maxiter": options.max_iterations,
                    "maxfun": options.max_evaluations,
                    "maxls": options.max_line_search_steps,
                    "gtol": options.force_tolerance / np.sqrt(3),
                    "ftol": options.energy_change_tolerance,
                },
            )
            message = str(optimized.message)
            optimizer_success = bool(optimized.success)
            if iterations >= options.max_iterations:
                reason = "maximum_iterations"
            elif "ABNORMAL" in message.upper() or "LINE SEARCH" in message.upper():
                reason = "line_search_failed"
            elif optimized.success:
                reason = "energy_stagnation"
            else:
                reason = "line_search_failed"
    except _Stop as stop:
        reason, message = stop.reason, str(stop)
        if reason == "force_converged":
            selected = latest
    except (RuntimeError, ValueError, FloatingPointError, TypeError) as error:
        reason, message = "evaluation_failed", f"SciPy failure: {error}"

    if selected is None:
        selected = best
    returned_state = (
        "accepted_iterate"
        if reason == "force_converged"
        else "best_admissible_evaluated_trial"
    )
    if selected is not None:
        try:
            final_result = evaluate(selected[0], final=True)
            selected = (selected[0], final_result)
            verified = True
        except _Stop as stop:
            reason, message = stop.reason, message + "; final check: " + str(stop)
            selected = best
            returned_state = "best_admissible_evaluated_trial"
    converged = bool(
        selected is not None
        and initial is not None
        and verified
        and reason
        not in ("evaluation_failed", "invalid_geometry", "stereochemistry_changed")
        and force_metrics(selected[1])[0] <= options.force_tolerance
        and selected[1].potential_energy
        <= initial.potential_energy + options.energy_increase_tolerance
    )
    if converged:
        reason = "force_converged"
    elif reason == "force_converged":
        reason = "energy_stagnation"
        message += "; independent final convergence checks not met"
    return MinimizationResult(
        initial_coordinates,
        initial_coordinates if selected is None else coordinates(selected[0]),
        initial,
        None if selected is None else selected[1],
        converged,
        reason,
        iterations,
        evaluations,
        message,
        optimizer_success,
        version,
        options,
        tuple(history),
        identity,
        "passed" if expected_cip else "not_assigned",
        "unevaluated_initial_coordinates" if selected is None else returned_state,
        verified,
    )
