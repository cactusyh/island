"""Reconstructed public results must be validated before application."""

from copy import deepcopy
from dataclasses import replace

import pytest

pytest.importorskip("scipy")
from test_minimization import Harmonic, molecule

from island.exceptions import MinimizationInputError
from island.minimization import minimize_geometry
from island.minimization.models import coordinate_hash


def test_missing_coordinate_reconstruction_rejected():
    system = molecule()
    before = deepcopy(system.to_dict())
    result = minimize_geometry(system, Harmonic())
    reduced = {17: result.coordinates[17]}
    reconstructed = replace(
        result,
        coordinates=reduced,
        final_evaluation=replace(
            result.final_evaluation,
            forces={17: result.forces[17]},
            coordinate_fingerprint=coordinate_hash(reduced),
        ),
    )
    with pytest.raises(MinimizationInputError):
        reconstructed.to_system(system)
    assert system.to_dict() == before


@pytest.fixture
def valid():
    system = molecule()
    return system, minimize_geometry(system, Harmonic())


@pytest.mark.parametrize(
    "changes",
    [
        {"coordinates": {17: (0, 0)}},
        {"coordinates": {17: (0, 0, 0, 0), 91: (0, 0, 0)}},
        {"coordinates": {17: (float("nan"), 0, 0), 91: (0, 0, 0)}},
        {"coordinates": {17: (float("inf"), 0, 0), 91: (0, 0, 0)}},
        {"coordinates": {17: ("0", 0, 0), 91: (0, 0, 0)}},
        {"coordinates": {"17": (0, 0, 0), 91: (0, 0, 0)}},
        {"coordinates": {True: (0, 0, 0), 91: (0, 0, 0)}},
        {"coordinates": {17: (10**400, 0, 0), 91: (0, 0, 0)}},
        {"coordinates": None},
        {"initial_coordinates": {}},
        {"initial_evaluation": None},
        {"final_evaluation": {}},
        {"options": {}},
        {"evaluations": 0},
        {"evaluations": True},
        {"iterations": -1},
        {"iterations": 10000},
        {"history": [{}]},
        {"converged": False},
        {"converged": 1},
        {"final_evaluation_verified": False},
        {"final_evaluation_verified": "yes"},
        {"termination_reason": "success"},
        {"optimizer_success": 1},
        {"optimizer_message": []},
        {"optimizer_version": []},
        {"returned_state": "unevaluated_initial_coordinates"},
        {"stereochemistry": "trusted"},
        {"system_fingerprint": []},
    ],
)
def test_malformed_public_record_uses_domain_exception(valid, changes):
    from island.exceptions import InvalidMinimizationResultError

    system, result = valid
    before = deepcopy(system.to_dict())
    with pytest.raises(InvalidMinimizationResultError):
        bad = replace(result, **changes)
        bad.validate_integrity()
    with pytest.raises(InvalidMinimizationResultError):
        bad = replace(result, **changes)
        bad.to_system(system, allow_unconverged=True)
    assert system.to_dict() == before


@pytest.mark.parametrize(
    "field,value",
    [
        ("coordinate_fingerprint", "wrong"),
        ("model_fingerprint", "different"),
        ("parameter_fingerprint", "different"),
        ("forces", {17: (0.0, 0.0, 0.0)}),
        ("energy_components", {"other": 10000.0}),
        ("platform", []),
    ],
)
def test_inconsistent_evaluation(valid, field, value):
    from island.exceptions import InvalidMinimizationResultError

    system, result = valid
    bad = replace(
        result, final_evaluation=replace(result.final_evaluation, **{field: value})
    )
    with pytest.raises(InvalidMinimizationResultError):
        bad.validate_integrity()
    with pytest.raises(InvalidMinimizationResultError):
        bad.to_system(system, allow_unconverged=True)


@pytest.mark.parametrize(
    "field,value",
    [
        ("potential_energy", float("nan")),
        ("force_unit", "newtons"),
        ("settings", {"mutable": []}),
    ],
)
def test_malformed_evaluation_fields(valid, field, value):
    from island.exceptions import InvalidMinimizationResultError

    _, result = valid
    evaluation = replace(result.final_evaluation)
    object.__setattr__(evaluation, field, value)
    bad = replace(result, final_evaluation=evaluation)
    with pytest.raises(InvalidMinimizationResultError):
        bad.validate_integrity()
    with pytest.raises(InvalidMinimizationResultError):
        deepcopy(bad)


def test_unknown_ids_even_with_consistent_reconstructed_frames(valid):
    from island.exceptions import InvalidMinimizationResultError

    system, result = valid

    def remap(mapping):
        return {site + 1000: value for site, value in mapping.items()}

    initial, final = remap(result.initial_coordinates), remap(result.coordinates)
    history = list(result.history)
    history[0] = replace(history[0], coordinate_fingerprint=coordinate_hash(initial))
    history[-1] = replace(history[-1], coordinate_fingerprint=coordinate_hash(final))
    bad = replace(
        result,
        initial_coordinates=initial,
        coordinates=final,
        history=history,
        initial_evaluation=replace(
            result.initial_evaluation,
            forces=remap(result.initial_evaluation.forces),
            coordinate_fingerprint=coordinate_hash(initial),
        ),
        final_evaluation=replace(
            result.final_evaluation,
            forces=remap(result.forces),
            coordinate_fingerprint=coordinate_hash(final),
        ),
    )
    bad.validate_integrity()  # internally coherent; incompatible with this target
    with pytest.raises(InvalidMinimizationResultError, match="coverage"):
        bad.to_system(system, allow_unconverged=True)


def test_extra_site_and_invalid_history(valid):
    from island.exceptions import InvalidMinimizationResultError

    system, result = valid
    coords = {**result.coordinates, 999: (0.0, 0.0, 0.0)}
    evaluation = replace(
        result.final_evaluation,
        forces={**result.forces, 999: (0.0, 0.0, 0.0)},
        coordinate_fingerprint=coordinate_hash(coords),
    )
    with pytest.raises(InvalidMinimizationResultError):
        replace(result, coordinates=coords, final_evaluation=evaluation).to_system(
            system, allow_unconverged=True
        )
    for step in (
        replace(result.history[0], energy=999),
        replace(result.history[0], evaluations=2),
        replace(result.history[0], fmax=float("nan")),
        replace(result.history[0], coordinate_fingerprint="wrong"),
    ):
        with pytest.raises(InvalidMinimizationResultError):
            replace(result, history=(step, *result.history[1:])).validate_integrity()


def test_valid_incomplete_diagnostics_and_owned_contents():
    from island.minimization import MinimizationOptions

    system = molecule()
    failure = minimize_geometry(system, Harmonic(fail_call=1))
    failure.validate_integrity()
    assert deepcopy(failure) is failure
    with pytest.raises(MinimizationInputError):
        failure.to_system(system, allow_unconverged=True)
    for budget in (1, 2, 3):
        diagnostic = minimize_geometry(
            system, Harmonic(), MinimizationOptions(max_evaluations=budget)
        )
        diagnostic.validate_integrity()
        assert not diagnostic.converged
        assert diagnostic.final_evaluation_verified == (budget > 1)
        with pytest.raises(MinimizationInputError, match="allow_unconverged"):
            diagnostic.to_system(system)
        diagnostic.to_system(system, allow_unconverged=True).validate()
    result = minimize_geometry(system, Harmonic())
    coordinates = {site: list(xyz) for site, xyz in result.coordinates.items()}
    history = list(result.history)
    copied = replace(result, coordinates=coordinates, history=history)
    coordinates[17][0] = 999
    history.clear()
    copied.validate_integrity()
    assert deepcopy(copied) is copied
    assert copied.coordinates == result.coordinates
    with pytest.raises(TypeError):
        copied.final_evaluation.settings["new"] = 1


@pytest.mark.parametrize(
    "mode", ["force", "energy", "options", "verified_calls", "termination_count"]
)
def test_convergence_and_verification_consistency(valid, mode):
    from island.exceptions import InvalidMinimizationResultError

    _, result = valid
    if mode == "force":
        final = replace(
            result.final_evaluation,
            forces={site: (999.0, 0.0, 0.0) for site in result.coordinates},
        )
        bad = replace(result, final_evaluation=final)
    elif mode == "energy":
        final = replace(
            result.final_evaluation,
            potential_energy=result.initial_energy + 100,
            energy_components={"harmonic": result.initial_energy + 100},
        )
        bad = replace(result, final_evaluation=final)
    elif mode == "options":
        options = replace(result.options)
        object.__setattr__(options, "force_tolerance", float("nan"))
        bad = replace(result, options=options)
    elif mode == "verified_calls":
        bad = replace(result, evaluations=result.history[-1].evaluations)
    else:
        bad = replace(
            result,
            converged=False,
            termination_reason="maximum_evaluations",
            returned_state="best_admissible_evaluated_trial",
        )
    with pytest.raises(InvalidMinimizationResultError):
        bad.validate_integrity()


def test_numeric_array_reconstruction_and_malformed_properties(valid):
    import numpy as np

    from island.exceptions import InvalidMinimizationResultError

    _, result = valid
    reconstructed = replace(
        result,
        coordinates={
            site: np.array(xyz, dtype=np.float64)
            for site, xyz in result.coordinates.items()
        },
    )
    reconstructed.validate_integrity()
    assert deepcopy(reconstructed) is reconstructed
    malformed = replace(result, final_evaluation={})
    for name in ("final_energy", "final_fmax", "forces", "coordinate_fingerprint"):
        with pytest.raises(InvalidMinimizationResultError):
            getattr(malformed, name)
