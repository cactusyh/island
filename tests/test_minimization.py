"""Independent optimizer acceptance with an analytical potential in angstroms."""

import copy
from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("scipy")

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.evaluation import EvaluationResult
from island.exceptions import (
    EvaluationError,
    EvaluationInputError,
    MinimizationInputError,
)
from island.minimization import MinimizationOptions, minimize_geometry
from island.minimization.models import coordinate_hash


class Harmonic:
    def __init__(self, stiffness=(1, 5, 11), fail_call=None, failure=None):
        self.stiffness = np.array(stiffness)
        self.calls = 0
        self.fail_call = fail_call
        self.failure = failure

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        self.calls += 1
        if self.calls == self.fail_call:
            raise (self.failure or EvaluationError("intentional evaluation failure"))
        assert coordinate_unit == "angstrom"
        energy = sum(
            0.5 * np.sum(self.stiffness * np.array(xyz) ** 2)
            for xyz in coordinates.values()
        )
        forces = {
            site: tuple(-self.stiffness * np.array(xyz))
            for site, xyz in coordinates.items()
        }
        return EvaluationResult(
            energy,
            {"harmonic": energy},
            forces,
            coordinate_hash(coordinates),
            "parameter-v1",
            "model-v1",
            "frame",
            "analytical",
            "1",
            "numpy",
            {},
        )


def molecule(reverse=False):
    graph = Topology()
    ids = (91, 17) if reverse else (17, 91)
    for site in ids:
        graph.add_site(
            AtomSite(
                site,
                "C",
                12.01,
                element="C",
                atomic_number=6,
                metadata={"repeat_index": site},
            )
        )
    return MolecularSystem(
        graph,
        Coordinates({17: (2.0, -1.0, 0.5), 91: (-1.0, 0.5, 1.0)}),
        metadata={
            "coordinate_source": "analytical_fixture",
            "ambertools_preparation": {"historical": "unchanged"},
        },
    )


def test_harmonic_sign_units_convergence_and_ownership():
    system = molecule()
    original = copy.deepcopy(system.to_dict())
    potential = Harmonic()
    options = MinimizationOptions(force_tolerance=1e-7)
    result = minimize_geometry(system, potential, options)
    assert result.converged and result.termination_reason == "force_converged"
    assert result.final_fmax <= options.force_tolerance
    assert result.final_energy < result.initial_energy
    np.testing.assert_allclose(list(result.coordinates.values()), 0, atol=1e-7)
    assert result.initial_fmax == pytest.approx(
        max(np.linalg.norm(f) for f in result.initial_evaluation.forces.values())
    )
    assert result.initial_rms_force == pytest.approx(
        np.sqrt(
            np.mean([np.dot(f, f) for f in result.initial_evaluation.forces.values()])
        )
    )
    assert potential.calls == result.evaluations
    assert result.evaluations <= options.max_evaluations
    assert result.history[-1].iteration == result.iterations
    independent = potential.evaluate(result.coordinates)
    assert independent.potential_energy == result.final_energy
    assert independent.forces == result.forces
    assert result.coordinate_unit == "angstrom"
    assert result.force_unit == "kJ/(mol*angstrom)"
    assert system.to_dict() == original
    with pytest.raises(TypeError):
        result.coordinates[17] = (0, 0, 0)
    assert copy.deepcopy(result) is result
    output = result.to_system(system)
    assert output is not system
    assert (
        output.metadata["ambertools_preparation"]
        == system.metadata["ambertools_preparation"]
    )
    assert output.metadata["coordinate_source"] == "local_minimization"
    assert output.topology.sites[17].metadata == system.topology.sites[17].metadata
    output.coordinates.set(17, (9, 9, 9))
    assert result.coordinates[17][0] != 9
    assert result.production_validated is False
    assert result.simulation_readiness == "not_established"


def test_already_converged_and_insertion_order():
    first = minimize_geometry(molecule(), Harmonic())
    second = minimize_geometry(molecule(True), Harmonic())
    assert dict(first.coordinates) == dict(second.coordinates)
    assert first.history == second.history
    system = molecule()
    for site in system.topology.sites:
        system.coordinates.set(site, (0, 0, 0))
    evaluator = Harmonic()
    result = minimize_geometry(system, evaluator)
    assert result.converged and result.iterations == 0
    assert result.evaluations == evaluator.calls == 2
    assert result.final_evaluation_verified


@pytest.mark.parametrize("budget", [1, 2, 3, 5])
def test_actual_evaluation_budget_exhaustion(budget):
    evaluator = Harmonic()
    result = minimize_geometry(
        molecule(),
        evaluator,
        MinimizationOptions(max_evaluations=budget, force_tolerance=1e-12),
    )
    assert result.evaluations == evaluator.calls == budget
    assert result.termination_reason == "maximum_evaluations"
    assert not result.converged
    assert result.final_energy <= result.initial_energy


def test_iteration_exhaustion_and_application_opt_in():
    system = molecule()
    result = minimize_geometry(
        system, Harmonic(), MinimizationOptions(max_iterations=1, force_tolerance=1e-12)
    )
    assert result.iterations == 1 and result.termination_reason == "maximum_iterations"
    assert not result.converged
    with pytest.raises(MinimizationInputError, match="allow_unconverged"):
        result.to_system(system)
    assert (
        result.to_system(system, allow_unconverged=True).metadata["coordinate_source"]
        == "local_minimization_diagnostic"
    )
    system.topology.sites[17].metadata["repeat_index"] = -1
    with pytest.raises(MinimizationInputError, match="provenance"):
        result.to_system(system, allow_unconverged=True)


def test_energy_stagnation_is_not_force_convergence():
    result = minimize_geometry(
        molecule(),
        Harmonic(),
        MinimizationOptions(force_tolerance=1e-14, energy_change_tolerance=0.9),
    )
    assert result.optimizer_success
    assert result.termination_reason == "energy_stagnation"
    assert not result.converged and result.final_fmax > result.options.force_tolerance


@pytest.mark.parametrize(
    "failure,reason",
    [
        (EvaluationError("bad backend"), "evaluation_failed"),
        (EvaluationInputError("singular trial"), "invalid_geometry"),
    ],
)
def test_failed_evaluation_preserves_best_actual_state(failure, reason):
    evaluator = Harmonic(fail_call=3, failure=failure)
    result = minimize_geometry(molecule(), evaluator)
    assert not result.converged and result.termination_reason == reason
    assert (
        result.evaluations == evaluator.calls == 4
    )  # includes failure and final recovery check
    assert result.returned_state == "best_admissible_evaluated_trial"
    assert result.final_energy <= result.initial_energy
    assert (
        Harmonic().evaluate(result.coordinates).potential_energy == result.final_energy
    )


def test_initial_failure_is_not_fabricated():
    result = minimize_geometry(molecule(), Harmonic(fail_call=1))
    assert result.evaluations == 1 and not result.converged
    assert result.initial_energy is None and result.final_energy is None
    assert result.final_fmax is None and not result.forces
    with pytest.raises(MinimizationInputError, match="No valid"):
        result.to_system(molecule(), allow_unconverged=True)


@pytest.mark.parametrize(
    "mode", ["units", "coverage", "model", "coordinate", "nonfinite"]
)
def test_inconsistent_evaluator_rejected(mode):
    class Bad(Harmonic):
        def evaluate(self, coordinates=None, **kwargs):
            result = super().evaluate(coordinates, **kwargs)
            if self.calls != 2:
                return result
            if mode == "units":
                object.__setattr__(result, "force_unit", "kJ/(mol*nm)")
            elif mode == "coverage":
                result = replace(result, forces={17: (0, 0, 0)})
            elif mode == "model":
                result = replace(result, model_fingerprint="changed")
            elif mode == "coordinate":
                result = replace(result, coordinate_fingerprint="wrong")
            else:
                object.__setattr__(result, "potential_energy", float("nan"))
            return result

    result = minimize_geometry(molecule(), Bad())
    assert not result.converged and result.termination_reason == "evaluation_failed"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_iterations": 0},
        {"max_evaluations": 0},
        {"max_line_search_steps": -1},
        {"force_tolerance": float("nan")},
        {"max_evaluations": True},
    ],
)
def test_invalid_options(kwargs):
    with pytest.raises(MinimizationInputError):
        MinimizationOptions(**kwargs)


def test_bounded_line_search_failure():
    class WrongGradient(Harmonic):
        def evaluate(self, *args, **kwargs):
            result = super().evaluate(*args, **kwargs)
            return replace(
                result,
                forces={
                    site: tuple(-np.array(force))
                    for site, force in result.forces.items()
                },
            )

    evaluator = WrongGradient()
    result = minimize_geometry(
        molecule(), evaluator, MinimizationOptions(max_line_search_steps=2)
    )
    assert not result.converged
    assert result.termination_reason == "line_search_failed"
    assert result.evaluations == evaluator.calls == 4
    assert result.final_energy == result.initial_energy


def test_changed_final_energy_cannot_be_force_converged():
    class ChangedFinal(Harmonic):
        def evaluate(self, *args, **kwargs):
            result = super().evaluate(*args, **kwargs)
            if self.calls == 2:
                result = replace(
                    result, potential_energy=10, energy_components={"harmonic": 10}
                )
            return result

    system = molecule()
    system.coordinates = Coordinates(
        {site: (0.0, 0.0, 0.0) for site in system.topology.sites}
    )
    result = minimize_geometry(system, ChangedFinal())
    assert not result.converged
    assert result.final_fmax == 0
    assert result.final_energy > result.initial_energy
