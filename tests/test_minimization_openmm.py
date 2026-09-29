"""Real Reference-platform minimization and independent original-prmtop checks."""

import copy

import pytest

pytest.importorskip("scipy")
pytest.importorskip("openmm")
pytest.importorskip("parmed")
from test_singlepoint import assert_reference, source_case

from island.evaluation import OpenMMSinglePointEvaluator
from island.exceptions import MinimizationInputError
from island.minimization import MinimizationOptions, minimize_geometry


def test_real_openmm_minimum(tmp_path):
    system, imported, path, mapping = source_case(tmp_path)
    before = copy.deepcopy(system.to_dict())
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    result = minimize_geometry(
        system, evaluator, MinimizationOptions(force_tolerance=0.01)
    )
    assert result.converged, result.optimizer_message
    assert result.initial_fmax > 1
    assert result.final_fmax <= 0.01
    assert result.final_energy < result.initial_energy
    checked = assert_reference(evaluator, path, mapping, result.coordinates)
    assert checked.potential_energy == pytest.approx(result.final_energy, abs=1e-10)
    assert system.to_dict() == before
    assert (
        evaluator.evaluate().coordinate_fingerprint
        == result.initial_evaluation.coordinate_fingerprint
    )
    applied = result.to_system(system)
    assert evaluator.evaluate_system(applied).potential_energy == pytest.approx(
        result.final_energy
    )
    changed = copy.deepcopy(system)
    changed.topology.sites[10].metadata["cip_label"] = "R"
    with pytest.raises(MinimizationInputError):
        minimize_geometry(changed, evaluator)


def test_actual_singular_start_is_diagnostic_failure(tmp_path):
    from island import Coordinates

    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    coordinates = {site: system.coordinates.get(site) for site in system.topology.sites}
    coordinates[20] = coordinates[10]
    system.coordinates = Coordinates(coordinates)
    result = minimize_geometry(system, evaluator)
    assert not result.converged
    assert result.termination_reason == "invalid_geometry"
    assert result.evaluations == 1
    assert result.final_energy is None


def test_improper_order_preserved_during_minimization(tmp_path):
    system, imported, path, mapping = source_case(tmp_path, improper=True)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    result = minimize_geometry(
        system, evaluator, MinimizationOptions(force_tolerance=0.1)
    )
    assert result.converged, result.optimizer_message
    assert result.final_energy < result.initial_energy
    assert_reference(evaluator, path, mapping, result.coordinates)
