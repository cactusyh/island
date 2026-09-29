"""Real NVE vs independent original-prmtop synchronized CustomIntegrator."""

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")
from test_singlepoint import source_case

from island.dynamics import DynamicsOptions, run_nve
from island.evaluation import OpenMMSinglePointEvaluator
from island.exceptions import DynamicsInputError


@pytest.mark.parametrize("improper", [False, True])
def test_real_openmm_independent_trajectory(tmp_path, monkeypatch, improper):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    from dynamics_reference import compare, reference_trajectory

    system, imported, path, mapping = source_case(tmp_path, improper=improper)
    before = deepcopy(system.to_dict())
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    velocities = {
        site: tuple(0.2 * np.sin(site + np.arange(3))) for site in system.topology.sites
    }
    result = run_nve(
        system,
        evaluator,
        velocities,
        DynamicsOptions(
            0.1, 20, 22, 5, recording_interval=5, max_energy_deviation=0.05
        ),
    )
    assert result.completed, result.message
    refs, _ = reference_trajectory(
        path,
        mapping,
        result.initial_state.coordinates,
        velocities,
        result.masses,
        0.1,
        [f.step for f in result.frames],
    )
    compare(result, refs)
    output = result.to_system(system)
    output.validate()
    assert evaluator.evaluate_system(output).potential_energy == pytest.approx(
        result.final_state.potential_energy, abs=1e-10
    )
    assert (
        evaluator.evaluate().coordinate_fingerprint
        == result.initial_state.coordinate_fingerprint
    )
    assert system.to_dict() == before
    assert result.evaluations == 22
    changed = deepcopy(system)
    changed.topology.sites[10].mass += 1
    with pytest.raises(DynamicsInputError):
        run_nve(changed, evaluator, velocities, result.options)


def test_real_singular_trial_returns_last_synchronized_state(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    starting = evaluator.evaluate()
    dt = 0.0001
    velocities = {
        site: tuple(
            -0.5
            * dt
            * 100
            * np.array(starting.forces[site])
            / system.topology.sites[site].mass
        )
        for site in system.topology.sites
    }
    velocities[20] = tuple(
        np.array(velocities[20])
        + (system.coordinates.get(10) - system.coordinates.get(20)) / dt
    )
    result = run_nve(system, evaluator, velocities, DynamicsOptions(0.1, 2, 4, 3))
    assert result.termination_reason == "invalid_geometry"
    assert result.completed_steps == 0 and result.evaluations == 3
    assert result.failed_trial_evaluations == 1 and result.final_evaluation_verified
    assert result.final_state.coordinates == result.initial_state.coordinates
    assert result.final_state.velocities == result.initial_state.velocities
