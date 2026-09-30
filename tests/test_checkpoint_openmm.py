"""Actual fresh/session OpenMM continuation with independently verified boundaries."""

from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")
from test_dynamics import options as nve_options
from test_langevin import options as thermal_options
from test_singlepoint import source_case

from island.dynamics import (
    DynamicsSegmentOptions,
    create_dynamics_checkpoint,
    initialize_velocities,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    save_dynamics_checkpoint,
)
from island.evaluation import OpenMMSinglePointEvaluator


@pytest.mark.parametrize("thermal", [False, True])
def test_fresh_sessions_split_original_records(tmp_path, thermal):
    system, imported, _, _ = source_case(tmp_path)
    opts = (
        thermal_options(12, timestep_fs=0.1)
        if thermal
        else nve_options(12, dt=0.1, max_energy_deviation=0.1)
    )
    velocities = initialize_velocities(
        system, temperature_kelvin=300, seed=771
    ).velocities
    with OpenMMSinglePointEvaluator(system, imported).open_session() as session:
        whole = run_dynamics_segment(system, session, velocities, opts)
    with OpenMMSinglePointEvaluator(system, imported).open_session() as session:
        first = run_dynamics_segment(
            system, session, velocities, replace(opts, steps=4, max_evaluations=6)
        )
    for i in range(2):
        path = tmp_path / f"cp-{i}.json"
        save_dynamics_checkpoint(create_dynamics_checkpoint(first), path)
        cp = load_dynamics_checkpoint(path)
        moved = first.to_system(system)
        with OpenMMSinglePointEvaluator(moved, imported).open_session() as session:
            first = resume_dynamics(cp, moved, session, DynamicsSegmentOptions(4, 6, 5))
    a, b = whole.final_state, first.final_state
    # Predetermined same-platform bound: much tighter than independent prmtop conversion.
    for site in velocities:
        for left, right in (
            (a.coordinates[site], b.coordinates[site]),
            (a.velocities[site], b.velocities[site]),
            (a.evaluation.forces[site], b.evaluation.forces[site]),
        ):
            np.testing.assert_allclose(left, right, rtol=1e-12, atol=1e-10)
    assert a.total_energy == pytest.approx(b.total_energy, abs=1e-10, rel=1e-12)
    assert whole.payload["rng"] == first.payload["rng"]
    assert (
        first.payload["diagnostic"]["final_verified"]
        and first.payload["diagnostic"]["startup_verified"]
    )


def test_session_failure_keeps_fresh_final_verification(tmp_path, monkeypatch):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    velocities = initialize_velocities(
        system, temperature_kelvin=300, seed=771
    ).velocities
    original = evaluator._evaluate_context
    calls = []

    def fail(xyz, context, mm, unit):
        calls.append(True)
        if len(calls) == 2:
            context.setPositions(xyz * 0.1 * unit.nanometer)
            raise RuntimeError("injected backend failure after installation")
        return original(xyz, context, mm, unit)

    monkeypatch.setattr(evaluator, "_evaluate_context", fail)
    with evaluator.open_session() as session:
        result = run_dynamics_segment(
            system, session, velocities, thermal_options(2, timestep_fs=0.1)
        )
        assert session.closed
        assert result.evaluations == 3 and result.completed_steps == 0
        assert result.payload["diagnostic"]["final_verified"]
        assert result.payload["counters"]["normal_draws"] == 12
    from island.exceptions import InvalidDynamicsCheckpointError

    with pytest.raises(InvalidDynamicsCheckpointError):
        create_dynamics_checkpoint(result)
