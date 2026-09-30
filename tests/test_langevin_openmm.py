"""Actual OpenMM force/model comparison and reusable-session integration."""

import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")
from test_langevin import options
from test_singlepoint import source_case

from island.dynamics import initialize_velocities, run_langevin
from island.evaluation import OpenMMSinglePointEvaluator
from island.exceptions import DynamicsInputError

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dynamics_reference import compare
from langevin_reference import reference_trajectory


def test_original_prmtop_baoab_with_prescribed_noise(tmp_path, monkeypatch):
    system, imported, prmtop, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    opts = options(40, timestep_fs=0.1)
    initial = initialize_velocities(system, temperature_kelvin=300, seed=771)
    noise = np.random.Generator(np.random.PCG64(opts.thermostat_seed)).standard_normal(
        (40, 4, 3)
    )
    increments = iter(noise)
    monkeypatch.setattr(
        "island.dynamics.langevin._normal_increments",
        lambda rng, shape: next(increments).copy(),
    )
    with evaluator.open_session() as session:
        result = run_langevin(system, session, initial.velocities, opts)
    refs, _ = reference_trajectory(
        prmtop,
        dict(imported.mapping),
        result.initial_state.coordinates,
        initial.velocities,
        result.masses,
        opts,
        [f.step for f in result.frames],
        noise,
    )
    compare(result, refs)
    assert result.completed and result.evaluations == 42 and result.normal_draws == 480
    result.to_system(system).validate()
    # Ordinary RNG path produces the same prescribed increments.
    from island.dynamics.thermal import _normal_increments

    monkeypatch.setattr(
        "island.dynamics.langevin._normal_increments", _normal_increments
    )
    with evaluator.open_session() as session:
        natural = run_langevin(system, session, initial.velocities, opts)
    assert result.frames == natural.frames


def test_binding_session_failure_and_fresh_verification(tmp_path, monkeypatch):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    initial = initialize_velocities(system, temperature_kelvin=300, seed=771)
    with evaluator.open_session() as session:
        for change in ("mass", "graph"):
            wrong = deepcopy(system)
            if change == "mass":
                wrong.topology.sites[10].mass += 1
            else:
                wrong.topology.sites[10].element = "N"
            with pytest.raises(DynamicsInputError):
                run_langevin(wrong, session, initial.velocities, options(2))
        original = evaluator._evaluate_context
        calls = []

        def failing(xyz, context, mm, unit):
            calls.append(True)
            if len(calls) == 2:
                context.setPositions(xyz * 0.1 * unit.nanometer)
                raise RuntimeError("backend failed after installation")
            return original(xyz, context, mm, unit)

        monkeypatch.setattr(evaluator, "_evaluate_context", failing)
        result = run_langevin(system, session, initial.velocities, options(2))
        assert session.closed
        assert result.final_evaluation_verified and result.evaluations == 3
        assert result.completed_steps == 0 and result.normal_draws == 12
        assert result.final_state.velocities == initial.velocities
        result.to_system(system, allow_incomplete=True).validate()
