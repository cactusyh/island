"""Real OpenMM session ownership, validation, and run integration."""

from copy import deepcopy

import numpy as np
import pytest

mm = pytest.importorskip("openmm")
pytest.importorskip("parmed")
from test_singlepoint import XYZ, source_case

from island.dynamics import DynamicsOptions, run_nve
from island.evaluation import OpenMMSinglePointEvaluator
from island.exceptions import (
    DynamicsInputError,
    EvaluationError,
    EvaluationInputError,
    EvaluationSessionBusyError,
    EvaluationSessionClosedError,
    MinimizationInputError,
)
from island.minimization import MinimizationOptions, minimize_geometry


def frame(xyz):
    return dict(zip((10, 20, 30, 40), xyz, strict=True))


def same(left, right):
    assert left.potential_energy == pytest.approx(
        right.potential_energy, rel=1e-12, abs=1e-10
    )
    assert left.energy_components == pytest.approx(
        right.energy_components, rel=1e-12, abs=1e-10
    )
    for site in left.forces:
        np.testing.assert_allclose(
            left.forces[site], right.forces[site], rtol=1e-12, atol=1e-10
        )
    assert left.coordinate_fingerprint == right.coordinate_fingerprint
    assert left.evaluation_fingerprint == right.evaluation_fingerprint
    assert left.model_fingerprint == right.model_fingerprint
    assert left.parameter_fingerprint == right.parameter_fingerprint


def test_real_frames_sessions_ownership_and_cleanup(tmp_path, monkeypatch):
    system, imported, _, _ = source_case(tmp_path)
    before = deepcopy(system.to_dict())
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    original_context = mm.Context
    counts = {"contexts": 0}

    def counted_context(*args, **kwargs):
        counts["contexts"] += 1
        return original_context(*args, **kwargs)

    monkeypatch.setattr(mm, "Context", counted_context)
    a = frame(XYZ)
    b = frame(XYZ + np.arange(12).reshape(4, 3) * 0.002)
    with evaluator.open_session() as one, evaluator.open_session() as two:
        assert counts["contexts"] == 2
        first = one.evaluate(dict(reversed(list(a.items()))))
        one_b = one.evaluate(b)
        other_b = two.evaluate(b)
        again = one.evaluate(a)
        same(one_b, other_b)
        same(first, again)
        same(first, evaluator.evaluate(a))
        same(one_b, evaluator.evaluate(b))
        assert first.coordinate_fingerprint != one_b.coordinate_fingerprint
        assert first.potential_energy == again.potential_energy
        assert counts["contexts"] == 4
    one.close()
    two.close()
    with pytest.raises(EvaluationSessionClosedError):
        one.evaluate(a)
    assert counts["contexts"] == 4
    assert system.to_dict() == before
    same(evaluator.evaluate(), first)


@pytest.mark.parametrize(
    "bad,unit",
    [
        ({10: (0, 0, 0)}, "angstrom"),
        (frame(XYZ), "nanometer"),
        ({**frame(XYZ), 10: (float("nan"), 0, 0)}, "angstrom"),
        ({**frame(XYZ), 10: (float("inf"), 0, 0)}, "angstrom"),
        ({**frame(XYZ), 20: XYZ[0]}, "angstrom"),
        ({**frame(XYZ), "10": XYZ[0]}, "angstrom"),
    ],
)
def test_invalid_frame_does_not_poison_context(tmp_path, bad, unit):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    with evaluator.open_session() as session:
        before = session.evaluate(frame(XYZ))
        with pytest.raises(EvaluationInputError):
            session.evaluate(bad, coordinate_unit=unit)
        same(before, session.evaluate(frame(XYZ)))


def test_exception_busy_and_backend_failure_close(tmp_path, monkeypatch):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    with pytest.raises(RuntimeError, match="body"), evaluator.open_session() as session:
        assert session._lock.acquire(blocking=False)
        try:
            with pytest.raises(EvaluationSessionBusyError):
                session.evaluate()
            with pytest.raises(EvaluationSessionBusyError):
                session.close()
        finally:
            session._lock.release()
        raise RuntimeError("body")
    with pytest.raises(EvaluationSessionClosedError):
        session.evaluate()
    session.close()
    session = evaluator.open_session()

    def broken(*_args):
        raise RuntimeError("backend failed after setPositions")

    monkeypatch.setattr(evaluator, "_evaluate_context", broken)
    with pytest.raises(EvaluationError, match="backend failed"):
        session.evaluate()
    with pytest.raises(EvaluationSessionClosedError):
        session.evaluate()
    session.close()


def test_bound_system_validation_nve_and_minimization(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    velocities = {site: (0.1, 0.2, -0.1) for site in system.topology.sites}
    options = DynamicsOptions(0.1, 2, 4, 3)
    with evaluator.open_session() as session:
        for change in ("mass", "graph"):
            altered = deepcopy(system)
            if change == "mass":
                altered.topology.sites[10].mass += 1
            else:
                altered.topology.bonds.pop(next(iter(altered.topology.bonds)))
            with pytest.raises(DynamicsInputError):
                run_nve(altered, session, velocities, options)
            with pytest.raises(MinimizationInputError):
                minimize_geometry(altered, session)


def test_nve_and_minimization_fresh_final_and_context_counts(tmp_path, monkeypatch):
    pytest.importorskip("scipy")
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    velocities = {
        site: tuple(0.2 * np.sin(site + np.arange(3))) for site in system.topology.sites
    }
    original_context = mm.Context
    contexts = []

    def counted_context(*args, **kwargs):
        result = original_context(*args, **kwargs)
        contexts.append(1)
        return result

    monkeypatch.setattr(mm, "Context", counted_context)
    options = DynamicsOptions(
        0.1, 5, 7, 4, recording_interval=2, max_energy_deviation=0.05
    )
    fresh = run_nve(system, evaluator, velocities, options)
    assert fresh.completed and fresh.evaluations == 7 and len(contexts) == 7
    with evaluator.open_session() as session:
        before = len(contexts)
        reused = run_nve(system, session, velocities, options)
        assert reused.completed and reused.evaluations == 7
        assert len(contexts) - before == 1
        assert reused.final_evaluation is not None
        assert (
            reused.final_evaluation.evaluation_fingerprint
            == fresh.final_evaluation.evaluation_fingerprint
        )
    assert len(contexts) == 9  # session plus independent fresh final
    for left, right in zip(fresh.frames, reused.frames, strict=True):
        assert left.step == right.step
        assert left.coordinates == right.coordinates
        assert left.velocities == right.velocities
        same(left.evaluation, right.evaluation)
        assert left.total_energy == right.total_energy
    min_options = MinimizationOptions(force_tolerance=0.1)
    fresh_min = minimize_geometry(system, evaluator, min_options)
    with evaluator.open_session() as session:
        before = len(contexts)
        reused_min = minimize_geometry(system, session, min_options)
        assert len(contexts) - before == 1
    assert len(contexts) - before == 1
    assert fresh_min.converged == reused_min.converged
    assert fresh_min.final_energy == pytest.approx(reused_min.final_energy, abs=1e-8)
    assert fresh_min.final_fmax == pytest.approx(reused_min.final_fmax, abs=1e-8)
    assert fresh_min.evaluations == reused_min.evaluations


def test_failed_trial_uses_fresh_last_state_verification(tmp_path, monkeypatch):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    velocities = {site: (0.1, 0, 0) for site in system.topology.sites}
    with evaluator.open_session() as session:
        original = session.evaluate
        calls = 0

        def fail_second(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                session.close()
                raise EvaluationError("injected backend failure")
            return original(*args, **kwargs)

        monkeypatch.setattr(session, "evaluate", fail_second)
        result = run_nve(system, session, velocities, DynamicsOptions(0.1, 2, 4, 3))
        assert result.termination_reason == "evaluation_failed"
        assert result.failed_trial_evaluations == 1
        assert result.evaluations == 3 and result.final_evaluation_verified
        assert result.final_state.coordinates == result.initial_state.coordinates
        assert result.final_state.velocities == result.initial_state.velocities


def test_owner_thread_nested_entry_no_steps_and_resource_release(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from copy import copy

    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)

    def forbidden(*args):
        raise AssertionError("Evaluator must never step its integrator")

    monkeypatch.setattr(mm.VerletIntegrator, "step", forbidden)
    with evaluator.open_session() as session:
        same(session.evaluate(), evaluator.evaluate())
        with pytest.raises(EvaluationSessionBusyError), session:
            pass
        with ThreadPoolExecutor(max_workers=1) as pool:
            for method in (session.evaluate, session.close, session.evaluate_fresh):
                with pytest.raises(EvaluationSessionBusyError):
                    pool.submit(method).result()
        with pytest.raises(EvaluationInputError):
            copy(session)
        with pytest.raises(EvaluationInputError):
            deepcopy(session)
        assert not session.closed
    assert session.closed
    assert session._resources is None
    session.close()
    with pytest.raises(EvaluationSessionClosedError):
        session.__enter__()


@pytest.mark.parametrize("mode", ["nve", "minimize"])
def test_actual_backend_failure_invalidates_session_and_fresh_verifies(
    tmp_path, monkeypatch, mode
):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    contexts = []
    create = mm.Context

    def counted(*args, **kwargs):
        contexts.append(1)
        return create(*args, **kwargs)

    monkeypatch.setattr(mm, "Context", counted)
    with evaluator.open_session() as session:
        evaluate = evaluator._evaluate_context
        calls = 0

        def fail_after_positions(xyz, context, backend, unit):
            nonlocal calls
            calls += 1
            if calls == 2:
                context.setPositions(xyz * 0.1 * unit.nanometer)
                raise RuntimeError("injected uncertain backend state")
            return evaluate(xyz, context, backend, unit)

        monkeypatch.setattr(evaluator, "_evaluate_context", fail_after_positions)
        if mode == "nve":
            result = run_nve(
                system,
                session,
                {site: (0.1, 0, 0) for site in system.topology.sites},
                DynamicsOptions(0.1, 5, 7, 6),
            )
            assert result.completed_steps == 0
            result.to_system(system, allow_incomplete=True).validate()
        else:
            pytest.importorskip("scipy")
            result = minimize_geometry(system, session)
            result.to_system(system, allow_unconverged=True).validate()
        assert session.closed and session._resources is None
        assert result.termination_reason == "evaluation_failed"
        assert result.evaluations == 3 and result.final_evaluation_verified
        assert len(contexts) == 2
        result.validate_integrity()
        same(evaluator.evaluate(), session.evaluate_fresh())


@pytest.mark.parametrize("budget", [1, 2, 3])
@pytest.mark.parametrize("mode", ["nve", "minimize"])
def test_session_budget_accounting(tmp_path, monkeypatch, budget, mode):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    count = []
    create = mm.Context

    def counted(*args, **kwargs):
        count.append(1)
        return create(*args, **kwargs)

    monkeypatch.setattr(mm, "Context", counted)
    with evaluator.open_session() as session:
        if mode == "nve":
            result = run_nve(
                system,
                session,
                {site: (0.1, 0, 0) for site in system.topology.sites},
                DynamicsOptions(0.1, 10, budget, 11),
            )
        else:
            pytest.importorskip("scipy")
            result = minimize_geometry(
                system,
                session,
                MinimizationOptions(max_evaluations=budget, force_tolerance=1e-12),
            )
        assert result.evaluations == budget
        assert result.termination_reason == "maximum_evaluations"
        assert result.final_evaluation_verified == (budget > 1)
        assert len(count) == (2 if budget > 1 else 1)
        result.validate_integrity()


def test_binding_stereo_and_parameter_integrity(tmp_path):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    with evaluator.open_session() as session:
        altered = deepcopy(system)
        altered.topology.sites[10].metadata["cip_label"] = "R"
        with pytest.raises(DynamicsInputError):
            run_nve(
                altered,
                session,
                {site: (0.1, 0, 0) for site in system.topology.sites},
                DynamicsOptions(0.1, 2, 4, 3),
            )
        with pytest.raises(MinimizationInputError):
            minimize_geometry(altered, session)
        from dataclasses import replace

        evaluator._imported = replace(evaluator._imported, result_signature="corrupt")
        with pytest.raises(DynamicsInputError):
            run_nve(
                system,
                session,
                {site: (0.1, 0, 0) for site in system.topology.sites},
                DynamicsOptions(0.1, 2, 4, 3),
            )


def test_reentrant_backend_callback_and_nonfinite_failure(tmp_path, monkeypatch):
    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    method = evaluator._evaluate_context
    reference = evaluator.evaluate()
    with evaluator.open_session() as session:

        def nested(xyz, context, backend, unit):
            with pytest.raises(EvaluationSessionBusyError):
                session.evaluate()
            return method(xyz, context, backend, unit)

        monkeypatch.setattr(evaluator, "_evaluate_context", nested)
        same(session.evaluate(), reference)

        def invalid(*args):
            raise EvaluationError("nonfinite backend forces")

        monkeypatch.setattr(evaluator, "_evaluate_context", invalid)
        with pytest.raises(EvaluationError, match="nonfinite"):
            session.evaluate()
        assert session.closed


def test_failed_context_construction_releases_partial_resources(tmp_path, monkeypatch):
    from island.evaluation.openmm import _OpenMMResources
    from island.exceptions import EvaluationUnavailableError

    system, imported, _, _ = source_case(tmp_path)
    evaluator = OpenMMSinglePointEvaluator(system, imported)
    original = _OpenMMResources.close
    cleaned = []

    def close(bundle):
        original(bundle)
        cleaned.append(bundle.system is bundle.integrator is bundle.context is None)

    def unavailable(*args, **kwargs):
        raise RuntimeError("injected Context construction failure")

    monkeypatch.setattr(_OpenMMResources, "close", close)
    monkeypatch.setattr(mm, "Context", unavailable)
    with pytest.raises(EvaluationUnavailableError, match="construction failure"):
        evaluator.open_session()
    assert cleaned == [True]
