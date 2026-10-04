"""Real OpenMM lifecycle; synthetic PCFF software fixture, not source acceptance."""
# ruff: noqa: F811 -- imported pytest fixture injection

from concurrent.futures import ThreadPoolExecutor
from copy import copy, deepcopy

import pytest
from test_pcff_evaluation import bound, parameters, synthetic  # noqa: F401

from island.exceptions import (
    EvaluationError,
    EvaluationInputError,
    EvaluationSessionBusyError,
    EvaluationSessionClosedError,
    MinimizationInputError,
)
from island.minimization import MinimizationOptions, minimize_geometry


def test_frames_and_lifecycle(bound, monkeypatch):
    import openmm as mm

    system, _, e = bound
    before = deepcopy(system.to_dict())
    a = {s: system.coordinates.get(s) for s in system.topology.sites}
    b = {
        s: v + [0.015 if s == 11 else 0, 0.02, 0.01]
        for s, v in reversed(list(a.items()))
    }
    first = e.evaluate(a)
    monkeypatch.setattr(mm.VerletIntegrator, "step", lambda *a: pytest.fail("Stepped"))
    with e.open_session() as one, e.open_session() as two:
        assert one._resources.context is not two._resources.context
        assert one.evaluate(a) == two.evaluate(a) == first
        prior = one.evaluate(b)
        assert prior.potential_energy != first.potential_energy
        assert two.evaluate(a) == one.evaluate(a) == first
        assert prior == e.evaluate(b) == one.evaluate_fresh(b)
        for copier in (copy, deepcopy):
            with pytest.raises(EvaluationInputError):
                copier(one)
        with pytest.raises(EvaluationSessionBusyError), one:
            pass
        one._lock.acquire()
        try:
            with pytest.raises(EvaluationSessionBusyError):
                one.evaluate()
        finally:
            one._lock.release()
        with (
            ThreadPoolExecutor(max_workers=1) as pool,
            pytest.raises(EvaluationSessionBusyError),
        ):
            pool.submit(one.evaluate).result()
    one.close()
    assert one.closed
    with pytest.raises(EvaluationSessionClosedError):
        one.evaluate()
    assert one.evaluate_fresh() == first
    assert system.to_dict() == before


@pytest.mark.parametrize("bad", [True, float("nan"), float("inf"), 1e308])
def test_prebackend_failure_keeps_session(bound, monkeypatch, bad):
    system, _, e = bound
    backend = e._evaluate_context
    calls = []

    def track(*args):
        calls.append(1)
        return backend(*args)

    monkeypatch.setattr(e, "_evaluate_context", track)
    coords = {s: system.coordinates.get(s).tolist() for s in system.topology.sites}
    coords[11][0] = bad
    with e.open_session() as session:
        with pytest.raises(EvaluationInputError):
            session.evaluate(coords)
        with pytest.raises(EvaluationInputError):
            session.evaluate(coordinate_unit="nm")
        assert not calls and not session.closed
        session.evaluate()
        assert calls == [1]


@pytest.mark.parametrize("defect", ["bond", "element", "mass", "stereo", "ids"])
@pytest.mark.parametrize("session_path", [False, True])
def test_binding_precedes_evaluation(bound, monkeypatch, defect, session_path):
    system, _, e = bound
    other = deepcopy(system)
    if defect == "bond":
        other.topology.bonds.pop(next(iter(other.topology.bonds)))
    elif defect == "element":
        other.topology.sites[11].element = "F"
        other.topology.sites[11].atomic_number = 9
    elif defect == "mass":
        other.topology.sites[11].mass = 13
    elif defect == "stereo":
        other.topology.sites[11].metadata["cip_label"] = "R"
    else:
        other.topology.sites.pop(11)
    before = deepcopy(other.to_dict())
    target = e.open_session() if session_path else e
    monkeypatch.setattr(
        target, "evaluate", lambda *a, **k: pytest.fail("Binding bypass")
    )
    try:
        with pytest.raises(MinimizationInputError):
            minimize_geometry(other, target)
    finally:
        if session_path:
            target.close()
    assert other.to_dict() == before


@pytest.mark.parametrize("budget", [1, 2, 20])
def test_budgets_fresh_context_application(bound, monkeypatch, budget):
    import openmm as mm

    system, spec, e = bound
    contexts = []
    calls = []
    ctor = mm.Context
    backend = e._evaluate_context

    def context(*a, **k):
        obj = ctor(*a, **k)
        contexts.append(obj)
        return obj

    def evaluate(*args):
        calls.append(args[1])
        return backend(*args)

    monkeypatch.setattr(mm, "Context", context)
    monkeypatch.setattr(e, "_evaluate_context", evaluate)
    before = deepcopy(system.to_dict())
    with e.open_session() as session:
        minimum = minimize_geometry(
            system,
            session,
            MinimizationOptions(
                max_evaluations=budget, force_tolerance=1e6 if budget == 20 else 0.1
            ),
        )
    minimum.validate_integrity()
    assert minimum.evaluations == len(calls) <= budget
    assert len(contexts) == (1 if budget == 1 else 2)
    if budget == 20:
        assert minimum.converged and minimum.final_evaluation_verified
        copied = minimum.to_system(system)
        copied.validate()
        spec.validate_integrity(copied)
    else:
        assert not minimum.converged
        with pytest.raises(MinimizationInputError):
            minimum.to_system(system)
    assert calls[-1] is contexts[-1]
    assert system.to_dict() == before


def test_backend_failure_recovery(bound, monkeypatch):
    system, _, e = bound
    calls = []
    backend = e._evaluate_context

    def fail(*args):
        calls.append(args[1])
        result = backend(*args)
        if len(calls) == 2:
            raise RuntimeError("after installation")
        return result

    monkeypatch.setattr(e, "_evaluate_context", fail)
    with e.open_session() as session:
        result = minimize_geometry(
            system, session, MinimizationOptions(max_evaluations=20)
        )
        assert session.closed
    assert result.evaluations == len(calls) == 3
    assert calls[0] is calls[1] and calls[-1] is not calls[0]
    assert result.termination_reason == "evaluation_failed" and not result.converged
    assert result.final_evaluation_verified
    result.validate_integrity()
    with pytest.raises(MinimizationInputError):
        result.to_system(system)


def test_final_failure_and_cleanup(bound, monkeypatch):
    system, _, e = bound
    with e.open_session() as session:

        def fail(*a, **k):
            raise EvaluationError("independent verification failure")

        monkeypatch.setattr(session, "evaluate_fresh", fail)
        result = minimize_geometry(
            system, session, MinimizationOptions(force_tolerance=1e6)
        )
    assert (
        result.evaluations == 2
        and not result.converged
        and not result.final_evaluation_verified
    )
    result.validate_integrity()
    with pytest.raises(MinimizationInputError):
        result.to_system(system)
    diagnostic = result.to_system(system, allow_unconverged=True)
    assert diagnostic.metadata["coordinate_source"] == "local_minimization_diagnostic"
    session = e.open_session()
    close = session._resources.close

    def bad_close():
        close()
        raise RuntimeError("secondary cleanup")

    monkeypatch.setattr(session._resources, "close", bad_close)
    with pytest.raises(ValueError, match="original") as caught, session:
        raise ValueError("original")
    assert "secondary cleanup" in caught.value.__notes__[0]
    assert session.closed and session._resources is None
    assert session.evaluate_fresh() == e.evaluate()


def test_lazy_session_import():
    import subprocess
    import sys

    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; sys.modules.update({k:None for k in ['openmm','scipy','rdkit','parmed','foyer']}); from island.evaluation.pcff_session import PCFFEvaluationSession",
        ],
        check=True,
    )
