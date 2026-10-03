"""Real Context lifecycle with explicitly synthetic spring parameters."""

from concurrent.futures import ThreadPoolExecutor
from copy import copy, deepcopy

import pytest
from test_oplsaa_binding import change
from test_oplsaa_parameters import synthetic as _synthetic

from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.exceptions import (
    EvaluationError,
    EvaluationInputError,
    EvaluationSessionBusyError,
    EvaluationSessionClosedError,
    MinimizationInputError,
)
from island.minimization import MinimizationOptions, minimize_geometry


@pytest.fixture
def bound(monkeypatch):
    pytest.importorskip("openmm")
    s, source, p = _synthetic.__wrapped__(monkeypatch)
    return s, source, p, OPLSSinglePointEvaluator(s, p, source)


def test_frames_identity_ownership_and_independent_sessions(bound):
    s, _, _, e = bound
    before = deepcopy(s.to_dict())
    a = {13: (0.0, 0.0, 0.0), 39: (1.2, 0.0, 0.0)}
    b = {39: (1.4, 0.1, 0.0), 13: (0.0, 0.0, 0.0)}
    first = e.evaluate(a)
    with e.open_session() as one, e.open_session() as two:
        assert one._resources.context is not two._resources.context
        assert one.evaluate(a) == two.evaluate(a) == first
        prior = one.evaluate(b)
        assert prior == e.evaluate(b)
        assert two.evaluate(a) == one.evaluate(a) == first
        assert one.evaluate_fresh(b) == prior
        assert one.model_fingerprint == e.model_fingerprint
        assert one.parameter_fingerprint == e.parameter_fingerprint
    assert e.evaluate() == first
    assert s.to_dict() == before
    assert prior == e.evaluate(b)


def test_lifecycle_and_no_stepping(bound, monkeypatch):
    import openmm as mm

    _, _, _, e = bound
    monkeypatch.setattr(
        mm.VerletIntegrator, "step", lambda *a: pytest.fail("Integrator stepped")
    )
    with e.open_session() as session:
        session.evaluate()
        for copier in (copy, deepcopy):
            with pytest.raises(EvaluationInputError):
                copier(session)
        with pytest.raises(EvaluationSessionBusyError), session:
            pass
        session._lock.acquire()
        try:
            with pytest.raises(EvaluationSessionBusyError):
                session.evaluate()
        finally:
            session._lock.release()
        with (
            ThreadPoolExecutor(max_workers=1) as pool,
            pytest.raises(EvaluationSessionBusyError),
        ):
            pool.submit(session.evaluate).result()
    session.close()
    assert session.closed and session._resources is None
    with pytest.raises(EvaluationSessionClosedError):
        session.evaluate()
    assert session.evaluate_fresh() == e.evaluate()


@pytest.mark.parametrize(
    "bad",
    [
        {13: (0, 0, 0)},
        {13: (True, 0, 0), 39: (1, 0, 0)},
        {13: (float("nan"), 0, 0), 39: (1, 0, 0)},
        {13: (1e308, 0, 0), 39: (-1e308, 0, 0)},
        {13: (0, 0, 0), 39: (0, 0, 0)},
    ],
)
def test_input_failure_before_backend_mutation(bound, monkeypatch, bad):
    _, _, _, e = bound
    original = e._evaluate_context
    calls = []

    def backend(*args):
        calls.append(1)
        return original(*args)

    monkeypatch.setattr(e, "_evaluate_context", backend)
    with e.open_session() as session:
        with pytest.raises(EvaluationInputError):
            session.evaluate(bad)
        with pytest.raises(EvaluationInputError):
            session.evaluate(coordinate_unit="nm")
        assert not calls and not session.closed
        session.evaluate()
        assert calls == [1]


@pytest.mark.parametrize("defect", ["graph", "element", "mass", "stereo"])
def test_minimizer_binding_before_evaluation(bound, monkeypatch, defect):
    s, _, _, e = bound
    other = change(s, defect)
    before = deepcopy(other.to_dict())
    with e.open_session() as session:
        monkeypatch.setattr(
            session, "evaluate", lambda *a, **k: pytest.fail("Binding bypass")
        )
        with pytest.raises(MinimizationInputError):
            minimize_geometry(other, session, MinimizationOptions())
    assert other.to_dict() == before


def test_backend_failure_closes_but_fresh_verification_survives(bound, monkeypatch):
    s, _, _, e = bound
    original = e._evaluate_context
    calls = []

    def backend(*args):
        calls.append(args[1])
        result = original(*args)
        if len(calls) == 2:
            raise RuntimeError("injected backend failure after state installation")
        return result

    monkeypatch.setattr(e, "_evaluate_context", backend)
    with e.open_session() as session:
        result = minimize_geometry(s, session, MinimizationOptions(max_evaluations=20))
        assert session.closed
    assert result.evaluations == len(calls) == 3
    assert calls[0] is calls[1] and calls[2] is not calls[1]
    assert not result.converged and result.termination_reason == "evaluation_failed"
    assert result.final_evaluation_verified
    result.validate_integrity()
    with pytest.raises(MinimizationInputError):
        result.to_system(s)
    result.to_system(s, allow_unconverged=True).validate()


def test_original_exception_survives_cleanup_failure(bound, monkeypatch):
    _, _, _, e = bound
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


@pytest.mark.parametrize("budget", [1, 2, 30])
def test_minimization_counts_final_reservation_and_application(
    bound, monkeypatch, budget
):
    import openmm as mm

    s, src, p, e = bound
    contexts, calls = [], []
    context = mm.Context
    backend = e._evaluate_context

    def construct(*a, **kw):
        obj = context(*a, **kw)
        contexts.append(obj)
        return obj

    def evaluate(*a):
        calls.append(a[1])
        return backend(*a)

    monkeypatch.setattr(mm, "Context", construct)
    monkeypatch.setattr(e, "_evaluate_context", evaluate)
    before = deepcopy(s.to_dict())
    with e.open_session() as session:
        result = minimize_geometry(
            s, session, MinimizationOptions(max_evaluations=budget)
        )
    assert result.evaluations == len(calls) <= budget
    assert len(contexts) == (1 if budget == 1 else 2)
    assert calls[-1] is contexts[-1]
    if budget == 1:
        assert not result.final_evaluation_verified
    assert all(c is contexts[0] for c in calls[:-1])
    if budget == 30:
        assert result.converged and result.final_evaluation_verified
        copied = result.to_system(s)
        copied.validate()
        p.validate_integrity(copied, src)
        assert result.final_fmax <= 0.1
    else:
        assert not result.converged
        with pytest.raises(MinimizationInputError):
            result.to_system(s)
    assert s.to_dict() == before


@pytest.mark.parametrize("session_path", [False, True])
def test_initialization_failure_preserves_original(bound, monkeypatch, session_path):
    import openmm as mm

    from island.evaluation.openmm import _OpenMMResources
    from island.exceptions import EvaluationUnavailableError

    _, _, _, e = bound
    close = _OpenMMResources.close

    def failed_context(*args, **kwargs):
        raise RuntimeError("primary initialization failure")

    def failed_cleanup(self):
        close(self)
        raise RuntimeError("secondary disposal failure")

    monkeypatch.setattr(mm, "Context", failed_context)
    monkeypatch.setattr(_OpenMMResources, "close", failed_cleanup)
    with pytest.raises(
        (EvaluationError, EvaluationUnavailableError), match="primary initialization"
    ) as caught:
        e.open_session() if session_path else e.evaluate()
    assert "secondary disposal" in caught.value.__cause__.__notes__[0]
