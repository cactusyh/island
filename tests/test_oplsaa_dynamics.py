"""OPLS session interoperability: synthetic spring, not scientific OPLS data."""

from copy import deepcopy
from dataclasses import replace

import pytest
from test_oplsaa_binding import change
from test_oplsaa_parameters import synthetic as _synthetic

from island.dynamics import (
    DynamicsOptions,
    DynamicsSegmentOptions,
    LangevinOptions,
    create_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    run_langevin,
    run_nve,
)
from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.exceptions import (
    DynamicsCheckpointCompatibilityError,
    DynamicsInputError,
    InvalidDynamicsCheckpointError,
)


@pytest.fixture
def bound(monkeypatch):
    pytest.importorskip("openmm")
    s, src, p = _synthetic.__wrapped__(monkeypatch)
    return s, OPLSSinglePointEvaluator(s, p, src)


def options(thermal, steps=2, budget=None):
    args = {
        "timestep_fs": 0.01,
        "steps": steps,
        "max_evaluations": steps + 2 if budget is None else budget,
        "max_frames": steps + 1,
    }
    return (
        LangevinOptions(
            **args, temperature_kelvin=0, friction_per_ps=5, thermostat_seed=37
        )
        if thermal
        else DynamicsOptions(**args)
    )


def velocities(s):
    return {i: (0.0, 0.0, 0.0) for i in s.topology.sites}


@pytest.mark.parametrize("thermal", [False, True])
@pytest.mark.parametrize("segmented", [False, True])
@pytest.mark.parametrize("defect", ["graph", "element", "mass", "stereo"])
def test_session_binding_before_dynamics(
    bound, monkeypatch, thermal, segmented, defect
):
    s, e = bound
    changed = change(s, defect)
    before = deepcopy(changed.to_dict())
    method = run_dynamics_segment if segmented else run_langevin if thermal else run_nve
    with e.open_session() as session:
        monkeypatch.setattr(
            e, "_evaluate_context", lambda *a: pytest.fail("Binding bypass")
        )
        with pytest.raises(DynamicsInputError):
            method(changed, session, velocities(s), options(thermal))
    assert changed.to_dict() == before


@pytest.mark.parametrize("thermal", [False, True])
def test_budget_and_split_contract(bound, thermal):
    s, e = bound
    with e.open_session() as session:
        direct = (run_langevin if thermal else run_nve)(
            s, session, velocities(s), options(thermal, 4)
        )
        whole = run_dynamics_segment(s, session, velocities(s), options(thermal, 4))
        first = run_dynamics_segment(s, session, velocities(s), options(thermal, 2))
        limited = run_dynamics_segment(
            s, session, velocities(s), options(thermal, 2, 2)
        )
    assert (
        direct.final_evaluation_verified
        and direct.evaluations == whole.evaluations == 6
    )
    assert limited.completed_steps == 0 and limited.evaluations == 2
    limited.validate_integrity()
    create_dynamics_checkpoint(limited)  # eligible verified pre-trial budget stop
    cp = create_dynamics_checkpoint(first)
    with e.open_session() as session:
        last = resume_dynamics(cp, s, session, DynamicsSegmentOptions(2, 4, 3))
    assert last.final_state == whole.final_state
    assert last.payload["rng"] == whole.payload["rng"]
    assert last.payload["origin"] == whole.payload["origin"]
    assert last.payload["counters"]["evaluations"] == 8


@pytest.mark.parametrize("stage", ["startup", "trial", "final"])
@pytest.mark.parametrize("thermal", [False, True])
def test_failure_accounting_and_ineligibility(bound, monkeypatch, stage, thermal):
    s, e = bound
    original = e._evaluate_context
    calls = []
    failure = {"startup": 1, "trial": 2, "final": 4}[stage]

    def backend(*args):
        calls.append(args[1])
        result = original(*args)
        if len(calls) == failure:
            raise RuntimeError("injected backend failure")
        return result

    monkeypatch.setattr(e, "_evaluate_context", backend)
    with e.open_session() as session:
        result = run_dynamics_segment(s, session, velocities(s), options(thermal))
        assert session.closed == (stage == "trial")
    result.validate_integrity()
    assert (
        result.evaluations
        == len(calls)
        == {"startup": 1, "trial": 3, "final": 4}[stage]
    )
    assert not result.completed
    if stage == "trial":
        assert result.payload["diagnostic"]["final_verified"]
        assert calls[-1] is not calls[1]
    with pytest.raises(InvalidDynamicsCheckpointError):
        create_dynamics_checkpoint(result)
    with pytest.raises(DynamicsInputError):
        result.to_system(s)
    if stage == "startup":
        with pytest.raises(DynamicsInputError):
            result.to_system(s, allow_incomplete=True)


@pytest.mark.parametrize("field", ["model_fingerprint", "parameter_fingerprint"])
def test_resume_identity_rejected_at_fresh_startup(bound, monkeypatch, field):
    s, e = bound
    with e.open_session() as session:
        cp = create_dynamics_checkpoint(
            run_dynamics_segment(s, session, velocities(s), options(True))
        )
    before = cp.payload
    original = e._evaluate_context
    calls = []

    def backend(*args):
        calls.append(1)
        return replace(original(*args), **{field: "incompatible"})

    monkeypatch.setattr(e, "_evaluate_context", backend)
    with (
        e.open_session() as session,
        pytest.raises(DynamicsCheckpointCompatibilityError) as error,
    ):
        resume_dynamics(cp, s, session, DynamicsSegmentOptions(2, 4, 3))
    assert error.value.evaluations == len(calls) == 1
    assert cp.payload == before


def test_session_resume_binding_after_valid_external_checkpoint(bound, monkeypatch):
    s, e = bound
    changed = change(s, "element")

    class Analytical:
        def evaluate(self, coordinates, *, coordinate_unit="angstrom"):
            return e.evaluate(coordinates, coordinate_unit=coordinate_unit)

    cp = create_dynamics_checkpoint(
        run_dynamics_segment(changed, Analytical(), velocities(s), options(False))
    )
    cp.validate_integrity()
    with e.open_session() as session:
        monkeypatch.setattr(
            e, "_evaluate_context", lambda *a: pytest.fail("Resume binding bypass")
        )
        with pytest.raises(DynamicsInputError):
            resume_dynamics(cp, changed, session, DynamicsSegmentOptions(2, 4, 3))
