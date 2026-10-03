"""Synthetic H2 contract tests, not OPLS scientific dynamics validation."""

from copy import deepcopy

import pytest
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
    DynamicsInputError,
    EvaluationInputError,
    MinimizationInputError,
    OPLSAssignmentError,
)
from island.minimization import MinimizationOptions, minimize_geometry


@pytest.fixture
def synthetic(monkeypatch):
    return _synthetic.__wrapped__(monkeypatch)


ENTRIES = ["minimize", "nve", "langevin", "segment_nve", "segment_langevin"]


def invoke(name, s, e):
    velocities = {i: (0.0, 0.0, 0.0) for i in s.topology.sites}
    nve = DynamicsOptions(0.01, 1, 3, 2)
    thermal = LangevinOptions(0.01, 0.0, 0.0, 123, 1, 3, 2)
    if name == "minimize":
        return minimize_geometry(
            s, e, MinimizationOptions(max_iterations=10, max_evaluations=30)
        )
    if name == "nve":
        return run_nve(s, e, velocities, nve)
    if name == "langevin":
        return run_langevin(s, e, velocities, thermal)
    return run_dynamics_segment(
        s, e, velocities, thermal if name == "segment_langevin" else nve
    )


def change(s, defect):
    s = deepcopy(s)
    if defect == "graph":
        s.topology.bonds.clear()
    elif defect == "element":
        for a in s.topology.sites.values():
            a.element, a.atomic_number, a.mass = "F", 9, 18.998
    elif defect == "mass":
        for a in s.topology.sites.values():
            a.mass = 2.0
    else:
        s.topology.sites[13].metadata["cip_label"] = "R"
    return s


@pytest.mark.parametrize("entry", ENTRIES)
@pytest.mark.parametrize("defect", ["graph", "element", "mass", "stereo"])
def test_bound_rejection_before_first_evaluation(synthetic, monkeypatch, entry, defect):
    pytest.importorskip("openmm")
    s, src, p = synthetic
    e = OPLSSinglePointEvaluator(s, p, src)
    previous = e.evaluate()
    previous_data = deepcopy(previous)
    altered = change(s, defect)
    before = deepcopy(altered.to_dict())
    with pytest.raises((EvaluationInputError, OPLSAssignmentError)):
        e.validate_system(altered)
    calls = []
    validate = e.validate_system

    def checked(system):
        calls.append("binding")
        return validate(system)

    def forbidden(*args, **kwargs):
        calls.append("energy")
        pytest.fail("Incompatible bound potential reached energy evaluation")

    monkeypatch.setattr(e, "validate_system", checked)
    monkeypatch.setattr(e, "evaluate", forbidden)
    with pytest.raises((MinimizationInputError, DynamicsInputError)):
        invoke(entry, altered, e)
    assert calls == ["binding"]
    assert altered.to_dict() == before
    assert previous == previous_data
    assert s.to_dict() == e._system.to_dict()


def test_resume_reaches_binding_with_valid_checkpoint(synthetic, monkeypatch):
    pytest.importorskip("openmm")
    s, src, p = synthetic
    e = OPLSSinglePointEvaluator(s, p, src)
    changed = change(s, "element")

    # An explicit analytical test potential is allowed to describe any graph.
    # Its numerical function happens to be this toy H2 spring. This creates a
    # fully valid F2 checkpoint with identical evaluation identities, without
    # changing checkpoint checksums, compatibility fields or numerical records.
    class Analytical:
        def evaluate(self, coordinates, *, coordinate_unit="angstrom"):
            return e.evaluate(coordinates, coordinate_unit=coordinate_unit)

    segment = invoke("segment_nve", changed, Analytical())
    cp = create_dynamics_checkpoint(segment)
    cp.validate_integrity()
    saved = cp.payload
    before = deepcopy(changed.to_dict())
    calls = []
    validate = e.validate_system

    def checked(system):
        calls.append("binding")
        return validate(system)

    def forbidden(*args, **kwargs):
        calls.append("energy")
        pytest.fail("Resume skipped bound validation")

    monkeypatch.setattr(e, "validate_system", checked)
    monkeypatch.setattr(e, "evaluate", forbidden)
    with pytest.raises(DynamicsInputError):
        resume_dynamics(cp, changed, e, DynamicsSegmentOptions(1, 3, 2))
    assert calls == ["binding"]
    assert cp.payload == saved
    assert changed.to_dict() == before


@pytest.mark.parametrize("entry", ENTRIES + ["resume"])
def test_compatible_coordinates_and_fresh_verification(synthetic, monkeypatch, entry):
    pytest.importorskip("openmm")
    s, src, p = synthetic
    e = OPLSSinglePointEvaluator(s, p, src)
    moved = deepcopy(s)
    moved.coordinates.set(39, (1.1, 0.0, 0.0))
    before = deepcopy(moved.to_dict())
    identity = e.model_fingerprint, e.parameter_fingerprint
    p.validate_integrity(moved, src)
    calls = []
    original = e.evaluate_fresh

    def fresh(*args, **kwargs):
        calls.append("fresh")
        return original(*args, **kwargs)

    monkeypatch.setattr(e, "evaluate_fresh", fresh)
    if entry == "resume":
        cp = create_dynamics_checkpoint(invoke("segment_nve", moved, e))
        calls.clear()
        result = resume_dynamics(cp, moved, e, DynamicsSegmentOptions(1, 3, 2))
        assert calls == ["fresh", "fresh"]  # startup and final
    else:
        result = invoke(entry, moved, e)
        assert calls == ["fresh"] * (2 if entry.startswith("segment_") else 1)
    result.validate_integrity()
    if entry == "minimize":
        assert result.converged and result.final_evaluation_verified
    assert moved.to_dict() == before
    assert identity == (e.model_fingerprint, e.parameter_fingerprint)


def test_fresh_really_constructs_new_contexts(synthetic, monkeypatch):
    mm = pytest.importorskip("openmm")
    s, src, p = synthetic
    e = OPLSSinglePointEvaluator(s, p, src)
    constructor = mm.Context
    calls = []

    def context(*args, **kwargs):
        calls.append("context")
        return constructor(*args, **kwargs)

    monkeypatch.setattr(mm, "Context", context)
    first = e.evaluate()
    second = e.evaluate_fresh()
    assert calls == ["context", "context"]
    assert first == second
