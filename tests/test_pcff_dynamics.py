"""PCFF interoperability: synthetic software fixture, not scientific charge data."""
# ruff: noqa: F811 -- imported pytest fixture injection

from copy import deepcopy
from dataclasses import replace

import pytest
from test_pcff_evaluation import bound as _bound  # noqa: F401
from test_pcff_evaluation import parameters, synthetic  # noqa: F401

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
from island.exceptions import (
    DynamicsCheckpointCompatibilityError,
    DynamicsInputError,
    InvalidDynamicsCheckpointError,
)


@pytest.fixture
def bound(_bound):
    s, _spec, evaluator = _bound
    return s, evaluator


def change(system, defect):
    s = deepcopy(system)
    first = next(iter(s.topology.sites))
    if defect == "graph":
        s.topology.bonds.pop(next(iter(s.topology.bonds)))
    elif defect == "element":
        for site in s.topology.sites.values():
            if site.element == "C":
                site.element = "F"
                site.atomic_number = 9
                site.mass = 18.998
    elif defect == "mass":
        s.topology.sites[first].mass += 0.01
    elif defect == "stereo":
        s.topology.sites[first].metadata["cip_label"] = "R"
    elif defect == "ids":
        s.topology.sites.pop(first)
    return s


def options(thermal, steps=2, budget=None):
    args = {
        "timestep_fs": 0.01,
        "steps": steps,
        "max_evaluations": steps + 2 if budget is None else budget,
        "max_frames": steps + 1,
    }
    return (
        LangevinOptions(
            **args, temperature_kelvin=300, friction_per_ps=5, thermostat_seed=37
        )
        if thermal
        else DynamicsOptions(**args)
    )


def velocities(s):
    return {i: (0.0, 0.0, 0.0) for i in s.topology.sites}


@pytest.mark.parametrize("thermal", [False, True])
@pytest.mark.parametrize("segmented", [False, True])
@pytest.mark.parametrize("defect", ["graph", "element", "mass", "stereo", "ids"])
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
    assert (
        last.payload["origin"]["initial_total_energy"]
        == whole.payload["origin"]["initial_total_energy"]
    )
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
        assert result.completed_steps == 0
        assert dict(result.final_state.velocities) == velocities(s)
        assert dict(result.final_state.coordinates) == {
            i: tuple(s.coordinates.get(i)) for i in s.topology.sites
        }
        assert result.payload["counters"]["random_steps"] == int(thermal)
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


@pytest.mark.parametrize("policy_change", [False, True])
def test_actual_changed_charge_or_pair_model_on_resume(
    _bound, monkeypatch, policy_change
):
    from test_pcff_class2 import assigned

    from island.evaluation import PCFFSinglePointEvaluator
    from island.forcefields.pcff import (
        automatic,
        define_pcff_model,
        model,
        source,
        special_pair_policy,
    )
    from island.forcefields.pcff.source import PCFFSource, digest

    s, spec, e = _bound
    with e.open_session() as session:
        cp = create_dynamics_checkpoint(
            run_dynamics_segment(s, session, velocities(s), options(True))
        )
    if policy_change:
        other = define_pcff_model(
            spec.assignment,
            special_pairs=special_pair_policy(lj=(0, 0, 0.5), coulomb=(0, 0, 1)),
        )
    else:
        raw = spec.assignment.source.raw.replace(b"c h -0.1 0.1", b"c h -0.11 0.11")
        assert raw != spec.assignment.source.raw
        sha = digest(raw)
        monkeypatch.setitem(source.PIN, "sha256", sha)
        monkeypatch.setitem(automatic.PROFILE, "source_sha256", sha)
        monkeypatch.setitem(model.PROFILE, "frc_sha256", sha)
        other = define_pcff_model(
            assigned(s, PCFFSource(raw, sha)),
            special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
        )
    evaluator = PCFFSinglePointEvaluator(s, other)
    with (
        evaluator.open_session() as session,
        pytest.raises(DynamicsCheckpointCompatibilityError) as caught,
    ):
        resume_dynamics(cp, s, session, DynamicsSegmentOptions(2, 4, 3))
    assert caught.value.evaluations == 1


@pytest.mark.parametrize("thermal", [False, True])
def test_frame_budget_no_trial_and_resume_startup_failure(bound, monkeypatch, thermal):
    s, e = bound
    with e.open_session() as session:
        limited = run_dynamics_segment(
            s, session, velocities(s), replace(options(thermal), max_frames=1)
        )
    assert limited.completed_steps == 0 and limited.evaluations == 2
    assert limited.termination_reason == "maximum_frames"
    assert limited.payload["counters"]["normal_draws"] == 0
    cp = create_dynamics_checkpoint(limited)
    before = cp.payload
    original = e._evaluate_context

    def changed(*args):
        record = original(*args)
        components = dict(record.energy_components)
        components[next(iter(components))] += 1
        return replace(
            record,
            potential_energy=record.potential_energy + 1,
            energy_components=components,
        )

    monkeypatch.setattr(e, "_evaluate_context", changed)
    with e.open_session() as session:
        result = resume_dynamics(cp, s, session, DynamicsSegmentOptions(2, 4, 3))
    result.validate_integrity()
    assert result.termination_reason == "startup_verification_failed"
    assert result.evaluations == 1 and result.completed_steps == 0
    assert result.final_state == cp.state
    assert (
        result.payload["rng"] == before["rng"]
        and result.payload["origin"] == before["origin"]
    )
    assert cp.payload == before
    with pytest.raises(InvalidDynamicsCheckpointError):
        create_dynamics_checkpoint(result)
    with pytest.raises(DynamicsInputError):
        result.to_system(s, allow_incomplete=True)


@pytest.fixture
def acceptance(monkeypatch):
    from pathlib import Path

    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    import validate_pcff_dynamics

    return validate_pcff_dynamics


@pytest.mark.parametrize("defect", ["missing", "altered", "inventory"])
def test_external_file_gates_before_backend(acceptance, tmp_path, monkeypatch, defect):
    from island.workflows import storage

    for name in acceptance.FILES:
        (tmp_path / name).write_text("{}")
    manifest = {
        "schema": "pcff_acceptance_inputs_v1",
        "files": {
            name: acceptance.digest(tmp_path / name) for name in acceptance.FILES
        },
    }
    if defect == "inventory":
        manifest["files"]["other"] = "0" * 64
    storage.publish(tmp_path / "inputs.json", storage.json_bytes(manifest))
    if defect == "missing":
        (tmp_path / "model.json").unlink()
    if defect == "altered":
        (tmp_path / "model.json").write_text("altered")
    monkeypatch.setattr(
        acceptance,
        "load_pcff_source",
        lambda *a: pytest.fail("Invalid input reached backend"),
    )
    with pytest.raises((ValueError, FileNotFoundError)):
        acceptance.reconstruct(tmp_path)


def test_valid_reconstruction_and_rechecksummed_mismatch(
    _bound, acceptance, tmp_path, monkeypatch
):
    import sys

    import numpy as np

    from island.dynamics import initialize_velocities
    from island.minimization import MinimizationOptions, minimize_geometry
    from island.workflows import storage
    from island.workflows.bundle import record

    s, spec, e = _bound
    with e.open_session() as session:
        minimum = minimize_geometry(
            s, session, MinimizationOptions(max_iterations=500, max_evaluations=1000)
        )
    assert minimum.converged and minimum.final_fmax <= 0.1
    optimized = minimum.to_system(s)
    init = initialize_velocities(optimized, temperature_kelvin=300, seed=78123)
    directory = tmp_path / "bundle"
    acceptance.bundle(
        directory,
        s,
        optimized,
        spec.assignment,
        spec,
        minimum,
        spec.assignment.source,
        init,
    )
    for module in ("rdkit", "parmed", "foyer", "scipy"):
        monkeypatch.setitem(sys.modules, module, None)
    before = np.random.get_state()
    restored, evaluator, loaded = acceptance.reconstruct(directory)
    assert restored.to_dict() == optimized.to_dict() and loaded == init
    assert evaluator.model_fingerprint == e.model_fingerprint
    after = np.random.get_state()
    assert all(np.array_equal(a, b) for a, b in zip(before, after))
    with evaluator.open_session() as session:
        cp = create_dynamics_checkpoint(
            run_dynamics_segment(restored, session, loaded.velocities, options(True))
        )
    acceptance.validate_origin(cp, restored, loaded)
    other_initialization = initialize_velocities(
        restored, temperature_kelvin=300, seed=17
    )
    with pytest.raises(ValueError, match="origin/initialization"):
        acceptance.validate_origin(cp, restored, other_initialization)
    # Independently valid initialization belonging to a different mass inventory.
    changed = deepcopy(optimized)
    changed.topology.sites[11].mass += 0.01
    bad = initialize_velocities(changed, temperature_kelvin=300, seed=78123)
    storage.publish(
        directory / "initialization.json",
        storage.json_bytes(storage.encode(record(bad))),
        replace=True,
    )
    manifest = storage.read_json(directory / "inputs.json")
    manifest["files"]["initialization.json"] = acceptance.digest(
        directory / "initialization.json"
    )
    storage.publish(
        directory / "inputs.json", storage.json_bytes(manifest), replace=True
    )
    with pytest.raises(ValueError, match="Initialization/system identity"):
        acceptance.reconstruct(directory)
