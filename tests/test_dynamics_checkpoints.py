"""Data-only boundaries: exact analytical continuation, corruption and failure contracts."""

from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest
from test_dynamics import VELOCITIES, Harmonic, molecule, options
from test_langevin import options as thermal_options

from island.dynamics import (
    DynamicsCheckpoint,
    DynamicsSegment,
    DynamicsSegmentOptions,
    create_dynamics_checkpoint,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    run_langevin,
    run_nve,
    save_dynamics_checkpoint,
)
from island.dynamics import _checkpoint_data as data
from island.exceptions import (
    DynamicsCheckpointCompatibilityError,
    DynamicsCheckpointIOError,
    DynamicsInputError,
    InvalidDynamicsCheckpointError,
)


def budgets(n, **kw):
    return DynamicsSegmentOptions(
        n, kw.pop("max_evaluations", n + 2), kw.pop("max_frames", n + 1), **kw
    )


def rebuild(cp, change):
    p = cp.payload
    change(p)
    return DynamicsCheckpoint(data.canonical(p), data.checksum(p))


@pytest.mark.parametrize("thermal", [False, True])
def test_exact_three_segments_lossless_state_and_legacy_equivalence(tmp_path, thermal):
    system = molecule()
    before = deepcopy(system.to_dict())
    opts = thermal_options(12) if thermal else options(12)
    whole = run_dynamics_segment(system, Harmonic(), VELOCITIES, opts)
    legacy = (run_langevin if thermal else run_nve)(
        system, Harmonic(), VELOCITIES, opts
    )
    for a, b in zip(whole.frames, legacy.frames, strict=True):
        assert a == b
    segment = run_dynamics_segment(
        molecule(True),
        Harmonic(),
        dict(reversed(list(VELOCITIES.items()))),
        replace(opts, steps=4, max_evaluations=6, max_frames=5),
    )
    np.random.seed(51)
    global_state = deepcopy(np.random.get_state())
    for index in range(3):
        checkpoint = create_dynamics_checkpoint(segment)
        assert deepcopy(checkpoint) is checkpoint
        path = tmp_path / f"{index}.json"
        save_dynamics_checkpoint(checkpoint, path)
        restored = load_dynamics_checkpoint(path)
        assert restored == checkpoint
        if index < 2:
            moved = segment.to_system(system)
            segment = resume_dynamics(restored, moved, Harmonic(), budgets(4))
            assert (
                segment.payload["lineage"][-1]["parent_checksum"]
                == restored.content_checksum
            )
    assert segment.final_state == whole.final_state
    assert segment.payload["rng"] == whole.payload["rng"]
    assert segment.payload["origin"] == whole.payload["origin"]
    assert segment.payload["counters"]["evaluations"] == 18 and whole.evaluations == 14
    assert (
        segment.payload["max_abs_energy_deviation"]
        == whole.payload["max_abs_energy_deviation"]
    )
    assert segment.final_state.step == 12 and segment.completed_steps == 4
    assert system.to_dict() == before
    now = np.random.get_state()
    assert global_state[0] == now[0]
    np.testing.assert_array_equal(global_state[1], now[1])
    assert global_state[2:] == now[2:]


@pytest.mark.parametrize("temp,gamma", [(0, 3), (300, 0), (0, 0)])
def test_zero_limits_preserve_actual_rng(tmp_path, temp, gamma):
    o = thermal_options(6, temperature_kelvin=temp, friction_per_ps=gamma)
    full = run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, o)
    first = run_dynamics_segment(
        molecule(), Harmonic(), VELOCITIES, replace(o, steps=2)
    )
    cp = create_dynamics_checkpoint(first)
    final = resume_dynamics(cp, molecule(), Harmonic(), budgets(4))
    assert (
        final.final_state == full.final_state
        and final.payload["rng"] == full.payload["rng"]
    )
    assert final.payload["counters"]["normal_draws"] == 36
    assert cp.payload["rng"]["state"]["state"]["state"] > 2**64


@pytest.mark.parametrize("budget", [1, 2, 4])
def test_budget_boundary_eligibility(budget):
    result = run_dynamics_segment(
        molecule(), Harmonic(), VELOCITIES, thermal_options(10, max_evaluations=budget)
    )
    assert (
        result.termination_reason == "maximum_evaluations"
        and result.evaluations == budget
    )
    if budget == 1:
        with pytest.raises(InvalidDynamicsCheckpointError):
            create_dynamics_checkpoint(result)
    else:
        cp = create_dynamics_checkpoint(result)
        assert cp.payload["counters"]["normal_draws"] == cp.absolute_step * 6
        continued = resume_dynamics(cp, molecule(), Harmonic(), budgets(1))
        assert continued.completed


def test_frame_stop_failure_and_old_results_ineligible():
    r = run_dynamics_segment(
        molecule(),
        Harmonic(),
        VELOCITIES,
        thermal_options(10, max_frames=2, recording_interval=3),
    )
    assert r.final_state.step == 3 and r.termination_reason == "maximum_frames"
    create_dynamics_checkpoint(r)
    for fail in (1, 2, 4):
        r = run_dynamics_segment(
            molecule(), Harmonic(fail_call=fail), VELOCITIES, thermal_options(2)
        )
        with pytest.raises(InvalidDynamicsCheckpointError):
            create_dynamics_checkpoint(r)
    old = run_langevin(molecule(), Harmonic(), VELOCITIES, thermal_options(2))
    with pytest.raises(InvalidDynamicsCheckpointError):
        create_dynamics_checkpoint(old)


def test_startup_reverification_before_noise_and_guard_reference():
    first = run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, thermal_options(2))
    cp = create_dynamics_checkpoint(first)
    r = resume_dynamics(cp, molecule(), Harmonic(fail_call=1), budgets(3))
    assert r.termination_reason == "startup_verification_failed" and r.evaluations == 1
    assert r.payload["rng"] == cp.payload["rng"] and r.final_state == first.final_state
    with pytest.raises(InvalidDynamicsCheckpointError):
        create_dynamics_checkpoint(r)
    # Guard chosen before execution: initially below limit, eventually crossed.
    o = options(200, dt=2, max_energy_deviation=1e-5)
    full = run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, o)
    initial = run_dynamics_segment(
        molecule(), Harmonic(), VELOCITIES, replace(o, steps=1)
    )
    assert initial.completed
    split = resume_dynamics(
        create_dynamics_checkpoint(initial), molecule(), Harmonic(), budgets(199)
    )
    assert (
        split.termination_reason == full.termination_reason == "energy_guard_exceeded"
    )
    assert split.final_state == full.final_state
    assert split.payload["origin"] == full.payload["origin"]


@pytest.mark.parametrize(
    "field",
    ["model_fingerprint", "parameter_fingerprint", "backend_version", "settings"],
)
def test_model_backend_changes_typed_before_propagation(field):
    cp = create_dynamics_checkpoint(
        run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, options(1))
    )

    class Changed(Harmonic):
        def evaluate(self, *a, **kw):
            r = super().evaluate(*a, **kw)
            return replace(
                r,
                **{
                    field: {"precision": "changed"}
                    if field == "settings"
                    else "changed"
                },
            )

    potential = Changed()
    with pytest.raises(DynamicsCheckpointCompatibilityError) as error:
        resume_dynamics(cp, molecule(), potential, budgets(3))
    assert potential.calls == error.value.evaluations == 1


def test_system_environment_and_physical_mismatch():
    cp = create_dynamics_checkpoint(
        run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, options(1))
    )
    system = molecule()
    system.topology.sites[17].mass += 1
    with pytest.raises(DynamicsCheckpointCompatibilityError):
        resume_dynamics(cp, system, Harmonic(), budgets(1))
    p = cp.payload
    p["environment"]["numpy"] = "unsupported"
    p["origin"]["trajectory_fingerprint"] = data.trajectory_identity(p)
    changed = DynamicsCheckpoint(data.canonical(p), data.checksum(p))
    with pytest.raises(DynamicsCheckpointCompatibilityError):
        resume_dynamics(changed, molecule(), Harmonic(), budgets(1))
    with pytest.raises(DynamicsInputError):
        resume_dynamics(cp, molecule(), Harmonic(), options(1, dt=1))
    with pytest.raises(InvalidDynamicsCheckpointError):
        rebuild(cp, lambda p: p["physical"].update(timestep_fs=8))


@pytest.mark.parametrize(
    "mutation",
    [
        lambda p: p.update(schema="v99"),
        lambda p: p["units"].update(coordinates="nm"),
        lambda p: p["state"]["coordinates"].append(p["state"]["coordinates"][0]),
        lambda p: p["state"]["coordinates"][0][1].__setitem__(0, float("inf")),
        lambda p: p.update(rng=None),
        lambda p: p["rng"]["state"]["state"].update(state=2**128),
        lambda p: p["rng"]["state"]["state"].update(inc=2),
        lambda p: p["rng"]["state"].update(has_uint32=True),
        lambda p: p["counters"].update(normal_draws=123),
        lambda p: p["state"].update(time_ps=88),
        lambda p: p.update(final_evaluation=None),
        lambda p: p["origin"].update(initial_total_energy=0),
    ],
)
def test_malformed_resigned_payload_rejected(mutation):
    cp = create_dynamics_checkpoint(
        run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, thermal_options(1))
    )
    with pytest.raises((InvalidDynamicsCheckpointError, ValueError)):
        rebuild(cp, mutation)


def test_persistence_corruption_no_overwrite_atomic_failure(tmp_path, monkeypatch):
    cp = create_dynamics_checkpoint(
        run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, thermal_options(1))
    )
    path = tmp_path / "checkpoint.json"
    save_dynamics_checkpoint(cp, path)
    original = path.read_bytes()
    with pytest.raises(DynamicsCheckpointIOError):
        save_dynamics_checkpoint(cp, path)

    def fail(*a):
        raise OSError("interrupted before replacement")

    monkeypatch.setattr("island.dynamics.checkpoints.os.replace", fail)
    with pytest.raises(DynamicsCheckpointIOError):
        save_dynamics_checkpoint(cp, path, replace=True)
    assert path.read_bytes() == original and load_dynamics_checkpoint(path) == cp
    assert len(list(tmp_path.iterdir())) == 1
    for text in (
        "{",
        '{"payload":null,"payload":null,"sha256":"x"}',
        original.decode().replace('"sha256":', '"sha256":"bad","other":', 1),
        "NaN",
        '{"a":Infinity}',
    ):
        path.write_text(text)
        with pytest.raises(InvalidDynamicsCheckpointError):
            load_dynamics_checkpoint(path)


def test_segment_reconstruction_and_input_ownership():
    r = run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, options(2))
    p = r.payload
    p["counters"]["evaluations"] = 999
    with pytest.raises(InvalidDynamicsCheckpointError):
        DynamicsSegment(data.canonical(p), data.checksum(p))
    assert r.evaluations == 4 and deepcopy(r) is r


def test_startup_force_disagreement_is_one_call_diagnostic():
    cp = create_dynamics_checkpoint(
        run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, thermal_options(2))
    )

    class WrongForce(Harmonic):
        def evaluate(self, *args, **kwargs):
            result = super().evaluate(*args, **kwargs)
            return replace(
                result,
                forces={s: tuple(np.array(f) + 0.1) for s, f in result.forces.items()},
            )

    potential = WrongForce()
    failed = resume_dynamics(cp, molecule(), potential, budgets(2))
    assert (
        failed.termination_reason == "startup_verification_failed"
        and potential.calls == 1
    )
    assert failed.payload["rng"]["state"] == cp.rng_state
    with pytest.raises(InvalidDynamicsCheckpointError):
        create_dynamics_checkpoint(failed)


def test_atomic_explicit_replace_and_no_aliases(tmp_path):
    one = create_dynamics_checkpoint(
        run_dynamics_segment(molecule(), Harmonic(), VELOCITIES, thermal_options(1))
    )
    two = create_dynamics_checkpoint(
        resume_dynamics(one, molecule(), Harmonic(), budgets(1))
    )
    path = tmp_path / "state.json"
    save_dynamics_checkpoint(one, path)
    save_dynamics_checkpoint(two, path, replace=True)
    assert load_dynamics_checkpoint(path) == two
    external = two.rng_state
    external["state"]["state"] = 0
    assert two.rng_state["state"]["state"] != 0
    with pytest.raises(InvalidDynamicsCheckpointError):
        replace(two, content_checksum="bad")


def test_nve_and_thermal_metadata_provenance_only_changes_allowed():
    system = molecule()
    system.metadata["ambertools_preparation"] = {"historical_signature": "unchanged"}
    segment = run_dynamics_segment(system, Harmonic(), VELOCITIES, thermal_options(1))
    moved = segment.to_system(system)
    assert (
        moved.metadata["ambertools_preparation"]
        == system.metadata["ambertools_preparation"]
    )
    cp = create_dynamics_checkpoint(segment)
    continued = resume_dynamics(cp, moved, Harmonic(), budgets(1))
    assert continued.completed
    moved.metadata["ambertools_preparation"]["historical_signature"] = "changed"
    with pytest.raises(DynamicsCheckpointCompatibilityError):
        resume_dynamics(cp, moved, Harmonic(), budgets(1))


@pytest.mark.parametrize("guard", [0.0, 1e-9])
@pytest.mark.parametrize("offset", [5e-9, -5e-9, 5e-7, -5e-7])
def test_startup_energy_rejection_preserves_accepted_boundary(guard, offset):
    system = molecule()
    zero_velocities = {site: (0.0, 0.0, 0.0) for site in VELOCITIES}
    cp = create_dynamics_checkpoint(
        run_dynamics_segment(
            system, Harmonic(0), zero_velocities,
            options(1, max_energy_deviation=guard),
        )
    )
    checkpoint_before = deepcopy(cp.payload)
    system_before = deepcopy(system.to_dict())
    rng_before = deepcopy(np.random.get_state())

    class ShiftedEnergy(Harmonic):
        def evaluate(self, *args, **kwargs):
            record = super().evaluate(*args, **kwargs)
            return replace(
                record,
                potential_energy=record.potential_energy + offset,
                energy_components={"harmonic": offset},
            )

    potential = ShiftedEnergy(0)
    failed = resume_dynamics(cp, system, potential, budgets(2))
    failed.validate_integrity()
    p = failed.payload
    assert failed.termination_reason == "startup_verification_failed"
    assert failed.completed_steps == 0
    assert failed.evaluations == potential.calls == 1
    assert failed.final_state == cp.state
    assert failed.frames == (cp.state,)
    for key in ("state", "origin", "max_abs_energy_deviation", "rng"):
        assert p[key] == checkpoint_before[key]
    assert p["counters"] == {
        **checkpoint_before["counters"],
        "evaluations": checkpoint_before["counters"]["evaluations"] + 1,
    }
    diagnostic = p["diagnostic"]
    for flag in ("completed", "startup_verified", "final_verified", "final_attempted"):
        assert diagnostic[flag] is False
    assert diagnostic["attempted_step"] == cp.absolute_step
    assert diagnostic["failure_stage"] == "startup"
    assert diagnostic["failed_trial_evaluations"] == 0
    assert diagnostic["message"]
    if abs(offset) < 1e-8:
        assert "guard" in diagnostic["message"]
    else:
        assert "Independent final energy differs" in diagnostic["message"]
    assert p["final_evaluation"] is None
    with pytest.raises(InvalidDynamicsCheckpointError):
        create_dynamics_checkpoint(failed)
    for allow in (False, True):
        with pytest.raises(DynamicsInputError):
            failed.to_system(system, allow_incomplete=allow)
    assert cp.payload == checkpoint_before
    cp.validate_integrity()
    assert system.to_dict() == system_before
    rng_after = np.random.get_state()
    assert rng_before[0] == rng_after[0]
    np.testing.assert_array_equal(rng_before[1], rng_after[1])
    assert rng_before[2:] == rng_after[2:]
    resumed = resume_dynamics(cp, system, Harmonic(0), budgets(2))
    resumed.validate_integrity()
    assert resumed.termination_reason == "completed"
    assert resumed.completed_steps == 2
