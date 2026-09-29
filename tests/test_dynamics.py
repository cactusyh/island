"""Independent analytical dynamics, units, ownership, and bounded failures."""

from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import DynamicsOptions, run_nve
from island.evaluation import EvaluationResult
from island.exceptions import DynamicsInputError, InvalidDynamicsResultError


def molecule(reverse=False):
    graph = Topology()
    for site, mass in [(91, 8.0), (17, 2.0)] if reverse else [(17, 2.0), (91, 8.0)]:
        graph.add_site(AtomSite(site, "C", mass, element="C", atomic_number=6))
    return MolecularSystem(
        graph,
        Coordinates({17: (1.0, 0.2, -0.3), 91: (-0.5, 0.4, 0.1)}),
        metadata={"coordinate_source": "analytical"},
    )


VELOCITIES = {17: (0.3, -0.4, 0.1), 91: (-0.2, 0.1, 0.5)}


class Harmonic:
    def __init__(self, k=2.0, fail_call=None):
        self.k, self.calls, self.fail_call = k, 0, fail_call

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        self.calls += 1
        if self.calls == self.fail_call:
            raise RuntimeError("intentional trial evaluation failure")
        assert coordinate_unit == "angstrom"
        energy = 0.5 * self.k * sum(np.dot(x, x) for x in coordinates.values())
        return EvaluationResult(
            float(energy),
            {"harmonic": float(energy)},
            {site: tuple(-self.k * np.array(x)) for site, x in coordinates.items()},
            coordinate_hash(coordinates),
            "params",
            "model",
            "evaluation",
            "analytical",
            "1",
            "numpy",
            {},
        )


def options(steps=20, dt=0.5, **kwargs):
    return DynamicsOptions(
        timestep_fs=dt,
        steps=steps,
        max_evaluations=kwargs.pop("max_evaluations", steps + 2),
        max_frames=kwargs.pop("max_frames", steps + 1),
        **kwargs,
    )


def test_force_free_exact_full_step_motion_and_units():
    system = molecule()
    result = run_nve(system, Harmonic(0), VELOCITIES, options(20, 1.0))
    assert result.completed and result.completed_steps == result.requested_steps == 20
    assert result.evaluations == 22 and len(result.frames) == 21
    expected_kinetic = 0.005 * (
        2 * (0.3**2 + 0.4**2 + 0.1**2) + 8 * (0.2**2 + 0.1**2 + 0.5**2)
    )
    for frame in result.frames:
        assert frame.time_ps == pytest.approx(frame.step * 0.001, abs=1e-15)
        assert frame.kinetic_energy == pytest.approx(expected_kinetic, abs=1e-15)
        assert frame.total_energy == frame.kinetic_energy
        for site, velocity in VELOCITIES.items():
            np.testing.assert_allclose(
                frame.coordinates[site],
                system.coordinates.get(site) + frame.time_ps * np.array(velocity),
                rtol=0,
                atol=2e-15,
            )
            assert frame.velocities[site] == velocity
    assert result.max_abs_energy_deviation == 0


def test_explicit_one_step_unequal_masses_acceleration_and_energy():
    system = molecule()
    result = run_nve(system, Harmonic(), VELOCITIES, options(1, 1.0))
    dt = 0.001
    for site, velocity in VELOCITIES.items():
        x, v, mass = (
            system.coordinates.get(site),
            np.array(velocity),
            system.topology.sites[site].mass,
        )
        # Independent closed one-step expansion for F=-2*x, a=-200*x/m.
        expected_x = x + dt * v - dt**2 * 100 * x / mass
        expected_v = v - dt * 100 * (x + expected_x) / mass
        np.testing.assert_allclose(
            result.final_state.coordinates[site], expected_x, rtol=0, atol=1e-15
        )
        np.testing.assert_allclose(
            result.final_state.velocities[site], expected_v, rtol=0, atol=1e-15
        )
    assert result.final_state.total_energy == pytest.approx(
        result.final_state.potential_energy + result.final_state.kinetic_energy,
        abs=1e-15,
    )


def analytic_state(system, elapsed):
    xyz, velocities = {}, {}
    for site, velocity in VELOCITIES.items():
        omega = np.sqrt(200 / system.topology.sites[site].mass)
        x, v = system.coordinates.get(site), np.array(velocity)
        xyz[site] = x * np.cos(omega * elapsed) + v / omega * np.sin(omega * elapsed)
        velocities[site] = v * np.cos(omega * elapsed) - omega * x * np.sin(
            omega * elapsed
        )
    return xyz, velocities


def test_harmonic_exact_solution_second_order_and_bounded_energy():
    system = molecule()
    errors = []
    for dt, steps in ((1.0, 100), (0.5, 200), (0.25, 400)):
        result = run_nve(
            system,
            Harmonic(),
            VELOCITIES,
            options(steps, dt, max_frames=11, recording_interval=steps // 10),
        )
        exact, speed = analytic_state(system, 0.1)
        errors.append(
            max(
                np.max(
                    np.abs(np.array(result.final_state.coordinates[site]) - exact[site])
                )
                for site in VELOCITIES
            )
        )
        for site in VELOCITIES:
            np.testing.assert_allclose(
                result.final_state.velocities[site], speed[site], rtol=2e-5, atol=2e-4
            )
        assert result.max_abs_energy_deviation < 3e-5
    assert errors[0] / errors[1] == pytest.approx(4, rel=0.005)
    assert errors[1] / errors[2] == pytest.approx(4, rel=0.005)


def test_time_reversal_and_insertion_order():
    system = molecule()
    first = run_nve(system, Harmonic(), VELOCITIES, options(100))
    same = run_nve(
        molecule(True),
        Harmonic(),
        dict(reversed(list(VELOCITIES.items()))),
        options(100),
    )
    assert first.dynamics_fingerprint == same.dynamics_fingerprint
    assert first.final_state.coordinates == same.final_state.coordinates
    end = first.to_system(system)
    reverse = {
        site: tuple(-np.array(v)) for site, v in first.final_state.velocities.items()
    }
    back = run_nve(end, Harmonic(), reverse, options(100))
    for site, velocity in VELOCITIES.items():
        np.testing.assert_allclose(
            back.final_state.coordinates[site],
            system.coordinates.get(site),
            rtol=0,
            atol=2e-14,
        )
        np.testing.assert_allclose(
            back.final_state.velocities[site],
            -np.array(velocity),
            rtol=0,
            atol=2e-13,
        )


@pytest.mark.parametrize(
    "budget,completed,calls", [(1, 0, 1), (2, 0, 2), (3, 1, 3), (6, 4, 6)]
)
def test_exact_evaluator_budget(budget, completed, calls):
    evaluator = Harmonic()
    result = run_nve(molecule(), evaluator, VELOCITIES, options(max_evaluations=budget))
    assert not result.completed and result.termination_reason == "maximum_evaluations"
    assert result.completed_steps == completed
    assert result.evaluations == evaluator.calls == calls
    assert result.final_evaluation_verified == (budget > 1)
    result.validate_integrity()
    result.to_system(molecule(), allow_incomplete=True).validate()


@pytest.mark.parametrize(
    "maximum,interval,completed,steps",
    [(1, 3, 0, [0]), (2, 3, 3, [0, 3]), (3, 3, 6, [0, 3, 6])],
)
def test_storage_budget(maximum, interval, completed, steps):
    evaluator = Harmonic()
    result = run_nve(
        molecule(),
        evaluator,
        VELOCITIES,
        options(max_frames=maximum, recording_interval=interval),
    )
    assert result.termination_reason == "maximum_frames"
    assert result.completed_steps == completed
    assert [frame.step for frame in result.frames] == steps
    assert evaluator.calls == completed + 2


def test_recording_interval_always_keeps_unique_final():
    result = run_nve(
        molecule(),
        Harmonic(),
        VELOCITIES,
        options(10, recording_interval=3, max_frames=5),
    )
    assert result.completed
    assert [frame.step for frame in result.frames] == [0, 3, 6, 9, 10]


def test_failure_keeps_last_complete_state_and_counts_final_check():
    good = run_nve(molecule(), Harmonic(), VELOCITIES, options(2))
    evaluator = Harmonic(fail_call=4)
    failed = run_nve(
        molecule(),
        evaluator,
        VELOCITIES,
        options(10, recording_interval=9, max_frames=3),
    )
    assert (
        failed.termination_reason == "evaluation_failed" and failed.attempted_step == 3
    )
    assert failed.completed_steps == 2 and failed.evaluations == evaluator.calls == 5
    assert failed.final_state.coordinates == good.final_state.coordinates
    assert failed.final_state.velocities == good.final_state.velocities
    assert failed.final_evaluation_verified and not failed.completed
    with pytest.raises(DynamicsInputError, match="allow_incomplete"):
        failed.to_system(molecule())
    failed.to_system(molecule(), allow_incomplete=True).validate()


def test_initial_failure_and_final_failure():
    failed = run_nve(molecule(), Harmonic(fail_call=1), VELOCITIES, options())
    failed.validate_integrity()
    assert failed.final_state.evaluation is None and failed.evaluations == 1
    with pytest.raises(DynamicsInputError):
        failed.to_system(molecule(), allow_incomplete=True)
    final = run_nve(molecule(), Harmonic(fail_call=4), VELOCITIES, options(2))
    assert final.completed_steps == 2 and not final.completed
    assert final.termination_reason == "final_evaluation_failed"
    assert final.failure_stage == "final_evaluation"
    assert final.final_evaluation_attempted and not final.final_evaluation_verified


def test_energy_guard_no_trial_commit_no_adaptation():
    result = run_nve(
        molecule(), Harmonic(), VELOCITIES, options(5, 20, max_energy_deviation=1e-12)
    )
    assert result.termination_reason == "energy_guard_exceeded"
    assert result.completed_steps == 0 and result.attempted_step == 1
    assert result.evaluations == 3
    assert result.final_state.coordinates == result.initial_state.coordinates
    assert result.final_state.velocities == VELOCITIES
    assert abs(result.failure_details["trial_energy_deviation"]) > 1e-12
    assert result.options.timestep_fs == 20


@pytest.mark.parametrize(
    "velocities",
    [
        {17: (0, 0, 0)},
        {**VELOCITIES, 999: (0, 0, 0)},
        {17: (0, 0), 91: (0, 0, 0)},
        {17: (float("nan"), 0, 0), 91: (0, 0, 0)},
        {17: (float("inf"), 0, 0), 91: (0, 0, 0)},
        None,
    ],
)
def test_bad_velocities(velocities):
    with pytest.raises(DynamicsInputError):
        run_nve(molecule(), Harmonic(), velocities, options())


@pytest.mark.parametrize("mass", [0, -1, float("nan"), float("inf")])
def test_bad_masses(mass):
    system = molecule()
    system.topology.sites[17].mass = mass
    with pytest.raises(DynamicsInputError):
        run_nve(system, Harmonic(), VELOCITIES, options())


@pytest.mark.parametrize(
    "mode", ["model", "parameter", "units", "coverage", "fingerprint", "components"]
)
def test_changed_evaluator_contract(mode):
    class Changing(Harmonic):
        def evaluate(self, coordinates=None, **kwargs):
            record = super().evaluate(coordinates, **kwargs)
            if self.calls == 2:
                if mode == "model":
                    record = replace(record, model_fingerprint="changed")
                elif mode == "parameter":
                    record = replace(record, parameter_fingerprint="changed")
                elif mode == "coverage":
                    record = replace(record, forces={17: (0.0, 0.0, 0.0)})
                elif mode == "fingerprint":
                    record = replace(record, coordinate_fingerprint="wrong")
                elif mode == "components":
                    record = replace(record, energy_components={"bad": 999.0})
                else:
                    object.__setattr__(record, "force_unit", "newtons")
            return record

    result = run_nve(molecule(), Changing(), VELOCITIES, options())
    assert result.completed_steps == 0 and result.evaluations == 3
    assert result.termination_reason == "evaluation_failed"


def test_ownership_reconstruction_and_identity():
    system = molecule()
    before = deepcopy(system.to_dict())
    velocities = dict(VELOCITIES)
    result = run_nve(system, Harmonic(), velocities, options())
    velocities[17] = (999, 0, 0)
    assert (
        system.to_dict() == before
        and result.initial_state.velocities[17] == VELOCITIES[17]
    )
    assert deepcopy(result) is result
    assert replace(result).dynamics_fingerprint == result.dynamics_fingerprint
    with pytest.raises(TypeError):
        result.final_state.coordinates[17] = (0, 0, 0)
    changed = run_nve(system, Harmonic(), {**VELOCITIES, 17: (0.2, 0, 0)}, options())
    assert changed.dynamics_fingerprint != result.dynamics_fingerprint
    result.to_system(system).validate()


@pytest.mark.parametrize(
    "changes",
    [
        {"completed_steps": 0},
        {"evaluations": 1},
        {"completed": False},
        {"final_evaluation_verified": False},
        {"masses": {17: 2.0}},
        {"dynamics_fingerprint": "wrong"},
        {"frames": []},
        {"failed_trial_evaluations": 7},
        {"failure_details": {"mutable": []}},
        {"options": {}},
    ],
)
def test_reconstructed_result_integrity(changes):
    system = molecule()
    before = deepcopy(system.to_dict())
    result = run_nve(system, Harmonic(), VELOCITIES, options(2))
    with pytest.raises(InvalidDynamicsResultError):
        bad = replace(result, **changes)
        bad.to_system(system, allow_incomplete=True)
    assert system.to_dict() == before


def test_zero_total_energy_uses_absolute_guard():
    kinetic = 0.005 * sum(
        molecule().topology.sites[site].mass * np.dot(v, v)
        for site, v in VELOCITIES.items()
    )

    class Offset(Harmonic):
        def evaluate(self, *args, **kwargs):
            record = super().evaluate(*args, **kwargs)
            return replace(
                record,
                potential_energy=-kinetic,
                energy_components={"offset": -kinetic},
            )

    result = run_nve(
        molecule(), Offset(0), VELOCITIES, options(5, max_energy_deviation=0)
    )
    assert result.completed and result.initial_state.total_energy == 0
    assert result.max_abs_energy_deviation == 0


def test_unsupported_velocity_units():
    with pytest.raises(DynamicsInputError, match="angstrom/ps"):
        run_nve(molecule(), Harmonic(), VELOCITIES, options(), velocity_unit="nm/ps")
