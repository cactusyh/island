"""BAOAB algebra and contracts using independent analytical potentials."""

from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest
from test_dynamics import VELOCITIES, Harmonic, molecule
from test_dynamics import options as nve_options

from island.dynamics import (
    LangevinOptions,
    initialize_velocities,
    run_langevin,
    run_nve,
)
from island.dynamics.thermal import GAS_CONSTANT as R
from island.exceptions import DynamicsInputError, InvalidDynamicsResultError


def options(steps=10, **kwargs):
    return LangevinOptions(
        timestep_fs=kwargs.pop("timestep_fs", 1.0),
        temperature_kelvin=kwargs.pop("temperature_kelvin", 300),
        friction_per_ps=kwargs.pop("friction_per_ps", 2),
        thermostat_seed=kwargs.pop("thermostat_seed", 314159),
        steps=steps,
        max_evaluations=kwargs.pop("max_evaluations", steps + 2),
        max_frames=kwargs.pop("max_frames", steps + 1),
        **kwargs,
    )


def test_initializer_identity_units_ownership_and_local_rng():
    system = molecule()
    before = deepcopy(system.to_dict())
    np.random.seed(817)
    expected = np.random.random(5)
    np.random.seed(817)
    result = initialize_velocities(system, temperature_kelvin=300, seed=18)
    np.testing.assert_array_equal(np.random.random(5), expected)
    result.validate_integrity()
    assert deepcopy(result) is result
    assert (
        result.velocities
        == initialize_velocities(
            molecule(True), temperature_kelvin=300, seed=18
        ).velocities
    )
    normal = np.random.Generator(np.random.PCG64(18)).standard_normal((2, 3))
    for i, site in enumerate(sorted(result.masses)):
        np.testing.assert_allclose(
            result.velocities[site],
            normal[i] * np.sqrt(100 * R * 300 / result.masses[site]),
            atol=1e-14,
        )
    assert result.instantaneous_temperature_kelvin == pytest.approx(
        2 * result.kinetic_energy / (6 * R)
    )
    assert result.instantaneous_temperature_kelvin != 300
    assert result.degrees_of_freedom == 6
    assert system.to_dict() == before
    with pytest.raises(TypeError):
        result.velocities[17] = (0, 0, 0)
    with pytest.raises(InvalidDynamicsResultError):
        replace(result, degrees_of_freedom=3).validate_integrity()
    zero = initialize_velocities(system, temperature_kelvin=0, seed=18)
    assert zero.kinetic_energy == zero.instantaneous_temperature_kelvin == 0


@pytest.mark.parametrize(
    "temperature,seed",
    [
        (-1, 0),
        (float("nan"), 1),
        (float("inf"), 2),
        (300, True),
        (300, -1),
        (300, 1.2),
        (300, 2**128),
    ],
)
def test_invalid_thermal_inputs(temperature, seed):
    with pytest.raises(DynamicsInputError):
        initialize_velocities(molecule(), temperature_kelvin=temperature, seed=seed)
    with pytest.raises(DynamicsInputError):
        options(temperature_kelvin=temperature, thermostat_seed=seed)


@pytest.mark.parametrize("mass", [0, -1, float("nan"), float("inf")])
def test_invalid_masses(mass):
    system = molecule()
    system.topology.sites[17].mass = mass
    with pytest.raises(DynamicsInputError):
        initialize_velocities(system, temperature_kelvin=300, seed=1)
    with pytest.raises(DynamicsInputError):
        run_langevin(system, Harmonic(), VELOCITIES, options())


def test_full_step_prescribed_noise_independent_algebra(monkeypatch):
    noise = np.array([[0.7, -0.2, 0.4], [-0.5, 1.2, 0.1]])
    monkeypatch.setattr(
        "island.dynamics.langevin._normal_increments", lambda rng, shape: noise.copy()
    )
    system = molecule()
    result = run_langevin(
        system, Harmonic(), VELOCITIES, options(1, timestep_fs=10, friction_per_ps=3)
    )
    dt = 0.01
    for i, site in enumerate(sorted(VELOCITIES)):
        mass = system.topology.sites[site].mass
        x = np.array(system.coordinates.get(site))
        v = np.array(VELOCITIES[site])
        c = np.exp(-3 * dt)
        # Closed expansion: force=-2*x; write x' and v' directly.
        w = (
            c * (v - 100 * dt * x / mass)
            + np.sqrt((1 - c * c) * 100 * R * 300 / mass) * noise[i]
        )
        xp = x + 0.5 * dt * (v - 100 * dt * x / mass + w)
        vp = w - 100 * dt * xp / mass
        np.testing.assert_allclose(result.final_state.coordinates[site], xp, atol=1e-14)
        np.testing.assert_allclose(result.final_state.velocities[site], vp, atol=1e-13)
    assert result.normal_draws == 6 and result.random_steps == 1
    assert result.evaluations == 3 and result.final_evaluation_verified
    assert result.final_state.instantaneous_temperature_kelvin == pytest.approx(
        2 * result.final_state.kinetic_energy / (6 * R)
    )


def test_zero_friction_is_nve_and_zero_temperature_damps():
    system = molecule()
    a = run_langevin(system, Harmonic(), VELOCITIES, options(40, friction_per_ps=0))
    b = run_nve(system, Harmonic(), VELOCITIES, nve_options(40, 1.0))
    for fa, fb in zip(a.frames, b.frames, strict=True):
        for site in VELOCITIES:
            np.testing.assert_allclose(
                fa.coordinates[site], fb.coordinates[site], atol=3e-15, rtol=1e-14
            )
            np.testing.assert_allclose(
                fa.velocities[site], fb.velocities[site], atol=3e-15, rtol=1e-14
            )
    z = run_langevin(system, Harmonic(0), VELOCITIES, options(10, temperature_kelvin=0))
    for site in VELOCITIES:
        np.testing.assert_allclose(
            z.final_state.velocities[site],
            np.exp(-2 * 0.01) * np.array(VELOCITIES[site]),
            atol=1e-15,
        )


def test_small_gamma_dt_stable_noise(monkeypatch):
    monkeypatch.setattr(
        "island.dynamics.langevin._normal_increments", lambda rng, shape: np.ones(shape)
    )
    speeds = {site: (0.0, 0.0, 0.0) for site in VELOCITIES}
    result = run_langevin(
        molecule(), Harmonic(0), speeds, options(1, friction_per_ps=1e-20)
    )
    for site, mass in result.masses.items():
        np.testing.assert_allclose(
            result.final_state.velocities[site],
            np.full(3, np.sqrt(2e-23 * 100 * R * 300 / mass)),
            rtol=1e-14,
            atol=0,
        )


def test_reproducible_ordering_and_rng_separation():
    a = run_langevin(molecule(), Harmonic(), VELOCITIES, options())
    initialize_velocities(molecule(), temperature_kelvin=500, seed=123)
    b = run_langevin(
        molecule(True), Harmonic(), dict(reversed(list(VELOCITIES.items()))), options()
    )
    assert a.dynamics_fingerprint == b.dynamics_fingerprint
    assert a.frames == b.frames
    assert (
        a.frames
        != run_langevin(
            molecule(), Harmonic(), VELOCITIES, options(thermostat_seed=44)
        ).frames
    )


@pytest.mark.parametrize(
    "budget,steps,draws", [(1, 0, 0), (2, 0, 0), (3, 1, 6), (5, 3, 18)]
)
def test_exact_evaluation_budget_before_noise(budget, steps, draws):
    evaluator = Harmonic()
    result = run_langevin(
        molecule(), evaluator, VELOCITIES, options(max_evaluations=budget)
    )
    assert result.termination_reason == "maximum_evaluations"
    assert result.evaluations == evaluator.calls == budget
    assert result.completed_steps == steps and result.normal_draws == draws
    assert result.final_evaluation_verified == (budget > 1)
    result.validate_integrity()


def test_frame_budget_and_failed_trial_keep_full_state():
    result = run_langevin(molecule(), Harmonic(), VELOCITIES, options(max_frames=1))
    assert result.normal_draws == 0 and result.evaluations == 2
    assert result.termination_reason == "maximum_frames"
    failing = run_langevin(molecule(), Harmonic(fail_call=3), VELOCITIES, options())
    accepted = run_langevin(molecule(), Harmonic(), VELOCITIES, options(1))
    assert failing.final_state == accepted.final_state
    assert failing.normal_draws == 12 and failing.random_steps == 2
    assert failing.evaluations == 4 and failing.failed_trial_evaluations == 1
    assert failing.final_evaluation_verified and not failing.completed
    with pytest.raises(DynamicsInputError):
        failing.to_system(molecule())
    failing.to_system(molecule(), allow_incomplete=True).validate()
    first = run_langevin(molecule(), Harmonic(fail_call=1), VELOCITIES, options())
    first.validate_integrity()
    assert first.normal_draws == 0 and first.final_state.evaluation is None
    with pytest.raises(DynamicsInputError):
        first.to_system(molecule(), allow_incomplete=True)


@pytest.mark.parametrize(
    "field,value",
    [
        ("evaluations", 0),
        ("normal_draws", 999),
        ("random_steps", -1),
        ("rng_algorithm", "unknown"),
        ("completed", False),
        ("final_evaluation_verified", False),
        ("dynamics_fingerprint", "wrong"),
        ("numpy_version", None),
    ],
)
def test_result_reconstruction_rejects(field, value):
    result = run_langevin(molecule(), Harmonic(), VELOCITIES, options(2))
    if field == "rng_algorithm":
        object.__setattr__(result, field, value)
    else:
        result = replace(result, **{field: value})
    with pytest.raises(InvalidDynamicsResultError):
        result.validate_integrity()


def test_frame_reconstruction_ownership_provenance():
    system = molecule()
    before = deepcopy(system.to_dict())
    result = run_langevin(system, Harmonic(), VELOCITIES, options(2))
    assert deepcopy(result) is result
    for changes in (
        {"degrees_of_freedom": 3},
        {"instantaneous_temperature_kelvin": float("nan")},
        {"velocities": {17: (0, 0, 0)}},
    ):
        broken = replace(
            result, frames=(*result.frames[:-1], replace(result.frames[-1], **changes))
        )
        with pytest.raises(InvalidDynamicsResultError):
            broken.to_system(system)
    output = result.to_system(system)
    assert output.metadata["coordinate_source"] == "langevin_dynamics"
    assert output.metadata["dynamics"]["molecular_system_is_complete_restart"] is False
    second = run_nve(output, Harmonic(), VELOCITIES, nve_options(1, 0.1)).to_system(
        output
    )
    assert second.metadata["coordinate_source"] == "nve_dynamics"
    assert (
        second.metadata["coordinate_history"][-1]["records"]["dynamics"][
            "coordinate_source"
        ]
        == "langevin_dynamics"
    )
    assert system.to_dict() == before


@pytest.mark.parametrize(
    "velocities",
    [
        {17: (0, 0, 0)},
        {17: (0, 0, 0), 91: (0, 0, 0), 12: (0, 0, 0)},
        {17: (0, 0), 91: (0, 0, 0)},
        {17: (float("nan"), 0, 0), 91: (0, 0, 0)},
        {17: (float("inf"), 0, 0), 91: (0, 0, 0)},
        {"17": (0, 0, 0), 91: (0, 0, 0)},
    ],
)
def test_malformed_velocity_inputs(velocities):
    with pytest.raises(DynamicsInputError):
        run_langevin(molecule(), Harmonic(), velocities, options())


@pytest.mark.parametrize(
    "change", ["model", "parameter", "units", "coverage", "components"]
)
def test_trial_contract_change_rejected(change):
    class Changing(Harmonic):
        def evaluate(self, *args, **kwargs):
            result = super().evaluate(*args, **kwargs)
            if self.calls == 2:
                if change == "model":
                    return replace(result, model_fingerprint="changed")
                if change == "parameter":
                    return replace(result, parameter_fingerprint="changed")
                if change == "coverage":
                    return replace(result, forces={17: (0, 0, 0)})
                if change == "components":
                    return replace(result, energy_components={"bad": 123.0})
                object.__setattr__(result, "force_unit", "wrong")
            return result

    result = run_langevin(molecule(), Changing(), VELOCITIES, options())
    assert result.termination_reason == "evaluation_failed"
    assert result.completed_steps == 0 and result.final_evaluation_verified
    assert result.evaluations == 3 and result.normal_draws == 6


def test_invalid_geometry_and_final_failure():
    system = molecule()
    enormous = {s: (1e150, 0, 0) for s in VELOCITIES}
    before = deepcopy(system.to_dict())
    result = run_langevin(
        system, Harmonic(), enormous, options(timestep_fs=1e160, friction_per_ps=0)
    )
    assert result.termination_reason == "invalid_geometry"
    assert result.completed_steps == 0 and result.random_steps == 1
    assert result.final_state.coordinates == result.initial_state.coordinates
    assert result.final_state.velocities == result.initial_state.velocities
    assert system.to_dict() == before
    failed = run_langevin(system, Harmonic(fail_call=4), VELOCITIES, options(2))
    assert failed.termination_reason == "final_evaluation_failed"
    assert failed.completed_steps == 2 and failed.random_steps == 2
    assert not failed.completed and not failed.final_evaluation_verified


def test_no_nve_conservation_guard_and_no_exact_temperature_claim():
    result = run_langevin(
        molecule(),
        Harmonic(0),
        {s: (0, 0, 0) for s in VELOCITIES},
        options(2, timestep_fs=100, temperature_kelvin=1000, friction_per_ps=100),
    )
    assert result.completed and result.max_abs_energy_deviation > 1
    assert (
        result.final_state.instantaneous_temperature_kelvin
        != result.options.temperature_kelvin
    )


def test_large_friction_small_time_product(monkeypatch):
    monkeypatch.setattr(
        "island.dynamics.langevin._normal_increments", lambda rng, shape: np.ones(shape)
    )
    result = run_langevin(
        molecule(),
        Harmonic(0),
        {s: (0, 0, 0) for s in VELOCITIES},
        options(1, timestep_fs=1e-305, friction_per_ps=1e308),
    )
    for site, mass in result.masses.items():
        np.testing.assert_allclose(
            result.final_state.velocities[site],
            np.full(3, np.sqrt(-np.expm1(-2) * 100 * R * 300 / mass)),
            rtol=1e-14,
        )
