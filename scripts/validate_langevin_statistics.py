"""Predetermined independent-ensemble tests of production thermal sampling/BAOAB.

No correlated time samples: independent Cartesian replicas at the final endpoint.
2048 independent atoms, half mass 2 and half mass 8; 3072 scalar samples/group.
Six-standard-error thresholds declared here before execution. No seed search.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import LangevinOptions, initialize_velocities, run_langevin
from island.dynamics.thermal import GAS_CONSTANT as R
from island.evaluation import EvaluationResult

N = 2048
TEMPERATURE = 300.0
VELOCITY_SEED = 52117
THERMOSTAT_SEED = 67231
SIGMA_LIMIT = 6.0


def ensemble():
    graph = Topology()
    for i in range(N):
        site = 7 * i + 11
        graph.add_site(
            AtomSite(
                site, "C", 2.0 if i < N // 2 else 8.0, element="C", atomic_number=6
            )
        )
    return MolecularSystem(
        graph, Coordinates({s: (0.0, 0.0, 0.0) for s in graph.sites})
    )


class Potential:
    def __init__(self, system, harmonic):
        self.constants = {
            s: atom.mass if harmonic else 0.0
            for s, atom in system.topology.sites.items()
        }

    def evaluate(self, coordinates, *, coordinate_unit="angstrom"):
        energy = float(
            sum(0.5 * self.constants[s] * np.dot(x, x) for s, x in coordinates.items())
        )
        return EvaluationResult(
            energy,
            {"harmonic": energy},
            {
                s: tuple(-self.constants[s] * np.array(x))
                for s, x in coordinates.items()
            },
            coordinate_hash(coordinates),
            "ensemble-parameters",
            "ensemble-model",
            "analytical",
            "analytical",
            "1",
            "numpy",
            {},
        )


def moments(values, mean, variance):
    values = np.asarray(values).ravel()
    n = len(values)
    actual_mean = float(np.mean(values))
    actual_variance = float(np.var(values, ddof=1))
    mean_z = abs(actual_mean - mean) / np.sqrt(variance / n)
    variance_z = abs(actual_variance - variance) / (variance * np.sqrt(2 / (n - 1)))
    assert mean_z < SIGMA_LIMIT and variance_z < SIGMA_LIMIT, (mean_z, variance_z)
    return {
        "samples": n,
        "mean": actual_mean,
        "expected_mean": float(mean),
        "variance": actual_variance,
        "expected_variance": float(variance),
        "mean_standard_errors": float(mean_z),
        "variance_standard_errors": float(variance_z),
    }


def discrete_covariance(mass, dt, gamma):
    # Independent 2x2 composition and discrete Lyapunov equation for k=mass.
    b = np.array([[1.0, 0.0], [-dt * 100 / 2, 1.0]])
    a = np.array([[1.0, dt / 2], [0.0, 1.0]])
    o = np.diag([1.0, np.exp(-gamma * dt)])
    transition = b @ a @ o @ a @ b
    injection = (
        b
        @ a
        @ np.array(
            [0.0, np.sqrt(-np.expm1(-2 * gamma * dt) * 100 * R * TEMPERATURE / mass)]
        )
    )
    cov = np.linalg.solve(
        np.eye(4) - np.kron(transition, transition),
        np.outer(injection, injection).ravel(),
    ).reshape(2, 2)
    np.testing.assert_allclose(
        cov,
        transition @ cov @ transition.T + np.outer(injection, injection),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        cov,
        np.diag(
            [
                R * TEMPERATURE / mass,
                100 * R * TEMPERATURE / mass * (1 - 100 * dt * dt / 4),
            ]
        ),
        rtol=1e-12,
        atol=1e-12,
    )
    assert np.linalg.norm(np.linalg.matrix_power(transition, 200)) < 1e-8
    return cov


def run():
    start = time.perf_counter()
    system = ensemble()
    ids = sorted(system.topology.sites)
    report = {
        "numpy": np.__version__,
        "atoms": N,
        "scalar_samples_per_mass": 3 * N // 2,
        "temperature_kelvin": TEMPERATURE,
        "velocity_seed": VELOCITY_SEED,
        "thermostat_seed": THERMOSTAT_SEED,
        "sigma_limit": SIGMA_LIMIT,
        "sampling": "independent Cartesian replicas; final endpoint only",
        "initialization": [],
        "ou": [],
        "harmonic": [],
    }
    initial = initialize_velocities(
        system, temperature_kelvin=TEMPERATURE, seed=VELOCITY_SEED
    )
    for j, mass in enumerate((2.0, 8.0)):
        speeds = np.array(
            [initial.velocities[s] for s in ids[j * N // 2 : (j + 1) * N // 2]]
        )
        entry = moments(speeds, 0.0, 100 * R * TEMPERATURE / mass)
        energy = 0.005 * mass * np.sum(speeds * speeds, axis=1)
        expected_mean = 1.5 * R * TEMPERATURE
        expected_variance = 1.5 * (R * TEMPERATURE) ** 2
        mean_z = abs(np.mean(energy) - expected_mean) / np.sqrt(
            expected_variance / len(energy)
        )
        # Chi-square(3): kurtosis=7, variance estimator SE ~ var*sqrt(6/n).
        variance_z = abs(np.var(energy, ddof=1) - expected_variance) / (
            expected_variance * np.sqrt(6 / len(energy))
        )
        assert mean_z < SIGMA_LIMIT and variance_z < SIGMA_LIMIT
        entry.update(
            mass=mass,
            kinetic_mean=float(np.mean(energy)),
            kinetic_variance=float(np.var(energy, ddof=1)),
            kinetic_mean_z=float(mean_z),
            kinetic_variance_z=float(variance_z),
        )
        report["initialization"].append(entry)
    opts = LangevinOptions(20.0, TEMPERATURE, 3.0, THERMOSTAT_SEED, 20, 22, 2, 20)
    result = run_langevin(
        system, Potential(system, False), {s: (1.0, 1.0, 1.0) for s in ids}, opts
    )
    assert result.completed
    final = result.final_state
    for j, mass in enumerate((2.0, 8.0)):
        values = [
            final.velocities[s] for s in ids[j * N // 2 : (j + 1) * N // 2]
        ]
        report["ou"].append(
            {
                "mass": mass,
                **moments(
                    values,
                    np.exp(-3 * 0.4),
                    100 * R * TEMPERATURE / mass * (-np.expm1(-2 * 3 * 0.4)),
                ),
            }
        )
    for dt in (0.1, 0.05):
        opts = LangevinOptions(
            dt * 1000, TEMPERATURE, 5.0, THERMOSTAT_SEED, 200, 202, 2, 200
        )
        result = run_langevin(
            system, Potential(system, True), {s: (0.0, 0.0, 0.0) for s in ids}, opts
        )
        assert result.completed
        final = result.final_state
        for j, mass in enumerate((2.0, 8.0)):
            covariance = discrete_covariance(mass, dt, 5.0)
            sites = ids[j * N // 2 : (j + 1) * N // 2]
            x = np.array([final.coordinates[s] for s in sites]).ravel()
            v = np.array([final.velocities[s] for s in sites]).ravel()
            cross = float(np.mean(x * v))
            cross_z = abs(cross) / np.sqrt(covariance[0, 0] * covariance[1, 1] / len(x))
            assert cross_z < SIGMA_LIMIT
            report["harmonic"].append(
                {
                    "mass": mass,
                    "timestep_ps": dt,
                    "burn_in_steps": 200,
                    "end_time_ps": 200 * dt,
                    "coordinates": moments(x, 0, covariance[0, 0]),
                    "velocities": moments(v, 0, covariance[1, 1]),
                    "cross_moment": cross,
                    "cross_standard_errors": float(cross_z),
                    "discrete_covariance": covariance.tolist(),
                    "kinetic_variance_fraction_of_continuous": 1 - 100 * dt * dt / 4,
                }
            )
    report.update(
        runtime_seconds=time.perf_counter() - start,
        production_validated=False,
        simulation_readiness="not_established",
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite report")
    report = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
