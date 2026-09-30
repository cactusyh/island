"""Run three separate Python processes, then compare with an uninterrupted run.

This analytical harmonic software example needs only NumPy/core ISLAND.
A real restart must reconstruct the same external chemical/parameter inputs.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import (
    DynamicsOptions,
    DynamicsSegmentOptions,
    LangevinOptions,
    create_dynamics_checkpoint,
    initialize_velocities,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    save_dynamics_checkpoint,
)
from island.evaluation import EvaluationResult


def external_inputs():
    graph = Topology()
    for site, mass in ((91, 8.0), (17, 2.0)):
        graph.add_site(AtomSite(site, "C", mass, element="C", atomic_number=6))
    return MolecularSystem(
        graph, Coordinates({17: (1.0, 0.2, -0.3), 91: (-0.5, 0.4, 0.1)})
    ), Harmonic()


class Harmonic:
    """Explicit test potential U=sum(x*x), force=-2*x; not a parameter assignment."""

    def evaluate(self, coordinates, *, coordinate_unit="angstrom"):
        assert coordinate_unit == "angstrom"
        energy = float(sum(np.dot(x, x) for x in coordinates.values()))
        return EvaluationResult(
            energy,
            {"harmonic": energy},
            {s: tuple(-2 * np.array(x)) for s, x in coordinates.items()},
            coordinate_hash(coordinates),
            "harmonic-k2",
            "harmonic-model",
            "point",
            "analytical",
            "1",
            "numpy",
            {},
        )


def options(kind, steps):
    return (
        LangevinOptions(0.5, 300.0, 2.0, 99181, steps, steps + 2, steps + 1)
        if kind == "langevin"
        else DynamicsOptions(0.5, steps, steps + 2, steps + 1, max_energy_deviation=1.0)
    )


def worker(directory, index, kind):
    system, evaluator = external_inputs()
    if index == 0:
        velocities = initialize_velocities(
            system, temperature_kelvin=300, seed=78123
        ).velocities
        segment = run_dynamics_segment(system, evaluator, velocities, options(kind, 5))
    else:
        checkpoint = load_dynamics_checkpoint(directory / f"{kind}-{index - 1}.json")
        segment = resume_dynamics(
            checkpoint, system, evaluator, DynamicsSegmentOptions(5, 7, 6)
        )
    save_dynamics_checkpoint(
        create_dynamics_checkpoint(segment), directory / f"{kind}-{index}.json"
    )
    print(kind, index, segment.final_state.step, segment.evaluations)


def demonstrate(directory):
    report = {}
    for kind in ("nve", "langevin"):
        for index in range(3):
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    str(directory),
                    "--index",
                    str(index),
                    "--kind",
                    kind,
                ],
                check=True,
            )
        system, evaluator = external_inputs()
        velocities = initialize_velocities(
            system, temperature_kelvin=300, seed=78123
        ).velocities
        whole = create_dynamics_checkpoint(
            run_dynamics_segment(system, evaluator, velocities, options(kind, 15))
        )
        split = load_dynamics_checkpoint(directory / f"{kind}-2.json")
        assert whole.payload["state"] == split.payload["state"]
        assert whole.payload["rng"] == split.payload["rng"]
        assert whole.payload["origin"] == split.payload["origin"]
        report[kind] = {
            "absolute_step": split.absolute_step,
            "time_ps": split.time_ps,
            "max_state_difference": 0.0,
            "rng_state_equal": True,
            "split_calls": split.payload["counters"]["evaluations"],
            "uninterrupted_calls": whole.payload["counters"]["evaluations"],
        }
    print(json.dumps(report, indent=2))
    print(
        "MolecularSystem coordinates and a reused seed are not a restart. No equilibration claim."
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--index", type=int)
    parser.add_argument("--kind", choices=["nve", "langevin"])
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.index, args.kind)
    else:
        with TemporaryDirectory() as directory:
            demonstrate(Path(directory))
