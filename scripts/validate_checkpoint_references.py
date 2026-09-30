"""Opt-in unchanged archived parameters: NVE/BAOAB fresh-session continuation."""

import argparse
import hashlib
import json
import time
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from validate_dynamics_references import MINIMIZATION
from validate_evaluation_sessions import instrumentation
from validate_singlepoint_references import load_archived_preparation

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
from island.evaluation import OpenMMSinglePointEvaluator
from island.minimization import minimize_geometry

# Declared before outcomes: 60 steps, 3x20 segments, 0.1 fs, both integrators.
ATOL = 1e-10
RTOL = 1e-12
CASES = {"pe_dp3_gaff2_provided", "halomethane_gaff2_provided"}


def compare(a, b):
    maximum = {"coordinates": 0.0, "velocities": 0.0, "forces": 0.0, "energy": 0.0}
    for name in ("coordinates", "velocities", "forces"):
        left = a.evaluation.forces if name == "forces" else getattr(a, name)
        right = b.evaluation.forces if name == "forces" else getattr(b, name)
        for site in left:
            np.testing.assert_allclose(left[site], right[site], atol=ATOL, rtol=RTOL)
            maximum[name] = max(
                maximum[name], float(np.max(np.abs(np.array(left[site]) - right[site])))
            )
    for x, y in (
        (a.total_energy, b.total_energy),
        (a.kinetic_energy, b.kinetic_energy),
        (a.potential_energy, b.potential_energy),
    ):
        np.testing.assert_allclose(x, y, atol=ATOL, rtol=RTOL)
        maximum["energy"] = max(maximum["energy"], abs(x - y))
    assert a.step == b.step and a.time_ps == b.time_ps
    return maximum


def run(manifest):
    results = []
    for row in json.loads(manifest.read_text())["cases"]:
        if row["case"] not in CASES:
            continue
        system, preparation, _ = load_archived_preparation(manifest.parent, row)
        imported = preparation.imported_result
        with OpenMMSinglePointEvaluator(system, imported).open_session() as session:
            minimized = minimize_geometry(system, session, MINIMIZATION)
        assert minimized.converged
        starting = minimized.to_system(system)
        for kind in ("nve", "langevin"):
            velocities = (
                initialize_velocities(
                    starting, temperature_kelvin=300, seed=78123
                ).velocities
                if kind == "langevin"
                else {
                    s: tuple(0.2 * np.sin(s + 1 + 0.7 * np.arange(3)))
                    for s in sorted(starting.topology.sites)
                }
            )
            opts = (
                LangevinOptions(0.1, 300, 5, 99181, 60, 62, 7, 10)
                if kind == "langevin"
                else DynamicsOptions(0.1, 60, 62, 7, 10, 0.01)
            )
            begin = time.perf_counter()
            with instrumentation() as counts:
                with OpenMMSinglePointEvaluator(
                    starting, imported
                ).open_session() as session:
                    whole = run_dynamics_segment(starting, session, velocities, opts)
                with OpenMMSinglePointEvaluator(
                    starting, imported
                ).open_session() as session:
                    part = run_dynamics_segment(
                        starting,
                        session,
                        velocities,
                        replace(opts, steps=20, max_evaluations=22, max_frames=3),
                    )
                differences = []
                with TemporaryDirectory() as directory:
                    for index in range(3):
                        assert part.completed
                        for frame in part.frames:
                            matching = next(
                                f for f in whole.frames if f.step == frame.step
                            )
                            differences.append(compare(frame, matching))
                        path = Path(directory) / f"{index}.json"
                        save_dynamics_checkpoint(create_dynamics_checkpoint(part), path)
                        cp = load_dynamics_checkpoint(path)
                        if index < 2:
                            moved = part.to_system(starting)
                            with OpenMMSinglePointEvaluator(
                                moved, imported
                            ).open_session() as session:
                                part = resume_dynamics(
                                    cp,
                                    moved,
                                    session,
                                    DynamicsSegmentOptions(20, 22, 3, 10),
                                )
                assert cp.payload["rng"] == whole.payload["rng"]
                assert cp.payload["origin"] == whole.payload["origin"]
                assert (
                    cp.payload["max_abs_energy_deviation"]
                    == whole.payload["max_abs_energy_deviation"]
                )
            preparation.validate_integrity(system)
            report = {
                "case": row["case"],
                "integrator": kind,
                "completed": part.completed,
                "absolute_step": cp.absolute_step,
                "time_ps": cp.time_ps,
                "uninterrupted_calls": whole.evaluations,
                "split_calls": cp.payload["counters"]["evaluations"],
                "contexts_all_four_runs": counts["contexts"],
                "max_errors": {
                    k: max(d[k] for d in differences) for k in differences[0]
                },
                "rng_state_equal": True,
                "runtime_seconds": time.perf_counter() - begin,
                "preparation_signature": preparation.record_signature,
                "source_sha256": imported.source_sha256,
                "checkpoint_payload": cp.payload,
                "checkpoint_checksum": cp.content_checksum,
            }
            results.append(report)
            print(row["case"], kind, report["max_errors"], flush=True)
    assert len(results) == 4
    return {
        "schema": "island_checkpoint_acceptance_v1",
        "atol": ATOL,
        "rtol": RTOL,
        "cases": results,
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite report")
    result = run(args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
