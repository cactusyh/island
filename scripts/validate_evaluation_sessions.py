"""Reproducible session correctness, archived trajectories and measured performance."""

import argparse
import cProfile
import hashlib
import json
import platform
import pstats
import statistics
import sys
import tempfile
import time
from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import numpy as np
from dynamics_reference import compare, reference_trajectory
from validate_dynamics_references import MINIMIZATION, OPTIONS
from validate_singlepoint_references import charge_report, load_archived_preparation

from island import Coordinates
from island.dynamics import DynamicsOptions, run_nve
from island.evaluation import OpenMMSinglePointEvaluator
from island.evaluation.openmm import _OpenMMResources
from island.minimization import minimize_geometry

# Declared before runs; reuse the Phase 4E3 independent reference tolerances.
LONG = DynamicsOptions(
    0.1, 1000, 1002, 11, recording_interval=100, max_energy_deviation=0.01
)
PAIR_ATOL, PAIR_RTOL = 1e-10, 1e-12


def same_evaluation(a, b):
    np.testing.assert_allclose(
        a.potential_energy, b.potential_energy, atol=PAIR_ATOL, rtol=PAIR_RTOL
    )
    assert set(a.energy_components) == set(b.energy_components)
    for key in a.energy_components:
        np.testing.assert_allclose(
            a.energy_components[key],
            b.energy_components[key],
            atol=PAIR_ATOL,
            rtol=PAIR_RTOL,
        )
    for site in a.forces:
        np.testing.assert_allclose(
            a.forces[site], b.forces[site], atol=PAIR_ATOL, rtol=PAIR_RTOL
        )
    assert a.model_fingerprint == b.model_fingerprint
    assert a.parameter_fingerprint == b.parameter_fingerprint
    assert a.settings == b.settings


def same_run(a, b):
    assert (a.completed, a.termination_reason, a.evaluations) == (
        b.completed,
        b.termination_reason,
        b.evaluations,
    )
    assert a.completed, a.message
    error = 0.0
    for left, right in zip(a.frames, b.frames, strict=True):
        assert left.step == right.step
        same_evaluation(left.evaluation, right.evaluation)
        for site in left.coordinates:
            for x, y in (
                (left.coordinates[site], right.coordinates[site]),
                (left.velocities[site], right.velocities[site]),
            ):
                np.testing.assert_allclose(x, y, atol=PAIR_ATOL, rtol=PAIR_RTOL)
                error = max(error, float(np.max(np.abs(np.array(x) - y))))
        np.testing.assert_allclose(
            [left.kinetic_energy, left.total_energy],
            [right.kinetic_energy, right.total_energy],
            atol=PAIR_ATOL,
            rtol=PAIR_RTOL,
        )
    return error


@contextmanager
def instrumentation():
    """Wrap real operations for timings/counts; never replace numerical output."""
    import openmm as mm

    metrics = {
        "contexts": 0,
        "resource_initialization_seconds": 0.0,
        "context_initialization_seconds": 0.0,
        "cleanup_seconds": 0.0,
        "backend_query_seconds": 0.0,
    }
    init, close = _OpenMMResources.__init__, _OpenMMResources.close
    context_init, query = (
        mm.Context.__init__,
        OpenMMSinglePointEvaluator._evaluate_context,
    )

    def wrapped_init(self, *args, **kwargs):
        start = time.perf_counter()
        try:
            return init(self, *args, **kwargs)
        finally:
            metrics["resource_initialization_seconds"] += time.perf_counter() - start

    def wrapped_close(self):
        start = time.perf_counter()
        try:
            return close(self)
        finally:
            metrics["cleanup_seconds"] += time.perf_counter() - start

    def wrapped_context(self, *args, **kwargs):
        metrics["contexts"] += 1
        start = time.perf_counter()
        try:
            return context_init(self, *args, **kwargs)
        finally:
            metrics["context_initialization_seconds"] += time.perf_counter() - start

    def wrapped_query(self, *args, **kwargs):
        start = time.perf_counter()
        try:
            return query(self, *args, **kwargs)
        finally:
            metrics["backend_query_seconds"] += time.perf_counter() - start

    with (
        patch.object(_OpenMMResources, "__init__", wrapped_init),
        patch.object(_OpenMMResources, "close", wrapped_close),
        patch.object(mm.Context, "__init__", wrapped_context),
        patch.object(OpenMMSinglePointEvaluator, "_evaluate_context", wrapped_query),
    ):
        yield metrics


def velocities(system):
    return {
        site: tuple(0.2 * np.sin(site + 1 + 0.7 * np.arange(3)))
        for site in sorted(system.topology.sites)
    }


def benchmark(system, imported, repeats=3):
    binding = []
    for _ in range(repeats):
        start = time.perf_counter()
        evaluator = OpenMMSinglePointEvaluator(system, imported)
        binding.append(time.perf_counter() - start)
    coords = {
        site: tuple(system.coordinates.get(site)) for site in system.topology.sites
    }
    speeds = velocities(system)
    report = {
        "atoms": len(coords),
        "repeats": repeats,
        "binding_seconds_median": statistics.median(binding),
        "backend_version": evaluator.evaluate().backend_version,
        "platform": "Reference",
        "precision": "double",
        "calls_per_singlepoint_batch": 50,
        "calls_per_nve_run": 22,
        "nve_steps": 20,
    }
    for workload in ("singlepoints", "nve"):
        for reuse in (False, True):

            def operation(potential, workload=workload):
                if workload == "singlepoints":
                    for _ in range(50):
                        potential.evaluate(coords)
                else:
                    result = run_nve(system, potential, speeds, OPTIONS)
                    assert result.completed, result.message

            # One complete unmeasured warmup in each mode, with identical work.
            if reuse:
                with evaluator.open_session() as session:
                    operation(session)
            else:
                operation(evaluator)
            samples = []
            for _ in range(repeats):
                with instrumentation() as metrics:
                    start = time.perf_counter()
                    if reuse:
                        opening = time.perf_counter()
                        session = evaluator.open_session()
                        metrics["session_open_seconds"] = time.perf_counter() - opening
                        try:
                            work_start = time.perf_counter()
                            operation(session)
                            metrics["work_seconds"] = time.perf_counter() - work_start
                        finally:
                            closing = time.perf_counter()
                            session.close()
                            metrics["session_close_seconds"] = (
                                time.perf_counter() - closing
                            )
                    else:
                        work_start = time.perf_counter()
                        operation(evaluator)
                        metrics["work_seconds"] = time.perf_counter() - work_start
                    metrics["total_seconds"] = time.perf_counter() - start
                expected = (
                    (1 if workload == "singlepoints" else 2)
                    if reuse
                    else (50 if workload == "singlepoints" else 22)
                )
                assert metrics["contexts"] == expected
                samples.append(dict(metrics))
            report[workload + ("_session" if reuse else "_fresh")] = {
                "samples": samples,
                "medians": {
                    key: statistics.median(row[key] for row in samples)
                    for key in samples[0]
                },
            }
    # Separate profile: do not fold profiler overhead into benchmark medians.
    profiler = cProfile.Profile()
    with evaluator.open_session() as session:
        profiler.enable()
        run_nve(system, session, speeds, OPTIONS)
        profiler.disable()
    stats = pstats.Stats(profiler)
    report["profile_top_cumulative"] = [
        {
            "function": f"{Path(key[0]).name}:{key[1]}:{key[2]}",
            "calls": value[1],
            "seconds": value[3],
        }
        for key, value in sorted(
            stats.stats.items(), key=lambda pair: pair[1][3], reverse=True
        )[:12]
    ]
    return report


def acceptance(original, preparation, directory, name):
    evaluator = OpenMMSinglePointEvaluator.from_preparation(original, preparation)
    system = deepcopy(original)
    rng = np.random.default_rng(20260929)
    system.coordinates = Coordinates(
        {
            site: original.coordinates.get(site) + rng.normal(0, 0.025, 3)
            for site in sorted(original.topology.sites)
        }
    )
    coords = {site: system.coordinates.get(site) for site in system.topology.sites}
    with evaluator.open_session() as session:
        for shift in (0.0, 0.0003, 0.0):
            frame = {site: xyz + shift * (site + 1) for site, xyz in coords.items()}
            same_evaluation(evaluator.evaluate(frame), session.evaluate(frame))
    fresh_min = minimize_geometry(system, evaluator, MINIMIZATION)
    with instrumentation() as counts, evaluator.open_session() as session:
        reused_min = minimize_geometry(system, session, MINIMIZATION)
    assert fresh_min.converged and reused_min.converged
    assert counts["contexts"] == 2
    same_evaluation(
        evaluator.evaluate(fresh_min.coordinates),
        evaluator.evaluate(reused_min.coordinates),
    )
    np.testing.assert_allclose(
        fresh_min.final_fmax, reused_min.final_fmax, atol=1e-8, rtol=1e-10
    )
    starting = reused_min.to_system(system)
    reports = []
    for options in (OPTIONS, LONG) if name == "phenol_gaff2_am1bcc" else (OPTIONS,):
        speed = velocities(starting)
        fresh = run_nve(starting, evaluator, speed, options)
        with instrumentation() as count:
            start = time.perf_counter()
            with evaluator.open_session() as session:
                reused = run_nve(starting, session, speed, options)
            elapsed = time.perf_counter() - start
        error = same_run(fresh, reused)
        assert count["contexts"] == 2
        refs, source_masses = reference_trajectory(
            directory / "result.prmtop",
            dict(preparation.imported_result.mapping),
            reused.initial_state.coordinates,
            speed,
            reused.masses,
            options.timestep_fs,
            [f.step for f in reused.frames],
        )
        independent_errors = compare(reused, refs)
        reused.to_system(starting).validate()
        reports.append(
            {
                "steps": reused.completed_steps,
                "duration_ps": reused.final_state.time_ps,
                "timestep_fs": options.timestep_fs,
                "recording_interval": options.recording_interval,
                "energy_guard_kj_mol": options.max_energy_deviation,
                "termination_reason": reused.termination_reason,
                "stereochemistry": reused.stereochemistry,
                "evaluations": reused.evaluations,
                "session_contexts": count["contexts"],
                "fresh_contexts": fresh.evaluations,
                "max_energy_deviation_kj_mol": reused.max_abs_energy_deviation,
                "elapsed_seconds": elapsed,
                "fresh_session_max_position_velocity_difference": error,
                "independent_reference_errors": independent_errors,
                "initial_coordinates": dict(reused.initial_state.coordinates),
                "initial_velocities": speed,
                "masses": dict(reused.masses),
                "source_masses": source_masses,
                "model_fingerprint": reused.model_fingerprint,
                "parameter_fingerprint": reused.parameter_fingerprint,
                "dynamics_fingerprint": reused.dynamics_fingerprint,
                "frames": [
                    {
                        "step": f.step,
                        "time_ps": f.time_ps,
                        "potential": f.potential_energy,
                        "kinetic": f.kinetic_energy,
                        "deviation": f.energy_deviation,
                        "coordinate_fingerprint": f.coordinate_fingerprint,
                        "velocity_fingerprint": f.velocity_fingerprint,
                    }
                    for f in reused.frames
                ],
            }
        )
        print(
            name,
            options.steps,
            "completed",
            reused.evaluations,
            "calls",
            count["contexts"],
            "contexts",
            "max_dE",
            reused.max_abs_energy_deviation,
            "reference",
            independent_errors,
            flush=True,
        )
    preparation.validate_integrity(original)
    return (
        {
            "case": name,
            "preparation_signature": preparation.record_signature,
            "source_prmtop_sha256": preparation.imported_result.source_sha256,
            "charge_report": charge_report(original, preparation, directory),
            "minimization": {
                "fresh_calls": fresh_min.evaluations,
                "session_calls": reused_min.evaluations,
                "session_contexts": 2,
                "fresh_final_energy": fresh_min.final_energy,
                "session_final_energy": reused_min.final_energy,
                "fresh_final_fmax": fresh_min.final_fmax,
                "session_final_fmax": reused_min.final_fmax,
            },
            "trajectories": reports,
        },
        starting,
        preparation.imported_result,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refusing to overwrite an existing report")
    report = {
        "schema": "island_evaluation_sessions_v1",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pair_atol": PAIR_ATOL,
        "pair_rtol": PAIR_RTOL,
        "benchmarks": {},
        "cases": [],
        "production_validated": False,
        "simulation_readiness": "not_established",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    # Repository synthetic fixture is an independently written prmtop with nonzero
    # charges and multi-term torsions, shared with the established numerical tests.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
    from test_singlepoint import source_case

    with tempfile.TemporaryDirectory() as directory:
        system, imported, _, _ = source_case(Path(directory))
        report["benchmarks"]["synthetic_four_site"] = benchmark(system, imported)
    if args.preparation_manifest:
        report["manifest_sha256"] = hashlib.sha256(
            args.preparation_manifest.read_bytes()
        ).hexdigest()
        for row in json.loads(args.preparation_manifest.read_text())["cases"]:
            original, prep, directory = load_archived_preparation(
                args.preparation_manifest.parent, row
            )
            record, minimized, imported = acceptance(
                original, prep, directory, row["case"]
            )
            report["cases"].append(record)
            if row["case"] == "pe_dp3_gaff2_provided":
                report["benchmarks"]["capped_PE_DP3"] = benchmark(minimized, imported)
    else:
        report["archive_skip"] = "No --preparation-manifest supplied"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    for name, bench in report["benchmarks"].items():
        for work in ("singlepoints", "nve"):
            fresh = bench[work + "_fresh"]["medians"]["total_seconds"]
            reused = bench[work + "_session"]["medians"]["total_seconds"]
            print(
                name,
                work,
                "median seconds fresh/session",
                fresh,
                reused,
                "ratio",
                fresh / reused,
            )


if __name__ == "__main__":
    main()
