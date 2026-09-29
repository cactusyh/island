"""Opt-in checksum-verified Amber archive NVE acceptance; no parameterization."""

import argparse
import hashlib
import json
import time
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path

import numpy as np
from dynamics_reference import TOLERANCES, compare, reference_trajectory
from validate_singlepoint_references import charge_report, load_archived_preparation

from island import Coordinates
from island.dynamics import DynamicsOptions, run_nve
from island.evaluation import OpenMMSinglePointEvaluator
from island.minimization import MinimizationOptions, minimize_geometry

# Predetermined before running any archived trajectory. Duration 0.002 ps.
OPTIONS = DynamicsOptions(
    timestep_fs=0.1,
    steps=20,
    max_evaluations=22,
    max_frames=5,
    recording_interval=5,
    max_energy_deviation=0.01,
)
MINIMIZATION = MinimizationOptions(
    force_tolerance=0.1, max_iterations=500, max_evaluations=2000
)


def run(manifest_path):
    reports = []
    for row in json.loads(manifest_path.read_text())["cases"]:
        original, preparation, directory = load_archived_preparation(
            manifest_path.parent, row
        )
        evaluator = OpenMMSinglePointEvaluator.from_preparation(original, preparation)
        system = deepcopy(original)
        rng = np.random.default_rng(20260929)
        system.coordinates = Coordinates(
            {
                site: original.coordinates.get(site) + rng.normal(0, 0.025, 3)
                for site in sorted(original.topology.sites)
            }
        )
        minimized = minimize_geometry(system, evaluator, MINIMIZATION)
        summary = {
            "case": row["case"],
            "minimization_converged": minimized.converged,
            "minimization_reason": minimized.termination_reason,
            "minimization_final_fmax": minimized.final_fmax,
            "preparation_signature": preparation.record_signature,
            "imported_signature": preparation.imported_result.result_signature,
            "source_prmtop_sha256": preparation.imported_result.source_sha256,
            "charge_report": charge_report(original, preparation, directory),
        }
        if not minimized.converged:
            reports.append(summary)
            print(
                row["case"],
                "minimization failed",
                minimized.termination_reason,
                flush=True,
            )
            continue
        starting = minimized.to_system(system)
        velocities = {
            site: tuple(0.2 * np.sin(site + 1 + 0.7 * np.arange(3)))
            for site in sorted(starting.topology.sites)
        }
        begin = time.perf_counter()
        result = run_nve(starting, evaluator, velocities, OPTIONS)
        runtime = time.perf_counter() - begin
        result.validate_integrity()
        refs, source_masses = reference_trajectory(
            directory / "result.prmtop",
            dict(preparation.imported_result.mapping),
            result.initial_state.coordinates,
            velocities,
            result.masses,
            OPTIONS.timestep_fs,
            [frame.step for frame in result.frames],
        )
        differences = compare(result, refs)
        result.to_system(starting, allow_incomplete=True).validate()
        preparation.validate_integrity(original)
        summary.update(
            {
                "completed": result.completed,
                "completed_steps": result.completed_steps,
                "time_ps": result.final_state.time_ps,
                "evaluations": result.evaluations,
                "termination_reason": result.termination_reason,
                "failure_stage": result.failure_stage,
                "attempted_step": result.attempted_step,
                "message": result.message,
                "failure_details": dict(result.failure_details),
                "stereochemistry": result.stereochemistry,
                "max_abs_energy_deviation_kj_mol": result.max_abs_energy_deviation,
                "runtime_seconds": runtime,
                "seconds_per_evaluator_call": runtime / result.evaluations,
                "independent_maximum_errors": differences,
                "initial_velocities_angstrom_ps": velocities,
                "masses_dalton": dict(result.masses),
                "original_prmtop_masses_dalton": source_masses,
                "model_fingerprint": result.model_fingerprint,
                "parameter_fingerprint": result.parameter_fingerprint,
                "dynamics_fingerprint": result.dynamics_fingerprint,
                "backend_version": result.initial_state.evaluation.backend_version,
                "platform": result.initial_state.evaluation.platform,
                "frames": [
                    {
                        "step": f.step,
                        "time_ps": f.time_ps,
                        "coordinates_angstrom": dict(f.coordinates),
                        "velocities_angstrom_ps": dict(f.velocities),
                        "potential_energy": f.potential_energy,
                        "kinetic_energy": f.kinetic_energy,
                        "total_energy": f.total_energy,
                        "energy_deviation": f.energy_deviation,
                        "fmax": f.fmax,
                        "coordinate_fingerprint": f.coordinate_fingerprint,
                        "velocity_fingerprint": f.velocity_fingerprint,
                    }
                    for f in result.frames
                ],
            }
        )
        reports.append(summary)
        print(
            row["case"],
            result.termination_reason,
            result.completed_steps,
            result.evaluations,
            result.max_abs_energy_deviation,
            "runtime",
            runtime,
            "reference",
            differences,
            flush=True,
        )
    return {
        "schema": "island_nve_acceptance_v1",
        "options": asdict(OPTIONS),
        "minimization_options": asdict(MINIMIZATION),
        "velocity_construction": "v(site,axis)=0.2*sin(site+1+0.7*axis) angstrom/ps; no COM removal or rescaling",
        "reference_tolerances": TOLERANCES,
        "preparation_manifest_sha256": hashlib.sha256(
            manifest_path.read_bytes()
        ).hexdigest(),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference_generator_sha256": hashlib.sha256(
            Path(__file__).with_name("dynamics_reference.py").read_bytes()
        ).hexdigest(),
        "cases": reports,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refusing to overwrite existing acceptance output")
    report = run(args.preparation_manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not all(row.get("completed", False) for row in report["cases"]):
        raise SystemExit("Incomplete trajectories retained in report")


if __name__ == "__main__":
    main()
