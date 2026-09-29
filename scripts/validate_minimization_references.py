"""Bounded minimization acceptance against checksum-verified Amber archives."""

import argparse
import hashlib
import json
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path

import numpy as np
from validate_singlepoint_references import (
    charge_report,
    load_archived_preparation,
    reference,
)

from island import Coordinates
from island.evaluation import OpenMMSinglePointEvaluator
from island.minimization import MinimizationOptions, minimize_geometry


def run(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    reports = []
    options = MinimizationOptions(
        force_tolerance=0.1,
        max_iterations=500,
        max_evaluations=2000,
        max_line_search_steps=20,
    )
    for row in manifest["cases"]:
        original, preparation, directory = load_archived_preparation(
            manifest_path.parent, row
        )
        evaluator = OpenMMSinglePointEvaluator.from_preparation(original, preparation)
        system = deepcopy(original)
        rng = np.random.default_rng(20260929)
        system.coordinates = Coordinates(
            {
                site: np.array(original.coordinates.get(site)) + rng.normal(0, 0.025, 3)
                for site in sorted(original.topology.sites)
            }
        )
        result = minimize_geometry(system, evaluator, options)
        checked = evaluator.evaluate(result.coordinates)
        energy, forces, _ = reference(
            directory / "result.prmtop",
            dict(preparation.imported_result.mapping),
            result.coordinates,
        )
        np.testing.assert_allclose(
            checked.potential_energy, result.final_energy, rtol=0, atol=1e-10
        )
        np.testing.assert_allclose(
            checked.potential_energy, energy, rtol=2e-10, atol=2e-7
        )
        errors = []
        for site in result.coordinates:
            np.testing.assert_allclose(
                checked.forces[site], result.forces[site], rtol=0, atol=1e-10
            )
            np.testing.assert_allclose(
                result.forces[site], forces[site], rtol=2e-9, atol=2e-6
            )
            errors.append(
                float(np.max(np.abs(np.array(result.forces[site]) - forces[site])))
            )
        reports.append(
            {
                "case": row["case"],
                "converged": result.converged,
                "termination_reason": result.termination_reason,
                "optimizer_message": result.optimizer_message,
                "iterations": result.iterations,
                "evaluations": result.evaluations,
                "initial_energy_kj_mol": result.initial_energy,
                "final_energy_kj_mol": result.final_energy,
                "initial_fmax_kj_mol_angstrom": result.initial_fmax,
                "final_fmax_kj_mol_angstrom": result.final_fmax,
                "initial_rms_force": result.initial_rms_force,
                "final_rms_force": result.final_rms_force,
                "reference_energy_error": abs(energy - result.final_energy),
                "reference_max_force_error": max(errors),
                "initial_coordinates_angstrom": dict(result.initial_coordinates),
                "returned_coordinates_angstrom": dict(result.coordinates),
                "coordinate_fingerprint": result.coordinate_fingerprint,
                "model_fingerprint": result.model_fingerprint,
                "parameter_fingerprint": result.parameter_fingerprint,
                "stereochemistry": result.stereochemistry,
                "preparation_signature": preparation.record_signature,
                "source_prmtop_sha256": preparation.imported_result.source_sha256,
                "charge_report": charge_report(original, preparation, directory),
                "optimizer_version": result.optimizer_version,
                "backend_version": checked.backend_version,
                "platform": checked.platform,
                "options": asdict(options),
            }
        )
        print(
            row["case"],
            result.termination_reason,
            result.initial_energy,
            result.final_energy,
            result.final_fmax,
            result.iterations,
            result.evaluations,
            flush=True,
        )
    return {
        "schema": "island_minimization_acceptance_v1",
        "cases": reports,
        "preparation_manifest_sha256": hashlib.sha256(
            manifest_path.read_bytes()
        ).hexdigest(),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refusing to overwrite existing report")
    report = run(args.preparation_manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
