"""Declared PE workflow acceptance: 60 x 0.1 fs, three processes, atol=1e-10/rtol=1e-12.

Live: explicit GAFF2/AM1-BCC tolerance 1e-4. Archive: original provided-charge
PE with its recorded charge tolerance. Never switch methods following failure.
"""

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

from island.dynamics import load_dynamics_checkpoint, run_dynamics_segment
from island.evaluation import OpenMMSinglePointEvaluator
from island.workflows import (
    WorkflowConfig,
    read_workflow_frames,
    start_prepared_workflow,
    start_workflow,
    storage,
    workflow_status,
)
from island.workflows.chain import _load_bundle


def run(output, manifest=None):
    config = WorkflowConfig.from_dict(
        json.loads(
            (
                Path(__file__).parents[1] / "examples/short_chain_workflow.json"
            ).read_text()
        )
    )
    config = replace(config, output_directory=str(output))
    if manifest is None:
        result = start_workflow(config)
    else:
        from validate_singlepoint_references import load_archived_preparation

        rows = json.loads(manifest.read_text())["cases"]
        row = next(r for r in rows if r["case"] == "pe_dp3_gaff2_provided")
        system, preparation, artifacts = load_archived_preparation(manifest.parent, row)
        config = replace(
            config,
            charge_method="provided",
            provided_charges={
                s: r.charge
                for s, r in preparation.imported_result.charge_result.assignments.items()
            },
            charge_tolerance=preparation.record["charge_validation_tolerance_e"],
        )
        result = start_prepared_workflow(
            config, system, preparation, artifacts, evidence="archived_parameters"
        )
    if result["status"] != "paused":
        return {"acceptance": "unmet", "manifest": result}
    moved = output.with_name(output.name + "-relocated")
    shutil.move(output, moved)
    for _ in range(2):
        subprocess.run(
            [sys.executable, "-m", "island.workflows", "resume", str(moved)],
            stdout=subprocess.DEVNULL,
            check=True,
        )
    result = workflow_status(moved)
    assert result["status"] == "completed"
    starting, preparation, initialized = _load_bundle(moved, result)
    with OpenMMSinglePointEvaluator(
        starting, preparation.imported_result
    ).open_session() as session:
        whole = run_dynamics_segment(
            starting,
            session,
            initialized.velocities,
            replace(config.langevin(60), max_evaluations=62, max_frames=7),
        )
    assert whole.completed
    cp = load_dynamics_checkpoint(moved / result["checkpoint"])
    errors = {}
    for name in ("coordinates", "velocities"):
        a = np.array(list(getattr(cp.state, name).values()))
        b = np.array(list(getattr(whole.final_state, name).values()))
        np.testing.assert_allclose(a, b, atol=1e-10, rtol=1e-12)
        errors[name] = float(np.max(np.abs(a - b)))
    a = np.array(list(cp.state.evaluation.forces.values()))
    b = np.array(list(whole.final_state.evaluation.forces.values()))
    np.testing.assert_allclose(a, b, atol=1e-10, rtol=1e-12)
    errors["forces"] = float(np.max(np.abs(a - b)))
    errors["energy"] = abs(cp.state.total_energy - whole.final_state.total_energy)
    np.testing.assert_allclose(
        cp.state.total_energy, whole.final_state.total_energy, atol=1e-10, rtol=1e-12
    )
    assert cp.rng_state == whole.payload["rng"]["state"]
    assert cp.payload["origin"] == whole.payload["origin"]
    assert (
        cp.payload["counters"]["normal_draws"]
        == whole.payload["counters"]["normal_draws"]
    )
    saved_bundle = storage.decode(storage.read_json(moved / result["bundle"]))
    minimum = saved_bundle["minimum"]
    from island.workflows.bundle import minimum_from

    validated_minimum = minimum_from(minimum)
    return {
        "acceptance": "passed",
        "evidence": result["evidence"],
        "directory": str(moved),
        "settings": config.to_dict(),
        "steps": cp.absolute_step,
        "time_ps": cp.time_ps,
        "max_errors": errors,
        "rng_state_equal": True,
        "split_calls": cp.payload["counters"]["evaluations"],
        "whole_calls": whole.evaluations,
        "frames": len(read_workflow_frames(moved)),
        "final_temperature_kelvin": cp.state.instantaneous_temperature_kelvin,
        "final_max_force": float(np.max(np.linalg.norm(a, axis=1))),
        "preparation_signature": preparation.record_signature,
        "imported_signature": preparation.imported_result.result_signature,
        "charge_residual_e": preparation.imported_result.charge_result.total_charge_residual,
        "minimization": {
            "initial_energy": validated_minimum.initial_energy,
            "final_energy": validated_minimum.final_energy,
            "final_fmax": validated_minimum.final_fmax,
            "evaluations": validated_minimum.evaluations,
            "final_verified": validated_minimum.final_evaluation_verified,
        },
        "backend": {
            "name": cp.state.evaluation.backend_name,
            "version": cp.state.evaluation.backend_version,
            "platform": cp.state.evaluation.platform,
        },
        "tool_package": preparation.record.get("tool_package"),
        "source_sha256": preparation.imported_result.source_sha256,
        "artifact_checksums": preparation.record["artifact_sha256"],
        "stages": result["stages"],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.output, args.manifest)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(report["acceptance"])
