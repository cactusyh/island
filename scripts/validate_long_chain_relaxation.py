"""Bounded reuse of Phase 4F1 PEO20 inputs; synthetic-charge software acceptance.

No construction, parameterization, or charge generation. Declare 5000/10000
minimization budgets, fmax .1, then 20 BAOAB steps in two relocated processes.
Continuation comparisons use atol=1e-10, rtol=1e-12, fixed before execution.
"""

import argparse
import shutil
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from island.analysis import analyze_workflow, export_analysis
from island.dynamics import load_dynamics_checkpoint, run_dynamics_segment
from island.evaluation import OpenMMSinglePointEvaluator
from island.minimization import MinimizationOptions
from island.workflows import (
    WorkflowConfig,
    read_workflow_frames,
    start_prepared_workflow,
    storage,
    workflow_status,
)
from island.workflows.bundle import preparation_from, system_from
from island.workflows.chain import _load_bundle

ATOL, RTOL = 1e-10, 1e-12


def inventory(root):
    return {
        str(p.relative_to(root)): storage.checksum(p.read_bytes())
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def load_retained(root):
    """Verify original signed preparation, supplied charges, and source artifacts."""
    case = root / "peo20"
    system = system_from(storage.decode(storage.read_json(case / "input-system.json")))
    payload = storage.decode(storage.read_json(case / "preparation.json"))
    r = payload["record"]
    artifacts = storage.child(case, Path(r["artifact_dir"]).name)
    for name, expected in {
        **r["artifact_sha256"],
        "input.mol2": r["input_mol2_sha256"],
        "lineage.json": r["input_lineage_sha256"],
    }.items():
        if storage.checksum(storage.child(artifacts, name).read_bytes()) != expected:
            raise ValueError(f"Retained artifact checksum mismatch: {name}")
    prepared = preparation_from(system, payload, artifacts / "result.prmtop")
    supplied = storage.read_json(case / "charges.json")
    historical = r["provided_charge_input"]
    if (
        supplied["charges"] != dict(historical["charges"])
        or supplied["source"] != historical["source"]
    ):
        raise ValueError("Saved charges/source contradict signed preparation")
    row = next(
        x
        for x in storage.read_json(root / "report.json")["cases"]
        if x["case"] == "peo20"
    )
    if (
        row["record_signature"] != prepared.record_signature
        or row["import_signature"] != prepared.imported_result.result_signature
    ):
        raise ValueError("Phase 4F1 report signatures disagree")
    if (
        system.number_of_sites != 142
        or r["requested_force_field"] != "gaff2"
        or r["charge_method"] != "provided"
    ):
        raise ValueError("Expected retained 142-site GAFF2 provided-charge case")
    return system, prepared, artifacts, supplied


def compare_states(a, b):
    errors = {}
    ids = sorted(a.coordinates)
    for name, left, right in (
        (
            "coordinates",
            [a.coordinates[s] for s in ids],
            [b.coordinates[s] for s in ids],
        ),
        ("velocities", [a.velocities[s] for s in ids], [b.velocities[s] for s in ids]),
        (
            "forces",
            [a.evaluation.forces[s] for s in ids],
            [b.evaluation.forces[s] for s in ids],
        ),
        (
            "potential_energy",
            a.evaluation.potential_energy,
            b.evaluation.potential_energy,
        ),
        ("kinetic_energy", a.kinetic_energy, b.kinetic_energy),
        ("total_energy", a.total_energy, b.total_energy),
    ):
        np.testing.assert_allclose(left, right, atol=ATOL, rtol=RTOL)
        errors[name] = float(np.max(np.abs(np.asarray(left) - np.asarray(right))))
    for key in a.evaluation.energy_components:
        np.testing.assert_allclose(
            a.evaluation.energy_components[key],
            b.evaluation.energy_components[key],
            atol=ATOL,
            rtol=RTOL,
        )
    return errors


def execute(source, output, report):
    system, preparation, artifacts, supplied = load_retained(source)
    config = WorkflowConfig(
        "[*]CCO[*]",
        20,
        str(output / "workflow"),
        "gaff2",
        "provided",
        provided_charges={int(k): v for k, v in supplied["charges"].items()},
        charge_source=supplied["source"],
        max_atoms=1000,
        charge_tolerance=preparation.record["charge_validation_tolerance_e"],
        minimization=MinimizationOptions(
            force_tolerance=0.1,
            max_iterations=5000,
            max_evaluations=10000,
            max_line_search_steps=20,
            energy_change_tolerance=1e-12,
            energy_increase_tolerance=1e-8,
        ),
        temperature_kelvin=300,
        friction_per_ps=5,
        timestep_fs=0.1,
        total_steps=20,
        segment_steps=10,
        recording_interval=5,
        velocity_seed=78123,
        thermostat_seed=99181,
        max_evaluations_per_segment=12,
        max_frames_per_segment=3,
    )
    storage.publish(
        output / "declared-settings.json",
        storage.json_bytes(
            {
                "config": config.to_dict(),
                "platform": "Reference",
                "atol": ATOL,
                "rtol": RTOL,
                "expected_retained_steps": [0, 5, 10, 15, 20],
                "uninterrupted_max_evaluations": 22,
                "uninterrupted_max_frames": 5,
                "charge_evidence": "synthetic_software_test",
            }
        ),
    )
    report.update(
        preparation_signature=preparation.record_signature,
        imported_signature=preparation.imported_result.result_signature,
        charge_source=supplied["source"],
        artifact_checksums=dict(preparation.record["artifact_sha256"]),
        settings=config.to_dict(),
    )
    print(
        "Declared 5000 iterations/10000 calls, fmax <=0.1; starting retained workflow",
        flush=True,
    )
    t = time.perf_counter()
    manifest = start_prepared_workflow(
        config,
        system,
        preparation,
        artifacts,
        evidence="synthetic_software_test",
        segments=1,
    )
    report["setup_and_first_segment_wall_seconds"] = time.perf_counter() - t
    report["stages"] = manifest["stages"]
    stage = next((s for s in manifest["stages"] if s["stage"] == "minimization"), None)
    if stage and "record" in stage:
        p = storage.decode(storage.read_json(output / "workflow" / stage["record"]))

        def metrics(e):
            f = np.array(list(e["forces"].values()))
            return {
                "energy": e["potential_energy"],
                "fmax": float(np.max(np.linalg.norm(f, axis=1))),
                "rms_force": float(np.sqrt(np.mean(np.sum(f * f, axis=1)))),
            }

        report["minimization"] = {
            k: p[k]
            for k in (
                "converged",
                "termination_reason",
                "iterations",
                "evaluations",
                "final_evaluation_verified",
                "optimizer_message",
                "options",
            )
        }
        report["minimization"].update(
            initial=metrics(p["initial_evaluation"]),
            final=metrics(p["final_evaluation"]),
        )
    if manifest["status"] != "paused":
        report.update(acceptance="unmet", workflow_status=manifest["status"])
        return
    moved = output / "workflow-relocated"
    shutil.move(output / "workflow", moved)
    t = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, "-m", "island.workflows", "resume", str(moved)],
        capture_output=True,
        text=True,
        check=False,
    )
    storage.publish(output / "resume.stdout", completed.stdout.encode())
    storage.publish(output / "resume.stderr", completed.stderr.encode())
    report["resume_wall_seconds"] = time.perf_counter() - t
    if completed.returncode:
        raise RuntimeError(f"Separate-process resume failed ({completed.returncode})")
    manifest = workflow_status(moved)
    if manifest["status"] != "completed":
        report.update(
            acceptance="unmet",
            workflow_status=manifest["status"],
            stages=manifest["stages"],
        )
        return
    starting, prepared, initialized = _load_bundle(moved, manifest)
    with OpenMMSinglePointEvaluator(
        starting, prepared.imported_result, platform="Reference"
    ).open_session() as session:
        whole = run_dynamics_segment(
            starting,
            session,
            initialized.velocities,
            replace(config.langevin(20), max_evaluations=22, max_frames=5),
        )
    storage.publish(output / "uninterrupted.json", storage.json_bytes(whole.payload))
    if not whole.completed:
        raise RuntimeError("Uninterrupted comparison did not complete")
    cp = load_dynamics_checkpoint(moved / manifest["checkpoint"])
    frames = read_workflow_frames(moved)
    if [f.step for f in frames] != [0, 5, 10, 15, 20]:
        raise ValueError("Unexpected retained frame schedule")
    errors = {
        str(a.step): compare_states(a, b)
        for a, b in zip(frames, whole.frames, strict=True)
    }
    compare_states(cp.state, whole.final_state)
    assert cp.absolute_step == 20 and cp.time_ps == whole.final_state.time_ps
    assert cp.rng_state == whole.payload["rng"]["state"]
    assert cp.payload["origin"] == whole.payload["origin"]
    assert (
        cp.payload["counters"]["normal_draws"]
        == whole.payload["counters"]["normal_draws"]
    )
    before = inventory(moved)
    analysis = analyze_workflow(moved)
    analysis.validate_integrity()
    export_analysis(analysis, output / "analysis")
    assert before == inventory(moved)
    report.update(
        acceptance="passed",
        workflow_status=manifest["status"],
        split_calls=cp.payload["counters"]["evaluations"],
        whole_calls=whole.evaluations,
        split_counters=cp.payload["counters"],
        whole_counters=whole.payload["counters"],
        max_errors_by_step=errors,
        rng_state_equal=True,
        origin_equal=True,
        accepted_steps=cp.absolute_step,
        time_ps=cp.time_ps,
        analyzed_samples=analysis.payload["sample_count"],
        analysis_source_unchanged=True,
        backend_version=cp.state.evaluation.backend_version,
        platform=cp.state.evaluation.platform,
    )


def run(source, output):
    source, output = source.resolve(), output.resolve()
    if output == source or source in output.parents or output in source.parents:
        raise ValueError("Use a new output directory separate from retained inputs")
    output.mkdir(parents=True, exist_ok=False)
    before = inventory(source)
    report = {
        "acceptance": "unmet",
        "source": str(source),
        "source_checksums": before,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        execute(source, output, report)
    except Exception as error:  # noqa: BLE001 -- durable failed acceptance diagnostics
        report["failure"] = f"{type(error).__name__}: {error}"
    finally:
        report["retained_source_unchanged"] = before == inventory(source)
        if not report["retained_source_unchanged"]:
            report["acceptance"] = "unmet"
        storage.publish(output / "report.json", storage.json_bytes(report))
    return report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args(argv)
    report = run(args.source, args.output)
    print(report["acceptance"], report.get("failure", ""), flush=True)
    return 0 if report["acceptance"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
