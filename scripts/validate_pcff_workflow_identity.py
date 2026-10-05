"""J1 identity correction: bounded workflow using a hash-checked retained bundle.

This is a new integration experiment, not independent force-model verification.
No construction, preparation, typing or QM occurs; original evidence stays intact.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from validate_oplsaa_dynamics import compare_frames
from validate_prepared_workflow import count, hashes, require

from island.dynamics import run_dynamics_segment
from island.evaluation.pcff_identity import pcff_evaluation_identity
from island.forcefields import (
    PreparedForceFieldSources,
    create_evaluator,
    load_prepared_forcefield,
)
from island.minimization import MinimizationOptions
from island.workflows import (
    PreparedWorkflowConfig,
    prepared_workflow_status,
    read_prepared_workflow_frames,
    resume_prepared_workflow,
    start_prepared_bundle_workflow,
    storage,
)
from island.workflows.chain import _load_segment
from island.workflows.prepared import _load_setup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j1.json")
    )
    parser.add_argument("--child", action="store_true")
    args = parser.parse_args()
    root = args.output.resolve()
    sources = PreparedForceFieldSources(pcff_frc=args.source.resolve())
    if args.child:
        import island.workflows.prepared as workflow

        def forbidden(*a, **k):
            raise AssertionError("setup during resume")

        with (
            patch.object(workflow, "minimize_geometry", forbidden),
            patch.object(workflow, "initialize_velocities", forbidden),
        ):
            result, contexts = count(
                lambda: resume_prepared_workflow(root / "relocated", sources=sources)
            )
        storage.publish(
            root / "child.json",
            storage.json_bytes(
                {"pid": os.getpid(), "status": result["status"], "contexts": contexts}
            ),
        )
        return 0 if result["status"] == "completed" else 1
    root.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    outcome = {
        "passed": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    before = None
    try:
        before = hashes(args.bundle)
        expected = {
            k.removeprefix("relocated/"): v
            for k, v in storage.read_json(args.evidence)["bundle_cli"][
                "artifacts"
            ].items()
            if k.startswith("relocated/")
        }
        require(before == expected, "Retained J1 bundle differs from committed hashes")
        loaded = load_prepared_forcefield(args.bundle, sources=sources)
        settings, parameter, model = pcff_evaluation_identity(
            loaded.prepared.native_result, system=loaded.system
        )
        require(
            settings["implementation"] == "island_pcff_source_singlepoint_v1",
            "Expanded model required",
        )
        config = PreparedWorkflowConfig(
            total_steps=4,
            segment_steps=2,
            recording_interval=1,
            max_evaluations_per_segment=4,
            max_frames_per_segment=3,
            minimization=MinimizationOptions(
                force_tolerance=0.1, max_iterations=5000, max_evaluations=10000
            ),
        )
        declaration = {
            "schema": "island_j1_workflow_identity_correction_v1",
            "experiment": "new bounded workflow integration from retained cyclohexane bundle; no independent force-model oracle",
            "construction": "retained J1 C1CCCCC1, seed 2026; original coordinates used unchanged",
            "input_hashes": before,
            "prior_evidence_sha256": storage.checksum(args.evidence.read_bytes()),
            "source_sha256": storage.checksum(args.source.read_bytes()),
            "prepared_identity": loaded.prepared.identity,
            "native_identity": loaded.prepared.native_result.identity,
            "parameter_fingerprint": parameter,
            "model_fingerprint": model,
            "evaluation_settings": settings,
            "config": config.to_dict(),
            "split_comparison": {
                "atol": 1e-10,
                "rtol": 1e-12,
                "units": "angstrom, angstrom/ps, kJ/mol, kJ/(mol*angstrom), ps",
            },
            "whole_budget": {"steps": 4, "evaluations": 6, "frames": 5},
            "command": sys.argv,
            "python": sys.version,
            "parent_pid": os.getpid(),
        }
        storage.publish(root / "declaration.json", storage.json_bytes(declaration))
        manifest, contexts = count(
            lambda: start_prepared_bundle_workflow(
                args.bundle, root / "run", config, sources=sources
            )
        )
        outcome.update(start_status=manifest["status"], start_contexts=contexts)
        require(
            manifest["status"] == "paused" and manifest["accepted_step"] == 2,
            "Minimum/first segment gate failed",
        )
        system, prepared, initialization, minimum = _load_setup(
            root / "run", manifest, sources, with_minimum=True
        )
        outcome["minimum"] = {
            "initial_energy": minimum.initial_energy,
            "final_energy": minimum.final_energy,
            "fmax": minimum.final_fmax,
            "rms": minimum.final_rms_force,
            "iterations": minimum.iterations,
            "evaluations": minimum.evaluations,
            "reason": minimum.termination_reason,
            "verified": minimum.final_evaluation_verified,
        }
        require(
            minimum.final_evaluation.model_fingerprint == model,
            "Minimum identity mismatch",
        )
        evaluator = create_evaluator(system, prepared)

        def whole_run():
            with evaluator.open_session() as session:
                return run_dynamics_segment(
                    system,
                    session,
                    initialization.velocities,
                    replace(config.langevin(4), max_evaluations=6, max_frames=5),
                )

        whole, whole_contexts = count(whole_run)
        storage.publish(
            root / "whole.json",
            storage.json_bytes(
                {"payload": whole.payload, "sha256": whole.content_checksum}
            ),
        )
        require(whole.completed, "Uninterrupted diagnostic comparison failed")
        shutil.move(root / "run", root / "relocated")
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--child",
            "--source",
            str(args.source.resolve()),
            "--output",
            str(root),
        ]
        child = subprocess.run(
            command, capture_output=True, text=True, check=False, timeout=600
        )
        storage.publish(root / "child.stdout", child.stdout.encode())
        storage.publish(root / "child.stderr", child.stderr.encode())
        require(child.returncode == 0, "Separate-process continuation failed")
        completed = prepared_workflow_status(root / "relocated", sources=sources)
        frames = read_prepared_workflow_frames(root / "relocated", sources=sources)
        final = _load_segment(root / "relocated", completed["segments"][-1])
        outcome["differences"] = compare_frames(whole.frames, frames)
        require(whole.payload["rng"] == final.payload["rng"], "RNG mismatch")
        require(
            [f.step for f in frames] == list(range(5)), "Retained schedule mismatch"
        )
        require(
            all(f.evaluation.model_fingerprint == model for f in frames),
            "Frame model mismatch",
        )
        require(
            completed["prepared_identity"] == loaded.prepared.identity,
            "Relocated facade mismatch",
        )
        require(
            resume_prepared_workflow(root / "relocated", sources=sources) == completed,
            "Completed resume mismatch",
        )
        outcome.update(
            passed=True,
            steps=4,
            time_ps=frames[-1].time_ps,
            frames=len(frames),
            whole_contexts=whole_contexts,
            whole_evaluations=whole.evaluations,
            split_evaluations=final.payload["counters"]["evaluations"],
            rng_equal=True,
            parameter_fingerprint=parameter,
            model_fingerprint=model,
            child_command=command,
            child=storage.read_json(root / "child.json"),
        )
    except Exception as error:  # noqa: BLE001 -- retain failed experiment
        outcome.update(error_type=type(error).__name__, error=str(error))
    if before is not None:
        outcome["inputs_unchanged"] = hashes(args.bundle) == before
        outcome["passed"] = outcome["passed"] and outcome["inputs_unchanged"]
    outcome["wall_seconds"] = time.monotonic() - start
    storage.publish(root / "outcome.json", storage.json_bytes(outcome))
    print(outcome)
    return 0 if outcome["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
