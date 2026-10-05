"""Bounded I3 retained-bundle acceptance; no preparation, typing or QM."""

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

from validate_oplsaa_dynamics import compare_frames

from island.dynamics import run_dynamics_segment
from island.forcefields import PreparedForceFieldSources, create_evaluator
from island.minimization import MinimizationOptions
from island.minimization.models import force_metrics
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


def sources(args):
    return PreparedForceFieldSources(opls_xml=args.xml, pcff_frc=args.frc)


def hashes(root):
    return {
        str(p.relative_to(root)): storage.checksum(p.read_bytes())
        for p in root.rglob("*")
        if p.is_file()
    }


def require(ok, message):
    if not ok:
        raise ValueError(message)


def count(call):
    import openmm

    original = openmm.Context
    contexts = []

    def construct(*a, **kw):
        contexts.append(1)
        return original(*a, **kw)

    with patch.object(openmm, "Context", construct):
        value = call()
    return value, len(contexts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--xml", type=Path)
    parser.add_argument("--frc", type=Path)
    parser.add_argument("--child", action="store_true")
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4i2.json")
    )
    args = parser.parse_args()
    source = sources(args)
    root = args.output.resolve()
    if args.child:
        # Resume must not initialize, minimize, type, or parameterize.
        import island.workflows.prepared as workflow

        def forbidden(*a, **kw):
            raise AssertionError("Scientific setup during resume")

        with (
            patch.object(workflow, "initialize_velocities", forbidden),
            patch.object(workflow, "minimize_geometry", forbidden),
        ):
            m, contexts = count(
                lambda: resume_prepared_workflow(root / "relocated", sources=source)
            )
        storage.publish(
            root / "child.json",
            storage.json_bytes(
                {"pid": os.getpid(), "contexts": contexts, "status": m["status"]}
            ),
        )
        return 0 if m["status"] == "completed" else 1
    root.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    outcome = {
        "passed": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        before = hashes(args.bundle)
        envelope = storage.read_json(args.bundle / "manifest.json")
        family = envelope["payload"]["family"]
        evidence = storage.read_json(args.evidence)["cases"][family]
        expected = {
            name.removeprefix("relocated/"): sha
            for name, sha in evidence["files"].items()
            if name.startswith("relocated/")
        }
        require(
            before == expected,
            "Retained I2 file identities differ from committed evidence",
        )
        config = PreparedWorkflowConfig(
            minimization=MinimizationOptions(max_iterations=5000, max_evaluations=10000)
        )
        declaration = {
            "schema": "island_prepared_workflow_acceptance_v1",
            "family": family,
            "config": config.to_dict(),
            "input_hashes": before,
            "prepared": envelope["payload"]["prepared"],
            "prepared_identity": envelope["payload"]["prepared_identity"],
            "split_atol": 1e-10,
            "split_rtol": 1e-12,
            "platform": "Reference",
            "units": {
                "coordinates": "angstrom",
                "velocities": "angstrom/ps",
                "energy": "kJ/mol",
                "forces": "kJ/(mol*angstrom)",
                "time": "ps",
            },
            "whole_budget": {"steps": 200, "evaluations": 202, "frames": 11},
            "charge_provenance": "synthetic provided +/-0.01e"
            if family in ("gaff", "gaff2")
            else "source native",
            "command": sys.argv,
            "python": sys.version,
            "parent_pid": os.getpid(),
        }
        storage.publish(root / "declaration.json", storage.json_bytes(declaration))
        m, contexts = count(
            lambda: start_prepared_bundle_workflow(
                args.bundle, root / "run", config, sources=source
            )
        )
        outcome.update(start_contexts=contexts, start_status=m["status"])
        require(
            m["status"] == "paused" and m["accepted_step"] == 100,
            "Minimum/first segment gate failed",
        )
        system, prepared, initialization, minimum = _load_setup(
            root / "run", m, source, with_minimum=True
        )
        outcome["minimum"] = {
            "initial_energy": minimum.initial_evaluation.potential_energy,
            "final_energy": minimum.final_evaluation.potential_energy,
            "fmax": force_metrics(minimum.final_evaluation)[0],
            "rms": force_metrics(minimum.final_evaluation)[1],
            "iterations": minimum.iterations,
            "evaluations": minimum.evaluations,
            "reason": minimum.termination_reason,
            "verified": minimum.final_evaluation_verified,
        }
        evaluator = create_evaluator(system, prepared)
        from dataclasses import replace

        def whole_run():
            with evaluator.open_session() as session:
                return run_dynamics_segment(
                    system,
                    session,
                    initialization.velocities,
                    replace(config.langevin(200), max_evaluations=202, max_frames=11),
                )

        whole, whole_contexts = count(whole_run)
        storage.publish(
            root / "whole.json",
            storage.json_bytes(
                {"payload": whole.payload, "sha256": whole.content_checksum}
            ),
        )
        require(whole.completed, "Uninterrupted gate failed")
        shutil.move(root / "run", root / "relocated")
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--output",
            str(root),
            "--child",
        ]
        for option, value in (("--xml", args.xml), ("--frc", args.frc)):
            if value is not None:
                command.extend([option, str(value.resolve())])
        child = subprocess.run(command, capture_output=True, text=True, check=False)
        storage.publish(root / "child.stdout", child.stdout.encode())
        storage.publish(root / "child.stderr", child.stderr.encode())
        require(child.returncode == 0, "Separate-process resume failed")
        completed = prepared_workflow_status(root / "relocated", sources=source)
        frames = read_prepared_workflow_frames(root / "relocated", sources=source)
        final = _load_segment(root / "relocated", completed["segments"][-1])
        outcome["differences"] = compare_frames(whole.frames, frames)
        require(whole.payload["rng"] == final.payload["rng"], "Complete RNG differs")
        require(
            whole.payload["origin"] == final.payload["origin"],
            "Trajectory origin differs",
        )
        require(
            [f.step for f in frames] == list(range(0, 201, 20)),
            "Retained schedule differs",
        )
        require(
            resume_prepared_workflow(root / "relocated", sources=source) == completed,
            "Completed no-op resume differs",
        )
        require(hashes(args.bundle) == before, "Input mutation")
        outcome.update(
            passed=True,
            whole_contexts=whole_contexts,
            whole_evaluations=whole.evaluations,
            split_evaluations=final.payload["counters"]["evaluations"],
            rng_equal=True,
            child=storage.read_json(root / "child.json"),
            child_command=command,
            frames=len(frames),
            steps=200,
            time_ps=frames[-1].time_ps,
            source_unchanged=True,
            model_identity=minimum.final_evaluation.model_fingerprint,
        )
    except Exception as error:  # noqa: BLE001 -- preserve failed declared experiment
        outcome.update(error_type=type(error).__name__, error=str(error))
    outcome["wall_seconds"] = time.monotonic() - started
    storage.publish(root / "outcome.json", storage.json_bytes(outcome))
    print(outcome)
    return 0 if outcome["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
