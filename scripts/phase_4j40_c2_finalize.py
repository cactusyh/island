"""Observe real C2 histories and compare independently continued partial runs."""

import argparse
import concurrent.futures
import sys
from pathlib import Path

from island.analysis.equilibration import require
from island.analysis.equilibration_protocol import CorrectiveProtocol
from island.analysis.npt import FLAGS, _digest, _signed, _write_json
from scripts.phase_4j40_c2_analysis import (
    components,
    observations,
    read,
    replica_analysis,
    trajectory,
)
from scripts.validate_phase_4j38_c3 import command
from scripts.validate_phase_4j40 import CASES
from scripts.validate_phase_4j40_c2 import BASE


def post(root):
    row = read(root / "result.json")
    restart = read(root / "relocated/restart.json")
    argv = [
        sys.executable,
        "-m",
        "scripts.phase_4j40_c2_resume",
        "--resume",
        str(root / "relocated"),
        "--destination",
        str(root / "resumed-attempt"),
        "--restart",
        restart["identity"],
        "--parent",
        restart["parent_workflow_identity"],
        "--config",
        restart["npt_config_identity"],
    ]
    if (root / "resumed-attempt/attempt.json").exists():
        # Idempotent publication may reuse an actually executed, relocated child
        # artifact. Its restart must equal this independently generated split.
        code = read(root / "resume-attempt.command.json")["exit_code"]
    else:
        code = command(argv, root / "resume-attempt")
    attempt = read(root / "resumed-attempt/attempt.json")
    require(attempt["parent_restart"] == restart, "Stale relocated continuation parent")
    require(
        attempt["config"]["identity"] == restart["npt_config_identity"],
        "Changed continuation configuration",
    )
    _, _, frames, complete = trajectory(root)
    same = frames == attempt["frames"]
    comparison = _signed(
        dict(
            schema="island.j40-c2-partial-comparison.v1",
            frames_identity=_digest(frames),
            resumed_frames_identity=_digest(attempt["frames"]),
            byte_exact_retained_history=same,
            completed_workflow=complete,
            completed_continuation=attempt["status"] == "completed",
            fields="all frame fields including coordinates velocities box energies temperature pressure volume checkpoint bytes and complete RNG",
            barostat_state="independent NPT prefix replay then boundary checkpoint restore",
            tolerances={
                "absolute": 0,
                "relative": 0,
                "checkpoint_and_rng": "byte_exact",
            },
            **FLAGS,
        )
    )
    _write_json(root / "partial-comparison.json", comparison)
    row.pop("identity", None)
    if (root / "uninterrupted/failure.json").exists():
        row["retained_failure_identity"] = read(root / "uninterrupted/failure.json")[
            "identity"
        ]
    _write_json(root / "result.json", _signed(row))
    _write_json(root / "observations.json", observations(root))
    _write_json(root / "component-energies.json", components(root))
    row.pop("identity", None)
    row["partial_continuation_comparison"] = comparison["identity"]
    row["fresh_process_attempt_exit_code"] = code
    row["retained_history_equality"] = same
    row["continuation_status"] = (
        "completed"
        if complete and same
        else "partial_history_equal_integration_failed"
        if same
        else "failed"
    )
    _write_json(root / "result.json", _signed(row))
    return read(root / "result.json")


def finalize(root):
    roots = [root / case / f"replica-{i}" for case in CASES for i in range(3)]
    require(
        all((r / "result.json").exists() for r in roots),
        "Missing actual replica attempts",
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
        list(pool.map(post, roots))
    matrix = _signed(
        dict(
            schema="island.j40-c2-matrix.v1",
            base_commit=BASE,
            protocol=CorrectiveProtocol().payload,
            execution_policy="execute full maximum budget; integration errors terminate a replica, scientific preflight failures do not suppress attempts",
            cases={
                case: [read(root / case / f"replica-{i}/result.json") for i in range(3)]
                for case in CASES
            },
            replica_diagnostics={
                case: replica_analysis([root / case / f"replica-{i}" for i in range(3)])
                for case in CASES
            },
            scientific_acceptance=False,
            blocked_gates=[
                "pressure_unavailable",
                "cross_engine_discrepancy",
                "integration_failure_or_insufficient_sampling",
            ],
            unsupported_origins={
                "native_aromatic_oplsaa": "source_incompatible",
                "native_am1bcc": "blocked_precision",
                "unsupported_pcff": "typing_rejection",
            },
            **FLAGS,
        )
    )
    _write_json(root / "matrix.json", matrix)
    return matrix


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finalize(a.output)
    return 1  # scientific acceptance remains blocked, independently of integrity


if __name__ == "__main__":
    raise SystemExit(main())
