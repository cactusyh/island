"""Reuse validated acceptance inputs; checkpoint, relocate, and resume in a child.

Example:
 python examples/oplsaa_checkpoint.py --inputs /run/pe3-baoab --output /new/demo
The inputs directory comes from scripts/validate_oplsaa_dynamics.py. It contains
external OPLS records and saved initialization, not a self-contained checkpoint.
"""

import argparse
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

# Share this repository's acceptance-only persistence helpers, not a new workflow API.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from validate_oplsaa_dynamics import (
    continuation_options,
    counted,
    digest,
    load_segment,
    options,
    reconstruct,
    save_segment,
)

from island.dynamics import (
    create_dynamics_checkpoint,
    run_dynamics_segment,
    save_dynamics_checkpoint,
)
from island.workflows import storage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    system, evaluator, initialization = reconstruct(args.inputs)
    manifest = storage.read_json(args.inputs / "inputs.json")
    args.output.mkdir(parents=True, exist_ok=False)
    original = args.output / "before-relocation"
    original.mkdir()
    for name in (*manifest["files"], "inputs.json"):
        shutil.copyfile(storage.child(args.inputs, name), original / name)
    segment, _ = counted(
        evaluator,
        lambda session: run_dynamics_segment(
            system, session, initialization.velocities, options("baoab", 100)
        ),
    )
    save_segment(segment, original / "first.json")
    if not segment.completed:
        print(segment.termination_reason)
        return 1
    save_dynamics_checkpoint(
        create_dynamics_checkpoint(segment), original / "checkpoint.json"
    )
    import os

    storage.publish(
        original / "resume-request.json",
        storage.json_bytes(
            {
                "inputs_sha256": digest(original / "inputs.json"),
                "checkpoint_sha256": digest(original / "checkpoint.json"),
                "parent_pid": os.getpid(),
                "options": asdict(continuation_options()),
            }
        ),
    )
    relocated = args.output / "relocated"
    shutil.move(str(original), str(relocated))
    command = [
        sys.executable,
        str(
            Path(__file__).resolve().parents[1] / "scripts/validate_oplsaa_dynamics.py"
        ),
        "--resume-directory",
        str(relocated.resolve()),
    ]
    subprocess.run(command, check=True)
    final = load_segment(relocated / "continued.json")
    print(
        "Completed:",
        final.completed,
        "absolute step:",
        final.final_state.step,
        "time (ps):",
        final.final_state.time_ps,
        "cumulative evaluations:",
        final.payload["counters"]["evaluations"],
    )
    print("Checkpoint:", relocated / "continued-checkpoint.json")
    print(
        "Not equilibration; production_validated=False; simulation_readiness=not_established"
    )
    return 0 if final.completed else 1


if __name__ == "__main__":
    raise SystemExit(main())
