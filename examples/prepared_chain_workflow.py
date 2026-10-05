"""Checked I2 bundle -> minimum -> BAOAB -> relocate -> child-process resume.

No construction, parameterization or charge assignment. Completion is not
scientific equilibration. External pinned sources must be resolved explicitly.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from island.forcefields import PreparedForceFieldSources
from island.minimization import MinimizationOptions
from island.workflows import (
    PreparedWorkflowConfig,
    prepared_workflow_status,
    read_prepared_workflow_frames,
    resume_prepared_workflow,
    start_prepared_bundle_workflow,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--xml", type=Path)
    parser.add_argument("--frc", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    sources = PreparedForceFieldSources(opls_xml=args.xml, pcff_frc=args.frc)
    if args.resume:
        result = resume_prepared_workflow(args.output, sources=sources)
        print(result["status"], result["accepted_step"])
        print(
            "Retained steps:",
            [
                f.step
                for f in read_prepared_workflow_frames(args.output, sources=sources)
            ],
        )
        assert prepared_workflow_status(args.output, sources=sources) == result
        return 0 if result["status"] == "completed" else 1
    if args.bundle is None:
        parser.error("--bundle required for start")
    args.output.mkdir(parents=True, exist_ok=False)
    config = PreparedWorkflowConfig(
        minimization=MinimizationOptions(max_iterations=5000, max_evaluations=10000)
    )
    result = start_prepared_bundle_workflow(
        args.bundle, args.output / "before", config, sources=sources
    )
    print(result["status"], result["accepted_step"])
    if result["status"] != "paused":
        return 1
    relocated = args.output / "relocated"
    shutil.move(args.output / "before", relocated)
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--output",
        str(relocated.resolve()),
        "--resume",
    ]
    for key in ("xml", "frc"):
        value = getattr(args, key)
        if value is not None:
            command.extend(["--" + key, str(value.resolve())])
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
