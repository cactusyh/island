"""Reproducible J32 Linux validation gates; no required test exclusions."""

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

PYTHON = "../island-validation/phase4e2-env/bin/python"
FOYER = "../island-validation/phase4g1-env/bin/python"
ANALYSIS = "../island-validation/phase4j32-analysis-env/bin/python"
SUITES = {
    "j32": ["tests/test_phase_4j32_condensed.py", "tests/test_phase_4j32_bundle.py"],
    "j31_c1": [
        "tests/test_phase_4j31_network.py",
        "tests/test_phase_4j31_execution.py",
        "tests/test_phase_4j31_receipt.py",
        "tests/test_phase_4j31_c1_outcomes.py",
    ],
    "j26_j30": [
        "tests/test_phase_4j26_integrity.py",
        "tests/test_phase_4j27_rg.py",
        "tests/test_phase_4j28_selection.py",
        "tests/test_phase_4j29_real_acceptance.py",
        "tests/test_phase_4j30_multichain.py",
        "tests/test_phase_4j30_c1_integrity.py",
        "tests/test_mdanalysis_layer.py",
        "tests/test_trajectory_store.py",
    ],
    "j22_j23": [
        "tests/test_phase_4j22_periodic.py",
        "tests/test_phase_4j23_periodic_md.py",
    ],
    "full_supported": [],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("gate")
    args = parser.parse_args()
    root, name = args.directory.resolve(), args.gate
    root.mkdir(parents=True, exist_ok=True)
    env = {
        **os.environ,
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONPATH": "src:.",
        "OPENMM_CPU_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "ISLAND_J32_TEST_RUN": str(root / "acceptance"),
        "PYTHONPYCACHEPREFIX": str(root / "pycache"),
    }
    checks = {
        "acceptance": [
            PYTHON,
            "scripts/validate_phase_4j32.py",
            "--output",
            str(root / "acceptance"),
        ],
        "ruff": ["ruff", "check", "--no-cache", "."],
        "compileall": [PYTHON, "-m", "compileall", "-q", "src", "tests", "scripts"],
        "pip_check": [PYTHON, "-m", "pip", "check"],
        "foyer_pip_check": [FOYER, "-m", "pip", "check"],
        "analysis_pip_check": [ANALYSIS, "-m", "pip", "check"],
        "git_diff_check": [
            "git",
            "diff",
            "d212fbaac72583e5a91f1519166c17a65e6e7218",
            "--check",
        ],
        "historical_receipt": [PYTHON, "scripts/phase_4j31_c1_receipt.py"],
        "retained_evidence": [
            PYTHON,
            "scripts/validate_phase_4j31_c1_evidence.py",
            "--output",
            str(root / "retained-evidence.json"),
        ],
        "j32_receipt": [ANALYSIS, "scripts/phase_4j32_receipt.py"],
    }
    if name in SUITES:
        command = [
            PYTHON,
            "-m",
            "pytest",
            "-p",
            "no:cacheprovider",
            "-q",
            "--tb=short",
            "--junitxml=" + str(root / (name + ".xml")),
            "--basetemp=" + str(root / (name + "-tmp")),
            *SUITES[name],
        ]
    else:
        command = checks[name]
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    start = time.monotonic()
    with (root / (name + ".log")).open("w") as log:
        p = subprocess.run(
            command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False
        )
    record = {
        "command": command,
        "cwd": str(Path.cwd()),
        "exit_code": p.returncode,
        "tested_commit": commit,
        "tested_commit_end": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "elapsed_seconds": time.monotonic() - start,
        "environment": {
            k: env[k]
            for k in (
                "PYTHONDONTWRITEBYTECODE",
                "PYTHONPATH",
                "OPENMM_CPU_THREADS",
                "OMP_NUM_THREADS",
                "ISLAND_J32_TEST_RUN",
                "PYTHONPYCACHEPREFIX",
            )
        },
    }
    (root / (name + "-command.json")).write_text(json.dumps(record, indent=2) + "\n")
    print(name, "exit", p.returncode, flush=True)
    print((root / (name + ".log")).read_text()[-4000:], flush=True)
    raise SystemExit(p.returncode)


if __name__ == "__main__":
    main()
