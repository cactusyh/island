"""Archive actual J31 acceptance artifacts and build a self-consistent receipt.

Run only after the implementation commit and all requested checks complete.
The output and artifact directory must be new or explicitly chosen by the caller.
"""

import argparse
import json
import shutil
import subprocess
from pathlib import Path

from island.analysis.network import FLAGS, _signed
from scripts.validate_phase_4j31_receipt import (
    BASE,
    BRANCH,
    SCHEMA,
    SUITES,
    digest_file,
    fixture_summaries,
    junit_counts,
    manifest_section,
)

IMPLEMENTATION = (
    "src/island/graph/__init__.py",
    "src/island/graph/network.py",
    "src/island/analysis/__init__.py",
    "src/island/analysis/network.py",
    "scripts/phase_4j31_network_fixture.py",
    "scripts/validate_phase_4j31_network.py",
    "scripts/validate_phase_4j31_receipt.py",
    "scripts/generate_phase_4j31_receipt.py",
    "tests/test_phase_4j31_network.py",
    "tests/test_phase_4j31_execution.py",
    "tests/test_phase_4j31_receipt.py",
    "docs/phase_4j31.md",
    "docs/paper/island_capability_matrix.json",
    "docs/paper/island_implementation_audit.md",
)


def generate(validation_directory, source_directory, artifacts, *, validated_commit):
    artifacts.mkdir(parents=True, exist_ok=False)
    acceptance = artifacts / "acceptance"
    acceptance.mkdir()
    shutil.copyfile(source_directory / "matrix.json", artifacts / "source_matrix.json")
    matrix = json.loads((source_directory / "matrix.json").read_text())
    for name, row in matrix["matrix"].items():
        original = source_directory / name.replace("/", "-")
        target = acceptance / name.replace("/", "-")
        target.mkdir()
        shutil.copyfile(original / "result.json", target / "result.json")
        if row["status"] == "complete":
            for subdirectory in ("full", "relocated", "analysis-relocated"):
                shutil.copytree(original / subdirectory, target / subdirectory)
            for directory in original.glob("island-ambertools-*"):
                shutil.copytree(directory, target / "amber-source-artifacts")
    suites = {}
    for name in sorted(SUITES):
        command_record = json.loads(
            (validation_directory / f"{name}-command.json").read_text()
        )
        if command_record["exit_code"]:
            raise ValueError(f"Cannot release failed suite: {name}")
        for suffix in ("xml", "log"):
            shutil.copyfile(
                validation_directory / f"{name}.{suffix}",
                artifacts / f"{name}.{suffix}",
            )
        suites[name] = {
            "command": command_record["command"],
            "counts": junit_counts(artifacts / f"{name}.xml"),
            "junit": str(artifacts / f"{name}.xml"),
            "log": str(artifacts / f"{name}.log"),
            "exit_code": 0,
        }
    checks = json.loads((validation_directory / "checks.json").read_text())
    for name, row in checks.items():
        target = artifacts / f"{name}.log"
        shutil.copyfile(row["log"], target)
        row["log"] = str(target)
    history = {
        str(path): digest_file(path)
        for path in sorted(Path("docs/evidence").glob("phase_4j*.json"))
        if path.name != "phase_4j31.json"
    }
    receipt = _signed(
        {
            "schema": SCHEMA,
            "branch": BRANCH,
            "base_commit": BASE,
            "j30_implementation_commit": "55de228de9247017ba3d706d2f6c7bdda63f6087",
            "validated_commit": validated_commit,
            "readiness": FLAGS,
            "fixtures": fixture_summaries(),
            "force_field_source_completeness": matrix,
            "source_matrix_path": str(artifacts / "source_matrix.json"),
            "acceptance_directory": str(acceptance),
            "validation": {"suites": suites, "checks": checks},
            "historical_receipts": history,
            "implementation_files": {
                name: digest_file(name) for name in IMPLEMENTATION
            },
            "artifact_files": {
                str(path): digest_file(path)
                for path in sorted(artifacts.rglob("*"))
                if path.is_file()
            },
            "blocked_claims": [
                "production equilibrium",
                "experimental agreement",
                "Tg",
                "mechanical properties",
                "AI prediction",
                "universal force-field coverage",
                "msi2lmp parity",
                "triclinic support",
                "LAMMPS KSPACE/Ewald/PPPM validation",
            ],
        }
    )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-directory", type=Path, required=True)
    parser.add_argument("--source-directory", type=Path, required=True)
    parser.add_argument(
        "--artifacts", type=Path, default=Path("docs/evidence/phase_4j31")
    )
    args = parser.parse_args()
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    receipt = generate(
        args.validation_directory,
        args.source_directory,
        args.artifacts,
        validated_commit=commit,
    )
    Path("docs/evidence/phase_4j31.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    )
    path = Path("docs/paper/reproducibility_manifest.json")
    # Keep every byte of the historical manifest prefix unchanged.
    original = subprocess.check_output(["git", "show", f"{BASE}:{path}"], text=True)
    closing = original.rfind("}")
    original_prefix = original[:closing].rstrip()
    path.write_text(
        original_prefix
        + ',\n  "phase_4j31": '
        + json.dumps(manifest_section(receipt), indent=2)
        + "\n}\n"
    )
    print("J31 receipt identity:", receipt["identity"])


if __name__ == "__main__":
    main()
