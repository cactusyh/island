"""Validate C2 source decisions and diagnostic receipts, not unsupported acceptance."""

import argparse
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from scripts.phase_4j38_c1_audit import require, sha
from scripts.phase_4j38_c2 import (
    BASE,
    SourceIncompatibleError,
    aromatic_decision,
    require_native_source,
    validate_matrix,
)
from scripts.phase_4j38_contract import validate_identity
from scripts.validate_phase_4j38_c1_receipt import validate_bundle as validate_c1_bundle

SOURCES = {
    "scripts/phase_4j38_c2.py",
    "scripts/validate_phase_4j38_c2.py",
    "scripts/validate_phase_4j38_c2_receipt.py",
    "tests/test_phase_4j38_c2.py",
}
COMMANDS = {
    "ruff",
    "compileall",
    "pip_supported",
    "pip_foyer",
    "pip_analysis",
    "pip_amber",
    "pip_quantum",
    "diff_check",
    "historical_receipts",
    "historical_hashes",
    "c1_diagnostic",
    "c1_acceptance",
    "raw_driver",
    "relocated_bundle",
}


def read(path):
    return json.loads(Path(path).read_text())


def validate_bundle(root):
    root = Path(root)
    raw = root / "c1-recheck"
    validate_c1_bundle(raw)
    decision = aromatic_decision(
        read(raw / "aromatic-audit.json"), read(raw / "aromatic-fixture.json")
    )
    require(
        decision == read(root / "native-source-decision.json"),
        "Changed native source compatibility, vector, graph or provenance",
    )
    try:
        require_native_source(decision)
    except SourceIncompatibleError as exc:
        rejection = {
            "exception": type(exc).__name__,
            "diagnostic": str(exc),
            "classification": "source_incompatible",
            "dynamics_executed": False,
        }
    else:
        raise ValueError("Unexpected native source acceptance")
    require(
        rejection == read(root / "native-rejection.json"),
        "Changed native source rejection",
    )
    validate_matrix(
        read(root / "matrix.json"),
        read(raw / "matrix.json"),
        decision,
        read(raw / "resp-inventory.json"),
    )
    require(
        read(root / "run-state.json")
        == {
            "status": "diagnostic_complete",
            "acceptance_met": False,
            "dynamics_executed": False,
        },
        "False execution/acceptance claim",
    )
    return True


def historical_hashes(base=BASE):
    names = subprocess.check_output(
        [
            "git",
            "ls-tree",
            "-r",
            "--name-only",
            base,
            "docs",
            ".gitattributes",
        ],
        text=True,
    ).splitlines()
    result = {}
    # Stream one blob at a time to avoid copying large historical archives into
    # a single process output buffer. No historical file is opened for writing.
    with subprocess.Popen(
        ["git", "cat-file", "--batch"], stdin=subprocess.PIPE, stdout=subprocess.PIPE
    ) as process:
        for name in names:
            process.stdin.write(f"{base}:{name}\n".encode())
            process.stdin.flush()
            header = process.stdout.readline().split()
            require(
                len(header) == 3 and header[1] == b"blob",
                "Missing historical blob: " + name,
            )
            blob = process.stdout.read(int(header[2]))
            require(process.stdout.read(1) == b"\n", "Invalid Git blob framing")
            expected = hashlib.sha256(blob).hexdigest()
            require(expected == sha(name), "Changed historical bytes: " + name)
            result[name] = expected
        process.stdin.close()
        require(process.wait() == 0, "Historical Git read failed")
    return result


def validate_receipt(path):
    root = Path(path).parent
    receipt, manifest = read(path), read(root / "manifest.json")
    for value in (receipt, manifest):
        validate_identity(value)
    require(
        receipt["schema"] == "island.j38-c2-receipt.v1"
        and manifest["schema"] == "island.j38-c2-manifest.v1",
        "Wrong receipt/manifest schema",
    )
    require(receipt["base_commit"] == manifest["base_commit"] == BASE, "Changed base")
    require(
        receipt["manifest_identity"] == manifest["identity"],
        "Manifest/receipt disagreement",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["acceptance_met"] is False
        and receipt["new_md_executions"] == 0,
        "Unsupported readiness, execution or acceptance claim",
    )
    require(
        receipt["unmet_cases"] == ["gaff2/aromatic", "gaff2/pe", "oplsaa/aromatic"],
        "Relabeled failed case",
    )
    files = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p.name not in {"receipt.json", "manifest.json"}
    }
    require(set(manifest["files"]) == files, "Unbound evidence artifact")
    for name, expected in manifest["files"].items():
        target = (root / name).resolve()
        require(
            target.is_relative_to(root.resolve()) and sha(target) == expected,
            "Changed artifact: " + name,
        )
    validate_bundle(root / "run")
    require(
        receipt["matrix_identity"] == read(root / "run/matrix.json")["identity"],
        "Changed matrix binding",
    )
    require(
        read(root / "historical-hashes.json") == historical_hashes(),
        "Incomplete historical byte inventory",
    )
    implementation = receipt["tested_implementation_commit"]
    require(
        implementation == manifest["tested_implementation_commit"],
        "Wrong tested implementation",
    )
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE, implementation], check=True
    )
    require(
        set(manifest["implementation_files"]) == SOURCES,
        "Missing implementation source",
    )
    for file, expected in manifest["implementation_files"].items():
        blob = subprocess.check_output(["git", "show", f"{implementation}:{file}"])
        require(
            hashlib.sha256(blob).hexdigest() == expected == sha(file),
            "Untested implementation file",
        )
    require(
        set(receipt["tests"]) == {"focused", "j22-j38-c1", "full_supported"},
        "Missing test group",
    )
    for name, counts in receipt["tests"].items():
        tree = ET.parse(root / "validation" / (name + ".xml")).getroot()
        actual = {
            k: sum(int(s.get(k, 0)) for s in tree.iter("testsuite"))
            for k in ("tests", "failures", "errors", "skipped")
        }
        require(
            counts == actual
            and counts["tests"] > counts["skipped"]
            and counts["failures"] == counts["errors"] == 0,
            "Invalid test counts",
        )
        command = read(root / "validation" / (name + "-command.json"))
        require(
            command["implementation_commit"] == implementation
            and command["exit_code"] == 0,
            "Wrong tested commit or test result",
        )
    require(set(receipt["validation_commands"]) == COMMANDS, "Missing validation gate")
    for name, expected in receipt["validation_commands"].items():
        require(
            expected == (1 if name in {"c1_acceptance", "raw_driver"} else 0),
            "Wrong expected gate semantics",
        )
        command = read(root / "validation" / (name + "-command.json"))
        require(
            command["exit_code"] == expected
            and command["implementation_commit"] == implementation,
            "Validation result mismatch",
        )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--bundle", action="store_true")
    parser.add_argument("--diagnostic-only", action="store_true")
    args = parser.parse_args()
    (validate_bundle if args.bundle else validate_receipt)(args.path)
    print("VALID DIAGNOSTIC; source compatibility/asset gates unmet; acceptance=False")
    return 0 if args.diagnostic_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
