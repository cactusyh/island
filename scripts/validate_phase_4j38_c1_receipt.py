"""Strict additive J38-C1 diagnostic receipt; failed cells never prove acceptance."""

import argparse
import json
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from scripts.phase_4j38_c1_audit import aromatic, require, resp_inventory, sha
from scripts.phase_4j38_c1_contract import BASE, FAILED, validate
from scripts.phase_4j38_contract import validate_identity


def read(path):
    return json.loads(Path(path).read_text())


def validate_bundle(root):
    from island.exceptions import OPLSAssignmentError
    from island.forcefields.oplsaa import load_oplsaa_source, type_atoms
    from island.periodic import periodic_chemical_view
    from scripts.phase_4j32_assignment import prepare
    from scripts.phase_4j32_real_fixture import XML
    from scripts.phase_4j38_fixture import build

    root = Path(root)
    fixtures = {}
    for key in ("pe", "aromatic"):
        fixture, evidence = build(key)
        require(
            read(root / f"{key}-fixture.json") == json.loads(json.dumps(evidence)),
            "Changed exact graph, coordinates, membership, mapping or construction",
        )
        fixtures[key] = fixture
    audit = aromatic(read(root / "aromatic-fixture.json"), XML)
    require(
        audit == read(root / "aromatic-audit.json"),
        "Changed independent source/charge audit",
    )
    typing = type_atoms(
        periodic_chemical_view(fixtures["aromatic"].packed.system),
        load_oplsaa_source(XML),
    )
    require(
        typing.payload["data"]["assignments"]
        == {r["atom_id"]: r["source_type"] for r in audit["rows"]},
        "Native/external typing disagreement",
    )
    with tempfile.TemporaryDirectory(prefix="j38-c1-native-rejection-") as temp:
        try:
            prepare(fixtures["aromatic"], "oplsaa", Path(temp))
        except OPLSAssignmentError as exc:
            rejection = {
                "exception": type(exc).__name__,
                "diagnostic": str(exc),
                "stage": "charge_completeness",
                "classification": "source_incomplete",
                "final_graph_identity": fixtures["aromatic"].packed.graph.identity,
                "dynamics_executed": False,
                "typing_identity": typing.identity,
            }
        else:
            raise ValueError("Unexpected aromatic preparation success")
    require(
        rejection == read(root / "native-rejection.json"), "Changed native rejection"
    )
    inventory = resp_inventory(
        {f"gaff2/{key}": f.packed.graph.identity for key, f in fixtures.items()}
    )
    require(
        inventory == read(root / "resp-inventory.json"),
        "Changed exact RESP source, executable, vector or mapping provenance",
    )
    validate(
        read(root / "matrix.json"),
        read("docs/evidence/phase_4j38/matrix.json"),
        audit,
        inventory,
    )
    require(
        read(root / "run-state.json")
        == {
            "status": "diagnostic_complete",
            "acceptance_met": False,
            "dynamics_executed": False,
            "unmet_cases": sorted(FAILED),
        },
        "Contradictory diagnostic completion",
    )
    return True


def validate_receipt(path):
    root = Path(path).parent
    receipt = read(path)
    manifest = read(root / "manifest.json")
    for value in (receipt, manifest):
        validate_identity(value)
    require(
        receipt["schema"] == "island.j38-c1-receipt.v1"
        and manifest["schema"] == "island.j38-c1-manifest.v1",
        "Wrong schema",
    )
    require(receipt["base_commit"] == manifest["base_commit"] == BASE, "Changed base")
    require(
        receipt["manifest_identity"] == manifest["identity"],
        "Manifest/receipt mismatch",
    )
    require(
        receipt["acceptance_met"] is False
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Unsupported readiness/acceptance claim",
    )
    require(
        receipt["new_md_executions"] == 0 and receipt["unmet_cases"] == sorted(FAILED),
        "Invented dynamics or acceptance",
    )
    actual_files = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p.name not in {"receipt.json", "manifest.json"}
    }
    require(
        set(manifest["files"]) == actual_files, "Missing or unbound receipt artifacts"
    )
    for name, expected in manifest["files"].items():
        target = (root / name).resolve()
        require(
            target.is_relative_to(root.resolve()) and sha(target) == expected,
            "Changed evidence artifact: " + name,
        )
    validate_bundle(root / "run")
    require(
        receipt["matrix_identity"] == read(root / "run/matrix.json")["identity"],
        "Changed matrix binding",
    )
    implementation = receipt["tested_implementation_commit"]
    require(
        implementation == manifest["tested_implementation_commit"],
        "Changed tested implementation",
    )
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE, implementation], check=True
    )
    required_sources = {
        "scripts/phase_4j38_c1_audit.py",
        "scripts/phase_4j38_c1_contract.py",
        "scripts/validate_phase_4j38_c1.py",
        "scripts/validate_phase_4j38_c1_receipt.py",
        "tests/test_phase_4j38_c1.py",
    }
    require(
        set(manifest["implementation_files"]) == required_sources,
        "Missing implementation source",
    )
    for file, expected in manifest["implementation_files"].items():
        blob = subprocess.check_output(["git", "show", f"{implementation}:{file}"])
        import hashlib

        require(
            hashlib.sha256(blob).hexdigest() == expected == sha(file),
            "Untested implementation: " + file,
        )
    historical = read(root / "historical-hashes.json")
    names = subprocess.check_output(
        [
            "git",
            "ls-tree",
            "-r",
            "--name-only",
            BASE,
            "docs/evidence",
            "docs/paper/reproducibility_manifest.json",
            ".gitattributes",
        ],
        text=True,
    ).splitlines()
    require(set(historical) == set(names), "Incomplete historical inventory")
    for file, expected in historical.items():
        blob = subprocess.check_output(["git", "show", f"{BASE}:{file}"])
        require(
            hashlib.sha256(blob).hexdigest() == expected == sha(file),
            "Changed historical bytes: " + file,
        )
    require(
        set(receipt["tests"]) == {"focused", "j22-j38", "full_supported"},
        "Missing test group",
    )
    for name, counts in receipt["tests"].items():
        suite = ET.parse(root / "validation" / (name + ".xml")).getroot()
        actual = {
            k: sum(int(s.get(k, 0)) for s in suite.iter("testsuite"))
            for k in ("tests", "failures", "errors", "skipped")
        }
        require(
            counts == actual
            and actual["tests"] > actual["skipped"]
            and actual["failures"] == actual["errors"] == 0,
            "Invalid JUnit counts",
        )
        command = read(root / "validation" / (name + "-command.json"))
        require(
            command["implementation_commit"] == implementation
            and command["exit_code"] == 0,
            "Tests did not pass on tested implementation",
        )
    required_commands = {
        "ruff",
        "compileall",
        "pip_supported",
        "pip_foyer",
        "pip_analysis",
        "pip_amber",
        "pip_quantum",
        "diff_check",
        "j38_diagnostic",
        "j38_acceptance",
        "historical_receipts",
        "historical_hashes",
        "raw_driver",
        "relocated_bundle",
    }
    require(
        set(receipt["validation_commands"]) == required_commands,
        "Missing verification gate",
    )
    require(
        all(
            code == (1 if name in {"raw_driver", "j38_acceptance"} else 0)
            for name, code in receipt["validation_commands"].items()
        ),
        "Wrong expected command exits",
    )
    for name, expected in receipt["validation_commands"].items():
        command = read(root / "validation" / (name + "-command.json"))
        require(
            command["exit_code"] == expected
            and command["implementation_commit"] == implementation,
            "Wrong validation result: " + name,
        )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--bundle", action="store_true")
    parser.add_argument("--diagnostic-only", action="store_true")
    args = parser.parse_args()
    if args.bundle:
        validate_bundle(args.path)
    else:
        validate_receipt(args.path)
    print("Diagnostic integrity verified; acceptance remains unmet")
    return 0 if args.diagnostic_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
