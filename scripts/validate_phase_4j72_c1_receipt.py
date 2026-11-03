"""Strict additive C1 receipts: typed semantics, quarantine and offline replay."""

import argparse
import lzma
from pathlib import Path
from xml.etree import ElementTree as ET

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.source_adjudication import verify_signed
from scripts.phase_4j72_c1_integrity import (
    BASE,
    BRANCH,
    DECLARATIONS,
    PARENT,
    REQUIRED,
    ROOT,
    read,
    reconstruct,
    typed_equal,
    verify_events,
)
from scripts.validate_phase_4j47_receipt import pin as historical_pin

FILES = (
    "scripts/phase_4j72_c1_fence.py",
    "scripts/phase_4j72_c1_supervisor.py",
    "scripts/run_phase_4j72_c1.py",
    "scripts/phase_4j72_c1_test_policy.py",
    "scripts/phase_4j72_c1_integrity.py",
    "scripts/validate_phase_4j72_c1_receipt.py",
    "scripts/verify_phase_4j72_c1.py",
    "tests/test_phase_4j72_c1.py",
    "docs/phase_4j72_c1.md",
)
RESULT_KEYS = {
    "software_counts",
    "terminal_unique_tests",
    "historical_files_preserved",
    "fence_verification",
    "offline_reconstruction",
}


def semantics(receipt, manifest):
    verify_signed(receipt)
    verify_signed(manifest)
    require(
        set(receipt)
        == set(DECLARATIONS)
        | RESULT_KEYS
        | {
            "identity",
            "schema",
            "base_commit",
            "required_J72",
            "branch",
            "manifest_sha256",
            "quarantine_manifest_sha256",
        },
        "C1_receipt_schema",
    )
    require(
        receipt["schema"] == "island.j72-c1-receipt.v1"
        and receipt["base_commit"] == BASE
        and receipt["required_J72"] == REQUIRED
        and receipt["branch"] == BRANCH,
        "C1_receipt_ancestry",
    )
    typed_equal(
        {k: receipt[k] for k in DECLARATIONS}, DECLARATIONS, "C1_receipt_semantics"
    )
    require(
        set(manifest)
        == {
            "schema",
            "identity",
            "base_commit",
            "implementation_files",
            "historical_files_sha256",
            "artifacts",
            "scientific_scope",
            "quarantine_manifest_sha256",
        }
        and manifest["schema"] == "island.j72-c1-manifest.v1"
        and manifest["base_commit"] == BASE,
        "C1_manifest_schema",
    )
    typed_equal(manifest["scientific_scope"], DECLARATIONS, "C1_manifest_semantics")
    require(
        receipt["quarantine_manifest_sha256"] == manifest["quarantine_manifest_sha256"]
        and len(receipt["quarantine_manifest_sha256"]) == 64,
        "C1_quarantine_binding",
    )
    typed_equal(
        receipt["offline_reconstruction"],
        {
            "offline": True,
            "unresolved_instances": 378,
            "duplicate_equivalences": 16,
            "unresolved_unequal_COMPASS_assignments": 2,
            **DECLARATIONS,
        },
        "C1_nested_scope_semantics",
    )
    counts = receipt["software_counts"]
    require(
        set(counts) == {"tests", "passed", "skipped", "failure", "error"}
        and all(type(v) is int and v >= 0 for v in counts.values())
        and counts["failure"] == counts["error"] == 0
        and counts["tests"] == counts["passed"] + counts["skipped"]
        and type(receipt["terminal_unique_tests"]) is int
        and receipt["terminal_unique_tests"] == counts["tests"],
        "C1_software_counts",
    )


def artifact_inventory(root):
    return {
        str(p.relative_to(root)): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file()
        and p.name not in {"manifest.json", "receipt.json"}
        and "__pycache__" not in p.parts
    }


def replay(repo, root):
    expected = reconstruct(repo, root)
    for name, value in expected.items():
        typed_equal(read(root / name), value, "C1_reconstruction:" + name)
    return {
        "offline": True,
        "unresolved_instances": 378,
        "duplicate_equivalences": 16,
        "unresolved_unequal_COMPASS_assignments": 2,
        **DECLARATIONS,
    }


def terminal(repo, root):
    nodes = read(root / "collection.json")["nodes"]
    # Later phases must not broaden this historical inventory. Preserve the
    # exact J72 collection and the executed C1 file, not a future test glob.
    parent_nodes = read(repo / PARENT / "collection.json")["nodes"]
    corrective = {n for n in nodes if n.startswith("tests/test_phase_4j72_c1.py::")}
    require(
        set(nodes) == set(parent_nodes) | corrective and len(corrective) == 145,
        "C1_terminal_collection_scope",
    )
    selection = read(root / "validation/terminal-selection.json")
    require(
        len(nodes) == len(set(nodes)) and set(selection) == set(nodes),
        "C1_terminal_unique_inventory",
    )
    counts = dict.fromkeys(("tests", "passed", "skipped", "failure", "error"), 0)
    all_cases = {}
    for source in set(selection.values()):
        with (
            lzma.open(root / source, "rb")
            if source.endswith(".xz")
            else (root / source).open("rb") as stream
        ):
            cases = list(ET.parse(stream).iter("testcase"))
        for case in cases:
            node = (
                case.attrib["classname"].replace(".", "/")
                + ".py::"
                + case.attrib["name"]
            )
            if selection.get(node) != source:
                continue
            require((source, node) not in all_cases, "C1_duplicate_raw_testcase")
            all_cases[source, node] = case
    for node in nodes:
        case = all_cases[selection[node], node]
        counts["tests"] += 1
        status = next(
            (k for k in ("error", "failure", "skipped") if case.find(k) is not None),
            "passed",
        )
        counts[status] += 1
    require(counts["failure"] == counts["error"] == 0, "C1_terminal_failures")
    return counts, len(nodes)


def audit(repo, root):
    result = replay(repo, root)
    cross = read(root / "available-quarantine-crosscheck.json")
    verify_signed(cross)
    typed_equal(
        {k: v for k, v in cross.items() if k != "identity"},
        {
            "schema": "island.j72-c1-available-quarantine-crosscheck.v1",
            "original_directory": "/tmp/island-j72-validation/pytest-terminal-0",
            "available_files": read(root / "quarantine-manifest.json")["members"],
            "file_count": 192,
            "unarchived_files": [],
            "missing_files": [],
            "all_original_available_files_match": True,
            "scientific_execution": False,
        },
        "C1_available_quarantine_crosscheck",
    )
    history = read(root / "historical-files.json")
    frozen = read(root / "frozen-inputs.json")
    require(
        pin(root / "historical-files.json") == frozen["historical_files_sha256"],
        "C1_historical_inventory_mutation",
    )
    for name, digest in history.items():
        require(
            historical_pin(repo / name) == digest,
            "C1_historical_byte_mutation",
            path=name,
        )
    counts, unique = terminal(repo, root)
    required = read(root / "validation/checks.json")
    expected_checks = {f"J{i}_receipt" for i in range(54, 73)} | {
        "pip-phase4j66",
        "pip-phase4e2",
        "pip-phase4g1",
        "pip-phase4j32",
        "pip-phase4j36",
        "pip-md_rg",
        "pip-analysis-solver-bridge",
        "ruff",
        "compileall",
        "diff",
        "offline-replay",
    }
    require(
        len(required) == len(expected_checks)
        and {c["name"] for c in required} == expected_checks
        and all(c["returncode"] == 0 for c in required),
        "C1_required_checks_failed",
    )
    fence = verify_events(
        [root / n for n in read(root / "validation/fenced-runs.json")]
    )
    return {
        "software_counts": counts,
        "terminal_unique_tests": unique,
        "historical_files_preserved": len(history),
        "fence_verification": fence,
        "offline_reconstruction": result,
    }


def seal(repo, root):
    result = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.j72-c1-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "historical_files_sha256": pin(root / "historical-files.json"),
            "artifacts": artifact_inventory(root),
            "scientific_scope": DECLARATIONS,
            "quarantine_manifest_sha256": pin(root / "quarantine-manifest.json"),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j72-c1-receipt.v1",
            **DECLARATIONS,
            "base_commit": BASE,
            "required_J72": REQUIRED,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            "quarantine_manifest_sha256": pin(root / "quarantine-manifest.json"),
            **result,
        }
    )
    semantics(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {"status": "sealed", **result}


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantics(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["quarantine_manifest_sha256"]
        == pin(root / "quarantine-manifest.json")
        and manifest["historical_files_sha256"] == pin(root / "historical-files.json")
        and manifest["implementation_files"] == {n: pin(repo / n) for n in FILES}
        and manifest["artifacts"] == artifact_inventory(root),
        "C1_manifest_inventory",
    )
    result = audit(repo, root)
    typed_equal(
        {k: receipt[k] for k in RESULT_KEYS}, result, "C1_receipt_reconstruction"
    )
    return {"status": "passed", "receipt_sha256": pin(root / "receipt.json"), **result}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path, nargs="?", default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--replay", action="store_true")
    p.add_argument("--seal", action="store_true")
    a = p.parse_args()
    function = replay if a.replay else seal if a.seal else validate
    print(function(a.repo.resolve(), a.directory.resolve()))
