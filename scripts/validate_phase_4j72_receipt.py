"""Strict J72 reconstruction, explicit inventories and scientific semantics."""

import argparse
from collections import Counter
from copy import deepcopy
from pathlib import Path
from xml.etree import ElementTree as ET

from island.pcff_completion import canonical, pin, require, signed, write
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.source_adjudication import verify_signed
from scripts.adjudicate_phase_4j72 import (
    BASE,
    BRANCH,
    DECLARATIONS,
    REQUIRED,
    ROOT,
    replay,
)
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts

FILES = (
    "src/island/pcff_qm_pilot/source_adjudication.py",
    "scripts/adjudicate_phase_4j72.py",
    "scripts/check_phase_4j72_sources.py",
    "scripts/phase_4j72_execution_fence.py",
    "scripts/phase_4j72_test_policy.py",
    "scripts/validate_phase_4j72_receipt.py",
    "tests/test_phase_4j72.py",
    "docs/phase_4j72.md",
)


def offline_declarations():
    return {
        "offline": True,
        **DECLARATIONS,
        "records": 467,
        "numerical_cells": 1470,
        "instances": 378,
        "resolved_instances": 0,
        "real_duplicate_equivalences": 16,
        "unequal_conflicts_resolved": 0,
        "converter_executed": False,
        "static_complete_closure": False,
    }


def exact_typed(value, expected, error):
    require(set(value) == set(expected), error)
    for key, reference in expected.items():
        require(
            type(value[key]) is type(reference) and value[key] == reference,
            error,
            field=key,
        )


def semantics(receipt, manifest):
    verify_signed(receipt)
    verify_signed(manifest)
    require(
        set(receipt)
        == set(DECLARATIONS)
        | {
            "identity",
            "schema",
            "base_commit",
            "required_J71",
            "branch",
            "manifest_sha256",
            "software_counts",
            "terminal_unique_tests",
            "offline_reconstruction",
            "historical_files_preserved",
        },
        "J72_receipt_schema",
    )
    require(
        receipt["schema"] == "island.j72-receipt.v1"
        and receipt["base_commit"] == BASE
        and receipt["required_J71"] == REQUIRED
        and receipt["branch"] == BRANCH,
        "J72_receipt_ancestry",
    )
    exact_typed(
        {k: receipt[k] for k in DECLARATIONS},
        DECLARATIONS,
        "J72_scientific_declarations",
    )
    require(
        set(manifest)
        == {
            "schema",
            "identity",
            "base_commit",
            "implementation_files",
            "historical_files",
            "artifacts",
            "scientific_scope",
        }
        and manifest["schema"] == "island.j72-evidence-manifest.v1"
        and manifest["base_commit"] == BASE,
        "J72_manifest_schema",
    )
    exact_typed(
        manifest["scientific_scope"], DECLARATIONS, "J72_manifest_scientific_semantics"
    )
    exact_typed(
        receipt["offline_reconstruction"],
        offline_declarations(),
        "J72_nested_semantics",
    )
    counts = receipt["software_counts"]
    require(
        set(counts) == {"tests", "skipped", "failure", "error"}
        and all(type(v) is int and v >= 0 for v in counts.values())
        and counts["failure"] == counts["error"] == 0
        and counts["skipped"] <= counts["tests"]
        and counts["tests"] == receipt["terminal_unique_tests"],
        "J72_software_count_schema",
    )
    require(
        type(receipt["terminal_unique_tests"]) is int
        and receipt["terminal_unique_tests"] > 0
        and type(receipt["historical_files_preserved"]) is int
        and receipt["historical_files_preserved"] > 0,
        "J72_verification_counts",
    )


def case_node(case):
    return case.attrib["classname"].replace(".", "/") + ".py::" + case.attrib["name"]


def terminal_sources(root, nodes):
    """Select the entire corrected file, retaining every original failed result."""
    raw = {}
    for i in range(8):
        name = f"validation/guarded-final/full-{i}.xml"
        for case in ET.parse(root / name).iter("testcase"):
            node = case_node(case)
            require(node not in raw, "J72_duplicate_baseline_test")
            raw[node] = (name, case)
    new = {n for n in nodes if n.startswith("tests/test_phase_4j72.py::")}
    affected = {
        n
        for n in nodes
        if n.startswith(("tests/test_phase_4j35_c1.py::", "tests/test_phase_4j37.py::"))
    }
    require(set(raw) == set(nodes) - new, "J72_baseline_test_inventory")
    name = "validation/guard-corrections-terminal.xml"
    corrected = {}
    for case in ET.parse(root / name).iter("testcase"):
        node = case_node(case)
        require(node not in corrected, "J72_duplicate_corrected_test")
        corrected[node] = (name, case)
    require(set(corrected) == new | affected, "J72_whole_file_correction_inventory")
    selected = {n: corrected[n] if n in corrected else raw[n] for n in nodes}
    require(
        all(
            c.find("failure") is None and c.find("error") is None
            for _, c in selected.values()
        ),
        "J72_terminal_test_failure",
    )
    return selected


def terminal_xml_closure(root, nodes):
    selected = terminal_sources(root, nodes)
    require(
        read(root / "validation/junit-files.json")
        == ["validation/terminal-results.xml"],
        "J72_terminal_xml_inventory",
    )
    require(
        read(root / "validation/terminal-selection.json")
        == {n: source for n, (source, _) in selected.items()},
        "J72_terminal_source_mapping",
    )
    cases = list(ET.parse(root / "validation/terminal-results.xml").iter("testcase"))
    require([case_node(c) for c in cases] == nodes, "J72_terminal_node_order")
    for case in cases:
        original = deepcopy(selected[case_node(case)][1])
        case.tail = original.tail = None
        require(
            ET.tostring(case) == ET.tostring(original), "J72_terminal_result_mutation"
        )


def audit(repo, root):
    result = replay(repo, root)
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == h for n, h in old.items()),
        "J72_historical_byte_mutation",
    )
    terminal_xml_closure(root, read(root / "collection.json")["nodes"])
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        len(nodes) == len(set(nodes))
        and Counter(nodes) == Counter(read(root / "collection.json")["nodes"]),
        "J72_unique_terminal_coverage",
    )
    require(
        counts["failure"] == counts["error"] == 0
        and all(c["returncode"] == 0 for c in read(root / "validation/checks.json")),
        "J72_required_checks",
    )
    return {
        "software_counts": counts,
        "terminal_unique_tests": len(nodes),
        "offline_reconstruction": result,
        "historical_files_preserved": len(old),
    }, old


def seal(repo, root):
    result, old = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.j72-evidence-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "historical_files": old,
            "artifacts": artifacts(root),
            "scientific_scope": DECLARATIONS,
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j72-receipt.v1",
            **DECLARATIONS,
            "base_commit": BASE,
            "required_J71": REQUIRED,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantics(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {"status": "sealed", "receipt_sha256": pin(root / "receipt.json"), **result}


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantics(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["implementation_files"] == {n: pin(repo / n) for n in FILES}
        and manifest["artifacts"] == artifacts(root),
        "J72_manifest_inventory",
    )
    result, old = audit(repo, root)
    require(
        manifest["historical_files"] == old
        and all(canonical(receipt[k]) == canonical(v) for k, v in result.items()),
        "J72_receipt_reconstruction",
    )
    return {"status": "passed", "receipt_sha256": pin(root / "receipt.json"), **result}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path, nargs="?", default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    a = p.parse_args()
    print((seal if a.seal else validate)(a.repo.resolve(), a.directory.resolve()))
