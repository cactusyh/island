"""Strict J69 receipt semantics, source replay and historical byte validation."""

import argparse
import tempfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from island.pcff_completion import inventory, pin, require, signed, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_augmentation import identity
from island.pcff_qm_pilot.hmds_precision import read
from scripts.audit_phase_4j69 import BASE, REQUIRED, ROOT, inputs
from scripts.replay_phase_4j69 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts

FILES = (
    "src/island/pcff_qm_pilot/hmds_source_audit.py",
    "scripts/audit_phase_4j69.py",
    "scripts/replay_phase_4j69.py",
    "scripts/validate_phase_4j69_receipt.py",
    "scripts/phase_4j69_test_policy.py",
    "tests/test_phase_4j69.py",
    "docs/phase_4j69.md",
)
DECLARATIONS = {
    "category": "source_incomplete",
    "source_complete": False,
    "available_source_subset_closed": True,
    "source_row_mismatch_demonstrated": False,
    "new_internal_convention_mismatch_demonstrated": False,
    "training_feasible_original_bounded": False,
    "training_performance_feasible_original_bound_free": False,
    "candidate_qualified": False,
    "accepted_coefficients": 0,
    "parameter_completeness": "incomplete",
    "scientific_validation": "not_established",
    "chemical_transfer": "not_established",
    "automatic_assignment": False,
    "new_QM": 0,
    "MD_steps": 0,
    "validation_QM": 0,
    "validation_activated": False,
    "validation_consumed": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "original_author_artifact_status": "not_yet_supplied",
    "author_confirmations": "28 user-attested",
    "hypotheses_frozen": [],
    "hypotheses_executed": 0,
    "research_solver_calls": 0,
    "implementation_correction_justified": False,
}
DECISION_EXTRA = {
    "historical_converter_phase_difference",
    "J66_original22_feasibility",
    "J67_expanded_incompatibility",
    "J68_diagnostic_failures",
    "J67_timing_deviation",
    "affected_rows",
    "affected_residuals",
    "next_step",
    "new_QM_informative_now",
    "blocked",
    "diagnostic_correction",
}


def declarations(record):
    identity(record)
    for k, v in DECLARATIONS.items():
        require(
            type(record.get(k)) is type(v) and record[k] == v,
            "J69_scientific_declaration",
            key=k,
        )


def decision_semantics(record):
    declarations(record)
    require(
        set(record) == set(DECLARATIONS) | DECISION_EXTRA | {"schema", "identity"}
        and record["schema"] == "island.j69-decision.v1",
        "J69_decision_schema",
    )
    require(
        record["new_QM_informative_now"] is False
        and all(
            record[k] == "preserved"
            for k in (
                "J66_original22_feasibility",
                "J67_expanded_incompatibility",
                "J68_diagnostic_failures",
            )
        ),
        "J69_historical_declarations",
    )


def semantics(receipt, manifest):
    declarations(receipt)
    identity(manifest)
    extra = {
        "schema",
        "identity",
        "base_commit",
        "required_J68",
        "branch",
        "manifest_sha256",
        "software_counts",
        "terminal_unique_tests",
        "numerical_reconstruction",
        "decision_identity",
        "historical_files_preserved",
    }
    require(
        set(receipt) == set(DECLARATIONS) | extra
        and receipt["schema"] == "island.j69-receipt.v1",
        "J69_receipt_schema",
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
        }
        and manifest["schema"] == "island.numerical-evidence-manifest.v1",
        "J69_manifest_schema",
    )
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["required_J68"] == REQUIRED
        and receipt["branch"] == "codex/phase-4j69-authoritative-sio-source-audit",
        "J69_ancestry",
    )
    n = receipt["numerical_reconstruction"]
    expected = {
        "offline": True,
        "observations": 26,
        "protected": 22,
        "physical_columns": 58,
        "energy_only_intercepts": 1,
        "source_available_subset_closed": True,
        "source_complete": False,
        "new_QM": 0,
        "MD_steps": 0,
        "validation_untouched": True,
        "hypotheses_executed": 0,
        "accepted_coefficients": 0,
    }
    require(
        n == expected and all(type(n[k]) is type(v) for k, v in expected.items()),
        "J69_nested_scientific_semantics",
    )
    require(
        all(
            isinstance(manifest[k], dict) and manifest[k]
            for k in ("implementation_files", "historical_files", "artifacts")
        ),
        "J69_manifest_inventory_schema",
    )


def terminal_selection(root):
    """Preserve full raw runs and replace only the overlapping J69 file once."""
    record = read(root / "validation/terminal-selection.json")
    identity(record)
    require(
        record["schema"] == "island.j69-terminal-selection.v1"
        and record["replacement"] == "validation/J69-terminal.xml",
        "J69_terminal_selection_schema",
    )
    normalize = lambda c: (
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
    )
    excluded = []
    with tempfile.TemporaryDirectory(prefix="J69-suite-integrity-") as temp:
        raw = Path(temp) / "raw"
        unpack(root / "validation/full-suite-raw.tar.xz", raw)
        require(
            inventory(raw) == read(root / "validation/full-suite-raw-inventory.json"),
            "J69_full_suite_archive",
        )
        for r in record["originals"]:
            require(pin(raw / r["source"]) == r["sha256"], "J69_raw_suite_XML_identity")
            old = list(ET.parse(raw / r["source"]).iter("testcase"))
            keep = [
                c
                for c in old
                if not normalize(c).startswith("tests/test_phase_4j69.py::")
            ]
            excluded += [
                normalize(c)
                for c in old
                if normalize(c).startswith("tests/test_phase_4j69.py::")
            ]
            require(
                [ET.tostring(c) for c in keep]
                == [
                    ET.tostring(c)
                    for c in ET.parse(root / r["selected"]).iter("testcase")
                ],
                "J69_unrelated_terminal_case_rewrite",
            )
    replacement = list(ET.parse(root / record["replacement"]).iter("testcase"))
    require(
        Counter(excluded) == Counter(map(normalize, replacement))
        and len(excluded) == len(set(excluded)) == 68
        and all(
            c.find("failure") is None and c.find("error") is None for c in replacement
        ),
        "J69_terminal_overlap_or_failure",
    )
    return True


def audit(repo, root):
    inputs(repo, root)
    decision = read(root / "decision.json")
    decision_semantics(decision)
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in old.items()),
        "J69_historical_bytes_changed",
    )
    numerical = replay(repo, root)
    terminal_selection(root)
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        len(nodes) == len(set(nodes))
        and Counter(nodes) == Counter(read(root / "collection.json")["nodes"]),
        "J69_unique_terminal_coverage",
    )
    require(
        all(r["returncode"] == 0 for r in read(root / "validation/checks.json")),
        "J69_required_checks",
    )
    return {
        "software_counts": counts,
        "terminal_unique_tests": len(nodes),
        "numerical_reconstruction": numerical,
        "decision_identity": decision["identity"],
        "historical_files_preserved": len(old),
    }, old


def seal(repo, root):
    result, old = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.numerical-evidence-manifest.v1",
            "base_commit": BASE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            **DECLARATIONS,
            "schema": "island.j69-receipt.v1",
            "base_commit": BASE,
            "required_J68": REQUIRED,
            "branch": "codex/phase-4j69-authoritative-sio-source-audit",
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantics(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {"receipt_sha256": pin(root / "receipt.json"), **result}


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantics(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["implementation_files"] == {n: pin(repo / n) for n in FILES}
        and manifest["artifacts"] == artifacts(root),
        "J69_manifest_inventory",
    )
    result, old = audit(repo, root)
    require(
        manifest["historical_files"] == old
        and all(receipt[k] == v for k, v in result.items()),
        "J69_receipt_result",
    )
    return {"status": "passed", "receipt_sha256": pin(root / "receipt.json"), **result}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path, nargs="?", default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    a = p.parse_args()
    print((seal if a.seal else validate)(a.repo.resolve(), a.directory.resolve()))
