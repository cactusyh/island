"""Strict J68 semantic declarations, explicit inventories and offline receipts."""

import argparse
import tempfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from island.pcff_completion import inventory, pin, require, signed, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_augmentation import identity
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j68 import FROZEN, HYPOTHESES
from scripts.replay_phase_4j68 import BASE, BRANCH, parent_semantics, replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts

FILES = (
    "src/island/pcff_qm_pilot/hmds_span.py",
    "src/island/pcff_qm_pilot/hmds_local_baseline.py",
    "scripts/analyze_phase_4j68.py",
    "scripts/diagnose_phase_4j68.py",
    "scripts/summarize_phase_4j68.py",
    "scripts/replay_phase_4j68.py",
    "scripts/validate_phase_4j68_receipt.py",
    "tests/test_phase_4j68.py",
    "docs/phase_4j68.md",
)
DECLARATIONS = {
    "schema": "island.j68-receipt.v1",
    "training_feasible_original_bounded": False,
    "training_performance_feasible_original_bound_free": False,
    "candidate_qualified": False,
    "accepted_coefficients": 0,
    "parameter_completeness": "incomplete",
    "scientific_validation": "not_established",
    "chemical_transfer": "not_established",
    "automatic_assignment": False,
    "new_QM": 0,
    "validation_QM": 0,
    "MD_steps": 0,
    "validation_activated": False,
    "validation_consumed": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "original_author_artifact_status": "not_yet_supplied",
    "author_confirmations": "28 user-attested",
}
RESULT_FIELDS = {
    "historical_files_preserved",
    "software_counts",
    "terminal_unique_tests",
    "focused_J68",
    "historical_J54_J67",
    "numerical_reconstruction",
    "decision_identity",
}
RECEIPT_FIELDS = (
    set(DECLARATIONS)
    | RESULT_FIELDS
    | {"identity", "base_commit", "required_J67", "branch", "manifest_sha256"}
)


def declarations(record, expected):
    for k, v in expected.items():
        require(
            type(record.get(k)) is type(v) and record[k] == v,
            "J68_scientific_declaration",
            key=k,
        )


def decision_semantics(decision):
    identity(decision)
    expected = {k: v for k, v in DECLARATIONS.items() if k != "schema"}
    expected.update(
        schema="island.j68-decision.v1",
        original_span_with_fixed_baseline="checked_insufficient_for_four_limits",
        unique_physical_cause="not_demonstrated",
        protocol_specific_target_reproduction="fails_simultaneous_training_limits",
        basis_sensitivity="historical_failed_or_unestablished_qualifications_unchanged",
        cross_engine_equivalence="not_claimed",
        hypothesis_results={
            "bridge_angle": "checked_incompatible",
            "SiO_stretch": "checked_incompatible; strong partial residual reduction",
        },
    )
    declarations(decision, expected)
    require(
        set(decision)
        == set(expected)
        | {
            "identity",
            "bounds_only_explanation",
            "next_step",
            "next_evidence_requirement",
            "minimum_bound_relaxation",
            "backend_correctness",
            "J66_original22_feasibility",
            "J67_expanded_incompatibility",
            "J67_timing_deviation",
            "scope",
            "comparison_sha256",
        },
        "J68_decision_schema",
    )


def semantics(receipt, manifest):
    identity(receipt)
    identity(manifest)
    require(set(receipt) == RECEIPT_FIELDS, "J68_receipt_schema")
    declarations(receipt, DECLARATIONS)
    require(
        set(manifest)
        == {
            "schema",
            "identity",
            "base_commit",
            "historical_files",
            "implementation_files",
            "artifacts",
        }
        and manifest["schema"] == "island.numerical-evidence-manifest.v1",
        "J68_manifest_schema",
    )
    require(
        all(
            isinstance(manifest[k], dict) and manifest[k]
            for k in ("historical_files", "implementation_files", "artifacts")
        ),
        "J68_manifest_inventory_schema",
    )
    require(
        receipt["base_commit"]
        == receipt["required_J67"]
        == manifest["base_commit"]
        == BASE
        and receipt["branch"] == BRANCH,
        "J68_ancestry",
    )
    n = receipt["numerical_reconstruction"]
    declarations(
        n,
        {
            "observations": 26,
            "protected": 22,
            "physical_columns": 58,
            "intercept_columns": 1,
            "research_solver_calls": 8,
            "failed_processes_retained": 1,
            "new_QM": 0,
            "MD_steps": 0,
            "validation_untouched": True,
            "accepted_coefficients": 0,
            "offline": True,
        },
    )
    require(
        len(n["classifications"]) == 8
        and set(n["classifications"].values()) == {"checked_lower_bound_above_one"},
        "J68_nested_classifications",
    )


def implementation(repo):
    return {n: pin(repo / n) for n in FILES}


def selected_cases(source, selected, replacement):
    """Resolve overlapping runs transparently, never erase an unrelated failure."""
    normalize = lambda c: (
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
    )
    old = list(source.iter("testcase"))
    excluded = [c for c in old if normalize(c).startswith("tests/test_phase_4j68.py::")]
    retained = [
        c for c in old if not normalize(c).startswith("tests/test_phase_4j68.py::")
    ]
    new = list(replacement.iter("testcase"))
    require(
        Counter(map(normalize, excluded)) == Counter(map(normalize, new)),
        "J68_overlap_node_identity",
    )
    require(
        all(
            c.find("failure") is None and c.find("error") is None
            for c in new + retained
        ),
        "J68_unresolved_terminal_failure",
    )
    require(
        [ET.tostring(c) for c in retained]
        == [ET.tostring(c) for c in selected.iter("testcase")],
        "J68_historical_case_rewrite",
    )
    return [normalize(c) for c in excluded]


def terminal_selection(root):
    record = read(root / "validation/terminal-selection.json")
    identity(record)
    require(
        record["schema"] == "island.j68-terminal-selection.v1"
        and record["source"] == "full-4.xml"
        and record["replacement"] == "validation/J68-terminal.xml"
        and record["selected"] == "validation/full-4-historical.xml",
        "J68_terminal_selection_schema",
    )
    with tempfile.TemporaryDirectory(prefix="J68-raw-suite-") as temp:
        target = Path(temp) / "raw"
        unpack(root / "validation/full-suite-raw.tar.xz", target)
        require(
            inventory(target)
            == read(root / "validation/full-suite-raw-inventory.json"),
            "J68_raw_suite_archive",
        )
        original = target / record["source"]
        require(pin(original) == record["source_sha256"], "J68_original_suite_identity")
        excluded = selected_cases(
            ET.parse(original),
            ET.parse(root / record["selected"]),
            ET.parse(root / record["replacement"]),
        )
        require(
            excluded == record["excluded_nodes"]
            and len(excluded) == 69
            and record["source_returncode"] == 1,
            "J68_terminal_overlap_accounting",
        )
    return True


def audit(repo, root):
    parent = repo / "docs/evidence/phase_4j67"
    parent_semantics(
        *(read(parent / n) for n in ("receipt.json", "manifest.json", "decision.json"))
    )
    require(
        pin(root / "plan/scientific-plan.json") == FROZEN
        and pin(root / "plan/baseline-hypotheses.json") == HYPOTHESES,
        "J68_frozen_contract",
    )
    historical = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in historical.items()),
        "J68_historical_bytes_changed",
    )
    for archive, listing in (
        ("campaign-logs.tar.xz", "campaign-log-file-inventory.json"),
        (
            "validation/preliminary-checks.tar.xz",
            "validation/preliminary-file-inventory.json",
        ),
        (
            "validation/subprocess-diagnostics.tar.xz",
            "validation/subprocess-diagnostic-inventory.json",
        ),
    ):
        with tempfile.TemporaryDirectory(prefix="J68-log-integrity-") as temp:
            target = Path(temp) / "archive"
            unpack(root / archive, target)
            require(
                inventory(target) == read(root / listing), "J68_retained_log_archive"
            )
    decision = read(root / "decision.json")
    decision_semantics(decision)
    require(
        decision["comparison_sha256"] == pin(root / "executed-comparison.json"),
        "J68_comparison_identity",
    )
    numerical = replay(repo, root)
    terminal_selection(root)
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        Counter(nodes) == Counter(read(root / "collection.json")["nodes"])
        and len(nodes) == len(set(nodes)),
        "J68_unique_terminal_coverage",
    )
    require(
        all(r["returncode"] == 0 for r in read(root / "validation/checks.json")),
        "J68_required_checks",
    )
    return {
        "historical_files": historical,
        "software_counts": counts,
        "terminal_unique_tests": len(nodes),
        "focused_J68": sum(n.startswith("tests/test_phase_4j68.py::") for n in nodes),
        "historical_J54_J67": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 68))
            for n in nodes
        ),
        "numerical_reconstruction": numerical,
        "decision_identity": decision["identity"],
    }


def seal(repo, root):
    result = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.numerical-evidence-manifest.v1",
            "base_commit": BASE,
            "historical_files": result["historical_files"],
            "implementation_files": implementation(repo),
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            **DECLARATIONS,
            "base_commit": BASE,
            "required_J67": BASE,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            **{k: v for k, v in result.items() if k != "historical_files"},
            "historical_files_preserved": len(result["historical_files"]),
        }
    )
    semantics(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "software_counts": result["software_counts"],
    }


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantics(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["implementation_files"] == implementation(repo)
        and manifest["artifacts"] == artifacts(root),
        "J68_manifest_inventory",
    )
    result = audit(repo, root)
    require(
        manifest["historical_files"] == result["historical_files"]
        and all(receipt[k] == v for k, v in result.items() if k != "historical_files"),
        "J68_receipt_result",
    )
    return {
        "status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        "software_counts": result["software_counts"],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    a = p.parse_args()
    print((seal if a.seal else validate)(a.repo.resolve(), a.directory.resolve()))
