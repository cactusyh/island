"""Strict J67 receipt semantics, explicit inventories, and offline numerical replay."""

import argparse
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_augmentation import identity, parent_semantics
from island.pcff_qm_pilot.hmds_precision import read
from scripts.phase_4j67 import BASE, BRANCH, PREREQUISITE
from scripts.replay_phase_4j67 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts

FILES = (
    "src/island/pcff_qm_pilot/hmds_augmentation.py",
    "scripts/phase_4j67.py",
    "scripts/analyze_phase_4j67.py",
    "scripts/replay_phase_4j67.py",
    "scripts/diagnose_phase_4j67.py",
    "scripts/validate_phase_4j67_receipt.py",
    "tests/test_phase_4j67.py",
    "docs/phase_4j67.md",
)
# Fixed before new native outcomes. The file binds every frozen request and draw.
FROZEN = "aaf39ed5a68e8a17ee7a1987b9228c8c0d652c080bb5e9cae34d2e66fc156208"
DECLARATIONS = {
    "schema": "island.j67-receipt.v1",
    "candidate_qualified": False,
    "accepted_coefficients": 0,
    "parameter_completeness": "incomplete",
    "scientific_validation": "not_established",
    "chemical_transfer": "not_established",
    "automatic_assignment": False,
    "new_QM": 4,
    "validation_QM": 0,
    "MD_steps": 0,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "original_author_artifact_status": "not_yet_supplied",
    "author_confirmations": "28 user-attested",
}
RECEIPT_FIELDS = set(DECLARATIONS) | {
    "base_commit",
    "required_J66",
    "branch",
    "manifest_sha256",
    "identity",
    "software_counts",
    "terminal_unique_tests",
    "focused_J67",
    "historical_J54_J66",
    "numerical_reconstruction",
    "training_feasible",
    "decision_identity",
    "historical_files_preserved",
}


def semantics(receipt, manifest, training_feasible=None):
    require(set(receipt) == RECEIPT_FIELDS, "J67_receipt_schema")
    require(
        type(receipt.get("training_feasible")) is bool
        and (
            training_feasible is None
            or receipt["training_feasible"] is training_feasible
        ),
        "J67_training_feasibility_declaration",
    )
    identity(receipt)
    identity(manifest)
    for key, value in DECLARATIONS.items():
        require(
            type(receipt.get(key)) is type(value) and receipt[key] == value,
            "J67_receipt_semantic_mutation",
            key=key,
        )
    require(
        manifest.get("schema") == "island.numerical-evidence-manifest.v1"
        and set(manifest)
        == {
            "schema",
            "identity",
            "base_commit",
            "historical_files",
            "implementation_files",
            "artifacts",
        },
        "J67_manifest_schema",
    )
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["required_J66"] == PREREQUISITE
        and receipt["branch"] == BRANCH,
        "J67_receipt_ancestry",
    )


def audit(repo, root):
    parent = repo / "docs/evidence/phase_4j66"
    parent_semantics(
        *(read(parent / n) for n in ("receipt.json", "manifest.json", "decision.json"))
    )
    require(
        pin(root / "plan/frozen-file-pins.json") == FROZEN,
        "J67_frozen_contract_rewritten",
    )
    frozen = read(root / "plan/frozen-file-pins.json")
    identity(frozen)
    require(
        all(pin(root / "plan" / n) == sha for n, sha in frozen["files"].items()),
        "J67_frozen_request_draw_mutation",
    )
    plan = read(root / "plan/scientific-plan.json")
    identity(plan)
    require(
        plan["base_commit"] == BASE
        and plan["required_J66"] == PREREQUISITE
        and plan["branch"] == BRANCH
        and plan["parent_receipt_sha256"] == pin(parent / "receipt.json"),
        "J67_plan_ancestry",
    )
    historical = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in historical.items()),
        "J67_historical_bytes_changed",
    )
    numerical = replay(repo, root)
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        Counter(nodes) == Counter(read(root / "collection.json")["nodes"])
        and len(nodes) == len(set(nodes)),
        "J67_terminal_unique_coverage",
    )
    require(
        all(r["returncode"] == 0 for r in read(root / "validation/checks.json")),
        "J67_required_check_failed",
    )
    decision = read(root / "decision.json")
    identity(decision)
    require(
        decision["training_feasible"] == numerical["training_feasible"]
        and decision["new_QM"] == 4
        and decision["validation_QM"] == decision["MD_steps"] == 0
        and decision["accepted_coefficients"] == 0
        and not decision["candidate_qualified"]
        and not decision["automatic_assignment"]
        and not decision["production_validated"]
        and decision["simulation_readiness"] == "not_established"
        and decision["original_author_artifact_status"] == "not_yet_supplied",
        "J67_decision_semantics",
    )
    return {
        "historical_files": historical,
        "software_counts": counts,
        "terminal_unique_tests": len(nodes),
        "focused_J67": sum(n.startswith("tests/test_phase_4j67.py::") for n in nodes),
        "historical_J54_J66": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 67))
            for n in nodes
        ),
        "numerical_reconstruction": numerical,
        "training_feasible": numerical["training_feasible"],
        "decision_identity": decision["identity"],
    }


def implementation(repo):
    return {n: pin(repo / n) for n in FILES}


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
            "required_J66": PREREQUISITE,
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
    receipt = read(root / "receipt.json")
    manifest = read(root / "manifest.json")
    semantics(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["implementation_files"] == implementation(repo)
        and manifest["artifacts"] == artifacts(root),
        "J67_manifest_hash_inventory",
    )
    result = audit(repo, root)
    require(
        manifest["historical_files"] == result["historical_files"]
        and all(receipt[k] == v for k, v in result.items() if k != "historical_files"),
        "J67_receipt_result_mutation",
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
