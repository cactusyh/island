"""Seal J66 with explicit implementation inventory and offline numerical replay."""

import argparse
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_precision import read
from scripts.phase_4j66 import BASE, BRANCH, PREREQUISITE
from scripts.replay_phase_4j66 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts, record

# Explicit inventory: future research modules must not invalidate this receipt.
FILES = (
    "src/island/pcff_qm_pilot/hmds_constrained.py",
    "scripts/phase_4j66.py",
    "scripts/analyze_phase_4j66.py",
    "scripts/fit_phase_4j66.py",
    "scripts/diagnose_phase_4j66.py",
    "scripts/replay_phase_4j66.py",
    "scripts/validate_phase_4j66_receipt.py",
    "tests/test_phase_4j66.py",
    "docs/phase_4j66.md",
)
FROZEN_PLANS = {
    "scientific-plan.json": "2628ccc395472953f942b5b31be6e2c5d712a5f6123782d4a36d950735ac9e33",
    "objective-epigraph.json": "fbab0fbf8e1557fb4e23d76364ee6bd23688a3b3c761709f75d555b4f8e81b23",
    "numerical-selection.json": "7fbf66ad2872ffedfdeec4fd5e36639b5fd33386f200412e9ca9ff66373d7ac5",
}


def implementation(repo):
    return {name: pin(repo / name) for name in FILES}


def audit(repo, root):
    require(
        all(pin(root / "plan" / name) == sha for name, sha in FROZEN_PLANS.items()),
        "J66_frozen_plan_mutation",
    )
    plan = record(root / "plan/scientific-plan.json")
    require(
        plan["base_commit"] == BASE
        and plan["J65_prerequisite"] == PREREQUISITE
        and plan["branch"] == BRANCH,
        "J66_base_mutation",
    )
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / name) == sha for name, sha in old.items()),
        "J66_historical_byte_mutation",
    )
    numerical = replay(repo, root)
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        Counter(nodes) == Counter(read(root / "collection.json")["nodes"])
        and len(nodes) == len(set(nodes)),
        "J66_unique_terminal_coverage",
    )
    require(
        all(r["returncode"] == 0 for r in read(root / "validation/checks.json")),
        "J66_required_check_failed",
    )
    fit = record(root / "constrained-epigraph.json")
    uncertainty = record(root / "block-uncertainty.json")
    decision = record(root / "decision.json")
    require(
        fit["original_unit_verification"]["original_limits_satisfied"]
        and decision["training_feasible"]
        and decision["accepted_candidates"] == []
        and not uncertainty["qualified"],
        "J66_qualification_mutation",
    )
    require(
        decision["new_QM"] == decision["validation_QM"] == decision["MD_steps"] == 0
        and not decision["production_validated"]
        and decision["simulation_readiness"] == "not_established",
        "J66_scope_mutation",
    )
    return {
        "historical_files": old,
        "software_counts": counts,
        "focused_tests": sum(n.startswith("tests/test_phase_4j66.py::") for n in nodes),
        "J54_J65_regressions": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 66))
            for n in nodes
        ),
        "numerical_reconstruction": numerical,
        "fit_identity": fit["identity"],
        "uncertainty_identity": uncertainty["identity"],
    }


def seal(repo, root):
    results = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.numerical-evidence-manifest.v1",
            "base_commit": BASE,
            "historical_files": results["historical_files"],
            "implementation_files": implementation(repo),
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    write(
        root / "receipt.json",
        signed(
            dict(
                schema="island.j66-receipt.v1",
                base_commit=BASE,
                J65_prerequisite=PREREQUISITE,
                branch=BRANCH,
                manifest_sha256=pin(root / "manifest.json"),
                **{k: v for k, v in results.items() if k != "historical_files"},
                historical_files_preserved=len(results["historical_files"]),
                training_feasible=True,
                candidate_qualified=False,
                parameter_completeness="incomplete",
                explicit_research_trial_vector_complete=True,
                validation_QM=0,
                new_QM=0,
                MD_steps=0,
                production_validated=False,
                simulation_readiness="not_established",
                scientific_validation="not_established",
                chemical_transfer="not_established",
                cross_engine_equivalence="not_claimed",
                automatic_assignment=False,
                original_author_artifact_status="not_yet_supplied",
            )
        ),
    )
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "software_counts": results["software_counts"],
    }


def validate(repo, root):
    receipt = record(root / "receipt.json")
    manifest = record(root / "manifest.json")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["J65_prerequisite"] == PREREQUISITE
        and receipt["branch"] == BRANCH,
        "J66_receipt_base_mutation",
    )
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["implementation_files"] == implementation(repo)
        and manifest["artifacts"] == artifacts(root),
        "J66_manifest_mutation",
    )
    result = audit(repo, root)
    require(
        manifest["historical_files"] == result["historical_files"]
        and all(receipt[k] == v for k, v in result.items() if k != "historical_files"),
        "J66_receipt_mutation",
    )
    require(
        not receipt["candidate_qualified"]
        and not receipt["production_validated"]
        and receipt["simulation_readiness"] == "not_established",
        "J66_readiness_promoted",
    )
    require(
        receipt["new_QM"] == receipt["validation_QM"] == receipt["MD_steps"] == 0
        and not receipt["automatic_assignment"]
        and receipt["parameter_completeness"] == "incomplete",
        "J66_receipt_scope_mutation",
    )
    return {
        "status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(result["historical_files"]),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--seal", action="store_true")
    args = parser.parse_args()
    print(
        (seal if args.seal else validate)(args.repo.resolve(), args.directory.resolve())
    )
