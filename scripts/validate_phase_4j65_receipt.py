"""Seal/replay J65 through the existing manifest, history and JUnit contracts."""

import argparse
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_precision import read
from scripts.phase_4j65 import BASE, BRANCH, PREREQUISITE
from scripts.replay_phase_4j65 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts, record

PLAN = "93043b3bf43ca1c5352d94d2a4822104e6b0ebd7202ba96e697bf0ebb8c9af88"


def implementation(repo):
    files = [
        *(repo / "scripts").glob("*phase_4j65*.py"),
        repo / "src/island/pcff_qm_pilot/hmds_fit.py",
        repo / "tests/test_phase_4j65.py",
        repo / "docs/phase_4j65.md",
    ]
    return {str(f.relative_to(repo)): pin(f) for f in sorted(files)}


def audit(root, repo):
    p = record(root / "plan/scientific-plan.json")
    require(
        pin(root / "plan/scientific-plan.json") == PLAN
        and p["base_commit"] == BASE
        and p["J64_prerequisite"] == PREREQUISITE,
        "J65_frozen_contract_mutation",
    )
    require(
        all(pin(repo / n) == sha for n, sha in p["execution_sources"].items()),
        "J65_execution_source_mutation",
    )
    binding = record(root / "plan/analysis-source-bindings.json")
    require(
        binding["bound_before_native_training"]
        and all(pin(repo / n) == sha for n, sha in binding["sources"].items()),
        "J65_fitting_implementation_mutation",
    )
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in old.items()),
        "J65_historical_bytes_changed",
    )
    numerical = replay(repo, root)
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        Counter(nodes) == Counter(read(root / "collection.json")["nodes"])
        and len(nodes) == len(set(nodes)),
        "J65_terminal_test_coverage",
    )
    require(
        all(r["returncode"] == 0 for r in read(root / "validation/checks.json")),
        "J65_required_check_failed",
    )
    fit = record(root / "trial-fit.json")
    require(fit["fit_executed"], "J65_executed_fit_not_supplied")
    controls = record(root / "targeted-controls.json")
    require(
        controls["historical_decisions_changed"] is False
        and controls["J64_complete_force_fit_authorization"] is False,
        "J65_historical_gate_promoted",
    )
    compute = record(root / "compute-accounting.json")
    require(
        compute["counted_attempts"] == numerical["counted_attempts"] <= 32
        and compute["new_SiH_QM"] == 0
        and compute["PCFF_MD"] is False
        and compute["new_registered_profiles"] == [],
        "J65_compute_or_scope_violation",
    )
    require(
        not fit["validation_authorized"]
        or (root / "untouched-validation.json").exists(),
        "J65_eligible_validation_incomplete",
    )
    return {
        "historical_files": old,
        "software_counts": counts,
        "focused_tests": sum(n.startswith("tests/test_phase_4j65.py::") for n in nodes),
        "J54_J64_regressions": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 65))
            for n in nodes
        ),
        "numerical_reconstruction": numerical,
        "fit_identity": fit["identity"],
        "compute_identity": compute["identity"],
    }


def seal(root, repo):
    results = audit(root, repo)
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
    receipt = signed(
        {
            "schema": "island.j65-receipt.v1",
            "base_commit": BASE,
            "J64_prerequisite": PREREQUISITE,
            "J64_in_main": True,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            **{k: v for k, v in results.items() if k != "historical_files"},
            "historical_files_preserved": len(results["historical_files"]),
            "parameter_completeness": "incomplete",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "scientific_validation": "not_established",
            "cross_engine_equivalence": "not_claimed",
            "PCFF_MD_executed": False,
            "new_registered_profiles": [],
            "original_author_artifact_status": "not_yet_supplied",
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "software_counts": results["software_counts"],
    }


def validate(root, repo):
    r, m = record(root / "receipt.json"), record(root / "manifest.json")
    require(
        r["base_commit"] == m["base_commit"] == BASE
        and r["J64_prerequisite"] == PREREQUISITE
        and r["branch"] == BRANCH,
        "J65_base_mutation",
    )
    require(
        pin(root / "manifest.json") == r["manifest_sha256"]
        and m["artifacts"] == artifacts(root)
        and m["implementation_files"] == implementation(repo),
        "J65_manifest_mutation",
    )
    result = audit(root, repo)
    require(
        m["historical_files"] == result["historical_files"]
        and all(r[k] == v for k, v in result.items() if k != "historical_files"),
        "J65_receipt_result_mutation",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established"
        and r["parameter_completeness"] == "incomplete"
        and r["new_registered_profiles"] == [],
        "J65_readiness_promoted",
    )
    return {
        "status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(result["historical_files"]),
    }


if __name__ == "__main__":
    import json

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    a = p.parse_args()
    print(
        json.dumps(
            (seal if a.seal else validate)(a.directory.resolve(), a.repo.resolve()),
            sort_keys=True,
        )
    )
