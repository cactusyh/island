"""Seal/replay additive J62 evidence with existing history and JUnit contracts."""

import argparse
import json
import lzma
from collections import Counter
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from scripts.replay_phase_4j62 import reconstruct
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

BASE = "1604a19f4f6245a5e35b302467c519aa62e7814d"
PREREQUISITE = "b4cae05e3666d6adb232e0bad711db85292994c7"
BRANCH = "codex/phase-4j62-angular-torsional-control-hmds-entry"
PLAN = "2fbb54e771e865d233b004153eb37799acca94f6aa003cfc956522b641f5e5d5"


def read(path):
    r = json.loads(Path(path).read_text())
    require(
        r["identity"] == digest({k: v for k, v in r.items() if k != "identity"}),
        "J62_record_identity_mutation",
        file=str(path),
    )
    return r


def implementation(repo):
    files = [
        *(repo / "scripts").glob("*phase_4j62*.py"),
        repo / "tests/test_phase_4j62.py",
        repo / "docs/phase_4j62.md",
    ]
    files.extend(
        repo / ("src/island/pcff_qm_pilot/" + n + ".py")
        for n in (
            "angular_baseline",
            "angular_fit",
            "molecular_control",
            "molecular_native",
            "control_engine",
        )
    )
    return {str(p.relative_to(repo)): pin(p) for p in sorted(files)}


def artifacts(root):
    return {
        str(p.relative_to(root)): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file() and p.name not in ("receipt.json", "manifest.json")
    }


def audit(root, repo):
    plan = read(root / "plan/scientific-plan.json")
    require(
        pin(root / "plan/scientific-plan.json") == PLAN
        and plan["base_commit"] == BASE
        and plan["J61_prerequisite"] == PREREQUISITE,
        "J62_frozen_plan_mutation",
    )
    for file in ("implementation-pins.json", "model-implementation-pins.json"):
        for path, sha in read(root / "plan" / file)["files"].items():
            require(
                pin(repo / path) == sha, "J62_frozen_implementation_mutation", file=path
            )
    coordinator = read(root / "plan/coordinator-pin.json")
    amendment = read(root / "plan/coordinator-maintenance-amendment.json")
    original = lzma.decompress(
        (root / amendment["preserved_executed_source"]).read_bytes()
    )
    require(
        sha256(original).hexdigest() == coordinator["sha256"]
        and pin(repo / coordinator["file"]) == amendment["current_source_sha256"],
        "J62_coordinator_maintenance_provenance_changed",
    )
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in old.items()),
        "J62_historical_bytes_changed",
    )
    decision = read(root / "fitting-decision.json")
    require(
        decision["accepted_coefficients"] == {}
        and decision["parameter_completeness"] == "incomplete"
        and decision["production_validated"] is False
        and decision["simulation_readiness"] == "not_established"
        and decision["PCFF_MD_executed"] is False
        and decision["new_registered_profiles"] == []
        and decision["original_author_artifact_status"] == "not_yet_supplied",
        "J62_scope_or_readiness_promoted",
    )
    control = read(root / "angular-torsional-control.json")
    require(
        control["reference_precision_passed"] == all(control["gates"].values()),
        "J62_precision_claim_changed",
    )
    selection = read(root / "nested-model-selection.json")
    for name in selection["models"]:
        trial = read(root / "nested" / f"{name}.json")
        gates = trial["training_gates"]
        n = trial["parameter_count"]
        ref = np.array(trial["nonstretch_reference_column_RMS"])
        err = np.array(trial["nonstretch_response_column_RMS_errors"])
        expected = {
            "optimizer_converged": trial["optimizer"]["success"],
            "full_rank_and_condition": trial["rank"] == n
            and trial["column_scaled_condition"] < 1000,
            "no_active_bounds": not trial["active_bounds"],
            "training_energy": trial["training"]["energy_RMSE"] <= 0.1,
            "training_forces": trial["training"]["force_RMS"] <= 0.5,
            "stationarity": trial["stationary_force_RMS"] <= 0.05,
            "combined_curvature": trial["projected_curvature_relative_error"] <= 0.05,
            "nonstretch_full_response": bool(np.all(err <= 0.5 + 0.05 * ref)),
        }
        require(
            gates == expected and trial["training_qualifies"] == all(expected.values()),
            "J62_nested_acceptance_mutation",
            model=name,
        )
    hmds = read(root / "HMDS-center.json")
    status = hmds["optimization_result"]["status"]
    stationary = (
        status == "completed"
        and hmds["gradient_max_hartree_bohr"] <= 2.5e-5
        and hmds["gradient_RMS_hartree_bohr"] <= 1e-5
    )
    require(
        hmds["stationary_candidate"] == stationary
        and not hmds["full_minimum_verified"],
        "J62_HMDS_stationarity_claim_changed",
    )
    compute = read(root / "compute-accounting.json")
    require(
        compute["attempted_energy_gradient_evaluations"] <= 48
        and compute["primary_control_gradients"] == 16,
        "J62_compute_budget_mutation",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(v["returncode"] == 0 for v in checks), "J62_required_check_failed")
    suites = json.loads((root / "validation/junit-files.json").read_text())
    counts, nodes = junit(root, suites)
    collection = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(collection) and len(nodes) == len(set(nodes)),
        "J62_suite_coverage_mismatch",
    )
    replay = reconstruct(repo, root)
    require(
        replay["native_samples_reconstructed"] == compute["completed_gradients"],
        "J62_native_count_changed",
    )
    return {
        "historical_files": old,
        "software_counts": counts,
        "focused_tests": sum(n.startswith("tests/test_phase_4j62.py::") for n in nodes),
        "J54_J61_regressions": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 62))
            for n in nodes
        ),
        "decision_identity": decision["identity"],
        "reconstruction": replay,
        "new_energy_gradient_evaluations": compute[
            "attempted_energy_gradient_evaluations"
        ],
    }


def seal(root, repo):
    result = audit(root, repo)
    write(
        root / "manifest.json",
        signed(
            {
                "schema": "island.j62-manifest.v1",
                "base_commit": BASE,
                "historical_files": result["historical_files"],
                "implementation_files": implementation(repo),
                "artifacts": artifacts(root),
            }
        ),
    )
    receipt = signed(
        {
            "schema": "island.j62-receipt.v1",
            "base_commit": BASE,
            "J61_prerequisite": PREREQUISITE,
            "J61_contained_in_base": True,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            **{k: v for k, v in result.items() if k != "historical_files"},
            "historical_files_preserved": len(result["historical_files"]),
            "accepted_coefficients": {},
            "parameter_completeness": "incomplete",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "PCFF_MD_executed": False,
            "new_registered_profiles": [],
            "user_attested_confirmations": 28,
            "original_author_artifact_status": "not_yet_supplied",
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "software_counts": result["software_counts"],
    }


def validate(root, repo, expected=None):
    if expected:
        require(pin(root / "receipt.json") == expected, "J62_receipt_hash_changed")
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["J61_prerequisite"] == PREREQUISITE
        and receipt["branch"] == BRANCH,
        "J62_ancestry_mutation",
    )
    require(
        pin(root / "manifest.json") == receipt["manifest_sha256"]
        and manifest["artifacts"] == artifacts(root)
        and manifest["implementation_files"] == implementation(repo),
        "J62_manifest_mutation",
    )
    result = audit(root, repo)
    require(
        result["historical_files"] == manifest["historical_files"]
        and all(receipt[k] == v for k, v in result.items() if k != "historical_files"),
        "J62_receipt_result_mutation",
    )
    require(
        receipt["accepted_coefficients"] == {}
        and receipt["parameter_completeness"] == "incomplete"
        and receipt["scientific_validation"] == "not_established"
        and receipt["cross_engine_equivalence"] == "not_claimed"
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["PCFF_MD_executed"] is False
        and receipt["new_registered_profiles"] == [],
        "J62_receipt_readiness_promoted",
    )
    return {
        "status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(result["historical_files"]),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    p.add_argument("--expected-sha256")
    a = p.parse_args()
    print(
        json.dumps(
            seal(a.directory.resolve(), a.repo.resolve())
            if a.seal
            else validate(a.directory.resolve(), a.repo.resolve(), a.expected_sha256),
            sort_keys=True,
        )
    )
