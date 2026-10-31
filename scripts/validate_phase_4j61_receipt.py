"""Seal/reconstruct J61 results with existing history, native and JUnit helpers."""

import argparse
import json
import lzma
from collections import Counter
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from scripts.replay_phase_4j61 import reconstruct
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

BASE = "3f7c1cf37fe1a28a12f4469409b3e2afabff35d8"
PREREQUISITE = "b441b295a0a23d39c67c62237acb292a564e5d99"
BRANCH = "codex/phase-4j61-pcff-stationary-baseline-joint-refit"
PLAN = "80cd44639bf18e5fab5afae20f310c51e751e370582b711075e2b94b9dd59534"


def read(path):
    record = json.loads(Path(path).read_text())
    require(
        record["identity"]
        == digest({k: v for k, v in record.items() if k != "identity"}),
        "J61_record_identity_mutation",
        file=str(path),
    )
    return record


def implementation(repo):
    files = [
        *(repo / "scripts").glob("*phase_4j61*.py"),
        repo / "tests/test_phase_4j61.py",
        repo / "docs/phase_4j61.md",
    ]
    files.extend(
        repo / ("src/island/pcff_qm_pilot/" + name + ".py")
        for name in ("stationary_engine", "stretch_baseline", "stretch_fit")
    )
    return {str(p.relative_to(repo)): pin(p) for p in sorted(files)}


def artifacts(root):
    return {
        str(p.relative_to(root)): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file() and p.name not in ("receipt.json", "manifest.json")
    }


def fit_gates(root, plan, fit, curvature):
    """Recompute acceptance from the retained numerical values, not claim flags."""
    joint = fit["baseline_adjustable_joint"]
    criteria = plan["fit"]["acceptance"]
    jacobian = np.array(read(root / "joint-parameter-jacobian.json")["matrix"])
    singular = np.linalg.svd(
        jacobian / np.linalg.norm(jacobian, axis=0), compute_uv=False
    )
    rank = int(sum(singular > singular[0] * 1e-10))
    condition = float(singular[0] / singular[-1])
    bootstrap = read(root / "training-bootstrap.json")
    samples = np.array(bootstrap["coefficient_samples"])
    names = bootstrap["names"]
    parameters = np.array([joint["parameters"][k] for k in names])
    sigma = np.std(samples, axis=0, ddof=1)
    require(
        np.allclose(sigma, bootstrap["sigma"], rtol=1e-10, atol=1e-10),
        "J61_bootstrap_sigma_changed",
    )
    require(
        len(samples) == bootstrap["valid_replicates"]
        and bootstrap["replicates"] == 100
        and bootstrap["seed"] == 4061,
        "J61_bootstrap_count_changed",
    )
    weak = [
        k for k, s, v in zip(names, sigma, parameters, strict=True) if s > 0.5 * abs(v)
    ]
    require(
        weak == bootstrap["weak_relative_sigma_parameters"],
        "J61_uncertainty_claim_changed",
    )
    require(
        rank == joint["rank"] and np.isclose(condition, joint["scaled_condition"]),
        "J61_identifiability_claim_changed",
    )
    restarted = read(root / "trials/source-centered-restart.json")["parameters"]
    difference = np.array([restarted[k] - joint["parameters"][k] for k in names])
    restart_error = float(
        np.max(abs(difference) / np.r_[[1, 1, 100, 100], np.ones(20)])
    )
    reduction = (
        1
        - joint["weighted_residual_SS"]
        / fit["coupling_only_refit_same_constraints"]["weighted_residual_SS"]
    )
    require(
        np.isclose(reduction, fit["training_weighted_SS_reduction_fraction"])
        and np.isclose(
            restart_error, fit["source_centered_restart_scaled_max_difference"]
        ),
        "J61_sensitivity_claim_changed",
    )
    return {
        "reference_curvature": curvature["reference_curvature_gate_passed"],
        "identifiable_full24": rank == 24 and condition < 1000,
        "no_active_bounds": not joint["active_bounds"],
        "fresh_energy": joint["fresh_validation"]["energy_RMSE"]
        <= criteria["validation_energy_RMSE_kcal_mol"],
        "fresh_forces": joint["fresh_validation"]["force_RMS"]
        <= criteria["validation_force_RMS_kcal_mol_A"],
        "stationarity": joint["stationary_force_RMS"]
        <= criteria["stationary_force_RMS_kcal_mol_A"],
        "curvature": joint["projected_curvature_relative_Frobenius"]
        <= criteria["projected_Hessian_relative_Frobenius"],
        "material_training_reduction": reduction
        >= criteria["material_training_weighted_SS_reduction_min_fraction"],
        "bootstrap": len(samples) >= 90 and not weak,
        "source_centered_restart": restart_error <= 1e-5,
    }


def audit(root, repo):
    require(pin(root / "plan/scientific-plan.json") == PLAN, "J61_plan_mutation")
    plan = read(root / "plan/scientific-plan.json")
    require(
        plan["budget"]["maximum_new_energy_gradient_evaluations"] == 120,
        "J61_budget_mutation",
    )
    for path, expected in plan["historical_inputs"].items():
        require(pin(repo / path) == expected, "J61_historical_input_changed", file=path)
    orchestration = json.loads((root / "plan/orchestration.json").read_text())
    require(
        pin(repo / "src/island/pcff_qm_pilot/stationary_engine.py")
        == orchestration["adapter_sha256"]
        and pin(repo / "scripts/phase_4j61.py") == orchestration["coordinator_sha256"],
        "J61_executed_implementation_mutation",
    )
    before = read(root / "plan/pre-validation-implementation.json")
    amendment = read(root / "plan/analysis-retention-amendment.json")
    require(
        all(
            pin(repo / path) == expected
            for path, expected in before["files"].items()
            if path != "scripts/analyze_phase_4j61.py"
        )
        and pin(repo / "scripts/analyze_phase_4j61.py")
        == amendment["current_analyzer_sha256"]
        and sha256(
            lzma.decompress(
                (root / amendment["preserved_previous_source"]).read_bytes()
            )
        ).hexdigest()
        == before["files"]["scripts/analyze_phase_4j61.py"],
        "J61_frozen_fit_implementation_mutation",
    )
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / name) == expected for name, expected in old.items()),
        "J61_historical_bytes_changed",
    )
    stationary = read(root / "stationary-center.json")
    curvature = read(root / "projected-Hessian.json")
    fit = read(root / "baseline-versus-joint-fit.json")
    decision = read(root / "fitting-decision.json")
    require(
        stationary["gradient_calls"] == 13
        and stationary["stationary_candidate"]
        and not stationary["all_direction_minimum_verified"],
        "J61_stationary_claim_mutation",
    )
    require(
        curvature["reference_curvature_gate_passed"]
        == all(curvature["gates"].values()),
        "J61_curvature_claim_mutation",
    )
    require(
        len(fit["baseline_adjustable_joint"]["parameters"]) == 24
        and fit["gates"] == fit_gates(root, plan, fit, curvature)
        and bool(fit["accepted_coefficients"]) == all(fit["gates"].values())
        and decision["accepted_coefficients"] == fit["accepted_coefficients"],
        "J61_fit_claim_mutation",
    )
    require(
        decision["parameter_completeness"] == "incomplete"
        and decision["scientific_validation"] == "not_established"
        and decision["production_validated"] is False
        and decision["simulation_readiness"] == "not_established"
        and decision["PCFF_MD_executed"] is False
        and decision["new_registered_profiles"] == []
        and decision["original_author_artifact_status"] == "not_yet_supplied",
        "J61_readiness_or_provenance_promoted",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(r["returncode"] == 0 for r in checks), "J61_required_check_failed")
    suites = json.loads((root / "validation/junit-files.json").read_text())
    counts, nodes = junit(root, suites)
    collection = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(collection) and len(set(nodes)) == len(nodes),
        "J61_suite_coverage_mismatch",
    )
    replay = reconstruct(repo, root)
    compute = read(root / "compute-accounting.json")
    require(
        compute["attempted_energy_gradient_evaluations"] <= 120
        and compute["attempted_energy_gradient_evaluations"]
        == replay["native_samples_reconstructed"],
        "J61_compute_claim_mutation",
    )
    return {
        "historical_files": old,
        "software_counts": counts,
        "focused_tests": sum(n.startswith("tests/test_phase_4j61.py::") for n in nodes),
        "J54_J60_regressions": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 61))
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
    manifest = signed(
        {
            "schema": "island.j61-manifest.v1",
            "base_commit": BASE,
            "historical_files": result["historical_files"],
            "implementation_files": implementation(repo),
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    fit = read(root / "baseline-versus-joint-fit.json")
    receipt = signed(
        {
            "schema": "island.j61-receipt.v1",
            "base_commit": BASE,
            "J60_prerequisite": PREREQUISITE,
            "J60_contained_in_base": True,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            **{k: v for k, v in result.items() if k != "historical_files"},
            "historical_files_preserved": len(result["historical_files"]),
            "trial_parameters": 24,
            "accepted_coefficients": fit["accepted_coefficients"],
            "fit_completed": True,
            "fit_validated": all(fit["gates"].values()),
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
        require(pin(root / "receipt.json") == expected, "J61_receipt_hash_changed")
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["J60_prerequisite"] == PREREQUISITE
        and receipt["branch"] == BRANCH,
        "J61_ancestry_mutation",
    )
    require(
        pin(root / "manifest.json") == receipt["manifest_sha256"]
        and manifest["artifacts"] == artifacts(root)
        and manifest["implementation_files"] == implementation(repo),
        "J61_manifest_mutation",
    )
    result = audit(root, repo)
    require(
        result["historical_files"] == manifest["historical_files"],
        "J61_history_mutation",
    )
    require(
        all(receipt[k] == v for k, v in result.items() if k != "historical_files"),
        "J61_receipt_result_mutation",
    )
    require(
        receipt["accepted_coefficients"]
        == read(root / "fitting-decision.json")["accepted_coefficients"]
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["parameter_completeness"] == "incomplete"
        and receipt["cross_engine_equivalence"] == "not_claimed"
        and receipt["scientific_validation"] == "not_established"
        and receipt["PCFF_MD_executed"] is False
        and receipt["new_registered_profiles"] == [],
        "J61_receipt_readiness_promoted",
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
    parser.add_argument("--expected-sha256")
    args = parser.parse_args()
    print(
        json.dumps(
            seal(args.directory.resolve(), args.repo.resolve())
            if args.seal
            else validate(
                args.directory.resolve(), args.repo.resolve(), args.expected_sha256
            ),
            sort_keys=True,
        )
    )
