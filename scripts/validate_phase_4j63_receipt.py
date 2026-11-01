"""Seal/replay J63 numerical evidence using the existing history/JUnit contracts."""

import argparse
import ast
import json
import lzma
from collections import Counter
from hashlib import sha256
from pathlib import Path

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.hmds_hamiltonian import accounting, load_model
from scripts.replay_phase_4j63 import reconstruct
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

BASE = "251e17db8a19b652196e559409ddd601fc20ceef"
PREREQUISITE = "197a4b8921cf8180124485d9f0ba840c2853f7af"
BRANCH = "codex/phase-4j63-hmds-stationary-coupled-mm"
PLAN = "ad92b223c2acb4fd3265a5c7800107f133fe8f648d7b3a2503cd6a40c92218a8"


def read(p):
    r = json.loads(Path(p).read_text())
    require(
        r["identity"] == digest({k: v for k, v in r.items() if k != "identity"}),
        "J63_record_identity_mutation",
        file=str(p),
    )
    return r


def execution_source(repo, root):
    pins = read(root / "plan/execution-pins.json")
    for name, original_hash in pins["files"].items():
        if pin(repo / name) == original_hash:
            continue
        amendment = read(root / "plan/engine-maintenance-amendment.json")
        original = lzma.decompress(
            (root / amendment["preserved_executed_source"]).read_bytes()
        )
        require(
            sha256(original).hexdigest() == original_hash
            and pin(repo / name) == amendment["current_source_sha256"],
            "J63_executed_source_provenance_changed",
        )
        before, after = ast.parse(original), ast.parse((repo / name).read_text())
        # Reconstruct only the declared import/status maintenance AST. Numerical
        # kernels, settings, callbacks, output capture and every other node must match.
        before.body = [
            n
            for n in before.body
            if not (isinstance(n, ast.Import) and [a.name for a in n.names] == ["os"])
        ]
        if "status_repair_helper" in amendment:
            helper_path = root / amendment["status_repair_helper"]
            require(
                pin(helper_path) == amendment["status_repair_helper_sha256"],
                "J63_status_repair_helper_changed",
            )
            helper = ast.parse(helper_path.read_text()).body
            require(
                len(helper) == 1
                and isinstance(helper[0], ast.FunctionDef)
                and helper[0].name == "optimization_budget_exception",
                "J63_status_repair_scope_changed",
            )
            position = next(
                i
                for i, n in enumerate(before.body)
                if isinstance(n, ast.FunctionDef) and n.name == "optimize"
            )
            before.body.insert(position, helper[0])
            for node in before.body:
                if (
                    isinstance(node, ast.Assign)
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == "VERSION"
                ):
                    node.value = ast.Constant("island.psi4-hmds-j63.v1.1")
            optimizer = before.body[position + 1]
            changed = 0
            for node in ast.walk(optimizer):
                if (
                    isinstance(node, ast.Assign)
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == "exhausted"
                ):
                    require(
                        isinstance(node.value, ast.BoolOp), "J63_status_repair_shape"
                    )
                    node.value.values.insert(
                        0,
                        ast.parse(
                            "optimization_budget_exception(exc)", mode="eval"
                        ).body,
                    )
                    changed += 1
            require(changed == 1, "J63_status_repair_scope_changed")
        require(
            ast.dump(before) == ast.dump(after),
            "J63_maintenance_changed_execution_logic",
        )
    amendment_path = root / "plan/engine-maintenance-amendment.json"
    if amendment_path.exists():
        amendment = read(amendment_path)
        coordinator = lzma.decompress(
            (root / amendment["preserved_executed_coordinator"]).read_bytes()
        )
        require(
            sha256(coordinator).hexdigest() == amendment["original_coordinator_sha256"]
            and pin(repo / "scripts/phase_4j63.py")
            == amendment["current_coordinator_sha256"],
            "J63_coordinator_maintenance_provenance_changed",
        )


def implementation(repo):
    files = [
        *(repo / "scripts").glob("*phase_4j63*.py"),
        *(repo / "src/island/pcff_qm_pilot").glob("hmds_*.py"),
        repo / "tests/test_phase_4j63.py",
        repo / "docs/phase_4j63.md",
    ]
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
        and plan["J62_prerequisite"] == PREREQUISITE,
        "J63_plan_or_ancestry_mutation",
    )
    execution_source(repo, root)
    for name, value in read(root / "plan/model-and-analysis-implementation-pins.json")[
        "files"
    ].items():
        require(pin(repo / name) == value, "J63_numerical_model_source_mutation")
    curvature_source = read(root / "plan/unknown-null-stretch-diagnostic.json")
    require(
        pin(repo / curvature_source["file"]) == curvature_source["sha256"],
        "J63_unknown_null_source_mutation",
    )
    trajectory_source = read(root / "plan/optimizer-development-model-diagnosis.json")
    require(
        pin(repo / trajectory_source["file"]) == trajectory_source["sha256"],
        "J63_trajectory_diagnosis_source_mutation",
    )
    if (root / "plan/path-implementation-pins.json").exists():
        for name, value in read(root / "plan/path-implementation-pins.json")[
            "files"
        ].items():
            require(pin(repo / name) == value, "J63_executed_path_source_mutation")
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in old.items()),
        "J63_historical_bytes_changed",
    )
    model = load_model(repo)
    require(
        model == read(root / "Hamiltonian-definition.json"), "J63_Hamiltonian_changed"
    )
    table = accounting(model, plan["initial_coordinates"])
    require(
        json.loads(json.dumps(table)) == read(root / "interaction-accounting.json")
        and table["independent_coefficient_count"] == 56
        and table["directed_missing_query_count"] == 26,
        "J63_interaction_accounting_changed",
    )
    decision = read(root / "fitting-decision.json")
    next_data = read(root / "next-data-requirements.json")
    require(
        decision["accepted_coefficients"] == {}
        and decision["parameter_completeness"] == "incomplete"
        and decision["production_validated"] is False
        and decision["simulation_readiness"] == "not_established"
        and decision["PCFF_MD_executed"] is False
        and decision["new_registered_profiles"] == []
        and decision["original_author_artifact_status"] == "not_yet_supplied"
        and decision["user_attested_confirmations"] == 28
        and decision["next_data_requirements_identity"] == next_data["identity"]
        and next_data["new_QM_or_baseline_parameters_activated"] is False,
        "J63_scope_promoted",
    )
    compute = read(root / "compute-accounting.json")
    trajectory = read(root / "optimizer-model-diagnosis.json")
    require(
        len(trajectory["samples"]) == compute["groups"].get("optimization", 0)
        and trajectory["all56_unknown_contributions_represented"] is True
        and trajectory["physical_coefficients_inferred"] is False
        and trajectory["validation_used"] is False,
        "J63_terminal_development_diagnosis_missing_or_promoted",
    )
    require(
        compute["attempted_energy_gradient_evaluations"] <= 80
        and compute["groups"].get("optimization", 0) <= 40
        and compute["new_SiH_QM_calls"] == 0
        and compute["numerical_Hessian_calls"] == 0,
        "J63_budget_or_domain_mutation",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(c["returncode"] == 0 for c in checks), "J63_required_check_failed")
    counts, nodes = junit(
        root, json.loads((root / "validation/junit-files.json").read_text())
    )
    collection = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(collection) and len(nodes) == len(set(nodes)),
        "J63_suite_coverage_mismatch",
    )
    replay = reconstruct(repo, root)
    require(
        replay["native_samples_reconstructed"] == compute["completed_gradients"]
        and replay["attempted_gradients"]
        == compute["attempted_energy_gradient_evaluations"]
        and replay["stationary_candidate"] == decision["stationary_candidate"]
        and replay["bounded_path_precision_passed"]
        == decision["bounded_reference_precision_passed"],
        "J63_numerical_status_changed",
    )
    return {
        "historical_files": old,
        "software_counts": counts,
        "focused_tests": sum(n.startswith("tests/test_phase_4j63.py::") for n in nodes),
        "J54_J62_regressions": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 63))
            for n in nodes
        ),
        "reconstruction": replay,
        "decision_identity": decision["identity"],
        "independent_coefficient_count": 56,
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
                "schema": "island.j63-manifest.v1",
                "base_commit": BASE,
                "historical_files": result["historical_files"],
                "implementation_files": implementation(repo),
                "artifacts": artifacts(root),
            }
        ),
    )
    write(
        root / "receipt.json",
        signed(
            {
                "schema": "island.j63-receipt.v1",
                "base_commit": BASE,
                "J62_prerequisite": PREREQUISITE,
                "J62_contained_in_base": True,
                "branch": BRANCH,
                "manifest_sha256": pin(root / "manifest.json"),
                **{k: v for k, v in result.items() if k != "historical_files"},
                "historical_files_preserved": len(result["historical_files"]),
                "accepted_coefficients": {},
                "parameter_completeness": "incomplete",
                "scientific_validation": "not_established",
                "cross_engine_equivalence": "not_claimed",
                "production_validated": False,
                "simulation_readiness": "not_established",
                "PCFF_MD_executed": False,
                "new_registered_profiles": [],
                "original_author_artifact_status": "not_yet_supplied",
            }
        ),
    )
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "software_counts": result["software_counts"],
    }


def validate(root, repo, expected=None):
    if expected:
        require(pin(root / "receipt.json") == expected, "J63_receipt_hash_changed")
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["J62_prerequisite"] == PREREQUISITE
        and receipt["branch"] == BRANCH,
        "J63_ancestry_changed",
    )
    require(
        pin(root / "manifest.json") == receipt["manifest_sha256"]
        and manifest["artifacts"] == artifacts(root)
        and manifest["implementation_files"] == implementation(repo),
        "J63_manifest_mutation",
    )
    result = audit(root, repo)
    require(
        result["historical_files"] == manifest["historical_files"]
        and all(receipt[k] == v for k, v in result.items() if k != "historical_files"),
        "J63_receipt_result_mutation",
    )
    require(
        receipt["accepted_coefficients"] == {}
        and receipt["parameter_completeness"] == "incomplete"
        and receipt["scientific_validation"] == "not_established"
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["PCFF_MD_executed"] is False
        and receipt["new_registered_profiles"] == [],
        "J63_receipt_readiness_promoted",
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
