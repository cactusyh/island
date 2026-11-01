"""Seal or independently replay J64 using existing manifests/history/JUnit contracts."""

import argparse
import ast
import json
import lzma
from collections import Counter
from hashlib import sha256
from pathlib import Path

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.hmds_precision import read
from scripts.phase_4j64 import BASE, PREREQUISITE
from scripts.replay_phase_4j64 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

PLAN = "6c324b0e8fdc9dae7955324f2d20b9cd3758bda4b5fe807073eaad407d160b1e"
BRANCH = "codex/phase-4j64-hmds-precision-baseline-calibration"


def record(path):
    r = read(path)
    require(
        r["identity"] == digest({k: v for k, v in r.items() if k != "identity"}),
        "J64_record_identity_mutation",
        file=str(path),
    )
    return r


def execution_binding(root, repo):
    plan = record(root / "plan/scientific-plan.json")
    amendment = record(root / "plan/startup-maintenance.json")
    for name, old_sha in plan["execution_sources"].items():
        current = pin(repo / name)
        if current == old_sha:
            continue
        b = amendment["sources"][name]
        archive = root / "plan" / b["original_archive"]
        original = lzma.decompress(archive.read_bytes())
        require(
            b["old_sha256"] == old_sha == sha256(original).hexdigest()
            and b["archive_sha256"] == pin(archive)
            and current == b["new_sha256"],
            "J64_execution_source_mutation",
        )
        before, after = ast.parse(original), ast.parse((repo / name).read_text())
        functions = lambda tree: {
            n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)
        }
        bf, af = functions(before), functions(after)
        if name.endswith("hmds_precision.py"):
            # Only getrusage reporting was added inside gradient's finally block.
            final = next(n for n in ast.walk(af["gradient"]) if isinstance(n, ast.Try))
            require(len(final.finalbody) == 3, "J64_resource_amendment_scope_changed")
            final.finalbody = final.finalbody[2:]
            for function in bf:
                require(
                    ast.dump(bf[function]) == ast.dump(af[function]),
                    "J64_numerical_adapter_logic_changed",
                    function=function,
                )
        else:
            # Requests, stage rules, and all scientific settings must be unchanged.
            for function in ("prepare", "execute"):
                require(
                    ast.dump(bf[function]) == ast.dump(af[function]),
                    "J64_stage_or_request_logic_changed",
                    function=function,
                )
    return plan


def implementation(repo):
    files = [
        *(repo / "scripts").glob("*phase_4j64*.py"),
        repo / "src/island/pcff_qm_pilot/hmds_precision.py",
        repo / "src/island/pcff_qm_pilot/hmds_baseline.py",
        repo / "tests/test_phase_4j64.py",
        repo / "docs/phase_4j64.md",
    ]
    return {str(p.relative_to(repo)): pin(p) for p in sorted(files)}


def artifacts(root):
    return {
        str(p.relative_to(root)): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file() and p.name not in ("receipt.json", "manifest.json")
    }


def audit(root, repo):
    plan = execution_binding(root, repo)
    require(
        pin(root / "plan/scientific-plan.json") == PLAN
        and plan["base_commit"] == BASE
        and plan["prerequisite"] == PREREQUISITE,
        "J64_frozen_plan_or_base_mutation",
    )
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in old.items()),
        "J64_historical_bytes_changed",
    )
    numerical = replay(repo, root)
    compute = record(root / "compute-accounting.json")
    require(
        compute["counted_process_attempts"] == numerical["counted_attempts"] <= 24
        and compute["successful_gradients"] == numerical["native_gradients"]
        and compute["new_SiH_QM"]
        == compute["joint_design_QM"]
        == compute["validation_QM"]
        == 0,
        "J64_budget_or_scope_changed",
    )
    baseline = record(root / "baseline-calibration.json")
    require(
        baseline["source_baseline_reference"]
        == "d9ae6279735fb219c7cdaf447971dde5995980292f2c019b71710213d5b3851c"
        and not baseline["equilibrium_and_higher_order_rows_changed"]
        and not baseline["unknown_coefficients_assigned"]
        and not baseline["full56_coefficient_fit_executed"]
        and baseline["accepted_for_automatic_assignment"] == [],
        "J64_parameter_authority_promoted",
    )
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        Counter(nodes) == Counter(read(root / "collection.json")["nodes"])
        and len(nodes) == len(set(nodes)),
        "J64_terminal_suite_coverage_changed",
    )
    checks = read(root / "validation/checks.json")
    require(all(c["returncode"] == 0 for c in checks), "J64_required_check_failed")
    require(
        (root / "independent-path-and-baseline-checks.json").is_file(),
        "J64_independent_numerical_controls_missing",
    )
    return {
        "historical_files": old,
        "software_counts": counts,
        "focused_tests": sum(n.startswith("tests/test_phase_4j64.py::") for n in nodes),
        "J54_J63_regressions": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 64))
            for n in nodes
        ),
        "reconstruction": numerical,
        "compute_identity": compute["identity"],
        "baseline_identity": baseline["identity"],
    }


def seal(root, repo):
    result = audit(root, repo)
    write(
        root / "manifest.json",
        signed(
            {
                "schema": "island.numerical-evidence-manifest.v1",
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
                "schema": "island.j64-receipt.v1",
                "base_commit": BASE,
                "J63_prerequisite": PREREQUISITE,
                "J63_in_main": True,
                "branch": BRANCH,
                "manifest_sha256": pin(root / "manifest.json"),
                **{k: v for k, v in result.items() if k != "historical_files"},
                "historical_files_preserved": len(result["historical_files"]),
                "parameter_completeness": "incomplete",
                "candidate_validation": "not_established",
                "production_validated": False,
                "simulation_readiness": "not_established",
                "scientific_validation": "not_established",
                "cross_engine_equivalence": "not_claimed",
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
        require(pin(root / "receipt.json") == expected, "J64_receipt_hash_mutation")
    receipt, manifest = record(root / "receipt.json"), record(root / "manifest.json")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["J63_prerequisite"] == PREREQUISITE
        and receipt["branch"] == BRANCH,
        "J64_base_mutation",
    )
    require(
        pin(root / "manifest.json") == receipt["manifest_sha256"]
        and manifest["artifacts"] == artifacts(root)
        and manifest["implementation_files"] == implementation(repo),
        "J64_manifest_mutation",
    )
    result = audit(root, repo)
    require(
        result["historical_files"] == manifest["historical_files"]
        and all(receipt[k] == v for k, v in result.items() if k != "historical_files"),
        "J64_receipt_result_mutation",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["parameter_completeness"] == "incomplete"
        and receipt["PCFF_MD_executed"] is False
        and receipt["new_registered_profiles"] == [],
        "J64_readiness_promoted",
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
