"""Validate J58 real QM evidence, strict intake, test membership and historical bytes."""

import argparse
import json
import tempfile
from collections import Counter
from pathlib import Path

from island.pcff_completion import digest, inventory, pin, require
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.evidence import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

BASE = "f471e06894275b5b6485461d2bc1dcd594f7dfa7"
BRANCH = "codex/phase-4j58-pcff-pdms-executed-qm-pilot"
IMPLEMENTATION = tuple(
    f"src/island/pcff_qm_pilot/{name}.py"
    for name in (
        "__init__",
        "engine",
        "campaign",
        "importer",
        "analysis",
        "recovery",
        "bundle",
        "evidence",
    )
) + (
    "scripts/phase_4j58.py",
    "scripts/independent_phase_4j58.py",
    "scripts/seal_phase_4j58_receipt.py",
    "scripts/validate_phase_4j58_receipt.py",
    "tests/test_phase_4j58.py",
    "docs/phase_4j58.md",
)
SUITES = [f"validation/full-{i}.xml" for i in range(8)] + [
    "validation/focused-final.xml"
]
REQUIRED_CHECKS = (
    {
        "ruff",
        "compileall",
        "diff",
        "historical_hashes",
        "independent_units",
        "relocation_fresh_process",
    }
    | {f"J{j}_receipt" for j in range(54, 58)}
    | {"pip-" + n for n in ("phase4e2", "phase4g1", "phase4j32", "phase4j36", "md_rg")}
)


def read(path):
    value = json.loads(Path(path).read_text())
    require(
        value.get("identity")
        == digest({k: v for k, v in value.items() if k != "identity"}),
        "record_checksum_mutation",
    )
    return value


def audit(root, *, repo, distribution):
    root, repo = Path(root), Path(repo)
    summary = json.loads((root / "summary.json").read_text())
    require(
        summary["QM_execution"] == "passed"
        and summary["QM_data_integrity"] == "passed",
        "executed_qm_acceptance_not_passed",
    )
    require(
        summary["scan_gradient_jobs"] == 50
        and summary["central_hessian_jobs"] == 2
        and summary["sensitivity_gradient_jobs"] == 6
        and summary["finite_difference_energy_jobs"] == 4,
        "pilot_budget_mutation",
    )
    require(
        summary["new_profiles"] == summary["fitted_parameters"] == []
        and not summary["PCFF_MD_executed"],
        "invented_profile_or_md",
    )
    require(
        summary["parameter_completeness"] == "incomplete"
        and summary["cross_engine_equivalence"] == "not_claimed"
        and summary["scientific_validation"]
        == summary["simulation_readiness"]
        == "not_established"
        and summary["production_validated"] is False,
        "promoted_readiness",
    )
    require(
        summary["author_confirmed_audited"] == 28
        and summary["primary_author_artifacts"] == "not_yet_supplied"
        and summary["remaining_qm_fit_required"] == 19,
        "author_or_term_provenance_mutation",
    )
    with tempfile.TemporaryDirectory(prefix="j58-receipt-replay-") as temp:
        temp = Path(temp)
        unpack(root / "native-outputs.tar.xz", temp / "native")
        unpack(root / "reference-datasets.tar.xz", temp / "normalized")
        require(
            inventory(temp / "native")
            == json.loads((root / "native-file-inventory.json").read_text()),
            "native_recursive_inventory_mutation",
        )
        require(
            inventory(temp / "normalized")
            == json.loads((root / "reference-file-inventory.json").read_text()),
            "reference_recursive_inventory_mutation",
        )
        result = replay(repo, distribution, root, temp / "native", temp / "normalized")
    require(
        result == json.loads((root / "relocation-replay.json").read_text()),
        "relocation_replay_mutation",
    )
    independent = json.loads((root / "independent-unit-checks.json").read_text())
    require(
        not independent["production_importer_used"]
        and len(independent["samples"]) == 50
        and independent["checker_sha256"]
        == pin(repo / "scripts/independent_phase_4j58.py"),
        "independent_importer_check_missing",
    )
    numerical = json.loads((root / "numerical-checks.json").read_text())
    require(
        len(numerical["finite_difference"]) == 2
        and all(
            r["passed"] and r["predeclared_tolerance"] == 0.01
            for r in numerical["finite_difference"]
        ),
        "finite_difference_not_passed",
    )
    require(
        all(
            r["rigid_mode_rank"] == 6
            and not r["strict_full_positive_definiteness_required"]
            for r in numerical["hessians"].values()
        ),
        "incorrect_hessian_stability_policy",
    )
    fit = json.loads((root / "identifiability.json").read_text())
    require(
        not fit["fit_completed"]
        and not fit["fit_validated"]
        and fit["created_parameters"] == []
        and fit["discrepancy_counterexample"]["nullity"] > 0
        and fit["next_data"],
        "unidentified_coefficient_promoted",
    )
    return summary


def validate(root, *, repo, distribution, expected_sha256=None):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    if expected_sha256:
        require(pin(root / "receipt.json") == expected_sha256, "receipt_hash_mutation")
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["branch"] == BRANCH,
        "base_or_branch_mutation",
    )
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json"),
        "manifest_hash_mutation",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "promoted_receipt_readiness",
    )
    historical = history(repo, BASE)
    require(
        manifest["historical_files"] == historical
        and all(historical_pin(repo / name) == sha for name, sha in historical.items()),
        "historical_bytes_mutation",
    )
    require(
        set(manifest["implementation_files"]) == set(IMPLEMENTATION),
        "implementation_inventory_mutation",
    )
    expected_artifacts = {
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if p.is_file() and p.name not in ("manifest.json", "receipt.json")
    }
    require(
        set(manifest["artifacts"]) == expected_artifacts, "evidence_inventory_mutation"
    )
    for base, files in (
        (repo, manifest["implementation_files"]),
        (root, manifest["artifacts"]),
    ):
        for name, sha in files.items():
            path = base / name
            require(
                path.resolve().is_relative_to(base) and pin(path) == sha,
                "file_hash_mutation",
                path=name,
            )
    summary = audit(root, repo=repo, distribution=distribution)
    require(
        receipt["summary_sha256"] == pin(root / "summary.json")
        and receipt["claims"] == summary,
        "receipt_claims_mutation",
    )
    counts, nodes = junit(root, SUITES)
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected) and len(nodes) == len(set(nodes)),
        "full_suite_membership_mutation",
    )
    require(receipt["full_suite_counts"] == counts, "full_suite_counts_mutation")
    checks = json.loads((root / "validation/checks.json").read_text())
    require(
        REQUIRED_CHECKS <= {r["name"] for r in checks}
        and all(r["returncode"] == 0 for r in checks),
        "required_check_failed_or_missing",
    )
    return {
        "status": "passed",
        "historical_files_preserved": len(historical),
        "full_suite_counts": counts,
        "QM_executed": True,
        "parameter_completeness": "incomplete",
        "receipt_sha256": pin(root / "receipt.json"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--distribution", type=Path, required=True)
    parser.add_argument("--receipt-sha256")
    args = parser.parse_args()
    print(
        json.dumps(
            validate(
                args.directory,
                repo=args.repo,
                distribution=args.distribution,
                expected_sha256=args.receipt_sha256,
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
