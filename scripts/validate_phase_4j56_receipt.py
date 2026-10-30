"""Validate J56 data absence, exact Class-II plan, test membership and history."""

import argparse
import json
import tempfile
from collections import Counter
from pathlib import Path

from island.pcff_completion import BASE, BRANCH, digest, pin, require
from scripts.phase_4j56 import FIXTURES, audit
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

IMPLEMENTATION = (
    "src/island/pcff_completion/__init__.py",
    "src/island/pcff_completion/__main__.py",
    "src/island/pcff_completion/schema.py",
    "scripts/phase_4j56.py",
    "scripts/seal_phase_4j56_receipt.py",
    "scripts/validate_phase_4j56_receipt.py",
    "tests/test_phase_4j56.py",
    "docs/phase_4j56.md",
)
SUITES = [f"validation/full-{i}.xml" for i in range(8)] + [
    "validation/focused-final.xml"
]


def read(path):
    p = json.loads(Path(path).read_text())
    require(
        p.get("identity") == digest({k: v for k, v in p.items() if k != "identity"}),
        "record_checksum_mutation",
        path=str(path),
    )
    return p


def validate_plan(root, *, repo, distribution):
    root, repo = Path(root), Path(repo)
    summary = read(root / "summary.json")
    require(
        summary["base_commit"] == BASE and summary["branch"] == BRANCH,
        "base_or_branch_mutation",
    )
    require(
        summary["reference_dataset_status"] == "not_supplied"
        and summary["registered_profiles"] == summary["created_terms"] == [],
        "invented_reference_or_profile",
    )
    require(
        summary["directed_source_queries"] == 47
        and summary["missing_site_requests"] == 2916
        and summary["query_classification_counts"]
        == {"requires_author_confirmation": 28, "requires_qm_fit": 19},
        "classification_count_mutation",
    )
    for k in (
        "production_validated",
        "MD_executed",
        "NPT_executed",
        "historical_source_fitting_assumptions_changed",
    ):
        require(summary[k] is False, "authority_mutation", field=k)
    require(
        summary["parameter_completeness"] == "incomplete"
        and summary["cross_engine_equivalence"] == "not_claimed"
        and summary["simulation_readiness"]
        == summary["scientific_validation"]
        == summary["production_readiness"]
        == "not_established",
        "readiness_mutation",
    )
    with tempfile.TemporaryDirectory(prefix="island-j56-plan-check-") as tmp:
        expected_root = Path(tmp) / "audit"
        expected = audit(repo, distribution, expected_root, fresh=False)
        actual_files = {p.name for p in root.iterdir() if p.is_file()}
        require(
            actual_files == {p.name for p in expected_root.iterdir()},
            "plan_membership_mutation",
        )
        for p in expected_root.iterdir():
            if p.name != "summary.json":
                require(
                    json.loads(p.read_text())
                    == json.loads((root / p.name).read_text()),
                    "plan_evidence_mutation",
                    file=p.name,
                )
        actual = json.loads(json.dumps(summary))
        actual.pop("identity")
        expected.pop("identity")
        require(
            [c["fixture"] for c in summary["fixtures"]] == list(FIXTURES),
            "fixture_membership_mutation",
        )
        for c in actual["fixtures"]:
            require(
                c["fresh_processes"]
                == [
                    {
                        "hash_seed": seed,
                        "returncode": 0,
                        "record_identity": c["classification_identity"],
                    }
                    for seed in (1, 99397)
                ],
                "fresh_process_evidence_mutation",
            )
            c["fresh_processes"] = []
        require(actual == expected, "summary_evidence_mutation")
    return summary


def validate(directory, *, repo, distribution, expected_sha256=None):
    root, repo = Path(directory), Path(repo)
    if expected_sha256 is not None:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "independent_receipt_hash_mutation",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "base_or_branch_mutation",
    )
    require(
        r["schema"] == "island.j56-receipt.v1"
        and m["schema"] == "island.j56-manifest.v1",
        "schema_mutation",
    )
    require(
        r["manifest_sha256"] == pin(root / "manifest.json"), "manifest_hash_mutation"
    )
    old = history(repo, BASE)
    require(
        m["historical_files"] == old and len(old) == r["historical_files_preserved"],
        "historical_membership_mutation",
    )
    require(
        all(historical_pin(repo / n) == v for n, v in old.items()),
        "historical_bytes_mutation",
    )
    require(
        set(m["implementation_files"]) == set(IMPLEMENTATION)
        and all(pin(repo / n) == v for n, v in m["implementation_files"].items()),
        "implementation_mutation",
    )
    artifacts = {
        str(p.relative_to(root)): pin(p)
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "receipt.json", root / "manifest.json")
    }
    require(artifacts == m["artifacts"], "artifact_mutation")
    summary = validate_plan(root / "audit", repo=repo, distribution=distribution)
    require(
        summary["identity"] == r["summary_identity"]
        and m["source_hashes"] == summary["source_hashes"],
        "source_or_summary_mutation",
    )
    counts, nodes = junit(root, SUITES)
    collected = json.loads((root / "collection.json").read_text())
    require(
        Counter(nodes) == Counter(collected["nodes"]) and len(nodes) == len(set(nodes)),
        "incomplete_or_repeated_suite_membership",
    )
    require(counts == r["full_suite_counts"], "full_suite_count_mutation")
    focused, _ = junit(root, ["validation/focused-final.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    require(
        r["focused_counts"] == focused
        and r["J54_J55_regression_counts"] == regressions,
        "focused_or_regression_count_mutation",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    required = {
        "ruff",
        "compileall",
        "diff",
        "historical_hashes",
        "J54_receipt",
        "J55_receipt",
    } | {
        "pip-" + n for n in ("phase4e2", "phase4g1", "phase4j32", "phase4j36", "md_rg")
    }
    require(
        required <= {c["name"] for c in checks}
        and all(c["returncode"] == 0 for c in checks)
        and r["check_returncodes"] == {c["name"]: c["returncode"] for c in checks},
        "required_check_failed_or_absent",
    )
    for k in (
        "parameter_completeness",
        "cross_engine_equivalence",
        "scientific_validation",
        "production_readiness",
        "production_validated",
        "simulation_readiness",
        "registered_profiles",
        "created_terms",
        "MD_executed",
        "NPT_executed",
        "reference_dataset_status",
    ):
        require(r[k] == summary[k], "receipt_authority_mutation", field=k)
    require(
        r["parameters_or_charges_changed"] is False
        and r["historical_claims_rewritten"] is False,
        "historical_authority_mutation",
    )
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(old),
        "full_suite_counts": counts,
        "parameter_completeness": "incomplete",
        "registered_profiles": [],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--distribution", type=Path, required=True)
    p.add_argument("--expected-sha256")
    a = p.parse_args()
    print(
        json.dumps(
            validate(
                a.directory,
                repo=a.repo,
                distribution=a.distribution,
                expected_sha256=a.expected_sha256,
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
