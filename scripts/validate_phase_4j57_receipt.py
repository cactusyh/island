"""Validate J57 immutable preparation evidence, suite membership and J17-J56 bytes."""

import argparse
import json
import tempfile
from collections import Counter
from pathlib import Path

from island.pcff_completion import digest, pin, require
from island.pcff_qm.__main__ import load_cases
from island.pcff_qm.records import BASE, BRANCH, PREREQUISITE
from scripts.phase_4j57 import audit, reconstruction
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

IMPLEMENTATION = tuple(
    f"src/island/pcff_qm/{name}.py"
    for name in (
        "__init__",
        "__main__",
        "records",
        "targets",
        "geometry",
        "inputs",
        "reference",
        "readiness",
        "schema",
        "supporting",
    )
) + (
    "scripts/phase_4j57.py",
    "scripts/seal_phase_4j57_receipt.py",
    "scripts/validate_phase_4j57_receipt.py",
    "tests/test_phase_4j57.py",
    "docs/phase_4j57.md",
)
SUITES = [f"validation/full-{i}.xml" for i in range(8)] + [
    "validation/focused-final.xml"
]
REQUIRED_CHECKS = {
    "ruff",
    "compileall",
    "diff",
    "historical_hashes",
    "J54_receipt",
    "J55_receipt",
    "J56_receipt",
} | {
    "pip-" + name
    for name in ("phase4e2", "phase4g1", "phase4j32", "phase4j36", "md_rg")
}


def read(path):
    p = json.loads(Path(path).read_text())
    require(
        p.get("identity") == digest({k: v for k, v in p.items() if k != "identity"}),
        "record_checksum_mutation",
        path=str(path),
    )
    return p


def validate_audit(root, *, repo, distribution):
    root = Path(root)
    s = read(root / "summary.json")
    require(
        s["base_commit"] == BASE
        and s["prerequisite_snapshot"] == PREREQUISITE
        and s["branch"] == BRANCH
        and not s["base_contains_merged_J56"],
        "base_or_prerequisite_mutation",
    )
    require(
        s["author_confirmed_queries"] == 28
        and s["qm_fit_required_queries"] == 19
        and s["missing_numerical_site_requests"] == 1068
        and s["author_attested_site_requests"] == 1848
        and s["historical_missing_site_requests"] == 2916,
        "query_inventory_mutation",
    )
    require(
        s["registered_profiles"] == s["created_terms"] == []
        and s["reference_dataset_status"] == "not_supplied",
        "invented_parameters_or_data",
    )
    for key in ("QM_executed", "MD_executed", "NPT_executed", "production_validated"):
        require(s[key] is False, "execution_authority_mutation", field=key)
    require(
        s["parameter_completeness"] == "incomplete"
        and s["cross_engine_equivalence"] == "not_claimed"
        and s["scientific_validation"]
        == s["production_readiness"]
        == s["simulation_readiness"]
        == "not_established",
        "readiness_mutation",
    )
    with tempfile.TemporaryDirectory(prefix="island-j57-reconstruct-") as tmp:
        expected_root = Path(tmp) / "audit"
        expected = audit(repo, distribution, expected_root, fresh=False)
        require(
            {p.name for p in root.iterdir()}
            == {p.name for p in expected_root.iterdir()},
            "audit_membership_mutation",
        )
        for file in expected_root.iterdir():
            if file.name != "summary.json":
                require(
                    file.read_bytes() == (root / file.name).read_bytes(),
                    "audit_reconstruction_mutation",
                    file=file.name,
                )
        actual = dict(s)
        actual.pop("identity")
        expected.pop("identity")
        geometry_checksums = reconstruction(load_cases(repo, distribution))[
            "fixture_geometry_checksums"
        ]
        require(
            actual["fresh_processes"]
            == [
                {
                    "hash_seed": seed,
                    "returncode": 0,
                    "inventory_identity": s["target_inventory_identity"],
                    "relocated_sources": True,
                    "all_five_fixture_graphs_reconstructed": True,
                    "fixture_geometry_checksums": geometry_checksums,
                }
                for seed in (1, 99397)
            ],
            "fresh_process_evidence_mutation",
        )
        actual["fresh_processes"] = []
        require(actual == expected, "summary_reconstruction_mutation")
    return s


def validate(directory, *, repo, distribution, expected_sha256=None):
    root, repo = Path(directory), Path(repo)
    if expected_sha256:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "independent_receipt_hash_mutation",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["schema"] == "island.j57-receipt.v1"
        and m["schema"] == "island.j57-manifest.v1",
        "schema_mutation",
    )
    require(
        r["manifest_sha256"] == pin(root / "manifest.json"), "manifest_hash_mutation"
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE
        and r["branch"] == BRANCH
        and r["prerequisite_snapshot"] == m["prerequisite_snapshot"] == PREREQUISITE,
        "base_or_branch_mutation",
    )
    old = history(repo, PREREQUISITE)
    require(
        m["historical_files"] == old and r["historical_files_preserved"] == len(old),
        "historical_membership_mutation",
    )
    require(
        all(historical_pin(repo / n) == sha for n, sha in old.items()),
        "historical_bytes_mutation",
    )
    require(
        m["implementation_files"] == {n: pin(repo / n) for n in IMPLEMENTATION},
        "implementation_mutation",
    )
    files = {
        p.relative_to(root).as_posix(): pin(p)
        for p in root.rglob("*")
        if p.is_file() and p.name not in ("manifest.json", "receipt.json")
    }
    require(files == m["artifacts"], "evidence_inventory_mutation")
    s = validate_audit(root / "audit", repo=repo, distribution=distribution)
    require(
        s["identity"] == r["summary_identity"]
        and s["source_hashes"] == m["source_hashes"],
        "summary_or_source_mutation",
    )
    counts, nodes = junit(root, SUITES)
    collected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(collected) and len(nodes) == len(set(nodes)),
        "suite_membership_mutation",
    )
    require(counts == r["full_suite_counts"], "suite_count_mutation")
    focused, _ = junit(root, ["validation/focused-final.xml"])
    regression, _ = junit(root, ["validation/regressions.xml"])
    require(
        focused == r["focused_counts"] and regression == r["J54_J56_regression_counts"],
        "focused_counts_mutation",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(
        REQUIRED_CHECKS <= {c["name"] for c in checks}
        and all(c["returncode"] == 0 for c in checks)
        and r["check_returncodes"] == {c["name"]: c["returncode"] for c in checks},
        "required_check_failed_or_absent",
    )
    for key in (
        "author_confirmation_status",
        "author_confirmed_queries",
        "primary_author_artifact_status",
        "qm_fit_required_queries",
        "parameter_completeness",
        "cross_engine_equivalence",
        "scientific_validation",
        "production_readiness",
        "production_validated",
        "simulation_readiness",
        "QM_executed",
        "MD_executed",
        "NPT_executed",
        "registered_profiles",
        "created_terms",
        "reference_dataset_status",
    ):
        require(r[key] == s[key], "receipt_claim_mutation", field=key)
    require(
        r["historical_claims_rewritten"] is False
        and r["parameters_or_charges_changed"] is False,
        "historical_claim_or_parameter_mutation",
    )
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(old),
        "full_suite_counts": counts,
        "parameter_completeness": "incomplete",
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
