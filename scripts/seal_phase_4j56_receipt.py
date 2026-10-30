"""Seal the exact J56 plan, validation membership and unchanged baseline files."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.pcff_completion import BASE, BRANCH, pin, require, signed, write
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j56_receipt import IMPLEMENTATION, SUITES, read


def seal(directory, repo):
    root, repo = Path(directory), Path(repo)
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == v for n, v in old.items()),
        "historical_bytes_mutation",
    )
    summary = read(root / "audit/summary.json")
    full, nodes = junit(root, SUITES)
    require(
        Counter(nodes)
        == Counter(json.loads((root / "collection.json").read_text())["nodes"])
        and len(nodes) == len(set(nodes)),
        "suite_membership_mutation",
    )
    focused, _ = junit(root, ["validation/focused-final.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(c["returncode"] == 0 for c in checks), "required_check_failed")
    manifest = signed(
        {
            "schema": "island.j56-manifest.v1",
            "base_commit": BASE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
            "source_hashes": summary["source_hashes"],
            "artifacts": {
                str(p.relative_to(root)): pin(p)
                for p in sorted(root.rglob("*"))
                if p.is_file()
                and p not in (root / "receipt.json", root / "manifest.json")
            },
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j56-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_identity": summary["identity"],
            **{
                k: summary[k]
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
                    "query_classification_counts",
                    "directed_source_queries",
                    "missing_site_requests",
                    "historical_missing_site_requests",
                    "author_confirmed_site_requests",
                    "author_confirmation_status",
                    "primary_author_artifact_status",
                    "author_confirmed_queries",
                    "qm_fit_required_queries",
                )
            },
            "parameters_or_charges_changed": False,
            "historical_claims_rewritten": False,
            "historical_files_preserved": len(old),
            "focused_counts": focused,
            "J54_J55_regression_counts": regressions,
            "full_suite_counts": full,
            "junit_files": SUITES,
            "check_returncodes": {c["name"]: c["returncode"] for c in checks},
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "full_suite_counts": full,
        "historical_files": len(old),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    a = p.parse_args()
    print(json.dumps(seal(a.directory, a.repo), sort_keys=True))


if __name__ == "__main__":
    main()
