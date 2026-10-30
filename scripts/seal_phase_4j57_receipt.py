"""Seal exact J57 input-preparation evidence only after completed checks."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm.records import BASE, BRANCH, PREREQUISITE
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j57_receipt import (
    IMPLEMENTATION,
    REQUIRED_CHECKS,
    SUITES,
    read,
)


def seal(directory, repo):
    root, repo = Path(directory), Path(repo)
    old = history(repo, PREREQUISITE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in old.items()),
        "historical_bytes_mutation",
    )
    s = read(root / "audit/summary.json")
    counts, nodes = junit(root, SUITES)
    collected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(collected) and len(nodes) == len(set(nodes)),
        "suite_membership_mutation",
    )
    focused, _ = junit(root, ["validation/focused-final.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    checks = json.loads((root / "validation/checks.json").read_text())
    require(
        REQUIRED_CHECKS <= {c["name"] for c in checks}
        and all(c["returncode"] == 0 for c in checks),
        "required_check_failed_or_absent",
    )
    m = signed(
        {
            "schema": "island.j57-manifest.v1",
            "base_commit": BASE,
            "prerequisite_snapshot": PREREQUISITE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
            "source_hashes": s["source_hashes"],
            "artifacts": {
                p.relative_to(root).as_posix(): pin(p)
                for p in sorted(root.rglob("*"))
                if p.is_file() and p.name not in ("manifest.json", "receipt.json")
            },
        }
    )
    write(root / "manifest.json", m)
    r = signed(
        {
            "schema": "island.j57-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "prerequisite_snapshot": PREREQUISITE,
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_identity": s["identity"],
            "historical_files_preserved": len(old),
            "historical_claims_rewritten": False,
            "parameters_or_charges_changed": False,
            "focused_counts": focused,
            "J54_J56_regression_counts": regressions,
            "full_suite_counts": counts,
            "check_returncodes": {c["name"]: c["returncode"] for c in checks},
            **{
                k: s[k]
                for k in (
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
                )
            },
        }
    )
    write(root / "receipt.json", r)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "full_suite_counts": counts,
        "historical_files_preserved": len(old),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    a = p.parse_args()
    print(json.dumps(seal(a.directory, a.repo), sort_keys=True))


if __name__ == "__main__":
    main()
