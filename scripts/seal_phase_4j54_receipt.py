"""Seal exact J54 source audit, test membership and preservation evidence."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.pcff_pdms import CANDIDATE
from scripts.phase_4j54 import BASE, BRANCH, signed, write
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history, pin, require
from scripts.validate_phase_4j54_receipt import IMPLEMENTATION, SUITES, read


def seal(root, repo):
    root, repo = Path(root), Path(repo)
    old = history(repo, BASE)
    require(all(pin(repo / n) == v for n, v in old.items()), "Historical bytes changed")
    summary = read(root / "audit/summary.json")
    counts, nodes = junit(root, SUITES)
    require(
        Counter(nodes)
        == Counter(json.loads((root / "collection.json").read_text())["nodes"])
        and len(nodes) == len(set(nodes)),
        "Full-suite membership changed",
    )
    focused, _ = junit(root, ["validation/focused.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(c["returncode"] == 0 for c in checks), "Required check failed")
    context = read(root / "audit/literature-context.json")
    manifest = signed(
        {
            "schema": "island.j54-manifest.v1",
            "base_commit": BASE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
            "source_hashes": summary["source_hashes"],
            "literature_sources": context["sources"],
            "candidate_profile": CANDIDATE,
            "registered_profiles": [],
            "artifacts": {
                str(p.relative_to(root)): pin(p)
                for p in sorted(root.rglob("*"))
                if p.is_file()
                and p not in (root / "manifest.json", root / "receipt.json")
            },
        }
    )
    write(root / "manifest.json", manifest)
    from scripts.select_phase_4j54_validation import selection

    receipt = signed(
        {
            "schema": "island.j54-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_identity": summary["identity"],
            "candidate_profile": CANDIDATE,
            "registered_profiles": [],
            "ordinary_forcefield": "pcff",
            "ordinary_profile_argument": None,
            "ordinary_charge_method_argument": None,
            "parameter_completeness": "incomplete",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "MD_executed": False,
            "NPT_executed": False,
            "parameters_charges_or_tolerances_changed": False,
            "historical_claims_rewritten": False,
            "historical_files_preserved": len(old),
            "focused_counts": focused,
            "J55_regression_counts": regressions,
            "full_suite_counts": counts,
            "ordinary_suite_retry": selection(root),
            "junit_files": SUITES,
            "check_returncodes": {c["name"]: c["returncode"] for c in checks},
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "full_suite_counts": counts,
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
