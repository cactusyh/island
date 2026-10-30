"""Seal J53 intake evidence without promoting a source or a scientific claim."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.forcefields.pdms_opls import PROFILE
from island.opls_conventions import signed
from island.opls_final_resolution import NINE_CONVENTIONS, corpus_audit
from scripts.phase_4j47 import write
from scripts.phase_4j53 import BASE, BRANCH
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history, pin, require

IMPLEMENTATION = (
    "src/island/opls_source_intake/__init__.py",
    "src/island/opls_source_intake/schema.py",
    "src/island/opls_source_execution.py",
    "src/island/intake_opls_source.py",
    "scripts/phase_4j53.py",
    "scripts/seal_phase_4j53_receipt.py",
    "scripts/validate_phase_4j53_receipt.py",
    "scripts/check_phase_4j53_receipt.py",
    "tests/test_phase_4j53.py",
    "docs/phase_4j53.md",
)


def seal(directory, repo):
    root, repo = Path(directory).resolve(), Path(repo).resolve()
    summary = json.loads((root / "audit/summary.json").read_text())
    require(
        summary["source_status"] == "not_supplied"
        and not summary["new_content_candidates"]
        and not summary["MD_executed"],
        "Unexpected new source or execution",
    )
    old = history(repo, BASE)
    require(all(pin(repo / n) == v for n, v in old.items()), "Historical bytes changed")
    files = [f"validation/full-{i}.xml" for i in range(8)] + ["validation/focused.xml"]
    counts, nodes = junit(root, files)
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected)
        and all(n == 1 for n in Counter(nodes).values()),
        "Incomplete full-suite union",
    )
    focused, _ = junit(root, ["validation/focused.xml"])
    regression, _ = junit(root, ["validation/j49-j52-regression.xml"])
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(c["returncode"] == 0 for c in checks), "Required check failed")
    audit = corpus_audit()
    manifest = signed(
        {
            "schema": "island.j53-manifest.v1",
            "base_commit": BASE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
            "profile_identities": [PROFILE],
            "new_profiles": [],
            "J52_source_audit_identity": audit["identity"],
            "source_hashes": {
                r["path"]: r["sha256"]
                for r in json.loads((root / "audit/source-delta.json").read_text())[
                    "files"
                ]
            },
            "artifacts": {
                str(p.relative_to(root)): pin(p)
                for p in sorted(root.rglob("*"))
                if p.is_file()
                and p not in (root / "manifest.json", root / "receipt.json")
            },
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j53-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            "public_family": "opls",
            "profile_identity": PROFILE,
            "source_status": "not_supplied",
            "new_profiles": [],
            "permanent_blockers": list(NINE_CONVENTIONS),
            "summary_identity": summary["identity"],
            "J52_source_audit_identity": audit["identity"],
            "new_source_executable_acceptance_gate": 1,
            "research_opt_in_acceptance_gate": 1,
            "production_candidate_acceptance_gate": 1,
            "MD_executed": False,
            "NPT_executed": False,
            "charges_parameters_or_tolerances_changed": False,
            "new_static_comparisons_executed": False,
            "checkpoint_RNG_continuation": summary["checkpoint_RNG_continuation"],
            "production_validated": False,
            "simulation_readiness": "not_established",
            "focused_tests_passed": focused["tests"] - focused["skipped"],
            "J49_J52_regressions_passed": regression["tests"] - regression["skipped"],
            "full_suite_counts": counts,
            "junit_files": files,
            "check_returncodes": {c["name"]: c["returncode"] for c in checks},
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(old),
        "full_suite_counts": counts,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    args = p.parse_args(argv)
    print(json.dumps(seal(args.directory, args.repo), sort_keys=True))


if __name__ == "__main__":
    main()
