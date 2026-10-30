"""Seal J58 after actual numerical acceptance and complete software verification."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j58_receipt import (
    BASE,
    BRANCH,
    IMPLEMENTATION,
    REQUIRED_CHECKS,
    SUITES,
    audit,
)


def seal(root, repo, distribution):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    summary = audit(root, repo=repo, distribution=distribution)
    historical = history(repo, BASE)
    require(
        all(historical_pin(repo / name) == sha for name, sha in historical.items()),
        "historical_bytes_mutation",
    )
    counts, nodes = junit(root, SUITES)
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected) and len(nodes) == len(set(nodes)),
        "full_suite_membership_mutation",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(
        REQUIRED_CHECKS <= {r["name"] for r in checks}
        and all(r["returncode"] == 0 for r in checks),
        "required_check_failed_or_missing",
    )
    manifest = signed(
        {
            "schema": "island.j58-manifest.v1",
            "base_commit": BASE,
            "historical_files": historical,
            "implementation_files": {name: pin(repo / name) for name in IMPLEMENTATION},
            "artifacts": {
                p.relative_to(root).as_posix(): pin(p)
                for p in sorted(root.rglob("*"))
                if p.is_file() and p.name not in ("manifest.json", "receipt.json")
            },
        }
    )
    write(root / "manifest.json", manifest)
    focused, _ = junit(root, ["validation/focused-final.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    receipt = signed(
        {
            "schema": "island.j58-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "prerequisites": json.loads((root / "plan/ancestry.json").read_text()),
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_sha256": pin(root / "summary.json"),
            "claims": summary,
            "historical_files_preserved": len(historical),
            "full_suite_counts": counts,
            "focused_counts": focused,
            "J54_J57_regression_counts": regressions,
            "check_returncodes": {r["name"]: r["returncode"] for r in checks},
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "full_suite_counts": counts,
        "historical_files_preserved": len(historical),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--distribution", type=Path, required=True)
    args = p.parse_args()
    print(
        json.dumps(seal(args.directory, args.repo, args.distribution), sort_keys=True)
    )


if __name__ == "__main__":
    main()
