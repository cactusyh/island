"""Seal completed J52 diagnostics, independent reconstruction and exact test union."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.forcefields.pdms_opls import PROFILE
from island.opls_conventions import AUTHOR_ENGINE, EXACT_ENGINE, signed
from island.opls_final_resolution import NINE_CONVENTIONS, corpus_audit
from scripts.phase_4j47 import write
from scripts.phase_4j52 import BASE, BRANCH
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history, pin, require

IMPLEMENTATION = (
    "src/island/opls_final_resolution/__init__.py",
    "src/island/opls_final_resolution/evidence.py",
    "src/island/opls_final_execution.py",
    "scripts/phase_4j52.py",
    "scripts/audit_phase_4j52_sources.py",
    "scripts/reconstruct_phase_4j52_evidence.py",
    "scripts/seal_phase_4j52_receipt.py",
    "scripts/validate_phase_4j52_receipt.py",
    "scripts/check_phase_4j52_receipt.py",
    "tests/test_phase_4j52.py",
    "docs/phase_4j52.md",
)


def seal(directory, repo):
    root, repo = Path(directory).resolve(), Path(repo).resolve()
    summary = json.loads((root / "audit/summary.json").read_text())
    independent = json.loads((root / "independent-reconstruction.json").read_text())
    require(
        independent["outside_production_parser"]
        and independent["numerical"]["raw_FD_samples_reconstructed"] == 288
        and independent["assets_independently_reconstructed"] == 64,
        "Independent J52 reconstruction incomplete",
    )
    old = history(repo, BASE)
    require(all(pin(repo / n) == v for n, v in old.items()), "Historical bytes changed")
    files = [f"validation/full-{i}.xml" for i in range(8)] + [
        "validation/focused-final.xml"
    ]
    counts, nodes = junit(root, files)
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected)
        and all(n == 1 for n in Counter(nodes).values()),
        "Incomplete or duplicated ordinary full suite",
    )
    focused, _ = junit(root, ["validation/focused-final.xml"])
    regressions, _ = junit(root, ["validation/j49-j51-regression.xml"])
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(v["returncode"] == 0 for v in checks), "Required check failed")
    audit = corpus_audit()
    m = signed(
        {
            "schema": "island.j52-manifest.v1",
            "base_commit": BASE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
            "profile_identities": [PROFILE],
            "source_audit_identity": audit["identity"],
            "source_assets": {
                a["source_file"]: a["source_sha256"] for a in audit["assets"]
            },
            "source_archive_sha256": audit["source_archive_sha256"],
            "engine_executable_sha256": {
                "source_replay_v1": AUTHOR_ENGINE,
                "exact_erfc_diagnostic_v1": EXACT_ENGINE,
            },
            "artifacts": {
                str(p.relative_to(root)): pin(p)
                for p in sorted(root.rglob("*"))
                if p.is_file()
                and p not in (root / "receipt.json", root / "manifest.json")
            },
        }
    )
    write(root / "manifest.json", m)
    r = signed(
        {
            "schema": "island.j52-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "profile_identity": PROFILE,
            "public_family": "opls",
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_identity": summary["identity"],
            "source_audit_identity": audit["identity"],
            "nine_conventions": list(NINE_CONVENTIONS),
            "permanently_unresolved_in_pinned_corpus": list(NINE_CONVENTIONS),
            "diagnostic_acceptance_gate": summary["diagnostic_acceptance_gate"],
            "research_opt_in_acceptance_gate": 1,
            "production_candidate_acceptance_gate": 1,
            "MD_executed": False,
            "NPT_executed": False,
            "checkpoint_RNG_continuation": "not run: source prerequisites remain blocked",
            "charges_parameters_or_tolerances_changed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "focused_tests_passed": focused["tests"] - focused["skipped"],
            "J49_J51_regressions_passed": regressions["tests"] - regressions["skipped"],
            "full_suite_counts": counts,
            "junit_files": files,
            "check_returncodes": {
                **{v["name"]: v["returncode"] for v in checks},
                "historical_hashes": 0,
                "independent_reconstruction": 0,
                "fresh_process_relocation": 0,
            },
        }
    )
    write(root / "receipt.json", r)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "full_suite_counts": counts,
        "historical_files": len(old),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    a = p.parse_args(argv)
    print(json.dumps(seal(a.directory, a.repo), sort_keys=True))


if __name__ == "__main__":
    main()
