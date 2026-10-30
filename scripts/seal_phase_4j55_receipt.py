"""Seal independent J55 parameter/backend claims and exact historical bytes."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.forcefields.opls_registry.registry import NATIVE
from island.forcefields.pdms_opls import PROFILE
from island.opls_conventions import signed
from scripts.phase_4j47 import write
from scripts.phase_4j55 import BASE, BRANCH
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history, pin, require

IMPLEMENTATION = (
    "src/island/forcefield_policy/__init__.py",
    "src/island/forcefield_policy/claims.py",
    "src/island/forcefield_policy/evidence.py",
    "src/island/inspect_parameters.py",
    "scripts/phase_4j55.py",
    "scripts/seal_phase_4j55_receipt.py",
    "scripts/validate_phase_4j55_receipt.py",
    "scripts/check_phase_4j55_receipt.py",
    "tests/test_phase_4j55.py",
    "docs/phase_4j55.md",
    "docs/parameter_completeness_policy.md",
    "docs/forcefield_readiness_matrix_j55.md",
)
SUITES = [f"validation/full-{i}.xml" for i in range(8)] + ["validation/focused.xml"]


def seal(root, repo):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    old = history(repo, BASE)
    require(
        all(pin(repo / n) == value for n, value in old.items()),
        "Historical bytes changed",
    )
    summary = json.loads((root / "audit/summary.json").read_text())
    counts, nodes = junit(root, SUITES)
    collection = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(collection) and len(nodes) == len(set(nodes)),
        "Incomplete suite membership",
    )
    focused, _ = junit(root, ["validation/focused.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(c["returncode"] == 0 for c in checks), "Required check failed")
    manifest = signed(
        {
            "schema": "island.j55-manifest.v1",
            "base_commit": BASE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
            "profile_identities": [NATIVE, PROFILE],
            "new_profiles": [],
            "source_hashes": summary["source_hashes"],
            "artifacts": {
                str(p.relative_to(root)): pin(p)
                for p in sorted(root.rglob("*"))
                if p.is_file()
                and p.name not in ("receipt.json",)
                and p not in (root / "manifest.json",)
            },
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j55-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_identity": summary["identity"],
            "profile_identities": [NATIVE, PROFILE],
            "new_profiles": [],
            "parameter_completeness_independent_of_equivalence": True,
            "ordinary_forcefield": "opls",
            "ordinary_profile_argument": None,
            "ordinary_charge_method_argument": None,
            "bounded_PDMS_parameter_completeness": "incomplete",
            "bounded_PDMS_missing_model_conventions": [
                "dielectric",
                "periodic_image_policy",
            ],
            "native_finite_control_parameter_completeness": "complete",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "MD_executed": False,
            "NPT_executed": False,
            "charges_parameters_or_tolerances_changed": False,
            "historical_claims_rewritten": False,
            "historical_files_preserved": len(old),
            "junit_files": SUITES,
            "focused_tests_passed": focused["tests"] - focused["skipped"],
            "J48_J53_regression_counts": regressions,
            "J54": "absent_from_exact_base",
            "full_suite_counts": counts,
            "check_returncodes": {c["name"]: c["returncode"] for c in checks},
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "counts": counts,
        "historical_files": len(old),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    args = p.parse_args(argv)
    print(json.dumps(seal(args.directory, args.repo), sort_keys=True))


if __name__ == "__main__":
    main()
