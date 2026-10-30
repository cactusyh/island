"""Seal or validate J59 with existing history, JUnit, manifest and native replay."""

import argparse
import json
import lzma
from collections import Counter
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET

from island.pcff_completion import pin, require, signed, write
from scripts.phase_4j59 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j58_receipt import read

BASE = "d361f7e0a7c16a3a3a3513987971db6127690386"
BRANCH = "codex/phase-4j59-pcff-residual-model-qm-convergence"
IMPLEMENTATION = tuple(
    "src/island/pcff_qm_pilot/" + n + ".py"
    for n in (
        "convergence",
        "convergence_fast_engine",
        "convergence_importer",
        "convergence_analysis",
        "residual",
        "residual_analysis",
    )
) + (
    "scripts/phase_4j59.py",
    "scripts/independent_phase_4j59.py",
    "scripts/validate_phase_4j59_receipt.py",
    "tests/test_phase_4j59.py",
    "docs/phase_4j59.md",
)
SUITES = [
    f"validation/full-{i}-final.xml" if i == 3 else f"validation/full-{i}.xml"
    for i in range(8)
] + ["validation/focused-final.xml"]
REQUIRED = (
    {
        "ruff",
        "compileall",
        "diff",
        "historical_hashes",
        "independent_units",
        "relocation_fresh_process",
    }
    | {f"J{j}_receipt" for j in range(54, 59)}
    | {"pip-" + n for n in ("phase4e2", "phase4g1", "phase4j32", "phase4j36", "md_rg")}
)


def verify_software_retry(root):
    """Keep the failed event; accept only its exact successful isolated rerun."""
    root = Path(root)
    original = lzma.decompress((root / "validation/full-3.xml.xz").read_bytes())
    retry_record = json.loads((root / "validation/software-retry.json").read_text())
    require(
        sha256(original).hexdigest() == retry_record["original_uncompressed_sha256"],
        "original_failed_JUnit_mutation",
    )
    raw = list(ET.fromstring(original).iter("testcase"))
    failed = [
        c for c in raw if c.find("failure") is not None or c.find("error") is not None
    ]
    retry = list(ET.parse(root / "validation/retry-j31-gaff.xml").iter("testcase"))
    expected = (
        "tests.test_phase_4j31_execution",
        "test_real_final_graph_source_completeness_and_exact_execution[crosslinked-gaff]",
    )
    require(
        len(failed) == len(retry) == 1
        and all(
            (c.attrib["classname"], c.attrib["name"]) == expected
            for c in failed + retry
        ),
        "unexpected_full_suite_retry",
    )
    require(
        retry[0].find("failure") is None
        and retry[0].find("error") is None
        and retry[0].find("skipped") is None,
        "isolated_retry_not_passed",
    )
    effective = list(ET.parse(root / "validation/full-3-final.xml").iter("testcase"))
    replacement = [retry[0] if c in failed else c for c in raw]
    require(
        [ET.tostring(c) for c in effective] == [ET.tostring(c) for c in replacement],
        "retry_aggregate_mutation",
    )
    return {
        "initial_failures": 1,
        "successful_isolated_retries": 1,
        "original_attempt_preserved": True,
    }


def audit(root, repo, distribution):
    root, repo = Path(root), Path(repo)
    verify_software_retry(root)
    require(
        replay(repo, root, distribution)
        == json.loads((root / "relocation-replay.json").read_text()),
        "native_source_replay_failed",
    )
    ident = read(root / "residual-identifiability.json")
    model = read(root / "interaction-table.json")
    qm = read(root / "QM-convergence.json")
    summary = json.loads((root / "summary.json").read_text())
    require(
        model["status_counts"] == {"fixed": 74, "unresolved": 115, "not_applicable": 8},
        "interaction_accounting_mutation",
    )
    require(
        ident["model_identity"] == model["identity"],
        "residual_source_model_binding_mutation",
    )
    require(
        ident["energy_design"]["rank"] == 22
        and ident["energy_force_design"]["rank"] == 31
        and not ident["physical_identifiability"],
        "conditional_rank_promoted_to_physical_authority",
    )
    require(
        not ident["fit_completed"]
        and not ident["fit_validated"]
        and ident["fitted_candidates"] == [],
        "unjustified_fit",
    )
    require(
        qm["final_QM_adequacy"] == "not_established" and not qm["fit_authorized"],
        "promoted_QM_adequacy",
    )
    require(
        summary["parameter_completeness"] == "incomplete"
        and summary["cross_engine_equivalence"] == "not_claimed"
        and summary["scientific_validation"]
        == summary["simulation_readiness"]
        == "not_established"
        and summary["production_validated"] is False
        and summary["new_profiles"] == []
        and not summary["PCFF_MD_executed"],
        "promoted_readiness",
    )
    require(
        summary["author_confirmed_audited"] == 28
        and summary["original_author_artifacts"] == "not_yet_supplied"
        and summary["remaining_QM_fit_required"] == 19,
        "invented_author_or_numerical_evidence",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(
        REQUIRED <= {r["name"] for r in checks}
        and all(r["returncode"] == 0 for r in checks),
        "missing_or_failed_check",
    )
    counts, nodes = junit(root, SUITES)
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected) and len(nodes) == len(set(nodes)),
        "full_suite_membership_mutation",
    )
    return summary, counts


def seal(root, repo, distribution):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    summary, counts = audit(root, repo, distribution)
    h = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == sha for n, sha in h.items()),
        "historical_bytes_mutation",
    )
    manifest = signed(
        {
            "schema": "island.j59-manifest.v1",
            "base_commit": BASE,
            "historical_files": h,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
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
            "schema": "island.j59-receipt.v1",
            "branch": BRANCH,
            "base_commit": BASE,
            "prerequisites": json.loads((root / "plan/ancestry.json").read_text()),
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_sha256": pin(root / "summary.json"),
            "claims": summary,
            "full_suite_counts": counts,
            "focused_counts": focused,
            "software_retry": verify_software_retry(root),
            "J54_J58_regression_counts": regressions,
            "historical_files_preserved": len(h),
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "status": "sealed",
        "receipt_sha256": pin(root / "receipt.json"),
        "full_suite_counts": counts,
    }


def validate(root, repo, distribution, expected_sha256=None):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    if expected_sha256:
        require(pin(root / "receipt.json") == expected_sha256, "receipt_pin_mutation")
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["branch"] == BRANCH,
        "base_branch_mutation",
    )
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json"),
        "manifest_hash_mutation",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "receipt_readiness_mutation",
    )
    h = history(repo, BASE)
    require(
        h == manifest["historical_files"]
        and all(historical_pin(repo / n) == sha for n, sha in h.items()),
        "historical_bytes_mutation",
    )
    artifacts = {
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if p.is_file() and p.name not in ("manifest.json", "receipt.json")
    }
    require(
        set(manifest["implementation_files"]) == set(IMPLEMENTATION)
        and set(manifest["artifacts"]) == artifacts,
        "recursive_inventory_mutation",
    )
    for base, files in (
        (repo, manifest["implementation_files"]),
        (root, manifest["artifacts"]),
    ):
        for n, sha in files.items():
            require(
                (base / n).resolve().is_relative_to(base) and pin(base / n) == sha,
                "artifact_hash_mutation",
                path=n,
            )
    summary, counts = audit(root, repo, distribution)
    require(
        summary == receipt["claims"]
        and counts == receipt["full_suite_counts"]
        and receipt["summary_sha256"] == pin(root / "summary.json"),
        "receipt_claims_mutation",
    )
    return {
        "status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files_preserved": len(h),
        "full_suite_counts": counts,
        "parameter_completeness": "incomplete",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--distribution", type=Path, required=True)
    p.add_argument("--seal", action="store_true")
    p.add_argument("--receipt-sha256")
    a = p.parse_args()
    print(
        json.dumps(
            seal(a.directory, a.repo, a.distribution)
            if a.seal
            else validate(a.directory, a.repo, a.distribution, a.receipt_sha256),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
