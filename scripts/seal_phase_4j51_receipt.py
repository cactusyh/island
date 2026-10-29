"""Seal a completed offline J51 audit and its independently collected checks."""

import argparse
import json
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.pdms_opls import PROFILE
from island.opls_conventions import AUTHOR_ENGINE, EXACT_ENGINE, signed
from island.opls_resolution import REMAINING_SOURCE_GAPS, TEN_CONVENTIONS
from scripts.phase_4j47 import write
from scripts.phase_4j51 import BASE, BRANCH
from scripts.validate_phase_4j47_receipt import history, pin, require

IMPLEMENTATION = (
    "src/island/opls_resolution/__init__.py",
    "src/island/opls_resolution/backend_evidence.py",
    "src/island/opls_resolved_execution.py",
    "scripts/phase_4j51.py",
    "scripts/audit_phase_4j51_sources.py",
    "scripts/validate_phase_4j51_receipt.py",
    "scripts/seal_phase_4j51_receipt.py",
    "tests/test_phase_4j51.py",
    "docs/phase_4j51.md",
)


def junit(directory, files):
    counts, nodes = Counter(), []
    for file in files:
        for case in ET.parse(directory / file).iter("testcase"):
            classname = case.attrib["classname"]
            require(classname.startswith("tests."), "Wrong test class")
            nodes.append(classname.replace(".", "/") + ".py::" + case.attrib["name"])
            counts["tests"] += 1
            for key in ("failure", "error", "skipped"):
                counts[key] += int(case.find(key) is not None)
    require(counts["failure"] == counts["error"] == 0, "Tests did not pass")
    return dict(counts), nodes


def seal(directory, repo):
    root, repo = Path(directory).resolve(), Path(repo).resolve()
    summary = json.loads((root / "audit/summary.json").read_text())
    independent = json.loads((root / "independent-reconstruction.json").read_text())
    require(
        independent["outside_production_parser"]
        and independent["numerical"]["raw_FD_samples_reconstructed"] == 288
        and independent["no_tail_rule_reconstructed"] == 8,
        "Independent reconstruction did not complete",
    )
    require(
        all(
            v["relocation"]["verified"]
            and all(p["returncode"] == 0 for p in v["relocation"]["fresh_processes"])
            for case in summary["cases"].values()
            for v in case["policies"].values()
        ),
        "Fresh process relocation did not complete",
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
        "Incomplete or repeated full-suite accounting",
    )
    focused, _ = junit(root, ["validation/focused-final.xml"])
    regressions, _ = junit(root, ["validation/j49-j50-regression.xml"])
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(v["returncode"] == 0 for v in checks), "A required check failed")
    m = signed(
        {
            "schema": "island.j51-manifest.v1",
            "base_commit": BASE,
            "historical_files": old,
            "implementation_files": {n: pin(repo / n) for n in IMPLEMENTATION},
            "profile_identities": [PROFILE],
            "source_assets": {
                **{
                    r["source_file"]: r["source_sha256"]
                    for r in independent["source_reviews"]
                },
                "PDMS_benchmark-main.zip": summary["source_archive_sha256"],
            },
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
            "schema": "island.j51-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "profile_identity": PROFILE,
            "manifest_sha256": pin(root / "manifest.json"),
            "summary_identity": summary["identity"],
            "ten_conventions": list(TEN_CONVENTIONS),
            "converter_resolved_from_ten": ["lj_tail_correction"],
            "remaining_source_gaps": list(REMAINING_SOURCE_GAPS),
            "diagnostic_acceptance_gate": summary["diagnostic_acceptance_gate"],
            "research_opt_in_acceptance_gate": 1,
            "production_candidate_acceptance_gate": 1,
            "MD_executed": False,
            "NPT_executed": False,
            "charges_parameters_or_tolerances_changed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "focused_tests_passed": focused["tests"] - focused["skipped"],
            "J49_J50_regressions_passed": regressions["tests"] - regressions["skipped"],
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
