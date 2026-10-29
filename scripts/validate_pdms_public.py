"""Validate an independently pinned PDMS research receipt; no MD execution."""

import argparse
import json
import re
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.external_pdms import LIMITS, digest
from island.forcefields.pdms_opls import PROFILE, check
from scripts.validate_pdms_execution_boundary import validate as validate_historical


def validate_coverage(root):
    """Check sequential observed outcomes plus the four native JUnit shards."""
    trace = root / "suite-trace"
    nodes = json.loads((trace / "collection.json").read_text())
    observed = json.loads((trace / "observed-prefix.json").read_text())
    schedule = json.loads((trace / "schedule.json").read_text())
    old = [
        n
        for n in nodes
        if "test_blocked_preparation_cannot_start_md" not in n
        and "test_workflow_resume_and_thermodynamics" not in n
    ]
    statuses = []
    for line in (trace / "original-progress.log").read_text().splitlines():
        match = re.fullmatch(r"([.s]+)(?:\s+\[\s*\d+%\])?", line.strip())
        if match:
            statuses.extend(match.group(1))
    check(
        len(old) == schedule["original_collected"]
        and len(statuses) == schedule["observed_status_characters"],
        "progress_inventory_changed",
    )
    reconstructed = {
        n: s
        for n, s in zip(old, statuses)
        if n.split("::")[0]
        not in (schedule["interrupted_module_rerun"], "tests/test_pdms_public.py")
    }
    check(reconstructed == observed, "sequential_progress_changed")

    def key(node):
        parts = node.split("::")
        return (
            parts[0][:-3].replace("/", ".") + "".join("." + p for p in parts[1:-1]),
            parts[-1],
        )

    expected_keys = {key(n) for n in nodes}
    full_cases = list(ET.parse(root / "full-suite.xml").getroot().iter("testcase"))
    full = {(c.get("classname"), c.get("name")): c for c in full_cases}
    check(
        len(full) == len(full_cases) == len(expected_keys)
        and set(full) == expected_keys,
        "full_collection_coverage_changed",
    )
    completed = set(observed)
    for i in range(4):
        check(
            json.loads((trace / f"exit-{i}.json").read_text())["returncode"] == 0,
            "complement_shard_failed",
        )
        selected = json.loads((trace / f"nodes-{i}.json").read_text())
        cases = list(ET.parse(trace / f"shard-{i}.xml").getroot().iter("testcase"))
        actual = {(c.get("classname"), c.get("name")) for c in cases}
        check(
            actual == {key(n) for n in selected} and not (completed & set(selected)),
            "shard_coverage_changed",
        )
        check(
            all(c.find("failure") is None and c.find("error") is None for c in cases),
            "complement_test_failure",
        )
        completed.update(selected)
    check(completed == set(nodes), "missing_full_suite_tests")
    for node, status in observed.items():
        check(
            (full[key(node)].find("skipped") is not None) == (status == "s"),
            "sequential_outcome_changed",
        )


def validate(directory, expected_sha256, repo):
    root, repo = Path(directory), Path(repo)
    raw = (root / "manifest.json").read_bytes()
    check(sha256(raw).hexdigest() == expected_sha256, "manifest_hash_mismatch")
    manifest = json.loads(raw)
    for name, pin in manifest["artifacts"].items():
        relative = Path(name)
        check(
            not relative.is_absolute() and ".." not in relative.parts,
            "invalid_artifact_path",
        )
        path = (
            root / name[len("docs/evidence/pdms_public_v1/") :]
            if name.startswith("docs/evidence/pdms_public_v1/")
            else repo / name
        )
        check(
            sha256(path.read_bytes()).hexdigest() == pin,
            "artifact_hash_mismatch",
            file=name,
        )
    receipt = json.loads((root / "receipt.json").read_text())
    check(
        receipt["identity"]
        == digest({k: v for k, v in receipt.items() if k != "identity"}),
        "receipt_identity_mismatch",
    )
    check(
        receipt["profile_identity"] == PROFILE and receipt["limitations"] == LIMITS,
        "profile_or_limits_changed",
    )
    check(
        receipt["software_acceptance_gate"] == 0
        and receipt["scientific_acceptance_gate"] == 1
        and receipt["native_pcff_acceptance_gate"] == 1,
        "forged_scientific_promotion",
    )
    check(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "readiness_changed",
    )
    qualification = json.loads((root / "qualification.json").read_text())
    check(
        qualification["passed"] and not qualification["scientifically_validated"],
        "numerical_audit_failed",
    )
    for case in qualification["cases"].values():
        check(
            case["passed"] and not case["lammps_exact_equation_energy_parity"],
            "raw_lammps_difference_hidden",
        )
    for filename in (
        "focused.xml",
        "physical.xml",
        "regressions.xml",
        "full-suite.xml",
    ):
        suites = ET.parse(root / filename).getroot()
        check(
            all(
                int(s.get("failures", 0)) == 0 and int(s.get("errors", 0)) == 0
                for s in suites.iter("testsuite")
            ),
            "test_failures",
            file=filename,
        )
    validate_coverage(root)
    historical = validate_historical(
        repo / "docs/evidence/opls_pdms_execution_boundary_v4",
        "d23c717a4d087a9903a6096f29157fd10225d8ccb9ebdfbc6da23984469fd89e",
        repo,
    )
    return {
        "valid": True,
        "profile_identity": PROFILE,
        "historical_files_unchanged": historical["historical_files_unchanged"],
        "software_acceptance_gate": 0,
        "scientific_acceptance_gate": 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--repo", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(validate(args.directory, args.sha256, args.repo)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
