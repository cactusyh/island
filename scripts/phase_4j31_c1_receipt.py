"""Generate/validate additive C1 evidence without rewriting any J17-J31 asset."""

import argparse
import json
import shutil
import subprocess
from hashlib import sha256
from pathlib import Path

from island.analysis.network import FLAGS, _body, _digest, _require, _signed
from scripts.phase_4j31_contract import CONTRACT
from scripts.validate_phase_4j31_c1_evidence import REVIEWED, historical_hashes
from scripts.validate_phase_4j31_receipt import (
    CHECKS,
    SUITES,
    junit_counts,
    validate_receipt,
)

BRANCH = "codex/phase-4j31-c1-strict-network-acceptance"
SCHEMA = "island.phase-4j31-c1-evidence.v1"
IMPLEMENTATION = (
    "scripts/phase_4j31_contract.py",
    "scripts/validate_phase_4j31_network.py",
    "scripts/validate_phase_4j31_receipt.py",
    "scripts/validate_phase_4j31_c1_evidence.py",
    "scripts/phase_4j31_c1_receipt.py",
    "tests/test_phase_4j31_execution.py",
    "tests/test_phase_4j31_receipt.py",
    "tests/test_phase_4j31_c1_outcomes.py",
    "docs/phase_4j31_c1.md",
)


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def generate(directory, destination):
    destination.mkdir(parents=True, exist_ok=False)
    source = json.loads((directory / "validation.json").read_text())
    artifacts = destination / "artifacts"
    artifacts.mkdir()
    for path in directory.iterdir():
        if path.is_file() and path.suffix in {".log", ".xml", ".json"}:
            shutil.copyfile(path, artifacts / path.name)
    suites = {}
    for name in sorted(SUITES):
        row = source["suites"][name]
        suites[name] = {
            "command": row["command"],
            "tested_commit": row["tested_commit"],
            "exit_code": row["exit_code"],
            "counts": junit_counts(artifacts / (name + ".xml")),
            "junit": str(artifacts / (name + ".xml")),
            "log": str(artifacts / (name + ".log")),
        }
    checks = {
        name: {**row, "log": str(artifacts / (name + ".log"))}
        for name, row in source["checks"].items()
    }
    receipt = _signed(
        {
            "schema": SCHEMA,
            "contract": CONTRACT,
            "branch": BRANCH,
            "base_commit": REVIEWED,
            "validated_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            "historical_receipt_identity": "32a5688d61cd367658c62dbab094794e969fac4046dc256c540f110ec9f91907",
            "historical_evidence_hashes": historical_hashes(),
            "historical_manifest_sha256": digest(
                "docs/paper/reproducibility_manifest.json"
            ),
            "historical_attributes_sha256": digest(".gitattributes"),
            "implementation_files": {name: digest(name) for name in IMPLEMENTATION},
            "validation": {
                "suites": suites,
                "checks": checks,
                "additional_gates": source["additional_gates"],
            },
            "environments": source["environments"],
            "retained_revalidation": json.loads(
                (directory / "retained-evidence.json").read_text()
            ),
            "before_negative_controls": json.loads(
                (directory / "before-reproductions.json").read_text()
            ),
            "before_regression_counts": junit_counts(
                directory / "before-regression.xml"
            ),
            "after_regression_counts": junit_counts(directory / "after-regression.xml"),
            "scientific_scope": {
                "negative_controls": "injected/mutated software controls; not scientific evidence",
                "retained_revalidation_new_dynamics": False,
                "j31_integration_tests": "focused and full suites each execute six native four-step uninterrupted versus two-plus-two-step network fixtures; PCFF rejects before dynamics",
                "other_integration_tests": "unchanged J22/J23 and historical analysis short integration tests also run in compatibility/regression/full suites",
                "new_scientific_campaign": False,
            },
            "remaining_unmet_gates": [],
            "artifact_files": {
                str(path): digest(path) for path in sorted(artifacts.iterdir())
            },
            **FLAGS,
        }
    )
    (destination / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    )
    return receipt


def validate(path):
    r = json.loads(path.read_text())
    _require(
        r.get("schema") == SCHEMA and r.get("identity") == _digest(_body(r)),
        "C1 receipt checksum/schema mismatch",
    )
    _require(
        r.get("base_commit") == REVIEWED
        and r.get("branch") == BRANCH
        and r.get("contract") == CONTRACT,
        "C1 branch/base/contract mismatch",
    )
    _require(
        all(r.get(k) == v and type(r.get(k)) is type(v) for k, v in FLAGS.items()),
        "C1 readiness mismatch",
    )
    commit = r["validated_commit"]
    _require(
        len(commit) == 40 and all(c in "0123456789abcdef" for c in commit),
        "C1 implementation SHA must be full",
    )
    subprocess.run(["git", "merge-base", "--is-ancestor", REVIEWED, commit], check=True)
    _require(
        set(r["implementation_files"]) == set(IMPLEMENTATION),
        "C1 implementation inventory mismatch",
    )
    for name, checksum in r["implementation_files"].items():
        original = subprocess.check_output(["git", "show", f"{commit}:{name}"])
        _require(
            digest(name) == sha256(original).hexdigest() == checksum,
            "C1 implementation source mismatch: " + name,
        )
    _require(
        historical_hashes() == r["historical_evidence_hashes"],
        "C1 historical evidence mismatch",
    )
    _require(
        digest("docs/paper/reproducibility_manifest.json")
        == r["historical_manifest_sha256"]
        and digest(".gitattributes") == r["historical_attributes_sha256"],
        "C1 historical metadata mismatch",
    )
    for name, checksum in r["artifact_files"].items():
        _require(digest(name) == checksum, "C1 artifact hash mismatch: " + name)
    _require(
        set(r["validation"]["suites"]) == SUITES
        and set(r["validation"]["checks"]) == CHECKS,
        "C1 validation inventory mismatch",
    )
    for name, row in r["validation"]["suites"].items():
        actual = junit_counts(row["junit"])
        _require(
            actual == row["counts"]
            and row["exit_code"] == 0
            and row["tested_commit"] == commit
            and actual["errors"] == actual["failures"] == 0
            and actual["passed"] > 0,
            "C1 suite evidence mismatch: " + name,
        )
    _require(
        all(row["exit_code"] == 0 for row in r["validation"]["checks"].values()),
        "C1 failed checks",
    )
    _require(
        set(r["validation"]["additional_gates"])
        == {"retained_evidence", "historical_receipt"}
        and all(
            row["exit_code"] == 0 and row["tested_commit"] == commit
            for row in r["validation"]["additional_gates"].values()
        ),
        "C1 additional gates incomplete",
    )
    _require(
        r["before_regression_counts"]["failures"] == 2
        and r["before_regression_counts"]["tests"] == 2,
        "C1 before regression controls missing",
    )
    _require(
        r["after_regression_counts"]["errors"]
        == r["after_regression_counts"]["failures"]
        == 0
        and r["after_regression_counts"]["passed"] >= 2,
        "C1 after regression controls missing",
    )
    before = r["before_negative_controls"]
    _require(
        before["defect_a"]["exit_code"] == 0
        and before["defect_a"]["result"]["status"] == "blocked"
        and before["defect_b"]["structural_acceptance"] == "ACCEPTED",
        "C1 original defect reproduction missing",
    )
    revalidation = r["retained_revalidation"]
    _require(
        revalidation["historical_file_hashes"] == r["historical_evidence_hashes"]
        and revalidation["new_dynamics_executed"] is False,
        "C1 retained revalidation mismatch",
    )
    expected = {
        f"{family}/{kind}"
        for family in ("oplsaa", "gaff", "gaff2")
        for kind in ("branched", "crosslinked")
    }
    _require(
        set(revalidation["revalidated_cases"]) == expected
        and all(
            row["observables_match"] is True
            and row["same_process_relocation"]
            == row["fresh_process_reconstruction"]
            == "PASS"
            for row in revalidation["revalidated_cases"].values()
        ),
        "C1 retained reconstruction incomplete",
    )
    _require(r["remaining_unmet_gates"] == [], "C1 has unmet gates")
    validate_receipt(
        Path("docs/evidence/phase_4j31.json"),
        Path("docs/paper/reproducibility_manifest.json"),
        expected_identity=r["historical_receipt_identity"],
    )
    return r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "receipt",
        type=Path,
        nargs="?",
        default=Path("docs/evidence/phase_4j31_c1/receipt.json"),
    )
    parser.add_argument("--generate", type=Path, metavar="VALIDATION_DIRECTORY")
    args = parser.parse_args()
    if args.generate:
        r = generate(args.generate, args.receipt.parent)
        print("C1 receipt generated:", r["identity"])
    else:
        validate(args.receipt)
        print("C1 receipt validation: PASS")


if __name__ == "__main__":
    main()
