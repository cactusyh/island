"""Verify the v2 evidence export, source pins and all preserved historical bytes."""

import argparse
import json
import lzma
from pathlib import Path

from island.charge_references.records import pack
from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_conventions import (
    DATA_ID,
    OPLSPDMSConventionsBlockedError,
    load_opls_pdms_conventions,
    source_evidence,
)
from island.forcefields.oplsaa.pdms_literature import build_opls_pdms_literature
from island.forcefields.pdms_literature import pinned_bytes
from island.workflows.storage import child
from scripts.assess_opls_pdms_conventions import fixtures, summary


def validate(directory, *, repo, references):
    root, repo = Path(directory), Path(repo)
    manifest = json.loads((root / "manifest.json").read_text())
    require(
        manifest["identity"]
        == digest({k: v for k, v in manifest.items() if k != "identity"}),
        "changed_manifest",
    )
    require(
        manifest["schema"] == "island_opls_pdms_conventions_manifest_v2",
        "wrong_manifest",
    )
    require(manifest["data_identity"] == DATA_ID, "wrong_data_identity")
    require(
        manifest["execution_profiles"] == [] and manifest["acceptance_gate"] == 1,
        "forged_authority",
    )
    for name, expected in manifest["artifacts"].items():
        pinned_bytes(child(repo, name), expected)
    baseline = json.loads(
        lzma.decompress((root / "protected-file-hashes.json.xz").read_bytes())
    )
    for name, expected in baseline.items():
        pinned_bytes(child(repo, name), expected)
    receipt = json.loads((root / "receipt.json").read_text())
    require(
        receipt["identity"]
        == digest({k: v for k, v in receipt.items() if k != "identity"}),
        "changed_receipt",
    )
    require(
        receipt["schema"] == "island_opls_pdms_conventions_receipt_v2"
        and receipt["data_identity"] == DATA_ID
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["profile_identity"] is None
        and receipt["limitations"] == LIMITS
        and receipt["acceptance_gate"] == 1
        and receipt["execution_profiles"] == receipt["backend_executions"] == [],
        "forged_receipt_authority",
    )
    supplement = build_opls_pdms_literature(
        DomainSelection(
            (repo / "docs/evidence/phase_4j45/selection.json").read_text(),
            (repo / "docs/evidence/phase_4j44/source-audit.json").read_text(),
        )
    )
    evidence = json.loads((root / "source-evidence.json").read_text())
    require(
        pack(evidence) == pack(source_evidence(supplement)), "changed_source_evidence"
    )
    report = json.loads((root / "report.json").read_text())
    expected = {}
    for name, system, oxygen in fixtures():
        case = report["cases"][name]
        result = load_opls_pdms_conventions(
            root / (name + ".json"),
            case["file_sha256"],
            supplement=supplement,
            references=references,
        )
        result.validate_integrity(system, references=references)
        require(
            pack(result.payload["requested_E_oxygen_site_id"]) == pack(oxygen),
            "changed_fixture_orientation",
        )
        expected[name] = {**summary(result), "file_sha256": case["file_sha256"]}
        try:
            result.require_executable()
        except OPLSPDMSConventionsBlockedError as error:
            require(error.result is result, "lost_blocked_result")
        else:
            raise ValueError("blocked_record_reached_execution")
    require(
        report
        == {
            "cases": expected,
            "acceptance_gate": 1,
            "execution_profiles": [],
            "backend_executions": [],
        },
        "changed_report",
    )
    return {
        "manifest_valid": True,
        "receipt_valid": True,
        "cases_reconstructed": len(expected),
        "historical_files_unchanged": len(baseline),
        "execution_profiles": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(validate(args.directory, repo=args.repo, references=args.references))
    )
