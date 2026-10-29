"""Offline validation of the independently pinned v3 author-audit manifest."""

import argparse
import json
import lzma
from hashlib import sha1
from pathlib import Path

from island.charge_references.records import pack
from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_author import (
    CANONICAL_REPOSITORY,
    COMMIT,
    DATA_ID,
    PINS,
    REPOSITORY,
    load_author_opls_audit,
)
from island.forcefields.oplsaa.pdms_literature import (
    OPLSPDMSLiteratureError,
    build_opls_pdms_literature,
)
from island.forcefields.pdms_literature import INPUT_ERRORS, pinned_bytes
from island.workflows.storage import child
from scripts.audit_opls_pdms_author_sources import audit as independent_audit


def validate(directory, expected_sha256, *, repo, sources, references, engine_docs):
    root, repo = Path(directory), Path(repo)
    manifest = json.loads(pinned_bytes(root / "manifest.json", expected_sha256))
    require(
        manifest["identity"]
        == digest({k: v for k, v in manifest.items() if k != "identity"})
        and manifest["schema"] == "island_opls_pdms_author_manifest_v3"
        and manifest["data_identity"] == DATA_ID
        and manifest["source_hashes"] == PINS
        and manifest["commit"] == COMMIT
        and manifest["acceptance_gate"] == 1
        and manifest["execution_profiles"] == [],
        "changed_manifest_or_authority",
    )
    for name, pin in manifest["artifacts"].items():
        pinned_bytes(child(repo, name), pin)
    baseline = json.loads(
        lzma.decompress((root / "protected-file-hashes.json.xz").read_bytes())
    )
    require(
        len(baseline) == manifest["historical_files_protected"],
        "changed_history_inventory",
    )
    for name, pin in baseline.items():
        pinned_bytes(child(repo, name), pin)
    receipt = json.loads((root / "receipt.json").read_text())
    require(
        receipt["identity"]
        == digest({k: v for k, v in receipt.items() if k != "identity"})
        and receipt["schema"] == "island_opls_pdms_author_receipt_v3"
        and receipt["data_identity"] == DATA_ID
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["limitations"] == LIMITS
        and receipt["execution_profiles"] == receipt["backend_executions"] == []
        and receipt["profile_identity"] is None
        and receipt["acceptance_gate"] == receipt["native_pcff_acceptance_gate"] == 1,
        "changed_receipt_or_authority",
    )
    for name, field in (
        ("record.json", "record_sha256"),
        ("audit.json", "audit_sha256"),
        ("report.json", "report_sha256"),
        ("validation.json", "validation_sha256"),
        ("source-provenance.json", "source_provenance_sha256"),
    ):
        pinned_bytes(root / name, receipt[field])
    supplement = build_opls_pdms_literature(
        DomainSelection(
            (repo / "docs/evidence/phase_4j45/selection.json").read_text(),
            (repo / "docs/evidence/phase_4j44/source-audit.json").read_text(),
        )
    )
    result = load_author_opls_audit(
        root / "record.json",
        receipt["record_sha256"],
        sources=sources,
        supplement=supplement,
        references=references,
        engine_docs=engine_docs,
    )
    p = result.payload
    require(
        receipt["record_identity"] == p["identity"]
        and receipt["status"] == p["status"] == "scientifically_unqualified"
        and receipt["stages"] == p["stages"]
        and receipt["missing_dependencies"] == p["missing_dependencies"]
        and pack(json.loads((root / "audit.json").read_text())) == pack(p),
        "changed_audit_or_receipt",
    )
    report = json.loads((root / "report.json").read_text())
    require(
        report.pop("record_sha256") == receipt["record_sha256"]
        and all(p[k] == v for k, v in report.items()),
        "changed_report",
    )
    provenance = json.loads((root / "source-provenance.json").read_text())
    require(
        provenance["repository"] == REPOSITORY
        and provenance["canonical_repository"] == CANONICAL_REPOSITORY
        and provenance["repository_id"] == 902160948
        and provenance["commit"] == COMMIT
        and provenance["tree"] == "5b3bcc406226a9f5af9fcaee4abf5a465e9abc26"
        and provenance["repository_license"] is None
        and {f["path"]: f["sha256"] for f in provenance["files"]} == PINS,
        "changed_author_provenance",
    )
    for record in provenance["files"]:
        raw = pinned_bytes(child(sources, record["path"]), record["sha256"])
        require(
            len(raw) == record["size"]
            and sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
            == record["git_blob_sha1"],
            "changed_git_blob",
        )
    independent = independent_audit(sources, root / "audit.json", receipt["audit_sha256"])
    require(
        independent
        == json.loads((root / "independent-source-reconstruction.json").read_text()),
        "changed_independent_reconstruction",
    )
    try:
        result.require_executable()
    except OPLSPDMSLiteratureError as error:
        require(
            error.status == "scientifically_unqualified" and error.result is result,
            "lost_classified_blocked_result",
        )
    else:
        raise ValueError("blocked_audit_reached_execution")
    return {
        "manifest_valid": True,
        "receipt_valid": True,
        "record_identity": p["identity"],
        "source_files_reconstructed": len(PINS),
        "selected_coefficient_rows": len(p["selected_coefficient_rows"]),
        "historical_files_unchanged": len(baseline),
        "execution_profiles": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--sha256", required=True, help="Independent manifest SHA-256")
    for name in ("repo", "sources", "references", "engine-docs"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    try:
        print(
            json.dumps(
                validate(
                    args.directory,
                    args.sha256,
                    repo=args.repo,
                    sources=args.sources,
                    references=args.references,
                    engine_docs=args.engine_docs,
                )
            )
        )
    except INPUT_ERRORS as error:
        print(json.dumps({"status": "source_incompatible", "reason": str(error)}))
        raise SystemExit(1) from None
