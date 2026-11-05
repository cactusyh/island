"""Fenced, explicit local PDMS source intake; no scientific execution path."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from island.pcff_pdms.source_ingestion.pipeline import SCOPE, capture, normalize
from island.pcff_pdms.source_ingestion.schema import VERSION, schema
from scripts.phase_4j73_source_intake import (
    digest,
    exact,
    pin,
    read,
    require,
    signed,
    write,
)
from scripts.phase_4j76_source_resolution import scope

BASE = "637957fc8167e2f14f983cc364512da3fede9c52"
BRANCH = "codex/phase-4j81-user-supplied-pdms-source-ingestion"
ROOT = Path("docs/evidence/phase_4j81")
# Pinned after the no-source capture; supplied bundles require their own external pin.
BLOCKED_INPUT_SHA = "86a079c14db6d3461534010a24489be5f18a2c4703ace1cf1d56392755270b5c"
J80_ROWS_SHA = "c2badc816b348050e29c2bf1cbcc3c2833facb3ea58cd15510adee1dfdb7a496"
J80_RECEIPT_SHA = "0de6c9c78a77112706e9f64f2c9b6b7f56daf8f2f08885421ea04ef55dbdaa0b"
IMPLEMENTATION = (
    "src/island/pcff_pdms/source_ingestion/__init__.py",
    "src/island/pcff_pdms/source_ingestion/schema.py",
    "src/island/pcff_pdms/source_ingestion/parsing.py",
    "src/island/pcff_pdms/source_ingestion/pipeline.py",
    "scripts/phase_4j81_ingest.py",
    "scripts/validate_phase_4j81_receipt.py",
    "tests/test_phase_4j81.py",
    "docs/phase_4j81.md",
)
REPORTS = (
    "source-bundle-manifest.json",
    "source-provenance.json",
    "candidate-rows.json",
    "convention-comparison.json",
    "unresolved-terms.json",
    "result.json",
    "intake-receipt.json",
    "source-bundle.schema.json",
)


def empty_bundle():
    return {
        "schema": VERSION,
        "artifacts": [],
        "candidates": [],
        "convention_contract": {"fields": {}},
    }


def prepare(repo, manifest=None, manifest_sha=None, audit=None, audit_sha=None):
    scope()  # Before even stat/open of any supplied source or bundle.
    require(
        pin(repo / "docs/evidence/phase_4j80/receipt.json") == J80_RECEIPT_SHA,
        "J80_receipt_mutation",
    )
    require(
        pin(repo / "docs/evidence/phase_4j80/candidate-row-classification.json")
        == J80_ROWS_SHA,
        "J80_row_matrix_mutation",
    )
    targets = read(repo / "docs/evidence/phase_4j80/candidate-row-classification.json")[
        "extension_targets"
    ]
    parent = read(repo / "docs/evidence/phase_4j80/source-provenance-manifest.json")
    require(
        (manifest is None) == (manifest_sha is None),
        "manifest_and_external_hash_required",
    )
    require(
        (audit is None) == (audit_sha is None), "audit_and_independent_hash_required"
    )
    require(audit is None or manifest is not None, "audit_without_source_bundle")
    if manifest is not None:
        require(pin(manifest) == manifest_sha, "bundle_manifest_checksum_mismatch")
    if audit is not None:
        require(pin(audit) == audit_sha, "authority_audit_checksum_mismatch")
    value = read(manifest) if manifest else empty_bundle()
    return {
        "schema": "island.j81-captured-input.v1",
        "J80_receipt_sha256": J80_RECEIPT_SHA,
        "expected_J79_candidate_ids": [
            r["candidate_id"] for r in parent["row_evidence"]
        ],
        "J80_extension_targets": targets,
        "bundle_manifest_sha256": manifest_sha,
        "independent_audit_sha256": audit_sha,
        "capture": capture(value, manifest.parent if manifest else repo),
        "trusted_audits": read(audit) if audit else {},
    }


def reconstruct(inputs):
    require(
        inputs["J80_receipt_sha256"] == J80_RECEIPT_SHA, "historical_parent_mutation"
    )
    result = normalize(
        inputs["capture"],
        inputs["trusted_audits"],
        inputs["expected_J79_candidate_ids"],
    )
    return {
        "source-bundle-manifest.json": inputs["capture"]["manifest"],
        "source-provenance.json": signed(
            {
                "artifacts": inputs["capture"]["artifacts"],
                "independent_audit_sha256": inputs["independent_audit_sha256"],
                "authority_note": "Hash proves bytes; independently reviewed original-artifact audits, not user attestations, are required for verification.",
            }
        ),
        "candidate-rows.json": signed(
            {
                "allowed_statuses": list(result["candidate_status_counts"]),
                "rows": result["candidate_rows"],
                "counts": result["candidate_status_counts"],
            }
        ),
        "convention-comparison.json": signed(
            {
                "source_gate": result["source_gate"],
                "convention_gate": result["convention_comparison"],
                "reused_contract": "J79 comparison via unchanged J80 twelve-group convention gate",
                "parameter_conversion_allowed": False,
            }
        ),
        "unresolved-terms.json": signed(
            {
                **{
                    k: result[k]
                    for k in (
                        "missing_cross_term_classes",
                        "historical_candidates_still_blocked",
                        "unresolved_cross_terms",
                        "unresolved_end_group_terms",
                        "J80_preserved",
                    )
                },
                "extension_targets": inputs["J80_extension_targets"],
                "cross_term_ids": [
                    r["J78_parameter_id"]
                    for r in inputs["J80_extension_targets"]
                    if r["cross_term"]
                ],
                "end_group_term_ids": [
                    r["J78_parameter_id"]
                    for r in inputs["J80_extension_targets"]
                    if r["end_group"]
                ],
                "automatic_target_resolution": False,
            }
        ),
        "intake-receipt.json": signed(
            {
                "schema": "island.j81-intake-status.v1",
                "base": BASE,
                "branch": BRANCH,
                "decision": result["decision"],
                "captured_input_identity": digest(inputs),
                "bundle_manifest_sha256": inputs["bundle_manifest_sha256"],
                "independent_audit_sha256": inputs["independent_audit_sha256"],
                "verified_literature_rows": result["verified_literature_rows"],
                "scope": SCOPE,
                "qm_authorization_granted": False,
            }
        ),
        "result.json": signed(
            {k: v for k, v in result.items() if k != "candidate_rows"}
        ),
        "source-bundle.schema.json": schema(),
    }


def emit(directory, inputs):
    require(not directory.exists(), "do_not_overwrite_ingestion_attempt")
    directory.mkdir(parents=True)
    write(directory / "captured-input.json", inputs)
    for name, value in reconstruct(inputs).items():
        write(directory / name, value)
    return {
        "snapshot_sha256": pin(directory / "captured-input.json"),
        "decision": read(directory / "result.json")["decision"],
    }


def replay(directory, snapshot_sha=None):
    scope()
    expected = snapshot_sha or BLOCKED_INPUT_SHA
    require(
        pin(directory / "captured-input.json") == expected,
        "captured_input_missing_or_mutated_external_pin",
    )
    inputs = read(directory / "captured-input.json")
    if snapshot_sha is None:
        require(
            inputs["capture"]["manifest"] == empty_bundle()
            and not inputs["trusted_audits"],
            "source_claims_in_blocked_snapshot",
        )
    reports = reconstruct(inputs)
    for name, value in reports.items():
        exact(read(directory / name), value, "resigned_semantic_mutation:" + name)
    return reports["result.json"]


def historical(repo):
    rows = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    hashes = {}
    for row in rows.splitlines():
        meta, name = row.split("\t", 1)
        p = repo / name
        raw = (
            str(p.readlink()).encode() if meta.startswith("120000") else p.read_bytes()
        )
        git_blob = hashlib.sha1(
            b"blob " + str(len(raw)).encode() + b"\0" + raw
        ).hexdigest()
        require(git_blob == meta.split()[2], "historical_byte_mutation:" + name)
        sha = hashlib.sha256(raw).hexdigest()
        hashes[name] = sha
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("ingest", "replay", "schema"))
    parser.add_argument("--directory", type=Path, default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--audit-sha256")
    parser.add_argument("--snapshot-sha256")
    a = parser.parse_args()
    try:
        scope()
        if a.mode == "schema":
            value = schema()
        elif a.mode == "replay":
            value = replay(a.directory, a.snapshot_sha256)
        else:
            try:
                inputs = prepare(
                    a.repo, a.manifest, a.manifest_sha256, a.audit, a.audit_sha256
                )
                value = emit(a.directory, inputs)
            except (ValueError, OSError, TypeError, KeyError) as e:
                require(not a.directory.exists(), "do_not_overwrite_failed_attempt")
                write(
                    a.directory / "blocked-receipt.json",
                    signed(
                        {
                            "schema": "island.j81-invalid-input.v1",
                            "decision": "source_ingestion_blocked",
                            "reason": str(e),
                            "scope": SCOPE,
                        }
                    ),
                )
                raise
        print(json.dumps(value, sort_keys=True))
        return 0
    except (ValueError, OSError, TypeError, KeyError, RuntimeError) as e:
        print(
            json.dumps(
                {
                    "decision": "source_ingestion_blocked",
                    "reason": str(e),
                    "scope": SCOPE,
                },
                sort_keys=True,
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
