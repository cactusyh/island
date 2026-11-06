"""J89: all-eight-graph, externally pinned OPLS qualification; never activation."""

import argparse
import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from island.forcefields.assignment_core.schema import data, is_hash
from island.forcefields.bundle_adapters.intake import _read
from island.forcefields.bundle_adapters.opls_authorized import (
    ALLOWED_FIXTURES,
    capture_request,
    external_path,
)
from island.forcefields.bundle_adapters.schema import keys, loads
from scripts import phase_4j88_opls as inherited
from scripts.phase_4j73_source_intake import exact, pin, read, require, signed, write
from scripts.phase_4j76_source_resolution import scope

BASE = "4ed7377a66c02362d69da0840b7701ba2b6a797f"
BRANCH = "codex/phase-4j89-first-real-authorized-opls-qualification"
ROOT = Path("docs/evidence/phase_4j89")
WHEN = "2026-11-06T00:00:00Z"
IMPLEMENTATION = (
    "scripts/phase_4j89_opls.py",
    "scripts/validate_phase_4j89_receipt.py",
    "tests/test_phase_4j89.py",
    "docs/phase_4j89.md",
)
SCOPE = {k.replace("J88_", "J89_"): v for k, v in inherited.SCOPE.items()}
READY = inherited.READY
INPUTS = ("request", "manifest", "audit")
REPORTS = (*inherited.REPORTS, "input-gate.json", "audit-references.json")


def supplied_inputs(external):
    """No IO: incomplete inputs cannot trigger inspection of any external file."""
    require(
        set(external) <= {k for n in INPUTS for k in (n, n + "_sha256")},
        "unknown_external_input",
    )
    return {
        n: {
            "path_supplied": external.get(n) is not None,
            "sha256": external.get(n + "_sha256"),
        }
        for n in INPUTS
    }


def input_blockers(provided):
    reasons = []
    for name in INPUTS:
        if not provided[name]["path_supplied"]:
            reasons.append("required_external_path_missing:" + name)
        if not is_hash(provided[name]["sha256"]):
            reasons.append("required_external_SHA256_pin_missing_or_invalid:" + name)
    return sorted(reasons)


def exact_eight_request(path, checksum):
    captured = capture_request(external_path(path), checksum)
    request = loads(captured["request_raw"])
    require(
        {f["id"] for f in request["fixtures"]} == set(ALLOWED_FIXTURES)
        and len(request["fixtures"]) == 8,
        "all_eight_immutable_J87_fixtures_required",
    )
    return request


def audit_identity(path, checksum):
    """Reviewer/provenance are caller-trusted through the independent audit pin."""
    audit = loads(_read(external_path(path), checksum))
    require(isinstance(audit, dict), "independent_audit_schema")
    invalid = {"unknown", "missing", "unresolved", "default", "n/a", "none"}
    reasons = []
    for name in ("reviewer", "authority_basis"):
        value = audit.get(name)
        if (
            not isinstance(value, str)
            or not value.strip()
            or value.strip().lower() in invalid
        ):
            reasons.append(
                "trusted_independent_audit_identity_or_provenance_required:" + name
            )
    if audit.get("policy") != "independent_source_review":
        reasons.append("independent_audit_policy_required")
    return audit, sorted(reasons)


def _audit_references(audit, checksum, authorized):
    if audit is None:
        return {
            "status": "not_supplied",
            "audit_sha256": None,
            "records": [],
            "trust_basis": "No independent audit supplied; no source authority established.",
        }
    rows = []
    for i, review in enumerate(audit["source_reviews"]):
        for j, claim in enumerate(review["claims"]):
            rows.append(
                {
                    "source_file": review["path"],
                    "source_sha256": review["sha256"],
                    "record_id": claim["record_id"],
                    "claim_sha256": claim["claim_sha256"],
                    "audit_pointer": f"/source_reviews/{i}/claims/{j}",
                    "audit_sha256": checksum,
                }
            )
    return {
        "status": "independently_audited_source"
        if authorized
        else "source_authorization_blocked",
        "audit_sha256": checksum,
        "reviewer": audit["reviewer"],
        "authority_basis": audit["authority_basis"],
        "policy": audit["policy"],
        "review_timestamp": audit["review_timestamp"],
        "manifest_identity": audit["manifest_identity"],
        "records": rows,
        "trust_basis": "Caller-supplied trusted independent audit SHA-256; hashes alone do not authenticate a reviewer.",
    }


def reconstruct(repo, *, capture_timestamp=WHEN, **external):
    scope()
    provided = supplied_inputs(external)
    blockers = input_blockers(provided)
    audit = None
    if not blockers:
        blockers = sorted(
            "required_external_file_unavailable:" + name
            for name in INPUTS
            if not external_path(external[name]).is_file()
        )
    if not blockers:
        # These are pinned JSON/fixture files. Parameter sources remain unopened
        # until the complete request scope and independent audit identity pass.
        exact_eight_request(external["request"], external["request_sha256"])
        audit, blockers = audit_identity(external["audit"], external["audit_sha256"])
    if blockers:
        # This constructs an explicit empty obligation report from repository
        # graph fixtures. It supplies no source facts or scientific coverage.
        _, reports = inherited.reconstruct(repo)
        capture_timestamp = WHEN
    else:
        _, reports = inherited.reconstruct(
            repo, capture_timestamp=capture_timestamp, **external
        )
    # Strip inherited signatures, replace only phase metadata, then sign J89 facts.
    reports = {
        n: {k: v for k, v in r.items() if k != "identity"} for n, r in reports.items()
    }
    for report in reports.values():
        if "scope" in report:
            report["scope"] = SCOPE
    result = reports["result.json"]
    if blockers:
        result["decision"] = "source_bundle_required"
        result["blocking_reasons"] = sorted(
            set(result["blocking_reasons"]) | set(blockers)
        )
        result["unique_blocking_reasons"] = len(result["blocking_reasons"])
        reports["source-qualification.json"].update(
            decision="source_bundle_required",
            blocking_reasons=result["blocking_reasons"],
        )
        for graph in reports["graph-coverage-report.json"]["reports"]:
            graph["blocking_reasons"] = sorted(
                set(graph["blocking_reasons"]) | set(blockers)
            )
        missing = reports["missing-terms-and-blockers.json"]
        missing["blocking_reasons"] = result["blocking_reasons"]
        for graph in missing["graphs"]:
            graph["blocking_reasons"] = sorted(
                set(graph["blocking_reasons"]) | set(blockers)
            )
            graph["blocker_categories"] = inherited.blocker_categories(
                graph["blocking_reasons"]
            )
    audit_refs = _audit_references(
        audit if not blockers else None,
        external.get("audit_sha256"),
        result["source_authorized"],
    )
    if audit is not None and blockers:
        audit_refs = {
            "status": "independent_audit_untrusted_or_incomplete",
            "audit_sha256": external.get("audit_sha256"),
            "records": [],
            "trust_basis": "Trusted independent reviewer identity/provenance not established.",
        }
    by_record = {(r["source_file"], r["record_id"]): r for r in audit_refs["records"]}
    if not blockers:

        def bind(row):
            reference = by_record.get((row["source_file"], row["record_id"]))
            require(
                reference is not None
                and reference["claim_sha256"] == row["claim_sha256"]
                and reference["source_sha256"] == row["source_sha256"],
                "independent_audit_source_claim_reference_mismatch",
            )
            row["audit_reference"] = reference

        for row in reports["source-locations.json"]["records"]:
            bind(row)
        for row in reports["source-evidence-graph.json"]["records"]:
            bind(row)
        for graph in reports["graph-coverage-report.json"]["reports"]:
            for row in graph["matched_source_records"]:
                bind(row)
    result.update(
        {
            "required_graph_count": 8,
            "parameter_rows_inspected": sum(
                r["kind"] == "parameter"
                for r in reports["source-locations.json"]["records"]
            ),
            "request_sha256": external.get("request_sha256") if not blockers else None,
            "manifest_sha256": external.get("manifest_sha256")
            if not blockers
            else None,
            "audit_sha256": external.get("audit_sha256") if not blockers else None,
            "source_bytes_rehashed": not blockers,
            "chemical_family_transfer": "not_established",
        }
    )
    reports["input-gate.json"] = {
        "status": "required_external_inputs_complete"
        if not blockers
        else "source_bundle_required",
        "provided": provided,
        "blocking_reasons": blockers,
        "external_source_inspection_performed": not blockers,
        "required_graph_count": 8,
    }
    reports["audit-references.json"] = audit_refs
    inputs = {
        "schema": "island.j89.external-input-pins.v1",
        "provided": provided,
        "capture_timestamp": capture_timestamp,
        "external_sources_used": not blockers,
        "raw_source_records_persisted": False,
    }
    return inputs, {n: data(signed(r)) for n, r in reports.items()}


def intake_receipt(root, outputs):
    return signed(
        {
            "schema": "island.j89-qualification-receipt.v1",
            "scope": SCOPE,
            "snapshot_sha256": pin(root / "captured-input.json"),
            "reports": {n: pin(root / n) for n in REPORTS},
            "decision": outputs["result.json"]["decision"],
        }
    )


def emit(repo, root, **external):
    scope()
    require(not root.exists(), "do_not_overwrite_qualification_attempt")
    inputs, reports = reconstruct(repo, **external)
    root.mkdir(parents=True)
    write(root / "captured-input.json", inputs)
    for name, report in reports.items():
        write(root / name, report)
    write(root / "qualification-receipt.json", intake_receipt(root, reports))
    return {
        "result": reports["result.json"],
        "snapshot_sha256": pin(root / "captured-input.json"),
    }


def replay(root, snapshot_sha=None, *, repo=None, **external):
    scope()
    inputs = loads((root / "captured-input.json").read_bytes())
    keys(
        inputs,
        (
            "schema",
            "provided",
            "capture_timestamp",
            "external_sources_used",
            "raw_source_records_persisted",
        ),
        "J89_snapshot_schema",
    )
    require(
        inputs["schema"] == "island.j89.external-input-pins.v1"
        and inputs["raw_source_records_persisted"] is False,
        "J89_snapshot_schema",
    )
    keys(inputs["provided"], INPUTS, "J89_external_input_schema")
    for item in inputs["provided"].values():
        keys(item, ("path_supplied", "sha256"), "J89_external_input_schema")
        require(type(item["path_supplied"]) is bool, "J89_external_path_flag")
    if inputs["provided"] != supplied_inputs({}):
        require(is_hash(snapshot_sha), "external_snapshot_pin_required")
    if snapshot_sha is not None:
        require(
            pin(root / "captured-input.json") == snapshot_sha, "snapshot_hash_mutation"
        )
    exact(
        inputs["provided"],
        supplied_inputs(external),
        "explicit_external_replay_paths_and_pins_required",
    )
    fresh, reports = reconstruct(
        repo or Path.cwd(), capture_timestamp=inputs["capture_timestamp"], **external
    )
    require(
        not inputs["external_sources_used"] or fresh["external_sources_used"],
        "external_originals_unavailable_or_independent_authority_blocked",
    )
    exact(inputs, fresh, "replayed_input_identity_mutation")
    for name, report in reports.items():
        exact(read(root / name), report, "resigned_semantic_mutation:" + name)
    exact(
        read(root / "qualification-receipt.json"),
        intake_receipt(root, reports),
        "qualification_receipt_mutation",
    )
    return reports["result.json"]


def historical(repo):
    entries = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    hashes = {}
    for row in entries.splitlines():
        meta, name = row.split("\t", 1)
        path = repo / name
        raw = (
            str(path.readlink()).encode()
            if meta.startswith("120000")
            else path.read_bytes()
        )
        blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        require(blob == meta.split()[2], "historical_byte_mutation:" + name)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    return hashes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("emit", "ingest", "replay"))
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    for name in INPUTS:
        p.add_argument("--" + name, type=Path)
        p.add_argument("--" + name + "-sha256")
    p.add_argument("--snapshot-sha256")
    p.add_argument("--capture-timestamp")
    a = p.parse_args()
    external = {
        k: getattr(a, k)
        for n in INPUTS
        for k in (n, n + "_sha256")
        if getattr(a, k) is not None
    }
    try:
        scope()
        if a.command == "replay":
            value = replay(a.directory, a.snapshot_sha256, repo=a.repo, **external)
        else:
            require(a.command == "ingest" or not external, "use_ingest_for_sources")
            value = emit(
                a.repo,
                a.directory,
                **external,
                capture_timestamp=a.capture_timestamp
                or (
                    datetime.now(UTC).isoformat().replace("+00:00", "Z")
                    if external
                    else WHEN
                ),
            )
        print(json.dumps(value, sort_keys=True))
        return 0
    except (
        ValueError,
        OSError,
        TypeError,
        KeyError,
        RuntimeError,
        AttributeError,
        IndexError,
        OverflowError,
    ) as error:
        reason = (
            str(error)
            if isinstance(error, ValueError)
            else "malformed_or_unavailable_input:" + type(error).__name__
        )
        print(
            json.dumps(
                {
                    "status": "blocked",
                    "decision": "source_bundle_required"
                    if reason.startswith(
                        ("source_file_unavailable", "audit_", "independent_audit_")
                    )
                    else "parameter_bundle_incomplete",
                    "blocking_reasons": list(getattr(error, "blockers", (reason,))),
                    "scope": SCOPE,
                },
                sort_keys=True,
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
