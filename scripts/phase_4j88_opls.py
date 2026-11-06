"""Fenced J88 source-bound qualification with metadata-only external-byte replay."""

import argparse
import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from island.forcefields.assignment_core.schema import data, is_hash
from island.forcefields.bundle_adapters.opls_authorized import (
    capture_opls,
    capture_request,
    external_path,
    qualify_metadata,
)
from island.forcefields.bundle_adapters.schema import keys, loads
from scripts.phase_4j73_source_intake import exact, pin, read, require, signed, write
from scripts.phase_4j76_source_resolution import scope
from scripts.phase_4j84_bundles import spec_from_data
from scripts.phase_4j87_opls import REQUEST_PATH, REQUEST_SHA
from scripts.phase_4j87_opls import SCOPE as J87_SCOPE

BASE = "7d6fd22f49c0a28a2bfece088fb8c3673241fa92"
BRANCH = "codex/phase-4j88-authorized-opls-bundle-qualification"
ROOT = Path("docs/evidence/phase_4j88")
WHEN = "2026-11-06T00:00:00Z"
IMPLEMENTATION = (
    "src/island/forcefields/bundle_adapters/opls_authorized.py",
    "scripts/phase_4j88_opls.py",
    "scripts/validate_phase_4j88_receipt.py",
    "tests/test_phase_4j88.py",
    "docs/phase_4j88.md",
)
REPORTS = (
    "qualification-request.json",
    "source-qualification.json",
    "source-evidence-graph.json",
    "source-locations.json",
    "source-manifest.json",
    "graph-coverage-report.json",
    "convention-comparison.json",
    "missing-terms-and-blockers.json",
    "result.json",
)
SCOPE = {k.replace("J87_", "J88_"): v for k, v in J87_SCOPE.items()}
READY = "family_scoped_qualification_ready_runtime_blocked"


def reconstruct(
    repo,
    *,
    request=None,
    request_sha256=None,
    manifest=None,
    manifest_sha256=None,
    audit=None,
    audit_sha256=None,
    capture_timestamp=WHEN,
):
    """Every invocation reopens the pinned external files, including opaque originals."""
    scope()
    require(
        manifest is None or request is not None,
        "explicit_request_required_for_source_bundle",
    )
    require(request is not None or request_sha256 is None, "request_path_required")
    # The built-in no-source case has a fixed replay identity, including when
    # invoked through ingest. Explicit external requests retain their timestamp.
    if request is None:
        capture_timestamp = WHEN
    request, manifest, audit = (external_path(p) for p in (request, manifest, audit))
    rc = capture_request(
        request or repo / REQUEST_PATH,
        request_sha256 if request is not None else REQUEST_SHA,
    )
    sc = capture_opls(
        manifest,
        manifest_sha256,
        audit,
        audit_sha256,
        capture_timestamp=capture_timestamp,
    )
    result = qualify_metadata(rc, sc, spec_from_data)
    q = result["base_qualification"]
    inputs = {
        "schema": "island.j88.external-input-pins.v1",
        "mode": "explicit_external" if request is not None else "no_source_default",
        "request_sha256": rc["request_sha256"],
        "manifest_sha256": sc["manifest_sha256"],
        "audit_sha256": sc["audit_sha256"],
        "capture_timestamp": capture_timestamp,
        "raw_source_records_persisted": False,
    }
    missing = [
        {
            "graph_identity": r["graph_identity"],
            "missing_terms": r["missing_terms"],
            "blocking_reasons": r["blocking_reasons"],
            "blocker_categories": blocker_categories(r["blocking_reasons"]),
        }
        for r in result["graph_reports"]
    ]
    counters = {
        "authorized_scientific_bundles_supplied": int(q["source_authorized"]),
        "exact_graph_fixtures": len(result["fixture_files"]),
        "qualified_graph_reports": sum(
            r["decision"] == READY for r in result["graph_reports"]
        ),
        "source_files": len(sc["files"]),
        "source_records_inspected": len(result["source_locations"]),
        "missing_term_instances": sum(len(r["missing_terms"]) for r in missing),
        "unique_blocking_reasons": len(result["blocking_reasons"]),
        "accepted_coefficients": 0,
    }
    outputs = {
        "qualification-request.json": {
            "request_sha256": rc["request_sha256"],
            "request": result["request"],
            "fixture_files": result["fixture_files"],
            "scope_of_qualification": result["scope_of_qualification"],
            "chemical_family_transfer": "not_established",
        },
        "source-qualification.json": {
            "profile": q,
            "decision": result["decision"],
            "blocking_reasons": result["blocking_reasons"],
            "scope": SCOPE,
        },
        "source-evidence-graph.json": result["source_evidence"],
        "source-locations.json": {
            "records": result["source_locations"],
            "source_hashes": result["source_hashes"],
        },
        "source-manifest.json": {
            "manifest": sc["manifest"],
            "manifest_sha256": sc["manifest_sha256"],
            "audit_sha256": sc["audit_sha256"],
            "file_manifest": sc["manifest"]["files"] if sc["manifest"] else [],
            "external_files_required_during_replay": manifest is not None,
            "raw_source_records_persisted": False,
        },
        "graph-coverage-report.json": {
            "family": "OPLS-AA",
            "functional_class": "I",
            "reports": result["graph_reports"],
            "scope_of_qualification": "declared_exact_graphs_only",
            "chemical_family_transfer": "not_established",
            "scope": SCOPE,
        },
        "convention-comparison.json": result["convention_comparison"],
        "missing-terms-and-blockers.json": {
            "graphs": missing,
            "blocking_reasons": result["blocking_reasons"],
            "no_zero_substitution": True,
            "no_fallback": True,
        },
        "result.json": {
            "decision": result["decision"],
            "family": "OPLS-AA",
            "profile_id": result["request"]["profile_id"],
            "graph_identities": result["request"]["graph_identities"],
            **counters,
            "source_authorized": q["source_authorized"],
            "bundle_identity": q["bundle_identity"],
            "source_hashes": result["source_hashes"],
            "blocking_reasons": result["blocking_reasons"],
            "scope_of_qualification": "declared_exact_graphs_only",
            "chemical_family_transfer": "not_established",
            "scope": SCOPE,
        },
    }
    return inputs, {n: data(signed(v)) for n, v in outputs.items()}


def blocker_categories(reasons):
    """Overlapping diagnostic views; retain the full original list alongside."""
    groups = {
        "authorization_blocked": (
            "authoriz",
            "audit",
            "license",
            "provenance",
            "source_bundle",
        ),
        "convention_blocked": (
            "convention",
            "equation",
            "unit",
            "phase",
            "policy",
            "version",
            "form_",
        ),
        "source_location_blocked": (
            "source_location",
            "source_hash",
            "source_record",
            "claim",
        ),
        "ambiguous": ("ambig", "conflict"),
        "graph_scope_blocked": ("graph", "topology", "scope"),
        "missing": ("missing", "unresolved", "required", "no_evidenced"),
    }
    return {
        name: [
            reason for reason in reasons if any(v in reason.lower() for v in patterns)
        ]
        for name, patterns in groups.items()
    }


def intake_receipt(root, outputs):
    return signed(
        {
            "schema": "island.j88-qualification-receipt.v1",
            "scope": SCOPE,
            "snapshot_sha256": pin(root / "captured-input.json"),
            "reports": {n: pin(root / n) for n in REPORTS},
            "decision": outputs["result.json"]["decision"],
        }
    )


def emit(repo, root, **external):
    scope()
    require(not root.exists(), "do_not_overwrite_qualification_attempt")
    inputs, outputs = reconstruct(repo, **external)
    root.mkdir(parents=True)
    write(root / "captured-input.json", inputs)
    for name, value in outputs.items():
        write(root / name, value)
    write(root / "qualification-receipt.json", intake_receipt(root, outputs))
    return {
        "result": outputs["result.json"],
        "snapshot_sha256": pin(root / "captured-input.json"),
    }


def replay(root, snapshot_sha=None, *, repo=None, **external):
    scope()
    repo = repo or Path.cwd()
    inputs = loads((root / "captured-input.json").read_bytes())
    keys(
        inputs,
        (
            "schema",
            "mode",
            "request_sha256",
            "manifest_sha256",
            "audit_sha256",
            "capture_timestamp",
            "raw_source_records_persisted",
        ),
        "J88_snapshot_schema",
    )
    require(
        inputs["schema"] == "island.j88.external-input-pins.v1"
        and inputs["raw_source_records_persisted"] is False,
        "J88_snapshot_schema",
    )
    if inputs["mode"] == "no_source_default":
        exact(
            inputs,
            {
                "schema": "island.j88.external-input-pins.v1",
                "mode": "no_source_default",
                "request_sha256": REQUEST_SHA,
                "manifest_sha256": None,
                "audit_sha256": None,
                "capture_timestamp": WHEN,
                "raw_source_records_persisted": False,
            },
            "source_claims_in_no_source_snapshot",
        )
        require(not external, "no_source_default_cannot_substitute_inputs")
    else:
        require(inputs["mode"] == "explicit_external", "J88_replay_mode_invalid")
        require(is_hash(snapshot_sha), "external_snapshot_pin_required")
        for name in ("request", "manifest", "audit"):
            expected = inputs[name + "_sha256"]
            exact(
                external.get(name + "_sha256"),
                expected,
                "external_replay_pin_mutation:" + name,
            )
            require(
                (external.get(name) is None) == (expected is None),
                "explicit_external_replay_path_required:" + name,
            )
    if snapshot_sha is not None:
        require(
            pin(root / "captured-input.json") == snapshot_sha, "snapshot_hash_mutation"
        )
    fresh, outputs = reconstruct(
        repo, **external, capture_timestamp=inputs["capture_timestamp"]
    )
    exact(inputs, fresh, "replayed_input_identity_mutation")
    for name, value in outputs.items():
        exact(read(root / name), value, "resigned_semantic_mutation:" + name)
    exact(
        read(root / "qualification-receipt.json"),
        intake_receipt(root, outputs),
        "qualification_receipt_mutation",
    )
    return outputs["result.json"]


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
    for name in ("request", "manifest", "audit"):
        p.add_argument("--" + name, type=Path)
        p.add_argument("--" + name + "-sha256")
    p.add_argument("--snapshot-sha256")
    p.add_argument("--capture-timestamp")
    args = p.parse_args()
    external = {
        name: getattr(args, name)
        for name in (
            "request",
            "request_sha256",
            "manifest",
            "manifest_sha256",
            "audit",
            "audit_sha256",
        )
        if getattr(args, name) is not None
    }
    try:
        scope()
        if args.command == "replay":
            value = replay(
                args.directory, args.snapshot_sha256, repo=args.repo, **external
            )
        else:
            require(args.command == "ingest" or not external, "use_ingest_for_sources")
            value = emit(
                args.repo,
                args.directory,
                **external,
                capture_timestamp=args.capture_timestamp
                or (
                    datetime.now(UTC).isoformat().replace("+00:00", "Z")
                    if args.command == "ingest"
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
        # Malformed nested source/request objects never escape as Python tracebacks.
        reason = (
            str(error)
            if isinstance(error, ValueError)
            else "malformed_or_unavailable_input:" + type(error).__name__
        )
        print(
            json.dumps(
                {
                    "decision": "source_bundle_required"
                    if not args.audit
                    else "parameter_bundle_incomplete",
                    "status": "blocked",
                    "blocking_reasons": list(getattr(error, "blockers", (reason,))),
                    "scope": SCOPE,
                },
                sort_keys=True,
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
