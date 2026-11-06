"""Fenced, single-profile OPLS-AA qualification. No source discovery or activation."""

import argparse
import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from island.forcefields.assignment_core.schema import READINESS, data
from island.forcefields.bundle_adapters.opls_qualification import (
    capture_opls,
    capture_request,
    qualify_opls,
)
from island.forcefields.bundle_adapters.schema import keys, loads
from scripts.phase_4j73_source_intake import exact, pin, read, require, signed, write
from scripts.phase_4j76_source_resolution import scope
from scripts.phase_4j84_bundles import spec_from_data

BASE = "fe4e82de3ff1b0704d600e226e79c66cdccec2be"
BRANCH = "codex/phase-4j87-first-authorized-opls-qualification"
ROOT = Path("docs/evidence/phase_4j87")
REQUEST_PATH = "tests/fixtures/phase_4j87/request.json"
REQUEST_SHA = "68f53de2c64a4d66af2cde976a3e48b7c09959137730eabbced307a63a61949c"
WHEN = "2026-11-06T00:00:00Z"
IMPLEMENTATION = (
    "src/island/forcefields/bundle_adapters/opls_qualification.py",
    "scripts/phase_4j87_opls.py",
    "scripts/validate_phase_4j87_receipt.py",
    "tests/test_phase_4j87.py",
    "docs/phase_4j87.md",
    REQUEST_PATH,
    *(
        "tests/fixtures/phase_4j87/graphs/" + n + ".json"
        for n in (
            "hydrocarbon",
            "aromatic",
            "ether",
            "ester",
            "carbonyl",
            "amide",
            "nitrile",
            "halogenated",
        )
    ),
)
REPORTS = (
    "qualification-request.json",
    "source-qualification.json",
    "source-evidence-graph.json",
    "source-locations.json",
    "graph-coverage-report.json",
    "convention-comparison.json",
    "result.json",
)
SCOPE = {
    **READINESS,
    "J87_QM_attempts": 0,
    "J87_MD_steps": 0,
    "converter_attempts": 0,
    "physical_fitting_executions": 0,
    "parameter_generation_executions": 0,
    "scientific_energy_force_evaluations": 0,
    "reserved_validation": 0,
    "production_executions": 0,
    "runtime_registration": False,
    "scientific_validation": "not_established",
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
}


def reconstruct(inputs):
    keys(inputs, ("schema", "request_capture", "source_capture"), "J87_snapshot_schema")
    require(inputs["schema"] == "island.j87.capture.v1", "J87_snapshot_schema")
    result = qualify_opls(
        inputs["request_capture"], inputs["source_capture"], spec_from_data
    )
    q = result["base_qualification"]
    outputs = {
        "qualification-request.json": {
            "request_sha256": inputs["request_capture"]["request_sha256"],
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
        "graph-coverage-report.json": {
            "functional_class": "I",
            "family": "OPLS-AA",
            "reports": result["graph_reports"],
            "scope_of_qualification": "declared_exact_graphs_only",
            "chemical_family_transfer": "not_established",
            "scope": SCOPE,
        },
        "convention-comparison.json": result["convention_comparison"],
        "result.json": {
            "decision": result["decision"],
            "family": "OPLS-AA",
            "profile_id": result["request"]["profile_id"],
            "graph_identities": result["request"]["graph_identities"],
            "exact_graph_fixtures": len(result["fixture_files"]),
            "source_authorized": q["source_authorized"],
            "bundle_identity": q["bundle_identity"],
            "source_hashes": result["source_hashes"],
            "qualified_graph_reports": sum(
                r["decision"] == "family_scoped_qualification_ready_runtime_blocked"
                for r in result["graph_reports"]
            ),
            "blocking_reasons": result["blocking_reasons"],
            "scope_of_qualification": "declared_exact_graphs_only",
            "chemical_family_transfer": "not_established",
            "scope": SCOPE,
        },
    }
    return {n: data(signed(v)) for n, v in outputs.items()}


def intake_receipt(root, outputs):
    return signed(
        {
            "schema": "island.j87-qualification-receipt.v1",
            "scope": SCOPE,
            "snapshot_sha256": pin(root / "captured-input.json"),
            "reports": {n: pin(root / n) for n in REPORTS},
            "decision": outputs["result.json"]["decision"],
        }
    )


def emit(
    repo,
    root,
    *,
    request=None,
    request_sha256=None,
    manifest=None,
    manifest_sha256=None,
    audit=None,
    audit_sha256=None,
    capture_timestamp=WHEN,
):
    scope()
    require(not root.exists(), "do_not_overwrite_qualification_attempt")
    require(
        manifest is None or request is not None,
        "explicit_request_required_for_source_bundle",
    )
    require(request is not None or request_sha256 is None, "request_path_required")
    request_cap = capture_request(
        request or repo / REQUEST_PATH,
        request_sha256 if request is not None else REQUEST_SHA,
    )
    source_cap = capture_opls(
        manifest,
        manifest_sha256,
        audit,
        audit_sha256,
        capture_timestamp=capture_timestamp,
    )
    inputs = {
        "schema": "island.j87.capture.v1",
        "request_capture": request_cap,
        "source_capture": source_cap,
    }
    outputs = reconstruct(inputs)
    root.mkdir(parents=True)
    write(root / "captured-input.json", inputs)
    for name, value in outputs.items():
        write(root / name, value)
    write(root / "qualification-receipt.json", intake_receipt(root, outputs))
    return {
        "result": outputs["result.json"],
        "snapshot_sha256": pin(root / "captured-input.json"),
    }


def replay(root, snapshot_sha=None):
    scope()
    inputs = loads((root / "captured-input.json").read_bytes())
    if snapshot_sha is not None:
        require(
            pin(root / "captured-input.json") == snapshot_sha, "snapshot_hash_mutation"
        )
    else:
        require(
            inputs["request_capture"]["request_sha256"] == REQUEST_SHA,
            "default_request_pin_mutation",
        )
        exact(
            inputs["source_capture"],
            capture_opls(capture_timestamp=WHEN),
            "source_claims_in_no_source_snapshot",
        )
    outputs = reconstruct(inputs)
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("emit", "ingest", "replay"))
    parser.add_argument("--directory", type=Path, default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    for name in ("request", "manifest", "audit"):
        parser.add_argument("--" + name, type=Path)
        parser.add_argument("--" + name + "-sha256")
    parser.add_argument("--snapshot-sha256")
    parser.add_argument("--capture-timestamp")
    args = parser.parse_args()
    try:
        scope()
        if args.command == "replay":
            value = replay(args.directory, args.snapshot_sha256)
        else:
            require(
                args.command == "ingest"
                or (args.manifest is None and args.audit is None),
                "use_ingest_for_sources",
            )
            value = emit(
                args.repo,
                args.directory,
                request=args.request,
                request_sha256=args.request_sha256,
                manifest=args.manifest,
                manifest_sha256=args.manifest_sha256,
                audit=args.audit,
                audit_sha256=args.audit_sha256,
                capture_timestamp=args.capture_timestamp
                or (
                    datetime.now(UTC).isoformat().replace("+00:00", "Z")
                    if args.command == "ingest"
                    else WHEN
                ),
            )
        print(json.dumps(value, sort_keys=True))
        return 0
    except (ValueError, OSError, TypeError, KeyError, RuntimeError) as error:
        print(
            json.dumps(
                {
                    "decision": "parameter_bundle_incomplete",
                    "reason": str(error),
                    "scope": SCOPE,
                }
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
