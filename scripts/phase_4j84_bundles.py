"""Fenced, reproducible software coverage audit. Not a scientific parameter assignment."""

import argparse
import hashlib
import json
import subprocess
from dataclasses import fields as dataclass_fields
from datetime import UTC
from pathlib import Path

from island.forcefields.assignment_core import builtin_contracts, normalize, type_graph
from island.forcefields.assignment_core.inventory import requests
from island.forcefields.assignment_core.profiles import (
    required_conventions,
    required_families,
)
from island.forcefields.assignment_core.schema import (
    CROSS_TERMS,
    FAMILIES,
    READINESS,
    Branch,
    Crosslink,
    PolymerSpec,
    data,
    unique,
)
from island.forcefields.bundle_adapters import ADAPTERS, adapt, capture
from island.forcefields.bundle_adapters.schema import schemas
from island.forcefields.coverage_registry import ForceFieldProfileRegistry
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

BASE = "8a1c92c992951edbf34dd73c776df399443e2f06"
BRANCH = "codex/phase-4j84-audited-parameter-bundle-adapters"
ROOT = Path("docs/evidence/phase_4j84")
FIXTURE = "tests/fixtures/phase_4j84/organic-motifs.json"
FIXTURE_SHA = "0c3a17d70dde6a99ecddcc3fbd6e235af26a40ddfcbd44ae4b48134992f0ad8a"
IMPLEMENTATION = (
    "src/island/forcefields/bundle_adapters/__init__.py",
    "src/island/forcefields/bundle_adapters/schema.py",
    "src/island/forcefields/bundle_adapters/intake.py",
    "src/island/forcefields/bundle_adapters/adapter.py",
    "scripts/phase_4j84_bundles.py",
    "scripts/validate_phase_4j84_receipt.py",
    "tests/test_phase_4j84.py",
    FIXTURE,
    "docs/phase_4j84.md",
)
REPORTS = (
    "bundle-ingestion.json",
    "bundle-schemas.json",
    "profile-adapters.json",
    "organic-coverage-report.json",
    "completeness-report.json",
    "result.json",
)

SCOPE = {
    **READINESS,
    "J84_QM_attempts": 0,
    "J84_MD_steps": 0,
    "converter_attempts": 0,
    "physical_fitting_executions": 0,
    "parameter_generation_executions": 0,
    "scientific_energy_force_evaluations": 0,
    "reserved_validation": 0,
    "runtime_registration": False,
    "production_executions": 0,
    "scientific_validation": "not_established",
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
}


def spec_from_data(row):
    allowed = {f.name for f in dataclass_fields(PolymerSpec)}
    require(set(row) - {"id"} <= allowed, "unknown_polymer_spec_field")
    fields = {k: v for k, v in row.items() if k != "id"}
    for k in ("repeat_units", "connection_points", "periodic_topology"):
        if k in fields:
            fields[k] = tuple(tuple(v) for v in fields[k])
    for k in ("sequence", "end_groups"):
        if k in fields:
            fields[k] = tuple(fields[k])
    for k, cls in (("branches", Branch), ("crosslinks", Crosslink)):
        if k in fields:
            fields[k] = tuple(cls(**v) for v in fields[k])
    return PolymerSpec(**fields)


def reconstruct(inputs):
    require(
        inputs["schema"] == "island.j84.organic-benchmark.v1"
        and inputs["synthetic_controls_only"] is True,
        "benchmark_schema_or_scientific_promotion",
    )
    unique([r["id"] for r in inputs["motifs"]], "duplicate_motif_identity")
    outcomes = [adapt(c) for c in inputs["bundle_captures"]]
    active = [o.profile for o in outcomes if o.profile is not None]
    unique([p.family for p in active], "multiple_bundles_for_family_no_overlay")
    unique([o.bundle_identity for o in outcomes], "duplicate_bundle_identity")
    profiles = tuple(active) + tuple(
        p for p in builtin_contracts() if p.family not in {p.family for p in active}
    )
    trusted = tuple(e for o in outcomes for e in o.trusted_evidence)
    registry = ForceFieldProfileRegistry(profiles, trusted)
    reports, completeness = [], []
    for row in sorted(inputs["motifs"], key=lambda r: r["id"]):
        graph = normalize(spec_from_data(row))
        report = registry.coverage(graph)
        reports.append({"motif": row["id"], "report": data(report)})
        for p in registry.profiles:
            core = next(r for r in report.profile_reports if r.profile_id == p.id)
            reqs = requests(graph, type_graph(graph, p, trusted), p)
            matched = dict(core.matched_records)
            terms = {
                f: {
                    "required": sum(r.family == f for r in reqs),
                    "covered": sum(r.family == f and r.id in matched for r in reqs),
                }
                for f in required_families(p)
            }
            completeness.append(
                {
                    "motif": row["id"],
                    "profile_id": p.id,
                    "family": p.family,
                    "functional_class": p.functional_class,
                    "terms": terms,
                    "required_conventions": required_conventions(p),
                    "blocking_reasons": core.blocking_reasons,
                    "missing_terms": data(core.missing_terms),
                    "source_hashes": core.source_hashes,
                    "replay_identity": core.replay_identity,
                    "out_of_plane_contract": "improper form and ordering must be explicitly evidenced",
                }
            )
    supplied = {
        c["manifest"]["family"] for c in inputs["bundle_captures"] if c["manifest"]
    }
    all_hashes = sorted({h for o in outcomes for h in o.source_hashes})
    output = {
        "bundle-ingestion.json": signed(
            {
                "outcomes": data(outcomes),
                "evidence_graphs": [
                    json.loads(o.evidence_graph_json) for o in outcomes
                ],
                "source_hashes": all_hashes,
                "blocked_missing_bundles": sorted(set(FAMILIES) - supplied),
                "source_authority_note": "No source discovery; external manifest and independent audit pins are required.",
            }
        ),
        "bundle-schemas.json": schemas(),
        "profile-adapters.json": signed(
            {
                "adapters": ADAPTERS,
                "profiles": data(registry.profiles),
                "registry_identity": registry.identity,
                "no_conversion": True,
            }
        ),
        "organic-coverage-report.json": signed(
            {
                "synthetic_geometry_fixtures_only": True,
                "scientific_validation": "not_established",
                "motifs": reports,
                "source_evidence": [
                    json.loads(o.evidence_graph_json) for o in outcomes
                ],
            }
        ),
        "completeness-report.json": signed(
            {
                "profiles": data(completeness),
                "class_i_cross_terms": "not_required_by_class_i_contract",
                "class_ii_cross_terms": CROSS_TERMS,
            }
        ),
        "result.json": signed(
            {
                "decision": "audited_bundle_ingestion_ready_runtime_blocked",
                "scope": SCOPE,
                "motifs": len(reports),
                "profile_contracts": len(profiles),
                "candidate_decisions": len(reports) * len(profiles),
                "class_i_reports": sum(
                    c["functional_class"] == "I" for c in completeness
                ),
                "class_ii_reports": sum(
                    c["functional_class"] == "II" for c in completeness
                ),
                "authorized_bundles_supplied": len(outcomes),
                "audited_bundle_candidates": len(active),
                "source_hashes": all_hashes,
                "accepted_coefficients": 0,
                "scientific_parameter_rows_added": 0,
                "replay_identity": digest(reports),
                "production_qualification_established": False,
                "mapped_source_parameter_rows": sum(len(p.parameters) for p in active),
                "coverage_decisions": {
                    k: sum(
                        c["coverage_decision"] == k
                        for r in reports
                        for c in r["report"]["capability_matrix"]["cells"]
                    )
                    for k in ("accepted", "provisional", "blocked")
                },
            }
        ),
    }
    return {k: data(v) for k, v in output.items()}


def emit(repo, root, bundle_capture=None):
    scope()
    require(pin(repo / FIXTURE) == FIXTURE_SHA, "fixture_pin_mismatch")
    require(not root.exists(), "do_not_overwrite_bundle_attempt")
    inputs = read(repo / FIXTURE)
    if bundle_capture is not None:
        inputs["bundle_captures"] = [bundle_capture]
    reports = reconstruct(inputs)
    root.mkdir(parents=True)
    if bundle_capture is None:
        (root / "captured-input.json").write_bytes((repo / FIXTURE).read_bytes())
    else:
        write(root / "captured-input.json", inputs)
    for name, value in reports.items():
        write(root / name, value)
    return {
        "result": reports["result.json"],
        "snapshot_sha256": pin(root / "captured-input.json"),
    }


def replay(root, snapshot_sha=None):
    scope()
    require(
        pin(root / "captured-input.json") == (snapshot_sha or FIXTURE_SHA),
        "captured_fixture_pin_mismatch",
    )
    inputs = read(root / "captured-input.json")
    if snapshot_sha is None:
        require(not inputs["bundle_captures"], "source_claims_in_no_source_snapshot")
    reports = reconstruct(inputs)
    for name, value in reports.items():
        exact(read(root / name), value, "resigned_semantic_mutation:" + name)
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
    p.add_argument("mode", choices=("emit", "ingest", "replay", "schema"))
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--manifest", type=Path)
    p.add_argument("--manifest-sha256")
    p.add_argument("--audit", type=Path)
    p.add_argument("--audit-sha256")
    p.add_argument("--capture-timestamp")
    p.add_argument("--snapshot-sha256")
    a = p.parse_args()
    try:
        scope()
        if a.mode == "schema":
            value = schemas()
        elif a.mode == "replay":
            value = replay(a.directory, a.snapshot_sha256)
        elif a.mode == "ingest":
            from datetime import datetime

            when = a.capture_timestamp or datetime.now(UTC).isoformat().replace(
                "+00:00", "Z"
            )
            cap = capture(
                a.manifest,
                a.manifest_sha256,
                a.audit,
                a.audit_sha256,
                capture_timestamp=when,
            )
            value = emit(a.repo, a.directory, cap)
        else:
            require(
                a.manifest is None and a.audit is None,
                "use_ingest_for_explicit_sources",
            )
            value = emit(a.repo, a.directory)
        print(json.dumps(value, sort_keys=True))
        return 0
    except (ValueError, OSError, TypeError, KeyError, RuntimeError) as e:
        print(json.dumps({"status": "blocked", "reason": str(e), "scope": SCOPE}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
