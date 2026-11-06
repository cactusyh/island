"""Fenced source qualification and organic coverage. No runtime assignment."""

import argparse
import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from island.forcefields.assignment_core import normalize
from island.forcefields.assignment_core.schema import FAMILIES, READINESS, data, unique
from island.forcefields.bundle_adapters import capture
from island.forcefields.bundle_adapters.qualification import (
    qualify,
    rank_qualifications,
)
from island.forcefields.bundle_adapters.schema import schemas
from scripts.phase_4j73_source_intake import exact, pin, read, require, signed, write
from scripts.phase_4j76_source_resolution import scope
from scripts.phase_4j84_bundles import spec_from_data

BASE = "f5b9abcf459f541134d4f1162fb3af12f454a1f5"
BRANCH = "codex/phase-4j85-authorized-bundle-qualification"
ROOT = Path("docs/evidence/phase_4j85")
FIXTURE = "tests/fixtures/phase_4j85/organic-motifs.json"
FIXTURE_SHA = "f97e6132f86b655788504316bb050b0602638d05b32447000066f3f774fb278a"
IMPLEMENTATION = (
    "src/island/forcefields/bundle_adapters/qualification.py",
    "scripts/phase_4j85_qualification.py",
    "scripts/validate_phase_4j85_receipt.py",
    "tests/test_phase_4j85.py",
    FIXTURE,
    "docs/phase_4j85.md",
)
REPORTS = (
    "bundle-requirements.json",
    "source-qualification.json",
    "source-evidence-graph.json",
    "organic-coverage-report.json",
    "class-i-completeness.json",
    "class-ii-completeness.json",
    "result.json",
)
SCOPE = {
    **READINESS,
    "J85_QM_attempts": 0,
    "J85_MD_steps": 0,
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
    require(inputs["schema"] == "island.j85.organic-benchmark.v1", "benchmark_schema")
    require(inputs["geometry_fixtures_only"] is True, "geometry_fixture_promotion")
    unique([r["id"] for r in inputs["motifs"]], "duplicate_motif_identity")
    graphs = {r["id"]: normalize(spec_from_data(r)) for r in inputs["motifs"]}
    require(graphs, "benchmark_requires_motifs")
    captures = inputs["bundle_captures"]
    unique([c["manifest"]["family"] for c in captures], "no_bundle_overlays_or_mixing")
    by_family = {c["manifest"]["family"]: c for c in captures}
    require(set(by_family) <= set(FAMILIES), "unknown_bundle_family")
    empty = capture(capture_timestamp=inputs["capture_timestamp"])
    qualified = tuple(
        qualify(by_family.get(f, empty), tuple(graphs.values()), family=f)
        for f in FAMILIES
    )
    rows = []
    names = {g.identity: name for name, g in graphs.items()}
    for q in qualified:
        for report in json.loads(q.graph_reports_json):
            rows.append(
                {
                    "motif": names[report["graph_identity"]],
                    "family": q.family,
                    "profile_id": q.profile_id,
                    "functional_class": q.functional_class,
                    **report,
                }
            )
    rows.sort(key=lambda r: (r["motif"], r["family"], r["profile_id"]))
    supplied = [q for q in qualified if q.bundle_identity is not None]
    incomplete = [
        q
        for q in supplied
        if q.decision != "authorized_bundle_qualified_runtime_blocked"
    ]
    decision = (
        "source_bundle_required"
        if not supplied
        else (
            "parameter_bundle_incomplete"
            if incomplete
            else "authorized_bundle_qualified_runtime_blocked"
        )
    )
    hashes = sorted({h for q in qualified for h in q.source_hashes})
    outputs = {
        "bundle-requirements.json": {
            "j84_schemas": schemas(),
            "additional_audited_convention": "energy_convention",
            "atom_type_references": "Every parameter type must name a declared source typing rule type.",
            "section_completeness": "All required final-graph instances; no generic section-presence acceptance.",
            "acceptance": "Exact declared graph coverage only; no runtime registration or coefficient acceptance.",
            "source_authority": "External manifest/audit pins and independent source review required; no discovery or fetching.",
        },
        "source-qualification.json": {
            "profiles": [
                {
                    k: v
                    for k, v in data(q).items()
                    if k not in ("graph_reports_json", "evidence_graph_json")
                }
                | {"replay_identity": q.replay_identity}
                for q in qualified
            ],
            "ranking": [
                {
                    "family": q.family,
                    "profile_id": q.profile_id,
                    "decision": q.decision,
                    "replay_identity": q.replay_identity,
                }
                for q in rank_qualifications(qualified)
            ],
        },
        "source-evidence-graph.json": {
            "source_hashes": hashes,
            "profiles": [
                {"family": q.family, "evidence": json.loads(q.evidence_graph_json)}
                for q in qualified
            ],
        },
        "organic-coverage-report.json": {
            "geometry_fixtures_only": True,
            "reports": rows,
            "scope": SCOPE,
        },
        "class-i-completeness.json": {
            "functional_class": "I",
            "reports": [r for r in rows if r["functional_class"] == "I"],
        },
        "class-ii-completeness.json": {
            "functional_class": "II",
            "reports": [r for r in rows if r["functional_class"] == "II"],
        },
        "result.json": {
            "decision": decision,
            "scope": SCOPE,
            "motifs": len(graphs),
            "profile_contracts": len(qualified),
            "candidate_decisions": len(rows),
            "authorized_bundles_supplied": sum(
                q.source_status == "audited_bundle_candidate"
                and not by_family[q.family]["manifest"]["synthetic"]
                for q in supplied
            ),
            "synthetic_bundles_supplied": sum(
                by_family[q.family]["manifest"]["synthetic"] for q in supplied
            ),
            "bundles_supplied": len(supplied),
            "qualified_profiles": sum(
                q.decision == "authorized_bundle_qualified_runtime_blocked"
                for q in qualified
            ),
            "qualified_graph_reports": sum(
                r["decision"] == "authorized_bundle_qualified_runtime_blocked"
                for r in rows
            ),
            "class_i_reports": sum(r["functional_class"] == "I" for r in rows),
            "class_ii_reports": sum(r["functional_class"] == "II" for r in rows),
            "accepted_coefficients": 0,
            "source_hashes": hashes,
            "replay_identity": hashlib.sha256(
                json.dumps(data(rows), sort_keys=True).encode()
            ).hexdigest(),
        },
    }
    return {name: data(signed(value)) for name, value in outputs.items()}


def emit(repo, root, bundle_capture=None):
    scope()
    require(pin(repo / FIXTURE) == FIXTURE_SHA, "fixture_pin_mismatch")
    require(not root.exists(), "do_not_overwrite_qualification_attempt")
    inputs = read(repo / FIXTURE)
    supplied = bundle_capture is not None and bundle_capture["manifest"] is not None
    if supplied:
        inputs["bundle_captures"] = [bundle_capture]
    reports = reconstruct(inputs)
    root.mkdir(parents=True)
    if not supplied:
        (root / "captured-input.json").write_bytes((repo / FIXTURE).read_bytes())
    else:
        write(root / "captured-input.json", inputs)
    for name, report in reports.items():
        write(root / name, report)
    # Every intake, including a blocked one, has a replay-bound receipt.
    write(
        root / "qualification-receipt.json",
        signed(
            {
                "schema": "island.j85-qualification-receipt.v1",
                "scope": SCOPE,
                "snapshot_sha256": pin(root / "captured-input.json"),
                "reports": {n: pin(root / n) for n in REPORTS},
                "decision": reports["result.json"]["decision"],
            }
        ),
    )
    return {
        "result": reports["result.json"],
        "snapshot_sha256": pin(root / "captured-input.json"),
    }


def replay(root, snapshot_sha=None):
    scope()
    sha = pin(root / "captured-input.json")
    require(sha == (snapshot_sha or FIXTURE_SHA), "captured_fixture_pin_mismatch")
    inputs = read(root / "captured-input.json")
    if snapshot_sha is None:
        require(not inputs["bundle_captures"], "source_claims_in_no_source_snapshot")
    outputs = reconstruct(inputs)
    for name, value in outputs.items():
        exact(read(root / name), value, "resigned_semantic_mutation:" + name)
    exact(
        read(root / "qualification-receipt.json"),
        signed(
            {
                "schema": "island.j85-qualification-receipt.v1",
                "scope": SCOPE,
                "snapshot_sha256": sha,
                "reports": {n: pin(root / n) for n in REPORTS},
                "decision": outputs["result.json"]["decision"],
            }
        ),
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
    p.add_argument("--manifest", type=Path)
    p.add_argument("--manifest-sha256")
    p.add_argument("--audit", type=Path)
    p.add_argument("--audit-sha256")
    p.add_argument("--snapshot-sha256")
    p.add_argument("--capture-timestamp")
    a = p.parse_args()
    try:
        scope()
        if a.command == "replay":
            value = replay(a.directory, a.snapshot_sha256)
        elif a.command == "ingest":
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
