"""Fenced source qualification and organic coverage. No runtime assignment."""

import argparse
import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from island.forcefields.assignment_core import normalize
from island.forcefields.assignment_core.schema import READINESS, data, unique
from island.forcefields.bundle_adapters import capture
from island.forcefields.bundle_adapters.scoped_qualification import (
    QualificationRequest,
    aggregate_qualifications,
    qualify_request,
    rank_qualifications,
)
from scripts.phase_4j73_source_intake import exact, pin, read, require, signed, write
from scripts.phase_4j76_source_resolution import scope
from scripts.phase_4j84_bundles import spec_from_data

BASE = "19fa2e76cd822987618cf03e42fca9dd9d884975"
BRANCH = "codex/phase-4j86-scoped-bundle-qualification"
ROOT = Path("docs/evidence/phase_4j86")
FIXTURE = "tests/fixtures/phase_4j86/qualification-set.json"
FIXTURE_SHA = "95f7e3e3957efa531c13a2acf38c0a2e4b15279fe0ec08773a9febe02b3f1d5c"
IMPLEMENTATION = (
    "src/island/forcefields/bundle_adapters/scoped_qualification.py",
    "scripts/phase_4j86_qualification.py",
    "scripts/validate_phase_4j86_receipt.py",
    "tests/test_phase_4j86.py",
    FIXTURE,
    "docs/phase_4j86.md",
)
REPORTS = (
    "qualification-set.json",
    "source-qualification.json",
    "source-evidence-graph.json",
    "graph-coverage-report.json",
    "class-i-completeness.json",
    "class-ii-completeness.json",
    "result.json",
)
SCOPE = {
    **READINESS,
    "J86_QM_attempts": 0,
    "J86_MD_steps": 0,
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
    require(
        inputs["schema"] == "island.j86.scoped-qualification.v1", "qualification_schema"
    )
    require(inputs["geometry_fixtures_only"] is True, "geometry_fixture_promotion")
    unique([r["id"] for r in inputs["motifs"]], "duplicate_fixture_identity")
    graphs = {r["id"]: normalize(spec_from_data(r)) for r in inputs["motifs"]}
    require(graphs, "qualification_requires_graphs")
    unique([g.identity for g in graphs.values()], "duplicate_fixture_graph")
    requests = []
    for r in inputs["required_profiles"]:
        require(
            set(r) == {"family", "profile_id", "graphs"}, "profile_declaration_fields"
        )
        require(set(r["graphs"]) <= graphs.keys(), "undeclared_graph_fixture")
        unique(r["graphs"], "duplicate_requested_graph")
        requests.append(
            QualificationRequest(
                r["family"],
                r["profile_id"],
                tuple(graphs[g].identity for g in r["graphs"]),
            )
        )
    require(
        set(graphs) == {g for r in inputs["required_profiles"] for g in r["graphs"]},
        "unused_graph_fixture",
    )
    captures = inputs["bundle_captures"]
    unique(
        [(c["manifest"]["family"], c["manifest"]["profile_id"]) for c in captures],
        "no_bundle_overlays_or_mixing",
    )
    by_profile = {
        (c["manifest"]["family"], c["manifest"]["profile_id"]): c for c in captures
    }
    require(set(by_profile) <= {r.key for r in requests}, "undeclared_bundle_profile")
    empty = capture(capture_timestamp=inputs["capture_timestamp"])
    by_identity = {g.identity: g for g in graphs.values()}
    scoped = tuple(
        qualify_request(
            by_profile.get(r.key, empty),
            tuple(by_identity[g] for g in r.graph_identities),
            request=r,
        )
        for r in requests
    )
    aggregate = aggregate_qualifications(scoped, qualification_set=requests)
    qualified = tuple(r.qualification for r in aggregate.profiles)
    rows = []
    names = {g.identity: name for name, g in graphs.items()}
    for result in aggregate.profiles:
        q = result.qualification
        for report in json.loads(q.graph_reports_json):
            rows.append(
                {
                    "fixture": names[report["graph_identity"]],
                    "family": q.family,
                    "declared_profile_id": result.request.profile_id,
                    "source_profile_id": q.profile_id,
                    "functional_class": q.functional_class,
                    **report,
                }
            )
    rows.sort(key=lambda r: (r["fixture"], r["family"], r["declared_profile_id"]))
    hashes = sorted({h for q in qualified for h in q.source_hashes})
    family_rankings = {}
    for family in sorted({q.family for q in qualified}):
        # Different exact scopes are intentionally incomparable, even within a family.
        groups = {}
        for q in qualified:
            if q.family == family and q.bundle_identity is not None:
                scope_id = tuple(
                    sorted(
                        r["graph_identity"] for r in json.loads(q.graph_reports_json)
                    )
                )
                groups.setdefault(scope_id, []).append(q)
        family_rankings[family] = [
            {
                "graph_identities": scope_id,
                "candidates": [
                    {
                        "family": q.family,
                        "profile_id": q.profile_id,
                        "decision": q.decision,
                        "replay_identity": q.replay_identity,
                    }
                    for q in rank_qualifications(
                        groups[scope_id], requested_family=family
                    )
                ],
            }
            for scope_id in sorted(groups)
        ]
    outputs = {
        "qualification-set.json": {
            "declarations": data(aggregate.qualification_set),
            "scope": "declared_exact_graphs_only",
            "chemical_family_transfer": "not_established",
        },
        "source-qualification.json": {
            "aggregate": data(aggregate),
            "replay_identity": aggregate.replay_identity,
            "rankings_by_requested_family_and_graph_scope": family_rankings,
        },
        "source-evidence-graph.json": {
            "source_hashes": hashes,
            "profiles": [
                {
                    "family": q.family,
                    "profile_id": q.profile_id,
                    "evidence": json.loads(q.evidence_graph_json),
                }
                for q in qualified
            ],
        },
        "graph-coverage-report.json": {
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
            "decision": aggregate.decision,
            "scope": SCOPE,
            "scope_of_qualification": "declared_exact_graphs_only",
            "chemical_family_transfer": "not_established",
            "graph_fixtures": len(graphs),
            "required_profiles": len(requests),
            "qualified_profiles": data(aggregate.qualified_profiles),
            "missing_or_unauthorized_profiles": data(
                aggregate.missing_or_unauthorized_profiles
            ),
            "incomplete_profiles": data(aggregate.incomplete_profiles),
            "per_profile_decisions": [
                {
                    "family": r.request.family,
                    "profile_id": r.request.profile_id,
                    "decision": r.qualification.decision,
                }
                for r in aggregate.profiles
            ],
            "bundles_supplied": len(captures),
            "authorized_bundles_supplied": sum(q.source_authorized for q in qualified),
            "source_hashes": hashes,
            "accepted_coefficients": 0,
            "replay_identity": aggregate.replay_identity,
        },
    }
    return {name: data(signed(value)) for name, value in outputs.items()}


def emit(repo, root, bundle_capture=None, declaration=None, declaration_sha=None):
    scope()
    require(pin(repo / FIXTURE) == FIXTURE_SHA, "fixture_pin_mismatch")
    require(not root.exists(), "do_not_overwrite_qualification_attempt")
    if declaration is not None:
        require(
            declaration_sha and pin(declaration) == declaration_sha,
            "qualification_declaration_pin_mismatch",
        )
        inputs = read(declaration)
        require(
            not inputs["bundle_captures"], "declaration_must_not_embed_source_captures"
        )
    else:
        require(declaration_sha is None, "declaration_path_required")
        inputs = read(repo / FIXTURE)
    supplied = bundle_capture is not None and bundle_capture["manifest"] is not None
    if supplied:
        inputs["bundle_captures"] = [bundle_capture]
    reports = reconstruct(inputs)
    root.mkdir(parents=True)
    if not supplied and declaration is None:
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
                "schema": "island.j86-qualification-receipt.v1",
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
                "schema": "island.j86-qualification-receipt.v1",
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
    p.add_argument("--qualification-set", type=Path)
    p.add_argument("--qualification-set-sha256")
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
            value = emit(
                a.repo,
                a.directory,
                cap,
                a.qualification_set,
                a.qualification_set_sha256,
            )
        else:
            require(
                a.manifest is None and a.audit is None,
                "use_ingest_for_explicit_sources",
            )
            value = emit(
                a.repo,
                a.directory,
                declaration=a.qualification_set,
                declaration_sha=a.qualification_set_sha256,
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
