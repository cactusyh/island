"""Fenced, reproducible software coverage audit. Not a scientific parameter assignment."""

import argparse
import hashlib
import json
import subprocess
from dataclasses import fields as dataclass_fields
from pathlib import Path

from island.forcefields.assignment_core import normalize
from island.forcefields.assignment_core.profiles import (
    required_conventions,
    required_families,
)
from island.forcefields.assignment_core.schema import (
    EVIDENCE_STATUSES,
    READINESS,
    Branch,
    Crosslink,
    PolymerSpec,
    data,
)
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

BASE = "c5f7c14b71b7bbcc1192c57d1f5dc02b744c1aa0"
BRANCH = "codex/phase-4j83-evidence-aware-coverage-registry"
ROOT = Path("docs/evidence/phase_4j83")
FIXTURE = "tests/fixtures/phase_4j83/motifs.json"
FIXTURE_SHA = "2387550345235ff68caece83bf0d3d43b34848e8657e8f3cbbccaa65428ef7d9"
IMPLEMENTATION = (
    "src/island/forcefields/coverage_registry/__init__.py",
    "src/island/forcefields/coverage_registry/catalog.py",
    "src/island/forcefields/coverage_registry/registry.py",
    "scripts/phase_4j83_registry.py",
    "scripts/validate_phase_4j83_receipt.py",
    "tests/test_phase_4j83.py",
    FIXTURE,
    "docs/phase_4j83.md",
)
REPORTS = (
    "motif-catalog.json",
    "profile-registry.json",
    "coverage-report.json",
    "capability-matrix.json",
    "result.json",
)

SCOPE = {
    **READINESS,
    "J83_QM_attempts": 0,
    "J83_MD_steps": 0,
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
        inputs["schema"] == "island.j83.synthetic-motifs.v1"
        and inputs["scientific_validation"] is False,
        "scientific_promotion_of_controls",
    )
    require(
        len(inputs["motifs"]) == len({r["id"] for r in inputs["motifs"]}),
        "duplicate_motif_identity",
    )
    registry = ForceFieldProfileRegistry.contracts()
    profiles = registry.profiles
    graphs, reports, matrices = [], [], []
    for row in sorted(inputs["motifs"], key=lambda r: r["id"]):
        graph = normalize(spec_from_data(row))
        result = registry.coverage(graph)
        graphs.append(
            {"motif": row["id"], "graph": data(graph), "identity": graph.identity}
        )
        reports.append({"fixture": row["id"], "report": data(result)})
        matrices.append(
            {"fixture": row["id"], "matrix": data(result.capability_matrix)}
        )
    contracts = [
        {
            "profile": data(p),
            "identity": p.identity,
            "required_terms": required_families(p),
            "required_conventions": required_conventions(p),
        }
        for p in profiles
    ]
    output = {
        "motif-catalog.json": signed(
            {
                "catalog": data(registry.catalog),
                "identity": registry.catalog.identity,
                "graphs": graphs,
            }
        ),
        "profile-registry.json": signed(
            {
                "profiles": data(contracts),
                "identity": registry.identity,
                "evidence_statuses": EVIDENCE_STATUSES,
            }
        ),
        "coverage-report.json": signed(
            {"software_controls_only": True, "motifs": reports}
        ),
        "capability-matrix.json": signed(
            {
                "matrices": matrices,
                "overlap_policy": "motif_cells_overlap_global_unique_inventory_authoritative",
            }
        ),
        "result.json": signed(
            {
                "decision": "coverage_registry_implemented_runtime_blocked",
                "scope": SCOPE,
                "motifs": len(graphs),
                "profile_contracts": len(profiles),
                "candidate_decisions": len(graphs) * len(profiles),
                "blocked_runtime_decisions": len(reports),
                "scientific_parameter_rows_added": 0,
                "accepted_coefficients": 0,
                "synthetic_fixture_validation_only": True,
                "capability_cells": sum(len(m["matrix"]["cells"]) for m in matrices),
                "cell_decisions": {
                    k: sum(
                        c["coverage_decision"] == k
                        for m in matrices
                        for c in m["matrix"]["cells"]
                    )
                    for k in ("accepted", "provisional", "blocked")
                },
                "class_i_reports": sum(
                    r["functional_class"] == "I"
                    for v in reports
                    for r in v["report"]["profile_reports"]
                ),
                "class_ii_reports": sum(
                    r["functional_class"] == "II"
                    for v in reports
                    for r in v["report"]["profile_reports"]
                ),
                "replay_identity": digest(reports),
            }
        ),
    }
    return {k: data(v) for k, v in output.items()}


def emit(repo, root):
    scope()
    require(pin(repo / FIXTURE) == FIXTURE_SHA, "fixture_pin_mismatch")
    require(not root.exists(), "do_not_overwrite_coverage_attempt")
    inputs = read(repo / FIXTURE)
    reports = reconstruct(inputs)
    root.mkdir(parents=True)
    (root / "captured-input.json").write_bytes((repo / FIXTURE).read_bytes())
    for name, value in reports.items():
        write(root / name, value)
    return reports["result.json"]


def replay(root):
    scope()
    require(
        pin(root / "captured-input.json") == FIXTURE_SHA,
        "captured_fixture_pin_mismatch",
    )
    reports = reconstruct(read(root / "captured-input.json"))
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
    p.add_argument("mode", choices=("emit", "replay"))
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    a = p.parse_args()
    try:
        scope()
        print(
            json.dumps(
                emit(a.repo, a.directory) if a.mode == "emit" else replay(a.directory),
                sort_keys=True,
            )
        )
        return 0
    except (ValueError, OSError, TypeError, KeyError, RuntimeError) as e:
        print(json.dumps({"status": "blocked", "reason": str(e), "scope": SCOPE}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
