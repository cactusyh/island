"""J90: offline OPLS provider intake, exact finite-oligomer reports and replay."""

import argparse
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path

from island.forcefields.assignment_core.schema import is_hash
from island.forcefields.bundle_adapters.schema import keys, loads
from island.forcefields.opls_providers import (
    FLAGS,
    OPLSProviderRegistry,
    benchmarks,
    finite_graph,
)
from island.forcefields.opls_providers.registry import pinned
from island.forcefields.opls_providers.schema import PROVIDERS, RESULTS
from scripts.phase_4j73_source_intake import exact, pin, read, require, signed, write
from scripts.phase_4j76_source_resolution import scope

BASE = "8a7933bd8088c23baeb75e9e9f3da61b4535dd67"
BRANCH = "codex/opls-provider-registry-provisional-polymers"
ROOT = Path("docs/evidence/phase_4j90")
PACKAGE = "src/island/forcefields/opls_providers/"
IMPLEMENTATION = (
    *(
        PACKAGE + n + ".py"
        for n in ("__init__", "schema", "polymer", "parsers", "registry", "export")
    ),
    "scripts/phase_4j90_providers.py",
    "scripts/validate_phase_4j90_receipt.py",
    "tests/test_phase_4j90.py",
    "docs/phase_4j90.md",
)
SCOPE = {
    "J90_QM_executions": 0,
    "J90_MD_steps": 0,
    "converter_executions": 0,
    "physical_fitting_executions": 0,
    "parameter_generation_executions": 0,
    "energy_force_evaluations": 0,
    "production_executions": 0,
    "reserved_validation_executions": 0,
    **FLAGS,
    "J72_historical_violation_and_C1_quarantine": "preserved_unchanged",
    "synthetic_controls_are_scientific_data": False,
}
REPORTS = (
    "graph-construction.json",
    "provider-results.json",
    "missing-terms-and-blockers.json",
    "result.json",
)
INPUTS = ("request", "manifest", "audit")


def provided_inputs(external):
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


def reconstruct(repo, **external):
    scope()
    provided = provided_inputs(external)
    registry = OPLSProviderRegistry()
    rows, graphs = [], []
    gate_blockers = []
    if external:
        for name in ("request", "manifest"):
            if not provided[name]["path_supplied"] or not is_hash(
                provided[name]["sha256"]
            ):
                gate_blockers.append("required_pinned_input_missing:" + name)
        if provided["audit"]["path_supplied"] != is_hash(provided["audit"]["sha256"]):
            gate_blockers.append("audit_path_and_pin_required")
    if external and not gate_blockers:
        request = loads(pinned(external["request"], external["request_sha256"]))
        keys(
            request,
            ("schema", "provider", "settings", "graph_identity"),
            "provider_request_schema",
        )
        require(
            request["schema"] == "island.opls-polymer-request.v1",
            "provider_request_schema",
        )
        graph = finite_graph(request["settings"])
        require(
            graph.identity == request["graph_identity"],
            "request_exact_graph_identity_mutation",
        )
        result = registry.assign(
            request["settings"],
            provider=request["provider"],
            **{k: v for k, v in external.items() if not k.startswith("request")},
        )
        rows.append({"benchmark": "explicit_user_graph", "report": result.report})
        graphs.append(
            {
                "benchmark": "explicit_user_graph",
                "settings": request["settings"],
                "graph_identity": graph.identity,
                "atom_count": len(graph.atom_order),
            }
        )
    else:
        for name, settings in benchmarks().items():
            graph = finite_graph(settings)
            graphs.append(
                {
                    "benchmark": name,
                    "settings": settings,
                    "graph_identity": graph.identity,
                    "atom_count": len(graph.atom_order),
                }
            )
            for provider in PROVIDERS:
                rows.append(
                    {
                        "benchmark": name,
                        "report": registry.assign(settings, provider=provider).report,
                    }
                )
    counts = {
        name: sum(r["report"]["decision"] == name for r in rows) for name in RESULTS
    }
    source_used = any(r["report"]["source_files_inspected"] for r in rows)
    decisions = Counter(r["report"]["decision"] for r in rows)
    # Benchmark aggregate never hides a blocked provider behind another provider.
    decision = next(
        (
            d
            for d in (
                "source_bundle_required",
                "parameter_bundle_incomplete",
                "provenance_only",
                "provisional_parameterization",
                "audited_exact_graph",
            )
            if decisions[d]
        ),
        "source_bundle_required",
    )
    snapshot = {
        "schema": "island.j90.external-input-pins.v1",
        "provided": provided,
        "external_sources_used": bool(source_used),
        "raw_source_records_persisted": False,
    }
    reports = {
        "graph-construction.json": signed(
            {
                "schema": "island.j90.graphs.v1",
                "graphs": graphs,
                "scope": "finite_capped_exact_graphs_only",
                "scientific_coverage_claimed_by_construction": False,
            }
        ),
        "provider-results.json": signed(
            {
                "schema": "island.j90.provider-reports.v1",
                "reports": rows,
                "provider_merging": False,
            }
        ),
        "missing-terms-and-blockers.json": signed(
            {
                "schema": "island.j90.blockers.v1",
                "input_blockers": gate_blockers,
                "reports": [
                    {
                        "benchmark": r["benchmark"],
                        "provider": r["report"]["provider"],
                        "missing_terms": r["report"]["missing_terms"],
                        "blocking_reasons": r["report"]["blocking_reasons"],
                    }
                    for r in rows
                ],
            }
        ),
        "result.json": signed(
            {
                "schema": "island.j90.result.v1",
                "decision": decision,
                "status_counts": counts,
                "declared_graphs": len(graphs),
                "provider_graph_reports": len(rows),
                "source_files_inspected": sum(
                    r["report"]["source_files_inspected"] for r in rows
                ),
                "source_rows_inspected": sum(
                    r["report"]["source_rows_inspected"] for r in rows
                ),
                "scientific_graphs_covered": sum(
                    r["report"]["scientific_graphs_covered"] for r in rows
                ),
                "exports_created": 0,
                "source_hashes": sorted(
                    {h for r in rows for h in r["report"]["source_hashes"]}
                ),
                "input_blockers": gate_blockers,
                "scope": SCOPE,
            }
        ),
    }
    return snapshot, reports


def intake_receipt(root, reports):
    return signed(
        {
            "schema": "island.j90.intake-receipt.v1",
            "base": BASE,
            "snapshot_sha256": pin(root / "captured-input.json"),
            "reports": {n: pin(root / n) for n in REPORTS},
            "decision": reports["result.json"]["decision"],
            "scope": SCOPE,
        }
    )


def emit(repo, root, **external):
    scope()
    require(not root.exists(), "do_not_overwrite_intake_attempt")
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
    inputs = read(root / "captured-input.json")
    if inputs["provided"] != provided_inputs({}):
        require(is_hash(snapshot_sha), "external_snapshot_pin_required")
    if snapshot_sha is not None:
        require(
            pin(root / "captured-input.json") == snapshot_sha, "snapshot_hash_mutation"
        )
    exact(
        inputs["provided"],
        provided_inputs(external),
        "explicit_external_originals_and_pins_required_for_replay",
    )
    fresh, reports = reconstruct(repo or Path.cwd(), **external)
    exact(inputs, fresh, "external_sources_unavailable_or_mutated")
    for name, report in reports.items():
        exact(read(root / name), report, "resigned_semantic_mutation:" + name)
    exact(
        read(root / "qualification-receipt.json"),
        intake_receipt(root, reports),
        "intake_receipt_mutation",
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
    a = p.parse_args()
    external = {
        k: getattr(a, k)
        for n in INPUTS
        for k in (n, n + "_sha256")
        if getattr(a, k) is not None
    }
    try:
        scope()
        value = (
            replay(a.directory, a.snapshot_sha256, repo=a.repo, **external)
            if a.command == "replay"
            else emit(a.repo, a.directory, **external)
        )
        print(json.dumps(value, sort_keys=True))
        return 0
    except (ValueError, TypeError, KeyError, OSError, RuntimeError) as error:
        print(
            json.dumps(
                {
                    "status": "blocked",
                    "decision": "parameter_bundle_incomplete",
                    "blocking_reasons": [str(error)],
                    "scope": SCOPE,
                },
                sort_keys=True,
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
