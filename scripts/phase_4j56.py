"""Reproduce the J56 completion plan from pinned J54 bundles, without new fits."""

import argparse
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

from island.pcff_completion import (
    BASE,
    BRANCH,
    VERSION,
    classify,
    derive,
    load_assessment,
    pin,
    request_template,
    require,
    signed,
    write,
)
from island.pcff_completion.schema import SCHEMA

FIXTURES = ("disiloxane", "hmds", "dp2", "dp4", "dp6")


def audit(repo, distribution, output, *, fresh=True):
    repo, distribution, out = (
        Path(repo).resolve(),
        Path(distribution).resolve(),
        Path(output).resolve(),
    )
    require(not out.exists(), "output_directory_already_exists")
    out.mkdir(parents=True)
    original = repo / "docs/evidence/phase_4j54/audit"
    grouped, cases, dependencies = defaultdict(list), [], []
    relocated = out.parent / (out.name + "-relocated")
    if fresh:
        require(not relocated.exists(), "relocation_directory_already_exists")
        relocated.mkdir()
        shutil.copytree(distribution, relocated / "sources")
    for name in FIXTURES:
        bundle = original / name
        manifest_hash = pin(bundle / "manifest.json")
        s, a = load_assessment(bundle, manifest_hash, distribution)
        p = a.payload
        independent_file = original / (name + "-independent.json")
        independent = json.loads(independent_file.read_text())
        require(
            independent["status"] == "passed"
            and independent["raw_source_sha256"] == p["source_hashes"]["pcff.frc"]
            and all(
                c["coverage_complete"] and Decimal(c["residual"]) == 0
                for c in independent["component_totals"]
            ),
            "independent_charge_conservation_not_verified",
        )
        classification = classify(s, a)
        write(out / f"{name}.json", classification)
        for row in classification["terms"]:
            grouped[row["query"]].append({"fixture": name, **row})
        dependencies.append(
            {
                "fixture": name,
                "graph_identity": p["graph_identity"],
                "typing_identity": p["typing_identity"],
                "charge_identity": p["charge_identity"],
                "assessment_identity": p["identity"],
                "parameter_completeness": "incomplete",
                "missing_rows": [
                    {k: r[k] for k in ("id", "family", "query", "dependencies")}
                    for r in p["missing_terms"]
                ],
                "missing_conventions": p["missing_conventions"],
                "typing_status": p["typing_status"],
                "charge_status": p["charge_status"],
                "charge_component_totals": p["charge_record"]["components"],
                "independent_charge_conservation": {
                    "J54_record_sha256": pin(independent_file),
                    "component_totals": independent["component_totals"],
                    "status": "passed_exact_decimal_without_redistribution",
                },
                "source_equilibrium_dependencies_complete": all(
                    not r["missing_dependencies"] for r in classification["terms"]
                ),
                "class2_dependency_closure": "incomplete",
                "cross_engine_equivalence": "not_claimed",
            }
        )
        processes = []
        if fresh:
            shutil.copytree(bundle, relocated / name)
            for seed in (1, 99397):
                target = relocated / f"{name}-{seed}.json"
                cmd = [
                    sys.executable,
                    "-m",
                    "island.pcff_completion",
                    "classify",
                    "--bundle",
                    str(relocated / name),
                    "--manifest-sha256",
                    manifest_hash,
                    "--distribution",
                    str(relocated / "sources"),
                    "--json",
                    str(target),
                ]
                run = subprocess.run(
                    cmd,
                    env={**os.environ, "PYTHONHASHSEED": str(seed)},
                    capture_output=True,
                    text=True,
                    check=False,
                )
                require(
                    run.returncode == 0,
                    "fresh_process_failed",
                    stderr=run.stderr,
                    stdout=run.stdout,
                )
                require(
                    json.loads(target.read_text()) == classification,
                    "relocation_or_hash_randomization_mutation",
                )
                processes.append(
                    {
                        "hash_seed": seed,
                        "returncode": run.returncode,
                        "record_identity": classification["identity"],
                    }
                )
        cases.append(
            {
                "fixture": name,
                "atoms": s.number_of_sites,
                "manifest_sha256": manifest_hash,
                "classification_identity": classification["identity"],
                "fresh_processes": processes,
            }
        )
    old_queries = json.loads((original / "missing-rows.json").read_text())["queries"]
    require(
        set(grouped) == {q["query"] for q in old_queries},
        "changed_J54_query_membership",
    )
    queries = []
    for q in old_queries:
        rows = grouped[q["query"]]
        require(len(rows) == len(q["requests"]), "changed_J54_request_count")
        require(
            len({r["classification"] for r in rows}) == 1, "inconsistent_classification"
        )
        queries.append(
            {
                "query": q["query"],
                "family": q["family"],
                "source_status": "source_missing",
                "classification": rows[0]["classification"],
                "reason": rows[0]["reason"],
                "requests": len(rows),
                "fixtures": sorted({r["fixture"] for r in rows}),
                "equation_identity": rows[0]["equation_evidence"],
                "reference_status": "not_supplied",
                "coefficient_created": False,
            }
        )
    write(
        out / "47-query-classifications.json",
        signed(
            {
                "schema": VERSION,
                "queries": queries,
                "J54_missing_rows_sha256": pin(original / "missing-rows.json"),
            }
        ),
    )
    write(
        out / "dependency-closure.json",
        signed(
            {
                "schema": VERSION,
                "fixtures": dependencies,
                "J54_dependency_matrix_sha256": pin(
                    original / "class2-dependency-matrix.json"
                ),
                "class2_dependency_closure": "incomplete",
            }
        ),
    )
    write(
        out / "provenance-and-uncertainty.json",
        signed(
            {
                "schema": VERSION,
                "reference_datasets": [],
                "derived_terms": [],
                "fitted_terms": [],
                "uncertainty": [
                    {
                        "query": q["query"],
                        "status": "not_identifiable_without_reference_data",
                        "value": None,
                    }
                    for q in queries
                ],
                "derivation_rule": "energy functional form does not determine independent missing coefficients",
                "fitting_rule": "typing/charges, source equations, licenses, hashed geometry/reference data, identifiable residual, train/held-out coverage required",
                "not_extrapolated_to": [
                    "cyclic",
                    "crosslinked",
                    "Si-H outside exact supplied disiloxane graph",
                    "Si-OH",
                    "Si-O-C",
                    "phosphorus",
                ],
            }
        ),
    )
    write(
        out / "derivation-blocked-example.json",
        derive(
            *load_assessment(
                original / "hmds", pin(original / "hmds/manifest.json"), distribution
            ),
            next(
                r["id"] for r in grouped[queries[0]["query"]] if r["fixture"] == "hmds"
            ),
        ),
    )
    s, a = load_assessment(
        original / "hmds", pin(original / "hmds/manifest.json"), distribution
    )
    write(
        out / "unfilled-request-template.json",
        request_template(
            s,
            a,
            next(
                r["id"]
                for r in a.payload["missing_terms"]
                if r["family"] == "angle-angle"
            ),
        ),
    )
    write(out / "missing-term-fitting.schema.json", SCHEMA)
    summary = signed(
        {
            "schema": VERSION,
            "base_commit": BASE,
            "branch": BRANCH,
            "fixtures": cases,
            "directed_source_queries": len(queries),
            "missing_site_requests": sum(q["requests"] for q in queries),
            "query_classification_counts": dict(
                sorted(Counter(q["classification"] for q in queries).items())
            ),
            "source_hashes": a.payload["source_hashes"],
            "registered_profiles": [],
            "created_terms": [],
            "reference_dataset_status": "not_supplied",
            "parameter_completeness": "incomplete",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "MD_executed": False,
            "NPT_executed": False,
            "historical_source_fitting_assumptions_changed": False,
        }
    )
    write(out / "summary.json", summary)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--distribution", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    summary = audit(args.repo, args.distribution, args.output)
    print(
        json.dumps(
            {"identity": summary["identity"], "parameter_completeness": "incomplete"}
        )
    )


if __name__ == "__main__":
    main()
