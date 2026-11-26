"""Paired EMC repeat/final-graph coverage with bounded workers and complete denominators."""

import argparse
import csv
import gzip
import json
import os
import signal
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
from pathlib import Path
from tempfile import NamedTemporaryFile


def write(path, value):
    with NamedTemporaryFile(
        mode="w", dir=path.parent, prefix=".report-", delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def cases(panel):
    with gzip.open(panel, "rt") as handle:
        records = [r for r in csv.DictReader(handle) if r["field"] == "opls"]
    if len(records) != 982 or len({r["id"] for r in records}) != 982:
        raise ValueError("Expected the fixed 982 unique OPLS panel inputs")
    return [
        {
            "id": r["id"],
            "psmiles": r["psmiles"],
            "panel": "fixed982",
            "end_groups": [None, None],
        }
        for r in records
    ]


def challenges():
    repeats = {
        "PDMS": "*O[Si](C)(C)*",
        "phosphoester": "*OCCOP(=O)(O)*",
        "polythioether": "*CCS*",
        "polysulfone": "*c1ccc(S(=O)(=O)c2ccc(*)cc2)cc1",
        "PVC": "*CC(Cl)*",
        "PVBr": "*CC(Br)*",
        "PTFE": "*C(F)(F)C(F)(F)*",
        "sulfonated_PS": "*CC(c1ccc(S(=O)(=O)O)cc1)*",
        "PEG": "*CCO*",
        "PA6": "*NCCCCCC(=O)*",
        "PET": "*OCCOC(=O)c1ccc(C(=O)*)cc1",
        "PU": "*OCCOC(=O)NCCCCCCNC(=O)*",
        "PP_chiral": "*C[C@H](C)*",
        "cis_diene": r"*C/C=C\C*",
        "ammonium": "*CC[N+](C)(C)CC*",
    }
    rows = [
        {
            "id": name + "-" + cap,
            "panel": "challenge",
            "psmiles": repeat,
            "end_groups": ends,
        }
        for name, repeat in repeats.items()
        for cap, ends in [("H", [None, None]), ("methyl", ["*C", "*C"])]
    ]
    rows.append(
        {
            "id": "PDMS-TMS",
            "panel": "challenge",
            "psmiles": "*O[Si](C)(C)*",
            "end_groups": ["*[Si](C)(C)C", "*O[Si](C)(C)C"],
        }
    )
    for name, smi in [
        ("ammonium-chloride", "C[N+](C)(C)C.[Cl-]"),
        ("sodium-acetate", "CC(=O)[O-].[Na+]"),
        ("lithium-ether", "COCCOC.[Li+].[Cl-]"),
    ]:
        rows.append({"id": name, "panel": "challenge", "smiles": smi})
    return rows


def run_case(job, progress_path=None):
    from rdkit import Chem

    from island.forcefields.assignment import assign_polymer_parameters
    from island.forcefields.assignment_core import PolymerSpec, normalize
    from island.forcefields.assignment_diagnostics import (
        assignment_failure,
        preflight_assignment,
    )
    from island.forcefields.emc_graph import _section

    start = time.monotonic()
    result = {**job, "assigned": False, "strict": False, "md": "unmeasured"}
    graph = None
    try:
        try:
            if "smiles" in job:
                graph = Chem.MolFromSmiles(job["smiles"])
                if graph is None:
                    raise ValueError("Invalid input SMILES")
            else:
                graph = normalize(
                    PolymerSpec(
                        (("A", job["psmiles"]),),
                        job["dp"],
                        end_groups=tuple(job["end_groups"]),
                    )
                )
            report = preflight_assignment(
                graph, job["field"], field_dir=job["field_dir"]
            )
            result["preflight"] = report
            result["graph_valid"] = True
        except Exception as exc:
            result["graph_valid"] = False
            result["preflight_failure"] = assignment_failure(exc, graph)
            if job["mode"] == "graph":
                raise
        if progress_path is not None:
            write(progress_path, result)
        if job["mode"] == "graph":
            assignment = assign_polymer_parameters(
                graph,
                job["field"],
                field_dir=job["field_dir"],
                seed=job["seed"],
                backend_timeout=job.get("backend_timeout"),
            )
        else:
            assignment = assign_polymer_parameters(
                job["psmiles"],
                job["field"],
                dp=job["dp"],
                emc_mode="repeat",
                nchains=1,
                density=0.03,
                stereo_policy="polyjarvis",
                field_dir=job["field_dir"],
                seed=job["seed"],
            )
        result.update(
            assigned=True,
            strict=not getattr(assignment, "bond_stereo_omitted", False),
            field_sha256=assignment.field_sha256,
            emc_version=assignment.emc_version,
        )
        result["type_counts"] = dict(
            Counter(
                r.split()[-1] for r in _section(assignment.data, "Atoms").splitlines()
            )
        )
        result["artifact_sha256"] = {
            k: sha256(getattr(assignment, k).encode()).hexdigest()
            for k in ("data", "parameters")
        }
        if job["mode"] == "graph":
            result["input_smiles"] = assignment.input_smiles
        else:
            result["input_smiles"] = assignment.emc_smiles
    except Exception as exc:  # noqa: BLE001 - every failure remains in the denominator
        result.update(assigned=False, strict=False)
        result["failure"] = assignment_failure(exc, graph)
    result["seconds"] = time.monotonic() - start
    return result


def summarize(rows, jobs):
    key = lambda row: tuple(row.get(k) for k in ("panel", "id", "field", "mode"))
    expected, observed = [key(j) for j in jobs], [key(r) for r in rows]
    result = {
        "expected_jobs": len(jobs),
        "recorded_jobs": len(rows),
        "complete": len(set(observed)) == len(rows) == len(jobs)
        and set(observed) == set(expected),
        "groups": {},
    }
    for panel in sorted({j["panel"] for j in jobs}):
        for field in ("pcff", "opls"):
            for mode in ("repeat", "graph"):
                selected = [
                    r
                    for r in rows
                    if (r["panel"], r["field"], r["mode"]) == (panel, field, mode)
                ]
                if not selected:
                    continue
                result["groups"][f"{panel}:{field}:{mode}"] = {
                    "attempted": len(selected),
                    "assigned": sum(r["assigned"] for r in selected),
                    "strict": sum(r["strict"] for r in selected),
                    "md": "unmeasured",
                    "failure_categories": dict(
                        Counter(
                            c
                            for r in selected
                            for c in r.get("failure", {}).get("categories", [])
                        )
                    ),
                }
    result["paired"] = {}
    for field in ("pcff", "opls"):
        pairs = {}
        for row in rows:
            if row["panel"] == "fixed982" and row["field"] == field:
                pairs.setdefault(row["id"], {})[row["mode"]] = row
        result["paired"][field] = dict(
            Counter(
                (
                    "both"
                    if pair["graph"]["assigned"] and pair["repeat"]["assigned"]
                    else "graph_only"
                    if pair["graph"]["assigned"]
                    else "repeat_only"
                    if pair["repeat"]["assigned"]
                    else "neither"
                )
                for pair in pairs.values()
                if set(pair) == {"graph", "repeat"}
            )
        )
    result["element_groups"] = {}
    result["element_classification_missing"] = {}
    result["reported_missing_terms"] = {}
    result["localized_environment_gaps"] = {}
    for field in ("pcff", "opls"):
        selected = [
            r
            for r in rows
            if r["panel"] == "fixed982" and r["mode"] == "graph" and r["field"] == field
        ]
        result["element_classification_missing"][field] = sum(
            "preflight" not in r for r in selected
        )
        elements = sorted(
            {e for r in selected for e in r.get("preflight", {}).get("elements", {})}
        )
        result["element_groups"][field] = {
            e: {
                "attempted": sum(
                    e in r.get("preflight", {}).get("elements", {}) for r in selected
                ),
                "assigned": sum(
                    r["assigned"] and e in r.get("preflight", {}).get("elements", {})
                    for r in selected
                ),
            }
            for e in elements
        }
        failed = [r for r in selected if not r["assigned"]]
        result["reported_missing_terms"][field] = dict(
            Counter(
                t
                for r in failed
                for t in r.get("failure", {}).get("parameter_tuples", [])
            )
        )
        result["localized_environment_gaps"][field] = dict(
            Counter(
                e["fragment"]
                for r in failed
                for e in r.get("failure", {}).get("environments", [])
            )
        )
    return result


def execute(job, root, timeout):
    name = f"{job['panel']}-{job['id']}-{job['field']}-{job['mode']}"
    # IDs come from a pinned CSV, but still must not escape the output root.
    name = sha256(name.encode()).hexdigest()[:24]
    request, response, log = (
        root / (name + suffix) for suffix in (".input.json", ".result.json", ".log")
    )
    if response.exists():
        row = json.loads(response.read_text())
        if all(row.get(k) == v for k, v in job.items()):
            return row
        raise ValueError("Resume input changed")
    write(request, job)
    with log.open("w") as handle:
        proc = subprocess.Popen(
            [sys.executable, __file__, "--case", str(request)],
            stdout=handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            code = "timeout"
    if code == 0 and response.exists():
        return json.loads(response.read_text())
    partial_path = request.with_name(
        request.name.replace(".input.json", ".progress.json")
    )
    partial = json.loads(partial_path.read_text()) if partial_path.exists() else {}
    if any(partial.get(k, v) != v for k, v in job.items()):
        raise ValueError("Progress input changed")
    row = {
        **partial,
        **job,
        "assigned": False,
        "strict": False,
        "md": "unmeasured",
        "failure": {
            "categories": ["timeout" if code == "timeout" else "worker_error"],
            "diagnostic": log.read_text()[-6000:],
            "exit_code": code,
        },
    }
    write(response, row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path)
    parser.add_argument(
        "--panel", type=Path, default=Path("docs/evidence/reference_assignment.csv.gz")
    )
    parser.add_argument("--pcff-field", type=Path)
    parser.add_argument("--opls-field", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--backend-timeout", type=float, default=30)
    parser.add_argument("--challenge-only", action="store_true")
    args = parser.parse_args()
    if args.case:
        result = run_case(
            json.loads(args.case.read_text()),
            args.case.with_name(
                args.case.name.replace(".input.json", ".progress.json")
            ),
        )
        write(
            args.case.with_name(args.case.name.replace(".input.json", ".result.json")),
            result,
        )
        return
    if args.out is None or args.pcff_field is None or args.opls_field is None:
        parser.error("--out, --pcff-field and --opls-field are required")
    if args.workers < 1 or args.timeout <= 0:
        parser.error("workers and timeout must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    root = args.out / "cases"
    root.mkdir(exist_ok=True)
    inputs = ([] if args.challenge_only else cases(args.panel)) + challenges()
    jobs = []
    for row in inputs:
        for field, directory in [("pcff", args.pcff_field), ("opls", args.opls_field)]:
            for mode in (
                ["graph", "repeat"] if row["panel"] == "fixed982" else ["graph"]
            ):
                jobs.append(
                    {
                        **row,
                        "field": field,
                        "field_dir": str(directory.resolve()),
                        "mode": mode,
                        "dp": 4,
                        "backend_timeout": args.backend_timeout,
                        "seed": 2026,
                    }
                )
    from island.forcefields._assignment_cache import file_hash, tree_hash
    from island.forcefields.emc import FIELDS, _field_sources

    sources = {}
    for field, directory in [("pcff", args.pcff_field), ("opls", args.opls_field)]:
        runtime, binary, _, hashes = _field_sources(FIELDS[field], None, directory)
        sources[field] = {
            "field_sha256": hashes,
            "binary_sha256": file_hash(binary),
            "scripts_sha256": tree_hash(runtime / "scripts"),
        }
    manifest = {
        "panel_sha256": sha256(args.panel.read_bytes()).hexdigest(),
        "timeout_seconds": args.timeout,
        "sources": sources,
        "jobs": jobs,
    }
    path = args.out / "inputs.json"
    if path.exists() and json.loads(path.read_text()) != manifest:
        raise ValueError("Resume manifest changed")
    write(path, manifest)
    rows = []
    with ThreadPoolExecutor(args.workers) as pool:
        futures = [pool.submit(execute, j, root, args.timeout) for j in jobs]
        for future in as_completed(futures):
            rows.append(future.result())
            if len(rows) % 50 == 0:
                write(args.out / "summary.json", summarize(rows, jobs))
                print(
                    len(rows),
                    len(jobs),
                    "assigned",
                    sum(r["assigned"] for r in rows),
                    flush=True,
                )
    # Classification may evolve; reclassify every raw diagnostic together so
    # resumed records never mix old and new category definitions.
    from island.forcefields.assignment_diagnostics import assignment_failure

    for row in rows:
        failure = row.get("failure")
        if failure and "error_type" in failure:
            classified = assignment_failure(RuntimeError(failure["diagnostic"]))
            for key in ("categories", "parameter_tuples", "backend_sites"):
                failure[key] = classified[key]
    rows.sort(key=lambda r: (r["panel"], r["id"], r["field"], r["mode"]))
    with gzip.open(args.out / "results.jsonl.gz", "wt") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")
    write(args.out / "summary.json", summarize(rows, jobs))


if __name__ == "__main__":
    main()
