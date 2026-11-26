"""Fixed 982 complete-graph parity with independently invoked native RadonPy.

This graph-only panel measures typing, coefficients and Gasteiger charges, not
conformer generation, packing, dynamics, or quantum-derived charge methods.
"""

import argparse
import gzip
import json
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from pathlib import Path

from scripts.benchmark_final_graphs import cases


def execute(job):
    from radonpy.core import calc, utils
    from radonpy.ff.gaff import GAFF
    from radonpy.ff.gaff2 import GAFF2
    from rdkit import Chem

    from island.forcefields.assignment_core import PolymerSpec, normalize
    from island.forcefields.assignment_core.topology import molecule
    from island.forcefields.radonpy import assign_radonpy

    row = {
        **job,
        "graph_valid": False,
        "island_assigned": False,
        "native_assigned": False,
        "parameter_match": False,
    }
    try:
        mol = molecule(normalize(PolymerSpec((("A", job["psmiles"]),), 4)))
        Chem.SanitizeMol(mol)
        row["graph_valid"] = True
    except Exception as exc:  # noqa: BLE001
        row["graph_error"] = str(exc)
        return row
    outputs = {}
    for backend in ("island", "native"):
        try:
            if backend == "island":
                assigned = assign_radonpy(mol, job["family"], charge="gasteiger")
                row["database_sha256"] = assigned.database_sha256
                product = assigned.mol
            else:
                product = Chem.Mol(mol)
                engine = {"gaff": GAFF, "gaff2": GAFF2}[job["family"]]()
                if not engine.ff_assign(product, useMDL=True, retryMDL=True):
                    raise ValueError("native typing/parameter assignment failed")
                if not calc.assign_charges(product, charge="gasteiger"):
                    raise ValueError("native Gasteiger charges failed")
            outputs[backend] = utils.MolToJSON_dict(product)
            row[backend + "_assigned"] = True
        except Exception as exc:  # noqa: BLE001
            row[backend + "_error"] = str(exc)
    if len(outputs) == 2:
        row["parameter_match"] = outputs["island"] == outputs["native"]
        row["parameters_sha256"] = {
            k: sha256(json.dumps(v, sort_keys=True).encode()).hexdigest()
            for k, v in outputs.items()
        }
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--panel", type=Path, default=Path("docs/evidence/reference_assignment.csv.gz")
    )
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    jobs = [
        {"id": c["id"], "psmiles": c["psmiles"], "family": f}
        for c in cases(a.panel)
        for f in ("gaff", "gaff2")
    ]
    rows = []
    with (
        ProcessPoolExecutor(a.workers) as pool,
        gzip.open(a.out / "results.jsonl.gz", "wt") as handle,
    ):
        for row in pool.map(execute, jobs):
            rows.append(row)
            handle.write(json.dumps(row) + "\n")
            handle.flush()
    summary = {
        "panel_sha256": sha256(a.panel.read_bytes()).hexdigest(),
        "dp": 4,
        "charge": "gasteiger",
        "groups": {},
    }
    for f in ("gaff", "gaff2"):
        selected = [r for r in rows if r["family"] == f]
        summary["groups"][f] = {
            "attempted": len(selected),
            **{
                k: sum(r[k] for r in selected)
                for k in (
                    "graph_valid",
                    "island_assigned",
                    "native_assigned",
                    "parameter_match",
                )
            },
            "status_mismatches": [
                r["id"]
                for r in selected
                if r["native_assigned"] != r["island_assigned"]
            ],
        }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
