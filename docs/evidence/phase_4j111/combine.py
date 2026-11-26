"""Combine the unchanged PCFF run and corrected OPLS run without replacing failures."""

import argparse
import gzip
import json
from pathlib import Path

from scripts.benchmark_final_graphs import summarize


def read(path):
    with gzip.open(path, "rt") as handle:
        return list(map(json.loads, handle))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pcff_run", type=Path)
    p.add_argument("opls_run", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--baseline",
        type=Path,
        default=Path("docs/evidence/phase_4j107b/paired-results.jsonl.gz"),
    )
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    records = []
    jobs = []
    for family, root in [("pcff", a.pcff_run), ("opls", a.opls_run)]:
        records.extend(
            r for r in read(root / "results.jsonl.gz") if r["field"] == family
        )
        jobs.extend(
            j
            for j in json.loads((root / "inputs.json").read_text())["jobs"]
            if j["field"] == family
        )
    records.sort(key=lambda r: (r["panel"], r["id"], r["field"]))
    summary = summarize(records, jobs)
    assert summary["complete"] and len(records) == 2032
    old = {(r["id"], r["field"]): r for r in read(a.baseline) if r["mode"] == "graph"}
    changes = []
    for r in records:
        previous = old[r["id"], r["field"]]
        if r["assigned"] != previous["assigned"]:
            changes.append(
                {
                    "id": r["id"],
                    "field": r["field"],
                    "panel": r["panel"],
                    "old_assigned": previous["assigned"],
                    "new_assigned": r["assigned"],
                    "failure": r.get("failure"),
                }
            )
    comparison = {
        "historical_rdkit_version": "not explicitly pinned in historical manifest",
        "current_benchmark_rdkit": "2024.09.6",
        "gains": sum(r["new_assigned"] for r in changes),
        "losses": sum(not r["new_assigned"] for r in changes),
        "records": changes,
        "type_count_changes": [
            {
                "id": r["id"],
                "field": r["field"],
                "old": old[r["id"], r["field"]]["type_counts"],
                "new": r["type_counts"],
            }
            for r in records
            if r["assigned"]
            and old[r["id"], r["field"]]["assigned"]
            and r["type_counts"] != old[r["id"], r["field"]]["type_counts"]
        ],
    }
    for name, value in [("summary.json", summary), ("comparison.json", comparison)]:
        (a.out / name).write_text(json.dumps(value, indent=2) + "\n")
    with gzip.open(a.out / "results.jsonl.gz", "wt") as handle:
        for r in records:
            handle.write(json.dumps(r) + "\n")
    print(json.dumps({k: comparison[k] for k in ("gains", "losses")}))


if __name__ == "__main__":
    main()
