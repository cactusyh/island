"""Run current-code old-profile controls for every observed coverage change."""

import argparse
import gzip
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from rdkit import rdBase

from island.forcefields.emc_profiles import install
from scripts.benchmark_final_graphs import execute, write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--pcff-run", type=Path, required=True)
    p.add_argument("--opls-run", type=Path, required=True)
    p.add_argument(
        "--baseline",
        type=Path,
        default=Path("docs/evidence/phase_4j107b/paired-results.jsonl.gz"),
    )
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "cases").mkdir(exist_ok=True)
    for family, profile in [
        ("pcff", "pcff-carboxamide-v1"),
        ("opls", "opls-polymer-v2"),
    ]:
        destination = a.out / family
        if not destination.exists():
            install(profile, destination)
    with gzip.open(a.baseline, "rt") as handle:
        old = {
            (r["id"], r["field"]): r
            for line in handle
            if (r := json.loads(line))["mode"] == "graph"
        }
    rows = []
    for family, root in [("pcff", a.pcff_run), ("opls", a.opls_run)]:
        for path in sorted((root / "cases").glob("*.result.json")):
            r = json.loads(path.read_text())
            if r["field"] == family:
                rows.append(r)
    changed = [r for r in rows if r["assigned"] != old[r["id"], r["field"]]["assigned"]]
    jobs = []
    for r in changed:
        job = {
            k: r[k]
            for k in ("id", "panel", "field", "mode", "dp", "seed", "backend_timeout")
        }
        for k in ("psmiles", "end_groups", "smiles"):
            if k in r:
                job[k] = r[k]
        job["field_dir"] = str((a.out / r["field"]).resolve())
        jobs.append(job)
    with ThreadPoolExecutor(2) as pool:
        controls = list(pool.map(lambda job: execute(job, a.out / "cases", 180), jobs))
    report = []
    for r, control in zip(changed, controls, strict=True):
        report.append(
            {
                "id": r["id"],
                "panel": r["panel"],
                "field": r["field"],
                "historical_assigned": old[r["id"], r["field"]]["assigned"],
                "new_profile_assigned": r["assigned"],
                "current_old_profile_assigned": control["assigned"],
                "profile_gain": r["assigned"] and not control["assigned"],
                "profile_regression": not r["assigned"] and control["assigned"],
                "new_field_sha256": r.get("field_sha256"),
                "control": control,
            }
        )
    write(
        a.out / "controls.json",
        {
            "rdkit_version": rdBase.rdkitVersion,
            "observed_graph_records": len(rows),
            "expected_graph_records": 2032,
            "profile_gains": sum(r["profile_gain"] for r in report),
            "profile_regressions": sum(r["profile_regression"] for r in report),
            "records": report,
        },
    )
    print(len(rows), "records,", len(report), "controls")


if __name__ == "__main__":
    main()
