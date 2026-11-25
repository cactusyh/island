"""Replay a pinned PCFF panel against an independently installed template overlay."""

import argparse
import csv
import gzip
import json
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from island.exceptions import ForceFieldError
from island.forcefields.assignment import assign_polymer_parameters


def evaluate(r, field_dir):
    report = {k: r[k] for k in ["id", "psmiles"]}
    report["stock_assigned"] = r["assigned"] == "True"
    try:
        a = assign_polymer_parameters(
            r["psmiles"],
            "pcff",
            field_dir=field_dir,
            dp=4,
            emc_mode="repeat",
            nchains=2,
            density=0.9,
            seed=2026,
            stereo_policy="polyjarvis",
        )
        report.update(
            assigned=True, strict=not a.bond_stereo_omitted, field_sha256=a.field_sha256
        )
    except (ValueError, ForceFieldError) as e:
        report.update(assigned=False, diagnostic=str(e))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", type=Path)
    p.add_argument("field_dir", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--workers", type=int, default=1)
    a = p.parse_args()
    with gzip.open(a.panel, "rt") as f:
        rows = [r for r in csv.DictReader(f) if r["field"] == "pcff"]
    with ProcessPoolExecutor(a.workers) as pool, gzip.open(a.out, "xt") as out:
        for r in pool.map(partial(evaluate, field_dir=a.field_dir), rows):
            out.write(json.dumps(r) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
