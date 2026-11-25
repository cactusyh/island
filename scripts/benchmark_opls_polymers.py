"""Replay the fixed 982-polymer OPLS panel with an explicitly selected field."""

import argparse
import csv
import gzip
import json
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from island.exceptions import ForceFieldError
from island.forcefields.assignment import assign_polymer_parameters


def check(row, field):
    result = {key: row[key] for key in ("id", "psmiles")}
    result["stock_assigned"] = row["assigned"] == "True"
    try:
        assignment = assign_polymer_parameters(
            row["psmiles"],
            "opls",
            dp=4,
            emc_mode="repeat",
            nchains=2,
            density=0.9,
            seed=2026,
            stereo_policy="polyjarvis",
            field_dir=field,
        )
        result.update(
            assigned=True,
            strict=not assignment.bond_stereo_omitted,
            field_sha256=assignment.field_sha256,
        )
    except (ValueError, ForceFieldError) as exc:
        result.update(assigned=False, diagnostic=str(exc))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", type=Path)
    parser.add_argument("field_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    with gzip.open(args.panel, "rt") as stream:
        records = [r for r in csv.DictReader(stream) if r["field"] == "opls"]
    with (
        ProcessPoolExecutor(args.workers) as pool,
        gzip.open(args.output, "xt") as output,
    ):
        for result in pool.map(partial(check, field=args.field_dir), records):
            output.write(json.dumps(result) + "\n")
            output.flush()


if __name__ == "__main__":
    main()
