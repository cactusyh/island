"""Real DP10/50/100 two-seed packing and MD with complete resource/failure records."""

import argparse
import json
import os
import resource
import time
from pathlib import Path

from island.workflows.recipe import run_recipe


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    rows = []
    for dp in (10, 50, 100):
        for seed, chains in ((173, 1), (719, 2)):
            cfg = {
                "architecture": "homopolymer",
                "polymer": {"repeat_units": {"A": "*CC*"}, "dp": dp},
                "forcefield": {"name": "gaff2", "charge": "gasteiger"},
                "chains": chains,
                "seed": seed,
                "density_g_cm3": 0.1,
                "md": {
                    "nvt_steps": 20,
                    "npt_steps": 20,
                    "production_steps": 20,
                    "sample_every": 10,
                    "segment_steps": 20,
                },
            }
            name = f"dp{dp}-seed{seed}-n{chains}"
            row = {"name": name, "recipe": cfg}
            start = time.monotonic()
            try:
                row["result"] = run_recipe(
                    cfg,
                    a.out / name,
                    lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                    packmol=os.environ["ISLAND_PACKMOL"],
                )
                row["chains"] = json.loads(
                    (a.out / name / "packing/chains.json").read_text()
                )
                row["passed"] = True
            except Exception as exc:  # noqa: BLE001
                row.update(passed=False, error_type=type(exc).__name__, error=str(exc))
            row["seconds"] = time.monotonic() - start
            row["process_peak_rss_kib_cumulative"] = resource.getrusage(
                resource.RUSAGE_SELF
            ).ru_maxrss
            rows.append(row)
            (a.out / "checks.json").write_text(json.dumps(rows, indent=2) + "\n")
            print(name, row["passed"], row["seconds"], flush=True)
    if not all(r["passed"] for r in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
