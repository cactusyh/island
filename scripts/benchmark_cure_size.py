"""Matched conversion finite-size and sampling benchmark for the PCFF cure."""

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np

from island.epoxy.cure import CureConfig, run_cure
from island.epoxy.properties import conditions


def case(job, *, out, library, lammps, packmol, steps):
    n, seed = job
    path = out / f"d{n}-p{n // 2}-s{seed}"
    start = time.monotonic()
    cfg = CureConfig(
        dgeba=n, pacm=n // 2, seed=seed, target=0.75, cutoff=5.5, max_idle=6
    )
    rec = run_cure(path / "cure", library, cfg, lammps=lammps, packmol=packmol)
    stat = rec[-1]["statistics"]
    result = {
        "dgeba": n,
        "pacm": n // 2,
        "seed": seed,
        "statistics": stat,
        "termination": rec[-1]["reason"],
    }
    if stat["epoxy_conversion"] == cfg.target:
        result["volumes"] = {}
        for label, idx in [("before", 0), ("after", rec[-1]["index"])]:
            stage = path / "cure" / f"stage-{idx:04d}"
            dest = path / label
            try:
                conditions(stage, dest, lammps, 300, steps)
            except RuntimeError as exc:
                result["conditioning_failure"] = {
                    "state": label,
                    "diagnostic": str(exc),
                }
                break
            values = np.loadtxt(dest / "samples.dat")
            result["volumes"][label] = {
                str(k): {
                    "volume_A3": float(values[k // 200 : k // 100, 3].mean()),
                    "density_g_cm3": float(values[k // 200 : k // 100, 2].mean()),
                }
                for k in [12000, steps]
            }
    result["elapsed_seconds"] = time.monotonic() - start
    (path / "measurement.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--library", type=Path, required=True)
    p.add_argument("--lammps", default="lmp")
    p.add_argument("--packmol", default="packmol")
    p.add_argument("--sizes", type=int, nargs="+", default=[4, 8, 16])
    p.add_argument("--seeds", type=int, nargs="+", default=[41, 173, 2026, 719])
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--steps", type=int, default=40000)
    a = p.parse_args()
    if any(n <= 0 or n % 4 for n in a.sizes) or a.steps <= 12000 or a.steps % 200:
        p.error(
            "sizes must be positive multiples of four and steps a multiple of 200 above 12000"
        )
    a.out.mkdir(parents=True, exist_ok=False)
    with ProcessPoolExecutor(a.workers) as pool:
        res = list(
            pool.map(
                partial(
                    case,
                    out=a.out,
                    library=a.library,
                    lammps=a.lammps,
                    packmol=a.packmol,
                    steps=a.steps,
                ),
                [(n, s) for n in a.sizes for s in a.seeds],
            )
        )
    (a.out / "cases.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
