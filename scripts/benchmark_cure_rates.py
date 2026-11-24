"""Matched-conversion temperature and strain-rate sensitivity on two box sizes."""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from island.epoxy.properties import conditions, cooling, tension


def run(job):
    fun, source, dst, lmp, args = job
    result = fun(source, dst, lmp, *args)
    (dst / "result.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--lammps", default="lmp")
    p.add_argument("--workers", type=int, default=1)
    a = p.parse_args()
    case = a.case
    r = json.loads((case / "measurement.json").read_text())
    if r["statistics"]["epoxy_conversion"] != 0.75:
        raise ValueError("Requires a 75%-cured sample")
    if "conditioning_failure" in r:
        raise ValueError("Rate comparisons require a completed 300 K conditioning run")
    a.out.mkdir(parents=True, exist_ok=False)
    run((conditions, case / "after", a.out / "hot", a.lammps, (550, 40000)))
    jobs = [
        (cooling, a.out / "hot", a.out / f"cool-{rate}", a.lammps, (rate,))
        for rate in (4, 8)
    ]
    jobs += [
        (tension, case / "after", a.out / f"pull-{rate}-{axis}", a.lammps, (rate, axis))
        for rate in (0.001, 0.002)
        for axis in "xyz"
    ]
    with ProcessPoolExecutor(a.workers) as pool:
        results = list(pool.map(run, jobs))
    (a.out / "fits.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
