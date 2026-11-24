"""Run a reproducible architecture/force-field matrix through the JSON recipe CLI."""

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--lammps", default="lmp")
    parser.add_argument("--packmol", default="packmol")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    recipes = Path(__file__).resolve().parents[1] / "examples/recipes"
    cases = []
    for base in ("pe-gaff2", "block-pcff", "branched-opls"):
        for ff in ("gaff", "gaff2", "pcff", "opls"):
            value = json.loads((recipes / (base + ".json")).read_text())
            value["forcefield"] = {"name": ff}
            if ff.startswith("gaff"):
                value["forcefield"]["charge"] = "gasteiger"
            elif ff == "opls" and base == "block-pcff":
                value["forcefield"]["profile"] = "opls-epoxy-amine-v1"
            cases.append((value["architecture"] + "-" + ff, value))
    for kind in ("pcff", "opls", "thiol"):
        value = json.loads((recipes / "epoxy-pcff.json").read_text())
        if kind == "opls":
            value["forcefield"] = {"name": "opls", "profile": "opls-epoxy-amine-v1"}
        if kind == "thiol":
            value["seed"] = 173
            value["cure"].pop("pacm")
            value["cure"].update(petmp=2, target=1.0)
        cases.append(("crosslinked-" + kind, value))

    def run(case):
        name, value = case
        input_path = args.directory / (name + ".json")
        input_path.write_text(json.dumps(value, indent=2) + "\n")
        output = args.directory / name
        command = [
            sys.executable,
            "-m",
            "island.workflows.recipe",
            str(input_path),
            str(output),
            "--lammps",
            args.lammps,
            "--packmol",
            args.packmol,
        ]
        with (args.directory / (name + ".log")).open("w") as log:
            proc = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        if proc.returncode:
            return {
                "case": name,
                "status": "failed",
                "diagnostic": (args.directory / (name + ".log")).read_text()[-4000:],
            }
        result = json.loads((output / "result.json").read_text())
        return {
            "case": name,
            "status": result["status"],
            "field": result["assignment"]["field"],
            "stages": result["stages"],
            "cure": result.get("cure"),
            "target_reached": result.get("target_reached"),
        }

    with ThreadPoolExecutor(args.workers) as pool:
        results = list(pool.map(run, cases))
    (args.directory / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    for row in results:
        print(row["case"], row["status"])
    return int(any(r["status"] != "completed" for r in results))


if __name__ == "__main__":
    raise SystemExit(main())
