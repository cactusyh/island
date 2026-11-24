"""Measure graph, assignment and periodic-MD coverage over a fixed chemistry panel."""

import argparse
import json
import subprocess
import sys
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_core.topology import molecule

REPEATS = {
    "PE": "*CC*",
    "PS": "*CC(c1ccccc1)*",
    "PEG": "*CCO*",
    "PCL": "*CCCCCC(=O)O*",
    "PA6": "*NCCCCCC(=O)*",
    "PMMA": "*CC(C)(C(=O)OC)*",
}

STAGES = ("graph", "packed", "assignment", "single_point", "md")


def summarize(records, expected):
    """Keep every attempted case in the denominator, including early failures."""
    groups = {"all": records}
    for key in ("chemistry", "architecture"):
        for value in sorted({r[key] for r in records}):
            groups[f"{key}:{value}"] = [r for r in records if r[key] == value]
    for field in sorted({r["field"]["name"] for r in records}):
        groups[f"field:{field}"] = [r for r in records if r["field"]["name"] == field]
    return {
        "expected_cases": expected,
        "recorded_cases": len(records),
        "complete": len(records) == expected,
        "equilibration": "unmeasured",
        "groups": {
            name: {
                "attempted": len(rows),
                **{stage: sum(bool(row[stage]) for row in rows) for stage in STAGES},
            }
            for name, rows in groups.items()
        },
    }


def finite_single_point(paths):
    for path in paths:
        try:
            if np.isfinite(float(path.read_text())):
                return True
        except (ValueError, OSError):
            continue
    return False


def cases():
    records = []
    for name, psmiles in REPEATS.items():
        for architecture in ("homopolymer", "block", "branched"):
            polymer = {"repeat_units": {"A": psmiles}, "dp": 6}
            if architecture == "block":
                polymer = {
                    "repeat_units": {
                        "A": psmiles,
                        "B": REPEATS["PEG" if name == "PE" else "PE"],
                    },
                    "blocks": [["A", 3], ["B", 3]],
                }
            if architecture == "branched":
                mol = molecule(normalize(PolymerSpec((("A", psmiles),), 6)))
                eligible = [
                    a.GetIdx() + 1
                    for a in mol.GetAtoms()
                    if a.GetAtomicNum() == 6
                    and not a.GetIsAromatic()
                    and sum(n.GetAtomicNum() > 1 for n in a.GetNeighbors()) == 2
                    and any(n.GetAtomicNum() == 1 for n in a.GetNeighbors())
                ]
                if not eligible:
                    raise ValueError(f"No replaceable backbone C-H site for {name}")
                polymer["branches"] = [
                    {"anchor_site": eligible[len(eligible) // 2], "psmiles": "*CC"}
                ]
            for field in ("gaff2", "pcff", "opls"):
                ff = {"name": field}
                if field == "gaff2":
                    ff["charge"] = "gasteiger"
                else:
                    ff["profile"] = (
                        "pcff-carboxamide-v1" if field == "pcff" else "opls-polymer-v1"
                    )
                records.append(
                    {
                        "name": f"{name}-{architecture}-{field}",
                        "chemistry": name,
                        "recipe": {
                            "architecture": architecture,
                            "polymer": polymer,
                            "forcefield": ff,
                            "chains": 4,
                            "density_g_cm3": 0.6,
                            "seed": 173,
                            "md": {
                                "nvt_steps": 1000,
                                "npt_steps": 8000,
                                "production_steps": 4000,
                                "sample_every": 200,
                                "segment_steps": 2000,
                            },
                        },
                    }
                )
    return records


def execute(case, root, lammps, packmol):
    name = case["name"]
    recipe_path, output = root / (name + ".json"), root / name
    recipe_path.write_text(json.dumps(case["recipe"], indent=2) + "\n")
    with (root / (name + ".log")).open("w") as log:
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "island.workflows.recipe",
                str(recipe_path),
                str(output),
                "--lammps",
                lammps,
                "--packmol",
                packmol,
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    graph = (output / "packing/chain.json").is_file() or bool(
        list(output.glob("failed-packing-*/chain.json"))
        or list(output.glob("failed-packing-*/attempt-*/chain.json"))
    )
    assignment = (output / "assignment/checkpoint.json").is_file()
    run0 = list(output.glob("*/start.txt"))
    single_point = finite_single_point(run0)
    record = {
        "name": name,
        "chemistry": case["chemistry"],
        "architecture": case["recipe"]["architecture"],
        "field": case["recipe"]["forcefield"],
        "density_g_cm3": case["recipe"]["density_g_cm3"],
        "graph": graph,
        "packed": (output / "packing/checkpoint.json").is_file(),
        "assignment": assignment,
        "single_point": single_point,
        "md": proc.returncode == 0,
        "equilibration": "unmeasured",
    }
    attempts = output / "packing/packing-attempts.json"
    if attempts.is_file():
        record["packing_attempts"] = json.loads(attempts.read_text())
    if proc.returncode == 0:
        result = json.loads((output / "result.json").read_text())
        record["stages"] = result["stages"]
    else:
        status = output / "status.json"
        record["failure"] = (
            json.loads(status.read_text())
            if status.is_file()
            else {"diagnostic": (root / (name + ".log")).read_text()[-4000:]}
        )
    return record


def save_reports(root, records, expected):
    ordered = sorted(records, key=lambda r: r["name"])
    for filename, value in (
        ("results.json", ordered),
        ("summary.json", summarize(ordered, expected)),
    ):
        temporary = root / (filename + ".tmp")
        temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
        temporary.replace(root / filename)


def archive_runs(root, destination):
    """Retain independent stage replay inputs and four representative trajectories."""
    names = {
        "panel.json",
        "results.json",
        "summary.json",
        "recipe.json",
        "runtime.json",
        "status.json",
        "system.data",
        "system.params",
        "assignment.json",
        "md.in",
        "screen.log",
        "end.txt",
        "start.txt",
        "final.dump",
        "start.restart",
        "pack.log",
        "pack.in",
        "monomer0.xyz",
        "packed.xyz",
        "packing-attempts.json",
        "chain.json",
        "site-map.json",
        "result.json",
        "hashes.json",
        "profile.json",
    }
    trajectories = {
        "PE-homopolymer-gaff2",
        "PEG-branched-opls",
        "PCL-homopolymer-opls",
        "PA6-homopolymer-pcff",
    }
    hashes = {}
    with zipfile.ZipFile(destination, "x", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(root.rglob("*")):
            relative = path.relative_to(root)
            if path.is_file() and (
                path.name in names
                or (
                    path.name == "trajectory.lammpstrj"
                    and relative.parts[0] in trajectories
                )
            ):
                data = path.read_bytes()
                hashes[str(relative)] = sha256(data).hexdigest()
                archive.writestr(str(relative), data)
        archive.writestr("manifest.json", json.dumps(hashes, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--lammps", default="lmp")
    parser.add_argument("--packmol", default="packmol")
    parser.add_argument(
        "--archive", type=Path, help="New ZIP for independent stage replay"
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    panel = cases()
    (args.out / "panel.json").write_text(json.dumps(panel, indent=2) + "\n")
    with ThreadPoolExecutor(args.workers) as pool:
        jobs = [
            pool.submit(execute, case, args.out, args.lammps, args.packmol)
            for case in panel
        ]
        records = []
        for job in as_completed(jobs):
            record = job.result()
            records.append(record)
            save_reports(args.out, records, len(panel))
            print(
                f"{len(records)}/{len(panel)} {record['name']}: MD={record['md']}",
                flush=True,
            )
    print(sum(r["md"] for r in records), "/", len(records), "completed MD")
    if args.archive:
        archive_runs(args.out, args.archive)


if __name__ == "__main__":
    main()
