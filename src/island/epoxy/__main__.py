"""Command line for compiling, verifying, curing and measuring epoxy networks."""

import argparse
import json
from pathlib import Path

from .cure import CureConfig, run_cure
from .parameters import compile_library
from .properties import conditions, cooling, tension
from .reacter import continuous
from .study import characterize, study
from .validation import opls_coverage, validate_library, validate_periodic_stage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in (
        "compile",
        "verify",
        "verify-periodic",
        "coverage",
        "cure",
        "study",
        "reacter",
        "condition",
        "characterize",
        "cool",
        "pull",
    ):
        p = sub.add_parser(name)
        p.add_argument("out", type=Path)
        if name in ("verify", "verify-periodic", "cure", "study", "reacter"):
            p.add_argument(
                "--library",
                type=Path,
                required=True,
                help="Directory produced by compile",
            )
        if name in ("cure", "study"):
            p.add_argument("--packmol", default="packmol")
        if name not in ("compile", "coverage"):
            p.add_argument("--lammps", default="lmp")
        if name == "compile":
            p.add_argument(
                "--field", choices=("pcff", "opls/2024/opls-aa"), default="pcff"
            )
            p.add_argument("--field-dir", type=Path)
            p.add_argument("--chemistry", choices=("pacm", "petmp"), default="pacm")
        if name == "cure":
            p.add_argument(
                "--config", type=Path, help="JSON object of CureConfig settings"
            )
            p.add_argument("--stop-after", type=int)
        if name == "verify-periodic":
            p.add_argument("--cure", required=True, type=Path)
        if name == "reacter":
            p.add_argument("--cure", type=Path, required=True)
            p.add_argument("--steps", type=int, default=60000)
        if name in ("condition", "pull", "cool"):
            p.add_argument("--source", required=True, type=Path)
        if name == "characterize":
            p.add_argument("--study", type=Path, required=True)
            p.add_argument("--extended", action="store_true")
        if name == "condition":
            p.add_argument("--temperature", type=float, default=300)
            p.add_argument("--steps", type=int, default=12000)
        if name in ("cool", "pull"):
            p.add_argument("--rate", type=float, required=True)
        if name == "pull":
            p.add_argument("--axis", choices=("x", "y", "z"), required=True)
    a = parser.parse_args()
    if a.command == "compile":
        compile_library(
            a.out, field=a.field, field_dir=a.field_dir, chemistry=a.chemistry
        )
        return
    if a.command == "coverage":
        opls_coverage(a.out)
        return
    if a.command == "verify-periodic":
        validate_periodic_stage(a.cure, a.library, a.out, a.lammps)
        return
    if a.command == "verify":
        validate_library(a.library, a.out, a.lammps)
        return
    if a.command == "cure":
        cfg = (
            CureConfig(**json.loads(a.config.read_text())) if a.config else CureConfig()
        )
        rec = run_cure(
            a.out,
            a.library / "library.json",
            cfg,
            lammps=a.lammps,
            packmol=a.packmol,
            stop_after=a.stop_after,
        )
        if rec:
            print(json.dumps(rec[-1]["statistics"]))
        return
    if a.command == "study":
        result = study(
            a.out, a.library / "library.json", lammps=a.lammps, packmol=a.packmol
        )
    elif a.command == "characterize":
        result = characterize(a.study, a.out, a.lammps, a.extended)
    elif a.command == "reacter":
        continuous(
            a.cure, a.library / "library.json", a.out, lammps=a.lammps, steps=a.steps
        )
        return
    elif a.command == "condition":
        result = conditions(a.source, a.out, a.lammps, a.temperature, a.steps)
    elif a.command == "cool":
        result = cooling(a.source, a.out, a.lammps, a.rate)
    else:
        result = tension(a.source, a.out, a.lammps, a.rate, a.axis)
    (a.out / "result.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
