"""Numerical heldouts and real periodic MD for J111 source-bound profiles."""

import argparse
import json
import os
from pathlib import Path

from island.epoxy.chemistry import monomer
from island.workflows.recipe import run_recipe
from scripts.validate_j107 import numerical_check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--fields", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    records = []

    def save():
        (args.out / "checks.json").write_text(json.dumps(records, indent=2) + "\n")

    for family, name, smiles in [
        ("pcff", "formyl-aryl", "O=CNc1ccc(CC)cc1"),
        ("pcff", "aryl-carbamate", "CCCOC(=O)Nc1ccc(C)cc1"),
        ("opls", "tertiary-amide", "CCCC(=O)N(CCC)CCC"),
        ("opls", "alkyl-sulfide", "CCCSC(C)CC"),
        ("opls", "alkyl-sulfone", "CCCS(=O)(=O)CCC"),
        ("opls", "cyclic-ether", "CCCOC1CCCCCC1"),
    ]:
        row = {
            "kind": "numerical_heldout",
            "family": family,
            "name": name,
            "smiles": smiles,
        }
        try:
            row.update(
                numerical_check(
                    monomer(smiles, 173),
                    family,
                    args.fields / (family + "-release"),
                    args.out / name,
                    os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                    ewald_accuracy=1e-10,
                )
            )
            row["passed"] = True
        except Exception as exc:  # noqa: BLE001 - preserve failed evidence
            row.update(passed=False, error=str(exc))
        records.append(row)
        save()
    for family, name, repeat, ends in [
        ("pcff", "aromatic-PU", "*OCCCCOC(=O)Nc1ccc(NC(=O)*)cc1", ["*C", "*OC"]),
        ("pcff", "formyl-polyamide", "*C(=O)Nc1ccc(*)cc1", [None, None]),
        ("opls", "amide-sidechain", "*CC(c1ccc(C(=O)N(CC)CC)cc1)*", [None, None]),
        ("opls", "polysulfone", "*CCCCS(=O)(=O)*", ["*C", "*C"]),
        ("opls", "cyclic-ether-sidechain", "*CC(OC1CCCCC1)*", [None, None]),
    ]:
        recipe = {
            "architecture": "homopolymer",
            "polymer": {"repeat_units": {"A": repeat}, "dp": 3, "end_groups": ends},
            "forcefield": {
                "name": family,
                "profile": {"pcff": "pcff-polymer-v2", "opls": "opls-polymer-v3"}[
                    family
                ],
            },
            "chains": 2,
            "density_g_cm3": 0.3,
            "seed": 173,
            "md": {
                "nvt_steps": 100,
                "npt_steps": 100,
                "production_steps": 100,
                "sample_every": 10,
                "segment_steps": 100,
            },
        }
        row = {"kind": "periodic_md", "family": family, "name": name, "recipe": recipe}
        try:
            result = run_recipe(
                recipe,
                args.out / ("md-" + name),
                lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                packmol=os.environ["ISLAND_PACKMOL"],
            )
            row.update(passed=result["status"] == "completed", result=result)
        except Exception as exc:  # noqa: BLE001 - preserve failed evidence
            row.update(passed=False, error=str(exc))
        records.append(row)
        save()
    if not all(r["passed"] for r in records):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
