"""Numerical invariance and periodic MD for the v2 profile extension families."""

import argparse
import json
import os
from pathlib import Path

from island.epoxy.chemistry import monomer
from island.workflows.recipe import run_recipe
from scripts.validate_j107 import numerical_check

# These longer/combined contexts were not used while creating the typing rules.
MODELS = [
    ("chloro", "CCC(Cl)CC"),
    ("bromo", "CCC(Br)CC"),
    ("amine", "CCCN(CCC)CCC"),
    ("aryl-ether", "CCC(C)Oc1ccccc1"),
]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("field", type=Path)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    rows = []

    def save():
        (args.out / "checks.json").write_text(json.dumps(rows, indent=2) + "\n")

    for name, smiles in MODELS:
        row = {"kind": "heldout_numerical", "name": name, "smiles": smiles}
        try:
            row.update(
                numerical_check(
                    monomer(smiles, 173),
                    "opls",
                    args.field,
                    args.out / name,
                    os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                )
            )
            row["passed"] = True
        except Exception as exc:  # noqa: BLE001
            row.update(passed=False, error=str(exc))
        rows.append(row)
        save()
    for name, smiles in [
        ("PVC", "*CC(Cl)*"),
        ("PVBr", "*CC(Br)*"),
        ("polyamine", "*CCN(C)CC*"),
        ("aryl-ether", "*CC(Oc1ccccc1)*"),
    ]:
        cfg = {
            "architecture": "homopolymer",
            "polymer": {"repeat_units": {"A": smiles}, "dp": 4},
            "forcefield": {"name": "opls", "profile": "opls-polymer-v2"},
            "chains": 2,
            "density_g_cm3": 0.3,
            "seed": 173,
            "md": {
                "nvt_steps": 200,
                "npt_steps": 200,
                "production_steps": 200,
                "sample_every": 20,
                "segment_steps": 200,
            },
        }
        row = {"kind": "periodic_md", "name": name, "recipe": cfg}
        try:
            result = run_recipe(
                cfg,
                args.out / ("md-" + name),
                lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                packmol=os.environ["ISLAND_PACKMOL"],
            )
            row.update(passed=result["total_md_time_fs"] == 300, result=result)
        except Exception as exc:  # noqa: BLE001
            row.update(passed=False, error=str(exc))
        rows.append(row)
        save()
    if not all(r["passed"] for r in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
