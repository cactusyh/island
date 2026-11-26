"""Seed, size and time contrasts for ordinary, polar and reacting recipes.

Finite smoke budgets intentionally retain the default 512-frame diagnostic
floor. A successful run may report budget_exhausted; that is not equilibrium.
"""

import argparse
import json
import os
from pathlib import Path

from island.workflows.recipe import run_recipe


def cases():
    for chemistry in ("PE", "PA6", "epoxy"):
        for contrast, seed, size, steps in [
            ("base", 173, 1, 200),
            ("seed", 719, 1, 200),
            ("size", 173, 2, 200),
            ("time", 173, 1, 400),
        ]:
            cfg = {
                "seed": seed,
                "density_g_cm3": 0.3,
                "md": {
                    "nvt_steps": 100,
                    "npt_steps": 100,
                    "production_steps": 100,
                    "sample_every": 10,
                    "segment_steps": 100,
                },
                "equilibration": {"max_segments": 2, "steps_per_segment": steps},
            }
            if chemistry == "epoxy":
                cfg.update(
                    architecture="crosslinked",
                    forcefield={"name": "gaff2", "charge": "gasteiger"},
                    cure={
                        "dgeba": 4 * size,
                        "pacm": 2 * size,
                        "pre_steps": 100,
                        "batch_steps": 100,
                        "post_steps": 100,
                        "max_idle": 2,
                        "cutoff": 5.5,
                        "target": 0.5,
                    },
                )
            else:
                cfg.update(
                    architecture="homopolymer",
                    chains=2 * size,
                    polymer={
                        "repeat_units": {
                            "A": "*CC*" if chemistry == "PE" else "*NCCCCCC(=O)*"
                        },
                        "dp": 4,
                    },
                    forcefield={"name": "opls", "profile": "opls-polymer-v2"}
                    if chemistry == "PE"
                    else {"name": "pcff", "profile": "pcff-carboxamide-v1"},
                )
            yield chemistry + "-" + contrast, cfg


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    rows = []
    for name, cfg in cases():
        row = {"name": name, "recipe": cfg}
        try:
            result = run_recipe(
                cfg,
                args.out / name,
                lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                packmol=os.environ["ISLAND_PACKMOL"],
            )
            assert (
                result["status"] == "budget_exhausted"
                and result["equilibrated"] is False
            )
            assert result["equilibration"]["segments"] == 2
            if name.startswith("epoxy"):
                assert result["cure"]["crosslinks"] > 0
            row.update(passed=True, result=result)
        except Exception as exc:  # noqa: BLE001
            row.update(passed=False, error_type=type(exc).__name__, diagnostic=str(exc))
        rows.append(row)
        (args.out / "checks.json").write_text(json.dumps(rows, indent=2) + "\n")
        print(name, row["passed"], flush=True)
    if not all(r["passed"] for r in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
