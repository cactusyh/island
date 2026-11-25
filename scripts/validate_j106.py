"""Real-engine J106 entry consistency and short periodic-MD regression evidence."""

import argparse
import json
import os
import platform
from pathlib import Path

import rdkit
from rdkit import Chem

from island.forcefields.assignment import assign_polymer_parameters
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.emc_profiles import install
from island.workflows.recipe import run_recipe
from island.workflows.recipe_config import PolymerRecipe


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    records = []

    def save():
        (root / "results.json").write_text(
            json.dumps(
                {
                    "python": platform.python_version(),
                    "rdkit": rdkit.__version__,
                    "equilibration": "not tested; short integration regression only",
                    "results": records,
                },
                indent=2,
            )
            + "\n"
        )

    cases = {
        "PEG": "*CCO*",
        "PA6": "*NCCCCCC(=O)*",
        "PET": "*OCCOC(=O)c1ccc(C(=O)*)cc1",
        "PU": "*OCCOC(=O)NCCCCCCNC(=O)*",
    }
    for field, profile in [
        ("pcff", "pcff-carboxamide-v1"),
        ("opls", "opls-polymer-v1"),
    ]:
        directory = root / (field + "-field")
        install(profile, directory)
        for name, repeat in cases.items():
            row = {"test": "entry_equivalence", "name": name, "field": field, "dp": 3}
            try:
                a = assign_polymer_parameters(repeat, field, dp=3, field_dir=directory)
                b = assign_polymer_parameters(
                    normalize(PolymerSpec((("A", repeat),), 3)),
                    field,
                    field_dir=directory,
                )

                # Same site order, chemistry, charges and coefficients. Coordinates
                # can depend on the graph adapter's embedding path.
                def atoms(result):
                    return [
                        r.split()[:4] + r.split()[7:]
                        for r in result.data.split("\nAtoms\n\n")[1]
                        .split("\n\n")[0]
                        .splitlines()
                    ]

                assert a.input_smiles == b.input_smiles
                assert a.site_groups == b.site_groups
                assert a.parameters == b.parameters
                assert atoms(a) == atoms(b)
                row.update(
                    passed=True,
                    field_sha256=a.field_sha256,
                    input_smiles=a.input_smiles,
                )
            except Exception as exc:  # noqa: BLE001 - retain every failed case in evidence
                row.update(passed=False, error=type(exc).__name__, diagnostic=str(exc))
            records.append(row)
            save()
            print(name, field, row["passed"], flush=True)

    for field in ("gaff", "gaff2"):
        for repeat in ("*C[C@H](C)*", r"*C/C=C\C*"):
            row = {"test": "gaff_entry_equivalence", "repeat": repeat, "field": field}
            try:
                a = assign_polymer_parameters(
                    repeat, field, dp=3, charge="gasteiger", atom_order="radonpy"
                )
                graph = normalize(PolymerSpec((("A", repeat),), 3))
                b = assign_polymer_parameters(graph, field, charge="gasteiger")
                expected = Chem.MolToSmiles(Chem.RemoveHs(molecule(graph)))
                assert Chem.MolToSmiles(Chem.RemoveHs(a.mol)) == expected
                assert Chem.MolToSmiles(Chem.RemoveHs(b.mol)) == expected
                assert sorted(x.GetProp("ff_type") for x in a.mol.GetAtoms()) == sorted(
                    x.GetProp("ff_type") for x in b.mol.GetAtoms()
                )
                row.update(passed=True, smiles=expected)
            except Exception as exc:  # noqa: BLE001 - retain every failed case in evidence
                row.update(passed=False, error=type(exc).__name__, diagnostic=str(exc))
            records.append(row)
            save()

    for field, repeat in [
        ("gaff2", "*C[C@H](C)*"),
        ("pcff", r"*C/C=C\C*"),
        ("opls", r"*C/C=C\C*"),
    ]:
        ff = (
            {"name": field, "charge": "gasteiger"}
            if field == "gaff2"
            else {
                "name": field,
                "profile": "pcff-carboxamide-v1"
                if field == "pcff"
                else "opls-polymer-v1",
            }
        )
        config = {
            "architecture": "homopolymer",
            "polymer": {"repeat_units": {"A": repeat}, "dp": 3},
            "forcefield": ff,
            "chains": 2,
            "density_g_cm3": 0.3,
            "seed": 173,
            "md": {
                "nvt_steps": 40,
                "npt_steps": 40,
                "production_steps": 40,
                "sample_every": 10,
                "segment_steps": 40,
            },
        }
        row = {"test": "periodic_md", "field": field, "recipe": config}
        try:
            output = root / (field + "-md")
            run_recipe(
                PolymerRecipe.from_dict(config),
                output,
                lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                packmol=os.environ["ISLAND_PACKMOL"],
            )
            status = json.loads((output / "status.json").read_text())
            assert status["status"] == "completed", status
            result = json.loads((output / "result.json").read_text())
            assert result["total_md_time_fs"] == 60.0
            assert len(result["stages"]) == 3
            assert (output / "production/final.restart").stat().st_size > 0
            row["stages"] = result["stages"]
            row["total_md_time_fs"] = result["total_md_time_fs"]
            row.update(passed=True, status=status)
        except Exception as exc:  # noqa: BLE001 - retain every failed case in evidence
            row.update(passed=False, error=type(exc).__name__, diagnostic=str(exc))
        records.append(row)
        save()
        print(field, "MD", row["passed"], flush=True)
    if not all(row["passed"] for row in records):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
