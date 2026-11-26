"""Source-residue, junction, force-gradient and periodic-MD evidence for J107."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from rdkit import Chem

from island.epoxy.chemistry import monomer
from island.epoxy.validation import evaluate
from island.forcefields.assignment import assign_emc_graph
from island.forcefields.emc_graph import _section
from island.workflows.recipe import run_recipe
from island.workflows.recipe_config import PolymerRecipe
from scripts.check_opls_polymer_profile import check as opls_check
from scripts.check_pcff_amides import check as pcff_check


def displaced(data, atom_id, axis, delta):
    block = _section(data, "Atoms")
    rows = [r.split() for r in block.splitlines()]
    for row in rows:
        if int(row[0]) == atom_id:
            row[4 + axis] = format(float(row[4 + axis]) + delta, ".16g")
    return data.replace(block, "\n".join(" ".join(r) for r in rows), 1)


def numerical_check(
    mol, field, directory, out, lammps, *, ewald_accuracy=1e-7, coulomb_table=True
):
    # Distortion prevents an equilibrium geometry from masking force errors.
    xyz = mol.GetConformer().GetPositions() + np.random.default_rng(193).normal(
        0, 0.02, (mol.GetNumAtoms(), 3)
    )
    for i, pos in enumerate(xyz):
        mol.GetConformer().SetAtomPosition(i, pos)
    assignment = assign_emc_graph(mol, field, field_dir=directory)
    _, _, force, energy = evaluate(
        out / "reference",
        assignment.data,
        assignment.parameters,
        lammps,
        direct_ewald=True,
        ewald_accuracy=ewald_accuracy,
        coulomb_table=coulomb_table,
    )
    order = np.random.default_rng(719).permutation(mol.GetNumAtoms()).tolist()
    permuted = assign_emc_graph(
        Chem.RenumberAtoms(mol, order), field, field_dir=directory
    )
    _, _, other_force, other_energy = evaluate(
        out / "permuted",
        permuted.data,
        permuted.parameters,
        lammps,
        direct_ewald=True,
        ewald_accuracy=ewald_accuracy,
        coulomb_table=coulomb_table,
    )
    first_order = [assignment.site_groups[i][0] for i in sorted(assignment.site_groups)]
    second_order = [
        order[permuted.site_groups[i][0]] for i in sorted(permuted.site_groups)
    ]
    df = float(
        np.max(
            np.abs(
                force[np.argsort(first_order)] - other_force[np.argsort(second_order)]
            )
        )
    )
    de = float(np.max(np.abs(energy - other_energy)))
    assert df < 1e-5 and de < 1e-6, (df, de)
    atom_id = assignment.site_to_atom_id[
        next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))
    ]
    errors = []
    for axis in range(3):
        energies = []
        for sign in (-1, 1):
            _, _, _, value = evaluate(
                out / f"fd-{axis}-{sign}",
                displaced(assignment.data, atom_id, axis, sign * 1e-4),
                assignment.parameters,
                lammps,
                direct_ewald=True,
                ewald_accuracy=ewald_accuracy,
                coulomb_table=coulomb_table,
            )
            energies.append(float(np.sum(value)))
        finite_difference = -(energies[1] - energies[0]) / 2e-4
        errors.append(abs(finite_difference - force[atom_id - 1, axis]))
    assert max(errors) < 5e-3, errors
    return {
        "permutation_force_max_error": df,
        "permutation_energy_component_max_error": de,
        "force_gradient_max_error": max(errors),
        "gradient_atom_id": atom_id,
        "gradient_step_A": 1e-4,
        "ewald_accuracy": ewald_accuracy,
        "coulomb_table": coulomb_table,
        "field_sha256": assignment.field_sha256,
        "energy_kcal_mol": float(np.sum(energy)),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--pcff-field", type=Path, required=True)
    p.add_argument("--opls-field", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    lammps = os.environ["ISLAND_ASSIGNMENT_LAMMPS"]
    records = []

    def save():
        (args.out / "checks.json").write_text(json.dumps(records, indent=2) + "\n")

    for field, checker, directory in [
        ("opls", opls_check, args.opls_field),
        ("pcff", pcff_check, args.pcff_field),
    ]:
        values = checker(directory, args.out / (field + "-source"), lammps)
        records.append(
            {
                "kind": "source_and_junction_checks",
                "field": field,
                "passed": True,
                "checks": values,
            }
        )
        save()
        for name, smi in [
            ("ether-junction", "CCOCCOCCOC"),
            ("amide-junction", "CCC(=O)NCCC"),
            ("ester-branch", "CCC(=O)OC(C)CC"),
        ]:
            result = numerical_check(
                monomer(smi, 173), field, directory, args.out / field / name, lammps
            )
            records.append(
                {
                    "kind": "force_gradient_and_permutation",
                    "field": field,
                    "name": name,
                    "smiles": smi,
                    "passed": True,
                    **result,
                }
            )
            save()
    recipes = [
        (
            "pcff-PDMS-TMS",
            "pcff",
            {
                "repeat_units": {"A": "*O[Si](C)(C)*"},
                "dp": 4,
                "end_groups": ["*[Si](C)(C)C", "*O[Si](C)(C)C"],
            },
        ),
        (
            "opls-PEG-block",
            "opls",
            {
                "repeat_units": {"A": "*CCO*", "B": "*CC*"},
                "blocks": [["A", 2], ["B", 2], ["A", 2]],
            },
        ),
        (
            "opls-PA6-methyl",
            "opls",
            {
                "repeat_units": {"A": "*NCCCCCC(=O)*"},
                "dp": 3,
                "end_groups": ["*C", "*C"],
            },
        ),
    ]
    for name, field, polymer in recipes:
        cfg = {
            "architecture": "block" if "blocks" in polymer else "homopolymer",
            "polymer": polymer,
            "forcefield": {
                "name": field,
                "profile": "pcff-carboxamide-v1"
                if field == "pcff"
                else "opls-polymer-v1",
            },
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
                PolymerRecipe.from_dict(cfg),
                args.out / name,
                lammps=lammps,
                packmol=os.environ["ISLAND_PACKMOL"],
            )
            assert result["total_md_time_fs"] == 300
            row.update(passed=True, result=result)
        except Exception as exc:  # noqa: BLE001 - retain failed MD separately from assignment
            row.update(passed=False, error_type=type(exc).__name__, diagnostic=str(exc))
        records.append(row)
        save()
    if not all(r["passed"] for r in records):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
