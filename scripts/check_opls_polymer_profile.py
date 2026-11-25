"""Validate source residue atom types/charges and held-out OPLS graph assignments."""

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from island.epoxy.chemistry import monomer
from island.epoxy.validation import evaluate
from island.forcefields.assignment import assign_emc_graph
from island.forcefields.emc_graph import _section

SOURCE_MODELS = [
    ("DMET", "COC"),
    ("EMET", "COCC"),
    ("MACE", "CC(=O)OC"),
    ("EACE", "CC(=O)OCC"),
    ("ACAM", "CC(=O)N"),
    ("NMAM", "CC(=O)NC"),
    ("DMAM", "CC(=O)N(C)C"),
]


def source_atoms(text, residue, molecule, elements=None):
    block = text.split("RESI\t" + residue + "\t")[1].split("\nRESI\t", 1)[0]
    graph = Chem.RWMol()
    names, expected = {}, {}
    for line in block.splitlines():
        if line.startswith("ATOM\t"):
            _, name, kind, q = line.split("!")[0].split()
            index = graph.AddAtom(
                Chem.Atom(elements[kind] if elements is not None else kind[0].upper())
            )
            graph.GetAtomWithIdx(index).SetNoImplicit(True)
            names[name] = index
            expected[index] = (kind, float(q))
    for line in block.splitlines():
        if line.startswith("CONECT\t"):
            _, first, *neighbors = line.split("!")[0].split()
            for neighbor in neighbors:
                i, j = names[first], names[neighbor.lstrip("=")]
                if graph.GetBondBetweenAtoms(i, j) is None:
                    graph.AddBond(
                        i,
                        j,
                        Chem.BondType.DOUBLE
                        if neighbor.startswith("=")
                        else Chem.BondType.SINGLE,
                    )
    graph = graph.GetMol()
    Chem.SanitizeMol(graph)
    mapping = molecule.GetSubstructMatch(graph)
    if len(mapping) != molecule.GetNumAtoms():
        raise ValueError(f"Source residue does not map to {residue}")
    return {mapping[i]: value for i, value in expected.items()}


def typed_atoms(assignment):
    return {
        assignment.site_groups[int(row[0])][0]: (row[-1], float(row[3]))
        for line in _section(assignment.data, "Atoms").splitlines()
        if (row := line.split())
    }


def check(field, out, lammps):
    text = (field / "opls/2024/opls-aa.define").read_text()
    results = []
    for residue, smiles in SOURCE_MODELS:
        mol = monomer(smiles, 41)
        expected = source_atoms(text, residue, mol)
        actual = typed_atoms(assign_emc_graph(mol, "opls", field_dir=field))
        for atom, (kind, charge) in expected.items():
            if actual[atom][0] != kind or abs(actual[atom][1] - charge) > 1e-6:
                raise ValueError(
                    f"{residue} source mismatch at {atom}: {actual[atom]} vs {expected[atom]}"
                )
        results.append(
            {
                "residue": residue,
                "smiles": smiles,
                "source_atoms_checked": len(expected),
                "source_match": True,
            }
        )
    heldouts = [
        ("secondary-amide", "CCC(=O)NCCC", "n3am1", -0.5),
        ("tertiary-amide", "CCC(=O)N(CC)CCC", "n3am0", -0.14),
        ("branched-ester", "CCC(=O)OC(C)CC", "o2s", -0.33),
        ("extended-ether", "COCCOCCOC", "o2e", -0.4),
    ]
    for name, smiles, kind, charge in heldouts:
        mol = monomer(smiles, 173)
        perturbation = np.random.default_rng(41).normal(0, 0.02, (mol.GetNumAtoms(), 3))
        coords = mol.GetConformer().GetPositions()
        for i, xyz in enumerate(coords + perturbation):
            mol.GetConformer().SetAtomPosition(i, xyz)
        permutation = np.random.default_rng(719).permutation(mol.GetNumAtoms()).tolist()
        frames = []
        for label, order in [
            ("original", list(range(mol.GetNumAtoms()))),
            ("permuted", permutation),
        ]:
            assignment = assign_emc_graph(
                Chem.RenumberAtoms(mol, order), "opls", field_dir=field
            )
            actual = typed_atoms(assignment)
            target = [q for t, q in actual.values() if t == kind]
            if not target or max(abs(q - charge) for q in target) > 1e-6:
                raise ValueError(f"{name} source-derived charge mismatch")
            _, _, forces, energies = evaluate(
                out / name / label,
                assignment.data,
                assignment.parameters,
                lammps,
                direct_ewald=True,
            )
            mapping = [
                order[assignment.site_groups[i][0]]
                for i in sorted(assignment.site_groups)
            ]
            frames.append((forces[np.argsort(mapping)], energies))
        df, de = (float(np.max(np.abs(a - b))) for a, b in zip(*frames))
        if df > 1e-5 or de > 1e-6:
            raise ValueError(f"{name} permutation error")
        results.append(
            {
                "name": name,
                "smiles": smiles,
                "type_checked": kind,
                "expected_charge": charge,
                "max_force_error": df,
                "max_energy_component_error": de,
            }
        )
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field_dir", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--lammps", default="lmp")
    args = parser.parse_args()
    result = check(args.field_dir, args.out, args.lammps)
    (args.out / "checks.json").write_text(json.dumps(result, indent=2) + "\n")
