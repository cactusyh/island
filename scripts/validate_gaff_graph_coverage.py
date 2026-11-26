"""Compare complete-graph assignment to native RadonPy, including junctions."""

import argparse
import json
from pathlib import Path

from rdkit import Chem

from island.forcefields.assignment import assign_polymer_parameters
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_core.schema import Branch
from island.forcefields.assignment_core.topology import molecule

CASES = [
    ("PE", "*CC*", 4),
    ("PP", "*C[C@H](C)*", 4),
    ("PEG", "*CCO*", 4),
    ("PET", "*OCCOC(=O)c1ccc(C(=O)*)cc1", 2),
    ("PA6", "*NCCCCCC(=O)*", 3),
    ("cis-diene", r"*C/C=C\C*", 4),
    ("side-substituted-ether", "*CC(C)(OCC*)C", 3),
    ("fluoro-chain", "*CC(F)C*", 4),
    ("chloro-chain", "*CC(Cl)C*", 4),
    ("aromatic-ether", "*CCOc1ccccc1*", 3),
]


def specifications():
    for name, smi, dp in CASES:
        yield name, PolymerSpec((("A", smi),), dp)
    yield (
        "PE-PEG-block",
        PolymerSpec(
            (("A", "*CC*"), ("B", "*CCO*")), 6, sequence=("A", "A", "B", "B", "A", "A")
        ),
    )
    yield (
        "PE-ether-branch",
        PolymerSpec((("A", "*CC*"),), 4, branches=(Branch(4, "*CCOCC"),)),
    )


def check(family):
    from radonpy.core import calc
    from radonpy.core.utils import Cell
    from radonpy.ff.gaff import GAFF
    from radonpy.ff.gaff2 import GAFF2
    from radonpy.sim.lammps import MolToLAMMPSdataBlock

    records = []
    for name, spec in specifications():
        row = {"name": name, "passed": False, "native_match": False}
        try:
            mol = molecule(normalize(spec))
            Chem.SanitizeMol(mol)
            assignment = assign_polymer_parameters(mol, family, charge="gasteiger")
            # Native engine sees the same explicit graph and coordinates, but no
            # assigned properties or topology objects from the ISLAND adapter.
            native = Chem.Mol(mol)
            native.RemoveAllConformers()
            native.AddConformer(Chem.Conformer(assignment.mol.GetConformer()))
            engine = {"gaff": GAFF, "gaff2": GAFF2}[family]()
            assert engine.ff_assign(native, useMDL=True, retryMDL=True)
            assert calc.assign_charges(native, charge="gasteiger")
            for m in (native, assignment.mol):
                m.cell = Cell(100, 0, 100, 0, 100, 0)
            # Exact native export compares every atom charge/type and every
            # bonded/nonbonded coefficient and connectivity, including junctions.
            expected = MolToLAMMPSdataBlock(native, velocity=False)
            actual = MolToLAMMPSdataBlock(assignment.mol, velocity=False)
            assert actual == expected, "Native topology/coefficients mismatch"
            row.update(
                passed=True,
                native_match=True,
                atoms=mol.GetNumAtoms(),
                database_sha256=assignment.database_sha256,
                lj14_scale=assignment.lj14_scale,
                coulomb14_scale=assignment.coulomb14_scale,
            )
        except Exception as exc:  # noqa: BLE001 - complete denominator
            row["error"] = str(exc)
        records.append(row)
    return {
        "family": family,
        "charge": "gasteiger",
        "attempted": len(records),
        "assigned": sum(r["passed"] for r in records),
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--families", nargs="+", default=["gaff", "gaff2"])
    args = parser.parse_args()
    result = {family: check(family) for family in args.families}
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {k: {x: v[x] for x in ("attempted", "assigned")} for k, v in result.items()}
        )
    )


if __name__ == "__main__":
    main()
