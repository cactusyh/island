"""Check native amide increments and permutation invariance on independent molecules."""

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from island.epoxy.chemistry import monomer
from island.epoxy.validation import evaluate
from island.forcefields.assignment import assign_emc_graph
from island.forcefields.emc_graph import _section


def check(fields, out, lammps):
    checks = []
    for name, smiles, expected in [
        ("primary", "CCC(=O)N", -0.8790),
        ("secondary", "CCC(=O)NCC", -0.6503),
        ("tertiary", "CCC(=O)N(C)CC", -0.4216),
        ("diamide", "CC(=O)NCCCNC(C)=O", -0.6503),
    ]:
        mol = monomer(smiles, 41)
        perm = np.random.default_rng(173).permutation(mol.GetNumAtoms()).tolist()
        frames = []
        for tag, order in [
            ("native", list(range(mol.GetNumAtoms()))),
            ("permuted", perm),
        ]:
            ref = assign_emc_graph(
                Chem.RenumberAtoms(mol, order), "pcff", field_dir=fields
            )
            atoms = [line.split() for line in _section(ref.data, "Atoms").splitlines()]
            n = [float(r[3]) for r in atoms if r[-1] == "n"]
            assert len(n) == sum(a.GetAtomicNum() == 7 for a in mol.GetAtoms())
            np.testing.assert_allclose(n, expected, atol=1e-8, rtol=0)
            _, _, f, e = evaluate(
                out / name / tag, ref.data, ref.parameters, lammps, direct_ewald=True
            )
            ids = [order[ref.site_groups[int(r[0])][0]] for r in atoms]
            frames.append((f[np.argsort(ids)], e))
        df, de = (float(np.max(np.abs(a - b))) for a, b in zip(*frames))
        if df > 1e-5 or de > 1e-6:
            raise ValueError(f"Permutation failure for {name}: {df}, {de}")
        checks.append(
            {
                "name": name,
                "smiles": smiles,
                "n_charge": expected,
                "max_force_error": df,
                "max_energy_component_error": de,
            }
        )
    return checks


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("field_dir", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--lammps", default="lmp")
    a = p.parse_args()
    print(json.dumps(check(a.field_dir, a.out, a.lammps), indent=2))
