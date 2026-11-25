"""Check recovered large-ring graph assignments against independent atom permutations."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from rdkit import Chem

from island.epoxy.validation import evaluate
from island.forcefields._emc_graph_input import graph_input
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.emc_graph import _section, assign_emc_graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--pcff-field", type=Path, required=True)
    parser.add_argument("--opls-field", type=Path, required=True)
    args = parser.parse_args()
    rows = json.loads(
        Path("docs/evidence/phase_4j107/repeat-only-cases.json").read_text()
    )
    root = args.out
    root.mkdir(parents=True, exist_ok=False)
    results = []
    for name in ["PI736", "PI936", "PI837", "PI1047"]:
        r = next(r for r in rows if r["id"] == name)
        r["field_dir"] = str(
            args.pcff_field if r["field"] == "pcff" else args.opls_field
        )
        g = normalize(PolymerSpec((("A", r["psmiles"]),), 4))
        source = graph_input(g, 2026)
        mol = source.mol
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, xyz in enumerate(source.xyz):
            conf.SetAtomPosition(i, xyz)
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        order = np.random.default_rng(41).permutation(mol.GetNumAtoms()).tolist()
        frames = []
        typed = []
        tries = []
        for label, mapping, m in [
            ("original", list(range(mol.GetNumAtoms())), mol),
            ("permuted", order, Chem.RenumberAtoms(mol, order)),
        ]:
            a = assign_emc_graph(
                m, r["field"], field_dir=r["field_dir"], backend_timeout=30
            )
            _, _, forces, energies = evaluate(
                root / name / label,
                a.data,
                a.parameters,
                os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
                direct_ewald=True,
            )
            indices = [mapping[a.site_groups[i][0]] for i in sorted(a.site_groups)]
            frames.append((forces[np.argsort(indices)], energies))
            typed.append(
                {
                    mapping[a.site_groups[int(row[0])][0]]: (row[-1], float(row[3]))
                    for line in _section(a.data, "Atoms").splitlines()
                    if (row := line.split())
                }
            )
            tries.append(a.assignment_attempts)
        df, de = [float(np.max(np.abs(a - b))) for a, b in zip(*frames)]
        result = {
            "id": name,
            "force_error": df,
            "energy_error": de,
            "types_charges_identical": typed[0] == typed[1],
            "passed": df < 1e-5 and de < 1e-6 and typed[0] == typed[1],
            "attempts": tries,
            "field_sha256": a.field_sha256,
        }
        results.append(result)
        print(name, result["passed"], df, de, flush=True)
        (root / "results.json").write_text(json.dumps(results, indent=2))

    assert len(results) == 4 and all(r["passed"] for r in results), results


if __name__ == "__main__":
    main()
