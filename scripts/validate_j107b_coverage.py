"""Validate source-bound OPLS profile extensions and development graph probes."""

import argparse
import json
from pathlib import Path

from rdkit import Chem

from island.forcefields.assignment import assign_polymer_parameters
from island.forcefields.emc_graph import _field_elements
from island.forcefields.emc_profiles import install
from scripts.check_opls_polymer_profile import source_atoms, typed_atoms

SOURCE_MODELS = [
    ("MEAN", "CN"),
    ("ETAN", "CCN"),
    ("DMAN", "CNC"),
    ("TMAN", "CN(C)C"),
    ("ETCL", "CCCl"),
    ("ETBR", "CCBr"),
    ("ETIO", "CCI"),
    ("ANIS", "COc1ccccc1"),
]

DEVELOPMENT_PROBES = [
    ("branched-chloroalkane", "CC(Cl)C"),
    ("branched-bromoalkane", "CC(Br)C"),
    ("tertiary-amine", "CCN(CC)CC"),
    ("branched-aryl-ether", "CC(C)Oc1ccccc1"),
    ("secondary-amide", "CC(=O)NCC"),
]


def check(field_dir, out):
    field_dir = Path(field_dir)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    records = []
    for residue, smiles in SOURCE_MODELS:
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        expected = source_atoms(
            (field_dir / "opls/2024/opls-aa.define").read_text(),
            residue,
            molecule,
            _field_elements(field_dir, "opls/2024/opls-aa"),
        )
        assignment = assign_polymer_parameters(
            molecule, "opls/2024/opls-aa", field_dir=field_dir
        )
        actual = typed_atoms(assignment)
        mismatches = {
            atom: [actual[atom], value]
            for atom, value in expected.items()
            if actual[atom][0] != value[0] or abs(actual[atom][1] - value[1]) > 1e-6
        }
        records.append(
            {
                "kind": "source_residue",
                "residue": residue,
                "smiles": smiles,
                "atoms_checked": len(expected),
                "mismatches": mismatches,
                "passed": not mismatches,
            }
        )
    for name, smiles in DEVELOPMENT_PROBES:
        molecule = Chem.MolFromSmiles(smiles)
        try:
            assignment = assign_polymer_parameters(
                molecule, "opls/2024/opls-aa", field_dir=field_dir
            )
        except Exception as exc:  # noqa: BLE001 - preserve every held-out result
            records.append(
                {
                    "kind": "development_probe",
                    "name": name,
                    "smiles": smiles,
                    "passed": False,
                    "error": str(exc),
                }
            )
        else:
            records.append(
                {
                    "kind": "development_probe",
                    "name": name,
                    "smiles": smiles,
                    "passed": True,
                    "atoms": len(assignment.site_groups),
                    "serialization": assignment.serialization,
                }
            )
    result = {
        "profile": "opls-polymer-v2",
        "field_dir": str(field_dir),
        "source_models": len(SOURCE_MODELS),
        "development_models": len(DEVELOPMENT_PROBES),
        "passed": all(row["passed"] for row in records),
        "records": records,
    }
    (out / "checks.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--field-dir", type=Path)
    parser.add_argument("--install", action="store_true")
    args = parser.parse_args()
    field_dir = args.field_dir
    if args.install:
        if field_dir is None:
            parser.error("--install requires --field-dir as the output directory")
        install("opls-polymer-v2", field_dir)
    if field_dir is None:
        parser.error("--field-dir is required")
    result = check(field_dir, args.out)
    print(
        json.dumps(
            {
                k: result[k]
                for k in ("profile", "passed", "source_models", "development_models")
            }
        )
    )


if __name__ == "__main__":
    main()
