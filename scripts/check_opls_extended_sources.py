"""Independent native residue checks for sulfur and halogen model compounds."""

import argparse
import json
from pathlib import Path

from island.epoxy.chemistry import monomer
from island.exceptions import ForceFieldError
from island.forcefields.assignment_diagnostics import assignment_failure
from island.forcefields.emc_graph import _field_elements, assign_emc_graph
from scripts.check_opls_polymer_profile import source_atoms, typed_atoms

MODELS = [
    ("DMES", "CSC"),
    ("EMES", "CSCC"),
    ("DMSO", "CS(C)=O"),
    ("MSME", "CS(C)(=O)=O"),
    ("ETBR", "CCBr"),
    ("ETCL", "CCCl"),
]


def check(field):
    text = (field / "opls/2024/opls-aa.define").read_text()
    elements = _field_elements(field, "opls/2024/opls-aa")
    records = []
    for residue, smiles in MODELS:
        mol = monomer(smiles, 173)
        expected = source_atoms(text, residue, mol, elements)
        record = {
            "residue": residue,
            "smiles": smiles,
            "source_atoms": {str(i): list(v) for i, v in expected.items()},
        }
        try:
            result = assign_emc_graph(mol, "opls", field_dir=field)
        except ForceFieldError as exc:
            record.update(
                assigned=False,
                source_match="unmeasured",
                failure=assignment_failure(exc, mol),
            )
        else:
            actual = typed_atoms(result)
            mismatch = {
                i: (actual.get(i), v)
                for i, v in expected.items()
                if i not in actual
                or actual[i][0] != v[0]
                or abs(actual[i][1] - v[1]) > 1e-6
            }
            if mismatch or actual.keys() != expected.keys():
                raise ValueError(f"{residue} native residue mismatch: {mismatch}")
            record.update(
                assigned=True, source_match=True, field_sha256=result.field_sha256
            )
        records.append(record)
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(check(args.field_dir), indent=2) + "\n")
