import runpy
from pathlib import Path

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("pyemc")
from rdkit import Chem

from island.forcefields.assignment import assign_emc_graph
from island.forcefields.emc_graph import _section
from island.forcefields.emc_profiles import install


def test_ether_hydrogens_and_amide_alkyl_atoms_match_source_residues(tmp_path):
    profile = install("opls-polymer-v1", tmp_path / "field")
    for name in ("opls/2024/opls-aa.define", "opls/2024/opls-aa.prm"):
        assert profile["installed"][name] == profile["base"][name]
    for smiles, indices, kind, charge in [
        ("COC", [0, 2], "c4oe", 0.11),
        ("CC(=O)NCC", [4], "c4", 0.08),
    ]:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assigned = assign_emc_graph(mol, "opls", field_dir=tmp_path / "field")
        actual = {
            assigned.site_groups[int(r[0])][0]: (r[-1], float(r[3]))
            for line in _section(assigned.data, "Atoms").splitlines()
            if (r := line.split())
        }
        for i in indices:
            assert actual[i] == (kind, charge)
            for neighbor in mol.GetAtomWithIdx(i).GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                    assert actual[neighbor.GetIdx()][0] == "h1"


@pytest.mark.parametrize(
    "residue,smiles",
    [
        ("DMET", "COC"),
        ("EMET", "COCC"),
        ("MACE", "CC(=O)OC"),
        ("EACE", "CC(=O)OCC"),
        ("ACAM", "CC(=O)N"),
        ("NMAM", "CC(=O)NC"),
        ("DMAM", "CC(=O)N(C)C"),
    ],
)
def test_every_source_residue_atom_retains_native_type_and_charge(
    tmp_path, residue, smiles
):
    # The independent .define residue graph supplies the expected values;
    # the changed .top rules are never used as the reference.
    checks = runpy.run_path(
        str(
            Path(__file__).resolve().parents[1]
            / "scripts/check_opls_polymer_profile.py"
        )
    )
    field = tmp_path / "field"
    install("opls-polymer-v1", field)
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    expected = checks["source_atoms"](
        (field / "opls/2024/opls-aa.define").read_text(), residue, mol
    )
    actual = checks["typed_atoms"](assign_emc_graph(mol, "opls", field_dir=field))
    assert set(actual) == set(expected) == set(range(mol.GetNumAtoms()))
    for atom, (kind, charge) in expected.items():
        assert actual[atom][0] == kind
        assert actual[atom][1] == pytest.approx(charge, abs=1e-6)
