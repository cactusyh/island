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
