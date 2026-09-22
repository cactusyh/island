import pytest

pytest.importorskip("rdkit")
from rdkit import Chem
from rdkit.Chem import AllChem

from island import AtomSite, BeadSite, Coordinates, MolecularSystem, Topology
from island.chemistry import from_rdkit, from_smiles, to_rdkit
from island.exceptions import MissingConformerError, UnsupportedRepresentationError


def embedded(smiles: str, *, hydrogens: bool = False) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if hydrogens:
        mol = Chem.AddHs(mol)
    assert AllChem.EmbedMolecule(mol, randomSeed=2026) == 0
    return mol


@pytest.mark.parametrize(
    "smiles", ["C", "CC", "CCO", "C=C", "C#C", "c1ccccc1", "CC(=O)C"]
)
def test_canonical_smiles_round_trip(smiles: str) -> None:
    source = embedded(smiles)
    system = from_rdkit(source).system
    restored = to_rdkit(system).mol
    assert Chem.MolToSmiles(restored) == Chem.MolToSmiles(source)


def test_smiles_helper_adds_hydrogens_coordinates_and_derived_topology() -> None:
    system = from_smiles("CCO", random_seed=17)
    assert system.representation == "atomistic"
    assert system.number_of_sites == 9
    assert len(system.coordinates) == 9
    assert system.topology.angles
    assert system.topology.dihedrals
    oxygen = [site for site in system.topology.sites.values() if site.element == "O"]
    assert oxygen[0].atomic_number == 8


@pytest.mark.parametrize(("smiles", "order"), [("CC", 1.0), ("C=C", 2.0), ("C#C", 3.0)])
def test_bond_orders_are_preserved(smiles: str, order: float) -> None:
    result = from_rdkit(embedded(smiles))
    assert next(iter(result.system.topology.bonds.values())).order == order
    assert to_rdkit(result.system).mol.GetBondWithIdx(0).GetBondTypeAsDouble() == order


def test_aromaticity_is_preserved() -> None:
    system = from_rdkit(embedded("c1ccccc1")).system
    assert all(
        bond.aromatic and bond.order == 1.5 for bond in system.topology.bonds.values()
    )
    assert all(bond.GetIsAromatic() for bond in to_rdkit(system).mol.GetBonds())


def test_formal_charge_is_preserved() -> None:
    system = from_rdkit(embedded("[NH4+]")).system
    site = next(iter(system.topology.sites.values()))
    assert site.charge == 1.0
    assert to_rdkit(system).mol.GetAtomWithIdx(0).GetFormalCharge() == 1


def test_noncontiguous_site_ids_are_explicitly_mapped() -> None:
    result = from_rdkit(embedded("CCO"), site_ids=[10, 20, 35])
    assert set(result.system.topology.sites) == {10, 20, 35}
    assert result.rdkit_index_to_site_id == {0: 10, 1: 20, 2: 35}
    assert to_rdkit(result.system).site_id_to_rdkit_index == {10: 0, 20: 1, 35: 2}


def test_missing_conformer_fails_clearly() -> None:
    with pytest.raises(MissingConformerError):
        from_rdkit(Chem.MolFromSmiles("CC"))


def test_bead_system_cannot_be_converted_to_rdkit() -> None:
    topology = Topology()
    topology.add_site(BeadSite(id=10, name="B", mass=44.0, bead_type="B"))
    system = MolecularSystem(
        topology, Coordinates({10: [0, 0, 0]}), representation="coarse_grained"
    )
    with pytest.raises(UnsupportedRepresentationError):
        to_rdkit(system)


def test_manual_atom_system_preserves_coordinates() -> None:
    topology = Topology()
    topology.add_site(
        AtomSite(id=35, name="Cl", mass=35.45, element="Cl", atomic_number=17)
    )
    result = to_rdkit(MolecularSystem(topology, Coordinates({35: [1, 2, 3]})))
    point = result.mol.GetConformer().GetAtomPosition(0)
    assert tuple(point) == pytest.approx((1, 2, 3))
