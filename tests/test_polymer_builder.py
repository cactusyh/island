import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from island.builders import build_linear_polymer
from island.chemistry import to_rdkit
from island.exceptions import InvalidRepeatUnitError, PolymerBuildError


def molecular_formula(system: object) -> str:
    return rdMolDescriptors.CalcMolFormula(to_rdkit(system).mol)


@pytest.mark.parametrize(
    ("dp", "formula"), [(1, "C2H6"), (2, "C4H10"), (3, "C6H14"), (5, "C10H22")]
)
def test_polyethylene_formula(dp: int, formula: str) -> None:
    system = build_linear_polymer("[*]CC[*]", dp=dp, random_seed=2026)
    assert molecular_formula(system) == formula
    assert system.number_of_sites == 6 * dp + 2


@pytest.mark.parametrize("dp", [1, 2, 3])
def test_connectivity_and_inter_repeat_bond_count(dp: int) -> None:
    system = build_linear_polymer("[*:1]CC[*:2]", dp=dp)
    assert len(system.topology.connected_components()) == 1
    cross_repeat_bonds = [
        bond
        for bond in system.topology.bonds.values()
        if system.topology.get_site(bond.site1).metadata["repeat_unit_index"]
        != system.topology.get_site(bond.site2).metadata["repeat_unit_index"]
    ]
    assert len(cross_repeat_bonds) == dp - 1
    assert system.metadata["polymer"]["number_of_inter_repeat_unit_bonds"] == dp - 1
    assert all(site.atomic_number != 0 for site in system.topology.sites.values())
    assert system.topology.angles
    assert system.topology.dihedrals


@pytest.mark.parametrize("psmiles", ["[*]CO[*]", "[*]CCO[*]", "[*]c1ccccc1[*]"])
def test_additional_repeat_unit_chemistries_build(psmiles: str) -> None:
    system = build_linear_polymer(psmiles, dp=2, random_seed=11)
    mol = to_rdkit(system).mol
    assert len(Chem.GetMolFrags(mol)) == 1
    assert not any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
    Chem.SanitizeMol(mol)


def test_atom_and_chain_provenance_is_complete() -> None:
    dp = 3
    system = build_linear_polymer("[*]CC[*]", dp=dp, chain_id="chain-7")
    repeat_indices = set()
    for site in system.topology.sites.values():
        assert site.metadata["chain_id"] == "chain-7"
        repeat_index = site.metadata["repeat_unit_index"]
        assert 0 <= repeat_index < dp
        repeat_indices.add(repeat_index)
        if site.element == "H":
            assert site.metadata["generated_hydrogen"] is True
            neighbor_id = next(iter(system.topology.neighbors(site.id)))
            neighbor = system.topology.get_site(neighbor_id)
            assert neighbor.metadata["repeat_unit_index"] == repeat_index
        else:
            assert "source_repeat_atom_index" in site.metadata
    assert repeat_indices == set(range(dp))

    metadata = system.metadata["polymer"]
    assert metadata["architecture"] == "linear"
    assert metadata["polymer_type"] == "homopolymer"
    assert metadata["source_psmiles"] == "[*]CC[*]"
    assert metadata["degree_of_polymerization"] == dp
    assert metadata["number_of_repeat_units"] == dp
    assert metadata["chain_id"] == "chain-7"
    assert metadata["head_site_id"] in system.topology.sites
    assert metadata["tail_site_id"] in system.topology.sites
    assert (
        system.topology.get_site(metadata["head_site_id"]).metadata["repeat_unit_index"]
        == 0
    )
    assert (
        system.topology.get_site(metadata["tail_site_id"]).metadata["repeat_unit_index"]
        == dp - 1
    )


def test_labeled_asymmetric_orientation_is_preserved() -> None:
    system = build_linear_polymer("[*:2]CO[*:1]", dp=2, generate_3d=False)
    polymer = system.metadata["polymer"]
    assert system.topology.get_site(polymer["head_site_id"]).element == "O"
    assert system.topology.get_site(polymer["tail_site_id"]).element == "C"
    cross_repeat = [
        bond
        for bond in system.topology.bonds.values()
        if system.topology.get_site(bond.site1).metadata["repeat_unit_index"]
        != system.topology.get_site(bond.site2).metadata["repeat_unit_index"]
    ]
    assert len(cross_repeat) == 1
    elements = {
        system.topology.get_site(cross_repeat[0].site1).element,
        system.topology.get_site(cross_repeat[0].site2).element,
    }
    assert elements == {"C", "O"}


def test_generate_3d_false_produces_explicit_2d_coordinates() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=2, generate_3d=False)
    positions = np.array(
        [system.coordinates.get(site_id) for site_id in system.topology.sites]
    )
    assert np.allclose(positions[:, 2], 0.0)
    assert system.metadata["polymer"]["coordinates"] == "rdkit_2d"


@pytest.mark.parametrize("dp", [0, -1, 2.5, True])
def test_invalid_degree_of_polymerization_is_rejected(dp: object) -> None:
    with pytest.raises(PolymerBuildError, match="integer >= 1"):
        build_linear_polymer("[*]CC[*]", dp=dp)


def test_unsupported_attachment_bond_order_is_rejected() -> None:
    with pytest.raises(InvalidRepeatUnitError, match="single-bond"):
        build_linear_polymer("[*]=CC=[*]", dp=2)
