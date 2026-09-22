import pytest

from island import AtomSite, Topology, TopologyError


def topology_with_sites(count: int = 4) -> Topology:
    topology = Topology()
    for site_id in range(10, 10 + count):
        topology.add_site(AtomSite(id=site_id, name="C", mass=12.0, element="C"))
    return topology


def test_bond_neighbors_and_removal() -> None:
    topology = topology_with_sites(2)
    topology.add_bond(10, 11, order=1)
    assert topology.neighbors(10) == {11}
    assert topology.degree(11) == 1
    assert topology.remove_bond(11, 10).order == 1
    assert topology.degree(10) == 0


def test_duplicate_and_unknown_bonds_are_rejected() -> None:
    topology = topology_with_sites(2)
    topology.add_bond(10, 11)
    with pytest.raises(TopologyError, match="Duplicate"):
        topology.add_bond(11, 10)
    with pytest.raises(TopologyError, match="Unknown"):
        topology.add_bond(10, 99)


def test_connected_components_use_site_ids_not_positions() -> None:
    topology = topology_with_sites(4)
    topology.add_bond(10, 11)
    topology.add_bond(12, 13)
    assert {frozenset(component) for component in topology.connected_components()} == {
        frozenset({10, 11}),
        frozenset({12, 13}),
    }


def test_rebuild_generates_unique_angles_and_proper_dihedrals() -> None:
    topology = topology_with_sites(4)
    topology.add_bond(10, 11)
    topology.add_bond(11, 12)
    topology.add_bond(12, 13)
    topology.rebuild_derived_interactions()
    assert len(topology.angles) == 2
    assert len(topology.dihedrals) == 1
