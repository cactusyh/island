import pytest

from island import Angle, AtomSite, Dihedral, Improper, Topology
from island.exceptions import UnsupportedParameterRequirementError
from island.forcefields.parameters import derive_interaction_inventory


def topology_with_edges(
    site_ids: tuple[int, ...], edges: tuple[tuple[int, int], ...]
) -> Topology:
    topology = Topology()
    for site_id in site_ids:
        topology.add_site(
            AtomSite(site_id, f"C{site_id}", 12.0, element="C", atomic_number=6)
        )
    for site1, site2 in edges:
        topology.add_bond(site1, site2, order=1)
    return topology


def test_linear_inventory_has_independently_expected_counts_and_keys() -> None:
    topology = topology_with_edges((40, 10, 30, 20), ((40, 30), (30, 20), (20, 10)))
    inventory = derive_interaction_inventory(topology)
    assert inventory.bonds == ((10, 20), (20, 30), (30, 40))
    assert inventory.angles == ((10, 20, 30), (20, 30, 40))
    assert inventory.proper_torsions == ((10, 20, 30, 40),)


def test_inventory_ignores_empty_or_stale_angle_and_dihedral_caches() -> None:
    topology = topology_with_edges((1, 2, 3, 4), ((1, 2), (2, 3), (3, 4)))
    topology.angles.clear()
    topology.dihedrals.clear()
    empty_cache = derive_interaction_inventory(topology)
    topology.angles[(1, 2, 999)] = Angle(1, 2, 999)
    topology.dihedrals[(1, 2, 999, 3)] = Dihedral(1, 2, 999, 3)
    stale_cache = derive_interaction_inventory(topology)
    assert empty_cache == stale_cache
    assert empty_cache.angles == ((1, 2, 3), (2, 3, 4))


def test_branched_and_small_ring_inventories_are_simple_and_unique() -> None:
    branched = topology_with_edges((1, 2, 3, 4), ((1, 2), (1, 3), (1, 4)))
    branch_inventory = derive_interaction_inventory(branched)
    assert len(branch_inventory.bonds) == 3
    assert len(branch_inventory.angles) == 3
    assert branch_inventory.proper_torsions == ()

    triangle = topology_with_edges((1, 2, 3), ((1, 2), (2, 3), (3, 1)))
    triangle_inventory = derive_interaction_inventory(triangle)
    assert len(triangle_inventory.bonds) == 3
    assert len(triangle_inventory.angles) == 3
    assert triangle_inventory.proper_torsions == ()

    square = topology_with_edges((1, 2, 3, 4), ((1, 2), (2, 3), (3, 4), (4, 1)))
    square_inventory = derive_interaction_inventory(square)
    assert len(square_inventory.bonds) == 4
    assert len(square_inventory.angles) == 4
    assert len(square_inventory.proper_torsions) == 4
    assert all(len(set(item)) == 4 for item in square_inventory.proper_torsions)


def test_required_impropers_are_rejected_explicitly() -> None:
    topology = topology_with_edges((1, 2, 3, 4), ((1, 2), (1, 3), (1, 4)))
    topology.impropers.append(Improper(1, 2, 3, 4))
    with pytest.raises(UnsupportedParameterRequirementError, match="improper"):
        derive_interaction_inventory(topology)
