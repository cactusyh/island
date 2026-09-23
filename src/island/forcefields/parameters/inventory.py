"""Interaction inventories derived directly from authoritative chemical bonds."""

from dataclasses import dataclass
from itertools import combinations

from island.core import Topology
from island.exceptions import UnsupportedParameterRequirementError


@dataclass(frozen=True)
class InteractionInventory:
    """Canonical stable-site interactions required by the supported scope."""

    bonds: tuple[tuple[int, int], ...]
    angles: tuple[tuple[int, int, int], ...]
    proper_torsions: tuple[tuple[int, int, int, int], ...]


def derive_interaction_inventory(topology: Topology) -> InteractionInventory:
    """Derive bonds, angles, and simple proper torsions without using caches.

    Angles and torsions are equivalent under complete reversal. Proper torsions
    are simple length-three graph paths and therefore contain four distinct sites.
    """
    topology.validate_bond_graph()
    if topology.impropers:
        raise UnsupportedParameterRequirementError(
            "Phase 4B does not support required improper interactions; "
            f"the topology contains {len(topology.impropers)} improper(s)"
        )

    adjacency = {site_id: set() for site_id in topology.sites}
    bonds: set[tuple[int, int]] = set()
    for bond in topology.bonds.values():
        key = tuple(sorted((bond.site1, bond.site2)))
        bonds.add(key)
        adjacency[bond.site1].add(bond.site2)
        adjacency[bond.site2].add(bond.site1)

    angles: set[tuple[int, int, int]] = set()
    for center in sorted(adjacency):
        for left, right in combinations(sorted(adjacency[center]), 2):
            forward = (left, center, right)
            angles.add(min(forward, tuple(reversed(forward))))

    torsions: set[tuple[int, int, int, int]] = set()
    for center1, center2 in sorted(bonds):
        for left in sorted(adjacency[center1] - {center2}):
            for right in sorted(adjacency[center2] - {center1}):
                forward = (left, center1, center2, right)
                if len(set(forward)) != 4:
                    continue
                torsions.add(min(forward, tuple(reversed(forward))))

    return InteractionInventory(
        bonds=tuple(sorted(bonds)),
        angles=tuple(sorted(angles)),
        proper_torsions=tuple(sorted(torsions)),
    )
