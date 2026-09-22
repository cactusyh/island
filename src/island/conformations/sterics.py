"""Geometric steric policies with no force-field parameter dependency."""

import math
from abc import ABC, abstractmethod
from collections.abc import Iterable, Mapping

import numpy as np
from numpy.typing import NDArray

from island.core import MolecularSystem

Pair = tuple[int, int]


class StericPolicy(ABC):
    """Decide minimum geometrically allowed distances between site pairs."""

    @abstractmethod
    def minimum_distance(
        self, system: MolecularSystem, site1: int, site2: int
    ) -> float:
        """Return the allowed minimum distance for a non-excluded pair."""


class FixedDistanceStericPolicy(StericPolicy):
    """Require the same minimum distance for every non-excluded pair."""

    def __init__(self, min_distance: float = 1.0) -> None:
        if not math.isfinite(min_distance) or min_distance <= 0:
            raise ValueError("min_distance must be positive")
        self.min_distance = float(min_distance)

    def minimum_distance(
        self, system: MolecularSystem, site1: int, site2: int
    ) -> float:
        del system, site1, site2
        return self.min_distance


class VanDerWaalsStericPolicy(StericPolicy):
    """Use scaled tabulated elemental van der Waals radii."""

    def __init__(self, scale: float = 0.6) -> None:
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError("scale must be positive")
        self.scale = float(scale)

    def minimum_distance(
        self, system: MolecularSystem, site1: int, site2: int
    ) -> float:
        first = system.topology.get_site(site1)
        second = system.topology.get_site(site2)
        return self.scale * (
            _vdw_radius(first.atomic_number) + _vdw_radius(second.atomic_number)
        )


def build_excluded_pairs(
    system: MolecularSystem, *, exclude_one_four: bool = False
) -> set[Pair]:
    """Build 1-2 and 1-3 exclusions; optionally exclude 1-4 pairs."""
    # Derive exclusions from bonds; callers need not have rebuilt cached angles.
    adjacency = {site_id: set() for site_id in system.topology.sites}
    for bond in system.topology.bonds.values():
        adjacency[bond.site1].add(bond.site2)
        adjacency[bond.site2].add(bond.site1)
    excluded: set[Pair] = set()
    for start in adjacency:
        seen = {start}
        frontier = {start}
        for _ in range(3 if exclude_one_four else 2):
            frontier = {n for site in frontier for n in adjacency[site]} - seen
            seen.update(frontier)
            excluded.update(tuple(sorted((start, end))) for end in frontier)
    return excluded


def find_clashes(
    system: MolecularSystem,
    positions: Mapping[int, NDArray[np.float64]],
    new_site_ids: Iterable[int],
    accepted_site_ids: Iterable[int],
    policy: StericPolicy,
    excluded_pairs: set[Pair],
) -> list[tuple[int, int, float, float]]:
    """Find new-versus-accepted geometric clashes incrementally."""
    clashes: list[tuple[int, int, float, float]] = []
    for site1 in new_site_ids:
        for site2 in accepted_site_ids:
            if tuple(sorted((site1, site2))) in excluded_pairs:
                continue
            distance = float(np.linalg.norm(positions[site1] - positions[site2]))
            threshold = policy.minimum_distance(system, site1, site2)
            if distance < threshold:
                clashes.append((site1, site2, distance, threshold))
    return clashes


def minimum_nonbonded_distance(
    positions: Mapping[int, NDArray[np.float64]], excluded_pairs: set[Pair]
) -> float | None:
    """Return the minimum distance over all non-excluded pairs."""
    site_ids = list(positions)
    minimum: float | None = None
    for index, site1 in enumerate(site_ids):
        for site2 in site_ids[index + 1 :]:
            if tuple(sorted((site1, site2))) in excluded_pairs:
                continue
            distance = float(np.linalg.norm(positions[site1] - positions[site2]))
            minimum = distance if minimum is None else min(minimum, distance)
    return minimum


def _vdw_radius(atomic_number: int | None) -> float:
    radii = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 16: 1.80, 17: 1.75}
    return radii.get(atomic_number, 1.70)
