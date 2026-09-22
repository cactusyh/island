"""Force-field-independent self-avoiding random-walk coordinates."""

import math
import random
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from island.conformations.base import ConformationGenerator
from island.conformations.result import ConformationResult
from island.conformations.sterics import (
    FixedDistanceStericPolicy,
    StericPolicy,
    build_excluded_pairs,
    find_clashes,
    minimum_nonbonded_distance,
)
from island.conformations.torsions import TorsionSampler, UniformTorsionSampler
from island.core import AtomSite, Coordinates, MolecularSystem
from island.exceptions import (
    ConformationGenerationError,
    UnsupportedConformationError,
)


@dataclass(frozen=True)
class RetryPolicy:
    """Retry and rollback limits for incremental repeat-unit placement."""

    attempts_per_unit: int = 100
    rollback_units: int = 5
    max_rollbacks: int = 20

    def __post_init__(self) -> None:
        for name in ("attempts_per_unit", "rollback_units"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError(f"{name} must be an integer >= 1")
        if (
            not isinstance(self.max_rollbacks, int)
            or isinstance(self.max_rollbacks, bool)
            or self.max_rollbacks < 0
        ):
            raise ValueError("max_rollbacks must be an integer >= 0")


@dataclass
class SelfAvoidingRandomWalkGenerator(ConformationGenerator):
    """Place rigid repeat-unit geometries incrementally with steric rejection."""

    seed: int = 2026
    torsion_sampler: TorsionSampler = field(default_factory=UniformTorsionSampler)
    steric_policy: StericPolicy = field(
        default_factory=lambda: FixedDistanceStericPolicy(0.8)
    )
    retry_policy: RetryPolicy = field(default_factory=RetryPolicy)
    bond_angle_degrees: float = 109.5
    direction_jitter_degrees: float = 25.0
    exclude_one_four: bool = False

    def generate(self, system: MolecularSystem) -> ConformationResult:
        """Generate new coordinates without modifying ``system``."""
        if not isinstance(self.seed, int) or isinstance(self.seed, bool):
            raise ValueError("seed must be an integer")
        layout = _validate_linear_polymer(system)
        source = {
            site_id: system.coordinates.get(site_id)
            for site_id in system.topology.sites
        }
        positions: dict[int, NDArray[np.float64]] = {}
        first_ids = layout.unit_site_ids[0]
        origin = source[layout.incoming_anchors[0]]
        for site_id in first_ids:
            positions[site_id] = source[site_id] - origin

        rng = random.Random(self.seed)
        excluded_pairs = build_excluded_pairs(
            system, exclude_one_four=self.exclude_one_four
        )
        attempts = 0
        rejected_trials = 0
        rollback_count = 0
        unit_attempts = [0] * layout.number_of_units
        accepted_torsions: dict[int, float] = {}
        connection_directions: dict[int, NDArray[np.float64]] = {}
        repeat_index = 1

        while repeat_index < layout.number_of_units:
            placed = False
            for trial in range(1, self.retry_policy.attempts_per_unit + 1):
                attempts += 1
                unit_attempts[repeat_index] += 1
                torsion = self.torsion_sampler.sample(
                    rng, repeat_index=repeat_index, trial=trial
                )
                proposed, direction = self._propose_unit(
                    system,
                    layout,
                    source,
                    positions,
                    repeat_index,
                    torsion,
                    rng,
                    connection_directions,
                )
                accepted_ids = [
                    site_id
                    for index in range(repeat_index)
                    for site_id in layout.unit_site_ids[index]
                ]
                clashes = find_clashes(
                    system,
                    {**positions, **proposed},
                    layout.unit_site_ids[repeat_index],
                    accepted_ids,
                    self.steric_policy,
                    excluded_pairs,
                )
                if clashes:
                    rejected_trials += 1
                    continue
                positions.update(proposed)
                accepted_torsions[repeat_index] = torsion
                connection_directions[repeat_index] = direction
                placed = True
                repeat_index += 1
                break

            if placed:
                continue
            if rollback_count >= self.retry_policy.max_rollbacks:
                raise ConformationGenerationError(
                    "Self-avoiding random walk exhausted retries at repeat "
                    f"{repeat_index} after {attempts} placement attempts and "
                    f"{rollback_count} rollbacks",
                    repeat_index=repeat_index,
                    attempts=attempts,
                    rejected_trials=rejected_trials,
                    rollback_count=rollback_count,
                )
            rollback_count += 1
            rollback_start = max(1, repeat_index - self.retry_policy.rollback_units)
            for index in range(rollback_start, repeat_index):
                for site_id in layout.unit_site_ids[index]:
                    positions.pop(site_id, None)
                accepted_torsions.pop(index, None)
                connection_directions.pop(index, None)
            repeat_index = rollback_start

        coordinates = Coordinates(positions)
        coordinates.validate(system.topology)
        minimum_distance = minimum_nonbonded_distance(positions, excluded_pairs)
        all_ids = list(system.topology.sites)
        position_array = np.array([positions[site_id] for site_id in all_ids])
        center = position_array.mean(axis=0)
        radius_of_gyration = float(
            np.sqrt(np.mean(np.sum((position_array - center) ** 2, axis=1)))
        )
        end_to_end_distance = float(
            np.linalg.norm(
                positions[layout.tail_site_id] - positions[layout.head_site_id]
            )
        )
        return ConformationResult(
            coordinates=coordinates,
            method="self_avoiding_random_walk",
            seed=self.seed,
            success=True,
            attempts=attempts,
            rejected_trials=rejected_trials,
            rollback_count=rollback_count,
            minimum_nonbonded_distance=minimum_distance,
            metadata={
                "attempts_per_repeat_unit": unit_attempts,
                "accepted_torsions_degrees": {
                    index: math.degrees(angle)
                    for index, angle in sorted(accepted_torsions.items())
                },
                "steric_policy": type(self.steric_policy).__name__,
                "excluded_pairs": "1-2 and 1-3",
                "one_four_pairs_checked": not self.exclude_one_four,
                "repeat_units_are_rigid": True,
                "ring_intersection_check": False,
                "end_to_end_distance": end_to_end_distance,
                "radius_of_gyration": radius_of_gyration,
            },
        )

    def _propose_unit(
        self,
        system: MolecularSystem,
        layout: "_PolymerLayout",
        source: dict[int, NDArray[np.float64]],
        positions: dict[int, NDArray[np.float64]],
        repeat_index: int,
        torsion: float,
        rng: random.Random,
        connection_directions: dict[int, NDArray[np.float64]],
    ) -> tuple[dict[int, NDArray[np.float64]], NDArray[np.float64]]:
        previous_anchor = layout.outgoing_anchors[repeat_index - 1]
        current_anchor = layout.incoming_anchors[repeat_index]
        previous_incoming = layout.incoming_anchors[repeat_index - 1]
        base = positions[previous_anchor] - positions[previous_incoming]
        if float(np.linalg.norm(base)) < 1e-10:
            base = connection_directions.get(repeat_index - 1, _random_unit_vector(rng))
        base = _normalize(base)
        bond_direction = _sample_cone_direction(
            base,
            rng,
            center_angle=math.pi - math.radians(self.bond_angle_degrees),
            half_width=math.radians(self.direction_jitter_degrees),
        )
        original_bond_length = float(
            np.linalg.norm(source[current_anchor] - source[previous_anchor])
        )
        if original_bond_length < 1e-6:
            raise ConformationGenerationError(
                f"Inter-repeat bond before repeat {repeat_index} has zero length",
                repeat_index=repeat_index,
            )
        target_anchor = (
            positions[previous_anchor] + bond_direction * original_bond_length
        )
        unit_ids = layout.unit_site_ids[repeat_index]
        local_anchor = source[current_anchor]
        reference = source[layout.outgoing_anchors[repeat_index]] - local_anchor
        if float(np.linalg.norm(reference)) < 1e-10:
            reference = max(
                (source[site_id] - local_anchor for site_id in unit_ids),
                key=lambda vector: float(np.linalg.norm(vector)),
            )
        bend_angle = math.pi - math.radians(self.bond_angle_degrees)
        perpendicular, second_perpendicular = _perpendicular_basis(bond_direction)
        desired_reference = math.cos(bend_angle) * bond_direction + math.sin(
            bend_angle
        ) * (
            math.cos(torsion) * perpendicular + math.sin(torsion) * second_perpendicular
        )
        rotation = _rotation_between(reference, desired_reference)
        proposed = {
            site_id: rotation @ (source[site_id] - local_anchor) + target_anchor
            for site_id in unit_ids
        }
        return proposed, bond_direction


@dataclass(frozen=True)
class _PolymerLayout:
    unit_site_ids: tuple[tuple[int, ...], ...]
    incoming_anchors: tuple[int, ...]
    outgoing_anchors: tuple[int, ...]
    head_site_id: int
    tail_site_id: int

    @property
    def number_of_units(self) -> int:
        return len(self.unit_site_ids)


def _validate_linear_polymer(system: MolecularSystem) -> _PolymerLayout:
    if system.representation != "atomistic" or any(
        not isinstance(site, AtomSite) for site in system.topology.sites.values()
    ):
        raise UnsupportedConformationError(
            "Random-walk generation supports atomistic systems only"
        )
    if len(system.topology.connected_components()) != 1:
        raise UnsupportedConformationError(
            "Random-walk generation requires one connected molecular component"
        )
    polymer = system.metadata.get("polymer")
    if not isinstance(polymer, dict) or polymer.get("architecture") != "linear":
        raise UnsupportedConformationError(
            "Random-walk generation requires linear-polymer metadata"
        )
    groups: dict[int, list[int]] = {}
    for site_id, site in system.topology.sites.items():
        repeat_index = site.metadata.get("repeat_unit_index")
        if not isinstance(repeat_index, int) or isinstance(repeat_index, bool):
            raise UnsupportedConformationError(
                f"Site {site_id} lacks an integer repeat_unit_index"
            )
        groups.setdefault(repeat_index, []).append(site_id)
    indices = sorted(groups)
    if indices != list(range(len(indices))):
        raise UnsupportedConformationError(
            "repeat_unit_index values must be contiguous and start at zero"
        )
    expected_units = polymer.get("number_of_repeat_units")
    if expected_units != len(indices):
        raise UnsupportedConformationError(
            "Polymer repeat count does not match atom provenance"
        )
    head_site_id = polymer.get("head_site_id")
    tail_site_id = polymer.get("tail_site_id")
    if (
        head_site_id not in system.topology.sites
        or tail_site_id not in system.topology.sites
    ):
        raise UnsupportedConformationError("Polymer head/tail site IDs are invalid")
    if (
        system.topology.get_site(head_site_id).metadata["repeat_unit_index"] != 0
        or system.topology.get_site(tail_site_id).metadata["repeat_unit_index"]
        != indices[-1]
    ):
        raise UnsupportedConformationError(
            "Polymer head/tail IDs do not identify terminal repeat units"
        )

    connections: dict[tuple[int, int], tuple[int, int]] = {}
    for bond in system.topology.bonds.values():
        first_index = system.topology.get_site(bond.site1).metadata["repeat_unit_index"]
        second_index = system.topology.get_site(bond.site2).metadata[
            "repeat_unit_index"
        ]
        if first_index == second_index:
            continue
        if abs(first_index - second_index) != 1:
            raise UnsupportedConformationError(
                "Non-adjacent repeat units are connected; networks are unsupported"
            )
        if first_index < second_index:
            key, endpoints = (first_index, second_index), (bond.site1, bond.site2)
        else:
            key, endpoints = (second_index, first_index), (bond.site2, bond.site1)
        if key in connections:
            raise UnsupportedConformationError(
                "Adjacent repeat units have multiple connections; branched or "
                "network polymers are unsupported"
            )
        connections[key] = endpoints
    expected_connections = {(index, index + 1) for index in range(len(indices) - 1)}
    if set(connections) != expected_connections:
        raise UnsupportedConformationError(
            "Repeat-unit connectivity is not one finite linear chain"
        )
    incoming = [head_site_id]
    incoming.extend(connections[(index - 1, index)][1] for index in indices[1:])
    outgoing = [connections[(index, index + 1)][0] for index in indices[:-1]]
    outgoing.append(tail_site_id)
    return _PolymerLayout(
        tuple(tuple(groups[index]) for index in indices),
        tuple(incoming),
        tuple(outgoing),
        head_site_id,
        tail_site_id,
    )


def _normalize(vector: NDArray[np.float64]) -> NDArray[np.float64]:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-12:
        raise ConformationGenerationError("Cannot normalize a zero-length vector")
    return vector / norm


def _random_unit_vector(rng: random.Random) -> NDArray[np.float64]:
    z = rng.uniform(-1.0, 1.0)
    azimuth = rng.uniform(-math.pi, math.pi)
    radius = math.sqrt(max(0.0, 1.0 - z * z))
    return np.array([radius * math.cos(azimuth), radius * math.sin(azimuth), z])


def _perpendicular_basis(
    axis: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, reference))) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    first = _normalize(np.cross(axis, reference))
    return first, _normalize(np.cross(axis, first))


def _sample_cone_direction(
    axis: NDArray[np.float64],
    rng: random.Random,
    *,
    center_angle: float,
    half_width: float,
) -> NDArray[np.float64]:
    first, second = _perpendicular_basis(axis)
    polar_angle = rng.uniform(
        max(0.0, center_angle - half_width),
        min(math.pi, center_angle + half_width),
    )
    cosine = math.cos(polar_angle)
    sine = math.sin(polar_angle)
    azimuth = rng.uniform(-math.pi, math.pi)
    return _normalize(
        cosine * axis + sine * (math.cos(azimuth) * first + math.sin(azimuth) * second)
    )


def _rotation_between(
    source: NDArray[np.float64], target: NDArray[np.float64]
) -> NDArray[np.float64]:
    first = _normalize(source)
    second = _normalize(target)
    cross = np.cross(first, second)
    sine = float(np.linalg.norm(cross))
    cosine = float(np.clip(np.dot(first, second), -1.0, 1.0))
    if sine < 1e-12:
        if cosine > 0:
            return np.eye(3)
        axis, _ = _perpendicular_basis(first)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    skew = np.array(
        [
            [0.0, -cross[2], cross[1]],
            [cross[2], 0.0, -cross[0]],
            [-cross[1], cross[0], 0.0],
        ]
    )
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / (sine * sine))
