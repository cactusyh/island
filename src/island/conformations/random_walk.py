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
    bond_angle_degrees: float | None = None
    direction_jitter_degrees: float | None = None
    exclude_one_four: bool = False

    def generate(self, system: MolecularSystem) -> ConformationResult:
        """Generate new coordinates without modifying ``system``."""
        if not isinstance(self.seed, int) or isinstance(self.seed, bool):
            raise ValueError("seed must be an integer")
        if (
            self.bond_angle_degrees is not None
            or self.direction_jitter_degrees is not None
        ):
            raise UnsupportedConformationError(
                "Bond-angle resampling is unsupported; source local geometry is preserved"
            )
        layout = _validate_linear_polymer(system)
        source = {
            site_id: system.coordinates.get(site_id)
            for site_id in system.topology.sites
        }
        _validate_source_geometry(system, source)
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
        unit_rotations: dict[int, NDArray[np.float64]] = {0: np.eye(3)}
        repeat_index = 1

        while repeat_index < layout.number_of_units:
            placed = False
            for trial in range(1, self.retry_policy.attempts_per_unit + 1):
                attempts += 1
                unit_attempts[repeat_index] += 1
                torsion = self.torsion_sampler.sample(
                    rng, repeat_index=repeat_index, trial=trial
                )
                if not math.isfinite(torsion):
                    raise ValueError("Torsion sampler must return a finite angle")
                proposed, rotation = self._propose_unit(
                    layout,
                    source,
                    positions,
                    repeat_index,
                    torsion,
                    unit_rotations,
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
                unit_rotations[repeat_index] = rotation
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
                unit_rotations.pop(index, None)
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
                "accepted_rotation_increments_degrees": {
                    index: math.degrees(angle)
                    for index, angle in sorted(accepted_torsions.items())
                },
                "steric_policy": type(self.steric_policy).__name__,
                "excluded_pairs": (
                    "1-2, 1-3 and 1-4" if self.exclude_one_four else "1-2 and 1-3"
                ),
                "one_four_pairs_checked": not self.exclude_one_four,
                "repeat_units_are_rigid": True,
                "source_local_geometry_preserved": True,
                "coordinate_units": "angstrom",
                "radius_of_gyration_weighting": "uniform_sites",
                "requires_valid_3d_source": True,
                "ring_intersection_check": False,
                "end_to_end_distance": end_to_end_distance,
                "radius_of_gyration": radius_of_gyration,
            },
        )

    def _propose_unit(
        self,
        layout: "_PolymerLayout",
        source: dict[int, NDArray[np.float64]],
        positions: dict[int, NDArray[np.float64]],
        repeat_index: int,
        torsion: float,
        unit_rotations: dict[int, NDArray[np.float64]],
    ) -> tuple[dict[int, NDArray[np.float64]], NDArray[np.float64]]:
        # The previous unit's full frame fixes its outgoing bond direction.
        # Rotating the next unit around that bond also preserves its incoming
        # bond direction. Thus *both* endpoint neighborhoods keep their geometry,
        # even when a stereocenter is itself an inter-repeat attachment atom.
        previous_anchor = layout.outgoing_anchors[repeat_index - 1]
        current_anchor = layout.incoming_anchors[repeat_index]
        previous_rotation = unit_rotations[repeat_index - 1]
        bond_vector = previous_rotation @ (
            source[current_anchor] - source[previous_anchor]
        )
        target_anchor = positions[previous_anchor] + bond_vector
        rotation = _axis_rotation(bond_vector, torsion) @ previous_rotation
        proposed = {
            site_id: rotation @ (source[site_id] - source[current_anchor])
            + target_anchor
            for site_id in layout.unit_site_ids[repeat_index]
        }
        return proposed, rotation


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
    system.validate()
    if system.box is not None and any(system.box.periodic):
        raise UnsupportedConformationError("Periodic systems are unsupported")
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
        if bond.order != 1 or bond.aromatic:
            raise UnsupportedConformationError(
                "Only single, non-aromatic inter-repeat bonds can be rotated"
            )
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


def _axis_rotation(axis: NDArray[np.float64], angle: float) -> NDArray[np.float64]:
    """A proper rotation around an existing bond, never a reflection."""
    axis = _normalize(axis)
    x, y, z = axis
    skew = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    cosine = math.cos(angle)
    return (
        cosine * np.eye(3)
        + (1 - cosine) * np.outer(axis, axis)
        + math.sin(angle) * skew
    )


def _validate_source_geometry(
    system: MolecularSystem, source: dict[int, NDArray[np.float64]]
) -> None:
    """Reject depiction coordinates and degenerate tetrahedral neighborhoods.

    This geometric check does not infer absolute CIP labels or certify energetic
    quality. Valid source stereochemistry is a precondition; rigid bond rotations
    preserve it. Coordinate-derived CIP verification belongs in chemistry tests.
    """
    if system.metadata["polymer"].get("coordinates") == "rdkit_2d":
        raise UnsupportedConformationError(
            "Random walk requires valid 3D source geometry, not a 2D depiction; "
            "build with generate_3d=True. Local 3D template generation is not implemented"
        )
    for site_id, point in source.items():
        if not np.all(np.isfinite(point)):
            raise UnsupportedConformationError(
                f"Non-finite source coordinates at site {site_id}"
            )
    for bond in system.topology.bonds.values():
        if np.linalg.norm(source[bond.site1] - source[bond.site2]) < 1e-6:
            raise UnsupportedConformationError(f"Zero-length source bond {bond.key}")
    for site_id, site in system.topology.sites.items():
        neighbors = sorted(system.topology.neighbors(site_id))
        tetrahedral = (
            site.metadata.get("hybridization") == "SP3"
            or site.metadata.get("chiral_tag")
            in {"CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW"}
            or site.metadata.get("controllable_stereocenter")
            or (site.atomic_number == 6 and len(neighbors) == 4)
        )
        if not tetrahedral or len(neighbors) < 3:
            continue
        points = np.array([source[neighbor] for neighbor in neighbors])
        if len(neighbors) == 4:
            vectors = points[:3] - points[3]
        elif len(neighbors) == 3:
            vectors = points - source[site_id]
        else:
            raise UnsupportedConformationError(
                f"Invalid tetrahedral coordination at site {site_id}"
            )
        norms = np.linalg.norm(vectors, axis=1)
        if np.any(norms < 1e-6) or abs(np.linalg.det(vectors / norms[:, None])) < 1e-3:
            raise UnsupportedConformationError(
                f"Degenerate tetrahedral source geometry at site {site_id}; "
                "rigid rotations cannot create a valid 3D stereocenter"
            )
