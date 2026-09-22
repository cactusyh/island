"""Initial polymer conformations assembled from short local 3D templates."""

import math
import random
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem

from island.chemistry import to_rdkit
from island.conformations.base import ConformationGenerator
from island.conformations.local_validation import (
    authoritative_cip_assignments,
    transfer_fragment_chirality,
    validate_coordinate_stereochemistry,
    validate_explicit_hydrogens,
)
from island.conformations.random_walk import (
    RetryPolicy,
    _axis_rotation,
    _normalize,
    _validate_linear_polymer,
)
from island.conformations.result import ConformationResult
from island.conformations.sterics import (
    FixedDistanceStericPolicy,
    StericPolicy,
    build_excluded_pairs,
    find_clashes,
    minimum_nonbonded_distance,
)
from island.conformations.torsions import TorsionSampler, UniformTorsionSampler
from island.core import Coordinates, MolecularSystem
from island.exceptions import (
    ConformationGenerationError,
    UnsupportedConformationError,
)

_FINAL_SITE_ID = "_island_final_site_id"
_CAP_ROLE = "_island_attachment_cap_role"


class ConnectionBondLengthPolicy:
    """Determine traceable inter-template connection-bond lengths."""

    def length(
        self, system: MolecularSystem, site1: int, site2: int, bond_order: float
    ) -> float:
        raise NotImplementedError


class CovalentRadiiBondLengthPolicy(ConnectionBondLengthPolicy):
    """Estimate bond lengths from RDKit covalent radii and bond order.

    Single bonds use the sum of elemental covalent radii. Higher bond orders are
    outside the current polymer-builder scope and are rejected explicitly.
    """

    def length(
        self, system: MolecularSystem, site1: int, site2: int, bond_order: float
    ) -> float:
        if bond_order != 1.0:
            raise UnsupportedConformationError(
                "Local-template assembly supports single inter-repeat bonds only"
            )
        first = system.topology.get_site(site1)
        second = system.topology.get_site(site2)
        if first.atomic_number is None or second.atomic_number is None:
            raise UnsupportedConformationError(
                "Connection atoms require atomic numbers for covalent-radius lengths"
            )
        table = Chem.GetPeriodicTable()
        length = table.GetRcovalent(first.atomic_number) + table.GetRcovalent(
            second.atomic_number
        )
        if not math.isfinite(length) or length <= 0:
            raise UnsupportedConformationError(
                f"No valid covalent-radius bond length for sites {site1}-{site2}"
            )
        return float(length)


@dataclass(frozen=True)
class LocalTemplate:
    """Rigid local 3D geometry and explicit attachment frames for one repeat."""

    repeat_index: int
    site_coordinates: dict[int, NDArray[np.float64]]
    head_site_id: int
    tail_site_id: int
    incoming_direction: NDArray[np.float64] | None
    outgoing_direction: NDArray[np.float64] | None
    embedded_atom_count: int
    template_smiles: str


@dataclass
class LocalTemplateConformationGenerator(ConformationGenerator):
    """Embed short repeat templates and assemble a self-avoiding 3D chain."""

    template_seed: int = 2026
    assembly_seed: int = 2026
    torsion_sampler: TorsionSampler = field(default_factory=UniformTorsionSampler)
    steric_policy: StericPolicy = field(
        default_factory=lambda: FixedDistanceStericPolicy(0.8)
    )
    retry_policy: RetryPolicy = field(default_factory=RetryPolicy)
    bond_length_policy: ConnectionBondLengthPolicy = field(
        default_factory=CovalentRadiiBondLengthPolicy
    )
    exclude_one_four: bool = False

    def generate(self, system: MolecularSystem) -> ConformationResult:
        """Return initial 3D coordinates without mutating the input system."""
        _validate_seed(self.template_seed, "template_seed")
        _validate_seed(self.assembly_seed, "assembly_seed")
        layout = _validate_linear_polymer(system)
        validate_explicit_hydrogens(system)
        expected_cip = authoritative_cip_assignments(system)
        templates = _generate_templates(system, layout, self.template_seed)
        positions, diagnostics = self._assemble(system, layout, templates)
        stereo_validation = validate_coordinate_stereochemistry(
            system, positions, expected_cip
        )
        coordinates = Coordinates(positions)
        coordinates.validate(system.topology)
        excluded = build_excluded_pairs(system, exclude_one_four=self.exclude_one_four)
        minimum_distance = minimum_nonbonded_distance(positions, excluded)
        array = np.array([positions[site_id] for site_id in system.topology.sites])
        center = array.mean(axis=0)
        radius_of_gyration = float(
            np.sqrt(np.mean(np.sum((array - center) ** 2, axis=1)))
        )
        end_to_end = float(
            np.linalg.norm(
                positions[layout.tail_site_id] - positions[layout.head_site_id]
            )
        )
        return ConformationResult(
            coordinates=coordinates,
            method="local_templates_self_avoiding_random_walk",
            seed=self.assembly_seed,
            success=True,
            attempts=diagnostics["attempts"],
            rejected_trials=diagnostics["rejected_trials"],
            rollback_count=diagnostics["rollback_count"],
            minimum_nonbonded_distance=minimum_distance,
            metadata={
                "coordinate_source": "local_templates_etkdg_incremental",
                "template_seed": self.template_seed,
                "assembly_seed": self.assembly_seed,
                "template_count": len(templates),
                "maximum_embedded_atom_count": max(
                    template.embedded_atom_count for template in templates
                ),
                "full_polymer_atom_count": system.number_of_sites,
                "full_polymer_embedded": False,
                "attachment_frames": "explicit_methyl_like_carbon_caps",
                "bond_length_policy": type(self.bond_length_policy).__name__,
                "accepted_rotation_increments_degrees": diagnostics[
                    "accepted_rotations"
                ],
                "attempts_per_repeat_unit": diagnostics["unit_attempts"],
                "steric_policy": type(self.steric_policy).__name__,
                "excluded_pairs": (
                    "1-2, 1-3 and 1-4" if self.exclude_one_four else "1-2 and 1-3"
                ),
                "one_four_pairs_checked": not self.exclude_one_four,
                "repeat_units_are_rigid": True,
                "source_local_geometry_preserved": True,
                "coordinate_units": "angstrom",
                "stereochemistry_validation": stereo_validation,
                "stereochemistry_validated_from_3d": (
                    stereo_validation["status"] == "validated"
                ),
                "ring_intersection_check": False,
                "energy_optimized": False,
                "end_to_end_distance": end_to_end,
                "radius_of_gyration": radius_of_gyration,
            },
        )

    def _assemble(
        self, system: MolecularSystem, layout: object, templates: list[LocalTemplate]
    ) -> tuple[dict[int, NDArray[np.float64]], dict[str, object]]:
        first = templates[0]
        origin = first.site_coordinates[first.head_site_id]
        positions = {
            site_id: point - origin for site_id, point in first.site_coordinates.items()
        }
        rotations: dict[int, NDArray[np.float64]] = {0: np.eye(3)}
        rng = random.Random(self.assembly_seed)
        excluded = build_excluded_pairs(system, exclude_one_four=self.exclude_one_four)
        attempts = rejected = rollbacks = 0
        unit_attempts = [0] * len(templates)
        accepted_rotations: dict[int, float] = {}
        repeat_index = 1
        while repeat_index < len(templates):
            placed = False
            for trial in range(1, self.retry_policy.attempts_per_unit + 1):
                attempts += 1
                unit_attempts[repeat_index] += 1
                angle = self.torsion_sampler.sample(
                    rng, repeat_index=repeat_index, trial=trial
                )
                if not math.isfinite(angle):
                    raise ValueError("Torsion sampler must return a finite angle")
                proposed, rotation = self._propose(
                    system, layout, templates, positions, rotations, repeat_index, angle
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
                    excluded,
                )
                if clashes:
                    rejected += 1
                    continue
                positions.update(proposed)
                rotations[repeat_index] = rotation
                accepted_rotations[repeat_index] = math.degrees(angle)
                repeat_index += 1
                placed = True
                break
            if placed:
                continue
            if rollbacks >= self.retry_policy.max_rollbacks:
                raise ConformationGenerationError(
                    "Local-template assembly exhausted retries at repeat "
                    f"{repeat_index} after {attempts} attempts and {rollbacks} rollbacks",
                    repeat_index=repeat_index,
                    attempts=attempts,
                    rejected_trials=rejected,
                    rollback_count=rollbacks,
                )
            rollbacks += 1
            rollback_start = max(1, repeat_index - self.retry_policy.rollback_units)
            for index in range(rollback_start, repeat_index):
                for site_id in layout.unit_site_ids[index]:
                    positions.pop(site_id, None)
                rotations.pop(index, None)
                accepted_rotations.pop(index, None)
            repeat_index = rollback_start
        return positions, {
            "attempts": attempts,
            "rejected_trials": rejected,
            "rollback_count": rollbacks,
            "unit_attempts": unit_attempts,
            "accepted_rotations": dict(sorted(accepted_rotations.items())),
        }

    def _propose(
        self,
        system: MolecularSystem,
        layout: object,
        templates: list[LocalTemplate],
        positions: dict[int, NDArray[np.float64]],
        rotations: dict[int, NDArray[np.float64]],
        repeat_index: int,
        angle: float,
    ) -> tuple[dict[int, NDArray[np.float64]], NDArray[np.float64]]:
        previous = templates[repeat_index - 1]
        current = templates[repeat_index]
        if previous.outgoing_direction is None or current.incoming_direction is None:
            raise UnsupportedConformationError(
                f"Missing attachment frame at repeat junction {repeat_index - 1}-{repeat_index}"
            )
        previous_rotation = rotations[repeat_index - 1]
        direction = _normalize(previous_rotation @ previous.outgoing_direction)
        alignment = _rotation_between(current.incoming_direction, direction)
        rotation = _axis_rotation(direction, angle) @ alignment
        previous_anchor = layout.outgoing_anchors[repeat_index - 1]
        current_anchor = layout.incoming_anchors[repeat_index]
        bond = system.topology.bonds[tuple(sorted((previous_anchor, current_anchor)))]
        length = self.bond_length_policy.length(
            system, previous_anchor, current_anchor, float(bond.order)
        )
        if not math.isfinite(length) or length <= 0:
            raise UnsupportedConformationError(
                f"Invalid connection-bond length {length!r} at junction "
                f"{repeat_index - 1}-{repeat_index}"
            )
        target = positions[previous_anchor] + direction * length
        current_origin = current.site_coordinates[current_anchor]
        proposed = {
            site_id: rotation @ (point - current_origin) + target
            for site_id, point in current.site_coordinates.items()
        }
        return proposed, rotation


def _generate_templates(
    system: MolecularSystem, layout: object, template_seed: int
) -> list[LocalTemplate]:
    converted = to_rdkit(system)
    source = converted.mol
    source.RemoveAllConformers()
    templates = []
    for repeat_index, unit_site_ids in enumerate(layout.unit_site_ids):
        templates.append(
            _embed_template(
                system,
                source,
                converted.site_id_to_rdkit_index,
                layout,
                repeat_index,
                unit_site_ids,
                template_seed + repeat_index,
            )
        )
    return templates


def _embed_template(
    system: MolecularSystem,
    source: Chem.Mol,
    site_to_source_index: dict[int, int],
    layout: object,
    repeat_index: int,
    unit_site_ids: tuple[int, ...],
    seed: int,
) -> LocalTemplate:
    editable = Chem.RWMol()
    site_to_template: dict[int, int] = {}
    for site_id in unit_site_ids:
        atom = Chem.Atom(source.GetAtomWithIdx(site_to_source_index[site_id]))
        atom.SetIntProp(_FINAL_SITE_ID, site_id)
        site_to_template[site_id] = editable.AddAtom(atom)
    for bond in system.topology.bonds.values():
        if bond.site1 not in site_to_template or bond.site2 not in site_to_template:
            continue
        editable.AddBond(
            site_to_template[bond.site1],
            site_to_template[bond.site2],
            _bond_type(bond.order, bond.aromatic),
        )
        copied = editable.GetBondBetweenAtoms(
            site_to_template[bond.site1], site_to_template[bond.site2]
        )
        copied.SetIsAromatic(bond.aromatic)
    template_to_source = {
        template_index: site_to_source_index[site_id]
        for site_id, template_index in site_to_template.items()
    }
    template_to_site = {
        template_index: site_id for site_id, template_index in site_to_template.items()
    }
    cap_to_source: dict[int, int] = {}

    incoming_cap = None
    if repeat_index > 0:
        incoming_cap = _add_cap(
            editable,
            site_to_template[layout.incoming_anchors[repeat_index]],
            "head",
            isotope=801,
        )
        cap_to_source[incoming_cap] = site_to_source_index[
            layout.outgoing_anchors[repeat_index - 1]
        ]
    outgoing_cap = None
    if repeat_index < len(layout.unit_site_ids) - 1:
        outgoing_cap = _add_cap(
            editable,
            site_to_template[layout.outgoing_anchors[repeat_index]],
            "tail",
            isotope=802,
        )
        cap_to_source[outgoing_cap] = site_to_source_index[
            layout.incoming_anchors[repeat_index + 1]
        ]
    fragment = editable.GetMol()
    transfer_fragment_chirality(
        source,
        fragment,
        template_to_source,
        cap_to_source,
        template_to_site,
        repeat_index,
    )
    try:
        Chem.SanitizeMol(fragment)
        fragment = Chem.AddHs(fragment)
    except Exception as error:
        raise UnsupportedConformationError(
            f"Cannot sanitize local template for repeat {repeat_index}"
        ) from error
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    try:
        status = AllChem.EmbedMolecule(fragment, parameters)
    except Exception as error:
        raise ConformationGenerationError(
            f"Local template embedding failed for repeat {repeat_index}",
            repeat_index=repeat_index,
        ) from error
    if status != 0:
        raise ConformationGenerationError(
            f"Local template embedding failed for repeat {repeat_index}",
            repeat_index=repeat_index,
        )
    conformer = fragment.GetConformer()
    for site_id, template_index in site_to_template.items():
        mapped = fragment.GetAtomWithIdx(template_index)
        if (
            not mapped.HasProp(_FINAL_SITE_ID)
            or mapped.GetIntProp(_FINAL_SITE_ID) != site_id
        ):
            raise ConformationGenerationError(
                f"Local template mapping changed for site {site_id}",
                repeat_index=repeat_index,
            )
    coordinates = {
        site_id: _point(conformer, template_index)
        for site_id, template_index in site_to_template.items()
    }
    incoming_direction = (
        _normalize(
            coordinates[layout.incoming_anchors[repeat_index]]
            - _point(conformer, incoming_cap)
        )
        if incoming_cap is not None
        else None
    )
    outgoing_direction = (
        _normalize(
            _point(conformer, outgoing_cap)
            - coordinates[layout.outgoing_anchors[repeat_index]]
        )
        if outgoing_cap is not None
        else None
    )
    return LocalTemplate(
        repeat_index=repeat_index,
        site_coordinates=coordinates,
        head_site_id=layout.incoming_anchors[repeat_index],
        tail_site_id=layout.outgoing_anchors[repeat_index],
        incoming_direction=incoming_direction,
        outgoing_direction=outgoing_direction,
        embedded_atom_count=fragment.GetNumAtoms(),
        template_smiles=Chem.MolToSmiles(fragment, isomericSmiles=True),
    )


def _add_cap(
    editable: Chem.RWMol, anchor_index: int, role: str, *, isotope: int
) -> int:
    cap = Chem.Atom(6)
    cap.SetIsotope(isotope)
    cap.SetProp(_CAP_ROLE, role)
    cap_index = editable.AddAtom(cap)
    editable.AddBond(anchor_index, cap_index, Chem.BondType.SINGLE)
    return cap_index


def _bond_type(order: float | None, aromatic: bool) -> Chem.BondType:
    if aromatic:
        return Chem.BondType.AROMATIC
    mapping = {
        1.0: Chem.BondType.SINGLE,
        2.0: Chem.BondType.DOUBLE,
        3.0: Chem.BondType.TRIPLE,
    }
    try:
        return mapping[float(order)]
    except (KeyError, TypeError, ValueError) as error:
        raise UnsupportedConformationError(
            f"Unsupported local-template bond order: {order!r}"
        ) from error


def _point(conformer: Chem.Conformer, atom_index: int) -> NDArray[np.float64]:
    point = conformer.GetAtomPosition(atom_index)
    return np.array([point.x, point.y, point.z], dtype=float)


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
        reference = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(first, reference))) > 0.9:
            reference = np.array([0.0, 1.0, 0.0])
        axis = _normalize(np.cross(first, reference))
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    skew = np.array(
        [
            [0.0, -cross[2], cross[1]],
            [cross[2], 0.0, -cross[0]],
            [-cross[1], cross[0], 0.0],
        ]
    )
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / (sine * sine))


def _validate_seed(seed: int, name: str) -> None:
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise ValueError(f"{name} must be an integer")
