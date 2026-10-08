"""Force-field-neutral multichain packing and periodic-box preparation.

This module moves complete molecular units as rigid bodies. It owns no
chemical reaction rules, atom typing, charges, force-field parameters or
periodic energy evaluation.
"""

import random
from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
from math import cos, pi, sin, sqrt
from pathlib import Path

import numpy as np

from island.core import BeadSite, Coordinates, MolecularSystem, SimulationBox, Topology
from island.exceptions import ValidationError
from island.graph import FinalChemicalGraph, final_graph
from island.graph.final import (
    TRANSFORMATION_SCHEMA as GRAPH_TRANSFORMATION_SCHEMA,
)
from island.graph.final import (
    pack,
    unpack,
)

CONFIG_SCHEMA = "island_periodic_packing_config_v1"
PLAN_SCHEMA = "island_periodic_packing_plan_v1"
PACKING_TRANSFORMATION_SCHEMA = "island_periodic_packing_transformation_v1"
PLACEMENT_MODE = "rigid_body_seeded_random"
BOX_SHAPE = "orthorhombic"
AVOGADRO = 6.02214076e23


def _require(condition, message):
    if not condition:
        raise ValidationError(message)


def _finite(value, label, *, positive=False):
    _require(
        type(value) in (int, float)
        and not isinstance(value, bool)
        and np.isfinite(value),
        f"{label} must be finite and nonboolean",
    )
    if positive:
        _require(value > 0, f"{label} must be positive")


def _text(value, label):
    _require(type(value) is str and value.strip(), f"{label} must be nonempty text")


def _digest(value):
    return sha256(pack(value).encode()).hexdigest()


def _identity_payload(value):
    return _digest(value)


def _system_identity(system):
    """Identity of coordinates/topology and source metadata, excluding runtime links."""
    system.validate()
    payload = system.to_dict()
    metadata = deepcopy(payload.get("metadata", {}))
    metadata.pop("packing", None)
    metadata.pop("final_graph_transformations", None)
    payload["metadata"] = metadata
    return _identity_payload(payload)


def _evidence(value):
    _require(
        type(value) in (tuple, list)
        and all(type(item) is str and item.strip() for item in value),
        "Evidence must contain nonempty references",
    )
    return tuple(value)


def _tuple3(value, label, *, boolean=False):
    _require(
        type(value) in (tuple, list) and len(value) == 3,
        f"{label} must have three values",
    )
    if boolean:
        _require(
            all(type(item) is bool for item in value), f"{label} must contain booleans"
        )
    else:
        for item in value:
            _finite(item, label, positive=True)
    return tuple(value)


def _config_payload(config):
    return {
        "schema": CONFIG_SCHEMA,
        "target_density": config.target_density,
        "box_lengths": None if config.box_lengths is None else list(config.box_lengths),
        "periodic": list(config.periodic),
        "seed": config.seed,
        "rotation_mode": config.rotation_mode,
        "minimum_interunit_distance": config.minimum_interunit_distance,
        "maximum_attempts": config.maximum_attempts,
        "density_tolerance": config.density_tolerance,
        "placement_mode": PLACEMENT_MODE,
        "box_shape": BOX_SHAPE,
    }


@dataclass(frozen=True)
class PeriodicPackingConfig:
    target_density: float | None = None
    box_lengths: tuple[float, float, float] | None = None
    periodic: tuple[bool, bool, bool] = (True, True, True)
    seed: int = 2026
    rotation_mode: str = "random_uniform"
    minimum_interunit_distance: float = 1.5
    maximum_attempts: int = 10000
    density_tolerance: float = 1e-8

    def __post_init__(self):
        _require(
            (self.target_density is None) != (self.box_lengths is None),
            "Specify exactly one of target_density or box_lengths",
        )
        if self.target_density is not None:
            _finite(self.target_density, "Target density", positive=True)
        if self.box_lengths is not None:
            object.__setattr__(
                self, "box_lengths", _tuple3(self.box_lengths, "Box lengths")
            )
        object.__setattr__(
            self, "periodic", _tuple3(self.periodic, "Periodic boundary", boolean=True)
        )
        _require(
            type(self.seed) is int and not isinstance(self.seed, bool),
            "Seed must be an integer",
        )
        _require(
            self.rotation_mode in ("random_uniform", "none"),
            "Unsupported rotation mode",
        )
        _finite(self.minimum_interunit_distance, "Minimum interunit distance")
        _require(
            self.minimum_interunit_distance >= 0,
            "Minimum interunit distance cannot be negative",
        )
        _require(
            type(self.maximum_attempts) is int and self.maximum_attempts > 0,
            "Maximum attempts must be positive",
        )
        _finite(self.density_tolerance, "Density tolerance", positive=True)

    @property
    def payload(self):
        return deepcopy(_config_payload(self))

    @property
    def identity(self):
        return _digest(self.payload)

    @property
    def json_text(self):
        payload = self.payload
        payload["identity"] = _digest(payload)
        return pack(payload)

    @classmethod
    def from_json(cls, json_text):
        return _config_from_payload(unpack(json_text))


def _config_from_payload(payload):
    required = {
        "schema",
        "target_density",
        "box_lengths",
        "periodic",
        "seed",
        "rotation_mode",
        "minimum_interunit_distance",
        "maximum_attempts",
        "density_tolerance",
        "placement_mode",
        "box_shape",
        "identity",
    }
    _require(
        type(payload) is dict and set(payload) == required, "Malformed packing config"
    )
    _require(
        payload["schema"] == CONFIG_SCHEMA
        and payload["placement_mode"] == PLACEMENT_MODE
        and payload["box_shape"] == BOX_SHAPE,
        "Unsupported packing config schema",
    )
    identity = payload.pop("identity")
    _require(identity == _digest(payload), "Packing config identity mismatch")
    return PeriodicPackingConfig(
        target_density=payload["target_density"],
        box_lengths=None
        if payload["box_lengths"] is None
        else tuple(payload["box_lengths"]),
        periodic=tuple(payload["periodic"]),
        seed=payload["seed"],
        rotation_mode=payload["rotation_mode"],
        minimum_interunit_distance=payload["minimum_interunit_distance"],
        maximum_attempts=payload["maximum_attempts"],
        density_tolerance=payload["density_tolerance"],
    )


def _rotation_matrix(rng, mode):
    if mode == "none":
        return np.eye(3)
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    q = (
        sqrt(1 - u1) * sin(2 * pi * u2),
        sqrt(1 - u1) * cos(2 * pi * u2),
        sqrt(u1) * sin(2 * pi * u3),
        sqrt(u1) * cos(2 * pi * u3),
    )
    x, y, z, w = q
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def _minimum_distance(candidate, existing, lengths, periodic, threshold):
    if not existing:
        return float("inf")
    for point in candidate:
        for other in existing:
            delta = np.abs(point - other)
            for axis in range(3):
                if periodic[axis]:
                    delta[axis] = min(delta[axis], lengths[axis] - delta[axis])
            if float(np.linalg.norm(delta)) < 0:
                return 0.0
            if float(np.linalg.norm(delta)) < threshold:
                return float(np.linalg.norm(delta))
    return float("inf")


def _unit_mass(system):
    return sum(float(site.mass) for site in system.topology.sites.values())


def _box_for(config, units):
    if config.box_lengths is not None:
        return config.box_lengths
    mass = sum(_unit_mass(system) for _, system, _ in units)
    _require(mass > 0, "Positive total mass is required for density packing")
    volume = mass / (config.target_density * AVOGADRO) * 1e24
    edge = volume ** (1 / 3)
    return (edge, edge, edge)


def _validate_units(units):
    _require(
        type(units) in (tuple, list) and units, "At least one packing unit is required"
    )
    seen = set()
    normalized = []
    for item in units:
        _require(
            type(item) in (tuple, list) and len(item) == 3,
            "Each unit must be (label, system, graph)",
        )
        label, system, graph = item
        _text(label, "Unit label")
        _require(label not in seen, f"Duplicate unit label: {label}")
        seen.add(label)
        _require(
            isinstance(system, MolecularSystem),
            "Packing unit system must be MolecularSystem",
        )
        _require(
            isinstance(graph, FinalChemicalGraph),
            "Packing unit graph must be FinalChemicalGraph",
        )
        system.validate()
        graph.validate_integrity(system)
        _require(
            not graph.periodic and system.box is None,
            "Periodic input units are unsupported; pack nonperiodic units",
        )
        _require(
            all(
                np.isfinite(system.coordinates.get(i)).all()
                for i in system.topology.sites
            ),
            "Unit coordinates must be finite",
        )
        normalized.append((label, system, graph))
    return tuple(normalized)


def _remap_metadata(metadata, mapping):
    metadata = deepcopy(metadata)
    crosslinks = []
    for row in metadata.get("crosslinks", []):
        _require(type(row) is dict and "sites" in row, "Malformed crosslink provenance")
        copied = deepcopy(row)
        copied["sites"] = [mapping[int(site)] for site in row["sites"]]
        crosslinks.append(copied)
    if "crosslinks" in metadata:
        metadata["crosslinks"] = crosslinks
    metadata.pop("final_graph_transformations", None)
    return metadata


def _build_output(units, config, placements, box_lengths, plan_identity=None):
    topology = Topology()
    coordinates = {}
    molecule_membership = {}
    graph_ids, system_ids, remaps, molecule_remaps, unit_records = [], [], [], [], []
    crosslinks = []
    polymer_units = []
    next_id = 1
    existing = set()
    for placement, (label, system, graph) in zip(placements, units):
        source_ids = sorted(system.topology.sites)
        mapping = {source_id: next_id + n for n, source_id in enumerate(source_ids)}
        next_id += len(source_ids)
        remaps.append({str(k): v for k, v in mapping.items()})
        graph_ids.append(graph.identity)
        system_ids.append(_system_identity(system))
        source_membership = {
            int(k): v for k, v in graph.payload["molecule_membership"].items()
        }
        member_map = {
            str(source_id): f"{label}:{source_membership[source_id]}"
            for source_id in source_ids
        }
        molecule_remaps.append(member_map)
        unit_records.append(
            {
                "label": label,
                "graph_identity": graph.identity,
                "system_identity": _system_identity(system),
            }
        )
        polymer_units.append(
            {
                "label": label,
                "repeat_provenance": deepcopy(graph.payload["repeat_provenance"]),
                "crosslink_provenance": deepcopy(graph.payload["crosslink_provenance"]),
            }
        )
        for source_id in source_ids:
            site = deepcopy(system.topology.sites[source_id])
            new_id = mapping[source_id]
            site.id = new_id
            if isinstance(site, BeadSite):
                site.mapped_atom_ids = tuple(
                    mapping.get(i, i) for i in site.mapped_atom_ids
                )
            topology.add_site(site)
            coordinates[new_id] = placement["coordinates"][str(source_id)]
            molecule_membership[new_id] = member_map[str(source_id)]
        for bond in system.topology.bonds.values():
            topology.add_bond(
                mapping[bond.site1],
                mapping[bond.site2],
                order=bond.order,
                aromatic=bond.aromatic,
            )
        for improper in system.topology.impropers:
            topology.impropers.append(
                type(improper)(
                    **{
                        field: mapping[getattr(improper, field)]
                        for field in improper.__dataclass_fields__
                    }
                )
            )
        for row in system.metadata.get("crosslinks", []):
            copied = deepcopy(row)
            copied["sites"] = [mapping[int(site)] for site in row["sites"]]
            crosslinks.append(copied)
        existing.update(mapping.values())
    metadata = {
        "polymer": {"packing_units": polymer_units},
        "crosslinks": crosslinks,
        "packing": {
            "plan_identity": plan_identity,
            "unit_labels": [item[0] for item in units],
        }
        if plan_identity
        else {"unit_labels": [item[0] for item in units]},
    }
    topology.rebuild_derived_interactions()
    result = MolecularSystem(
        topology,
        Coordinates(coordinates),
        box=SimulationBox(*box_lengths, periodic=config.periodic),
        metadata=metadata,
    )
    graph = final_graph(result, molecule_membership=molecule_membership)
    return result, graph, graph_ids, system_ids, remaps, molecule_remaps, unit_records


def _place_units(units, config, box_lengths):
    rng = random.Random(config.seed)
    lengths = np.asarray(box_lengths, dtype=float)
    existing = []
    placements = []
    rejected = 0
    for label, system, graph in units:
        ids = sorted(system.topology.sites)
        source = np.array([system.coordinates.get(i) for i in ids], dtype=float)
        center = source.mean(axis=0)
        centered = source - center
        placed = None
        attempts = 0
        while attempts < config.maximum_attempts:
            attempts += 1
            rotation = _rotation_matrix(rng, config.rotation_mode)
            translation = np.array([rng.random() * lengths[i] for i in range(3)])
            transformed = (centered @ rotation.T) + translation
            wrapped = np.mod(transformed, lengths)
            if (
                _minimum_distance(
                    wrapped,
                    existing,
                    lengths,
                    config.periodic,
                    config.minimum_interunit_distance,
                )
                >= config.minimum_interunit_distance
            ):
                placed = {
                    "unit_label": label,
                    "source_ids": ids,
                    "rotation": rotation.tolist(),
                    "translation": translation.tolist(),
                    "coordinates": {
                        str(i): point.tolist() for i, point in zip(ids, wrapped)
                    },
                    "attempts": attempts - 1,
                }
                existing.extend(wrapped)
                break
            rejected += 1
        _require(
            placed is not None,
            f"Unable to place unit {label} after {config.maximum_attempts} attempts",
        )
        placements.append(placed)
    return placements, rejected


def _validate_placements(units, payload, config):
    """Recompute serialized rigid placements from the source coordinates."""
    lengths = np.asarray(payload["box_lengths"], dtype=float)
    existing = []
    next_id = 1
    for index, ((label, system, graph), placement) in enumerate(
        zip(units, payload["accepted_placements"])
    ):
        ids = sorted(system.topology.sites)
        _require(placement["unit_label"] == label, "Placement unit order changed")
        _require(
            placement["source_ids"] == ids, f"Placement source IDs changed for {label}"
        )
        expected_mapping = {
            str(source_id): next_id + offset for offset, source_id in enumerate(ids)
        }
        _require(
            payload["site_id_remapping"][index] == expected_mapping,
            f"Site-ID remapping changed for {label}",
        )
        next_id += len(ids)
        source_membership = {
            int(key): value
            for key, value in graph.payload["molecule_membership"].items()
        }
        expected_membership = {
            str(source_id): f"{label}:{source_membership[source_id]}"
            for source_id in ids
        }
        _require(
            payload["molecule_membership_remapping"][index] == expected_membership,
            f"Molecule membership remapping changed for {label}",
        )
        source = np.array([system.coordinates.get(i) for i in ids], dtype=float)
        center = source.mean(axis=0)
        rotation = np.asarray(placement["rotation"], dtype=float)
        translation = np.asarray(placement["translation"], dtype=float)
        _require(
            rotation.shape == (3, 3) and translation.shape == (3,),
            "Malformed rigid placement",
        )
        _require(
            np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10)
            and np.isclose(np.linalg.det(rotation), 1, atol=1e-10),
            "Rigid rotation is not orthonormal",
        )
        expected = np.mod((source - center) @ rotation.T + translation, lengths)
        coordinates = np.array(
            [placement["coordinates"][str(i)] for i in ids], dtype=float
        )
        _require(
            coordinates.shape == (len(ids), 3) and np.all(np.isfinite(coordinates)),
            "Malformed packed coordinates",
        )
        _require(
            np.allclose(expected, coordinates, atol=1e-10, rtol=0),
            f"Placement coordinates disagree for {label}",
        )
        _require(
            _minimum_distance(
                coordinates,
                existing,
                lengths,
                config.periodic,
                config.minimum_interunit_distance,
            )
            >= config.minimum_interunit_distance,
            f"Minimum interunit distance violated by {label}",
        )
        existing.extend(coordinates)


def _plan_payload(plan):
    return unpack(plan.json_text)


@dataclass(frozen=True)
class PeriodicPackingPlan:
    json_text: str

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        payload = self.payload
        _require(
            payload["identity"]
            == _digest({k: v for k, v in payload.items() if k != "identity"}),
            "Packing plan identity mismatch",
        )
        return payload["identity"]

    def validate_integrity(self, units=None):
        payload = self.payload
        required = {
            "schema",
            "config",
            "config_identity",
            "input_graph_identities",
            "input_system_identities",
            "unit_labels",
            "site_id_remapping",
            "molecule_membership_remapping",
            "box_lengths",
            "box_shape",
            "periodic",
            "placement_mode",
            "minimum_interunit_distance",
            "seed",
            "accepted_placements",
            "rejected_attempt_count",
            "provenance",
            "evidence",
            "output_system_identity",
            "output_graph_identity",
            "identity",
        }
        _require(
            type(payload) is dict
            and set(payload) == required
            and payload["schema"] == PLAN_SCHEMA,
            "Malformed packing plan",
        )
        _require(
            payload["identity"]
            == _digest({k: v for k, v in payload.items() if k != "identity"}),
            "Packing plan identity mismatch",
        )
        config_payload = deepcopy(payload["config"])
        config = _config_from_payload(config_payload)
        _require(
            config.identity == payload["config_identity"],
            "Packing config identity mismatch",
        )
        _require(
            payload["placement_mode"] == PLACEMENT_MODE
            and payload["seed"] == config.seed,
            "Packing placement policy mismatch",
        )
        _tuple3(payload["box_lengths"], "Box lengths")
        _require(payload["box_shape"] == BOX_SHAPE, "Unsupported packing box shape")
        _tuple3(payload["periodic"], "Periodic boundary", boolean=True)
        _require(
            tuple(payload["periodic"]) == config.periodic,
            "Periodic boundary differs from config",
        )
        if config.box_lengths is not None:
            _require(
                tuple(payload["box_lengths"]) == config.box_lengths,
                "Box lengths differ from config",
            )
        _finite(payload["minimum_interunit_distance"], "Minimum interunit distance")
        _require(
            payload["minimum_interunit_distance"] == config.minimum_interunit_distance,
            "Minimum-distance policy differs from config",
        )
        _require(
            type(payload["rejected_attempt_count"]) is int
            and payload["rejected_attempt_count"] >= 0,
            "Invalid rejected attempt count",
        )
        _text(payload["provenance"], "Packing provenance")
        _evidence(payload["evidence"])
        labels = payload["unit_labels"]
        _require(
            type(labels) is list and labels and len(labels) == len(set(labels)),
            "Packing unit labels must be ordered and unique",
        )
        _require(
            len(payload["input_graph_identities"])
            == len(labels)
            == len(payload["input_system_identities"])
            == len(payload["site_id_remapping"])
            == len(payload["molecule_membership_remapping"])
            == len(payload["accepted_placements"]),
            "Packing unit record length mismatch",
        )
        for placement in payload["accepted_placements"]:
            _require(
                set(placement)
                == {
                    "unit_label",
                    "source_ids",
                    "rotation",
                    "translation",
                    "coordinates",
                    "attempts",
                },
                "Malformed accepted placement",
            )
            _require(placement["unit_label"] in labels, "Unknown placement unit label")
            _require(
                type(placement["source_ids"]) is list
                and placement["source_ids"] == sorted(set(placement["source_ids"]))
                and all(
                    type(item) is int and not isinstance(item, bool)
                    for item in placement["source_ids"]
                ),
                "Invalid placement source IDs",
            )
            _require(
                type(placement["attempts"]) is int and placement["attempts"] >= 0,
                "Invalid rejected placement count",
            )
            _require(
                np.asarray(placement["rotation"]).shape == (3, 3)
                and np.all(np.isfinite(placement["rotation"])),
                "Invalid rotation",
            )
            _require(
                np.asarray(placement["translation"]).shape == (3,)
                and np.all(np.isfinite(placement["translation"])),
                "Invalid translation",
            )
            _require(
                set(placement["coordinates"])
                == {str(item) for item in placement["source_ids"]},
                "Placement coordinate ID coverage mismatch",
            )
        _require(
            sum(item["attempts"] for item in payload["accepted_placements"])
            == payload["rejected_attempt_count"],
            "Rejected-attempt accounting mismatch",
        )
        if units is not None:
            normalized = _validate_units(units)
            _require(
                [item[0] for item in normalized] == labels, "Packing unit order changed"
            )
            _require(
                [item[2].identity for item in normalized]
                == payload["input_graph_identities"],
                "Packing graph identity changed",
            )
            _require(
                [_system_identity(item[1]) for item in normalized]
                == payload["input_system_identities"],
                "Packing system identity changed",
            )
            _validate_placements(normalized, payload, config)
        return config


@dataclass(frozen=True)
class PeriodicPackingTransformation:
    json_text: str

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        payload = self.payload
        _require(
            payload["identity"]
            == _digest({k: v for k, v in payload.items() if k != "identity"}),
            "Packing transformation identity mismatch",
        )
        return payload["identity"]

    def validate_integrity(self, output_graph=None):
        payload = self.payload
        required = {
            "schema",
            "operation",
            "input_graph_identity",
            "output_graph_identity",
            "parameters",
            "seed",
            "provenance",
            "evidence",
            "identity",
        }
        _require(
            type(payload) is dict
            and set(payload) == required
            and payload["schema"] == GRAPH_TRANSFORMATION_SCHEMA,
            "Malformed packing transformation",
        )
        _require(
            payload["identity"]
            == _digest({k: v for k, v in payload.items() if k != "identity"}),
            "Packing transformation identity mismatch",
        )
        _require(
            payload["operation"] == "periodic_multichain_packing",
            "Unsupported packing transformation operation",
        )
        params = payload["parameters"]
        _require(
            type(params) is dict
            and set(params)
            == {
                "schema",
                "input_graph_identities",
                "output_system_identity",
                "plan_identity",
            }
            and params["schema"] == PACKING_TRANSFORMATION_SCHEMA,
            "Malformed packing transformation parameters",
        )
        _require(
            type(params["input_graph_identities"]) is list
            and params["input_graph_identities"],
            "Input graph identities required",
        )
        _text(payload["output_graph_identity"], "Output graph identity")
        _text(params["output_system_identity"], "Output system identity")
        _text(params["plan_identity"], "Packing plan identity")
        _text(payload["provenance"], "Transformation provenance")
        _evidence(payload["evidence"])
        if output_graph is not None:
            _require(
                payload["output_graph_identity"] == output_graph.identity,
                "Packing transformation output mismatch",
            )


@dataclass(frozen=True)
class PeriodicPackingResult:
    system: MolecularSystem
    graph: FinalChemicalGraph
    config: PeriodicPackingConfig
    plan: PeriodicPackingPlan
    transformation: PeriodicPackingTransformation


def _make_plan(
    units,
    config,
    placements,
    rejected,
    box_lengths,
    system,
    graph,
    remaps,
    molecule_remaps,
    graph_ids,
    system_ids,
    unit_records,
    provenance,
    evidence,
):
    payload = {
        "schema": PLAN_SCHEMA,
        "config": config.payload | {"identity": config.identity},
        "config_identity": config.identity,
        "input_graph_identities": graph_ids,
        "input_system_identities": system_ids,
        "unit_labels": [item[0] for item in units],
        "site_id_remapping": remaps,
        "molecule_membership_remapping": molecule_remaps,
        "box_lengths": list(box_lengths),
        "box_shape": BOX_SHAPE,
        "periodic": list(config.periodic),
        "placement_mode": PLACEMENT_MODE,
        "minimum_interunit_distance": config.minimum_interunit_distance,
        "seed": config.seed,
        "accepted_placements": placements,
        "rejected_attempt_count": rejected,
        "provenance": provenance,
        "evidence": list(evidence),
        "output_system_identity": _system_identity(system),
        "output_graph_identity": graph.identity,
    }
    payload["identity"] = _digest(payload)
    return PeriodicPackingPlan(pack(payload))


def _build_transformation(plan, graph):
    payload = {
        "schema": GRAPH_TRANSFORMATION_SCHEMA,
        "operation": "periodic_multichain_packing",
        "input_graph_identity": _digest(plan.payload["input_graph_identities"]),
        "output_graph_identity": graph.identity,
        "parameters": {
            "schema": PACKING_TRANSFORMATION_SCHEMA,
            "input_graph_identities": plan.payload["input_graph_identities"],
            "output_system_identity": plan.payload["output_system_identity"],
            "plan_identity": plan.identity,
        },
        "seed": plan.payload["seed"],
        "provenance": plan.payload["provenance"],
        "evidence": plan.payload["evidence"],
    }
    payload["identity"] = _digest(payload)
    return PeriodicPackingTransformation(pack(payload))


def pack_multichain_periodic(
    units,
    *,
    config,
    provenance="Deterministic force-field-neutral periodic packing",
    evidence=("caller supplied units and final graph identities",),
):
    _require(
        isinstance(config, PeriodicPackingConfig), "PeriodicPackingConfig required"
    )
    normalized = _validate_units(units)
    box_lengths = _box_for(config, normalized)
    placements, rejected = _place_units(normalized, config, box_lengths)
    system, graph, graph_ids, system_ids, remaps, molecule_remaps, unit_records = (
        _build_output(normalized, config, placements, box_lengths)
    )
    plan = _make_plan(
        normalized,
        config,
        placements,
        rejected,
        box_lengths,
        system,
        graph,
        remaps,
        molecule_remaps,
        graph_ids,
        system_ids,
        unit_records,
        provenance,
        _evidence(evidence),
    )
    system.metadata["packing"]["plan_identity"] = plan.identity
    system, graph, *_ = _build_output(
        normalized, config, placements, box_lengths, plan.identity
    )
    graph = final_graph(
        system,
        molecule_membership={
            int(i): v for i, v in graph.payload["molecule_membership"].items()
        },
    )
    transformation = _build_transformation(plan, graph)
    graph = final_graph(
        system,
        transformations=[transformation.payload],
        molecule_membership={
            int(i): v for i, v in graph.payload["molecule_membership"].items()
        },
    )
    system.metadata["final_graph_transformations"] = deepcopy(
        graph.payload["transformations"]
    )
    graph.validate_integrity(system)
    transformation.validate_integrity(graph)
    return PeriodicPackingResult(system, graph, config, plan, transformation)


def apply_periodic_packing_plan(units, plan):
    _require(isinstance(plan, PeriodicPackingPlan), "PeriodicPackingPlan required")
    normalized = _validate_units(units)
    config = plan.validate_integrity(normalized)
    box_lengths = tuple(plan.payload["box_lengths"])
    system, graph, *_ = _build_output(
        normalized,
        config,
        plan.payload["accepted_placements"],
        box_lengths,
        plan.identity,
    )
    _require(
        _system_identity(system) == plan.payload["output_system_identity"],
        "Packed system identity mismatch",
    )
    _require(
        graph.identity == plan.payload["output_graph_identity"],
        "Packed graph identity mismatch",
    )
    transformation = _build_transformation(plan, graph)
    graph = final_graph(
        system,
        transformations=[transformation.payload],
        molecule_membership={
            int(i): v for i, v in graph.payload["molecule_membership"].items()
        },
    )
    system.metadata["final_graph_transformations"] = deepcopy(
        graph.payload["transformations"]
    )
    graph.validate_integrity(system)
    transformation.validate_integrity(graph)
    return PeriodicPackingResult(system, graph, config, plan, transformation)


def save_periodic_packing_plan(plan, path):
    _require(isinstance(plan, PeriodicPackingPlan), "PeriodicPackingPlan required")
    plan.validate_integrity()
    target = Path(path)
    _require(
        not target.exists() and not target.is_symlink(),
        "Packing plan destination exists",
    )
    _require(target.parent.is_dir(), "Packing plan parent must exist")
    target.write_text(plan.json_text, encoding="utf-8")
    return target


def load_periodic_packing_plan(path, units=None):
    plan = PeriodicPackingPlan(Path(path).read_text(encoding="utf-8"))
    plan.validate_integrity(units)
    return plan
