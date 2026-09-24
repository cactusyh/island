"""Owned, signed resolved Amber assignments; no SMARTS typing is implied."""

import hashlib
import json
from copy import deepcopy
from dataclasses import asdict, dataclass, field, replace
from itertools import pairwise
from types import MappingProxyType
from typing import Any

from island.core import MolecularSystem
from island.exceptions import (
    InvalidAmberImportResultError,
    InvalidChargeAssignmentResultError,
    NonbondedPolicyError,
)
from island.forcefields.charges.models import ChargeAssignmentResult
from island.forcefields.nonbonded import NonbondedPolicy, nonbonded_policy_signature
from island.forcefields.parameterized import ParameterizedSystem
from island.forcefields.parameters.models import (
    HarmonicAngleParameter,
    HarmonicBondParameter,
    LennardJonesParameter,
    PeriodicTorsionTerm,
    ProperTorsionParameter,
)
from island.forcefields.typing.signatures import graph_signature


def _digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


@dataclass(frozen=True)
class PeriodicImproperParameter:
    """Amber periodic improper; atom 3 is central and ordering is preserved."""

    parameter_id: str
    atom_types: tuple[str, str, str, str]
    terms: tuple[PeriodicTorsionTerm, ...]
    source: str
    library_name: str
    library_version: str
    central_atom_position: int = field(default=3, init=False)
    family: str = field(default="periodic_improper", init=False)
    functional_form: str = field(default="periodic_torsion", init=False)

    def __post_init__(self) -> None:
        if len(self.atom_types) != 4 or not self.terms or any(
            not isinstance(term, PeriodicTorsionTerm) for term in self.terms
        ):
            raise InvalidAmberImportResultError("Invalid periodic improper record")
        if any(not value for value in (
            self.parameter_id, self.source, self.library_name, self.library_version,
            *self.atom_types,
        )):
            raise InvalidAmberImportResultError("Improper identity fields are required")


AmberParameter = (
    LennardJonesParameter | HarmonicBondParameter | HarmonicAngleParameter
    | ProperTorsionParameter | PeriodicImproperParameter
)


@dataclass(frozen=True)
class ImportedSelection:
    """One resolved source record mapped to stable site IDs."""

    site_ids: tuple[int, ...]
    atom_types: tuple[str, ...]
    parameter: AmberParameter
    source_reference: str
    source_type_indices: tuple[int, ...]


@dataclass(frozen=True)
class ImportedAmberResult:
    """Versioned imported source scope, with independent charge and policy results."""

    graph_signature: str
    source_sha256: str
    source: str
    parser_version: str
    adapter_version: str
    mapping: dict[int, int]
    mapping_signature: str
    atom_types: dict[int, str]
    source_type_indices: dict[int, int]
    site_assignments: dict[int, ImportedSelection]
    bond_assignments: dict[tuple[int, int], ImportedSelection]
    angle_assignments: dict[tuple[int, int, int], ImportedSelection]
    proper_torsion_assignments: dict[tuple[int, int, int, int], ImportedSelection]
    improper_assignments: dict[tuple[int, int, int, int], ImportedSelection]
    source_exclusions: tuple[tuple[int, int], ...]
    source_14_pairs: tuple[tuple[int, int], ...]
    charge_result: ChargeAssignmentResult
    nonbonded_policy: NonbondedPolicy
    provenance: dict[str, Any]
    result_signature: str
    schema: str = "island_amber_resolved_v1"

    def __post_init__(self) -> None:
        for name in (
            "mapping", "atom_types", "source_type_indices", "site_assignments",
            "bond_assignments", "angle_assignments", "proper_torsion_assignments",
            "improper_assignments",
        ):
            object.__setattr__(self, name, MappingProxyType(dict(getattr(self, name))))
        object.__setattr__(self, "provenance", MappingProxyType(deepcopy(dict(self.provenance))))
        object.__setattr__(self, "source_exclusions", tuple(self.source_exclusions))
        object.__setattr__(self, "source_14_pairs", tuple(self.source_14_pairs))

    def __deepcopy__(self, memo: dict[int, object]) -> "ImportedAmberResult":
        copied = replace(
            self,
            **{
                name: deepcopy(dict(getattr(self, name)), memo)
                for name in (
                    "mapping", "atom_types", "source_type_indices", "site_assignments",
                    "bond_assignments", "angle_assignments", "proper_torsion_assignments",
                    "improper_assignments",
                )
            },
            charge_result=deepcopy(self.charge_result, memo),
            provenance=deepcopy(dict(self.provenance), memo),
        )
        memo[id(self)] = copied
        return copied

    def _content(self) -> dict[str, Any]:
        def selections(name: str) -> list[dict[str, Any]]:
            items = getattr(self, name)
            return [
                {"key": list(key) if isinstance(key, tuple) else key,
                 "selection": asdict(value)}
                for key, value in sorted(items.items())
            ]

        return {
            "schema": self.schema, "graph": self.graph_signature,
            "source_sha256": self.source_sha256, "source": self.source,
            "parser_version": self.parser_version, "adapter_version": self.adapter_version,
            "mapping": sorted(self.mapping.items()),
            "mapping_signature": self.mapping_signature,
            "atom_types": sorted(self.atom_types.items()),
            "source_type_indices": sorted(self.source_type_indices.items()),
            "assignments": {name: selections(name) for name in (
                "site_assignments", "bond_assignments", "angle_assignments",
                "proper_torsion_assignments", "improper_assignments",
            )},
            "source_exclusions": sorted(self.source_exclusions),
            "source_14_pairs": sorted(self.source_14_pairs),
            "charge_signature": self.charge_result.result_signature,
            "policy_signature": nonbonded_policy_signature(self.nonbonded_policy),
            "provenance": dict(self.provenance),
        }

    def content_signature(self) -> str:
        return _digest(self._content())

    def validate_integrity(self, system: MolecularSystem) -> None:
        """Reject stale or malformed imported contents before snapshot creation."""
        try:
            if self.schema != "island_amber_resolved_v1":
                raise ValueError("unsupported result schema")
            if system.representation != "atomistic":
                raise ValueError("system must be atomistic")
            if graph_signature(system.topology) != self.graph_signature:
                raise ValueError("authoritative chemical graph changed")
            ids = set(system.topology.sites)
            if (set(self.mapping.values()) != ids or len(self.mapping) != len(ids)
                or set(self.mapping) != set(range(len(ids)))
                or set(self.atom_types) != ids or set(self.source_type_indices) != ids):
                raise ValueError("source mapping or atom-type coverage changed")
            if self.mapping_signature != _digest(sorted(self.mapping.items())):
                raise ValueError("source mapping signature changed")
            expected_bonds = set(system.topology.bonds)
            if set(self.bond_assignments) != expected_bonds:
                raise ValueError("bond inventory/assignment mismatch")
            adjacency = {i: system.topology.neighbors(i) for i in ids}
            expected_angles = {
                min((a, c, b), (b, c, a))
                for c in ids for a in adjacency[c] for b in adjacency[c] if a < b
            }
            if set(self.angle_assignments) != expected_angles:
                raise ValueError("angle inventory/assignment mismatch")
            for name, size, record_type, family in (
                ("site_assignments", 1, LennardJonesParameter, "site"),
                ("bond_assignments", 2, HarmonicBondParameter, "bond"),
                ("angle_assignments", 3, HarmonicAngleParameter, "angle"),
                ("proper_torsion_assignments", 4, ProperTorsionParameter, "proper_torsion"),
                ("improper_assignments", 4, PeriodicImproperParameter, "periodic_improper"),
            ):
                for key, selection in getattr(self, name).items():
                    sites = (key,) if size == 1 else key
                    if (not isinstance(selection, ImportedSelection)
                        or selection.site_ids != sites
                        or len(sites) != size or len(set(sites)) != size
                        or any(i not in ids for i in sites)
                        or not isinstance(selection.parameter, record_type)
                        or selection.parameter.family != family
                        or selection.atom_types != tuple(self.atom_types[i] for i in sites)
                        or selection.parameter.atom_types != selection.atom_types
                        or selection.source_type_indices != tuple(
                            self.source_type_indices[i] for i in sites
                        )
                        or not isinstance(selection.source_reference, str)
                        or not selection.source_reference
                        or selection.parameter.source != self.source
                        or selection.parameter.library_name != "amber_prmtop_resolved"
                        or selection.parameter.library_version != "1"):
                        raise ValueError(f"{name} selection {key} is inconsistent")
                    if size == 2 and frozenset(sites) not in map(
                        frozenset, expected_bonds
                    ):
                        raise ValueError(f"bond {key} is not authoritative")
                    if size == 3 and not (
                        sites[0] in adjacency[sites[1]]
                        and sites[2] in adjacency[sites[1]]
                    ):
                        raise ValueError(f"angle {key} is not connected")
                    if size == 4 and family == "proper_torsion" and not all(
                        b in adjacency[a] for a, b in pairwise(sites)
                    ):
                        raise ValueError(f"proper torsion {key} is not connected")
                    if size == 4 and family == "periodic_improper" and not all(
                        i in adjacency[sites[2]] for i in (sites[0], sites[1], sites[3])
                    ):
                        raise ValueError(f"improper {key} has invalid central atom")
            if set(self.site_assignments) != ids:
                raise ValueError("site LJ coverage changed")
            expected_propers = {
                min((a, b, c, d), (d, c, b, a))
                for b, c in expected_bonds
                for a in adjacency[b] - {c}
                for d in adjacency[c] - {b}
                if len({a, b, c, d}) == 4
            }
            if set(self.proper_torsion_assignments) != expected_propers:
                raise ValueError("proper torsion inventory/assignment mismatch")
            if system.topology.impropers and set(self.improper_assignments) != {
                (i.site1, i.site2, i.site3, i.site4)
                for i in system.topology.impropers
            }:
                raise ValueError("authoritative improper inventory mismatch")
            for pair in (*self.source_exclusions, *self.source_14_pairs):
                if len(pair) != 2 or pair[0] >= pair[1] or any(i not in ids for i in pair):
                    raise ValueError(f"invalid source pair {pair}")
            from island.forcefields.amber.importer import _local_pair_distances

            distances = _local_pair_distances(system.topology)
            if set(self.source_exclusions) != set(distances):
                raise ValueError("source exclusion inventory changed")
            if set(self.source_14_pairs) != {
                pair for pair, distance in distances.items() if distance == 3
            }:
                raise ValueError("source 1-4 inventory changed")
            self.charge_result.validate_integrity(system)
            if self.charge_result.graph_signature != self.graph_signature or not self.charge_result.complete:
                raise ValueError("source charges are incomplete or incompatible")
            if self.result_signature != self.content_signature():
                raise ValueError("imported result content signature changed")
        except (
            TypeError, AttributeError, KeyError, ValueError,
            InvalidChargeAssignmentResultError, NonbondedPolicyError,
        ) as error:
            raise InvalidAmberImportResultError(str(error)) from error

    def is_compatible_with(self, system: MolecularSystem) -> bool:
        try:
            self.validate_integrity(system)
        except InvalidAmberImportResultError:
            return False
        return True

    def to_parameterized_system(self, system: MolecularSystem) -> ParameterizedSystem:
        """Create an owned snapshot; not a claim of MD readiness."""
        self.validate_integrity(system)
        return ParameterizedSystem(
            system=deepcopy(system), backend_name="amber_resolved_import",
            site_assignments=deepcopy(dict(self.site_assignments)),
            interaction_assignments={name: deepcopy(dict(getattr(self, name))) for name in (
                "bond_assignments", "angle_assignments", "proper_torsion_assignments",
                "improper_assignments",
            )},
            charge_assignments=deepcopy(dict(self.charge_result.assignments)),
            nonbonded_policy=self.nonbonded_policy,
            aggregate_signature=self.result_signature,
            metadata={
                "amber_import": {"result_signature": self.result_signature,
                                 "source_sha256": self.source_sha256,
                                 "provenance": deepcopy(dict(self.provenance))},
                "aggregate": {"parameter_coverage_complete": True,
                              "charge_assignment_complete": True,
                              "nonbonded_policy_available": True,
                              "production_validated": False,
                              "simulation_readiness": "not_established"},
            },
        )
