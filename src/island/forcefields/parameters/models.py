"""Typed, RDKit-independent force-field parameter definitions and results."""

from copy import deepcopy
from dataclasses import dataclass, field
from math import isfinite
from types import MappingProxyType
from typing import Any, Literal

from island.exceptions import InvalidParameterDefinitionError

ParameterFamily = Literal["site", "bond", "angle", "proper_torsion"]
DiagnosticReason = Literal["missing", "ambiguous"]

ENERGY_UNIT = "kJ/mol"
LENGTH_UNIT = "nm"
ANGLE_UNIT = "degree"
BOND_FORCE_UNIT = "kJ/(mol*nm^2)"
ANGLE_FORCE_UNIT = "kJ/(mol*rad^2)"


def _require_text(value: object, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise InvalidParameterDefinitionError(f"{label} must be a non-empty string")


def _positive(value: object, label: str) -> None:
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not isfinite(value)
        or value <= 0
    ):
        raise InvalidParameterDefinitionError(
            f"{label} must be a positive finite number"
        )


def _validate_common(record: object) -> None:
    for name in ("parameter_id", "source", "library_name", "library_version"):
        _require_text(getattr(record, name), name)


def _validate_types(atom_types: tuple[str, ...], size: int) -> None:
    if len(atom_types) != size or any(
        not isinstance(atom_type, str) or not atom_type.strip()
        for atom_type in atom_types
    ):
        raise InvalidParameterDefinitionError(
            f"atom_types must contain exactly {size} non-empty names"
        )
    if any("*" in atom_type for atom_type in atom_types):
        raise InvalidParameterDefinitionError(
            "Phase 4B atom-type patterns are exact and cannot contain wildcards"
        )


@dataclass(frozen=True)
class LennardJonesParameter:
    """12-6 LJ parameters using ``4*epsilon*((sigma/r)^12-(sigma/r)^6)``."""

    parameter_id: str
    atom_type: str
    epsilon: float
    sigma: float
    source: str
    library_name: str
    library_version: str
    epsilon_unit: str = ENERGY_UNIT
    sigma_unit: str = LENGTH_UNIT
    family: ParameterFamily = field(default="site", init=False)
    functional_form: str = field(default="lj_12_6", init=False)

    def __post_init__(self) -> None:
        _validate_common(self)
        _validate_types((self.atom_type,), 1)
        _positive(self.epsilon, "epsilon")
        _positive(self.sigma, "sigma")
        if self.epsilon_unit != ENERGY_UNIT or self.sigma_unit != LENGTH_UNIT:
            raise InvalidParameterDefinitionError(
                f"LJ units must be {ENERGY_UNIT!r} and {LENGTH_UNIT!r}"
            )

    @property
    def atom_types(self) -> tuple[str]:
        return (self.atom_type,)


@dataclass(frozen=True)
class HarmonicBondParameter:
    """Harmonic bond parameters for ``0.5*k*(r-r0)^2``."""

    parameter_id: str
    atom_types: tuple[str, str]
    force_constant: float
    equilibrium_length: float
    source: str
    library_name: str
    library_version: str
    force_constant_unit: str = BOND_FORCE_UNIT
    length_unit: str = LENGTH_UNIT
    family: ParameterFamily = field(default="bond", init=False)
    functional_form: str = field(default="harmonic_bond", init=False)

    def __post_init__(self) -> None:
        _validate_common(self)
        object.__setattr__(self, "atom_types", tuple(self.atom_types))
        _validate_types(self.atom_types, 2)
        _positive(self.force_constant, "force_constant")
        _positive(self.equilibrium_length, "equilibrium_length")
        if (
            self.force_constant_unit != BOND_FORCE_UNIT
            or self.length_unit != LENGTH_UNIT
        ):
            raise InvalidParameterDefinitionError(
                f"Bond units must be {BOND_FORCE_UNIT!r} and {LENGTH_UNIT!r}"
            )


@dataclass(frozen=True)
class HarmonicAngleParameter:
    """Harmonic angle parameters for ``0.5*k*(theta-theta0)^2``."""

    parameter_id: str
    atom_types: tuple[str, str, str]
    force_constant: float
    equilibrium_angle: float
    source: str
    library_name: str
    library_version: str
    force_constant_unit: str = ANGLE_FORCE_UNIT
    angle_unit: str = ANGLE_UNIT
    family: ParameterFamily = field(default="angle", init=False)
    functional_form: str = field(default="harmonic_angle", init=False)

    def __post_init__(self) -> None:
        _validate_common(self)
        object.__setattr__(self, "atom_types", tuple(self.atom_types))
        _validate_types(self.atom_types, 3)
        _positive(self.force_constant, "force_constant")
        _positive(self.equilibrium_angle, "equilibrium_angle")
        if self.equilibrium_angle > 180:
            raise InvalidParameterDefinitionError(
                "equilibrium_angle must not exceed 180 degrees"
            )
        if (
            self.force_constant_unit != ANGLE_FORCE_UNIT
            or self.angle_unit != ANGLE_UNIT
        ):
            raise InvalidParameterDefinitionError(
                f"Angle units must be {ANGLE_FORCE_UNIT!r} and {ANGLE_UNIT!r}"
            )


@dataclass(frozen=True)
class PeriodicTorsionTerm:
    """One term ``k*(1+cos(n*phi-phase))`` in a proper torsion record."""

    force_constant: float
    periodicity: int
    phase: float
    force_constant_unit: str = ENERGY_UNIT
    phase_unit: str = ANGLE_UNIT

    def __post_init__(self) -> None:
        _positive(self.force_constant, "torsion force_constant")
        if (
            not isinstance(self.periodicity, int)
            or isinstance(self.periodicity, bool)
            or self.periodicity <= 0
        ):
            raise InvalidParameterDefinitionError(
                "torsion periodicity must be a positive integer"
            )
        if (
            not isinstance(self.phase, (int, float))
            or isinstance(self.phase, bool)
            or not isfinite(self.phase)
            or not 0 <= self.phase < 360
        ):
            raise InvalidParameterDefinitionError(
                "torsion phase must be finite in [0, 360) degrees"
            )
        if self.force_constant_unit != ENERGY_UNIT or self.phase_unit != ANGLE_UNIT:
            raise InvalidParameterDefinitionError(
                f"Torsion units must be {ENERGY_UNIT!r} and {ANGLE_UNIT!r}"
            )


@dataclass(frozen=True)
class ProperTorsionParameter:
    """One exact type record containing one or more periodic terms."""

    parameter_id: str
    atom_types: tuple[str, str, str, str]
    terms: tuple[PeriodicTorsionTerm, ...]
    source: str
    library_name: str
    library_version: str
    family: ParameterFamily = field(default="proper_torsion", init=False)
    functional_form: str = field(default="periodic_torsion", init=False)

    def __post_init__(self) -> None:
        _validate_common(self)
        object.__setattr__(self, "atom_types", tuple(self.atom_types))
        object.__setattr__(self, "terms", tuple(self.terms))
        _validate_types(self.atom_types, 4)
        if not self.terms or any(
            not isinstance(term, PeriodicTorsionTerm) for term in self.terms
        ):
            raise InvalidParameterDefinitionError(
                "A proper torsion record requires PeriodicTorsionTerm entries"
            )


ParameterRecord = (
    LennardJonesParameter
    | HarmonicBondParameter
    | HarmonicAngleParameter
    | ProperTorsionParameter
)


@dataclass(frozen=True)
class ParameterLibrary:
    """Versioned exact-type parameter records for one typing ruleset."""

    name: str
    version: str
    required_ruleset_name: str
    required_ruleset_version: str
    required_ruleset_signature: str
    records: tuple[ParameterRecord, ...]
    supported_representation: str = "atomistic"
    description: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        for name in (
            "name",
            "version",
            "required_ruleset_name",
            "required_ruleset_version",
            "required_ruleset_signature",
            "description",
            "source",
        ):
            _require_text(getattr(self, name), f"library {name}")
        if self.supported_representation != "atomistic":
            raise InvalidParameterDefinitionError(
                "Phase 4B libraries support atomistic representation only"
            )
        records = tuple(self.records)
        if not records:
            raise InvalidParameterDefinitionError("Parameter library cannot be empty")
        if any(
            not isinstance(
                record,
                (
                    LennardJonesParameter,
                    HarmonicBondParameter,
                    HarmonicAngleParameter,
                    ProperTorsionParameter,
                ),
            )
            for record in records
        ):
            raise InvalidParameterDefinitionError(
                "Library contains an unsupported parameter record"
            )
        object.__setattr__(self, "records", records)
        ids = [record.parameter_id for record in records]
        if len(ids) != len(set(ids)):
            duplicates = sorted({item for item in ids if ids.count(item) > 1})
            raise InvalidParameterDefinitionError(
                f"Duplicate parameter IDs: {duplicates}"
            )
        for record in records:
            if (record.library_name, record.library_version) != (
                self.name,
                self.version,
            ):
                raise InvalidParameterDefinitionError(
                    f"Parameter {record.parameter_id!r} declares library "
                    f"{record.library_name!r} {record.library_version!r}, expected "
                    f"{self.name!r} {self.version!r}"
                )

    def records_for_family(
        self, family: ParameterFamily
    ) -> tuple[ParameterRecord, ...]:
        return tuple(record for record in self.records if record.family == family)


@dataclass(frozen=True)
class ParameterSelection:
    """Selected immutable record and provenance for one required item."""

    family: ParameterFamily
    site_ids: tuple[int, ...]
    atom_types: tuple[str, ...]
    parameter_id: str
    source: str
    parameter: ParameterRecord


@dataclass(frozen=True)
class ParameterAssignmentDiagnostic:
    """A missing or conflicting exact-type parameter assignment."""

    family: ParameterFamily
    site_ids: tuple[int, ...]
    atom_types: tuple[str, ...]
    candidate_parameter_ids: tuple[str, ...]
    reason: DiagnosticReason


@dataclass(frozen=True)
class FamilyCoverage:
    """Coverage counts for one parameter family."""

    required: int
    assigned: int
    missing: int
    ambiguous: int

    def __post_init__(self) -> None:
        values = (self.required, self.assigned, self.missing, self.ambiguous)
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in values
        ):
            raise InvalidParameterDefinitionError(
                "Coverage counts must be non-negative integers"
            )
        if self.assigned + self.missing + self.ambiguous != self.required:
            raise InvalidParameterDefinitionError(
                "Assigned, missing, and ambiguous coverage must sum to required"
            )

    @property
    def complete(self) -> bool:
        return (
            self.required == self.assigned and not self.missing and not self.ambiguous
        )


@dataclass(frozen=True)
class ParameterAssignmentResult:
    """Parameter selections, coverage, provenance, and reuse signatures."""

    library_name: str
    library_version: str
    representation: str
    site_assignments: dict[int, ParameterSelection]
    bond_assignments: dict[tuple[int, int], ParameterSelection]
    angle_assignments: dict[tuple[int, int, int], ParameterSelection]
    proper_torsion_assignments: dict[tuple[int, int, int, int], ParameterSelection]
    diagnostics: tuple[ParameterAssignmentDiagnostic, ...]
    coverage: dict[ParameterFamily, FamilyCoverage]
    complete_supported_scope: bool
    charges_status: str
    production_validated: bool
    simulation_readiness: str
    graph_signature: str
    typing_signature: str
    typing_assignment_signature: str
    ruleset_signature: str
    library_signature: str
    assignment_signature: str
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "site_assignments", MappingProxyType(dict(self.site_assignments))
        )
        object.__setattr__(
            self, "bond_assignments", MappingProxyType(dict(self.bond_assignments))
        )
        object.__setattr__(
            self, "angle_assignments", MappingProxyType(dict(self.angle_assignments))
        )
        object.__setattr__(
            self,
            "proper_torsion_assignments",
            MappingProxyType(dict(self.proper_torsion_assignments)),
        )
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))
        object.__setattr__(self, "coverage", MappingProxyType(dict(self.coverage)))
        object.__setattr__(self, "metadata", MappingProxyType(deepcopy(self.metadata)))

    def is_compatible_with(
        self,
        system: object,
        typing_result: object,
        library: ParameterLibrary,
    ) -> bool:
        """Check graph, typing, and library content without inspecting coordinates."""
        from island.forcefields.parameters.signatures import (
            parameter_library_signature,
            typing_assignment_content_signature,
        )
        from island.forcefields.typing.signatures import graph_signature

        topology = getattr(system, "topology", system)
        representation = getattr(system, "representation", None)
        if representation is not None and (
            representation != library.supported_representation
            or representation != self.representation
        ):
            return False
        return (
            self.graph_signature == graph_signature(topology)
            and self.typing_signature
            == getattr(typing_result, "typing_signature", None)
            and self.typing_assignment_signature
            == typing_assignment_content_signature(typing_result)
            and self.ruleset_signature
            == getattr(typing_result, "ruleset_signature", None)
            and self.library_signature == parameter_library_signature(library)
        )

    def to_parameterized_system(self, system: object) -> object:
        """Create an owned snapshot; subsequent caller mutations are isolated."""
        from island.core import MolecularSystem
        from island.exceptions import InvalidTypingResultError
        from island.forcefields.parameterized import ParameterizedSystem
        from island.forcefields.parameters.inventory import (
            derive_interaction_inventory,
        )
        from island.forcefields.typing.signatures import graph_signature

        if not isinstance(system, MolecularSystem):
            raise TypeError("A MolecularSystem is required for parameterized output")
        if system.representation != self.representation:
            raise InvalidTypingResultError(
                "System representation does not match parameter assignments"
            )
        expected_families = {"site", "bond", "angle", "proper_torsion"}
        if (
            not self.complete_supported_scope
            or set(self.coverage) != expected_families
            or not all(item.complete for item in self.coverage.values())
            or self.diagnostics
        ):
            raise InvalidParameterDefinitionError(
                "Cannot create ParameterizedSystem from incomplete assignments"
            )
        topology = system.topology
        if graph_signature(topology) != self.graph_signature:
            raise InvalidTypingResultError(
                "System graph does not match the validated parameter assignments"
            )
        inventory = derive_interaction_inventory(topology)
        if set(self.site_assignments) != set(topology.sites) or (
            set(self.bond_assignments) != set(inventory.bonds)
            or set(self.angle_assignments) != set(inventory.angles)
            or set(self.proper_torsion_assignments) != set(inventory.proper_torsions)
        ):
            raise InvalidParameterDefinitionError(
                "Assignment keys do not cover the authoritative interaction inventory"
            )
        return ParameterizedSystem(
            system=deepcopy(system),
            backend_name=f"{self.library_name}:{self.library_version}",
            site_assignments=deepcopy(dict(self.site_assignments)),
            interaction_assignments={
                "bond": deepcopy(dict(self.bond_assignments)),
                "angle": deepcopy(dict(self.angle_assignments)),
                "proper_torsion": deepcopy(dict(self.proper_torsion_assignments)),
            },
            metadata={
                "parameter_assignment": {
                    "assignment_signature": self.assignment_signature,
                    "typing_signature": self.typing_signature,
                    "library_signature": self.library_signature,
                    "complete_supported_scope": True,
                    "charges_status": self.charges_status,
                    "production_validated": self.production_validated,
                    "simulation_readiness": self.simulation_readiness,
                }
            },
        )
