"""RDKit-independent partial-charge definitions and assignment results."""

from copy import deepcopy
from dataclasses import dataclass, field, replace
from math import isfinite
from types import MappingProxyType
from typing import Any, Literal

from island.exceptions import InvalidChargeDefinitionError

CHARGE_UNIT = "elementary_charge"
ChargeDiagnosticReason = Literal["missing", "ambiguous", "charge_mismatch"]


def validate_charge_value(value: object, label: str = "charge") -> float:
    """Return a finite real charge while rejecting booleans."""
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not isfinite(value)
    ):
        raise InvalidChargeDefinitionError(f"{label} must be a finite real number")
    return float(value)


@dataclass(frozen=True)
class AtomTypeChargeEntry:
    """One versioned synthetic atom-type partial-charge lookup entry."""

    entry_id: str
    atom_type: str
    charge: float
    source: str
    table_name: str
    table_version: str
    unit: str = CHARGE_UNIT

    def __post_init__(self) -> None:
        for name in (
            "entry_id",
            "atom_type",
            "source",
            "table_name",
            "table_version",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InvalidChargeDefinitionError(
                    f"Charge entry {name} must be a non-empty string"
                )
        object.__setattr__(self, "charge", validate_charge_value(self.charge))
        if self.unit != CHARGE_UNIT:
            raise InvalidChargeDefinitionError(f"Charge unit must be {CHARGE_UNIT!r}")


@dataclass(frozen=True)
class AtomTypeChargeTable:
    """Versioned exact atom-type charge entries for one typing ruleset."""

    name: str
    version: str
    required_ruleset_name: str
    required_ruleset_version: str
    required_ruleset_signature: str
    entries: tuple[AtomTypeChargeEntry, ...]
    source: str
    supported_representation: str = "atomistic"
    unit: str = CHARGE_UNIT

    def __post_init__(self) -> None:
        for name in (
            "name",
            "version",
            "required_ruleset_name",
            "required_ruleset_version",
            "required_ruleset_signature",
            "source",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InvalidChargeDefinitionError(
                    f"Charge table {name} must be a non-empty string"
                )
        if self.supported_representation != "atomistic":
            raise InvalidChargeDefinitionError(
                "Phase 4C charge tables support atomistic systems only"
            )
        if self.unit != CHARGE_UNIT:
            raise InvalidChargeDefinitionError(
                f"Charge table unit must be {CHARGE_UNIT!r}"
            )
        entries = tuple(self.entries)
        if not entries or any(
            not isinstance(entry, AtomTypeChargeEntry) for entry in entries
        ):
            raise InvalidChargeDefinitionError(
                "Charge table requires AtomTypeChargeEntry records"
            )
        object.__setattr__(self, "entries", entries)
        entry_ids = [entry.entry_id for entry in entries]
        if len(entry_ids) != len(set(entry_ids)):
            raise InvalidChargeDefinitionError("Charge entry IDs must be unique")
        for entry in entries:
            if (entry.table_name, entry.table_version, entry.unit) != (
                self.name,
                self.version,
                self.unit,
            ):
                raise InvalidChargeDefinitionError(
                    f"Charge entry {entry.entry_id!r} has inconsistent table identity"
                )

    def entries_for_type(self, atom_type: str) -> tuple[AtomTypeChargeEntry, ...]:
        return tuple(entry for entry in self.entries if entry.atom_type == atom_type)


@dataclass(frozen=True)
class ChargeAssignment:
    """One stable-site partial-charge assignment with provenance."""

    site_id: int
    charge: float
    method: str
    source: str
    entry_id: str | None = None
    atom_type: str | None = None
    unit: str = CHARGE_UNIT

    def __post_init__(self) -> None:
        if not isinstance(self.site_id, int) or isinstance(self.site_id, bool):
            raise InvalidChargeDefinitionError("Charge site_id must be an integer")
        object.__setattr__(self, "charge", validate_charge_value(self.charge))
        for name in ("method", "source"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise InvalidChargeDefinitionError(
                    f"Charge assignment {name} must be non-empty"
                )
        if self.unit != CHARGE_UNIT:
            raise InvalidChargeDefinitionError(
                f"Charge assignment unit must be {CHARGE_UNIT!r}"
            )


@dataclass(frozen=True)
class ChargeAssignmentDiagnostic:
    """Missing, ambiguous, or total-charge inconsistency diagnostic."""

    reason: ChargeDiagnosticReason
    site_ids: tuple[int, ...]
    atom_type: str | None = None
    candidate_entry_ids: tuple[str, ...] = ()
    expected_charge: float | None = None
    observed_charge: float | None = None
    residual: float | None = None


@dataclass(frozen=True)
class ChargeCoverage:
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
            raise InvalidChargeDefinitionError(
                "Charge coverage counts must be non-negative integers"
            )
        if self.assigned + self.missing + self.ambiguous != self.required:
            raise InvalidChargeDefinitionError(
                "Assigned, missing, and ambiguous charge counts must sum to required"
            )

    @property
    def complete(self) -> bool:
        return (
            self.required == self.assigned and not self.missing and not self.ambiguous
        )


@dataclass(frozen=True)
class ComponentChargeDiagnostic:
    """Stable component charge comparison using authoritative formal charge."""

    site_ids: tuple[int, ...]
    assigned_charge: float | None
    expected_charge: float
    residual: float | None
    within_tolerance: bool
    assignments_complete: bool


@dataclass(frozen=True)
class ChargeAssignmentResult:
    """Versioned partial charges, diagnostics, signatures, and integrity methods."""

    method: str
    method_version: str
    source: str
    representation: str
    assignments: dict[int, ChargeAssignment]
    diagnostics: tuple[ChargeAssignmentDiagnostic, ...]
    coverage: ChargeCoverage
    component_diagnostics: tuple[ComponentChargeDiagnostic, ...]
    complete: bool
    unit: str
    tolerance: float
    total_assigned_charge: float
    total_formal_charge: float
    total_charge_residual: float | None
    total_within_tolerance: bool
    target_charge: float | None
    graph_signature: str
    typing_signature: str | None
    typing_assignment_signature: str | None
    table_signature: str | None
    input_signature: str
    result_signature: str
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "assignments", MappingProxyType(dict(self.assignments))
        )
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))
        object.__setattr__(
            self, "component_diagnostics", tuple(self.component_diagnostics)
        )
        object.__setattr__(
            self, "metadata", MappingProxyType(deepcopy(dict(self.metadata)))
        )

    def __deepcopy__(self, memo: dict[int, object]) -> "ChargeAssignmentResult":
        copied = replace(
            self,
            assignments=deepcopy(dict(self.assignments), memo),
            diagnostics=deepcopy(tuple(self.diagnostics), memo),
            component_diagnostics=deepcopy(tuple(self.component_diagnostics), memo),
            metadata=deepcopy(dict(self.metadata), memo),
        )
        memo[id(self)] = copied
        return copied

    def validate_integrity(
        self,
        topology_or_system: object,
        *,
        typing_result: object | None = None,
        table: AtomTypeChargeTable | None = None,
    ) -> None:
        from island.forcefields.charges.validation import validate_charge_result

        validate_charge_result(
            self, topology_or_system, typing_result=typing_result, table=table
        )

    def is_input_compatible_with(
        self,
        topology_or_system: object,
        *,
        typing_result: object | None = None,
        table: AtomTypeChargeTable | None = None,
        input_signature: str | None = None,
    ) -> bool:
        from island.forcefields.charges.validation import charge_inputs_compatible

        return charge_inputs_compatible(
            self,
            topology_or_system,
            typing_result=typing_result,
            table=table,
            input_signature=input_signature,
        )

    def is_compatible_with(
        self,
        topology_or_system: object,
        *,
        typing_result: object | None = None,
        table: AtomTypeChargeTable | None = None,
        input_signature: str | None = None,
    ) -> bool:
        from island.exceptions import InvalidChargeAssignmentResultError

        if not self.is_input_compatible_with(
            topology_or_system,
            typing_result=typing_result,
            table=table,
            input_signature=input_signature,
        ):
            return False
        try:
            self.validate_integrity(
                topology_or_system, typing_result=typing_result, table=table
            )
        except InvalidChargeAssignmentResultError:
            return False
        return True
