"""Graph-based provided and atom-type partial-charge assignment engines."""

from abc import ABC, abstractmethod
from dataclasses import replace
from math import fsum, isfinite

from island.core import AtomSite, Topology
from island.exceptions import (
    IncompleteChargeAssignmentError,
    InvalidChargeDefinitionError,
    InvalidTypingResultError,
)
from island.forcefields.charges.models import (
    CHARGE_UNIT,
    AtomTypeChargeTable,
    ChargeAssignment,
    ChargeAssignmentDiagnostic,
    ChargeAssignmentResult,
    ChargeCoverage,
    ComponentChargeDiagnostic,
    validate_charge_value,
)
from island.forcefields.charges.signatures import (
    charge_input_signature,
    charge_result_signature,
    charge_table_signature,
)
from island.forcefields.charges.validation import (
    formal_charge,
    ordered_components,
)
from island.forcefields.parameters.signatures import typing_assignment_content_signature
from island.forcefields.typing import AtomTypingResult, AtomTypingRuleSet
from island.forcefields.typing.signatures import (
    graph_signature,
    ruleset_signature,
    typing_signature,
)


class ChargeAssignmentEngine(ABC):
    """Extensible contract for partial-charge methods.

    Phase 4C implementations are graph-only. Future subclasses may explicitly
    declare conformer requirements without changing the result model.
    """

    method: str
    method_version: str
    requires_coordinates: bool = False

    @abstractmethod
    def assign(
        self, topology_or_system: object, *args: object, **kwargs: object
    ) -> ChargeAssignmentResult:
        """Assign partial charges without mutating the supplied system."""


class ProvidedChargeEngine(ChargeAssignmentEngine):
    """Assign an exact caller-provided stable-site charge mapping."""

    method = "provided_stable_site_charges"
    method_version = "1"

    def assign(
        self,
        topology_or_system: object,
        charges: dict[int, float],
        *,
        source: str,
        target_charge: float | None = None,
        tolerance: float = 1.0e-6,
        strict: bool = True,
    ) -> ChargeAssignmentResult:
        topology, representation = _validate_system(topology_or_system)
        _validate_source_tolerance(source, tolerance)
        if not isinstance(charges, dict):
            raise InvalidChargeDefinitionError(
                "Provided charges must be a stable-site-ID dictionary"
            )
        if any(not isinstance(key, int) or isinstance(key, bool) for key in charges):
            raise InvalidChargeDefinitionError("Provided charge keys must be integers")
        expected_ids = set(topology.sites)
        provided_ids = set(charges)
        if provided_ids != expected_ids:
            missing = sorted(expected_ids - provided_ids)
            extra = sorted(provided_ids - expected_ids)
            raise InvalidChargeDefinitionError(
                f"Provided charges require exact site coverage; missing={missing}, "
                f"unknown={extra}"
            )
        normalized = {
            site_id: validate_charge_value(value, f"charge for site {site_id}")
            for site_id, value in charges.items()
        }
        formal_total = formal_charge(topology, tuple(sorted(topology.sites)))
        normalized_target = _validate_target(target_charge, formal_total, tolerance)
        assignments = {
            site_id: ChargeAssignment(
                site_id=site_id,
                charge=normalized[site_id],
                method=self.method,
                source=source,
            )
            for site_id in sorted(normalized)
        }
        graph_digest = graph_signature(topology)
        input_digest = charge_input_signature(
            {
                "method": self.method,
                "method_version": self.method_version,
                "graph_signature": graph_digest,
                "charges": sorted(normalized.items()),
                "source": source,
                "target_charge": normalized_target,
                "tolerance": tolerance,
            }
        )
        return _build_result(
            topology_or_system,
            assignments=assignments,
            diagnostics=[],
            method=self.method,
            method_version=self.method_version,
            source=source,
            representation=representation,
            tolerance=tolerance,
            target_charge=normalized_target,
            typing_result=None,
            table=None,
            input_signature=input_digest,
            strict=strict,
        )

    def input_signature(
        self,
        topology_or_system: object,
        charges: dict[int, float],
        *,
        source: str,
        target_charge: float | None = None,
        tolerance: float = 1.0e-6,
    ) -> str:
        topology, _ = _validate_system(topology_or_system)
        _validate_source_tolerance(source, tolerance)
        if not isinstance(charges, dict) or any(
            not isinstance(key, int) or isinstance(key, bool) for key in charges
        ):
            raise InvalidChargeDefinitionError(
                "Provided charge signature requires an integer-keyed dictionary"
            )
        if set(charges) != set(topology.sites):
            raise InvalidChargeDefinitionError(
                "Provided charge signature requires exact site coverage"
            )
        normalized = sorted(
            (site_id, validate_charge_value(value))
            for site_id, value in charges.items()
        )
        formal_total = formal_charge(topology, tuple(sorted(topology.sites)))
        target = _validate_target(target_charge, formal_total, tolerance)
        return charge_input_signature(
            {
                "method": self.method,
                "method_version": self.method_version,
                "graph_signature": graph_signature(topology),
                "charges": normalized,
                "source": source,
                "target_charge": target,
                "tolerance": tolerance,
            }
        )


class AtomTypeChargeEngine(ChargeAssignmentEngine):
    """Assign exact atom-type table charges with explicit diagnostics."""

    method = "exact_atom_type_charge_table"
    method_version = "1"

    def assign(
        self,
        topology_or_system: object,
        typing_result: AtomTypingResult,
        ruleset: AtomTypingRuleSet,
        table: AtomTypeChargeTable,
        *,
        target_charge: float | None = None,
        tolerance: float = 1.0e-6,
        strict: bool = True,
    ) -> ChargeAssignmentResult:
        topology, representation = _validate_system(topology_or_system)
        if not isinstance(table, AtomTypeChargeTable):
            raise InvalidChargeDefinitionError(
                "Atom-type charge assignment requires an AtomTypeChargeTable"
            )
        if not isinstance(typing_result, AtomTypingResult) or not isinstance(
            ruleset, AtomTypingRuleSet
        ):
            raise InvalidTypingResultError(
                "Atom-type charge assignment requires a typing result and ruleset"
            )
        _validate_source_tolerance(table.source, tolerance)
        _validate_typing(topology, representation, typing_result, ruleset, table)
        formal_total = formal_charge(topology, tuple(sorted(topology.sites)))
        normalized_target = _validate_target(target_charge, formal_total, tolerance)
        assignments: dict[int, ChargeAssignment] = {}
        diagnostics: list[ChargeAssignmentDiagnostic] = []
        for site_id in sorted(topology.sites):
            atom_type = typing_result.assignments[site_id].atom_type
            candidates = tuple(
                sorted(
                    table.entries_for_type(atom_type), key=lambda item: item.entry_id
                )
            )
            if len(candidates) != 1:
                diagnostics.append(
                    ChargeAssignmentDiagnostic(
                        reason="missing" if not candidates else "ambiguous",
                        site_ids=(site_id,),
                        atom_type=atom_type,
                        candidate_entry_ids=tuple(
                            entry.entry_id for entry in candidates
                        ),
                    )
                )
                continue
            entry = candidates[0]
            assignments[site_id] = ChargeAssignment(
                site_id=site_id,
                charge=entry.charge,
                method=self.method,
                source=entry.source,
                entry_id=entry.entry_id,
                atom_type=atom_type,
            )
        graph_digest = graph_signature(topology)
        typing_content = typing_assignment_content_signature(typing_result)
        table_digest = charge_table_signature(table)
        input_digest = charge_input_signature(
            {
                "method": self.method,
                "method_version": self.method_version,
                "graph_signature": graph_digest,
                "typing_signature": typing_result.typing_signature,
                "typing_assignment_signature": typing_content,
                "table_signature": table_digest,
                "target_charge": normalized_target,
                "tolerance": tolerance,
            }
        )
        return _build_result(
            topology_or_system,
            assignments=assignments,
            diagnostics=diagnostics,
            method=self.method,
            method_version=self.method_version,
            source=table.source,
            representation=representation,
            tolerance=tolerance,
            target_charge=normalized_target,
            typing_result=typing_result,
            table=table,
            input_signature=input_digest,
            strict=strict,
        )


def _build_result(
    topology_or_system: object,
    *,
    assignments: dict[int, ChargeAssignment],
    diagnostics: list[ChargeAssignmentDiagnostic],
    method: str,
    method_version: str,
    source: str,
    representation: str,
    tolerance: float,
    target_charge: float | None,
    typing_result: AtomTypingResult | None,
    table: AtomTypeChargeTable | None,
    input_signature: str,
    strict: bool,
) -> ChargeAssignmentResult:
    topology = getattr(topology_or_system, "topology", topology_or_system)
    components = ordered_components(topology)
    component_diagnostics = []
    for component in components:
        complete = all(site_id in assignments for site_id in component)
        observed = (
            fsum(assignments[site_id].charge for site_id in component)
            if complete
            else None
        )
        expected = formal_charge(topology, component)
        residual = observed - expected if observed is not None else None
        within = residual is not None and abs(residual) <= tolerance
        component_diagnostics.append(
            ComponentChargeDiagnostic(
                component,
                observed,
                expected,
                residual,
                within,
                complete,
            )
        )
        if complete and not within:
            diagnostics.append(
                ChargeAssignmentDiagnostic(
                    reason="charge_mismatch",
                    site_ids=component,
                    expected_charge=expected,
                    observed_charge=observed,
                    residual=residual,
                )
            )
    missing = sum(item.reason == "missing" for item in diagnostics)
    ambiguous = sum(item.reason == "ambiguous" for item in diagnostics)
    coverage = ChargeCoverage(len(topology.sites), len(assignments), missing, ambiguous)
    total_assigned = fsum(item.charge for item in assignments.values())
    total_formal = formal_charge(topology, tuple(sorted(topology.sites)))
    total_residual = (
        total_assigned - total_formal
        if len(assignments) == len(topology.sites)
        else None
    )
    total_within = total_residual is not None and abs(total_residual) <= tolerance
    complete = coverage.complete and all(
        item.within_tolerance for item in component_diagnostics
    ) and total_within
    result = ChargeAssignmentResult(
        method=method,
        method_version=method_version,
        source=source,
        representation=representation,
        assignments=assignments,
        diagnostics=tuple(diagnostics),
        coverage=coverage,
        component_diagnostics=tuple(component_diagnostics),
        complete=complete,
        unit=CHARGE_UNIT,
        tolerance=tolerance,
        total_assigned_charge=total_assigned,
        total_formal_charge=total_formal,
        total_charge_residual=total_residual,
        total_within_tolerance=total_within,
        target_charge=target_charge,
        graph_signature=graph_signature(topology),
        typing_signature=(
            typing_result.typing_signature if typing_result is not None else None
        ),
        typing_assignment_signature=(
            typing_assignment_content_signature(typing_result)
            if typing_result is not None
            else None
        ),
        table_signature=charge_table_signature(table) if table is not None else None,
        input_signature=input_signature,
        result_signature="",
        metadata={
            "requires_coordinates": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "signature_schema": "island_charge_assignment_result_v1",
        },
    )
    result = replace(result, result_signature=charge_result_signature(result))
    result.validate_integrity(
        topology_or_system, typing_result=typing_result, table=table
    )
    if strict and not result.complete:
        raise IncompleteChargeAssignmentError(
            "Charge assignment is incomplete or inconsistent with component formal "
            "charges",
            result=result,
        )
    return result


def _validate_system(topology_or_system: object) -> tuple[Topology, str]:
    topology = getattr(topology_or_system, "topology", topology_or_system)
    if not isinstance(topology, Topology):
        raise InvalidChargeDefinitionError(
            "Charge assignment requires a Topology or MolecularSystem"
        )
    representation = getattr(topology_or_system, "representation", "atomistic")
    if representation != "atomistic":
        raise InvalidChargeDefinitionError(
            "Phase 4C charge assignment supports atomistic systems only"
        )
    topology.validate_bond_graph()
    if any(not isinstance(site, AtomSite) for site in topology.sites.values()):
        raise InvalidChargeDefinitionError(
            "Phase 4C charge assignment supports AtomSite objects only"
        )
    return topology, representation


def _validate_source_tolerance(source: str, tolerance: float) -> None:
    if not isinstance(source, str) or not source.strip():
        raise InvalidChargeDefinitionError("Charge source must be a non-empty string")
    if (
        not isinstance(tolerance, (int, float))
        or isinstance(tolerance, bool)
        or not isfinite(tolerance)
        or tolerance < 0
    ):
        raise InvalidChargeDefinitionError(
            "Charge tolerance must be a non-negative finite number"
        )


def _validate_target(
    target_charge: float | None, formal_total: float, tolerance: float
) -> float | None:
    if target_charge is None:
        return None
    normalized = validate_charge_value(target_charge, "target_charge")
    if abs(normalized - formal_total) > tolerance:
        raise InvalidChargeDefinitionError(
            f"Explicit target charge {normalized} disagrees with authoritative "
            f"formal charge {formal_total}"
        )
    return normalized


def _validate_typing(
    topology: Topology,
    representation: str,
    result: AtomTypingResult,
    ruleset: AtomTypingRuleSet,
    table: AtomTypeChargeTable,
) -> None:
    problems = []
    if representation != table.supported_representation:
        problems.append("charge table representation does not match system")
    current_ruleset = ruleset_signature(ruleset)
    if (
        table.required_ruleset_name != ruleset.name
        or table.required_ruleset_version != ruleset.version
        or table.required_ruleset_signature != current_ruleset
    ):
        problems.append("charge table typing-ruleset requirement does not match")
    if not result.is_compatible_with(topology, ruleset):
        problems.append("atom typing result is stale")
    expected_typing_signature = typing_signature(
        graph_signature(topology),
        current_ruleset,
        engine_name=result.engine_name,
        engine_version=result.engine_version,
    )
    if result.typing_signature != expected_typing_signature:
        problems.append("atom typing combined signature is inconsistent")
    if (result.ruleset_name, result.ruleset_version) != (
        ruleset.name,
        ruleset.version,
    ):
        problems.append("atom typing ruleset identity does not match")
    site_ids = set(topology.sites)
    if set(result.assignments) != site_ids or set(result.diagnostics) != site_ids:
        problems.append("atom typing result does not exactly cover sites")
    if not result.complete or result.untyped_site_ids or result.ambiguous_site_ids:
        problems.append("atom typing result is incomplete")
    for site_id in sorted(site_ids & set(result.assignments)):
        assignment = result.assignments[site_id]
        diagnostic = result.diagnostics.get(site_id)
        if assignment.site_id != site_id:
            problems.append(f"atom typing assignment key mismatch at {site_id}")
        if not assignment.selected_rule_ids or not set(
            assignment.selected_rule_ids
        ) <= set(assignment.matched_rule_ids):
            problems.append(f"atom typing selected rules are invalid at {site_id}")
        selected_types = {
            ruleset.rules_by_id[rule_id].atom_type
            for rule_id in assignment.selected_rule_ids
            if rule_id in ruleset.rules_by_id
        }
        if selected_types != {assignment.atom_type}:
            problems.append(f"atom typing selected rules disagree at {site_id}")
        if diagnostic is None or diagnostic.status != "assigned":
            problems.append(f"atom typing diagnostic is unresolved at {site_id}")
        elif (
            diagnostic.surviving_atom_types != (assignment.atom_type,)
            or diagnostic.surviving_rule_ids != assignment.selected_rule_ids
            or diagnostic.matched_rule_ids != assignment.matched_rule_ids
        ):
            problems.append(f"atom typing diagnostic disagrees at {site_id}")
    if problems:
        raise InvalidTypingResultError(
            "Invalid atom typing for charge assignment: "
            + "; ".join(dict.fromkeys(problems))
        )
