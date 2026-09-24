"""Integrity and compatibility validation for charge-assignment results."""

from collections.abc import Mapping
from math import fsum, isclose, isfinite

from island.core import AtomSite, Topology
from island.exceptions import InvalidChargeAssignmentResultError, IslandError
from island.forcefields.charges.models import (
    CHARGE_UNIT,
    AtomTypeChargeTable,
    ChargeAssignment,
    ChargeAssignmentDiagnostic,
    ChargeAssignmentResult,
    ChargeCoverage,
    ComponentChargeDiagnostic,
)
from island.forcefields.charges.signatures import (
    charge_input_signature,
    charge_result_signature,
    charge_table_signature,
)
from island.forcefields.parameters.signatures import typing_assignment_content_signature
from island.forcefields.typing.signatures import (
    graph_signature,
)


def ordered_components(topology: Topology) -> tuple[tuple[int, ...], ...]:
    """Return deterministic connected components from authoritative bonds."""
    remaining = set(topology.sites)
    components = []
    while remaining:
        pending = [min(remaining)]
        component: set[int] = set()
        while pending:
            site_id = pending.pop()
            if site_id in component:
                continue
            component.add(site_id)
            pending.extend(
                sorted(topology.neighbors(site_id) - component, reverse=True)
            )
        remaining -= component
        components.append(tuple(sorted(component)))
    return tuple(components)


def formal_charge(topology: Topology, site_ids: tuple[int, ...]) -> float:
    values = []
    for site_id in site_ids:
        site = topology.sites[site_id]
        if not isinstance(site, AtomSite):
            raise InvalidChargeAssignmentResultError(
                "Charge assignment supports AtomSite objects only"
            )
        values.append(float(site.formal_charge))
    return fsum(values)


def _finite_number(value: object, *, optional: bool = False) -> bool:
    return (optional and value is None) or (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and isfinite(value)
    )


def _site_tuple(value: object) -> bool:
    return isinstance(value, tuple) and all(
        isinstance(site_id, int) and not isinstance(site_id, bool) for site_id in value
    )


def _validate_result_structure(result: ChargeAssignmentResult) -> None:
    """Reject malformed records and scalars before arithmetic or dereferencing."""
    problems: list[str] = []
    if not isinstance(result.assignments, Mapping):
        problems.append("assignments must be a mapping")
    else:
        for site_id, assignment in result.assignments.items():
            if not isinstance(site_id, int) or isinstance(site_id, bool):
                problems.append(
                    f"assignment key {site_id!r} must be an integer site ID"
                )
            if not isinstance(assignment, ChargeAssignment):
                problems.append(
                    f"site {site_id!r} has an invalid charge assignment record"
                )
                continue
            if assignment.site_id != site_id:
                problems.append(
                    f"site {site_id!r} assignment site_id disagrees with key"
                )
            if not _finite_number(assignment.charge):
                problems.append(
                    f"site {site_id!r} charge must be finite and non-boolean"
                )
            for name in ("method", "source"):
                value = getattr(assignment, name)
                if not isinstance(value, str) or not value.strip():
                    problems.append(f"site {site_id!r} assignment {name} is invalid")
            for name in ("entry_id", "atom_type"):
                value = getattr(assignment, name)
                if value is not None and (
                    not isinstance(value, str) or not value.strip()
                ):
                    problems.append(f"site {site_id!r} assignment {name} is invalid")
    if not isinstance(result.diagnostics, tuple):
        problems.append("diagnostics must be a tuple")
    else:
        for index, diagnostic in enumerate(result.diagnostics):
            if not isinstance(diagnostic, ChargeAssignmentDiagnostic):
                problems.append(f"diagnostic {index} has an invalid record type")
                continue
            if diagnostic.reason not in {"missing", "ambiguous", "charge_mismatch"}:
                problems.append(f"diagnostic {index} has an invalid reason")
            if not _site_tuple(diagnostic.site_ids):
                problems.append(f"diagnostic {index} has invalid site_ids")
            if diagnostic.atom_type is not None and (
                not isinstance(diagnostic.atom_type, str) or not diagnostic.atom_type
            ):
                problems.append(f"diagnostic {index} has invalid atom_type")
            if not isinstance(diagnostic.candidate_entry_ids, tuple) or any(
                not isinstance(entry_id, str) or not entry_id
                for entry_id in diagnostic.candidate_entry_ids
            ):
                problems.append(f"diagnostic {index} has invalid candidate_entry_ids")
            for name in ("expected_charge", "observed_charge", "residual"):
                if not _finite_number(getattr(diagnostic, name), optional=True):
                    problems.append(f"diagnostic {index} has invalid {name}")
    if not isinstance(result.coverage, ChargeCoverage):
        problems.append("coverage must be a ChargeCoverage record")
    else:
        counts = (
            result.coverage.required,
            result.coverage.assigned,
            result.coverage.missing,
            result.coverage.ambiguous,
        )
        if (
            any(
                not isinstance(value, int) or isinstance(value, bool) or value < 0
                for value in counts
            )
            or sum(counts[1:]) != counts[0]
        ):
            problems.append("coverage counts are invalid")
    if not isinstance(result.component_diagnostics, tuple):
        problems.append("component_diagnostics must be a tuple")
    else:
        for index, diagnostic in enumerate(result.component_diagnostics):
            if not isinstance(diagnostic, ComponentChargeDiagnostic):
                problems.append(
                    f"component diagnostic {index} has an invalid record type"
                )
                continue
            if not _site_tuple(diagnostic.site_ids):
                problems.append(f"component diagnostic {index} has invalid site_ids")
            if not _finite_number(diagnostic.expected_charge):
                problems.append(
                    f"component diagnostic {index} has invalid expected_charge"
                )
            for name in ("assigned_charge", "residual"):
                if not _finite_number(getattr(diagnostic, name), optional=True):
                    problems.append(f"component diagnostic {index} has invalid {name}")
            for name in ("within_tolerance", "assignments_complete"):
                if type(getattr(diagnostic, name)) is not bool:
                    problems.append(f"component diagnostic {index} has invalid {name}")
    if not _finite_number(result.tolerance) or result.tolerance < 0:
        problems.append("tolerance must be a non-negative finite number")
    if not _finite_number(result.target_charge, optional=True):
        problems.append("target_charge must be finite or None")
    for name in ("total_assigned_charge", "total_formal_charge"):
        if not _finite_number(getattr(result, name)):
            problems.append(f"{name} must be a finite number")
    if not _finite_number(result.total_charge_residual, optional=True):
        problems.append("total_charge_residual must be finite or None")
    for name in ("total_within_tolerance", "complete"):
        if type(getattr(result, name)) is not bool:
            problems.append(f"{name} must be a boolean")
    if problems:
        raise InvalidChargeAssignmentResultError(
            "Invalid charge-assignment result: " + "; ".join(dict.fromkeys(problems))
        )


def charge_inputs_compatible(
    result: ChargeAssignmentResult,
    topology_or_system: object,
    *,
    typing_result: object | None = None,
    table: AtomTypeChargeTable | None = None,
    input_signature: str | None = None,
) -> bool:
    try:
        topology = getattr(topology_or_system, "topology", topology_or_system)
        if not isinstance(topology, Topology):
            return False
        representation = getattr(topology_or_system, "representation", None)
        if representation is not None and representation != result.representation:
            return False
        if result.graph_signature != graph_signature(topology):
            return False
        if input_signature is not None and result.input_signature != input_signature:
            return False
        if result.typing_signature is None:
            if typing_result is not None:
                return False
        elif (
            typing_result is None
            or result.typing_signature
            != getattr(typing_result, "typing_signature", None)
            or result.typing_assignment_signature
            != typing_assignment_content_signature(typing_result)
        ):
            return False
        if result.table_signature is None:
            return table is None
        return table is not None and result.table_signature == charge_table_signature(
            table
        )
    except (AttributeError, IslandError, TypeError, ValueError):
        return False


def validate_charge_result(
    result: ChargeAssignmentResult,
    topology_or_system: object,
    *,
    typing_result: object | None = None,
    table: AtomTypeChargeTable | None = None,
) -> None:
    problems: list[str] = []
    topology = getattr(topology_or_system, "topology", topology_or_system)
    if not isinstance(topology, Topology):
        raise InvalidChargeAssignmentResultError(
            "Charge-result validation requires a Topology or MolecularSystem"
        )
    try:
        topology.validate_bond_graph()
        components = ordered_components(topology)
        current_graph_signature = graph_signature(topology)
    except (AttributeError, IslandError, TypeError, ValueError) as error:
        raise InvalidChargeAssignmentResultError(
            f"Cannot validate charge result against graph: {error}"
        ) from error
    _validate_result_structure(result)
    representation = getattr(topology_or_system, "representation", None)
    if representation is not None and representation != result.representation:
        problems.append("result representation does not match system")
    if result.graph_signature != current_graph_signature:
        problems.append("charge result graph signature is stale")
    if result.unit != CHARGE_UNIT:
        problems.append("charge result has unsupported unit")
    if (
        not isinstance(result.tolerance, (int, float))
        or isinstance(result.tolerance, bool)
        or not isfinite(result.tolerance)
        or result.tolerance < 0
    ):
        problems.append("charge result tolerance is invalid")

    site_ids = set(topology.sites)
    assigned_ids = set(result.assignments)
    missing_diagnostics = [
        item for item in result.diagnostics if item.reason == "missing"
    ]
    ambiguous_diagnostics = [
        item for item in result.diagnostics if item.reason == "ambiguous"
    ]
    unresolved_ids = {
        item.site_ids[0]
        for item in (*missing_diagnostics, *ambiguous_diagnostics)
        if len(item.site_ids) == 1
    }
    unresolved_seen: set[int] = set()
    typing_assignments = getattr(typing_result, "assignments", {})
    table_entries = table.entries if table is not None else ()
    for diagnostic in (*missing_diagnostics, *ambiguous_diagnostics):
        if len(diagnostic.site_ids) != 1:
            problems.append("missing/ambiguous charge diagnostic must target one site")
            continue
        site_id = diagnostic.site_ids[0]
        if site_id in unresolved_seen:
            problems.append(f"duplicate charge diagnostic for site {site_id}")
        unresolved_seen.add(site_id)
        if site_id not in site_ids:
            problems.append(f"charge diagnostic references unknown site {site_id}")
        if tuple(sorted(set(diagnostic.candidate_entry_ids))) != (
            diagnostic.candidate_entry_ids
        ):
            problems.append(f"site {site_id} charge candidates are not canonical")
        if typing_result is not None and site_id in typing_assignments:
            expected_type = typing_assignments[site_id].atom_type
            if diagnostic.atom_type != expected_type:
                problems.append(f"site {site_id} charge diagnostic type is stale")
        if table is not None and diagnostic.atom_type is not None:
            candidates = tuple(
                sorted(
                    entry.entry_id
                    for entry in table_entries
                    if entry.atom_type == diagnostic.atom_type
                )
            )
            if candidates != diagnostic.candidate_entry_ids:
                problems.append(f"site {site_id} charge candidates disagree with table")
            expected_reason = "missing" if not candidates else "ambiguous"
            if len(candidates) == 1 or diagnostic.reason != expected_reason:
                problems.append(f"site {site_id} charge diagnostic reason is invalid")
    for key, assignment in result.assignments.items():
        if not isinstance(assignment, ChargeAssignment):
            problems.append(f"site {key} charge assignment has wrong type")
            continue
        if assignment.site_id != key:
            problems.append(f"charge mapping key/site mismatch at {key}")
        if key not in site_ids:
            problems.append(f"charge assignment references unknown site {key}")
        if assignment.unit != CHARGE_UNIT:
            problems.append(f"site {key} charge has unsupported unit")
        if assignment.method != result.method:
            problems.append(f"site {key} charge method disagrees with result")
        if result.method == "provided_stable_site_charges" and (
            assignment.source != result.source
        ):
            problems.append(f"site {key} provided charge source disagrees with result")
        if not isinstance(assignment.charge, float) or not isfinite(assignment.charge):
            problems.append(f"site {key} charge is not a finite normalized float")
    if assigned_ids & unresolved_ids:
        problems.append("sites are both charge-assigned and unresolved")
    if assigned_ids | unresolved_ids != site_ids:
        problems.append("charge assignments/diagnostics do not cover exact site set")
    expected_coverage = (
        len(site_ids),
        len(assigned_ids),
        len(missing_diagnostics),
        len(ambiguous_diagnostics),
    )
    actual_coverage = (
        result.coverage.required,
        result.coverage.assigned,
        result.coverage.missing,
        result.coverage.ambiguous,
    )
    if actual_coverage != expected_coverage:
        problems.append("charge coverage disagrees with assignments/diagnostics")

    component_by_ids = {item.site_ids: item for item in result.component_diagnostics}
    if len(component_by_ids) != len(result.component_diagnostics):
        problems.append("component charge diagnostics contain duplicates")
    if set(component_by_ids) != set(components):
        problems.append("component charge diagnostics disagree with graph components")
    expected_component_complete = True
    mismatch_by_component = {
        item.site_ids: item
        for item in result.diagnostics
        if item.reason == "charge_mismatch"
    }
    if set(mismatch_by_component) - set(components):
        problems.append("charge mismatch diagnostic targets a non-component")
    for component in components:
        diagnostic = component_by_ids.get(component)
        if diagnostic is None:
            expected_component_complete = False
            continue
        expected = formal_charge(topology, component)
        component_complete = all(site_id in result.assignments for site_id in component)
        observed = (
            fsum(result.assignments[site_id].charge for site_id in component)
            if component_complete
            else None
        )
        residual = observed - expected if observed is not None else None
        within = residual is not None and abs(residual) <= result.tolerance
        if (
            diagnostic.expected_charge != expected
            or diagnostic.assignments_complete != component_complete
            or diagnostic.assigned_charge != observed
            or diagnostic.residual != residual
            or diagnostic.within_tolerance != within
        ):
            problems.append(f"component charge diagnostic {component} is inconsistent")
        mismatch = mismatch_by_component.get(component)
        if within and mismatch is not None:
            problems.append(f"component {component} has a spurious mismatch diagnostic")
        if not within and component_complete:
            if mismatch is None:
                problems.append(f"component {component} lacks a mismatch diagnostic")
            elif (
                mismatch.expected_charge != expected
                or mismatch.observed_charge != observed
                or mismatch.residual != residual
            ):
                problems.append(
                    f"component {component} mismatch diagnostic is inconsistent"
                )
        expected_component_complete &= within

    assigned_total = fsum(item.charge for item in result.assignments.values())
    formal_total = formal_charge(topology, tuple(sorted(site_ids)))
    total_residual = assigned_total - formal_total if assigned_ids == site_ids else None
    total_within = (
        total_residual is not None and abs(total_residual) <= result.tolerance
    )
    if not isclose(result.total_assigned_charge, assigned_total, abs_tol=0.0):
        problems.append("total assigned charge disagrees with site assignments")
    if result.total_formal_charge != formal_total:
        problems.append("total formal charge disagrees with chemical graph")
    if (
        result.total_charge_residual != total_residual
        or result.total_within_tolerance != total_within
    ):
        problems.append("total charge diagnostics disagree with assignments")
    if result.target_charge is not None and (
        abs(result.target_charge - formal_total) > result.tolerance
    ):
        problems.append("explicit target charge disagrees with formal charge")
    expected_complete = (
        result.coverage.complete and expected_component_complete and total_within
    )
    if result.complete != expected_complete:
        problems.append("charge completeness flag disagrees with result content")

    if typing_result is not None:
        if result.typing_signature != getattr(typing_result, "typing_signature", None):
            problems.append("charge typing signature is stale")
        if result.typing_assignment_signature != typing_assignment_content_signature(
            typing_result
        ):
            problems.append("charge typing-content signature is stale")
        for site_id, assignment in result.assignments.items():
            if site_id in typing_assignments and assignment.atom_type != (
                typing_assignments[site_id].atom_type
            ):
                problems.append(f"site {site_id} charge atom type is stale")
    if table is not None:
        if result.table_signature != charge_table_signature(table):
            problems.append("charge table signature is stale")
        table_by_id = {entry.entry_id: entry for entry in table.entries}
        for site_id, assignment in result.assignments.items():
            entry = table_by_id.get(assignment.entry_id)
            if entry is None:
                problems.append(f"site {site_id} charge entry is absent from table")
            elif (
                entry.charge != assignment.charge
                or entry.atom_type != assignment.atom_type
                or entry.source != assignment.source
            ):
                problems.append(f"site {site_id} charge entry differs from table")
    elif result.method == "provided_stable_site_charges":
        for site_id, assignment in result.assignments.items():
            if assignment.entry_id is not None or assignment.atom_type is not None:
                problems.append(
                    f"provided charge at site {site_id} has table-only provenance"
                )

    if result.method == "provided_stable_site_charges":
        expected_input = charge_input_signature(
            {
                "method": result.method,
                "method_version": result.method_version,
                "graph_signature": result.graph_signature,
                "charges": sorted(
                    (site_id, assignment.charge)
                    for site_id, assignment in result.assignments.items()
                ),
                "source": result.source,
                "target_charge": result.target_charge,
                "tolerance": result.tolerance,
            }
        )
        if result.input_signature != expected_input:
            problems.append("provided charge input signature is stale")
    elif result.method == "exact_atom_type_charge_table":
        expected_input = charge_input_signature(
            {
                "method": result.method,
                "method_version": result.method_version,
                "graph_signature": result.graph_signature,
                "typing_signature": result.typing_signature,
                "typing_assignment_signature": result.typing_assignment_signature,
                "table_signature": result.table_signature,
                "target_charge": result.target_charge,
                "tolerance": result.tolerance,
            }
        )
        if result.input_signature != expected_input:
            problems.append("atom-type charge input signature is stale")
    else:
        problems.append(f"unsupported charge method {result.method!r}")

    try:
        expected_signature = charge_result_signature(result)
    except (AttributeError, TypeError, ValueError) as error:
        problems.append(f"charge result cannot be signed: {error}")
    else:
        if result.result_signature != expected_signature:
            problems.append("charge result content signature is stale")
    if problems:
        raise InvalidChargeAssignmentResultError(
            "Invalid charge-assignment result: " + "; ".join(dict.fromkeys(problems))
        )
