"""Reusable integrity checks for immutable parameter-assignment results."""

from collections.abc import Mapping

from island.core import Topology
from island.exceptions import InvalidParameterAssignmentResultError, IslandError
from island.forcefields.parameters.inventory import derive_interaction_inventory
from island.forcefields.parameters.models import (
    HarmonicAngleParameter,
    HarmonicBondParameter,
    LennardJonesParameter,
    ParameterAssignmentResult,
    ParameterFamily,
    ParameterLibrary,
    ParameterSelection,
    ProperTorsionParameter,
)
from island.forcefields.parameters.signatures import (
    parameter_library_signature,
    parameter_result_content_signature,
    typing_assignment_content_signature,
)
from island.forcefields.typing.signatures import graph_signature

_FAMILIES: tuple[ParameterFamily, ...] = (
    "site",
    "bond",
    "angle",
    "proper_torsion",
)
_ARITY: dict[ParameterFamily, int] = {
    "site": 1,
    "bond": 2,
    "angle": 3,
    "proper_torsion": 4,
}
_RECORD_TYPES = (
    LennardJonesParameter,
    HarmonicBondParameter,
    HarmonicAngleParameter,
    ProperTorsionParameter,
)


def validate_parameter_assignment_result(
    result: ParameterAssignmentResult,
    topology_or_system: object,
    *,
    typing_result: object | None = None,
    library: ParameterLibrary | None = None,
) -> None:
    """Validate result content against its graph and optional authoritative inputs."""
    problems: list[str] = []
    topology = getattr(topology_or_system, "topology", topology_or_system)
    if not isinstance(topology, Topology):
        raise InvalidParameterAssignmentResultError(
            "Result integrity validation requires a Topology or MolecularSystem"
        )
    representation = getattr(topology_or_system, "representation", None)
    if representation is not None and representation != result.representation:
        problems.append("result representation does not match the system")
    try:
        current_graph_signature = graph_signature(topology)
        inventory = derive_interaction_inventory(topology)
    except (AttributeError, IslandError, TypeError, ValueError) as error:
        raise InvalidParameterAssignmentResultError(
            f"Cannot validate result against the authoritative graph: {error}"
        ) from error
    if result.graph_signature != current_graph_signature:
        problems.append("result graph signature is stale")

    required: dict[ParameterFamily, set[tuple[int, ...]]] = {
        "site": {(site_id,) for site_id in topology.sites},
        "bond": set(inventory.bonds),
        "angle": set(inventory.angles),
        "proper_torsion": set(inventory.proper_torsions),
    }
    mappings: dict[ParameterFamily, Mapping[object, ParameterSelection]] = {
        "site": result.site_assignments,
        "bond": result.bond_assignments,
        "angle": result.angle_assignments,
        "proper_torsion": result.proper_torsion_assignments,
    }

    assigned_keys: dict[ParameterFamily, set[tuple[int, ...]]] = {
        family: set() for family in _FAMILIES
    }
    typing_assignments = (
        getattr(typing_result, "assignments", {}) if typing_result is not None else {}
    )
    library_by_id = (
        {record.parameter_id: record for record in library.records}
        if library is not None
        else {}
    )
    resolved_site_types: dict[int, str] = {}
    for raw_key, selection in result.site_assignments.items():
        if (
            isinstance(raw_key, int)
            and isinstance(selection, ParameterSelection)
            and selection.site_ids == (raw_key,)
            and len(selection.atom_types) == 1
        ):
            resolved_site_types[raw_key] = selection.atom_types[0]

    for family in _FAMILIES:
        for raw_key, selection in mappings[family].items():
            try:
                key = (raw_key,) if family == "site" else tuple(raw_key)
            except TypeError:
                problems.append(f"{family} assignment has a non-iterable key")
                continue
            assigned_keys[family].add(key)
            if not isinstance(selection, ParameterSelection):
                problems.append(f"{family} assignment {key} is not a selection")
                continue
            if key != selection.site_ids:
                problems.append(
                    f"{family} mapping key {key} disagrees with selection.site_ids "
                    f"{selection.site_ids}"
                )
            if selection.family != family:
                problems.append(f"{family} assignment {key} has wrong selection family")
            if (
                len(key) != _ARITY[family]
                or len(selection.atom_types) != _ARITY[family]
            ):
                problems.append(f"{family} assignment {key} has wrong arity")
            if key not in required[family]:
                problems.append(
                    f"{family} assignment {key} is not in the authoritative inventory"
                )
            parameter = selection.parameter
            if not isinstance(parameter, _RECORD_TYPES):
                problems.append(f"{family} assignment {key} has unsupported record")
                continue
            if parameter.family != family:
                problems.append(f"{family} assignment {key} has wrong record family")
            if selection.parameter_id != parameter.parameter_id:
                problems.append(f"{family} assignment {key} has wrong parameter ID")
            if selection.source != parameter.source:
                problems.append(f"{family} assignment {key} has wrong source")
            if (parameter.library_name, parameter.library_version) != (
                result.library_name,
                result.library_version,
            ):
                problems.append(f"{family} assignment {key} has wrong record library")
            if not exact_type_pattern_matches(
                parameter.atom_types, selection.atom_types
            ):
                problems.append(f"{family} assignment {key} has wrong type pattern")
            if all(site_id in resolved_site_types for site_id in key):
                expected_result_types = tuple(
                    resolved_site_types[site_id] for site_id in key
                )
                if selection.atom_types != expected_result_types:
                    problems.append(
                        f"{family} assignment {key} types disagree with site "
                        "assignments"
                    )
            if typing_result is not None and all(
                site_id in typing_assignments for site_id in key
            ):
                expected_types = tuple(
                    typing_assignments[site_id].atom_type for site_id in key
                )
                if selection.atom_types != expected_types:
                    problems.append(
                        f"{family} assignment {key} types disagree with atom typing"
                    )
            if library is not None:
                authoritative = library_by_id.get(selection.parameter_id)
                if authoritative is None:
                    problems.append(
                        f"{family} assignment {key} parameter is absent from library"
                    )
                elif authoritative != parameter:
                    problems.append(
                        f"{family} assignment {key} record content differs from library"
                    )

    diagnostic_keys: dict[ParameterFamily, set[tuple[int, ...]]] = {
        family: set() for family in _FAMILIES
    }
    missing_counts = {family: 0 for family in _FAMILIES}
    ambiguous_counts = {family: 0 for family in _FAMILIES}
    for diagnostic in result.diagnostics:
        if diagnostic.family not in _FAMILIES:
            problems.append(f"Unknown diagnostic family {diagnostic.family!r}")
            continue
        family = diagnostic.family
        key = tuple(diagnostic.site_ids)
        if key in diagnostic_keys[family]:
            problems.append(f"Duplicate {family} diagnostic for {key}")
        diagnostic_keys[family].add(key)
        if len(key) != _ARITY[family] or len(diagnostic.atom_types) != _ARITY[family]:
            problems.append(f"{family} diagnostic {key} has wrong arity")
        if key not in required[family]:
            problems.append(
                f"{family} diagnostic {key} is not in the authoritative inventory"
            )
        if diagnostic.reason == "missing":
            missing_counts[family] += 1
        elif diagnostic.reason == "ambiguous":
            ambiguous_counts[family] += 1
        else:
            problems.append(f"{family} diagnostic {key} has invalid reason")
        if tuple(sorted(set(diagnostic.candidate_parameter_ids))) != (
            diagnostic.candidate_parameter_ids
        ):
            problems.append(f"{family} diagnostic {key} candidates are not canonical")
        if typing_result is not None and all(
            site_id in typing_assignments for site_id in key
        ):
            expected_types = tuple(
                typing_assignments[site_id].atom_type for site_id in key
            )
            if diagnostic.atom_types != expected_types:
                problems.append(
                    f"{family} diagnostic {key} types disagree with atom typing"
                )
        if all(site_id in resolved_site_types for site_id in key):
            expected_result_types = tuple(
                resolved_site_types[site_id] for site_id in key
            )
            if diagnostic.atom_types != expected_result_types:
                problems.append(
                    f"{family} diagnostic {key} types disagree with site assignments"
                )
        if library is not None:
            candidates = tuple(
                sorted(
                    record.parameter_id
                    for record in library.records_for_family(family)
                    if exact_type_pattern_matches(
                        record.atom_types, diagnostic.atom_types
                    )
                )
            )
            if candidates != diagnostic.candidate_parameter_ids:
                problems.append(
                    f"{family} diagnostic {key} candidates disagree with library"
                )
            if len(candidates) == 1:
                problems.append(
                    f"{family} diagnostic {key} has one assignable candidate"
                )
            expected_reason = "missing" if not candidates else "ambiguous"
            if diagnostic.reason != expected_reason:
                problems.append(f"{family} diagnostic {key} reason is inconsistent")

    if set(result.coverage) != set(_FAMILIES):
        problems.append("coverage families are incomplete or unsupported")
    for family in _FAMILIES:
        if assigned_keys[family] & diagnostic_keys[family]:
            problems.append(f"{family} items are both assigned and diagnosed")
        if assigned_keys[family] | diagnostic_keys[family] != required[family]:
            problems.append(f"{family} assignments/diagnostics do not cover inventory")
        expected_counts = (
            len(required[family]),
            len(assigned_keys[family]),
            missing_counts[family],
            ambiguous_counts[family],
        )
        actual_coverage = result.coverage.get(family)
        actual_counts = (
            getattr(actual_coverage, "required", None),
            getattr(actual_coverage, "assigned", None),
            getattr(actual_coverage, "missing", None),
            getattr(actual_coverage, "ambiguous", None),
        )
        if actual_counts != expected_counts:
            problems.append(f"{family} coverage disagrees with result content")

    expected_complete = all(
        assigned_keys[family] == required[family] and not diagnostic_keys[family]
        for family in _FAMILIES
    )
    if result.complete_supported_scope != expected_complete:
        problems.append("complete_supported_scope disagrees with result content")
    if result.charges_status != "unassigned":
        problems.append("Phase 4B result must report charges as unassigned")
    if result.production_validated:
        problems.append("Phase 4B result cannot claim production validation")
    if result.simulation_readiness != "not_established":
        problems.append("Phase 4B result cannot claim simulation readiness")

    if typing_result is not None:
        if result.typing_signature != getattr(typing_result, "typing_signature", None):
            problems.append("typing signature disagrees with supplied typing result")
        if result.typing_assignment_signature != typing_assignment_content_signature(
            typing_result
        ):
            problems.append("typing content signature is stale")
        if result.ruleset_signature != getattr(
            typing_result, "ruleset_signature", None
        ):
            problems.append("ruleset signature disagrees with supplied typing result")
    if library is not None:
        if (result.library_name, result.library_version) != (
            library.name,
            library.version,
        ):
            problems.append("library identity disagrees with supplied library")
        if result.library_signature != parameter_library_signature(library):
            problems.append("library signature disagrees with supplied library")

    try:
        expected_signature = parameter_result_content_signature(result)
    except (AttributeError, TypeError, ValueError) as error:
        problems.append(f"result content cannot be signed: {error}")
    else:
        if result.assignment_signature != expected_signature:
            problems.append("assignment result content signature is stale")

    if problems:
        raise InvalidParameterAssignmentResultError(
            "Invalid parameter-assignment result: " + "; ".join(dict.fromkeys(problems))
        )


def exact_type_pattern_matches(
    record_pattern: tuple[str, ...], required_pattern: tuple[str, ...]
) -> bool:
    """Return exact equality under documented complete reversal symmetry."""
    if len(record_pattern) == 1:
        return record_pattern == required_pattern
    return record_pattern == required_pattern or tuple(reversed(record_pattern)) == (
        required_pattern
    )
