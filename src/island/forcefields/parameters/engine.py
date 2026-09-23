"""Deterministic exact-type numerical parameter assignment."""

from dataclasses import replace
from importlib.metadata import PackageNotFoundError, version

from island.core import Topology
from island.exceptions import (
    IncompleteParameterAssignmentError,
    InvalidTypingResultError,
    UnsupportedParameterRequirementError,
)
from island.forcefields.parameters.inventory import derive_interaction_inventory
from island.forcefields.parameters.models import (
    FamilyCoverage,
    ParameterAssignmentDiagnostic,
    ParameterAssignmentResult,
    ParameterFamily,
    ParameterLibrary,
    ParameterSelection,
)
from island.forcefields.parameters.signatures import (
    parameter_library_signature,
    parameter_result_content_signature,
    typing_assignment_content_signature,
)
from island.forcefields.parameters.validation import exact_type_pattern_matches
from island.forcefields.typing import AtomTypingResult, AtomTypingRuleSet
from island.forcefields.typing.signatures import (
    graph_signature,
    ruleset_signature,
    typing_signature,
)

_FAMILY_ORDER: tuple[ParameterFamily, ...] = (
    "site",
    "bond",
    "angle",
    "proper_torsion",
)


class ParameterAssignmentEngine:
    """Assign exact atom-type parameter patterns with reversal symmetry."""

    engine_name = "exact_type_parameter_assignment"
    engine_version = "2"

    def assign(
        self,
        topology_or_system: object,
        typing_result: AtomTypingResult,
        ruleset: AtomTypingRuleSet,
        library: ParameterLibrary,
        *,
        strict: bool = True,
    ) -> ParameterAssignmentResult:
        """Assign all supported interactions without consulting coordinates."""
        if not isinstance(typing_result, AtomTypingResult):
            raise InvalidTypingResultError("typing_result must be an AtomTypingResult")
        if not isinstance(ruleset, AtomTypingRuleSet):
            raise InvalidTypingResultError("ruleset must be an AtomTypingRuleSet")
        if not isinstance(library, ParameterLibrary):
            raise UnsupportedParameterRequirementError(
                "library must be a ParameterLibrary"
            )
        topology = getattr(topology_or_system, "topology", topology_or_system)
        if not isinstance(topology, Topology):
            raise UnsupportedParameterRequirementError(
                "Parameter assignment requires a Topology or MolecularSystem"
            )
        representation = getattr(topology_or_system, "representation", "atomistic")
        if representation != library.supported_representation:
            raise UnsupportedParameterRequirementError(
                f"Library {library.name!r} supports "
                f"{library.supported_representation!r}, not {representation!r}"
            )
        self._validate_typing(topology, typing_result, ruleset, library)
        inventory = derive_interaction_inventory(topology)
        atom_types = {
            site_id: typing_result.assignments[site_id].atom_type
            for site_id in sorted(topology.sites)
        }
        required: dict[ParameterFamily, tuple[tuple[int, ...], ...]] = {
            "site": tuple((site_id,) for site_id in sorted(topology.sites)),
            "bond": inventory.bonds,
            "angle": inventory.angles,
            "proper_torsion": inventory.proper_torsions,
        }
        assignment_maps: dict[
            ParameterFamily, dict[tuple[int, ...], ParameterSelection]
        ] = {family: {} for family in _FAMILY_ORDER}
        diagnostics: list[ParameterAssignmentDiagnostic] = []

        for family in _FAMILY_ORDER:
            records = library.records_for_family(family)
            for site_ids in required[family]:
                pattern = tuple(atom_types[site_id] for site_id in site_ids)
                candidates = tuple(
                    sorted(
                        (
                            record
                            for record in records
                            if exact_type_pattern_matches(record.atom_types, pattern)
                        ),
                        key=lambda record: record.parameter_id,
                    )
                )
                if len(candidates) != 1:
                    diagnostics.append(
                        ParameterAssignmentDiagnostic(
                            family=family,
                            site_ids=site_ids,
                            atom_types=pattern,
                            candidate_parameter_ids=tuple(
                                record.parameter_id for record in candidates
                            ),
                            reason="missing" if not candidates else "ambiguous",
                        )
                    )
                    continue
                record = candidates[0]
                assignment_maps[family][site_ids] = ParameterSelection(
                    family=family,
                    site_ids=site_ids,
                    atom_types=pattern,
                    parameter_id=record.parameter_id,
                    source=record.source,
                    parameter=record,
                )

        coverage = {}
        for family in _FAMILY_ORDER:
            family_diagnostics = [
                diagnostic for diagnostic in diagnostics if diagnostic.family == family
            ]
            coverage[family] = FamilyCoverage(
                required=len(required[family]),
                assigned=len(assignment_maps[family]),
                missing=sum(item.reason == "missing" for item in family_diagnostics),
                ambiguous=sum(
                    item.reason == "ambiguous" for item in family_diagnostics
                ),
            )
        complete = all(item.complete for item in coverage.values())
        library_digest = parameter_library_signature(library)
        typing_assignment_digest = typing_assignment_content_signature(typing_result)
        result = ParameterAssignmentResult(
            library_name=library.name,
            library_version=library.version,
            representation=representation,
            site_assignments={
                key[0]: value for key, value in assignment_maps["site"].items()
            },
            bond_assignments=assignment_maps["bond"],
            angle_assignments=assignment_maps["angle"],
            proper_torsion_assignments=assignment_maps["proper_torsion"],
            diagnostics=tuple(
                sorted(
                    diagnostics,
                    key=lambda item: (
                        _FAMILY_ORDER.index(item.family),
                        item.site_ids,
                        item.atom_types,
                    ),
                )
            ),
            coverage=coverage,
            complete_supported_scope=complete,
            charges_status="unassigned",
            production_validated=False,
            simulation_readiness="not_established",
            graph_signature=typing_result.graph_signature,
            typing_signature=typing_result.typing_signature,
            typing_assignment_signature=typing_assignment_digest,
            ruleset_signature=typing_result.ruleset_signature,
            library_signature=library_digest,
            engine_name=self.engine_name,
            engine_version=self.engine_version,
            assignment_signature="",
            metadata={
                "engine_name": self.engine_name,
                "engine_version": self.engine_version,
                "assignment_signature_schema": (
                    "island_parameter_assignment_result_v2"
                ),
                "island_version": _package_version("island"),
                "matching": "exact_atom_types_with_reversal_symmetry",
                "coordinates_consulted": False,
                "charges_assigned": False,
                "production_validated": False,
                "simulation_readiness": "not_established",
                "unsupported_families": ("improper", "class_ii_cross_terms"),
            },
        )
        result = replace(
            result,
            assignment_signature=parameter_result_content_signature(result),
        )
        result.validate_integrity(
            topology_or_system, typing_result=typing_result, library=library
        )
        if strict and not complete:
            raise IncompleteParameterAssignmentError(
                "Parameter assignment is incomplete: "
                + ", ".join(
                    f"{family} missing={coverage[family].missing} "
                    f"ambiguous={coverage[family].ambiguous}"
                    for family in _FAMILY_ORDER
                    if not coverage[family].complete
                ),
                result=result,
            )
        return result

    @staticmethod
    def _validate_typing(
        topology: Topology,
        result: AtomTypingResult,
        ruleset: AtomTypingRuleSet,
        library: ParameterLibrary,
    ) -> None:
        """Validate mutable typing-result content rather than trusting flags."""
        current_graph = graph_signature(topology)
        current_ruleset = ruleset_signature(ruleset)
        problems: list[str] = []
        if library.required_ruleset_name != ruleset.name or (
            library.required_ruleset_version != ruleset.version
        ):
            problems.append("library requires a different typing ruleset identity")
        if library.required_ruleset_signature != current_ruleset:
            problems.append("library requires different typing ruleset content")
        if (result.ruleset_name, result.ruleset_version) != (
            ruleset.name,
            ruleset.version,
        ):
            problems.append("typing result ruleset identity does not match")
        if result.graph_signature != current_graph:
            problems.append("typing result graph signature is stale")
        if result.ruleset_signature != current_ruleset:
            problems.append("typing result ruleset signature is stale")
        expected_typing_signature = typing_signature(
            current_graph,
            current_ruleset,
            engine_name=result.engine_name,
            engine_version=result.engine_version,
        )
        if result.typing_signature != expected_typing_signature:
            problems.append("typing result combined signature is inconsistent")

        site_ids = set(topology.sites)
        if set(result.assignments) != site_ids:
            problems.append("typing assignment keys do not exactly cover sites")
        if set(result.diagnostics) != site_ids:
            problems.append("typing diagnostic keys do not exactly cover sites")
        if result.untyped_site_ids or result.ambiguous_site_ids:
            problems.append("typing result reports unresolved sites")
        if not result.complete:
            problems.append("typing result is not marked complete")

        rules = ruleset.rules_by_id
        for site_id in sorted(site_ids & set(result.assignments)):
            assignment = result.assignments[site_id]
            diagnostic = result.diagnostics.get(site_id)
            if assignment.site_id != site_id:
                problems.append(f"assignment key/site mismatch at site {site_id}")
                continue
            if not assignment.selected_rule_ids or not assignment.matched_rule_ids:
                problems.append(f"site {site_id} has empty selected or matched rules")
            if not set(assignment.selected_rule_ids) <= set(
                assignment.matched_rule_ids
            ):
                problems.append(f"site {site_id} selected rules are not matched rules")
            if any(rule_id not in rules for rule_id in assignment.matched_rule_ids):
                problems.append(f"site {site_id} references unknown typing rules")
            selected_types = {
                rules[rule_id].atom_type
                for rule_id in assignment.selected_rule_ids
                if rule_id in rules
            }
            if selected_types != {assignment.atom_type}:
                problems.append(f"site {site_id} selected rules disagree with type")
            if diagnostic is None:
                continue
            if diagnostic.site_id != site_id or diagnostic.status != "assigned":
                problems.append(f"site {site_id} diagnostic is not assigned")
            if diagnostic.matched_rule_ids != assignment.matched_rule_ids:
                problems.append(f"site {site_id} match diagnostics disagree")
            if diagnostic.surviving_rule_ids != assignment.selected_rule_ids:
                problems.append(f"site {site_id} selected diagnostics disagree")
            if diagnostic.surviving_atom_types != (assignment.atom_type,):
                problems.append(f"site {site_id} surviving type diagnostics disagree")
            if set(diagnostic.eliminated_rule_ids) & set(diagnostic.surviving_rule_ids):
                problems.append(f"site {site_id} eliminated and surviving overlap")
            if set(diagnostic.eliminated_rule_ids) | set(
                diagnostic.surviving_rule_ids
            ) != set(diagnostic.matched_rule_ids):
                problems.append(f"site {site_id} diagnostic partition is inconsistent")
            if set(diagnostic.eliminated_by) != set(
                diagnostic.eliminated_rule_ids
            ) or any(
                not winners or not set(winners) <= set(diagnostic.matched_rule_ids)
                for winners in diagnostic.eliminated_by.values()
            ):
                problems.append(f"site {site_id} elimination diagnostics disagree")
            for values, label in (
                (assignment.selected_rule_ids, "selected rules"),
                (assignment.matched_rule_ids, "matched rules"),
                (diagnostic.eliminated_rule_ids, "eliminated rules"),
                (diagnostic.surviving_rule_ids, "surviving rules"),
            ):
                if tuple(sorted(set(values))) != values:
                    problems.append(f"site {site_id} {label} are not canonical")

        if problems:
            raise InvalidTypingResultError(
                "Invalid atom-typing result: " + "; ".join(dict.fromkeys(problems))
            )


def _package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unknown"
