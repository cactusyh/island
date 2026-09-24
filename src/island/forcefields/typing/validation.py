"""RDKit-free structural validation for precomputed atom-typing results."""

from collections.abc import Mapping

from island.core import Topology
from island.exceptions import InvalidTypingResultError
from island.forcefields.typing.models import (
    AtomTypeAssignment,
    AtomTypingResult,
    AtomTypingRuleSet,
    SiteTypingDiagnostic,
)
from island.forcefields.typing.signatures import (
    graph_signature,
    ruleset_signature,
    typing_signature,
)


def validate_complete_typing_result(
    topology: Topology, result: AtomTypingResult, ruleset: AtomTypingRuleSet
) -> None:
    """Check signatures and records without rerunning SMARTS matching."""
    if (
        not isinstance(topology, Topology)
        or not isinstance(result, AtomTypingResult)
        or not isinstance(ruleset, AtomTypingRuleSet)
    ):
        raise InvalidTypingResultError(
            "Typing validation requires topology, result, and ruleset"
        )
    problems: list[str] = []
    graph_digest = graph_signature(topology)
    rules_digest = ruleset_signature(ruleset)
    if (result.ruleset_name, result.ruleset_version) != (ruleset.name, ruleset.version):
        problems.append("typing result ruleset identity does not match")
    if result.graph_signature != graph_digest:
        problems.append("typing result graph signature is stale")
    if result.ruleset_signature != rules_digest:
        problems.append("typing result ruleset signature is stale")
    if result.typing_signature != typing_signature(
        graph_digest,
        rules_digest,
        engine_name=result.engine_name,
        engine_version=result.engine_version,
    ):
        problems.append("typing result combined signature is inconsistent")
    if not isinstance(result.assignments, Mapping) or not isinstance(
        result.diagnostics, Mapping
    ):
        raise InvalidTypingResultError(
            "Typing assignments and diagnostics must be mappings"
        )
    site_ids = set(topology.sites)
    if set(result.assignments) != site_ids:
        problems.append("typing assignment keys do not exactly cover sites")
    if set(result.diagnostics) != site_ids:
        problems.append("typing diagnostic keys do not exactly cover sites")
    if result.untyped_site_ids or result.ambiguous_site_ids:
        problems.append("typing result reports unresolved sites")
    if result.complete is not True:
        problems.append("typing result is not marked complete")

    rules = ruleset.rules_by_id
    for site_id in sorted(site_ids):
        assignment = result.assignments.get(site_id)
        diagnostic = result.diagnostics.get(site_id)
        if not isinstance(assignment, AtomTypeAssignment):
            problems.append(f"site {site_id} has an invalid typing assignment record")
            continue
        if not isinstance(diagnostic, SiteTypingDiagnostic):
            problems.append(f"site {site_id} has an invalid typing diagnostic record")
            continue
        if assignment.site_id != site_id:
            problems.append(f"assignment key/site mismatch at site {site_id}")
        if not isinstance(assignment.atom_type, str) or not assignment.atom_type:
            problems.append(f"site {site_id} assignment atom_type is invalid")
        if diagnostic.site_id != site_id:
            problems.append(f"diagnostic key/site mismatch at site {site_id}")
        if diagnostic.status != "assigned":
            problems.append(f"site {site_id} diagnostic is not assigned")
        if not isinstance(diagnostic.surviving_atom_types, tuple) or any(
            not isinstance(atom_type, str) or not atom_type
            for atom_type in diagnostic.surviving_atom_types
        ):
            problems.append(f"site {site_id} surviving atom types are invalid")
        fields = (
            (assignment.selected_rule_ids, "selected rules"),
            (assignment.matched_rule_ids, "matched rules"),
            (diagnostic.eliminated_rule_ids, "eliminated rules"),
            (diagnostic.surviving_rule_ids, "surviving rules"),
            (diagnostic.matched_rule_ids, "diagnostic matched rules"),
        )
        valid_fields = True
        for values, label in fields:
            if not isinstance(values, tuple) or any(
                not isinstance(item, str) or not item for item in values
            ):
                problems.append(f"site {site_id} {label} must be rule ID tuples")
                valid_fields = False
                continue
            if tuple(sorted(set(values))) != values:
                problems.append(f"site {site_id} {label} are not canonical")
            unknown = set(values) - set(rules)
            if unknown:
                problems.append(
                    f"site {site_id} {label} reference unknown rules {sorted(unknown)}"
                )
        if not valid_fields:
            continue
        selected = assignment.selected_rule_ids
        matched = assignment.matched_rule_ids
        if not selected or not matched:
            problems.append(f"site {site_id} has empty selected or matched rules")
        if not set(selected) <= set(matched):
            problems.append(f"site {site_id} selected rules are not matched rules")
        selected_types = {rules[item].atom_type for item in selected if item in rules}
        if selected_types != (
            {assignment.atom_type} if isinstance(assignment.atom_type, str) else set()
        ):
            problems.append(f"site {site_id} selected rules disagree with type")
        if diagnostic.matched_rule_ids != matched:
            problems.append(f"site {site_id} match diagnostics disagree")
        if diagnostic.surviving_rule_ids != selected:
            problems.append(f"site {site_id} selected diagnostics disagree")
        if diagnostic.surviving_atom_types != (assignment.atom_type,):
            problems.append(f"site {site_id} surviving type diagnostics disagree")
        eliminated = diagnostic.eliminated_rule_ids
        surviving = diagnostic.surviving_rule_ids
        if set(eliminated) & set(surviving) or set(eliminated) | set(surviving) != set(
            matched
        ):
            problems.append(f"site {site_id} diagnostic partition is inconsistent")
        eliminated_by = diagnostic.eliminated_by
        if not isinstance(eliminated_by, Mapping):
            problems.append(f"site {site_id} elimination diagnostics must be a mapping")
        elif set(eliminated_by) != set(eliminated) or any(
            not isinstance(winners, tuple)
            or not winners
            or any(not isinstance(winner, str) for winner in winners)
            or set(winners) - set(matched)
            for winners in eliminated_by.values()
        ):
            problems.append(f"site {site_id} elimination diagnostics disagree")
    if problems:
        raise InvalidTypingResultError(
            "Invalid atom-typing result: " + "; ".join(dict.fromkeys(problems))
        )
