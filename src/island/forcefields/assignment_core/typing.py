"""Evidence-aware priority adapter around the existing SMARTS typing engine."""

from island.forcefields.assignment_core.schema import (
    TypedGraph,
    claim_identity,
    require,
    unique,
)
from island.forcefields.assignment_core.topology import as_topology


def evidence_problem(profile, evidence_id, trusted_evidence, claim=None):
    evidence = {e.id: e for e in profile.evidence}.get(evidence_id)
    if evidence is None:
        return "missing_evidence:" + evidence_id
    if evidence.synthetic:
        return (
            None
            if profile.purpose == "synthetic_control"
            else "synthetic_source_forbidden"
        )
    if evidence.status not in ("verified_stock", "verified_extension"):
        return "unqualified_evidence:" + evidence.status
    if evidence not in trusted_evidence:
        return "independent_evidence_audit_missing:" + evidence.id
    if claim is not None and claim not in evidence.audited_claims:
        return "claim_not_independently_audited:" + evidence.id
    return None


def type_graph(graph, profile, trusted_evidence=()):
    from island.forcefields.typing.models import AtomTypingRule, AtomTypingRuleSet
    from island.forcefields.typing.smarts import RDKitSmartsAtomTypingEngine

    unique([e.id for e in trusted_evidence], "duplicate_trusted_evidence_identity")
    diagnostics = []
    usable = []
    by_id = {r.id: r for r in profile.rules}
    for r in profile.rules:
        require(type(r.priority) is int, "noninteger_typing_priority")
        require(set(r.overrides) <= set(by_id), "unknown_typing_override")
        issue = evidence_problem(
            profile, r.evidence_id, trusted_evidence, claim_identity("typing_rule", r)
        )
        if issue:
            diagnostics.append("typing_rule_blocked:" + r.id + ":" + issue)
        elif "@" in r.smarts:
            # Historical SMARTS engine ignores chirality; never imply stereospecific typing.
            diagnostics.append("stereospecific_rule_requires_chiral_engine:" + r.id)
        else:
            usable.append(r)
    if not usable:
        return TypedGraph(
            graph.identity,
            profile.identity,
            (),
            tuple(sorted({*diagnostics, "no_evidenced_typing_rules"})),
            (),
        )
    ids = {r.id for r in usable}
    rules = []
    for r in usable:
        overrides = set(r.overrides) & ids
        overrides |= {o.id for o in usable if o.priority < r.priority}
        rules.append(
            AtomTypingRule(r.id, r.atom_type, r.smarts, tuple(sorted(overrides)))
        )
    # AtomTypingRuleSet rejects cycles, including priority/override contradictions.
    from island.exceptions import InvalidAtomTypingRuleError

    try:
        result = RDKitSmartsAtomTypingEngine().type_topology(
            as_topology(graph),
            AtomTypingRuleSet(profile.id, profile.identity, tuple(rules)),
            strict=False,
        )
    except InvalidAtomTypingRuleError as error:
        raise ValueError("typing_contract_invalid:" + str(error)) from error
    diagnostics.extend("untyped_site:" + str(i) for i in result.untyped_site_ids)
    diagnostics.extend("ambiguous_site:" + str(i) for i in result.ambiguous_site_ids)
    return TypedGraph(
        graph.identity,
        profile.identity,
        tuple((i, a.atom_type) for i, a in sorted(result.assignments.items())),
        tuple(sorted(diagnostics)),
        tuple((i, d.surviving_rule_ids) for i, d in sorted(result.diagnostics.items())),
    )
