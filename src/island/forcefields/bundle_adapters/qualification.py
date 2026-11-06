"""Audited, graph-specific source qualification; never runtime parameter acceptance."""

import json
from dataclasses import dataclass, field
from decimal import Decimal

from island.forcefields.assignment_core import type_graph
from island.forcefields.assignment_core.inventory import requests
from island.forcefields.assignment_core.profiles import (
    profile_contract,
    required_families,
)
from island.forcefields.assignment_core.schema import (
    CROSS_TERMS,
    Coverage,
    data,
    identity,
    require,
    unique,
)
from island.forcefields.assignment_core.typing import evidence_problem
from island.forcefields.coverage_registry import ForceFieldProfileRegistry

from .adapter import adapt

REQUIRED_ENERGY_CONVENTION = "energy_convention"
INVALID = {"unknown", "unresolved", "default", "missing"}


@dataclass(frozen=True)
class BundleQualification:
    family: str
    functional_class: str
    decision: str
    source_status: str
    bundle_identity: str | None
    profile_id: str
    source_checks: tuple[tuple[str, str], ...]
    graph_reports_json: str
    evidence_graph_json: str
    source_hashes: tuple[str, ...]
    blocking_reasons: tuple[str, ...]
    automatic_assignment_enabled: bool = field(default=False, init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)
    accepted_coefficients: int = field(default=0, init=False)
    runtime_registration: bool = field(default=False, init=False)

    @property
    def replay_identity(self):
        return identity(self)


def _explicit(value):
    return bool(value.strip()) and value.strip().lower() not in INVALID


def _source_checks(profile, trusted):
    """Inspect all source rows, including records unused by the benchmark graphs."""
    checks, reasons = {}, []
    types = {r.atom_type for r in profile.rules}
    dangling = sorted({t for p in profile.parameters for t in p.atom_types} - types)
    checks["atom_type_references"] = "complete" if types and not dangling else "blocked"
    reasons += ["undefined_atom_type:" + t for t in dangling]
    if not types:
        reasons.append("atom_type_section_missing")
    forms = {f.family: f for f in profile.forms}
    keys = {
        identity((p.family, p.atom_types, p.environment)) for p in profile.parameters
    }
    for p in profile.parameters:
        if p.form != forms.get(p.family):
            reasons.append("parameter_equation_units_or_ordering_mismatch:" + p.id)
        if not set(p.dependency_keys) <= keys:
            reasons.append("parameter_dependency_reference_missing:" + p.id)
        if {k for k, _ in p.equilibrium_references} != set(p.dependency_keys):
            reasons.append("parameter_equilibrium_reference_keys_mismatch:" + p.id)
    checks["parameter_references_and_forms"] = (
        "blocked" if any(r.startswith("parameter_") for r in reasons) else "complete"
    )
    conventions = {c.name: c for c in profile.conventions}
    c = conventions.get(REQUIRED_ENERGY_CONVENTION)
    issue = (
        "energy_convention_missing"
        if c is None or not _explicit(c.value)
        else evidence_problem(profile, c.evidence_id, trusted, _convention_claim(c))
    )
    checks["energy_convention"] = "blocked" if issue else "complete"
    if issue:
        reasons.append("energy_convention_blocked:" + issue)
    # Audit-compatible equations/units are compared exactly, without numerical conversion.
    for f in profile.forms:
        if any(not _explicit(v) for v in (f.equation, f.ordering, *f.units)):
            reasons.append("incomplete_functional_form_declaration:" + f.family)
    checks["functional_form_declarations"] = (
        "blocked"
        if any(r.startswith("incomplete_functional_form") for r in reasons)
        else "complete"
    )
    return checks, reasons


def _convention_claim(c):
    from island.forcefields.assignment_core.schema import claim_identity

    return claim_identity("convention", c)


def _coverage(graph, profile, trusted, core, source_reasons, source_present):
    typed = type_graph(graph, profile, trusted)
    rows = requests(graph, typed, profile)
    matched = dict(core.matched_records)
    dimensions = {}
    for name, families in (
        ("bond", ("bond",)),
        ("angle", ("angle",)),
        ("torsion", ("torsion",)),
        ("improper", ("improper",)),
        ("charge", ("charge",)),
        ("nonbonded", ("nonbonded",)),
        ("class_ii_cross_terms", CROSS_TERMS),
    ):
        selected = {r.id for r in rows if r.family in families}
        covered = len(selected & matched.keys())
        status = (
            "not_required_by_class_i_contract"
            if name == "class_ii_cross_terms" and profile.functional_class == "I"
            else "complete"
            if len(selected) == covered
            else "incomplete"
        )
        dimensions[name] = Coverage(len(selected), covered, status)
    dimensions["atom_typing"] = dict(core.coverage)["atom_typing"]
    conventions = {c.name: c for c in profile.conventions}
    policy_names = ("scaling_14_lj", "scaling_14_coulomb")
    supported = sum(
        name in conventions
        and _explicit(conventions[name].value)
        and evidence_problem(
            profile,
            conventions[name].evidence_id,
            trusted,
            _convention_claim(conventions[name]),
        )
        is None
        for name in policy_names
    )
    dimensions["interaction_14"] = Coverage(
        2, supported, "complete" if supported == 2 else "incomplete"
    )
    scope_ok = not graph.blockers and graph.identity in profile.graph_scope
    scope_ok = scope_ok and graph.architecture in profile.topology_support
    scope_ok = (
        scope_ok
        and "graph_scope_not_independently_audited" not in core.blocking_reasons
    )
    dimensions["topology"] = Coverage(
        1, int(scope_ok), "complete" if scope_ok else "incomplete"
    )
    reasons = {*source_reasons, *core.blocking_reasons}
    if not source_present:
        reasons.add("authorized_source_bundle_required")
    if profile.purpose == "synthetic_control":
        reasons.add("synthetic_software_control_not_scientific_evidence")
    from island.forcefields.coverage_registry import ChemicalMotifCatalog

    motifs = ChemicalMotifCatalog().detect(graph)
    if any(m.motif == "unclassified" for m in motifs):
        reasons.add("unsupported_motif_not_in_catalog")
    qualified = not reasons and core.status == "coverage_complete_not_activated"
    return {
        "graph_identity": graph.identity,
        "topology_class": graph.architecture,
        "motifs": data(motifs),
        "decision": "authorized_bundle_qualified_runtime_blocked"
        if qualified
        else (
            "parameter_bundle_incomplete"
            if source_present
            else "source_bundle_required"
        ),
        "coverage": data(dimensions),
        "sections": {
            f: {
                "source_records": sum(p.family == f for p in profile.parameters),
                "required_instances": sum(r.family == f for r in rows),
                "covered_instances": sum(
                    r.family == f and r.id in matched for r in rows
                ),
            }
            for f in required_families(profile)
        },
        "missing_terms": data(core.missing_terms),
        "blocking_reasons": sorted(reasons),
        "source_hashes": core.source_hashes,
        "core_report": data(core),
        "scope_of_qualification": "declared_exact_graph_only",
    }


def qualify(captured, graphs, *, family):
    """No IO/execution. J84 pins and an independent audit precede qualification."""
    contract = profile_contract(family)
    require(graphs, "qualification_requires_declared_graphs")
    unique([g.identity for g in graphs], "duplicate_qualification_graph")
    outcome = adapt(captured)
    m = captured["manifest"]
    require(m is None or m["family"] == family, "cross_family_qualification_forbidden")
    p = outcome.profile or contract
    trusted = outcome.trusted_evidence
    checks, reasons = _source_checks(p, trusted) if outcome.profile else ({}, [])
    reasons += list(outcome.blocking_reasons)
    registry = ForceFieldProfileRegistry((p,), trusted)
    reports = []
    for g in sorted(graphs, key=lambda g: g.identity):
        core = registry.coverage(g).profile_reports[0]
        reports.append(_coverage(g, p, trusted, core, reasons, m is not None))
    all_reasons = sorted({r for report in reports for r in report["blocking_reasons"]})
    complete = all(
        r["decision"] == "authorized_bundle_qualified_runtime_blocked" for r in reports
    )
    decision = (
        "source_bundle_required"
        if m is None
        else "authorized_bundle_qualified_runtime_blocked"
        if complete
        else "parameter_bundle_incomplete"
    )
    evidence = {
        "adapter": json.loads(outcome.evidence_graph_json),
        "adapter_status": outcome.status,
        "qualification_checks": checks,
        "source_hashes": outcome.source_hashes,
        "source_authority": "Caller-supplied independent review and external pins; hashes alone prove no authorship.",
    }
    return BundleQualification(
        family,
        contract.functional_class,
        decision,
        outcome.status,
        outcome.bundle_identity,
        p.id,
        tuple(sorted(checks.items())),
        json.dumps(data(reports), sort_keys=True, separators=(",", ":")),
        json.dumps(data(evidence), sort_keys=True, separators=(",", ":")),
        outcome.source_hashes,
        tuple(all_reasons),
    )


def rank_qualifications(results):
    """Rank independent candidates, with stable identifiers breaking coverage ties."""
    unique(
        [(r.family, r.bundle_identity) for r in results], "duplicate_candidate_bundle"
    )

    def key(result):
        reports = json.loads(result.graph_reports_json)
        required = sum(c["required"] for r in reports for c in r["coverage"].values())
        covered = sum(c["covered"] for r in reports for c in r["coverage"].values())
        return (
            result.decision != "authorized_bundle_qualified_runtime_blocked",
            -Decimal(covered) / max(1, required),
            len(result.blocking_reasons),
            result.family,
            result.profile_id,
            result.bundle_identity or "",
        )

    return tuple(sorted(results, key=key))
