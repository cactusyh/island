"""Pure evidence/coverage resolution. No evaluator, exporter or runtime activation."""

from collections import defaultdict
from decimal import Decimal

from island.forcefields.assignment_core.inventory import requests
from island.forcefields.assignment_core.profiles import (
    required_conventions,
    required_families,
)
from island.forcefields.assignment_core.schema import (
    CROSS_TERMS,
    AssignmentDecision,
    CompletenessReport,
    Coverage,
    MissingTerm,
    claim_identity,
    identity,
    require,
    unique,
)
from island.forcefields.assignment_core.topology import as_topology, molecule
from island.forcefields.assignment_core.typing import evidence_problem, type_graph


def _record_problem(row, req, profile, trusted):
    if row.form != {f.family: f for f in profile.forms}.get(req.family):
        return "equation_units_or_ordering_mismatch"
    if row.convention_identity != profile.convention_identity:
        return "convention_identity_mismatch"
    issue = evidence_problem(
        profile, row.evidence_id, trusted, claim_identity("parameter", row)
    )
    if issue:
        return issue
    # Dependencies are parameter identities, not site aliases. Checked below.
    return None


def _validate(graph, typed, profile, trusted):
    typed.validate(graph, profile)
    require(typed == type_graph(graph, profile, trusted), "typing_snapshot_mutation")
    rows = requests(graph, typed, profile)
    by_id = {r.id: r for r in rows}
    blockers = [*graph.blockers, *typed.diagnostics]
    if profile.purpose != "synthetic_control":
        claim = claim_identity("graph_scope", graph.identity)
        if not any(
            evidence_problem(profile, e.id, trusted, claim) is None
            for e in profile.evidence
        ):
            blockers.append("graph_scope_not_independently_audited")
    if graph.identity not in profile.graph_scope:
        blockers.append("exact_graph_scope_not_evidenced")
    if graph.architecture not in profile.topology_support:
        blockers.append("topology_architecture_unsupported:" + graph.architecture)
    conventions = {c.name: c for c in profile.conventions}
    for name in required_conventions(profile):
        c = conventions.get(name)
        if (
            c is None
            or not c.value.strip()
            or c.value.lower() in ("unknown", "unresolved", "default")
        ):
            blockers.append("convention_missing:" + name)
        else:
            issue = evidence_problem(
                profile, c.evidence_id, trusted, claim_identity("convention", c)
            )
            if issue:
                blockers.append("convention_blocked:" + name + ":" + issue)
    forms = {f.family: f for f in profile.forms}
    for family in required_families(profile):
        if family not in forms:
            blockers.append("functional_form_missing:" + family)
        elif profile.purpose != "synthetic_control":
            claim = claim_identity("form", forms[family])
            if not any(
                evidence_problem(profile, e.id, trusted, claim) is None
                for e in profile.evidence
            ):
                blockers.append("functional_form_not_independently_audited:" + family)
    for family in ("bond", "angle"):
        if family in forms and not forms[family].equilibrium_fields:
            blockers.append("equilibrium_fields_missing:" + family)
    # The nonbonded equation and advertised convention must identify the same form.
    if (
        "nonbonded" in forms
        and "nonbonded_form" in conventions
        and conventions["nonbonded_form"].value != forms["nonbonded"].equation
    ):
        blockers.append("nonbonded_9_6_vs_12_6_or_equation_mismatch")
    parameters = defaultdict(list)
    for p in profile.parameters:
        parameters[(p.family, p.atom_types, p.environment)].append(p)
    matched = {}
    missing = {}
    for req in rows:
        candidates = parameters.get((req.family, req.atom_types, req.environment), [])
        reason = None
        if not all(req.atom_types):
            reason = "untyped_or_ambiguous_site"
        elif not candidates:
            reason = "exact_parameter_missing_no_fallback"
        else:
            issues = [_record_problem(p, req, profile, trusted) for p in candidates]
            if any(issues):
                reason = "candidate_blocked:" + ";".join(
                    sorted({i for i in issues if i})
                )
            elif (
                len(
                    {
                        identity(
                            (
                                p.values,
                                p.form,
                                p.applicability,
                                p.absence_reason,
                                p.dependency_keys,
                                p.equilibrium_references,
                            )
                        )
                        for p in candidates
                    }
                )
                != 1
            ):
                reason = "conflicting_exact_parameters_no_selection"
            else:
                needed = tuple(
                    sorted({by_id[d].parameter_key for d in req.dependencies})
                )
                if any(tuple(sorted(p.dependency_keys)) != needed for p in candidates):
                    reason = "class_ii_equilibrium_dependency_mismatch"
        if reason:
            missing[req.id] = MissingTerm(
                req.id,
                req.family,
                req.sites,
                reason,
                tuple(sorted(p.id for p in candidates)),
            )
        else:
            matched[req.id] = tuple(sorted(p.id for p in candidates))
    param_by_id = {p.id: p for p in profile.parameters}
    for req in rows:
        if (
            req.id in matched
            and req.dependencies
            and all(d in matched for d in req.dependencies)
        ):
            references = {}
            for d in req.dependencies:
                dep = param_by_id[matched[d][0]]
                values = dict(zip(dep.form.fields, dep.values, strict=True))
                references[by_id[d].parameter_key] = tuple(
                    (k, values[k]) for k in dep.form.equilibrium_fields
                )
            expected = tuple(sorted(references.items()))
            if any(
                param_by_id[p].equilibrium_references != expected
                for p in matched[req.id]
            ):
                missing[req.id] = MissingTerm(
                    req.id,
                    req.family,
                    req.sites,
                    "equilibrium_reference_value_mismatch",
                    matched.pop(req.id),
                )
    # All dependencies must themselves resolve. No implicit zero or generic reference.
    changed = True
    while changed:
        changed = False
        for req in rows:
            if req.id in matched and any(d not in matched for d in req.dependencies):
                missing[req.id] = MissingTerm(
                    req.id,
                    req.family,
                    req.sites,
                    "required_coupling_dependency_unresolved",
                    matched.pop(req.id),
                )
                changed = True
    charges = {}
    for req in rows:
        if req.family == "charge" and req.id in matched:
            p = param_by_id[matched[req.id][0]]
            if p.form.fields != ("q",) or p.form.units != ("e",):
                blockers.append("charge_form_requires_q_in_elementary_charge_units")
            else:
                charges[req.sites[0]] = Decimal(p.values[0])
    if len(charges) == len(graph.atom_order):
        try:
            tolerance = Decimal(conventions["charge_tolerance"].value)
            require(
                tolerance.is_finite()
                and Decimal(0) <= tolerance <= Decimal("0.000001"),
                "charge_tolerance_invalid_or_relaxed",
            )
            topology = as_topology(graph)
            for component in topology.connected_components():
                target = sum(topology.sites[i].formal_charge for i in component)
                if abs(sum(charges[i] for i in component) - target) > tolerance:
                    blockers.append(
                        "component_charge_not_conserved:" + str(min(component))
                    )
        except (KeyError, ArithmeticError, ValueError):
            blockers.append("charge_conservation_policy_unresolved")
    else:
        blockers.append("charge_conservation_not_checkable")
    # No implied generic chemical scope, including radicals/charges/stereochemistry.
    mol = molecule(graph)
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        blockers.append("radical_scope_requires_dedicated_adapter")
    if missing:
        blockers.append("required_terms_unresolved:" + str(len(missing)))
    coverage = {}
    categories = {
        "bonded_terms": {"bond", "angle", "improper"},
        "torsion": {"torsion"},
        "nonbonded": {"nonbonded"},
        "charge": {"charge"},
        "cross_terms": set(CROSS_TERMS),
    }
    for name, families in categories.items():
        required = {r.id for r in rows if r.family in families}
        covered = len(required & matched.keys())
        status = (
            "not_required_by_class_i_contract"
            if not required
            and name == "cross_terms"
            and profile.functional_class == "I"
            else "complete"
            if len(required) == covered
            else "incomplete"
        )
        coverage[name] = Coverage(len(required), covered, status)
    coverage["atom_typing"] = Coverage(
        len(graph.atom_order),
        len(typed.assignments),
        "complete" if len(typed.assignments) == len(graph.atom_order) else "incomplete",
    )
    topology_ok = not graph.blockers and graph.architecture in profile.topology_support
    coverage["topology"] = Coverage(
        1, int(topology_ok), "complete" if topology_ok else "incomplete"
    )
    blockers = sorted(set(blockers))
    status = (
        "incomplete"
        if blockers
        else "synthetic_complete"
        if profile.purpose == "synthetic_control"
        else "coverage_complete_not_activated"
    )
    # Convert nested dataclasses explicitly before hashing JSON dictionaries.
    from island.forcefields.assignment_core.schema import data

    replay = identity(
        (
            graph.identity,
            profile.identity,
            data(typed),
            [data(e) for e in trusted],
            {k: data(v) for k, v in coverage.items()},
            {k: data(v) for k, v in missing.items()},
            matched,
            blockers,
        )
    )
    return CompletenessReport(
        profile.id,
        profile.identity,
        profile.functional_class,
        graph.identity,
        status,
        tuple(sorted(coverage.items())),
        tuple(missing[k] for k in sorted(missing)),
        tuple(blockers),
        tuple(sorted({e.source_sha256 for e in profile.evidence})),
        tuple(sorted(matched.items())),
        replay,
    )


def validate_class_i(graph, typed, profile, trusted_evidence=()):
    require(profile.functional_class == "I", "class_i_validator_requires_class_i")
    return _validate(graph, typed, profile, trusted_evidence)


def validate_class_ii(graph, typed, profile, trusted_evidence=()):
    require(profile.functional_class == "II", "class_ii_validator_requires_class_ii")
    return _validate(graph, typed, profile, trusted_evidence)


def resolve(graph, profiles, trusted_evidence=()):
    """Rank coverage candidates deterministically. J82 always blocks runtime assignment."""
    unique([p.id for p in profiles], "duplicate_profile_identity")
    trusted = tuple(sorted(trusted_evidence, key=lambda e: e.id))
    reports = []
    for profile in sorted(profiles, key=lambda p: p.id):
        typed = type_graph(graph, profile, trusted)
        validator = (
            validate_class_ii if profile.functional_class == "II" else validate_class_i
        )
        reports.append(validator(graph, typed, profile, trusted))

    def rank(report):
        required = sum(v.required for _, v in report.coverage)
        covered = sum(v.covered for _, v in report.coverage)
        # Ratio before raw counts prevents Class II's larger inventory distorting ranking.
        return (
            report.status == "incomplete",
            -Decimal(covered) / max(required, 1),
            len(report.blocking_reasons),
            report.profile_id,
        )

    reports.sort(key=rank)
    blockers = ("J82_runtime_assignment_disabled",)
    if not reports:
        blockers += ("no_candidate_profiles",)
    return AssignmentDecision(
        "blocked",
        tuple(r.profile_id for r in reports),
        tuple(reports),
        blockers,
        identity([r.replay_identity for r in reports]),
    )
