"""Opt-in source supplementation; no missing-cross-term zero convention.

This policy is deliberately distinct from historical exact/ordinary resolution.
Numbered wildcard labels are patterns only when they START with an asterisk.
The center of Wilson and angle-angle interactions never moves.
"""

from decimal import Decimal
from itertools import permutations

from .catalog import AUTO_ROLES
from .class2 import FAMILIES, reverse_values
from .source import records, require, select

POLICY = "island_pcff_positional_fallbacks_v1"
DOMAIN_POLICY = "island_pcff_positional_fallbacks_v2"
MSI_POLICY = "island_pcff_msi_guarded_source_v3"
SCHEMA = "island_pcff_source_class2_assignment_v2"
MODEL_SCHEMA = "island_pcff_source_model_v2"
LOWER = {
    "quartic_bond": "quadratic_bond",
    "quartic_angle": "quadratic_angle",
    "torsion_3": "torsion_1",
    "wilson_out_of_plane": "wilson_out_of_plane",
}
POLICY_EVIDENCE = {
    "name": POLICY,
    "precedence": "ordinary direct exact/wildcard then ordinary family equivalence exact/wildcard; only missing base families enter automatic direct then positional equivalence",
    "wildcards": "leading * including numbered patterns; most constrained pattern first; equal-specificity conflicting oriented coefficients fail",
    "versions": "highest version per ordered source key, never across distinct patterns",
    "cross_terms": "all existing Class II coupling requirements retained; lower-order substitution does not imply absent coupling",
    "equilibrium": "use selected base equilibrium irrespective of polynomial order; retain dependency row identity",
    "torsion_1": "source plus-cos equation, cis=0; no class2 phase shift, no coefficient negation",
    "wilson": "same signed mean-of-three LAMMPS Class II convention; nonzero equilibrium still unsupported",
    "charge": "zero base; direct then ordinary bond then automatic bond_increment column; no wildcard increments or normalization",
}


def orientations(family, n):
    if family == "angle-angle":
        return [(0, 1, 2, 3), (3, 1, 2, 0)]
    if family == "wilson_out_of_plane":
        return [(p[0], 1, p[1], p[2]) for p in permutations((0, 2, 3))]
    return [tuple(range(n)), tuple(reversed(range(n)))]


def lookup(
    family,
    supplied,
    namespace,
    catalog,
    eqrows,
    *,
    roles_override=None,
    policy=POLICY,
    path_override=None,
    trace=False,
):
    """Owned provenance including every candidate at the winning lookup tier."""
    guarded = policy == MSI_POLICY
    trace = trace or guarded
    roles = roles_override or (
        AUTO_ROLES[family]
        if namespace == "cff91_auto"
        else [FAMILIES[family][1]] * len(supplied)
    )
    evidence = []

    def paths_in_order():
        # The guarded policy does not consult lower-priority equivalence rows
        # when the direct lookup succeeds. Preserve historical eager validation.
        if guarded:
            yield ("direct", supplied, [])
        evidence.extend(
            select([r for r in eqrows if r["data"]["type"] == t]) for t in supplied
        )
        if not guarded:
            yield ("direct", supplied, [])
        if all(evidence):
            yield (
                (path_override or "automatic_position_equivalence")
                if namespace == "cff91_auto"
                else (path_override or "ordinary_family_equivalence"),
                [
                    r["record"]["data"]["families"][role]
                    for r, role in zip(evidence, roles)
                ],
                evidence,
            )

    rows = [
        r
        for r in catalog.values()
        if r["section"] == family and r["namespace"] == namespace
    ]
    searches = []
    for path, labels, eq in paths_in_order():
        candidates = []
        for r in rows:
            for perm in orientations(family, len(labels)):
                query = [labels[k] for k in perm]
                if all(a == b or a.startswith("*") for a, b in zip(r["types"], query)):
                    values = r["normalized_values"][:]
                    if perm == tuple(reversed(range(len(labels)))):
                        values = reverse_values(family, values)
                    candidates.append(
                        {
                            "record_id": r["id"],
                            "version": r["version"],
                            "permutation": list(perm),
                            "source_types": r["types"],
                            "specificity": sum(
                                not t.startswith("*") for t in r["types"]
                            ),
                            "values": values,
                            **({"source_line": r["line"]} if guarded else {}),
                        }
                    )
        if trace:
            searches.append(
                {
                    "path": path,
                    "family": family,
                    "namespace": namespace,
                    "resolved_types": list(labels),
                    "position_roles": list(roles),
                    "equivalence_evidence": eq,
                    "legal_permutations": [
                        list(p) for p in orientations(family, len(labels))
                    ],
                    "candidate_ids": sorted({c["record_id"] for c in candidates}),
                    "records_examined": len(rows),
                    "status": "candidate_matches"
                    if candidates
                    else "no_matching_source_row",
                }
            )
        if not candidates:
            continue
        highest = {}
        for row in rows:
            key = tuple(row["types"])
            highest[key] = max(highest.get(key, Decimal(-1)), Decimal(row["version"]))
        current = [
            c
            for c in candidates
            if Decimal(c["version"]) == highest[tuple(c["source_types"])]
        ]
        specificity = max(c["specificity"] for c in current)
        # find_match tries exact rows before all wildcard rows. Numeric wildcard
        # suffixes are not priorities. Unlike its first-file-row choice, require
        # coefficient agreement at the winning tier before authorizing physics.
        eligible = (
            [c for c in current if c["specificity"] == len(labels)] or current
            if guarded
            else [c for c in current if c["specificity"] == specificity]
        )
        chosen = sorted(
            eligible,
            key=lambda c: (c["record_id"], c["permutation"]),
        )
        common = {
            **({"searches": searches} if trace else {}),
            "namespace": namespace,
            "supplied_types": list(supplied),
            "resolved_types": list(labels),
            "path": path,
            "position_roles": list(roles),
            "equivalence_evidence": eq,
            "candidates": candidates,
            "selection_policy": policy,
            "selected_family": family,
        }
        if guarded:
            common["decision_trace"] = {
                "tier": "exact" if specificity == len(labels) else "wildcard",
                "authority": "GetParameters.c:1055-1169; conflict rejection is an explicit ISLAND safeguard",
                "role_safeguard": "AA center/shared arm fixed; Wilson center fixed; no generic improper reversal",
                "candidates": [
                    {
                        "record_id": c["record_id"],
                        "permutation": c["permutation"],
                        "decision": "superseded_source_version"
                        if c not in current
                        else "lower_match_tier"
                        if c not in eligible
                        else "coefficient_agreement_required",
                    }
                    for c in candidates
                ],
            }
        if any(c["values"] != chosen[0]["values"] for c in chosen):
            return dict(
                common,
                status="ambiguous",
                reason="conflicting oriented candidates in msi2lmp match tier; file order is not physical authority"
                if guarded
                else "conflicting equally specific oriented candidates",
            )
        return dict(
            common,
            status="assigned",
            selected=chosen,
            normalized_values=chosen[0]["values"],
        )
    return {
        "status": "missing",
        **({"selection_policy": policy} if guarded else {}),
        **(
            {
                "searches": searches,
                "missing_equivalence_types": [
                    t for t, r in zip(supplied, evidence) if r is None
                ],
            }
            if trace
            else {}
        ),
        "supplied_types": list(supplied),
        "namespace": namespace,
        "selected_family": family,
        "candidates": [],
        "reason": "no source row after declared exact/equivalent/wildcard searches",
    }


def resolve(family, labels, catalog, inventory, *, policy=POLICY, trace=False):
    validate_policy(policy)
    ordinary = lookup(
        family,
        labels,
        "cff91",
        catalog,
        records(inventory, "equivalence"),
        policy=policy,
        trace=trace,
    )
    if (
        ordinary["status"] == "missing"
        and family == "nonbond(9-6)"
        and policy in (DOMAIN_POLICY, MSI_POLICY)
    ):
        fallback = lookup(
            family,
            labels,
            "cff91",
            catalog,
            records(inventory, "auto_equivalence"),
            roles_override=("nonbond",),
            policy=policy,
            trace=trace,
            path_override="auto_equivalence.nonbond",
        )
        return dict(fallback, prior_search=ordinary)
    if ordinary["status"] != "missing" or family not in LOWER:
        return ordinary
    fallback = lookup(
        LOWER[family],
        labels,
        "cff91_auto",
        catalog,
        records(inventory, "auto_equivalence"),
        policy=policy,
        trace=trace,
    )
    return dict(fallback, prior_search=ordinary)


def validate_policy(policy):
    require(
        policy in (POLICY, DOMAIN_POLICY, MSI_POLICY),
        "Unsupported PCFF resolution policy",
    )


def policy_evidence(policy):
    from copy import deepcopy

    validate_policy(policy)
    result = deepcopy(POLICY_EVIDENCE)
    if policy in (DOMAIN_POLICY, MSI_POLICY):
        result.update(
            name=policy,
            nonbonded="ordinary direct/family then automatic nonbond column into unchanged cff91 9-6 rows; no automatic cross-term equivalence invented",
        )
    if policy == MSI_POLICY:
        result.update(
            precedence="direct exact then wildcard; ordinary family exact then wildcard; only missing base terms enter separate automatic positional supplementation",
            wildcards="all highest-version rows at winning exact/wildcard tier must agree in legal physical orientation; no file-order or numeric-suffix conflict resolution",
            automatic="ISLAND extension using FRC positional columns; not implemented by pinned msi2lmp",
            cross_terms="every active coupling requires a source row and assigned equilibrium dependencies; no policy-derived BB13 zeros",
            charge="unchanged v2 zero-base oriented increment policy; not implemented by pinned msi2lmp; no correction",
            reference_revision="e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
            reference_routines="GetParameters.c find_match/match_types:1055-1169; get_equivs:1241 onward",
        )
    return result
