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


def lookup(family, supplied, namespace, catalog, eqrows):
    """Owned provenance including every candidate at the winning lookup tier."""
    roles = (
        AUTO_ROLES[family]
        if namespace == "cff91_auto"
        else [FAMILIES[family][1]] * len(supplied)
    )
    evidence = [select([r for r in eqrows if r["data"]["type"] == t]) for t in supplied]
    paths = [("direct", supplied, [])]
    if all(evidence):
        paths.append(
            (
                "automatic_position_equivalence"
                if namespace == "cff91_auto"
                else "ordinary_family_equivalence",
                [
                    r["record"]["data"]["families"][role]
                    for r, role in zip(evidence, roles)
                ],
                evidence,
            )
        )
    rows = [
        r
        for r in catalog.values()
        if r["section"] == family and r["namespace"] == namespace
    ]
    for path, labels, eq in paths:
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
        chosen = sorted(
            [c for c in current if c["specificity"] == specificity],
            key=lambda c: (c["record_id"], c["permutation"]),
        )
        common = {
            "namespace": namespace,
            "supplied_types": list(supplied),
            "resolved_types": list(labels),
            "path": path,
            "position_roles": list(roles),
            "equivalence_evidence": eq,
            "candidates": candidates,
            "selection_policy": POLICY,
            "selected_family": family,
        }
        if any(c["values"] != chosen[0]["values"] for c in chosen):
            return dict(
                common,
                status="ambiguous",
                reason="conflicting equally specific oriented candidates",
            )
        return dict(
            common,
            status="assigned",
            selected=chosen,
            normalized_values=chosen[0]["values"],
        )
    return {
        "status": "missing",
        "supplied_types": list(supplied),
        "namespace": namespace,
        "selected_family": family,
        "candidates": [],
        "reason": "no source row after declared exact/equivalent/wildcard searches",
    }


def resolve(family, labels, catalog, inventory):
    ordinary = lookup(
        family, labels, "cff91", catalog, records(inventory, "equivalence")
    )
    if ordinary["status"] != "missing" or family not in LOWER:
        return ordinary
    fallback = lookup(
        LOWER[family],
        labels,
        "cff91_auto",
        catalog,
        records(inventory, "auto_equivalence"),
    )
    return dict(fallback, prior_search=ordinary)


def validate_policy(policy):
    require(policy == POLICY, "Unsupported PCFF resolution policy")
