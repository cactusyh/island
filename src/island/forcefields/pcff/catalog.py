"""Full pinned-source numerical inventory, independent of executable coverage.

Historical H1/H3 catalogs are not changed. Namespace and source rows are never
collapsed. Automatic forms are interpreted as records, not implicit fallbacks.
"""

from collections import Counter
from decimal import Decimal
from math import isfinite

from .class2 import FAMILIES, FIELDS, conversion
from .source import boundary, require

EXTRA = {
    "quadratic_bond": (2, ["angstrom", "E/L^2"], ["r0", "K2"]),
    "quadratic_angle": (3, ["degree", "E/A^2"], ["theta0", "K2"]),
    "torsion_1": (4, ["E", "integer", "degree"], ["Kphi", "n", "phase"]),
}


@boundary
def inspect_pcff_full_source(source):
    """Return owned complete section inventory, numerical rows and conflicts.

    An empty/unknown family remains explicit. Interpretation is not evidence of
    applicability, a typing rule, or a verified model.
    """
    source.require_assignment()
    inventory = source.inventory
    sections, rows = [], []
    for section in inventory["sections"]:
        name, namespace = section["name"], section["namespace"]
        semantic = section["records"]
        numerical = []
        spec = EXTRA.get(name)
        if name in FAMILIES:
            arity, _, units = FAMILIES[name]
            spec = arity, units, FIELDS[name]
        for line in section["lines"]:
            text = line["raw"].strip()
            if not text or text.startswith(("!", ">", "@")):
                continue
            if not spec:
                continue
            fields = text.split("!", 1)[0].split()
            arity, units, names = spec
            version = Decimal(fields[0])
            require(
                version.is_finite() and version >= 0 and fields[1].isdigit(),
                f"{name}:{line['line']}: invalid version/reference",
            )
            labels = fields[2 : 2 + arity]
            values = list(map(float, fields[2 + arity :]))
            require(
                len(labels) == arity and all(isfinite(v) for v in values),
                f"{name}:{line['line']}: invalid numeric row",
            )
            expanded = values
            if name in ("bond-angle", "end_bond-torsion_3", "angle-torsion_3") and len(
                values
            ) * 2 == len(units):
                expanded = values * 2
            require(
                len(expanded) == len(units), f"{name}:{line['line']}: coefficient count"
            )
            if name == "torsion_1":
                require(
                    expanded[1] >= 0 and expanded[1].is_integer(),
                    "Invalid torsion multiplicity",
                )
            numerical.append(
                {
                    "id": f"{name}:{namespace}:{line['line']}",
                    "line": line["line"],
                    "section": name,
                    "namespace": namespace,
                    "raw": line["raw"],
                    "version": fields[0],
                    "reference": fields[1],
                    "types": labels,
                    "original_values": values,
                    "expanded_values": expanded,
                    "fields": list(names),
                    "source_units": list(units),
                    "conversion_factors": [conversion(u) for u in units],
                    "normalized_values": [
                        v * conversion(u) for v, u in zip(expanded, units)
                    ],
                }
            )
        keys = Counter(tuple(r["types"]) for r in numerical)
        sections.append(
            {
                "name": name,
                "namespace": namespace,
                "line": section["line"],
                "numerical_records": len(numerical),
                "semantic_records": len(semantic),
                "repeated_ordered_keys": [list(k) for k, n in keys.items() if n > 1],
                "record_ids": [r["id"] for r in numerical or semantic],
                "interpretation": "numerical_records"
                if spec
                else section["interpretation"],
                "raw_lines": section["lines"],
            }
        )
        rows.extend(numerical)
    return {
        "schema": "island_pcff_full_source_catalog_v1",
        "source": source.identity,
        "declarations": inventory["declarations"],
        "sections": sections,
        "records": rows,
        "counts": {
            name: sum(
                len(s["records"]) for s in inventory["sections"] if s["name"] == name
            )
            for name in (
                "atom_types",
                "equivalence",
                "auto_equivalence",
                "bond_increments",
            )
        },
    }


AUTO_ROLES = {
    "quadratic_bond": ("bond", "bond"),
    "quadratic_angle": ("angle_end", "angle_apex", "angle_end"),
    "torsion_1": ("torsion_end", "torsion_center", "torsion_center", "torsion_end"),
    "wilson_out_of_plane": (
        "out_of_plane_end",
        "out_of_plane_center",
        "out_of_plane_end",
        "out_of_plane_end",
    ),
}


@boundary
def resolve_pcff_source_record(source, family, types, *, namespace="cff91"):
    """Resolve a named family/namespace; never substitute another potential form.

    Automatic end/apex/center mappings are positional. Wildcard ties with
    differing coefficients are unresolved, not resolved by incidental file order.
    The automatic query is explicit and does not activate charge fallback.
    """
    from itertools import permutations

    from .class2 import resolve
    from .source import records, select

    require(
        type(types) in (list, tuple) and all(type(t) is str for t in types),
        "Exact type sequence required",
    )
    supplied = list(types)
    catalog = inspect_pcff_full_source(source)
    rows = [
        r
        for r in catalog["records"]
        if r["section"] == family and r["namespace"] == namespace
    ]
    if namespace == "cff91":
        require(
            family in FAMILIES and len(supplied) == FAMILIES[family][0],
            "Unsupported ordinary family/arity",
        )
        return resolve(
            family,
            supplied,
            {r["id"]: r for r in rows},
            records(source.inventory, "equivalence"),
        )
    require(
        namespace == "cff91_auto" and family in AUTO_ROLES,
        "Unsupported family/namespace",
    )
    roles = AUTO_ROLES[family]
    require(len(supplied) == len(roles), "Incorrect automatic interaction arity")
    eqrows = records(source.inventory, "auto_equivalence")
    evidence = [select([r for r in eqrows if r["data"]["type"] == t]) for t in supplied]
    paths = [("direct", supplied, [])]
    if all(evidence):
        paths.append(
            (
                "automatic_position_equivalence",
                [
                    r["record"]["data"]["families"][role]
                    for r, role in zip(evidence, roles)
                ],
                evidence,
            )
        )
    for path, labels, eq in paths:
        queries = [
            (list(range(len(labels))), labels),
            (list(reversed(range(len(labels)))), labels[::-1]),
        ]
        if family == "wilson_out_of_plane":
            queries = [
                (
                    [p[0], 1, p[1], p[2]],
                    [labels[p[0]], labels[1], labels[p[1]], labels[p[2]]],
                )
                for p in permutations((0, 2, 3))
            ]
        for wildcards in (False, True):
            candidates = []
            for row in rows:
                for permutation, query in queries:
                    if all(
                        a == b or (wildcards and a == "*")
                        for a, b in zip(row["types"], query)
                    ):
                        candidates.append(
                            {
                                "record_id": row["id"],
                                "permutation": permutation,
                                "version": row["version"],
                                "values": row["normalized_values"],
                            }
                        )
            if not candidates:
                continue
            # Versions apply to a source key, not across unrelated wildcard patterns.
            selected = []
            for key in sorted({r["record_id"] for r in candidates}):
                row = next(r for r in rows if r["id"] == key)
                latest = max(
                    Decimal(r["version"]) for r in rows if r["types"] == row["types"]
                )
                selected.extend(
                    c
                    for c in candidates
                    if c["record_id"] == key and Decimal(c["version"]) == latest
                )
            common = {
                "family": family,
                "namespace": namespace,
                "supplied_types": supplied,
                "resolved_types": labels,
                "position_roles": list(roles),
                "path": path,
                "equivalence_evidence": eq,
                "wildcard_policy": "exact_then_consistent_wildcard_candidates_v1",
                "candidates": candidates,
            }
            if any(r["values"] != selected[0]["values"] for r in selected):
                return {
                    **common,
                    "status": "ambiguous",
                    "reason": "conflicting oriented/wildcard candidates; no order-based precedence",
                }
            return {
                **common,
                "status": "assigned",
                "selected": selected,
                "normalized_values": selected[0]["values"],
            }
    return {
        "status": "missing",
        "family": family,
        "namespace": namespace,
        "supplied_types": supplied,
        "reason": "No direct or positional automatic-equivalence row",
        "candidates": [],
    }
