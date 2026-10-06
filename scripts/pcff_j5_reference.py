"""Independent raw-FRC oracle for J5 inspection, never a runtime resolver.

Uses no production parser, equivalence selector, conversion or inventory helper.
Conflicting numbered wildcard patterns remain unresolved: no priority is invented.
"""

from collections import defaultdict
from decimal import Decimal
from itertools import combinations, permutations
from math import pi

ARITY = {
    "quartic_bond": 2,
    "quartic_angle": 3,
    "torsion_3": 4,
    "wilson_out_of_plane": 4,
    "nonbond(9-6)": 1,
    "bond-bond": 3,
    "bond-angle": 3,
    "bond-bond_1_3": 4,
    "end_bond-torsion_3": 4,
    "middle_bond-torsion_3": 4,
    "angle-torsion_3": 4,
    "angle-angle-torsion_1": 4,
    "angle-angle": 4,
    "quadratic_bond": 2,
    "quadratic_angle": 3,
    "torsion_1": 4,
    "bond_increments": 2,
}
ORDINARY = ("nonbond", "bond", "angle", "torsion", "out_of_plane")
AUTO = (
    "nonbond",
    "bond_increment",
    "bond",
    "angle_end",
    "angle_apex",
    "torsion_end",
    "torsion_center",
    "out_of_plane_end",
    "out_of_plane_center",
)
ROLE = {
    f: (
        "bond"
        if f == "quartic_bond"
        else "angle"
        if f in ("quartic_angle", "bond-bond", "bond-angle")
        else "out_of_plane"
        if f in ("angle-angle", "wilson_out_of_plane")
        else "nonbond"
        if f == "nonbond(9-6)"
        else "torsion"
    )
    for f in ARITY
}
LOWER = {
    "quartic_bond": "quadratic_bond",
    "quartic_angle": "quadratic_angle",
    "torsion_3": "torsion_1",
    "wilson_out_of_plane": "wilson_out_of_plane",
}
POSITION = {
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


def read_source(raw):
    rows = defaultdict(list)
    section = namespace = None
    for n, line in enumerate(raw.decode().splitlines(), 1):
        w = line.split("!", 1)[0].split()
        if not w:
            continue
        if w[0].startswith("#"):
            section = w[0][1:]
            namespace = " ".join(w[1:])
            continue
        if not w[0][0].isdigit():
            continue
        if section in ARITY:
            k = ARITY[section]
            rows[section, namespace].append(
                {
                    "id": f"{section}:{namespace}:{n}",
                    "line": n,
                    "version": Decimal(w[0]),
                    "types": w[2 : 2 + k],
                    "values": [Decimal(x) for x in w[2 + k :]],
                    "raw": line,
                }
            )
        elif section in ("equivalence", "auto_equivalence"):
            cols = ORDINARY if section == "equivalence" else AUTO
            rows[section, namespace].append(
                {
                    "id": f"{section}:{namespace}:{n}",
                    "line": n,
                    "version": Decimal(w[0]),
                    "type": w[2],
                    "map": dict(zip(cols, w[3:], strict=True)),
                    "raw": line,
                }
            )
    return rows


def equivalents(rows, types, kind, roles):
    available = rows[kind, "cff91" if kind == "equivalence" else "cff91_auto"]
    out = []
    ids = []
    for t, role in zip(types, roles, strict=True):
        found = [r for r in available if r["type"] == t]
        if not found:
            return None, []
        version = max(r["version"] for r in found)
        found = [r for r in found if r["version"] == version]
        assert len({tuple(sorted(r["map"].items())) for r in found}) == 1
        out.append(found[0]["map"][role])
        ids.extend(r["id"] for r in found)
    return out, ids


def oriented_values(row, family, perm):
    v = [float(x) for x in row["values"]]
    if family in ("bond-angle", "end_bond-torsion_3", "angle-torsion_3"):
        required = 2 if family == "bond-angle" else 6
        if len(v) * 2 == required:
            v *= 2
    if tuple(perm) == tuple(reversed(range(len(row["types"])))):
        if family == "bond-angle":
            v.reverse()
        elif family in ("end_bond-torsion_3", "angle-torsion_3"):
            v = v[3:] + v[:3]
    for i, x in enumerate(v):
        if family in ("quartic_bond", "quadratic_bond", "nonbond(9-6)") and i == 0:
            continue
        if (
            (family in ("quartic_angle", "quadratic_angle") and i == 0)
            or (family == "torsion_3" and i % 2 == 1)
            or (family == "torsion_1" and i == 2)
            or (family == "wilson_out_of_plane" and i == 1)
        ):
            v[i] = x * pi / 180
        elif family == "torsion_1" and i == 1:
            continue
        else:
            v[i] = x * 4.184
    return v


def match(rows, family, namespace, types):
    n = len(types)
    perms = [tuple(range(n)), tuple(reversed(range(n)))]
    if family == "angle-angle":
        perms = [(0, 1, 2, 3), (3, 1, 2, 0)]
    elif family == "wilson_out_of_plane":
        perms = [(p[0], 1, p[1], p[2]) for p in permutations((0, 2, 3))]
    candidates = []
    for row in rows[family, namespace]:
        for p in perms:
            if all(
                t.startswith("*") or t == types[i]
                for t, i in zip(row["types"], p, strict=True)
            ):
                candidates.append((row, p))
    if not candidates:
        return {"status": "missing", "candidate_ids": []}
    versions = {}
    for r, p in candidates:
        versions[tuple(r["types"])] = max(
            versions.get(tuple(r["types"]), Decimal(-1)), r["version"]
        )
    current = [
        (r, p) for r, p in candidates if r["version"] == versions[tuple(r["types"])]
    ]
    specificity = max(
        sum(not t.startswith("*") for t in r["types"]) for r, p in current
    )
    winners = [
        (r, p)
        for r, p in current
        if sum(not t.startswith("*") for t in r["types"]) == specificity
    ]
    values = [oriented_values(r, family, p) for r, p in winners]
    return {
        "status": "assigned" if all(v == values[0] for v in values) else "ambiguous",
        "candidate_ids": sorted({r["id"] for r, p in candidates}),
        "selected_ids": sorted({r["id"] for r, p in winners}),
        "values": values[0] if all(v == values[0] for v in values) else None,
        "permutations": [list(p) for r, p in winners],
    }


def resolve_raw(rows, family, types):
    searches = []
    stages = [(family, "cff91", "direct", types, [])]
    ordinary, e = equivalents(rows, types, "equivalence", [ROLE[family]] * len(types))
    if ordinary:
        stages.append((family, "cff91", "ordinary_family_equivalence", ordinary, e))
    if family in LOWER:
        f = LOWER[family]
        stages.append((f, "cff91_auto", "direct", types, []))
        automatic, e = equivalents(rows, types, "auto_equivalence", POSITION[f])
        if automatic:
            stages.append(
                (f, "cff91_auto", "automatic_position_equivalence", automatic, e)
            )
    elif family == "nonbond(9-6)":
        stages.append((family, "cff91", "direct", types, []))
        auto, e = equivalents(rows, types, "auto_equivalence", ["nonbond"])
        if auto:
            stages.append((family, "cff91", "auto_equivalence.nonbond", auto, e))
    for f, ns, path, query, eq in stages:
        result = match(rows, f, ns, query)
        searches.append(
            {
                "family": f,
                "namespace": ns,
                "path": path,
                "types": list(query),
                "equivalence_ids": eq,
                "candidate_ids": result["candidate_ids"],
            }
        )
        if result["status"] != "missing":
            return dict(result, family=f, searches=searches)
    return {
        "status": "missing",
        "family": family,
        "searches": searches,
        "candidate_ids": [],
    }


def charge_raw(rows, system, labels):
    values = {i: Decimal(0) for i in labels}
    contributions = []
    missing = []
    for bond in system.topology.bonds.values():
        sites = (bond.site1, bond.site2)
        types = [labels[i] for i in sites]
        paths = [("direct", types, [])]
        for kind, role in [
            ("equivalence", "bond"),
            ("auto_equivalence", "bond_increment"),
        ]:
            eq, ev = equivalents(rows, types, kind, [role, role])
            if eq:
                paths.append((kind + "." + role, eq, ev))
        found = False
        searches = []
        for path, query, eq in paths:
            matches = []
            for (f, ns), rs in rows.items():
                if f != "bond_increments":
                    continue
                for r in rs:
                    for p in ((0, 1), (1, 0)):
                        if r["types"] == [query[i] for i in p]:
                            matches.append((r, p))
            searches.append(
                {
                    "path": path,
                    "types": query,
                    "candidate_ids": sorted({r["id"] for r, p in matches}),
                    "equivalence_ids": eq,
                }
            )
            if not matches:
                continue
            top = max(r["version"] for r, p in matches)
            matches = [(r, p) for r, p in matches if r["version"] == top]
            increments = [
                tuple(r["values"][p.index(i)] for i in range(2)) for r, p in matches
            ]
            if len(set(increments)) != 1:
                missing.append(
                    {
                        "sites": sites,
                        "reason": "conflicting_increment_candidates",
                        "searches": searches,
                    }
                )
                found = True
                break
            for i, x in zip(sites, increments[0]):
                values[i] += x
            contributions.append(
                {
                    "sites": sites,
                    "types": types,
                    "path": path,
                    "row_ids": sorted({r["id"] for r, p in matches}),
                    "endpoint_values": [str(x) for x in increments[0]],
                    "searches": searches,
                }
            )
            found = True
            break
        if not found:
            missing.append(
                {
                    "sites": sites,
                    "reason": "source_parameter_missing",
                    "searches": searches,
                }
            )
    adj = {i: set() for i in labels}
    for b in system.topology.bonds.values():
        adj[b.site1].add(b.site2)
        adj[b.site2].add(b.site1)
    seen = set()
    components = []
    for i in sorted(adj):
        if i in seen:
            continue
        todo = [i]
        c = set()
        while todo:
            j = todo.pop()
            if j in c:
                continue
            c.add(j)
            todo.extend(adj[j] - c)
        seen |= c
        formal = sum(system.topology.sites[k].formal_charge for k in c)
        total = sum((values[k] for k in c), Decimal(0))
        components.append(
            {
                "sites": sorted(c),
                "formal": formal,
                "known_contribution_total": str(total),
                "residual": str(total - formal),
                "coverage_complete": not any(set(m["sites"]) & c for m in missing),
            }
        )
    return {
        "partial_charges": {i: str(v) for i, v in values.items()},
        "contributions": contributions,
        "missing": missing,
        "components": components,
        "complete": not missing
        and all(abs(Decimal(c["residual"])) <= Decimal("1e-12") for c in components),
    }


def raw_inventory(system):
    adj = {i: set() for i in system.topology.sites}
    edges = set()
    for b in system.topology.bonds.values():
        edges.add(tuple(sorted((b.site1, b.site2))))
        adj[b.site1].add(b.site2)
        adj[b.site2].add(b.site1)
    result = []

    def emit(f, ids):
        result.append((f, tuple(ids)))

    for i in sorted(adj):
        emit("nonbond(9-6)", [i])
    for b in sorted(edges):
        emit("quartic_bond", b)
    angles = [(a, j, b) for j in adj for a, b in combinations(sorted(adj[j]), 2)]
    for f in ("quartic_angle", "bond-bond", "bond-angle"):
        for a in angles:
            emit(f, a)
    torsions = set()
    for b, c in edges:
        for a in adj[b] - {c}:
            for d in adj[c] - {b}:
                if a != d:
                    t = (a, b, c, d)
                    torsions.add(min(t, t[::-1]))
    for f in (
        "torsion_3",
        "bond-bond_1_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-torsion_3",
        "angle-angle-torsion_1",
    ):
        for t in sorted(torsions):
            emit(f, t)
    for j, ns in adj.items():
        for arms in combinations(sorted(ns), 3):
            emit("wilson_out_of_plane", (arms[0], j, arms[1], arms[2]))
            for k in arms:
                a, b = sorted(set(arms) - {k})
                emit("angle-angle", (a, j, k, b))
    return result
