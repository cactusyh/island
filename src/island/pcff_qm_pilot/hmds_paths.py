"""Graph-defined HMDS fragment paths and scalable exact duplicate screening."""

from collections import Counter
from math import factorial

import numpy as np

from island.pcff_completion import require

from .hamiltonian import bend, length
from .molecular_control import cartesian, torsion

GRAPH = "d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd"


def adjacency(model):
    require(
        model["graph_identity"] == GRAPH and model["atom_order"] == list(range(1, 28)),
        "unsupported_HMDS_graph",
    )
    a = {i: [] for i in model["atom_order"]}
    for i, j, order in model["bonds"]:
        require(order == 1 and i in a and j in a, "unsupported_HMDS_bond")
        a[i].append(j)
        a[j].append(i)
    require(len(model["bonds"]) == 26, "HMDS_tree_required")
    return a


def fragment(model, bond, moving_root):
    a = adjacency(model)
    i, j = bond
    require(j in a[i] and moving_root in bond, "invalid_fragment_cut")
    seen = set()
    todo = [moving_root]
    while todo:
        n = todo.pop()
        if n in seen:
            continue
        seen.add(n)
        todo.extend(k for k in a[n] if {k, n} != {i, j})
    require(len(seen) < 27 and (i in seen) != (j in seen), "cyclic_or_disconnected_cut")
    return sorted(seen)


def path(model, center, spec, s):
    x = np.asarray(center, float)
    require(x.shape == (27, 3) and np.isfinite(x).all(), "HMDS_path_mapping")
    i, j = spec["bond"]
    moving = fragment(model, [i, j], spec["moving_root"])
    require(spec["moving_root"] == j, "explicit_path_axis_orientation")
    axis = x[j - 1] - x[i - 1]
    axis /= np.linalg.norm(axis)
    indices = np.array(moving) - 1
    y = x.copy()
    t = np.zeros_like(x)
    a = np.zeros_like(x)
    if spec["kind"] == "stretch":
        require(np.linalg.norm(x[j - 1] - x[i - 1]) + s > 0.5, "invalid_stretch_range")
        y[indices] += s * axis
        t[indices] = axis
    else:
        require(spec["kind"] == "rotation", "unsupported_path_kind")
        v = x[indices] - x[i - 1]
        parallel = np.outer(v @ axis, axis)
        perp = v - parallel
        c = np.cross(axis, perp)
        y[indices] = x[i - 1] + parallel + np.cos(s) * perp + np.sin(s) * c
        t[indices] = -np.sin(s) * perp + np.cos(s) * c
        a[indices] = -np.cos(s) * perp - np.sin(s) * c
    return {
        "coordinates": y,
        "tangent": t,
        "acceleration": a,
        "moving_ids": moving,
        "parameter": float(s),
        "parameter_units": "angstrom" if spec["kind"] == "stretch" else "radian",
    }


def affected(model, center, displaced):
    x, y = cartesian(center), cartesian(displaced)
    changes = []
    for row in model["interaction_inventory"]:
        f, s = row["family"], row["sites"]
        if f not in ("quartic_bond", "quartic_angle", "torsion_3"):
            continue
        fun = {"quartic_bond": length, "quartic_angle": bend, "torsion_3": torsion}[f]
        d = fun(y, s).v - fun(x, s).v
        if f == "torsion_3":
            d = (d + np.pi) % (2 * np.pi) - np.pi
        if abs(d) > 1e-9:
            changes.append(
                {"id": row["id"], "sites": s, "family": f, "change": float(d)}
            )
    return changes


def derivative(sample, p):
    return float(-np.sum(np.asarray(sample["forces"]) * p["tangent"]))


def response(center, plus, minus, p0, pp, pm, h):
    gp, gm, g0 = [-np.array(r["forces"]) for r in (plus, minus, center)]
    dg = (gp - gm) / (2 * h)
    contraction = float(np.sum(p0["tangent"] * dg))
    acceleration = float(np.sum(g0 * p0["acceleration"]))
    return {
        "Cartesian_gradient_response": dg.tolist(),
        "path_first_derivative_plus": derivative(plus, pp),
        "path_first_derivative_minus": derivative(minus, pm),
        "curvature_from_path_gradient": (derivative(plus, pp) - derivative(minus, pm))
        / (2 * h),
        "curvature_from_energies": (
            plus["energy"] - 2 * center["energy"] + minus["energy"]
        )
        / h**2,
        "t_H_t_finite_step": contraction,
        "gradient_path_acceleration": acceleration,
        "chain_rule_sum": contraction + acceleration,
        "interpretation": "Total molecular path response, not an isolated source coefficient.",
    }


def automorphism_count(model):
    a = adjacency(model)

    def rooted(i, parent):
        children = [rooted(j, i) for j in a[i] if j != parent]
        counts = Counter(s for s, n in children)
        number = 1
        for s, n in children:
            number *= n
        for n in counts.values():
            number *= factorial(n)
        return (model["typing"][str(i)], tuple(sorted(s for s, n in children))), number

    root = next(i for i in a if model["typing"][str(i)] == "osi")
    return rooted(root, None)[1]


def equivalent(model, left, right, tolerance=1e-6, maximum_states=100000):
    """Distance-colored graph backtracking covers all allowed H/methyl/end maps.

    Generic geometry fingerprints prune the 3,359,232 graph automorphisms without
    enumerating them. Exhaustion is fail-closed, never evidence of independence.
    Final alignment uses proper rotations, not reflections.
    """
    a = adjacency(model)
    x = np.asarray(left)
    y = np.asarray(right)
    require(
        x.shape == y.shape == (27, 3) and np.isfinite(x).all() and np.isfinite(y).all(),
        "leakage_mapping",
    )
    dx = np.linalg.norm(x[:, None] - x[None, :], axis=2)
    dy = np.linalg.norm(y[:, None] - y[None, :], axis=2)
    labels = [model["typing"][str(i)] for i in range(1, 28)]
    kinds = sorted(set(labels))
    # Alignment accepts coordinate RMS over81 components. A necessary pair-
    # distance bound must allow the largest individual residual compatible
    # with that RMS, rather than pruning with an unjustified4*RMS cutoff.
    tol = 2 * np.sqrt(x.size) * tolerance
    fingerprints = lambda d: [
        np.concatenate(
            [np.sort(d[i, [j for j, v in enumerate(labels) if v == k]]) for k in kinds]
        )
        for i in range(27)
    ]
    fx, fy = fingerprints(dx), fingerprints(dy)
    candidates = {
        i: [
            j
            for j in range(27)
            if labels[i] == labels[j] and np.max(abs(fx[i] - fy[j])) <= tol
        ]
        for i in range(27)
    }
    if any(not v for v in candidates.values()):
        return False
    order = sorted(candidates, key=lambda i: (len(candidates[i]), -len(a[i + 1])))
    mapping = {}
    used = set()
    states = 0

    def visit(depth):
        nonlocal states
        states += 1
        require(states <= maximum_states, "HMDS_leakage_search_inconclusive")
        if depth == 27:
            z = y[[mapping[i] for i in range(27)]]
            xx = x - x.mean(axis=0)
            zz = z - z.mean(axis=0)
            u, _, vt = np.linalg.svd(xx.T @ zz)
            r = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
            return np.sqrt(np.mean((xx @ r - zz) ** 2)) <= tolerance
        i = order[depth]
        for j in candidates[i]:
            if j in used:
                continue
            if any(
                (b + 1 in a[i + 1]) != (c + 1 in a[j + 1])
                or abs(dx[i, b] - dy[j, c]) > tol
                for b, c in mapping.items()
            ):
                continue
            mapping[i] = j
            used.add(j)
            if visit(depth + 1):
                return True
            used.remove(j)
            del mapping[i]
        return False

    return visit(0)


def leakage(model, training, validation):
    require(len(validation) > 0, "validation_reservation_required")
    for i, v in enumerate(validation):
        for x in [*training, *validation[:i]]:
            require(not equivalent(model, x, v), "HMDS_validation_symmetry_leakage")
    return {
        "graph_automorphisms_covered": automorphism_count(model),
        "algorithm": "distance-colored graph backtracking plus proper Kabsch; search cap fails closed",
        "training_count": len(training),
        "reserved_count": len(validation),
    }
