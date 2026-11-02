"""Additive research diagnostics of an affine Hamiltonian span; no native engine.

Physical performance, original physical bounds, and numerical certification are
separate. A derived computational radius is never a new physical constraint.
"""

import numpy as np
from scipy import sparse

from island.pcff_completion import require

from .hmds_augmentation import protected_system
from .hmds_constrained import canonical, verify


def system(problem, lo, hi, protected=None, bounded=False):
    if bounded:
        return protected_system(problem, lo, hi, protected=protected)
    # Historical canonical validates finite shapes even with bounded=False.
    # These placeholders create ZERO constraint rows; added variables remain free.
    absent = problem["e"].shape[1] - 1 - len(lo)
    shape_lo, shape_hi = np.r_[lo, -np.ones(absent)], np.r_[hi, np.ones(absent)]
    out = canonical(problem, shape_lo, shape_hi, bounded=False)
    if protected is None:
        return out
    extra = canonical(protected, shape_lo, shape_hi, bounded=False)
    transform = sparse.diags(np.r_[extra["scale"] / out["scale"], 1.0])
    a = extra["A"] @ transform
    l, m = out["linear"], extra["linear"]
    out["A"] = sparse.vstack([out["A"][:l], a[:m], out["A"][l:], a[m:]], format="csc")
    out["b"] = np.r_[out["b"][:l], extra["b"][:m], out["b"][l:], extra["b"][m:]]
    out["linear"] += m
    out["soc"] += extra["soc"]
    out["protected"] = True
    out["labels"] = {"full": out["labels"], "protected": extra["labels"]}
    return out


def check(problem, protected, values, lo, hi, rho=1):
    checks = {}
    for name, p in [("full26", problem), ("protected22", protected)]:
        if p is None:
            continue
        extra = np.asarray(values)[len(lo) : -1]
        vlo, vhi = np.r_[lo, extra - 1], np.r_[hi, extra + 1]
        c = verify(p, values, vlo, vhi, rho)
        for key in ("parameter_bound_distances", "active_bound_indices"):
            c[key] = (
                c[key][: len(lo)]
                if key == "parameter_bound_distances"
                else [i for i in c[key] if i < len(lo)]
            )
        for key in ("lower_bound", "upper_bound"):
            c["slacks"][key] = c["slacks"][key][: len(lo)]
        c["performance_limits_satisfied"] = c["required_rho"] <= 1
        c["original_bounds_satisfied"] = bool(
            np.all(np.array(c["parameter_bound_distances"]) >= 0)
        )
        c["performance_rho_satisfied"] = c["required_rho"] <= rho
        c["original_limits_satisfied"] = (
            c["performance_limits_satisfied"] and c["original_bounds_satisfied"]
        )
        c["rho_constraints_satisfied"] = (
            c["performance_rho_satisfied"] and c["original_bounds_satisfied"]
        )
        c["original_limit_slacks"] = verify(p, values, vlo, vhi)["slacks"]
        for key in ("lower_bound", "upper_bound"):
            c["original_limit_slacks"][key] = c["original_limit_slacks"][key][: len(lo)]
        checks[name] = c
    x = np.array(values[: len(lo)])
    return {
        "performance_limits_satisfied": all(
            c["performance_limits_satisfied"] for c in checks.values()
        ),
        "original_bounds_satisfied": all(
            c["original_bounds_satisfied"] for c in checks.values()
        ),
        "required_rho": max(c["required_rho"] for c in checks.values()),
        "subsets": checks,
        "outside_original_bounds": [
            {
                "index": int(j),
                "value": float(x[j]),
                "lower": float(lo[j]),
                "upper": float(hi[j]),
                "signed_violation": float(x[j] - np.clip(x[j], lo[j], hi[j])),
                "width_fraction": float(
                    abs(x[j] - np.clip(x[j], lo[j], hi[j])) / (hi[j] - lo[j])
                ),
            }
            for j in np.flatnonzero((x < lo) | (x > hi))
        ],
    }


def radius(problem, scale, upper):
    """Checked left inverse bounds the smallest singular value, including offset.

    B is merely an approximate left inverse. Its construction need not be exact:
    ||B M-I||<1 is verified afresh in extended precision. Frobenius bounds
    dominate spectral norms. Conservative error allowances cover matrix division,
    dot products and norm rounding. Rank-deficient cases are inconclusive.
    """
    m = problem["w"] / scale
    u, s, vt = np.linalg.svd(m, full_matrices=False)
    rank = int(np.sum(s > max(m.shape) * np.finfo(float).eps * s[0]))
    if rank != m.shape[1]:
        return {
            "certified": False,
            "rank": rank,
            "reason": "no full-column-rank radius; nullspace not ignored",
        }
    b = (vt.T / s) @ u.T
    d = np.longdouble
    ml, bl = m.astype(d), b.astype(d)
    defect = bl @ ml - np.eye(m.shape[1], dtype=d)
    # Dominates long-double accumulation and double matrix division roundoff.
    allowance = d("1e-12") * max(d(1), np.linalg.norm(b) * np.linalg.norm(m))
    eps = np.sqrt(np.sum(defect**2, dtype=d)) + allowance
    bn = np.sqrt(np.sum(bl**2, dtype=d)) * (1 + d("1e-12"))
    if eps >= 1:
        return {"certified": False, "rank": rank, "reason": "left inverse error >=1"}
    sigma = (1 - eps) / bn
    yn = np.sqrt(np.sum(np.asarray(problem["y"], dtype=d) ** 2, dtype=d))
    r = (yn + np.sqrt(d(2)) * d(upper)) / sigma * (1 + d("1e-12"))
    return {
        "certified": True,
        "rank": rank,
        "singular_values": s.tolist(),
        "scaled_condition": float(s[0] / s[-1]),
        "left_inverse": b.tolist(),
        "left_inverse_defect_F": float(np.sqrt(np.sum(defect**2, dtype=d))),
        "left_inverse_roundoff_allowance": float(allowance),
        "left_inverse_norm_upper": float(bn),
        "sigma_lower": float(sigma),
        "scaled_coefficient_radius": float(r),
        "upper_rho": float(upper),
        "target_norm": float(yn),
        "includes_energy_only_intercept": True,
    }


def certificate(sys, result, problem, protected, lo, hi):
    """Bound-free support certificate with a derived finite computational ball."""
    require(not sys["bounded"], "bound_free_certificate_only")
    checked = check(problem, protected, result["values"], lo, hi)
    upper = checked["required_rho"] + 1e-6
    upper_check = check(problem, protected, result["values"], lo, hi, upper)
    require(
        all(c["performance_rho_satisfied"] for c in upper_check["subsets"].values()),
        "both_subset_upper_required",
    )
    region = radius(problem, sys["scale"], upper)
    if not region["certified"]:
        return {"certified": False, "region": region, "lower_bound": None}
    d = np.longdouble
    dual = np.asarray(result["dual"], dtype=d).copy()
    dual[: sys["linear"]] = np.maximum(dual[: sys["linear"]], 0)
    offset = sys["linear"]
    margins = []
    for size in sys["soc"]:
        norm = np.sqrt(np.sum(dual[offset + 1 : offset + size] ** 2, dtype=d))
        dual[offset] = max(dual[offset], norm + d("1e-14") * (1 + norm))
        margins.append(float(dual[offset] - norm))
        offset += size
    a, b, q = (np.asarray(v, dtype=d) for v in (sys["A"].toarray(), sys["b"], sys["q"]))
    g = q + a.T @ dual
    ball_support = -d(region["scaled_coefficient_radius"]) * np.sqrt(
        np.sum(g[:-1] ** 2, dtype=d)
    )
    rho_support = min(d(0), d(upper) * g[-1])
    dual_value = -np.sum(b * dual, dtype=d)
    magnitude = np.sum(abs(b * dual), dtype=d) + np.sum(
        abs(a * dual[:, None]), dtype=d
    ) * max(d(region["scaled_coefficient_radius"]), d(upper), d(1))
    allowance = d("1e-12") * max(d(1), magnitude)
    bound = dual_value + ball_support + rho_support - allowance
    return {
        "certified": True,
        "lower_bound": float(bound),
        "primal_upper_bound": upper,
        "upper_bound_checks": upper_check,
        "region": region,
        "projected_dual": dual.astype(float).tolist(),
        "dual_cone_margins": margins,
        "stationarity_residual": g.astype(float).tolist(),
        "ball_support": float(ball_support),
        "rho_support": float(rho_support),
        "raw_dual_value": float(dual_value),
        "arithmetic_allowance": float(allowance),
        "sufficient_margin": bool(bound > 1 + 1e-5),
        "derivation": "For both-subset checked U, ||M z-y||<=sqrt(2)U; full column rank gives ||z||<=R. Projected self-dual lambda gives rho>=-b.lambda+g.z. Minimize all stationarity residuals over ||z||<=R and 0<=rho<=U. No physical coefficient box imposed.",
    }


def classification(checks, cert):
    if checks["performance_limits_satisfied"]:
        return "independently_verified_performance_feasible_witness"
    if cert.get("certified") and cert["lower_bound"] > 1 + 1e-5:
        return "checked_lower_bound_above_one"
    if abs(checks["required_rho"] - 1) < 1e-5:
        return "boundary_or_accuracy_limited"
    return "numerically_inconclusive"


def projection(problem):
    scale = np.maximum(np.linalg.norm(problem["w"], axis=0), 1e-30)
    m = problem["w"] / scale
    z, _, rank, singular = np.linalg.lstsq(m, problem["y"], rcond=None)
    residual = m @ z - problem["y"]
    return {
        "values": (z / scale).tolist(),
        "rank": int(rank),
        "scale": scale.tolist(),
        "singular_values": singular.tolist(),
        "scaled_condition": float(singular[0] / singular[-1])
        if rank == m.shape[1]
        else None,
        "weighted_residual": residual.tolist(),
        "outside_span_target": (-residual).tolist(),
        "weighted_residual_norm": float(np.linalg.norm(residual)),
        "normal_equation_residual": (m.T @ residual).tolist(),
        "scope": "Weighted projection floor, not alone a four-limit incompatibility certificate or physical identifiability.",
    }


def relaxation_system(problem, lo, hi, protected=None):
    """Minimum normalized interval relaxation, fixed physical rho=1.

    eta is auxiliary, never a second energy intercept. Bounds are exactly
    [lo-eta*(hi-lo), hi+eta*(hi-lo)]; no arbitrary coefficient box.
    """
    out = system(problem, lo, hi, protected)
    a = out["A"].toarray()
    out["b"] = out["b"] - a[:, -1]
    a[:, -1] = 0
    rows, rhs = [], []
    for j in range(len(lo)):
        for sign, limit in ((1, hi[j]), (-1, -lo[j])):
            r = np.zeros(a.shape[1])
            r[j] = sign / (out["scale"][j] * (hi[j] - lo[j]))
            r[-1] = -1
            rows.append(r)
            rhs.append(limit / (hi[j] - lo[j]))
    r = np.zeros(a.shape[1])
    r[-1] = -1
    rows.append(r)
    rhs.append(0)
    l = out["linear"]
    out["A"] = sparse.csc_matrix(np.vstack([a[:l], rows, a[l:]]))
    out["b"] = np.r_[out["b"][:l], rhs, out["b"][l:]]
    out["linear"] += len(rows)
    out["auxiliary_role"] = "eta_normalized_bound_relaxation; physical rho fixed1"
    return out


def intake(
    model,
    rows,
    data,
    specifications,
    original_rows,
    original_arrays,
    original_specifications,
):
    """An outer seal does not authorize changed graph, mapping or physical units.

    The external immutable J67 records are the semantic reference, including
    every producer. This layer adds no replacement historical execution record.
    """
    from .hmds_hamiltonian import validate

    validate(model)
    require(
        rows == original_rows and len(rows) == 26,
        "J68_observation_protocol_producer_mapping",
    )
    require(
        specifications == original_specifications
        and [r["name"] for r in specifications] == data["names"],
        "J68_units_bounds_parameter_order",
    )
    require(
        data["names"] == original_arrays["names"]
        and all(
            np.array_equal(data[k], original_arrays[k]) for k in data if k != "names"
        ),
        "J68_physical_matrix_mutation",
    )
    return True


def authorize_attempt(ledger, name, backend, bounded, recover=False):
    require(
        len(ledger) <= 12 and len({r["name"] for r in ledger}) == len(ledger),
        "J68_solver_budget",
    )
    matches = [r for r in ledger if r["name"] == name]
    require(len(matches) == 1, "J68_count_before_worker")
    row = matches[0]
    require(
        row["backend"] == backend
        and row["bounded"] == bounded
        and row["status"] == ("failed" if recover else "started"),
        "J68_exact_attempt_authorization",
    )
    return True
