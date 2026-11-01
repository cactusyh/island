"""Research-only scaled conic fitting; native QM execution is not available here.

The physical residual arrays remain authoritative. Solver tolerances never
modify an acceptance limit. Optional conic backends are imported only on use.
"""

import time

import numpy as np
from scipy import sparse

from island.pcff_completion import require


def physical_problem(data, indices=None):
    ids = (
        np.arange(len(data["energy"])) if indices is None else np.asarray(indices, int)
    )
    require(ids.ndim == 1 and len(ids) > 0, "nonempty_conic_selection")
    n, atoms = len(ids), data["force"].shape[1]
    p = len(data["names"])
    e = np.c_[data["energy"][ids], np.ones(n)]
    f = np.c_[data["force"][ids].reshape(3 * atoms * n, p), np.zeros(3 * atoms * n)]
    ey = data["QM_energy"][ids] - data["fixed_energy"][ids]
    fy = (data["QM_force"][ids] - data["fixed_force"][ids]).ravel()
    require(
        e.shape == (n, p + 1) and f.shape == (3 * atoms * n, p + 1), "conic_mapping"
    )
    require(all(np.isfinite(v).all() for v in (e, f, ey, fy)), "finite_conic_arrays")
    w = np.vstack([e / (0.1 * np.sqrt(n)), f / (0.5 * np.sqrt(3 * atoms * n))])
    y = np.r_[ey / (0.1 * np.sqrt(n)), fy / (0.5 * np.sqrt(3 * atoms * n))]
    return {
        "e": e,
        "f": f,
        "ey": ey,
        "fy": fy,
        "w": w,
        "y": y,
        "n": n,
        "atoms": atoms,
        "indices": ids,
    }


def verify(problem, values, lower, upper, rho=1.0):
    """Independent original-unit checks, without solver scaling/canonical rows."""
    x = np.asarray(values, float)
    require(
        x.shape == (problem["e"].shape[1],) and np.isfinite(x).all(),
        "complete_conic_vector",
    )
    de = problem["e"] @ x - problem["ey"]
    df = (problem["f"] @ x - problem["fy"]).reshape(problem["n"], problem["atoms"], 3)
    local = np.linalg.norm(df, axis=2)
    scores = {
        "energy_RMSE": float(np.sqrt(np.mean(de**2))),
        "force_RMS": float(np.sqrt(np.mean(df**2))),
        "maximum_absolute_energy_error": float(abs(de).max()),
        "maximum_local_force_error": float(local.max()),
    }
    ratios = [
        scores["energy_RMSE"] / 0.1,
        scores["force_RMS"] / 0.5,
        scores["maximum_absolute_energy_error"] / 0.3,
        scores["maximum_local_force_error"] / 2,
    ]
    distances = np.minimum(x[:-1] - lower, upper - x[:-1])
    slacks = {
        "energy_L2": 0.1 * np.sqrt(problem["n"]) * rho - np.linalg.norm(de),
        "force_L2": 0.5 * np.sqrt(df.size) * rho - np.linalg.norm(df),
        "configuration_energy": (0.3 * rho - abs(de)).tolist(),
        "local_vector_force": (2 * rho - local).tolist(),
        "lower_bound": (x[:-1] - lower).tolist(),
        "upper_bound": (upper - x[:-1]).tolist(),
    }
    return dict(
        **scores,
        required_rho=float(max(ratios)),
        original_limits_satisfied=bool(max(ratios) <= 1 and np.all(distances >= 0)),
        rho_constraints_satisfied=bool(max(ratios) <= rho and np.all(distances >= 0)),
        slacks=slacks,
        parameter_bound_distances=distances.tolist(),
        active_bound_indices=np.flatnonzero(
            distances <= 1e-6 * (upper - lower)
        ).tolist(),
        signed_energy_errors=de.tolist(),
        signed_Cartesian_force_errors=df.tolist(),
        local_vector_force_errors=local.tolist(),
        weighted_objective=float(np.linalg.norm(problem["w"] @ x - problem["y"]) ** 2),
    )


def canonical(
    problem, lower, upper, kind="rho", interior=0.0, bounded=True, force_only=False
):
    """A x+s=b: linear nonnegative cone then self-dual Lorentz cones.

    Physical coefficients p=z/column_norm; offset is scaled the same way.
    Each residual constraint is divided by its frozen physical limit.
    """
    require(kind in ("rho", "least_squares"), "unknown_conic_objective")
    lo, hi = np.array(lower, float), np.array(upper, float)
    require(
        lo.shape == hi.shape == (problem["e"].shape[1] - 1,)
        and np.all(lo < hi)
        and np.isfinite(lo).all()
        and np.isfinite(hi).all()
        and 0 <= interior < 0.5,
        "conic_bounds",
    )
    scale = np.maximum(np.linalg.norm(problem["w"], axis=0), 1e-30)
    dim = len(scale) + 1
    rows, rhs, labels = [], [], []
    lower_i, upper_i = lo + interior * (hi - lo), hi - interior * (hi - lo)
    if bounded:
        for j in range(len(lo)):
            for sign, limit, label in (
                (1, upper_i[j], "upper"),
                (-1, -lower_i[j], "lower"),
            ):
                row = np.zeros(dim)
                row[j] = sign / (scale[j] * (hi[j] - lo[j]))
                rows.append(row)
                rhs.append(limit / (hi[j] - lo[j]))
                labels.append([label, j])
    if kind == "rho":
        row = np.zeros(dim)
        row[-1] = -1
        rows.append(row)
        rhs.append(0.0)
        labels.append(["rho_nonnegative"])
    e, f = problem["e"] / scale, problem["f"] / scale
    if not force_only:
        for i in range(problem["n"]):
            for sign in (1, -1):
                row = np.zeros(dim)
                row[: len(scale)] = sign * e[i] / 0.3
                if kind == "rho":
                    row[-1] = -1
                rows.append(row)
                rhs.append(sign * problem["ey"][i] / 0.3 + (kind != "rho"))
                labels.append(["energy_absolute", i, sign])
    linear = len(rows)
    cone_sizes = []

    def soc(matrix, target, limit, label):
        block = np.zeros((len(target) + 1, dim))
        if kind == "rho":
            block[0, -1] = -1
        block[1:, : len(scale)] = -matrix / limit
        rows.extend(block)
        rhs.extend(np.r_[float(kind != "rho"), -target / limit])
        cone_sizes.append(len(target) + 1)
        labels.append([label, len(rows) - len(block), len(rows)])

    if not force_only:
        soc(e, problem["ey"], 0.1 * np.sqrt(problem["n"]), "energy_L2")
    soc(f, problem["fy"], 0.5 * np.sqrt(len(f)), "force_L2")
    for i in range(len(f) // 3):
        soc(
            f[3 * i : 3 * i + 3],
            problem["fy"][3 * i : 3 * i + 3],
            2.0,
            ["local_vector", i // problem["atoms"], i % problem["atoms"] + 1],
        )
    if kind == "least_squares":
        # ||[2r,t-1]|| <= t+1 is exactly ||r||^2<=t; avoids
        # subtraction of a large QP objective constant in gap termination.
        w = problem["w"] / scale
        block = np.zeros((len(w) + 2, dim))
        block[0, -1] = -1
        block[1:-1, : len(scale)] = -2 * w
        block[-1, -1] = -1
        start = len(rows)
        rows.extend(block)
        rhs.extend(np.r_[1.0, -2 * problem["y"], -1.0])
        cone_sizes.append(len(block))
        labels.append(["objective_epigraph", start, len(rows)])
    a = sparse.csc_matrix(np.array(rows))
    q = np.zeros(dim)
    q[-1] = 1.0
    p = sparse.csc_matrix((dim, dim))
    constant = 0.0
    return {
        "A": a,
        "b": np.array(rhs),
        "P": p,
        "q": q,
        "linear": linear,
        "soc": cone_sizes,
        "scale": scale,
        "labels": labels,
        "objective_constant": constant,
        "kind": kind,
        "interior": interior,
        "bounded": bounded,
        "force_only": force_only,
    }


def solve(system, backend="clarabel", settings=None):
    """Retain full primal/dual/slack state; failed statuses remain inconclusive."""
    started = time.monotonic()
    a, b, p, q = (system[k] for k in ("A", "b", "P", "q"))
    if backend == "clarabel":
        import clarabel

        config = clarabel.DefaultSettings()
        defaults = {
            "verbose": False,
            "max_iter": 300,
            "time_limit": 120.0,
            "tol_gap_abs": 1e-9,
            "tol_gap_rel": 1e-9,
            "tol_feas": 1e-10,
        }
        defaults.update(settings or {})
        for k, v in defaults.items():
            setattr(config, k, v)
        cones = [clarabel.NonnegativeConeT(system["linear"])] + [
            clarabel.SecondOrderConeT(n) for n in system["soc"]
        ]
        optimizer = clarabel.DefaultSolver(
            sparse.triu(p).tocsc(), q, a, b, cones, config
        )
        result = optimizer.solve()
        info = optimizer.get_info()
        record = {
            "status": str(result.status),
            "x": list(result.x),
            "dual": list(result.z),
            "slack": list(result.s),
            "primal_objective": float(result.obj_val),
            "dual_objective": float(result.obj_val_dual),
            "iterations": result.iterations,
            "primal_residual": float(result.r_prim),
            "dual_residual": float(result.r_dual),
            "absolute_gap": float(info.gap_abs),
            "relative_gap": float(info.gap_rel),
            "version": clarabel.__version__,
        }
    elif backend == "scs":
        import scs

        defaults = {
            "verbose": False,
            "max_iters": 100000,
            "time_limit_secs": 120.0,
            "eps_abs": 1e-7,
            "eps_rel": 1e-7,
        }
        defaults.update(settings or {})
        result = scs.solve(
            {"A": a, "b": b, "c": q, "P": sparse.triu(p).tocsc()},
            {"l": system["linear"], "q": system["soc"]},
            **defaults,
        )
        info = result["info"]
        record = {
            "status": info["status"],
            "x": result["x"].tolist(),
            "dual": result["y"].tolist(),
            "slack": result["s"].tolist(),
            "primal_objective": info["pobj"],
            "dual_objective": info["dobj"],
            "iterations": info["iter"],
            "primal_residual": info["res_pri"],
            "dual_residual": info["res_dual"],
            "absolute_gap": info["gap"],
            "version": scs.__version__,
        }
    else:
        raise ValueError("unsupported conic backend")
    record.update(
        backend=backend, settings=defaults, runtime_seconds=time.monotonic() - started
    )
    x = np.asarray(record["x"])
    record["values"] = (x[: len(system["scale"])] / system["scale"]).tolist()
    record["rho"] = float(x[-1]) if system["kind"] == "rho" else 1.0
    record["objective_including_constant"] = (
        record["primal_objective"] + system["objective_constant"]
    )
    return record


def classify(result, verification, checked_lower_bound=None):
    """Physical witness precedes optimality, failure never proves incompatibility."""
    if verification is not None and verification["original_limits_satisfied"]:
        return "verified_feasible_witness"
    if checked_lower_bound is not None and checked_lower_bound > 1 + 1e-5:
        return "checked_lower_bound_above_one"
    if result["status"] in ("Solved", "solved"):
        return (
            "boundary_or_accuracy_limited"
            if abs(result.get("rho", 1) - 1) < 1e-5
            else "no_certificate_inconclusive"
        )
    return "solver_failure_inconclusive"


def save_system(file, system):
    values = {k: system[k] for k in ("b", "q", "scale")}
    for name in ("A", "P"):
        m = system[name]
        values.update(
            {
                name + "_data": m.data,
                name + "_indices": m.indices,
                name + "_indptr": m.indptr,
                name + "_shape": m.shape,
            }
        )
    np.savez_compressed(file, **values)


def dual_lower_bound(system, result, problem, lower, upper):
    """Check a conic-dual/support-function bound in extended precision.

    No unbounded-variable stationarity residual is ignored. A physical primal
    witness supplies U; rho<=U then bounds the energy intercept via row zero.
    The global optimum lies in this box, since a feasible rho=U witness exists.
    Projected cone duals give L(x)=q.x+z.(Ax-b)<=q.x for every feasible x.
    Minimize this affine L over that box, including its stationarity error.
    """
    require(
        system["kind"] == "rho" and system["bounded"], "finite_box_certificate_required"
    )
    dtype = np.longdouble
    lo, hi = np.array(lower, float), np.array(upper, float)
    values = np.array(result["values"], float)
    values[:-1] = np.clip(values[:-1], lo, hi)
    check = verify(problem, values, lo, hi)
    upper_rho = check["required_rho"] + 1e-6
    z = np.maximum(np.asarray(result["dual"], dtype=dtype), 0)
    # Only the initial block is a nonnegative cone; signed SOC tails retained.
    raw = np.asarray(result["dual"], dtype=dtype)
    offset = system["linear"]
    z[offset:] = raw[offset:]
    margins = []
    for size in system["soc"]:
        tail_norm = np.sqrt(np.sum(z[offset + 1 : offset + size] ** 2, dtype=dtype))
        z[offset] = max(z[offset], tail_norm + dtype("1e-14") * (1 + tail_norm))
        margins.append(float(z[offset] - tail_norm))
        offset += size
    a = np.asarray(system["A"].toarray(), dtype=dtype)
    b, q = np.asarray(system["b"], dtype=dtype), np.asarray(system["q"], dtype=dtype)
    gradient = q + a.T @ z
    scale = np.asarray(system["scale"], dtype=dtype)
    e0 = np.asarray(problem["e"][0, :-1], dtype=dtype)
    ey0 = dtype(problem["ey"][0])
    low_e = np.minimum(e0 * lo, e0 * hi).sum(dtype=dtype)
    high_e = np.maximum(e0 * lo, e0 * hi).sum(dtype=dtype)
    intercept_lo, intercept_hi = (
        ey0 - high_e - dtype(".3") * upper_rho,
        ey0 - low_e + dtype(".3") * upper_rho,
    )
    box_lo = np.r_[lo * scale[:-1], intercept_lo * scale[-1], dtype(0)]
    box_hi = np.r_[hi * scale[:-1], intercept_hi * scale[-1], dtype(upper_rho)]
    if system["force_only"]:
        # Intercept is absent; gradient is exactly zero in this column.
        require(gradient[-2] == 0, "force_only_intercept_not_zero")
        box_lo[-2] = box_hi[-2] = 0
    support = np.minimum(gradient * box_lo, gradient * box_hi)
    dual_value = -(b * z).sum(dtype=dtype)
    raw_bound = dual_value + support.sum(dtype=dtype)
    # Conservative absolute numerical safety allowance, much larger than the
    # long-double dot-product roundoff for these finite arrays.
    magnitude = abs(b * z).sum(dtype=dtype) + abs(a * z[:, None]).sum(
        dtype=dtype
    ) * max(np.max(abs(box_lo)), np.max(abs(box_hi)), dtype(1))
    allowance = dtype("1e-12") * max(magnitude, dtype(1))
    return {
        "lower_bound": float(raw_bound - allowance),
        "raw_support_bound": float(raw_bound),
        "arithmetic_allowance": float(allowance),
        "primal_upper_bound": upper_rho,
        "upper_bound_witness": values.tolist(),
        "dual_cone_margins": margins,
        "projected_dual": z.astype(float).tolist(),
        "stationarity_residual": gradient.astype(float).tolist(),
        "support_correction": float(support.sum(dtype=dtype)),
        "raw_dual_value": float(dual_value),
        "box_lower": box_lo.astype(float).tolist(),
        "box_upper": box_hi.astype(float).tolist(),
        "original_finite_bounds": True,
        "intercept_error_bounded": True,
        "sufficient_margin": bool(raw_bound - allowance > 1 + 1e-5),
        "derivation": "For feasible rho<=primalU, all variables lie in this derived box. Project dual into self-dual cones, then -b.z+min_box((q+A.Tz).x) bounds rho from below. Include every stationarity residual and conservative arithmetic allowance. rho>U already exceeds this bound.",
    }
