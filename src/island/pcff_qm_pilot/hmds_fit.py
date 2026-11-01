"""Explicit linear coupled-model fitting, with a training-only energy intercept.

Numerical fitting to a specified protocol is distinct from physical validation.
All fixed contributions come from the source-bound J64 model; no missing
coefficient is defaulted to zero in a complete physical evaluation.
"""

import numpy as np
from scipy.optimize import lsq_linear

from island.pcff_completion import require

from .hmds_baseline import extended_parts
from .hmds_design import metrics


def features(model, observations):
    require(observations, "nonempty_fit_observations_required")
    rows = []
    for item in observations:
        x = np.asarray(item["coordinates"], float)
        forces = np.asarray(item["forces"], float)
        require(
            x.shape == forces.shape == (27, 3)
            and np.isfinite(x).all()
            and np.isfinite(forces).all()
            and np.isfinite(item["relative_energy_kcal_mol"]),
            "finite_complete_HMDS_fit_observation_required",
        )
        rows.append(extended_parts(model, x))
    names = rows[0]["columns"]
    require(
        len(names) == 58 and all(r["columns"] == names for r in rows),
        "complete58_fit_mapping_required",
    )
    return {
        "names": names,
        "energy": np.array([r["energy_jacobian"] for r in rows]),
        "force": np.array([r["force_jacobian"] for r in rows]),
        "fixed_energy": np.array([r["fixed_energy"] for r in rows]),
        "fixed_force": np.array([r["fixed_forces"] for r in rows]),
        "QM_energy": np.array([r["relative_energy_kcal_mol"] for r in observations]),
        "QM_force": np.array([r["forces"] for r in observations]),
    }


def weighted_problem(data, contract, indices=None):
    ids = np.arange(len(data["energy"])) if indices is None else np.array(indices)
    require(len(ids) > 0, "nonempty_training_selection_required")
    n = len(ids)
    e = np.c_[data["energy"][ids], np.ones(n)]
    f = np.c_[data["force"][ids].reshape(81 * n, 58), np.zeros(81 * n)]
    ey = data["QM_energy"][ids] - data["fixed_energy"][ids]
    fy = (data["QM_force"][ids] - data["fixed_force"][ids]).ravel()
    ew = contract["energy_RMSE_scale_kcal_mol"] * np.sqrt(n)
    fw = contract["force_RMS_scale_kcal_mol_A"] * np.sqrt(81 * n)
    return np.vstack([e / ew, f / fw]), np.r_[ey / ew, fy / fw]


def infer(data, specifications, contract, indices=None):
    require([s["name"] for s in specifications] == data["names"], "fit_bound_mapping")
    a, y = weighted_problem(data, contract, indices)
    report = metrics(a)
    norms = np.linalg.norm(a, axis=0)
    safe = np.maximum(norms, 1e-30)
    scaled = a / safe
    lo = np.r_[[s["lower"] for s in specifications], -np.inf]
    hi = np.r_[[s["upper"] for s in specifications], np.inf]
    require(np.all(lo < hi), "invalid_fit_bounds")
    solve = lsq_linear(
        scaled,
        y,
        bounds=(lo * safe, hi * safe),
        method="trf",
        tol=1e-12,
        max_iter=1000,
        lsmr_tol=1e-12,
    )
    values = solve.x / safe
    distance = np.minimum(values[:-1] - lo[:-1], hi[:-1] - values[:-1])
    active = distance <= contract["bound_activity_fraction"] * (hi[:-1] - lo[:-1])
    unconstrained = np.linalg.lstsq(
        scaled, y, rcond=contract["rank_relative_tolerance"]
    )[0]
    lower = float(np.linalg.norm(scaled @ unconstrained - y) ** 2)
    return {
        "values": values.tolist(),
        "parameter_names": data["names"] + ["energy_offset"],
        "solver_success": bool(solve.success),
        "solver_status": int(solve.status),
        "solver_message": solve.message,
        "iterations": int(solve.nit or 0),
        "optimality": float(solve.optimality),
        "actual_weighted_design": report,
        "active_bound_parameters": [
            n for n, flag in zip(data["names"], active, strict=True) if flag
        ],
        "weighted_objective": float(np.linalg.norm(a @ values - y) ** 2),
        "unconstrained_objective_lower_bound": lower,
        "simultaneous_energy_force_RMSE_targets_impossible_in_linear_model": lower
        > 2 + 1e-8,
        "lower_bound_scope": "Global unbounded linear least-squares lower bound on mean(Eerror/0.1)^2+mean(Ferror/0.5)^2; both RMSE gates imply objective<=2. Does not prove infeasibility when lower bound<=2.",
        "identifiable": report["rank"] == len(values)
        and report["scaled_condition"] <= contract["maximum_scaled_condition"],
    }


def predictions(data, values):
    p = np.asarray(values, float)
    require(
        p.shape == (59,) and np.isfinite(p).all(), "all58_and_training_offset_required"
    )
    return data["fixed_energy"] + data["energy"] @ p[:-1] + p[-1], data[
        "fixed_force"
    ] + data["force"] @ p[:-1]


def performance(data, values, ids):
    e, f = predictions(data, values)
    de, df = e - data["QM_energy"], f - data["QM_force"]
    records = []
    for index, name in enumerate(ids):
        local = np.linalg.norm(df[index], axis=1)
        records.append(
            {
                "id": name,
                "energy_error_kcal_mol": float(de[index]),
                "force_RMS_kcal_mol_A": float(np.sqrt(np.mean(df[index] ** 2))),
                "maximum_force_component": float(abs(df[index]).max()),
                "maximum_local_force": float(local.max()),
                "worst_stable_atom_id": int(local.argmax()) + 1,
                "atom_wise_force_errors": local.tolist(),
                "Cartesian_force_errors": df[index].tolist(),
            }
        )
    return {
        "energy_RMSE": float(np.sqrt(np.mean(de**2))),
        "force_RMS": float(np.sqrt(np.mean(df**2))),
        "maximum_absolute_energy_error": float(abs(de).max()),
        "maximum_local_force_error": float(np.linalg.norm(df, axis=2).max()),
        "per_configuration": records,
        "energy_offset_from_training": float(values[-1]),
        "validation_offsets_refitted": False,
    }


def accepts(report, contract):
    return all(
        report[key] <= contract[limit]
        for key, limit in (
            ("energy_RMSE", "energy_RMSE_max"),
            ("force_RMS", "force_RMS_max"),
            ("maximum_absolute_energy_error", "maximum_absolute_energy_error"),
            ("maximum_local_force_error", "maximum_local_force_error"),
        )
    )


def prediction_uncertainty(data, values, nominal, ids):
    require(len(values) >= 2, "two_bootstrap_predictions_required")
    ep, fp = predictions(data, nominal)
    er, fr = [], []
    for v in values:
        e, f = predictions(data, v)
        er.append(e - ep)
        fr.append(f - fp)
    er, fr = np.array(er), np.array(fr)
    erms = np.sqrt(np.mean(er**2, axis=1))
    frms = np.sqrt(np.mean(fr**2, axis=(1, 2, 3)))
    return {
        "energy_RMS95": float(np.quantile(erms, 0.95)),
        "force_RMS95": float(np.quantile(frms, 0.95)),
        "per_configuration": [
            {
                "id": name,
                "energy_absolute95": float(np.quantile(abs(er[:, i]), 0.95)),
                "force_RMS95": float(
                    np.quantile(np.sqrt(np.mean(fr[:, i] ** 2, axis=(1, 2))), 0.95)
                ),
            }
            for i, name in enumerate(ids)
        ],
        "includes_model_discrepancy_or_basis_uncertainty": False,
    }


def uncertainty(data, observations, nominal, specifications, contract):
    spec = contract["uncertainty"]
    blocks = sorted({r["block"] for r in observations})
    grouped = {
        b: [i for i, r in enumerate(observations) if r["block"] == b] for b in blocks
    }
    rng = np.random.default_rng(spec["seed"])
    attempts, values = [], []
    for i in range(spec["replicates"]):
        selected = rng.choice(blocks, len(blocks), replace=True).tolist()
        indices = [j for b in selected for j in grouped[b]]
        result = infer(data, specifications, contract, indices)
        eligible = (
            result["identifiable"]
            and result["solver_success"]
            and not result["active_bound_parameters"]
        )
        attempts.append(
            {
                "replicate": i,
                "sampled_blocks": selected,
                "rank": result["actual_weighted_design"]["rank"],
                "scaled_condition": result["actual_weighted_design"][
                    "scaled_condition"
                ],
                "solver_success": result["solver_success"],
                "bound_activity": result["active_bound_parameters"],
                "eligible": bool(eligible),
                "values": result["values"],
            }
        )
        if eligible:
            values.append(result["values"])
    fraction = len(values) / spec["replicates"]
    answer = {
        "blocks": blocks,
        "attempts": attempts,
        "replacement_resamples": 0,
        "eligible_replicates": len(values),
        "eligible_fraction": fraction,
        "failed_replicates_retained": spec["replicates"] - len(values),
        "qualified": False,
        "parameter_intervals": None,
        "prediction_spread": None,
        "scope": spec["interval"],
    }
    if len(values) >= 2:
        array = np.array(values)
        intervals = np.quantile(array, [0.025, 0.975], axis=0)
        covariance = np.cov(array, rowvar=False)
        sd = np.sqrt(np.maximum(np.diag(covariance), 0))
        correlation = np.divide(
            covariance,
            np.outer(sd, sd),
            out=np.zeros_like(covariance),
            where=np.outer(sd, sd) > 0,
        )
        answer["parameter_intervals"] = [
            {
                "name": n,
                "nominal": float(v),
                "absolute_standard_deviation": float(s),
                "lower": float(l),
                "upper": float(u),
            }
            for n, v, s, l, u in zip(
                data["names"] + ["energy_offset"], nominal, sd, *intervals, strict=True
            )
        ]
        answer["covariance"] = covariance.tolist()
        answer["correlations"] = correlation.tolist()
        spread = prediction_uncertainty(
            data, values, nominal, [r["id"] for r in observations]
        )
        e95, f95 = spread["energy_RMS95"], spread["force_RMS95"]
        answer["prediction_spread"] = {**spread, "evaluation_domain": "training"}
        answer["qualified"] = bool(
            fraction >= spec["minimum_identifiable_fraction"]
            and e95 <= spec["prediction_energy_RMS95_max"]
            and f95 <= spec["prediction_force_RMS95_max"]
        )
    return answer
