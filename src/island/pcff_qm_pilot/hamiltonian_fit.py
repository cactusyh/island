"""Conditional complete-model joint least squares; no automatic physical fit."""

from itertools import permutations, product

import numpy as np

from island.pcff_completion import digest, require

from .hamiltonian import parts


def dataset(model, samples):
    rows = []
    for sample in samples:
        require(
            sample["graph_identity"] == model["graph_identity"]
            and sample["atom_order"] == model["atom_order"],
            "fit_mapping_mutation",
        )
        require(sample["split"] in ("train", "validation_only"), "fit_split_required")
        require(sample.get("conformer_identity"), "conformer_identity_required")
        p = parts(model, sample["coordinates"])
        f = np.asarray(sample["forces"], float)
        require(
            f.shape == (9, 3)
            and np.isfinite(f).all()
            and np.isfinite(sample["energy"]),
            "invalid_reference_data",
        )
        rows.append((sample, p))
    return rows


def conformer_equivalent(left, right):
    """Graph-symmetry/distance identity, including rigid transformations.

    This is an identity rejection tolerance, not a tuned fit-error criterion.
    Slightly strained same-basin geometries still need explicit conformer labels.
    """
    a, b = np.asarray(left, float), np.asarray(right, float)
    da = np.linalg.norm(a[:, None] - a[None, :], axis=2)
    db = np.linalg.norm(b[:, None] - b[None, :], axis=2)
    for first, second in product(permutations((3, 4, 5)), permutations((6, 7, 8))):
        for order in ([0, 1, 2, *first, *second], [2, 1, 0, *second, *first]):
            if np.max(abs(da - db[np.ix_(order, order)])) < 1e-7:
                return True
    return False


def joint_fit(model, training, held_out, *, reference_gate_passed, criteria):
    """Fit QM-fixed with all unknowns explicit and untouched validation.

    Bounds/weights must be the predeclared J60 plan. Energy zero is anchored
    once at the first training sample, never refit independently per conformer.
    """
    require(reference_gate_passed is True, "reference_accuracy_gate_failed")
    require(training and len(held_out) >= 2, "independent_validation_required")
    hashes = lambda samples: {digest(s["coordinates"]) for s in samples}
    require(not hashes(training) & hashes(held_out), "held_out_geometry_leakage")
    require(
        not {s["conformer_identity"] for s in training}
        & {s["conformer_identity"] for s in held_out},
        "held_out_conformer_leakage",
    )
    require(
        criteria["scaled_condition_max"] == 1000
        and criteria["rank_relative_threshold"] == 1e-10,
        "unregistered_fit_criteria",
    )
    require(
        not any(
            conformer_equivalent(t["coordinates"], v["coordinates"])
            for t in training
            for v in held_out
        ),
        "held_out_rigid_or_symmetry_leakage",
    )
    tr = dataset(model, training)
    va = dataset(model, held_out)
    columns = tr[0][1]["columns"]
    ref_s, ref_p = tr[0]
    a = []
    b = []
    for s, p in tr:
        a.append((p["energy_jacobian"] - ref_p["energy_jacobian"])[None, :] / 0.05)
        b.append(
            np.array(
                [
                    (s["energy"] - p["fixed_energy"])
                    - (ref_s["energy"] - ref_p["fixed_energy"])
                ]
            )
            / 0.05
        )
        a.append(p["force_jacobian"].reshape(27, -1) / 0.5)
        b.append((np.array(s["forces"]) - p["fixed_forces"]).reshape(27) / 0.5)
    a = np.vstack(a)
    b = np.concatenate(b)
    scale = np.linalg.norm(a, axis=0)
    require(np.all(scale > 1e-12), "unidentifiable_zero_column")
    singular = np.linalg.svd(a / scale, compute_uv=False)
    rank = int(sum(singular > singular[0] * criteria["rank_relative_threshold"]))
    condition = float(singular[0] / singular[-1])
    require(
        rank == len(columns) and condition < criteria["scaled_condition_max"],
        "joint_model_unidentifiable",
    )
    from scipy.optimize import lsq_linear

    result = lsq_linear(a, b, bounds=(-1000, 1000), tol=1e-12, max_iter=500)
    require(result.success, "fit_optimizer_failed")
    values = result.x
    offset = ref_s["energy"] - ref_p["fixed_energy"] - ref_p["energy_jacobian"] @ values

    def errors(rows):
        es = []
        fs = []
        records = []
        for s, p in rows:
            de = float(
                p["fixed_energy"] + p["energy_jacobian"] @ values + offset - s["energy"]
            )
            df = (
                p["fixed_forces"] + p["force_jacobian"] @ values - np.array(s["forces"])
            )
            es.append(de)
            fs.extend(df.ravel())
            records.append(
                {
                    "conformer_identity": s["conformer_identity"],
                    "geometry_identity": digest(s["coordinates"]),
                    "energy_error": de,
                    "force_errors": df.tolist(),
                    "force_RMS": float(np.sqrt(np.mean(df * df))),
                    "force_max": float(np.max(abs(df))),
                }
            )
        return {
            "residual_structure": records,
            "energy_errors": es,
            "energy_RMSE": float(np.sqrt(np.mean(np.square(es)))),
            "force_RMS": float(np.sqrt(np.mean(np.square(fs)))),
            "force_max": float(np.max(np.abs(fs))),
        }

    train_error = errors(tr)
    validation_error = errors(va)
    covariance = np.linalg.pinv(a.T @ a) * float(
        np.sum((a @ values - b) ** 2) / max(len(b) - len(columns), 1)
    )
    uncertainty = np.sqrt(np.maximum(np.diag(covariance), 0))
    denom = np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    correlation = np.divide(
        covariance, denom, out=np.zeros_like(covariance), where=denom > 0
    )
    return {
        "origin": "fitted_external",
        "scope": "experimental research Si-H disiloxane extension only",
        "coefficients": dict(zip(columns, values.tolist(), strict=True)),
        "rank": rank,
        "condition": condition,
        "training": train_error,
        "untouched_validation": validation_error,
        "uncertainty_1sigma_residual_only": uncertainty.tolist(),
        "covariance": covariance.tolist(),
        "parameter_correlation": correlation.tolist(),
        "energy_offset": float(offset),
        "held_out_thresholds_passed": validation_error["energy_RMSE"]
        <= criteria["validation_energy_RMSE_kcal_mol"]
        and validation_error["force_RMS"]
        <= criteria["validation_force_RMS_kcal_mol_A"],
        "validation_status": "requires_protocol_sensitivity_and_conformer_bootstrap",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def bootstrap_training(model, training, *, replicates=100, seed=4060):
    """Predeclared whole-training-geometry resampling; no held-out input."""
    rows = dataset(model, training)
    require(
        len(rows) == 18 and replicates == 100 and seed == 4060,
        "bootstrap_plan_mutation",
    )
    s0, p0 = rows[0]
    blocks = []
    targets = []
    for s, p in rows:
        blocks.append(
            np.vstack(
                [
                    (p["energy_jacobian"] - p0["energy_jacobian"])[None, :] / 0.05,
                    p["force_jacobian"].reshape(27, -1) / 0.5,
                ]
            )
        )
        targets.append(
            np.concatenate(
                [
                    [
                        (
                            s["energy"]
                            - p["fixed_energy"]
                            - s0["energy"]
                            + p0["fixed_energy"]
                        )
                        / 0.05
                    ],
                    (np.array(s["forces"]) - p["fixed_forces"]).ravel() / 0.5,
                ]
            )
        )
    from scipy.optimize import lsq_linear

    rng = np.random.default_rng(seed)
    values = []
    failed = []
    for k in range(replicates):
        indices = [0, *rng.integers(1, len(rows), len(rows) - 1).tolist()]
        a = np.vstack([blocks[i] for i in indices])
        b = np.concatenate([targets[i] for i in indices])
        norms = np.linalg.norm(a, axis=0)
        if np.any(norms <= 1e-12):
            failed.append({"replicate": k, "reason": "zero_column"})
            continue
        singular = np.linalg.svd(a / norms, compute_uv=False)
        rank = int(sum(singular > singular[0] * 1e-10))
        condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else None
        if rank != len(norms) or condition is None or condition >= 1000:
            failed.append({"replicate": k, "rank": rank, "condition": condition})
            continue
        solution = lsq_linear(a, b, bounds=(-1000, 1000), tol=1e-12, max_iter=500)
        if not solution.success:
            failed.append({"replicate": k, "reason": "optimizer_not_converged"})
            continue
        values.append(solution.x)
    require(len(values) >= 2, "insufficient_valid_bootstrap_replicates")
    v = np.array(values)
    covariance = np.cov(v, rowvar=False)
    denom = np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    correlation = np.divide(
        covariance, denom, out=np.zeros_like(covariance), where=denom > 0
    )
    return {
        "seed": seed,
        "replicates": replicates,
        "valid_replicates": len(values),
        "failed_replicates": failed,
        "columns": rows[0][1]["columns"],
        "coefficient_samples": v.tolist(),
        "sigma": np.std(v, axis=0, ddof=1).tolist(),
        "percentiles_2_5_50_97_5": np.quantile(v, [0.025, 0.5, 0.975], axis=0).tolist(),
        "correlation": correlation.tolist(),
        "held_out_used": False,
        "scope": "training-conformer sampling only; source/model/reference-systematic uncertainty is not estimated by this bootstrap",
    }


def discrepancy_diagnostics(model, training, fit):
    """Training-only residual projections; no new parameter is inferred/adopted.

    Probe source-fixed stretch/bend baseline directions after projecting out
    all20 unknown columns. This is not a held-out-guided second fit.
    """
    from .hamiltonian import Jet, bend, geometry, length

    rows = dataset(model, training)
    s0, p0 = rows[0]
    blocks, targets, probes = [], [], []
    for sample, p in rows:
        blocks.extend(
            [
                (p["energy_jacobian"] - p0["energy_jacobian"])[None, :] / 0.05,
                p["force_jacobian"].reshape(27, -1) / 0.5,
            ]
        )
        targets.extend(
            [
                np.array(
                    [
                        (sample["energy"] - p["fixed_energy"])
                        - (s0["energy"] - p0["fixed_energy"])
                    ]
                )
                / 0.05,
                (np.array(sample["forces"]) - p["fixed_forces"]).ravel() / 0.5,
            ]
        )
        x = geometry(sample["coordinates"])
        features = {}
        for r in model["interaction_inventory"]:
            if r["family"] not in ("quartic_bond", "quartic_angle"):
                continue
            f = r["family"]
            types = tuple(model["typing"][str(i)] for i in r["sites"])
            labels = min(types, types[::-1])
            delta = (
                length(x, r["sites"]) if f == "quartic_bond" else bend(x, r["sites"])
            ) - r["assignment"]["values"][0]
            for power in (1, 2):
                name = f + ":" + ",".join(labels) + f":delta^{power}"
                features.setdefault(name, Jet(0, np.zeros(27)))
                features[name] += delta**power
        probes.append(features)
    a, b = np.vstack(blocks), np.concatenate(targets)
    values = np.array([fit["coefficients"][k] for k in rows[0][1]["columns"]])
    residual = b - a @ values
    orthogonal = residual - a @ np.linalg.lstsq(a, residual, rcond=1e-10)[0]
    results = []
    for name in sorted(probes[0]):
        z = np.concatenate(
            [
                v
                for p in probes
                for v in (
                    np.array([(p[name].v - probes[0][name].v) / 0.05]),
                    -p[name].d / 0.5,
                )
            ]
        )
        outside = z - a @ np.linalg.lstsq(a, z, rcond=1e-10)[0]
        correlation = (
            float(
                outside @ residual / np.linalg.norm(outside) / np.linalg.norm(residual)
            )
            if np.linalg.norm(outside) > 1e-10
            else 0.0
        )
        results.append(
            {
                "probe": name,
                "absolute_residual_cosine": abs(correlation),
                "fraction_of_training_weighted_residual_SS_projected": correlation**2,
                "outside_current_Jacobian_norm": float(np.linalg.norm(outside)),
                "not_a_fitted_source_parameter": True,
            }
        )
    return {
        "training_weighted_residual_SS": float(residual @ residual),
        "fraction_residual_SS_outside_all20_columns": float(
            orthogonal @ orthogonal / (residual @ residual)
        ),
        "active_bound_coefficients": [
            k for k, v in fit["coefficients"].items() if abs(v) >= 1000 - 1e-6
        ],
        "probes": sorted(
            results,
            key=lambda r: r["fraction_of_training_weighted_residual_SS_projected"],
            reverse=True,
        ),
        "probe_meaning": "Unit polynomial energy features sum(delta) or sum(delta^2) for source-fixed type families. Linear probes diagnose equilibrium-force bias; they are outside this Hamiltonian and are not inserted or fitted as PCFF rows. All20 unknown columns are removed before projection.",
        "held_out_used_for_probe_selection": False,
        "additional_coefficients_created": [],
    }
