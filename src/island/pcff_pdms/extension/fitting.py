"""Staged adapter over existing pilot solver/diagnostics, on supplied matrices.

No second Hamiltonian/evaluator or engine runner is introduced. Existing pilot
infer()/joint_fit() own fixed 58-column/27-atom or 9-atom historical contracts;
this adapter must not mutate those receipts to handle new arbitrary matrices.
Only explicitly synthetic controls may execute without a scientific profile.
"""

import numpy as np

from . import digest, require, signed

STAGES = ("primary", "torsions", "cross_terms", "charges_LJ", "optional_refinement")
LIMITS = {
    "energy_RMSE_kcal_mol": 0.1,
    "force_RMS_kcal_mol_A": 0.5,
    "maximum_energy_error_kcal_mol": 0.3,
    "maximum_local_force_error_kcal_mol_A": 2.0,
    "equilibrium_bond_error_A": 0.02,
    "equilibrium_angle_error_rad": 0.04,
    "stationarity_RMS_kcal_mol_A": 0.05,
    "torsion_RMSE_kcal_mol": 0.1,
    "relative_curvature_error": 0.1,
    "scaled_condition_max": 1000.0,
    "rank_relative_tolerance": 1e-10,
    "uncertainty_fraction_of_bound_width_max": 0.1,
    "parameter_correlation_max": 0.999,
}
UNCERTAINTY = {
    "method": "training_only_Rademacher_target_perturbations",
    "seed": 4078,
    "replicates": 16,
    "noise_scale": 1e-4,
    "scope": "conditional numerical sensitivity only; excludes QM/source/model systematics",
}


def contract(inv):
    order = [t["id"] + "/" + f for t in inv["terms"] for f in t["fields"]]
    return signed(
        {
            "schema": "island.pdms-extension-fitting-contract.v1",
            "parameter_order": order,
            "parameter_bounds": [
                {
                    "parameter": p,
                    "lower": None,
                    "upper": None,
                    "status": "finite_physical_bounds_required_before_real_fit",
                }
                for p in order
            ],
            "limits": LIMITS,
            "uncertainty": UNCERTAINTY,
            "stages": [
                {
                    "stage": s,
                    "freeze_predecessors": list(STAGES[:i]),
                    "required_design": "orthogonal_signed_deformations"
                    if s == "cross_terms"
                    else "independent_reference_blocks",
                    "status": "blocked_no_execution_profile_or_reference_data",
                }
                for i, s in enumerate(STAGES)
            ],
            "solver": "island.pcff_qm_pilot.hmds_fit:lsq_linear/scipy.optimize.lsq_linear(trf,tol=1e-12,max_iter=1000)",
            "diagnostics": "island.pcff_qm_pilot.hmds_design:metrics",
            "energy_offset": "one training-only energy intercept per independent chemical composition; no force offset; frozen for heldout",
            "design_matrix_identity": None,
            "execution_profile": None,
            "joint_all_terms_fit_allowed": False,
            "optional_refinement_allowed": False,
            "zero_fill": False,
            "wildcard_fallback": False,
            "analog_substitution": False,
        }
    )


def check_bounds(parameters, bounds):
    require(
        len(parameters) == len(set(parameters)) and len(parameters) > 0,
        "duplicate_parameter_alias",
    )
    require(
        [b["parameter"] for b in bounds] == parameters, "fitting_bound_order_mutation"
    )
    for b in bounds:
        lo, hi = b["lower"], b["upper"]
        require(
            type(lo) in (float, int)
            and type(hi) in (float, int)
            and np.isfinite([lo, hi]).all()
            and lo < hi,
            "finite_explicit_bounds_required",
        )


def split_check(training, validation, test=()):
    require(training and validation, "separate_training_and_validation_required")
    for samples, role in (
        (training, "training"),
        (validation, "validation"),
        (test, "test"),
    ):
        require(all(s["role"] == role for s in samples), "split_role_mutation")
    all_samples = [*training, *validation, *test]
    for key in ("id", "geometry_identity", "conformer_group"):
        values = [s[key] for s in all_samples]
        require(
            all(values) and len(values) == len(set(values)),
            "training_validation_leakage",
            identity=key,
        )


def synthetic_stage(
    a,
    y,
    validation_a,
    validation_y,
    *,
    parameters,
    bounds,
    metadata,
    stage="primary",
    predecessors=(),
    uncertainty=UNCERTAINTY,
):
    """Manufactured matrix control. NEVER produces an executable coefficient."""
    require(
        stage in STAGES
        and tuple(p["stage"] for p in predecessors) == STAGES[: STAGES.index(stage)],
        "stage_order_or_frozen_predecessor_missing",
    )
    require(all(p["software_gates_passed"] for p in predecessors), "predecessor_failed")
    require(
        metadata["evidence_kind"] == "synthetic_contract_control",
        "scientific_execution_profile_absent",
    )
    check_bounds(parameters, bounds)
    require(
        metadata["parameter_order"] == parameters
        and metadata["bounds_identity"] == digest(bounds),
        "fitting_bound_mutation",
    )
    require(uncertainty == UNCERTAINTY, "uncertainty_protocol_mutation")
    split_check(metadata["training"], metadata["validation"], metadata.get("test", []))
    arrays = [np.asarray(z, float) for z in (a, y, validation_a, validation_y)]
    a, y, va, vy = arrays
    require(
        a.ndim == va.ndim == 2
        and a.shape[1] == va.shape[1] == len(parameters)
        and y.shape == (len(a),)
        and vy.shape == (len(va),)
        and len(a) > len(parameters)
        and all(np.isfinite(z).all() for z in arrays),
        "finite_overdetermined_design_required",
    )
    require(
        metadata["design_identity"] == digest([z.tolist() for z in arrays]),
        "design_matrix_identity_mutation",
    )
    from island.pcff_qm_pilot.hmds_design import metrics

    diagnostics = metrics(a)
    blockers = []
    if diagnostics["rank"] != len(parameters):
        blockers.append("rank_deficient_design")
    if diagnostics["scaled_condition"] >= LIMITS["scaled_condition_max"]:
        blockers.append("ill_conditioned_design")
    if (
        abs(diagnostics["strongest_column_correlation"]["cosine"])
        > LIMITS["parameter_correlation_max"]
    ):
        blockers.append("parameter_correlation_gate_failed")
    result = {
        "stage": stage,
        "diagnostics": diagnostics,
        "matrix_identity": metadata["design_identity"],
        "bounds_identity": digest(bounds),
        "parameter_order": parameters,
        "evidence_kind": "synthetic_contract_control",
        "scientific_qualification": False,
        "accepted_coefficients": 0,
        "executable": False,
        "frozen_predecessor_identities": [p["identity"] for p in predecessors],
    }
    if blockers:
        return signed(
            {
                **result,
                "software_gates_passed": False,
                "blockers": blockers,
                "coefficients": None,
                "uncertainty": None,
            }
        )
    # Same backend used by the established pilot. No MM features are evaluated.
    import scipy

    from island.pcff_qm_pilot.hmds_fit import lsq_linear

    lo = np.array([b["lower"] for b in bounds])
    hi = np.array([b["upper"] for b in bounds])
    norms = np.linalg.norm(a, axis=0)

    def solve(target):
        r = lsq_linear(
            a / norms,
            target,
            bounds=(lo * norms, hi * norms),
            method="trf",
            tol=1e-12,
            max_iter=1000,
            lsmr_tol=1e-12,
        )
        return r, r.x / norms

    raw, values = solve(y)
    tr = float(np.sqrt(np.mean((a @ values - y) ** 2)))
    held = float(np.sqrt(np.mean((va @ values - vy) ** 2)))
    rng = np.random.default_rng(uncertainty["seed"])
    attempts, coefficients = [], []
    for i in range(uncertainty["replicates"]):
        r, v = solve(y + rng.choice([-1, 1], len(y)) * uncertainty["noise_scale"])
        attempts.append(
            {
                "index": i,
                "success": bool(r.success),
                "status": int(r.status),
                "values": v.tolist(),
            }
        )
        coefficients.append(v)
    fraction = float(
        np.max(
            np.quantile(abs(np.array(coefficients) - values), 0.95, axis=0) / (hi - lo)
        )
    )
    if not raw.success or not all(r["success"] for r in attempts):
        blockers.append("solver_failure")
    if tr > LIMITS["energy_RMSE_kcal_mol"]:
        blockers.append("training_error_gate_failed")
    if held > LIMITS["energy_RMSE_kcal_mol"]:
        blockers.append("heldout_error_gate_failed")
    if fraction > LIMITS["uncertainty_fraction_of_bound_width_max"]:
        blockers.append("uncertainty_gate_failed")
    if np.any(np.minimum(values - lo, hi - values) <= 1e-8 * (hi - lo)):
        blockers.append("active_bound_requires_review")
    return signed(
        {
            **result,
            "software_gates_passed": not blockers,
            "blockers": blockers,
            "coefficients": values.tolist(),
            "training_RMSE": tr,
            "heldout_RMSE": held,
            "raw_solver_state": {
                "status": int(raw.status),
                "message": raw.message,
                "success": bool(raw.success),
                "iterations": int(raw.nit or 0),
                "optimality": float(raw.optimality),
                "x_scaled": raw.x.tolist(),
            },
            "solver_identity": {
                "backend": "scipy.optimize.lsq_linear",
                "version": scipy.__version__,
            },
            "uncertainty": {
                "fraction95_of_bound_width": fraction,
                "attempts": attempts,
                **uncertainty,
            },
            "data_used_for_uncertainty": "training_only",
        }
    )
