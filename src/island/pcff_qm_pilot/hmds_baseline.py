"""Two explicit quadratic research corrections, never replacement source rows."""

import numpy as np

from island.pcff_completion import require

from .hmds_hamiltonian import evaluate as original_evaluate
from .hmds_hamiltonian import parts, validate
from .residual import bond

IDENTITY = "island_research_class2_hmds_j64_quadratic_v1"
NAMES = ("delta_K2:SiC", "delta_K2:CH")
UNITS = "kcal/mol/angstrom^2"


def corrections(model, xyz):
    """E_corr = delta_K2 * sum (r-r0)^2, with no extra 1/2 prefactor."""
    validate(model)
    x = np.asarray(xyz, float)
    require(x.shape == (27, 3) and np.isfinite(x).all(), "HMDS_correction_mapping")
    energy = np.zeros(2)
    forces = np.zeros((27, 3, 2))
    rows = []
    for row in model["interaction_inventory"]:
        if row["family"] != "quartic_bond":
            continue
        sites = row["sites"]
        types = sorted(model["typing"][str(i)] for i in sites)
        column = 0 if types == ["c3", "sio"] else 1 if types == ["c3", "hc"] else None
        if column is None:
            continue
        require(row["assignment"] is not None, "source_bond_equilibrium_required")
        r, dr = bond(x, sites)
        r0 = row["assignment"]["values"][0]
        energy[column] += (r - r0) ** 2
        forces[:, :, column] -= 2 * (r - r0) * dr
        rows.append({"id": row["id"], "column": NAMES[column], "r0_A": r0})
    require(len(rows) == 24, "HMDS_quadratic_source_domain")
    return {
        "columns": list(NAMES),
        "energy_jacobian": energy,
        "force_jacobian": forces,
        "source_dependencies": rows,
    }


def extended_parts(model, xyz):
    p = parts(model, xyz)
    b = corrections(model, xyz)
    return {
        **p,
        "columns": p["columns"] + b["columns"],
        "energy_jacobian": np.r_[p["energy_jacobian"], b["energy_jacobian"]],
        "force_jacobian": np.concatenate(
            [p["force_jacobian"], b["force_jacobian"]], axis=-1
        ),
    }


def evaluate(model, xyz, coefficients=None):
    p = extended_parts(model, xyz)
    require(
        coefficients is not None and set(coefficients) == set(p["columns"]),
        "explicit_all58_HMDS_coefficients_required",
    )
    require(
        np.isfinite([coefficients[n] for n in p["columns"]]).all(),
        "nonfinite_HMDS_coefficient",
    )
    old = original_evaluate(model, xyz, {n: coefficients[n] for n in p["columns"][:-2]})
    c = corrections(model, xyz)
    v = np.array([coefficients[n] for n in NAMES])
    return {
        "energy": float(old["energy"] + c["energy_jacobian"] @ v),
        "forces": old["forces"] + c["force_jacobian"] @ v,
        "Hamiltonian": IDENTITY,
        "parameter_origin": "explicit_research_inputs",
    }


def calibrate(qm_curvature, fixed_curvature, sensitivity):
    s = np.asarray(sensitivity, float)
    q, f = np.asarray(qm_curvature, float), np.asarray(fixed_curvature, float)
    require(
        s.shape == (2, 2)
        and q.shape == f.shape == (2,)
        and np.isfinite(s).all()
        and np.isfinite(q).all()
        and np.isfinite(f).all(),
        "two_finite_unknown_null_observables_required",
    )
    require(np.linalg.matrix_rank(s) == 2, "baseline_curvature_not_identifiable")
    return np.linalg.solve(s, q - f)


def objective_status(
    native_success, stretching_precision, torsion_precision, calibration_supplied
):
    return {
        "native_calculation_success": bool(native_success),
        "stretching_reference_precision": bool(stretching_precision),
        "soft_torsion_reference_precision": bool(torsion_precision),
        "baseline_calibration_completed": bool(
            stretching_precision and calibration_supplied
        ),
        "candidate_validation": "not_established",
        "complete_force_fit_authorized": False,
        "parameter_completeness": "incomplete",
    }
