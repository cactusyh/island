"""Source-baseline curvature diagnosis in verified unknown-null stretch paths."""

import argparse
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed, write
from island.pcff_qm_pilot.hmds_hamiltonian import load_model, parts
from island.pcff_qm_pilot.hmds_independent import scalar
from island.pcff_qm_pilot.hmds_paths import path, response
from scripts.analyze_phase_4j63 import read


def unknown_null_check(model, center, spec, h):
    geometries = [path(model, center, spec, s * h) for s in (0, 1, -1)]
    p = [parts(model, g["coordinates"]) for g in geometries]
    q = [scalar(model, g["coordinates"]) for g in geometries]
    names = p[0]["columns"]
    energy_curvature = (
        p[1]["energy_jacobian"] - 2 * p[0]["energy_jacobian"] + p[2]["energy_jacobian"]
    ) / h**2
    independent = np.array(
        [
            (q[1]["columns"][n] - 2 * q[0]["columns"][n] + q[2]["columns"][n]) / h**2
            for n in names
        ]
    )
    # dE_unknown/ds is constant on a pure single-bond stretch. Check its
    # derivative from all Cartesian-force parameter columns as well.
    d = [
        -np.einsum("iap,ia->p", row["force_jacobian"], geo["tangent"])
        for row, geo in zip(p, geometries, strict=True)
    ]
    gradient_curvature = (d[1] - d[2]) / (2 * h)
    require(
        max(
            np.max(abs(energy_curvature)),
            np.max(abs(independent)),
            np.max(abs(gradient_curvature)),
        )
        < 1e-7,
        "HMDS_stretch_path_not_unknown_null",
    )
    fixed = [
        {"energy": row["fixed_energy"], "forces": row["fixed_forces"]} for row in p
    ]
    return {
        "unknown_columns": names,
        "unknown_energy_curvature_maximum": float(np.max(abs(energy_curvature))),
        "independent_scalar_unknown_curvature_maximum": float(np.max(abs(independent))),
        "unknown_force_chain_rule_curvature_maximum": float(
            np.max(abs(gradient_curvature))
        ),
        "verified_fixed_response": response(*fixed, *geometries, h),
    }


def audit(repo, root):
    model = load_model(repo)
    plan = read(root / "plan/scientific-plan.json")
    reference = read(root / "deformation-reference-results.json")
    x = np.array(read(root / "stationary-center.json")["center"]["coordinates"])
    results = {}
    for spec in plan["paths"]:
        if spec["kind"] != "stretch":
            continue
        rows = {}
        for h in spec["steps"]:
            control = unknown_null_check(model, x, spec, h)
            qm = reference["primary"][spec["name"]]["steps"][str(h)][
                "curvature_from_path_gradient"
            ]
            mm = control["verified_fixed_response"]["curvature_from_path_gradient"]
            rows[str(h)] = {
                **control,
                "QM_total_path_curvature": qm,
                "model_total_path_curvature_independent_of56_unknowns": mm,
                "QM_minus_model_curvature": qm - mm,
                "relative_curvature_difference": abs(qm - mm) / max(abs(qm), 1e-30),
            }
        fine = rows[str(min(spec["steps"]))]
        basis = reference["controls"]["basis"]["paths"][spec["name"]][
            "control_response"
        ]["curvature_from_path_gradient"]
        results[spec["name"]] = {
            "steps": rows,
            "basis_polarization_curvature_shift": basis
            - fine["QM_total_path_curvature"],
            "comparison_scope": "Total fixed molecular curvature includes nonbonded pairs and all source-complete interactions; not an isolated bond constant. Unknown couplings can change first derivatives and other Cartesian responses, but not this controlled path curvature.",
        }
    result = signed(
        {
            "schema": "island.j63-unknown-null-source-baseline.v1",
            "paths": results,
            "units": "kcal/mol/angstrom^2",
            "all56_unknowns_represented_in_null_check": True,
            "mathematical_reason": "Only one graph bond changes, linearly; all graph angles and proper torsions remain fixed. Every selected missing Class-II kernel contains at most one factor of that bond displacement. Its second path derivative is identically zero. Verified both through independent scalar energies and production Cartesian-force parameter derivatives.",
            "physical_coefficients_inferred": False,
            "source_baselines_released": False,
            "uncertainty_scope": "Finite-step and TZVP/TZVPP sensitivity are recorded as sensitivities, not a statistical or complete-basis uncertainty bound.",
        }
    )
    write(root / "source-baseline-path-curvature.json", result)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    a = p.parse_args()
    print(audit(a.repo, a.evidence)["paths"].keys())
