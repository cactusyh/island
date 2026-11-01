"""Complete HMDS Jacobian design and source-baseline compatibility diagnostics."""

import argparse
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed, write
from island.pcff_qm_pilot.hmds_design import matrix, metrics, screen
from island.pcff_qm_pilot.hmds_hamiltonian import load_model, parts
from island.pcff_qm_pilot.hmds_paths import leakage
from scripts.analyze_phase_4j63 import read


def design(repo, root):
    model = load_model(repo)
    center = read(root / "stationary-center.json")["center"]
    native = read(root / "normalized-reference.json")["samples"]
    # One protocol only. Optimizer trajectories are development data, not
    # stationary Hessian samples or independent conformer validation.
    keys = sorted(k for k in native if k.startswith("paths/"))
    samples = [center, *[native[k] for k in keys]]
    x = np.array(center["coordinates"])
    geometries = [np.array(s["coordinates"]) for s in samples]
    proposal = screen(model, x, geometries)
    reserved = read(root / "plan/reserved-validation.json")["coordinates"]
    chosen = proposal["selected_additional_configurations"]
    leakage_check = leakage(
        model, [*geometries, *[r["coordinates"] for r in chosen]], reserved
    )
    write(root / "joint-design.json", proposal)
    a, names = matrix(model, geometries, anchor=x)
    base = parts(model, x)
    target = []
    details = []
    for key, s, xyz in zip(["center", *keys], samples, geometries, strict=True):
        p = parts(model, xyz)
        de = (s["energy"] - center["energy"]) - (
            p["fixed_energy"] - base["fixed_energy"]
        )
        df = np.array(s["forces"]) - p["fixed_forces"]
        target.extend([de / 0.05, *(df.reshape(-1) / 0.5)])
        details.append(
            {
                "id": key,
                "QM_relative_energy": s["energy"] - center["energy"],
                "verified_fixed_relative_energy": p["fixed_energy"]
                - base["fixed_energy"],
                "QM_minus_fixed_force_RMS": float(np.sqrt(np.mean(df**2))),
            }
        )
    target = np.array(target)
    norms = np.linalg.norm(a, axis=0)
    u, sv, _ = np.linalg.svd(a / np.maximum(norms, 1e-30), full_matrices=False)
    keep = sv > sv[0] * 1e-10
    # No coefficient values are selected/reported. This is the optimistic
    # unconstrained column-space lower bound, not a physical fit with zero
    # nuisance terms. All56 unknown contributions are represented in a.
    orthogonal = target - u[:, keep] @ (u[:, keep].T @ target)
    blocks = orthogonal.reshape(len(samples), 82)
    projected = {
        "energy_RMS_kcal_mol": float(np.sqrt(np.mean((blocks[:, 0] * 0.05) ** 2))),
        "force_RMS_kcal_mol_A": float(np.sqrt(np.mean((blocks[:, 1:] * 0.5) ** 2))),
        "weighted_squared_lower_bound": float(orthogonal @ orthogonal),
        "target_weighted_squared_norm": float(target @ target),
        "unexplained_weighted_fraction": float(
            (orthogonal @ orthogonal) / (target @ target)
        ),
        "per_configuration": [
            {
                "id": key,
                "relative_energy_residual": float(r[0] * 0.05),
                "force_RMS": float(np.sqrt(np.mean((r[1:] * 0.5) ** 2))),
                "maximum_atom_force_error": float(
                    np.linalg.norm(r[1:].reshape(27, 3) * 0.5, axis=1).max()
                ),
            }
            for key, r in zip(["center", *keys], blocks, strict=True)
        ],
    }
    require(len(names) == 56, "HMDS_complete_unknown_design_required")
    diagnostic = signed(
        {
            "schema": "island.j63-source-baseline-compatibility.v1",
            "names": names,
            "sample_keys": ["center", *keys],
            "actual_design": metrics(a),
            "relative_energy_reference": "one common stationary center; no per-validation offset",
            "QM_minus_verified_fixed": details,
            "optimistic_joint_column_space_lower_bound": projected,
            "coefficient_values_inferred": False,
            "all_unknowns_explicitly_represented": True,
            "bounds_or_regularization_applied": False,
            "scope": "Training diagnostic under fixed source baseline and full56-column research model. Projection is not an accepted fit, statistical uncertainty or proof that a particular source row is wrong. Missing coefficients are never assigned zero.",
            "validation_leakage_check": leakage_check,
            "validation_QM_executed": False,
            "validation_used_for_selection": False,
            "physical_parameters_fitted": False,
        }
    )
    write(root / "source-baseline-compatibility.json", diagnostic)
    return proposal, diagnostic


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    a = p.parse_args()
    d, b = design(a.repo, a.evidence)
    print(
        d["current"]["rank"],
        d["current"]["scaled_condition"],
        d["after_selected"]["scaled_condition"],
    )
