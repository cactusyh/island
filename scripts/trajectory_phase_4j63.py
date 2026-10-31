"""Complete-model numerical diagnosis on actual HMDS optimizer development data."""

import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import qr

from island.pcff_completion import require, signed, write
from island.pcff_qm_pilot.hmds_design import matrix, metrics
from island.pcff_qm_pilot.hmds_hamiltonian import load_model, parts
from scripts.analyze_phase_4j63 import read


def diagnose(repo, root):
    model = load_model(repo)
    samples = read(root / "normalized-reference.json")["samples"]
    keys = sorted(k for k in samples if k.startswith("optimization/"))
    require(keys, "HMDS_optimizer_development_data_missing")
    selected = [samples[k] for k in keys]
    protocol = selected[0]["protocol"]
    require(
        all(s["protocol"] == protocol for s in selected),
        "mixed_optimizer_reference_protocols",
    )
    geometries = [np.array(s["coordinates"]) for s in selected]
    a, names = matrix(model, geometries, anchor=geometries[0])
    fixed = [parts(model, x) for x in geometries]
    target = []
    for s, p in zip(selected, fixed, strict=True):
        de = (s["energy"] - selected[0]["energy"]) - (
            p["fixed_energy"] - fixed[0]["fixed_energy"]
        )
        df = np.array(s["forces"]) - p["fixed_forces"]
        target.extend([de / 0.05, *(df.reshape(-1) / 0.5)])
    target = np.array(target)
    scaling = np.maximum(np.linalg.norm(a, axis=0), 1e-30)
    scaled = a / scaling
    u, singular, _ = np.linalg.svd(scaled, full_matrices=False)
    rank = int(sum(singular > singular[0] * 1e-10))
    residual = target - u[:, :rank] @ (u[:, :rank].T @ target)
    qr_check = {
        "available": False,
        "reason": "truncated rank requires a separately resolved null-space convention",
    }
    if rank == len(names):
        q, _, permutation = qr(scaled, pivoting=True, mode="economic")
        other = target - q @ (q.T @ target)
        error = float(np.linalg.norm(other - residual))
        require(
            error < 1e-7 * max(np.linalg.norm(target), 1),
            "HMDS_residual_projection_numerical_disagreement",
        )
        qr_check = {
            "available": True,
            "SVD_pivoted_QR_residual_difference_norm": error,
            "column_pivot_order": permutation.tolist(),
        }
    blocks = residual.reshape(len(keys), 82)
    energy = blocks[:, 0] * 0.05
    force = blocks[:, 1:].reshape(len(keys), 27, 3) * 0.5
    result = signed(
        {
            "schema": "island.j63-actual-optimizer-model-diagnosis.v1",
            "model_identity": model["identity"],
            "scope": "All successfully retained new HMDS optimizer snapshots, one exact protocol; development trajectory, not independent conformers or held-out validation.",
            "samples": keys,
            "protocol": protocol,
            "parameter_columns": names,
            "actual_energy_force_Jacobian": metrics(a),
            "independent_projection_check": qr_check,
            "all56_unknown_contributions_represented": len(names) == 56,
            "omitted_unknown_contributions_assigned_zero": False,
            "optimistic_unconstrained_lower_bound": {
                "weighted_squared_residual": float(residual @ residual),
                "weighted_target_squared_norm": float(target @ target),
                "unexplained_weighted_fraction": float(
                    (residual @ residual) / (target @ target)
                ),
                "relative_energy_RMS_kcal_mol": float(np.sqrt(np.mean(energy**2))),
                "Cartesian_force_RMS_kcal_mol_A": float(np.sqrt(np.mean(force**2))),
                "maximum_atom_force_error_kcal_mol_A": float(
                    np.linalg.norm(force, axis=2).max()
                ),
                "per_configuration": [
                    {
                        "key": key,
                        "relative_energy_residual": float(energy[i]),
                        "Cartesian_force_RMS": float(np.sqrt(np.mean(force[i] ** 2))),
                        "maximum_atom_force_error": float(
                            np.linalg.norm(force[i], axis=1).max()
                        ),
                        "Cartesian_force_residual": force[i].tolist(),
                    }
                    for i, key in enumerate(keys)
                ],
            },
            "energy_anchor": "first optimizer snapshot; one shared relative-energy origin",
            "weights": "frozen diagnostic energy/0.05 kcal/mol and Cartesian force/0.5 kcal/mol/A",
            "interpretation": "Orthogonal residual of QM minus verified fixed MM against the complete unknown space. An optimistic unlimited-coefficient representation bound, not an accepted physical fit, coefficient uncertainty, proof of a uniquely wrong source row, or evidence of global infeasibility. Optimization snapshots are correlated. No baseline release or parameter/weight selection uses this diagnostic.",
            "physical_coefficients_inferred": False,
            "validation_used": False,
            "new_QM_invoked": False,
        }
    )
    write(root / "optimizer-model-diagnosis.json", result)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    a = p.parse_args()
    print(diagnose(a.repo, a.evidence)["optimistic_unconstrained_lower_bound"])
