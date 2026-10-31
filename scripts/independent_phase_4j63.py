"""Numerical, manufactured-only checks of all HMDS Hamiltonian columns."""

import argparse
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, require, signed, write
from island.pcff_qm_pilot.hmds_hamiltonian import evaluate, load_model, parts
from island.pcff_qm_pilot.hmds_independent import energy, scalar
from island.pcff_qm_pilot.hmds_paths import path, response
from scripts.analyze_phase_4j63 import read


def check(repo, root):
    model = load_model(repo)
    plan = read(root / "plan/scientific-plan.json")
    x = np.array(plan["initial_coordinates"])
    rng = np.random.default_rng(4063)
    p = parts(model, x)
    values = dict(zip(p["columns"], rng.uniform(0.02, 0.5, 56), strict=True))
    q = scalar(model, x)
    e = evaluate(model, x, values)
    errors = {}
    for h in (2e-5, 1e-5):
        fd = np.zeros((27, 3))
        jac = np.zeros((27, 3, 56))
        for i, j in np.ndindex(27, 3):
            a, b = x.copy(), x.copy()
            a[i, j] += h
            b[i, j] -= h
            sa, sb = scalar(model, a), scalar(model, b)
            fd[i, j] = -(energy(model, a, values) - energy(model, b, values)) / (2 * h)
            jac[i, j] = [
                -(sa["columns"][n] - sb["columns"][n]) / (2 * h) for n in p["columns"]
            ]
        force_error = float(np.max(abs(fd - e["forces"])))
        jac_error = float(np.max(abs(jac - p["force_jacobian"])))
        require(
            force_error < 3e-5 and jac_error < 3e-7,
            "HMDS_independent_derivative_failed",
        )
        errors[str(h)] = {
            "all81_force_components_maximum_error": force_error,
            "all4536_force_parameter_components_maximum_error": jac_error,
        }
    direction = []
    for spec in plan["paths"]:
        h = min(spec["steps"])
        a, b, c = [path(model, x, spec, s * h) for s in (1, -1, 0)]
        samples = [evaluate(model, d["coordinates"], values) for d in (c, a, b)]
        r = response(*samples, c, a, b, h)
        direction.append({"name": spec["name"], "response": r})
    matrices, targets = [], []
    v = np.array([values[n] for n in p["columns"]])
    for _ in range(6):
        y = x + rng.normal(0, 0.025, x.shape)
        b = parts(model, y)
        a = np.vstack([b["energy_jacobian"], b["force_jacobian"].reshape(81, 56)])
        matrices.append(a)
        # Independent scalar energies plus independent finite-difference forces
        # above establish kernel correctness; manufactured observations here
        # exercise parameter mapping and complete model reconstruction only.
        targets.append(a @ v)
    a = np.vstack(matrices)
    recovered, _, rank, _ = np.linalg.lstsq(a, np.concatenate(targets), rcond=1e-12)
    require(
        rank == 56 and np.max(abs(recovered - v)) < 1e-7,
        "HMDS_manufactured_recovery_failed",
    )
    output = signed(
        {
            "schema": "island.j63-independent-HMDS-kernels.v1",
            "model_identity": model["identity"],
            "geometry_identity": digest(x.tolist()),
            "geometry_scope": "retained J62 starting geometry, software control only",
            "manufactured_coefficients": values,
            "physical_coefficients": False,
            "independent_scalar_component_errors": {
                k: abs(v - q["components"][k]) for k, v in p["components"].items()
            },
            "finite_differences": errors,
            "path_controls": direction,
            "manufactured_recovery_rank": int(rank),
            "manufactured_parameter_maximum_error": float(np.max(abs(recovered - v))),
            "unknown_coefficients_assigned_to_profile": False,
        }
    )
    write(root / "independent-kernel-checks.json", output)
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    a = p.parse_args()
    r = check(a.repo, a.evidence)
    print(r["finite_differences"], r["manufactured_parameter_maximum_error"])
