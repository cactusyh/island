"""Training-only weighting/tradeoff audit; never rewrites J61 conclusions."""

import copy
import json
from pathlib import Path

import numpy as np

from island.pcff_completion import CompletionError, signed, write
from island.pcff_qm_pilot.hamiltonian import load_model
from island.pcff_qm_pilot.stretch_fit import Problem


def build(repo):
    repo = Path(repo)
    old = repo / "docs/evidence/phase_4j60"
    new = repo / "docs/evidence/phase_4j61"
    read = lambda p: json.loads(p.read_text())
    design = read(old / "joint-target-design.json")
    prior = read(old / "stages/tight-d-joint-fit/normalized-reference.json")["samples"]
    fresh = read(new / "normalized-reference.json")["samples"]
    training = [
        {**prior["cc5/" + s["id"]], "key": "J60/" + s["id"]}
        for s in design["selected_training"]
    ] + [
        {**fresh[f"fresh/training-{i}"], "key": f"J61/fresh/training-{i}"}
        for i in range(4)
    ]
    p = Problem(
        load_model(repo),
        training,
        read(new / "stationary-center.json")["center"],
        read(new / "stretching-subspace.json")["Q"],
        read(new / "projected-Hessian.json")["steps"]["0.002"][
            "projected_symmetric_kcal_mol_A2"
        ],
    )
    trial = read(new / "trials/baseline-joint.json")
    values = np.array([trial["parameters"][k] for k in p.names])
    return p, values


def audit(repo, root):
    problem, initial = build(repo)
    r, _j = problem.system(initial)
    rows = []
    blocks = {
        "relative_energies": r[:616].reshape(22, 28)[:, 0],
        "Cartesian_forces": r[:616].reshape(22, 28)[:, 1:],
        "stationarity": r[-63:-36],
        "projected_curvature": r[-36:],
    }
    for factor in (1, 10, 100, 1000):
        p = copy.copy(problem)
        p._last = None
        original = Problem.system.__get__(p)

        def weighted(values, original=original, factor=factor):
            a, b = original(values)
            a, b = a.copy(), b.copy()
            a[-63:-36] *= factor
            b[-63:-36] *= factor
            return a, b

        p.system = weighted
        try:
            fit = p.fit(initial, maximum_evaluations=160)
            v = np.array([fit["parameters"][k] for k in p.names])
            rr, _ = problem.system(v)
            rows.append(
                {
                    "stationarity_penalty_multiplier": factor,
                    "optimizer": fit["optimizer"],
                    "parameters": fit["parameters"],
                    "stationary_force_RMS": fit["stationary_force_RMS"],
                    "curvature_relative_error": fit[
                        "projected_curvature_relative_Frobenius"
                    ],
                    "unmodified_objective_SS": float(rr @ rr),
                    "training": {
                        k: fit["training"][k] for k in ("energy_RMSE", "force_RMS")
                    },
                    "stationarity_and_curvature_feasible": fit["stationary_force_RMS"]
                    <= 0.05
                    and fit["projected_curvature_relative_Frobenius"] <= 0.05,
                }
            )
        except CompletionError as exc:
            rows.append(
                {"stationarity_penalty_multiplier": factor, "failure": str(exc)}
            )
    report = signed(
        {
            "schema": "island.j62-training-tradeoff.v1",
            "J61_soft_optimum_reconstructed_SS": float(r @ r),
            "weighted_SS_blocks": {k: float(np.sum(v * v)) for k, v in blocks.items()},
            "diagnostic_penalty_multipliers": [1, 10, 100, 1000],
            "maximum_evaluations_per_trial": 160,
            "trials": rows,
            "stationarity_and_curvature_feasible_witness_found": any(
                s.get("stationarity_and_curvature_feasible") for s in rows
            ),
            "global_infeasibility_demonstrated": False,
            "held_out_used": False,
            "parameter_or_acceptance_weights_selected_from_validation": False,
            "interpretation": "Diagnostic reweighting is not a change to the J61 scientific fit. A soft optimum failing a constraint is not proof of global infeasibility. Nonlinear bounded searches can provide feasible witnesses but cannot prove nonexistence.",
        }
    )
    write(Path(root) / "training-weighting-tradeoff.json", report)
    return report


if __name__ == "__main__":
    result = audit(Path.cwd(), Path("docs/evidence/phase_4j62"))
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "J61_soft_optimum_reconstructed_SS",
                    "weighted_SS_blocks",
                    "stationarity_and_curvature_feasible_witness_found",
                ]
            }
        )
    )
