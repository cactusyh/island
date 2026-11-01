"""Reconstruct J65 and execute J66 conic problems using retained observations."""

import argparse
import resource
import time
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_baseline import evaluate
from island.pcff_qm_pilot.hmds_constrained import (
    canonical,
    classify,
    physical_problem,
    save_system,
    solve,
    verify,
)
from island.pcff_qm_pilot.hmds_fit import features, infer, performance, weighted_problem
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import read
from scripts.phase_4j66 import PARENTS


def reconstruct(repo, root):
    parent = repo / "docs/evidence/phase_4j65"
    plan = read(root / "plan/scientific-plan.json")
    require(
        all(pin(parent / n) == plan["parent_files"][n] for n in PARENTS),
        "J66_parent_mutation",
    )
    old = read(parent / "plan/scientific-plan.json")
    require(
        plan["parameter_specifications"] == old["parameter_specifications"]
        and plan["protocol"] == old["protocol"]
        and plan["graph_identity"] == old["graph_identity"]
        and plan["atom_order"] == list(range(1, 28))
        and plan["inherited_fit"] == old["fit"],
        "J66_contract_mutation",
    )
    dataset = read(parent / "fitting-dataset.json")
    observations = dataset["observations"]
    require(
        len(observations) == 22
        and [r["id"] for r in observations] == plan["configuration_order"],
        "J66_training_mapping",
    )
    model = load_model(repo)
    data = features(model, observations)
    require(
        [s["name"] for s in plan["parameter_specifications"]] == data["names"],
        "J66_parameter_order",
    )
    matrix, target = weighted_problem(data, old["fit"])
    with np.load(parent / "actual-fitting-system.npz", allow_pickle=False) as stored:
        require(
            np.max(abs(matrix - stored["weighted_Jacobian"])) < 1e-10
            and np.max(abs(target - stored["QM_minus_verified_fixed"])) < 1e-10
            and stored["names"].tolist() == data["names"] + ["energy_offset"],
            "J66_actual_matrix_mutation",
        )
    fit = infer(data, old["parameter_specifications"], old["fit"])
    metrics = performance(data, fit["values"], plan["configuration_order"])
    original = read(parent / "trial-fit.json")
    errors = {
        k: abs(metrics[k] - original["training"][k])
        for k in (
            "energy_RMSE",
            "force_RMS",
            "maximum_absolute_energy_error",
            "maximum_local_force_error",
        )
    }
    require(max(errors.values()) <= 1e-9, "J66_nominal_reconstruction")
    bootstrap = read(parent / "block-uncertainty.json")
    require(len(bootstrap["attempts"]) == 64, "J66_all_draws_required")
    write(
        root / "reconstruction.json",
        signed(
            {
                "parent_pins": plan["parent_files"],
                "errors": errors,
                "metrics": {k: metrics[k] for k in errors},
                "weighted_shape": list(matrix.shape),
                "parameters": data["names"],
                "fixed_MM_reconstructed": True,
                "all58_columns_reconstructed": True,
                "original_producers_preserved": True,
                "graph_identity": plan["graph_identity"],
                "native_replay": "validation/J65-reconstruction.log",
                "bootstrap_draws": bootstrap["attempts"],
                "validation_QM": 0,
                "new_QM": 0,
                "MD_steps": 0,
            }
        ),
    )
    np.savez_compressed(
        root / "physical-arrays.npz",
        **{k: v for k, v in data.items() if k != "names"},
        names=np.array(data["names"]),
    )
    return plan, model, observations, data, bootstrap


def full_evaluator_check(model, observations, data, values):
    energy, force = [], []
    for row in observations:
        actual = evaluate(
            model,
            row["coordinates"],
            dict(zip(data["names"], values[:-1], strict=True)),
        )
        energy.append(actual["energy"] + values[-1])
        force.append(actual["forces"])
    energy, force = np.array(energy), np.array(force)
    expected_e = (
        data["fixed_energy"] + data["energy"] @ np.array(values[:-1]) + values[-1]
    )
    expected_f = data["fixed_force"] + data["force"] @ np.array(values[:-1])
    require(
        np.max(abs(energy - expected_e)) < 1e-9
        and np.max(abs(force - expected_f)) < 1e-9,
        "J66_complete_evaluator_disagreement",
    )
    de, df = energy - data["QM_energy"], force - data["QM_force"]
    return {
        "maximum_absolute_energy_error": float(abs(de).max()),
        "original_four_limits_satisfied": bool(
            np.sqrt(np.mean(de**2)) <= 0.1
            and np.sqrt(np.mean(df**2)) <= 0.5
            and abs(de).max() <= 0.3
            and np.linalg.norm(df, axis=2).max() <= 2
        ),
        "maximum_energy_array_difference": float(np.max(abs(energy - expected_e))),
        "maximum_force_array_difference": float(np.max(abs(force - expected_f))),
        "complete_evaluator_all58_supplied": True,
        "energy_RMSE": float(np.sqrt(np.mean((energy - data["QM_energy"]) ** 2))),
        "force_RMS": float(np.sqrt(np.mean((force - data["QM_force"]) ** 2))),
        "maximum_local_force_error": float(
            np.linalg.norm(force - data["QM_force"], axis=2).max()
        ),
    }


def execute(repo, root):
    require(
        not (root / "nominal-rho-solver.json").exists(),
        "additive_solve_attempt_identity_required",
    )
    plan, model, observations, data, _bootstrap = reconstruct(repo, root)
    lo = np.array([s["lower"] for s in plan["parameter_specifications"]])
    hi = np.array([s["upper"] for s in plan["parameter_specifications"]])
    problem = physical_problem(data)
    for name, options, backend in (
        ("nominal-rho", {}, "clarabel"),
        ("nominal-rho-SCS", {}, "scs"),
        ("unbounded-rho", {"bounded": False}, "clarabel"),
        ("force-only-rho", {"force_only": True}, "clarabel"),
    ):
        start = time.monotonic()
        system = canonical(problem, lo, hi, **options)
        save_system(root / (name + "-canonical.npz"), system)
        answer = solve(system, backend)
        write(root / (name + "-solver.json"), signed(answer))
        checked = verify(problem, answer["values"], lo, hi, answer["rho"])
        record = {
            "solution": answer,
            "original_unit_verification": checked,
            "complete_evaluator_check": full_evaluator_check(
                model, observations, data, answer["values"]
            ),
            "canonical_labels": system["labels"],
            "linear_rows": system["linear"],
            "soc_dimensions": system["soc"],
            "scaling": system["scale"].tolist(),
            "classification": classify(answer, checked),
            "wall_seconds": time.monotonic() - start,
            "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        }
        write(root / (name + ".json"), signed(record))
        print(
            name, answer["status"], answer["rho"], checked["required_rho"], flush=True
        )
    nominal = read(root / "nominal-rho.json")
    if nominal["original_unit_verification"]["original_limits_satisfied"]:
        system = canonical(problem, lo, hi, kind="least_squares")
        save_system(root / "constrained-fit-canonical.npz", system)
        result = solve(system)
        checked = verify(problem, result["values"], lo, hi)
        write(
            root / "constrained-fit.json",
            signed(
                {
                    "solution": result,
                    "original_unit_verification": checked,
                    "complete_evaluator_check": full_evaluator_check(
                        model, observations, data, result["values"]
                    ),
                    "performance": performance(
                        data, result["values"], plan["configuration_order"]
                    ),
                }
            ),
        )
    else:
        write(
            root / "constrained-fit-activation.json",
            signed(
                {
                    "executed": False,
                    "reason": "Original-limit feasible witness not found; certificate analysis separate. Do not execute an ineligible constrained training trial or bootstrap rule.",
                    "reserved_validation_QM": 0,
                }
            ),
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j66"))
    a = p.parse_args()
    execute(a.repo.resolve(), a.evidence.resolve())
