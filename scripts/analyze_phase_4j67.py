"""Execute protected/full constrained analysis of the four actual HMDS observations."""

import argparse
import resource
import shutil
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, inventory, pin, require, signed, write
from island.pcff_qm_pilot.bundle import pack
from island.pcff_qm_pilot.campaign import verify_completed
from island.pcff_qm_pilot.hmds_augmentation import (
    bootstrap_draws,
    candidate_gates,
    exact_intake,
    exact_request,
    physically_feasible,
    protected_system,
    solver_qualified,
    solver_state_diagnostics,
    validation_activation,
    verified_subsets,
)
from island.pcff_qm_pilot.hmds_constrained import (
    dual_lower_bound,
    physical_problem,
    save_system,
    solve,
    verify,
)
from island.pcff_qm_pilot.hmds_design import metrics
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_fit import (
    features,
    performance,
    prediction_uncertainty,
    predictions,
)
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import decimal_energy, read
from scripts.analyze_phase_4j63 import independent_native, sample
from scripts.analyze_phase_4j65 import independent_value, relative
from scripts.analyze_phase_4j66 import full_evaluator_check
from scripts.fit_phase_4j66 import intervals


def assemble(repo, root, work):
    plan = read(root / "plan/scientific-plan.json")
    rows = read(repo / "docs/evidence/phase_4j65/fitting-dataset.json")["observations"]
    model = load_model(repo)
    proposed = read(repo / "docs/evidence/phase_4j66/geometry-only-next-design.json")[
        "selected"
    ]
    # Strip only the block added during freeze; coordinates/actions are never regenerated.
    exact_intake(
        model,
        [
            {k: v for k, v in r.items() if k != "block"}
            for r in plan["training_requests"]
        ],
        proposed,
        rows,
        plan["validation_coordinates"],
    )
    center = rows[0]["native_center_energy_hartree_decimal"]
    checks = {}
    for row in plan["training_requests"]:
        key = "training/" + row["id"]
        r = read(work / "requests" / (key + ".json"))
        exact_request(r, request(plan, row["coordinates"]))
        protocol_keys = (
            "method",
            "basis",
            "grid",
            "thresholds",
            "scf_type",
            "settings",
        )
        require(
            all(r[k] == rows[0]["original_request"][k] for k in protocol_keys),
            "J67_protocol_pooling",
        )
        verify_completed(work / "outputs" / key, digest(r))
        value = sample(work, key)
        result = read(work / "outputs" / key / "result.json")
        require(
            result["basis_functions"] == 399 and result["version"] == "1.9.1",
            "J67_pinned_native_basis_engine",
        )
        checked = {
            **independent_native(work, key, value),
            "independent_decimal_ASCII": independent_value(work, key, value),
        }
        energy = str(decimal_energy(work / "outputs" / key)[0])
        observation = {
            "id": row["id"],
            "block": row["block"],
            "coordinates": row["coordinates"],
            "forces": value["forces"],
            "relative_energy_kcal_mol": relative(energy, center),
            "native_energy_hartree_decimal": energy,
            "native_center_energy_hartree_decimal": center,
            "original_normalized_sample": value,
            "original_request": r,
            "producer": {
                "kind": "new_J67_execution",
                "key": key,
                "producer_result_sha256": pin(work / "outputs" / key / "result.json"),
            },
            "independent_native_check": checked,
        }
        rows.append(observation)
        checks[row["id"]] = checked
    require(
        len(rows) == 26 and [r["id"] for r in rows] == plan["configuration_order"],
        "J67_dataset_mapping",
    )
    return plan, model, rows, checks


def retain(root, work):
    attempts = read(work / "attempts.json")["calls"]
    require(len(attempts) <= 8, "J67_attempt_budget")
    with tempfile.TemporaryDirectory(prefix="j67-native-retention-") as temp:
        target = Path(temp)
        for name in ("requests", "outputs", "processes", "resources"):
            if (work / name).exists():
                shutil.copytree(work / name, target / name)
        for name in (
            "scientific-plan.json",
            "jobs.json",
            "attempts.json",
            "ledger.json",
        ):
            if (work / name).exists():
                shutil.copyfile(work / name, target / name)
        write(root / "native-file-inventory.json", inventory(target))
        pack(target, root / "native-outputs.tar.xz")
    native = read(work / "ledger.json")["calls"]
    write(
        root / "compute-accounting.json",
        signed(
            {
                "counted_attempts": len(attempts),
                "cap": 8,
                "successful_new_gradients": sum(
                    r["status"] == "completed" for r in native.values()
                ),
                "attempts": attempts,
                "native_calls": native,
                "validation_attempts": sum(
                    k.startswith("validation/") for k in attempts
                ),
                "new_MD_steps": 0,
                "new_SiH_QM": 0,
                "new_DP2_QM": 0,
                "geometry_optimizations": 0,
                "Hessians": 0,
            }
        ),
    )


def attempt(
    problem,
    lo,
    hi,
    directory,
    name,
    plan,
    kind="rho",
    protected=None,
    witness=None,
    interior=0,
    backend="clarabel",
):
    file = directory / (name + ".json")
    require(
        not any(
            (directory / (name + suffix)).exists()
            for suffix in (".json", "-solver.json", "-canonical.npz")
        ),
        "J67_additive_solve_required",
    )
    system = protected_system(problem, lo, hi, kind, protected, interior)
    save_system(directory / (name + "-canonical.npz"), system)
    result = solve(system, backend=backend)
    # Always save raw solver states before independent evaluation.
    write(directory / (name + "-solver.json"), signed(result))
    raw = np.array(result["values"])
    checks = verified_subsets(problem, protected, raw, lo, hi, result["rho"])
    selection = {"rule": "raw solver vector", "blend_fraction": 0.0}
    selected = raw.copy()
    if (
        kind == "least_squares"
        and not physically_feasible(checks)
        and witness is not None
    ):
        reference = np.array(witness)
        require(
            physically_feasible(
                verified_subsets(problem, protected, reference, lo, hi)
            ),
            "J67_selection_witness",
        )
        raw_objective = checks["full"]["weighted_objective"]
        for alpha in plan["selection"]["alphas"]:
            candidate = (1 - alpha) * raw + alpha * reference
            new = verified_subsets(problem, protected, candidate, lo, hi)
            if physically_feasible(new):
                selected = candidate
                selection = {
                    "rule": plan["selection"]["rule"],
                    "blend_fraction": alpha,
                    "raw_values": raw.tolist(),
                    "objective_change": new["full"]["weighted_objective"]
                    - raw_objective,
                }
                checks = new
                break
    record = {
        "solution": result,
        "raw_state_diagnostics": solver_state_diagnostics(system, result),
        "selected_values": selected.tolist(),
        "checks": checks,
        "selection": selection,
        "solver_accuracy_qualified": solver_qualified(result),
        "old_status_qualified": result["status"] == "Solved",
        "raw_active_original_bounds": verify(problem, raw, lo, hi)[
            "active_bound_indices"
        ],
        "canonical_labels": system["labels"],
        "soc": system["soc"],
        "linear": system["linear"],
        "scale": system["scale"].tolist(),
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    if kind == "rho":
        # The largest physical ratio supplies a feasible rho upper box for ALL
        # cones; the existing support-certificate implementation includes all rows.
        clipped = raw.copy()
        clipped[:-1] = np.clip(clipped[:-1], lo, hi)
        problems = [problem] + ([protected] if protected is not None else [])
        worst = max(problems, key=lambda p: verify(p, clipped, lo, hi)["required_rho"])
        record["checked_dual_bound"] = dual_lower_bound(system, result, worst, lo, hi)
        record["classification"] = (
            "verified_training_feasible"
            if physically_feasible(checks)
            else "checked_incompatible"
            if record["checked_dual_bound"]["sufficient_margin"]
            else "numerically_inconclusive"
        )
    write(file, signed(record))
    return record


def subset_data(data, ids):
    return {k: (v if k == "names" else v[ids]) for k, v in data.items()}


def report_predictions(data, values, ids):
    return {
        key: performance(subset_data(data, index), values, [ids[i] for i in index])
        for key, index in {
            "protected22": np.arange(22),
            "new4": np.arange(22, 26),
            "full26": np.arange(26),
        }.items()
    }


def weak_effects(data, design, ids):
    vector = np.array(design["weakest_scaled_combination"]) / np.array(
        design["column_norms"]
    )
    de = data["energy"] @ vector[:-1] + vector[-1]
    df = data["force"] @ vector[:-1]
    return {
        "scaled_combination": design["weakest_scaled_combination"],
        "physical_per_unit_scaled_direction": vector.tolist(),
        "unit_effect_per_configuration": [
            {
                "id": name,
                "energy": float(e),
                "maximum_local_force": float(np.linalg.norm(f, axis=1).max()),
            }
            for name, e, f in zip(ids, de, df, strict=True)
        ],
    }


def run(repo, root, work):
    plan, model, observations, native_checks = assemble(repo, root, work)
    write(
        root / "fitting-dataset.json",
        signed(
            {
                "observations": observations,
                "protected_original22": plan["original22_ids"],
                "new4": [r["id"] for r in plan["training_requests"]],
                "original_producers_preserved": True,
                "consistent_protocol": plan["protocol"],
                "energy_offset_count": 1,
            }
        ),
    )
    write(root / "independent-native-checks.json", signed(native_checks))
    retain(root, work)
    data = features(model, observations)
    require(
        data["names"] == [r["name"] for r in plan["parameter_specifications"]],
        "J67_58column_mapping",
    )
    np.savez_compressed(
        root / "physical-arrays.npz",
        **{k: v for k, v in data.items() if k != "names"},
        names=np.array(data["names"]),
    )
    lo = np.array([r["lower"] for r in plan["parameter_specifications"]])
    hi = np.array([r["upper"] for r in plan["parameter_specifications"]])
    problem = physical_problem(data)
    protected = physical_problem(data, np.arange(22))
    design = metrics(problem["w"])
    write(
        root / "actual-weighted-design.json",
        signed(
            {
                "design": design,
                "weak_direction_effects": weak_effects(
                    data, design, plan["configuration_order"]
                ),
            }
        ),
    )
    witness = attempt(problem, lo, hi, root, "nominal-rho", plan, protected=protected)
    witness["complete_evaluator_checks"] = {
        key: full_evaluator_check(
            model,
            [observations[i] for i in index],
            subset_data(data, index),
            witness["selected_values"],
        )
        for key, index in {
            "protected22": np.arange(22),
            "full26": np.arange(26),
        }.items()
    }
    write(root / "nominal-rho.json", signed(witness))
    witness_feasible = physically_feasible(witness["checks"]) and all(
        c["original_four_limits_satisfied"]
        for c in witness["complete_evaluator_checks"].values()
    )
    fit = None
    if witness_feasible:
        fit = attempt(
            problem,
            lo,
            hi,
            root,
            "constrained-trial",
            plan,
            kind="least_squares",
            protected=protected,
            witness=witness["selected_values"],
        )
        fit["predictions"] = report_predictions(
            data, fit["selected_values"], plan["configuration_order"]
        )
        fit["complete_evaluator_checks"] = {
            key: full_evaluator_check(
                model,
                [observations[i] for i in index],
                subset_data(data, index),
                fit["selected_values"],
            )
            for key, index in {
                "protected22": np.arange(22),
                "full26": np.arange(26),
                "new4": np.arange(22, 26),
            }.items()
        }
        fit["change_from_J66"] = (
            np.array(fit["selected_values"])
            - read(repo / "docs/evidence/phase_4j66/constrained-epigraph.json")[
                "solution"
            ]["values"]
        ).tolist()
        write(root / "constrained-trial.json", signed(fit))
    else:
        cross = attempt(
            problem,
            lo,
            hi,
            root,
            "nominal-rho-SCS",
            plan,
            protected=protected,
            backend="scs",
        )
        write(
            root / "constrained-trial-nonactivation.json",
            signed(
                {
                    "executed": False,
                    "reason": witness["classification"],
                    "crosscheck": cross["classification"],
                }
            ),
        )
    draws = read(root / "plan/expanded-bootstrap-draws.json")["draws"]
    require(draws == bootstrap_draws(observations), "J67_frozen_draws_mutation")
    directory = root / "bootstrap"
    directory.mkdir(exist_ok=True)
    outcomes = []
    # Each draw constrains only its sampled rows. Post-fit evaluations NEVER
    # participate in optimization or selection/restoration.
    for draw in draws:
        indices = [
            i
            for b in draw["sampled_blocks"]
            for i, r in enumerate(observations)
            if r["block"] == b
        ]
        sample_problem = physical_problem(data, indices)
        d = metrics(sample_problem["w"])
        name = f"replicate-{draw['replicate']:02d}"
        rho = attempt(sample_problem, lo, hi, directory, name + "-rho", plan)
        out = {
            **draw,
            "indices": indices,
            "design": d,
            "feasibility_classification": rho["classification"],
            "rho": rho["solution"]["rho"],
            "rho_checked_lower_bound": rho["checked_dual_bound"]["lower_bound"],
            "failure_categories": [],
            "eligible": False,
            "fit_executed": False,
            "values": None,
        }
        if d["rank"] != 59:
            out["failure_categories"].append("rank")
        if d["scaled_condition"] > 10000:
            out["failure_categories"].append("conditioning")
        if physically_feasible(rho["checks"]):
            trial = attempt(
                sample_problem,
                lo,
                hi,
                directory,
                name + "-fit",
                plan,
                kind="least_squares",
                witness=rho["selected_values"],
            )
            out.update(
                fit_executed=True,
                values=trial["selected_values"],
                solver_status=trial["solution"]["status"],
                solver_accuracy_qualified=trial["solver_accuracy_qualified"],
                bound_active_indices=sorted(
                    set(trial["raw_active_original_bounds"])
                    | set(trial["checks"]["full"]["active_bound_indices"])
                ),
                original_unit_metrics=trial["checks"]["full"],
                selection=trial["selection"],
            )
            if not trial["solver_accuracy_qualified"]:
                out["failure_categories"].append("solver_accuracy")
            if not physically_feasible(trial["checks"]):
                out["failure_categories"].append("physical_constraints")
            if out["bound_active_indices"]:
                out["failure_categories"].append("original_bound_activity")
            out["complete_evaluator_sample_check"] = full_evaluator_check(
                model,
                [observations[i] for i in indices],
                subset_data(data, np.array(indices)),
                out["values"],
            )
            if (
                not out["complete_evaluator_sample_check"][
                    "original_four_limits_satisfied"
                ]
                and "physical_constraints" not in out["failure_categories"]
            ):
                out["failure_categories"].append("physical_constraints")
            out["postfit_predictions"] = report_predictions(
                data, out["values"], plan["configuration_order"]
            )
        else:
            out["failure_categories"].append("no_verified_sample_feasibility")
        out["eligible"] = not out["failure_categories"]
        write(directory / (name + "-outcome.json"), signed(out))
        outcomes.append(out)
        print(draw["replicate"], out["failure_categories"], flush=True)
    survivors = [r["values"] for r in outcomes if r["eligible"]]
    summary = {
        "attempts": outcomes,
        "replacement_resamples": 0,
        "eligible_replicates": len(survivors),
        "eligible_fraction": len(survivors) / 64,
        "failure_counts": dict(
            Counter(c for r in outcomes for c in r["failure_categories"])
        ),
        "qualified": False,
        "parameter_intervals": None,
        "prediction_spread": None,
        "draws_are_genuine_expanded_block_omission": True,
        "J66_guaranteed_addition_is_bootstrap": False,
    }
    if len(survivors) >= 2 and fit is not None:
        summary.update(
            intervals(
                survivors, fit["selected_values"], data["names"] + ["energy_offset"]
            )
        )
        spreads = {
            key: prediction_uncertainty(
                subset_data(data, index),
                survivors,
                fit["selected_values"],
                [plan["configuration_order"][i] for i in index],
            )
            for key, index in {
                "protected22": np.arange(22),
                "full26": np.arange(26),
            }.items()
        }
        summary["prediction_spread"] = spreads
        summary["qualified"] = bool(
            len(survivors) / 64 >= 0.8
            and all(
                s["energy_RMS95"] <= 0.1 and s["force_RMS95"] <= 0.5
                for s in spreads.values()
            )
        )
        _ep, fp = predictions(data, fit["selected_values"])
        allforces = np.array([predictions(data, v)[1] - fp for v in survivors])
        summary["survivor_conditional_atom_vector_prediction95"] = np.quantile(
            np.linalg.norm(allforces, axis=3), 0.95, axis=0
        ).tolist()
    write(root / "block-uncertainty.json", signed(summary))
    gates = (
        candidate_gates(fit, summary, design)
        if fit is not None
        else {
            "training": False,
            "numerical": False,
            "bounds": False,
            "identifiability": design["rank"] == 59
            and design["scaled_condition"] <= 10000,
            "uncertainty": summary["qualified"],
        }
    )
    if fit is not None and fit["raw_active_original_bounds"]:
        gates["bounds"] = False
    if fit is not None and not all(
        fit["complete_evaluator_checks"][k]["original_four_limits_satisfied"]
        for k in ("protected22", "full26")
    ):
        gates["training"] = False

    active = validation_activation(gates)
    write(
        root / "validation-activation.json",
        signed(
            {
                "activated": active,
                "gates": gates,
                "reserved_configurations": 3,
                "candidate_values_sha256": digest(fit["selected_values"])
                if active
                else None,
                "validation_consumed": False,
                "validation_native_attempts": 0,
                "reason": "All frozen prerequisites pass; explicit request/one-time activation still required."
                if active
                else "At least one frozen prerequisite fails; all three remain untouched.",
            }
        ),
    )
    write(
        root / "decision.json",
        signed(
            {
                "training_feasible": witness_feasible,
                "constrained_trial_executed": fit is not None,
                "candidate_qualified": False,
                "accepted_coefficients": 0,
                "accepted_candidates": [],
                "validation_activated": active,
                "validation_QM": 0,
                "new_QM": 4,
                "MD_steps": 0,
                "parameter_completeness": "incomplete",
                "explicit_trial_vector_complete": fit is not None,
                "uncertainty_qualified": summary["qualified"],
                "automatic_assignment": False,
                "scientific_validation": "not_established",
                "chemical_transfer": "not_established",
                "original_author_artifact_status": "not_yet_supplied",
                "author_confirmations": "28 user-attested",
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
    )
    write(root / "bootstrap-file-inventory.json", inventory(directory))
    pack(directory, root / "bootstrap-solves.tar.xz")
    print("summary", len(survivors), gates, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j67"))
    p.add_argument(
        "--work", type=Path, default=Path("../island-validation/phase4j67-qm")
    )
    a = p.parse_args()
    run(a.repo.resolve(), a.evidence.resolve(), a.work.resolve())
