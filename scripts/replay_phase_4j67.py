"""Offline, relocated native/complete-MM/constraint and omission-bootstrap replay."""

import argparse
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from island.pcff_completion import inventory, require
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_augmentation import (
    bootstrap_draws,
    candidate_gates,
    parent_semantics,
    physically_feasible,
    proposal_sensitivities,
    protected_system,
    solver_qualified,
    solver_state_diagnostics,
    validation_activation,
    verified_subsets,
)
from island.pcff_qm_pilot.hmds_constrained import (
    dual_lower_bound,
    physical_problem,
    verify,
)
from island.pcff_qm_pilot.hmds_design import metrics
from island.pcff_qm_pilot.hmds_fit import features, prediction_uncertainty
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j66 import full_evaluator_check
from scripts.analyze_phase_4j67 import assemble, report_predictions, subset_data
from scripts.replay_phase_4j66 import compare_system
from scripts.replay_phase_4j66 import replay as replay_parent


def verify_record(record, problem, protected, lo, hi, kind="rho"):
    closure = solver_state_diagnostics(
        protected_system(problem, lo, hi, kind=kind, protected=protected),
        record["solution"],
    )
    for key, value in closure.items():
        if isinstance(value, float):
            require(
                abs(value - record["raw_state_diagnostics"][key]) < 1e-9,
                "J67_solver_raw_state_replay",
            )
        else:
            require(
                value == record["raw_state_diagnostics"][key],
                "J67_solver_cone_state_replay",
            )
    actual = verified_subsets(
        problem, protected, record["selected_values"], lo, hi, record["solution"]["rho"]
    )
    for subset, check in actual.items():
        previous = record["checks"][subset]
        for key in (
            "energy_RMSE",
            "force_RMS",
            "maximum_absolute_energy_error",
            "maximum_local_force_error",
            "weighted_objective",
        ):
            require(
                abs(previous[key] - check[key]) < 1e-9, "J67_constraint_reconstruction"
            )
        for family in ("slacks", "original_limit_slacks"):
            require(
                all(
                    np.allclose(
                        np.asarray(value),
                        np.asarray(previous[family][name]),
                        atol=1e-9,
                        rtol=0,
                    )
                    for name, value in check[family].items()
                ),
                "J67_physical_slack_replay",
            )
        for key in (
            "signed_energy_errors",
            "signed_Cartesian_force_errors",
            "local_vector_force_errors",
            "parameter_bound_distances",
        ):
            require(
                np.allclose(check[key], previous[key], atol=1e-9, rtol=0),
                "J67_original_residual_array_replay",
            )
        require(
            check["original_limits_satisfied"] == previous["original_limits_satisfied"]
            and check["active_bound_indices"] == previous["active_bound_indices"],
            "J67_physical_classification",
        )
    require(
        record["solver_accuracy_qualified"] == solver_qualified(record["solution"]),
        "J67_accuracy_promoted",
    )
    if kind == "rho":
        raw = np.array(record["solution"]["values"])
        raw[:-1] = np.clip(raw[:-1], lo, hi)
        choices = [problem] + ([protected] if protected is not None else [])
        worst = max(choices, key=lambda p: verify(p, raw, lo, hi)["required_rho"])
        bound = dual_lower_bound(
            protected_system(problem, lo, hi, protected=protected),
            record["solution"],
            worst,
            lo,
            hi,
        )
        require(
            abs(bound["lower_bound"] - record["checked_dual_bound"]["lower_bound"])
            < 1e-10
            and bound["sufficient_margin"]
            == record["checked_dual_bound"]["sufficient_margin"],
            "J67_checked_lower_bound_replay",
        )
    selection = record["selection"]
    if selection["blend_fraction"]:
        raw = np.array(record["solution"]["values"])
        require(selection["raw_values"] == raw.tolist(), "J67_raw_solver_vector_lost")
    else:
        require(
            record["solution"]["values"] == record["selected_values"],
            "J67_unrecorded_selection",
        )
    return actual


def replay(repo, root):
    parent = repo / "docs/evidence/phase_4j66"
    parent_semantics(
        *(read(parent / n) for n in ("receipt.json", "manifest.json", "decision.json"))
    )
    parent_reconstruction = replay_parent(repo, parent)
    require(
        parent_reconstruction["bootstrap_draws"] == 64
        and parent_reconstruction["native_reconstructed"],
        "J67_full_parent_reconstruction",
    )
    with tempfile.TemporaryDirectory(prefix="j67-offline-relocated-") as temp:
        target = Path(temp)
        native = target / "native"
        unpack(root / "native-outputs.tar.xz", native)
        require(
            inventory(native) == read(root / "native-file-inventory.json"),
            "J67_native_archive_mutation",
        )
        plan, model, observations, checks = assemble(repo, root, native)
        require(
            observations == read(root / "fitting-dataset.json")["observations"],
            "J67_observation_producer_mutation",
        )
        require(
            checks
            == {
                k: v
                for k, v in read(root / "independent-native-checks.json").items()
                if k != "identity"
            },
            "J67_independent_native_mutation",
        )
        reconstructed_proposals = proposal_sensitivities(
            model,
            read(parent / "geometry-only-next-design.json"),
            read(parent / "block-identifiability-diagnosis.json"),
        )
        require(
            reconstructed_proposals
            == read(root / "proposal-sensitivity-reconstruction.json")["rows"],
            "J67_proposal_replay",
        )
        data = features(model, observations)
        with np.load(root / "physical-arrays.npz", allow_pickle=False) as saved:
            require(
                saved["names"].tolist() == data["names"]
                and all(
                    np.array_equal(saved[k], v) for k, v in data.items() if k != "names"
                ),
                "J67_physical_matrix_unit_order_mutation",
            )
        lo = np.array([s["lower"] for s in plan["parameter_specifications"]])
        hi = np.array([s["upper"] for s in plan["parameter_specifications"]])
        problem = physical_problem(data)
        protected = physical_problem(data, np.arange(22))
        require(
            plan["parameter_specifications"]
            == read(repo / "docs/evidence/phase_4j65/plan/scientific-plan.json")[
                "parameter_specifications"
            ],
            "J67_bounds_units_order_mutation",
        )
        rho = read(root / "nominal-rho.json")
        compare_system(
            root / "nominal-rho-canonical.npz",
            protected_system(problem, lo, hi, protected=protected),
        )
        verify_record(rho, problem, protected, lo, hi)
        witness_checks = {
            k: full_evaluator_check(
                model,
                [observations[i] for i in index],
                subset_data(data, index),
                rho["selected_values"],
            )
            for k, index in {
                "protected22": np.arange(22),
                "full26": np.arange(26),
            }.items()
        }
        require(
            witness_checks == rho["complete_evaluator_checks"],
            "J67_witness_complete_MM_replay",
        )
        witness_feasible = physically_feasible(rho["checks"]) and all(
            c["original_four_limits_satisfied"] for c in witness_checks.values()
        )

        if (root / "nominal-rho-SCS.json").exists():
            cross = read(root / "nominal-rho-SCS.json")
            compare_system(
                root / "nominal-rho-SCS-canonical.npz",
                protected_system(problem, lo, hi, protected=protected),
            )
            verify_record(cross, problem, protected, lo, hi)
            direct_cross = full_evaluator_check(
                model, observations, data, cross["selected_values"]
            )
            require(
                abs(direct_cross["force_RMS"] - cross["checks"]["full"]["force_RMS"])
                < 1e-9,
                "J67_independent_solver_MM_replay",
            )
        trial = None
        if (root / "constrained-trial.json").exists():
            require(physically_feasible(rho["checks"]), "J67_trial_activation")
            trial = read(root / "constrained-trial.json")
            verify_record(trial, problem, protected, lo, hi, kind="least_squares")
            compare_system(
                root / "constrained-trial-canonical.npz",
                protected_system(
                    problem, lo, hi, kind="least_squares", protected=protected
                ),
            )
            require(
                report_predictions(
                    data, trial["selected_values"], plan["configuration_order"]
                )
                == trial["predictions"],
                "J67_postfit_metric_mutation",
            )
            for key, index in {
                "protected22": np.arange(22),
                "full26": np.arange(26),
                "new4": np.arange(22, 26),
            }.items():
                direct = full_evaluator_check(
                    model,
                    [observations[i] for i in index],
                    subset_data(data, index),
                    trial["selected_values"],
                )
                require(
                    direct == trial["complete_evaluator_checks"][key],
                    "J67_complete_evaluator_disagreement",
                )
        directory = target / "bootstrap"
        unpack(root / "bootstrap-solves.tar.xz", directory)
        require(
            inventory(directory) == read(root / "bootstrap-file-inventory.json"),
            "J67_bootstrap_archive_mutation",
        )
        summary = read(root / "block-uncertainty.json")
        draws = bootstrap_draws(observations)
        require(
            draws == read(root / "plan/expanded-bootstrap-draws.json")["draws"],
            "J67_bootstrap_draw_identity",
        )
        outcomes = []
        for draw, original in zip(draws, summary["attempts"], strict=True):
            name = f"replicate-{draw['replicate']:02d}"
            out = read(directory / (name + "-outcome.json"))
            require(
                {k: v for k, v in out.items() if k != "identity"} == original,
                "J67_failed_outcome_retention",
            )
            require(
                out["sampled_blocks"] == draw["sampled_blocks"],
                "J67_bootstrap_replacement",
            )
            ids = [
                i
                for b in draw["sampled_blocks"]
                for i, r in enumerate(observations)
                if r["block"] == b
            ]
            require(out["indices"] == ids, "J67_omitted_block_contamination")
            sample = physical_problem(data, ids)
            d = metrics(sample["w"])
            feasibility = read(directory / (name + "-rho.json"))
            compare_system(
                directory / (name + "-rho-canonical.npz"),
                protected_system(sample, lo, hi),
            )
            verify_record(feasibility, sample, None, lo, hi)
            failures = []
            if d["rank"] != 59:
                failures.append("rank")
            if d["scaled_condition"] > 10000:
                failures.append("conditioning")
            if out["fit_executed"]:
                require(
                    physically_feasible(feasibility["checks"]),
                    "J67_sample_fit_activation",
                )
                fit = read(directory / (name + "-fit.json"))
                compare_system(
                    directory / (name + "-fit-canonical.npz"),
                    protected_system(sample, lo, hi, kind="least_squares"),
                )
                fc = verify_record(fit, sample, None, lo, hi, kind="least_squares")
                require(
                    out["values"] == fit["selected_values"], "J67_draw_selected_vector"
                )
                if not solver_qualified(fit["solution"]):
                    failures.append("solver_accuracy")
                if not physically_feasible(fc):
                    failures.append("physical_constraints")
                if set(fit["raw_active_original_bounds"]) | set(
                    fc["full"]["active_bound_indices"]
                ):
                    failures.append("original_bound_activity")
                if (
                    not out["complete_evaluator_sample_check"][
                        "original_four_limits_satisfied"
                    ]
                    and "physical_constraints" not in failures
                ):
                    failures.append("physical_constraints")
                require(
                    out["postfit_predictions"]
                    == report_predictions(
                        data, out["values"], plan["configuration_order"]
                    ),
                    "J67_prediction_mutation",
                )
            else:
                require(
                    not physically_feasible(feasibility["checks"]),
                    "J67_feasible_draw_dropped",
                )
                failures.append("no_verified_sample_feasibility")
            require(
                out["failure_categories"] == failures
                and out["eligible"] == (not failures),
                "J67_qualification_failure_erased",
            )
            outcomes.append(out)
        survivors = [r["values"] for r in outcomes if r["eligible"]]
        require(
            len(outcomes) == 64
            and summary["eligible_replicates"] == len(survivors)
            and summary["failure_counts"]
            == dict(Counter(c for r in outcomes for c in r["failure_categories"]))
            and summary["replacement_resamples"] == 0,
            "J67_bootstrap_accounting",
        )
        diagnosis = read(root / "block-identifiability-diagnosis.json")
        require(
            diagnosis["draw_prediction_roles"]
            == [
                {
                    "replicate": row["replicate"],
                    "role": "constrained_fit"
                    if row["fit_executed"]
                    else "incompatible_sample_minimax_only",
                }
                for row in outcomes
            ],
            "J67_draw_prediction_role_mutation",
        )
        incompatible = [row for row in outcomes if not row["fit_executed"]]
        for row, prediction in zip(
            incompatible, diagnosis["incompatible_draw_minimax_predictions"], strict=True
        ):
            minimax = read(directory / f"replicate-{row['replicate']:02d}-rho.json")
            require(
                prediction["replicate"] == row["replicate"]
                and prediction["predictions"]
                == report_predictions(
                    data, minimax["selected_values"], plan["configuration_order"]
                ),
                "J67_incompatible_draw_prediction_mutation",
            )
        qualified = False
        if len(survivors) >= 2 and trial:
            for key, index in {
                "protected22": np.arange(22),
                "full26": np.arange(26),
            }.items():
                s = prediction_uncertainty(
                    subset_data(data, index),
                    survivors,
                    trial["selected_values"],
                    [plan["configuration_order"][i] for i in index],
                )
                require(
                    s == summary["prediction_spread"][key],
                    "J67_prediction_spread_mutation",
                )
            qualified = len(survivors) / 64 >= 0.8 and all(
                s["energy_RMS95"] <= 0.1 and s["force_RMS95"] <= 0.5
                for s in summary["prediction_spread"].values()
            )
        require(summary["qualified"] == qualified, "J67_uncertainty_promoted")
        activation = read(root / "validation-activation.json")
        if trial:
            gates = candidate_gates(trial, summary, metrics(problem["w"]))
            if trial["raw_active_original_bounds"]:
                gates["bounds"] = False
            if not all(
                trial["complete_evaluator_checks"][k]["original_four_limits_satisfied"]
                for k in ("protected22", "full26")
            ):
                gates["training"] = False
            require(gates == activation["gates"], "J67_activation_gate_mutation")
        require(
            activation["activated"] == validation_activation(activation["gates"])
            and not activation["activated"]
            and not activation["validation_consumed"]
            and activation["validation_native_attempts"] == 0,
            "J67_reserved_validation_activation",
        )
        accounting = read(root / "compute-accounting.json")
        require(
            len(read(native / "attempts.json")["calls"])
            == accounting["counted_attempts"]
            == 4
            and accounting["successful_new_gradients"] == 4
            and accounting["validation_attempts"] == 0,
            "J67_native_attempt_accounting",
        )
        return {
            "native_attempts": 4,
            "new_gradients": 4,
            "training_configurations": 26,
            "protected_configurations": 22,
            "bootstrap_draws": 64,
            "eligible": len(survivors),
            "training_feasible": witness_feasible,
            "candidate_qualified": False,
            "validation_untouched": True,
            "offline": True,
            "MD_steps": 0,
        }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j67"))
    a = p.parse_args()
    print(replay(a.repo.resolve(), a.evidence.resolve()))
