"""Summarize executed J61 calculations and plot the frozen-model comparison."""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, signed, write
from island.pcff_qm_pilot.hamiltonian import load_model
from island.pcff_qm_pilot.stretch_baseline import Configuration, source_baseline
from scripts.phase_4j61 import read


def finish(work, root, repo):
    work, root = Path(work), Path(root)
    ledger = read(work / "ledger.json")
    calls = ledger["calls"]
    optimization = read(work / "optimization/result.json")
    fit = read(root / "baseline-versus-joint-fit.json")
    curvature = read(root / "projected-Hessian.json")
    joint = fit["baseline_adjustable_joint"]
    compute = signed(
        {
            "schema": "island.j61-compute-accounting.v1",
            "absolute_initial_cap": 120,
            "attempted_energy_gradient_evaluations": len(calls),
            "status_counts": dict(Counter(c["status"] for c in calls.values())),
            "group_counts": dict(Counter(k.split("/")[0] for k in calls)),
            "optimization_gradient_calls": len(optimization["calls"]),
            "optimization_wall_seconds": optimization["elapsed_seconds"],
            "sum_native_gradient_wall_seconds": sum(
                c.get("elapsed_seconds", 0) for c in calls.values()
            ),
            "maximum_native_gradient_wall_seconds": max(
                c.get("elapsed_seconds", 0) for c in calls.values()
            ),
            "numerical_Hessian_internal_calls": "Every displaced gradient is a separately reserved call; no hidden engine Hessian.",
            "failed_attempts_preserved": [
                k for k, c in calls.items() if c["status"] != "completed"
            ],
            "ledger_sha256": pin(work / "ledger.json"),
            "plan_sha256": pin(root / "plan/scientific-plan.json"),
            "MD_executed": False,
        }
    )
    write(root / "compute-accounting.json", compute)
    accepted = fit["accepted_coefficients"]
    original_training = read(
        Path(repo) / "docs/evidence/phase_4j60/experimental-joint-fit.json"
    )["training"]
    statuses = {
        "parameter_completeness": "incomplete",
        "QM_execution": "passed",
        "QM_data_integrity": "passed",
        "stationary_candidate": "passed",
        "stretch_curvature_precision": "passed"
        if curvature["reference_curvature_gate_passed"]
        else "failed",
        "fit_readiness": "passed",
        "fit_completion": "passed",
        "fit_validation": "passed" if accepted else "failed",
        "typing_status": "complete",
        "charge_status": "complete",
        "lammps_export_status": "not_run",
        "lammps_runtime_status": "not_run",
        "openmm_build_status": "not_run",
        "cross_engine_equivalence": "not_claimed",
        "scientific_validation": "not_established",
        "production_readiness": "not_established",
    }
    decision = signed(
        {
            "schema": "island.j61-fitting-decision.v1",
            "decision": "bounded_method_specific_candidate_accepted"
            if accepted
            else "baseline_joint_trial_rejected_by_predeclared_gates",
            "trial_parameters": joint["parameters"],
            "accepted_coefficients": accepted,
            "parameter_origin": "fitted_external",
            "rank": joint["rank"],
            "column_L2_scaled_condition": joint["scaled_condition"],
            "active_bounds": joint["active_bounds"],
            "origin": "fitted_external research-extension candidates; not native source rows",
            "fit_completed": True,
            "fit_validated": bool(accepted),
            "gates": fit["gates"],
            "weighted_training_SS_reduction_fraction": fit[
                "training_weighted_SS_reduction_fraction"
            ],
            "historical_18_configuration_training_comparison": {
                "unchanged_J60": {
                    k: original_training[k] for k in ("energy_RMSE", "force_RMS")
                },
                "baseline_joint": {
                    k: joint["J60_training_only"][k]
                    for k in ("energy_RMSE", "force_RMS")
                },
                "energy_units": "kcal/mol",
                "force_units": "kcal/mol/angstrom",
                "inspected_development_data_not_fresh_validation": True,
            },
            "training_energy_RMSE_kcal_mol": joint["training"]["energy_RMSE"],
            "training_force_RMS_kcal_mol_A": joint["training"]["force_RMS"],
            "fresh_validation_energy_RMSE_kcal_mol": joint["fresh_validation"][
                "energy_RMSE"
            ],
            "fresh_validation_force_RMS_kcal_mol_A": joint["fresh_validation"][
                "force_RMS"
            ],
            "scope": fit["scope"],
            "typing_and_charge_scope": "exact ordered nine-atom Si-H control only",
            "parameter_checksum": digest(joint["parameters"]),
            "independent_conformational_basins_validated": False,
            "next_experiment": {
                "record": "plan/next-control-experiment.json",
                "sha256": pin(root / "plan/next-control-experiment.json"),
                "basis": "Training/stationary/curvature diagnostics only; no fresh validation error used to select another parameter family.",
                "additional_QM_executed": False,
            },
            "cause_separation": {
                "reference_precision": "Declared stretch-curvature and representative basis/grid/SCF gates pass; full alternate-protocol Hessians and QM-method systematics are not established.",
                "identifiability": "Local combined24-column rank/conditioning reported separately from empirical bootstrap uncertainty; neither optimizer convergence nor rank implies model adequacy.",
                "fixed_baseline": "Changing stretch r0/K2 materially reduces training discrepancy. Remaining training-only torsional/angular probes and stationary forces outside stretching motivate testing fixed baselines, without proving a unique physical cause.",
                "functional_form": "This fit cannot distinguish every fixed-baseline mismatch from restrictions of the selected functional form. It does not authorize zero couplings, altered pairs/phases or transfer to Si-C.",
            },
            "user_attested_confirmations": 28,
            "original_author_artifact_status": "not_yet_supplied",
            "new_registered_profiles": [],
            "PCFF_MD_executed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            **statuses,
        }
    )
    write(root / "fitting-decision.json", decision)
    directional_comparison(root, repo, fit, curvature)
    predictive_uncertainty(root, repo, fit)
    plot(root, fit, curvature)
    return decision


def directional_comparison(root, repo, fit, curvature):
    """Report the retained full response too; no new post-hoc acceptance gate."""
    model = load_model(repo)
    center = read(root / "stationary-center.json")["center"]
    x = np.array(center["coordinates"])
    q = np.array(read(root / "stretching-subspace.json")["Q"])
    reference = np.array(
        curvature["steps"]["0.002"]["Cartesian_gradient_response_kcal_mol_A2"]
    )
    configurations = [
        Configuration(model, x + (sign * 0.002 * q[:, i]).reshape(9, 3))
        for sign in (1, -1)
        for i in range(8)
    ]
    original = read(Path(repo) / "docs/evidence/phase_4j60/experimental-joint-fit.json")
    vectors = {
        label: np.array([fit[label]["parameters"][k] for k in configurations[0].names])
        for label in (
            "coupling_only_refit_same_constraints",
            "baseline_adjustable_joint",
        )
    }
    vectors["original_J60"] = np.r_[
        source_baseline(model),
        [original["coefficients"][k] for k in configurations[0].columns],
    ]
    rows = {}
    for label, parameters in vectors.items():
        forces = [c.evaluate(parameters)["forces"] for c in configurations]
        response = -np.column_stack(
            [
                (a - b).ravel() / 0.004
                for a, b in zip(forces[:8], forces[8:], strict=True)
            ]
        )
        delta = response - reference
        outside = delta - q @ (q.T @ delta)
        stationary = Configuration(model, x).evaluate(parameters)["forces"].ravel()
        rows[label] = {
            "Cartesian_gradient_response_kcal_mol_A2": response.tolist(),
            "full_response_relative_Frobenius": float(
                np.linalg.norm(delta) / np.linalg.norm(reference)
            ),
            "off_stretch_response_error_norm_kcal_mol_A2": float(
                np.linalg.norm(outside)
            ),
            "stationary_force_stretch_projection_norm": float(
                np.linalg.norm(q.T @ stationary)
            ),
            "stationary_force_complement_projection_norm": float(
                np.linalg.norm(stationary - q @ (q.T @ stationary))
            ),
        }
    write(
        root / "directional-response-comparison.json",
        signed(
            {
                "models": rows,
                "units": "kcal/mol/angstrom^2 for responses; kcal/mol/angstrom for forces",
                "post_fit_diagnostic_not_new_acceptance_gate": True,
                "interpretation": "Orthogonal coordinate projections, not attribution of errors to particular physical terms; no isolated-bond coefficient interpretation.",
            }
        ),
    )


def predictive_uncertainty(root, repo, fit):
    model = load_model(repo)
    bootstrap = read(root / "training-bootstrap.json")
    names = bootstrap["names"]
    parameters = np.array(bootstrap["coefficient_samples"])
    old = Path(repo) / "docs/evidence/phase_4j60"
    anchor_name = read(old / "joint-target-design.json")["selected_training"][0]["id"]
    anchor = read(old / "stages/tight-d-joint-fit/normalized-reference.json")[
        "samples"
    ]["cc5/" + anchor_name]
    configuration = Configuration(model, anchor["coordinates"])
    offsets = [
        anchor["energy"] - configuration.evaluate(p)["energy"] for p in parameters
    ]
    data = read(root / "normalized-reference.json")["samples"]
    rows = []
    for i in range(4):
        key = f"fresh/validation-{i}"
        sample = data[key]
        configuration = Configuration(model, sample["coordinates"])
        predictions = [configuration.evaluate(p) for p in parameters]
        de = np.array(
            [
                v["energy"] + offset - sample["energy"]
                for v, offset in zip(predictions, offsets, strict=True)
            ]
        )
        df = np.array([v["forces"] - np.array(sample["forces"]) for v in predictions])
        rms = np.sqrt(np.mean(df**2, axis=(1, 2)))
        rows.append(
            {
                "key": key,
                "energy_error_resampling_std_kcal_mol": float(np.std(de, ddof=1)),
                "energy_error_resampling_percentiles_2p5_97p5": np.percentile(
                    de, [2.5, 97.5]
                ).tolist(),
                "force_RMS_error_resampling_percentiles_2p5_97p5": np.percentile(
                    rms, [2.5, 97.5]
                ).tolist(),
                "force_component_resampling_std_kcal_mol_A": np.std(
                    df, axis=0, ddof=1
                ).tolist(),
            }
        )
    write(
        root / "uncertainty-propagation.json",
        signed(
            {
                "parameter_names": names,
                "parameter_resampling_percentiles_2p5_97p5": np.percentile(
                    parameters, [2.5, 97.5], axis=0
                ).tolist(),
                "fresh_validation": rows,
                "valid_replicates": len(parameters),
                "training_anchor": anchor_name,
                "validation_offsets_independently_shifted": False,
                "fit_parameter_checksum": digest(
                    fit["baseline_adjustable_joint"]["parameters"]
                ),
                "scope": "Empirical whole-training-configuration resampling dispersion only. Correlated stationary/curvature block fixed. These percentiles exclude source, method/basis, model discrepancy and chemical-transfer systematics; they are not calibrated confidence intervals.",
                "used_to_select_weights_bounds_or_parameters": False,
            }
        ),
    )


def plot(root, fit, curvature):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    for column, split in enumerate(("training", "fresh_validation")):
        for label, color in (
            ("coupling_only_refit_same_constraints", "#b4522c"),
            ("baseline_adjustable_joint", "#176a8b"),
        ):
            rows = fit[label][split]["configurations"]
            x = np.arange(len(rows))
            text = (
                "Couplings only"
                if label.startswith("coupling_only")
                else "Stretch baseline + couplings"
            )
            axes[0, column].plot(
                x,
                [r["energy_error_kcal_mol"] for r in rows],
                "o-",
                color=color,
                label=text,
                markersize=3,
            )
            axes[1, column].plot(
                x,
                [r["force_RMS_kcal_mol_A"] for r in rows],
                "o-",
                color=color,
                label=text,
                markersize=3,
            )
        axes[0, column].axhline(0, color="gray", lw=0.7)
        axes[1, column].axhline(
            0.5, color="gray", lw=0.7, ls="--", label="RMS acceptance threshold"
        )
        axes[0, column].set_title(
            "Training (22 configurations)"
            if column == 0
            else "Fresh deformed validation (4 configurations)"
        )
        axes[0, column].set_ylabel("Relative energy error (kcal/mol)")
        axes[1, column].set_ylabel("Cartesian force RMS error (kcal/mol/Å)")
        for ax in axes[:, column]:
            ax.set_xlabel("Configuration index")
            ax.grid(alpha=0.15)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Nine-atom Si–H control: fixed thresholds, one training energy offset")
    fig.savefig(root / "baseline-versus-joint-fit.png", dpi=160)
    plt.close(fig)
    ref = np.array(curvature["steps"]["0.002"]["projected_symmetric_kcal_mol_A2"])
    diff = np.array(fit["baseline_adjustable_joint"]["projected_curvature"]) - ref
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, matrix, title in zip(
        axes,
        (ref, diff),
        ("Total QM projected Hessian", "Joint-model minus QM"),
        strict=True,
    ):
        image = ax.imshow(
            matrix, cmap="RdBu_r", vmin=-np.max(abs(matrix)), vmax=np.max(abs(matrix))
        )
        ax.set_title(title)
        ax.set_xlabel("Orthonormal stretch direction")
        ax.set_ylabel("Orthonormal stretch direction")
        fig.colorbar(image, ax=ax, label="kcal/mol/Å²")
    fig.savefig(root / "projected-curvature-comparison.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j61")
    )
    args = parser.parse_args()
    print(finish(args.work, args.evidence, args.repo)["decision"])
