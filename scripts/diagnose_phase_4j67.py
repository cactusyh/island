"""Localize retained constrained-fit and omitted-block failures without new QM."""

import argparse
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from island.pcff_completion import signed, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_baseline import evaluate
from island.pcff_qm_pilot.hmds_fit import prediction_uncertainty, predictions
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_paths import path
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j67 import report_predictions
from scripts.fit_phase_4j66 import intervals


def run(repo, root):
    plan = read(root / "plan/scientific-plan.json")
    summary = read(root / "block-uncertainty.json")
    outcomes = summary["attempts"]
    specs = plan["parameter_specifications"]
    names = [r["name"] for r in specs] + ["energy_offset"]
    blocks = sorted({b for r in outcomes for b in r["sampled_blocks"]})
    failures = Counter(
        names[i] for r in outcomes for i in r.get("bound_active_indices", [])
    )
    omission = []
    for block in blocks:
        absent = [r for r in outcomes if block not in r["sampled_blocks"]]
        present = [r for r in outcomes if block in r["sampled_blocks"]]

        def scores(rows):
            return {
                "replicates": len(rows),
                "eligible": sum(r["eligible"] for r in rows),
                "failure_categories": dict(
                    Counter(c for r in rows for c in r["failure_categories"])
                ),
                "maximum_condition": max(r["design"]["scaled_condition"] for r in rows)
                if rows
                else None,
                "replicate_ids": [r["replicate"] for r in rows],
            }

        omission.append(
            {
                "block": block,
                "omitted": scores(absent),
                "present": scores(present),
                "interpretation": "Correlated inclusion associations; not an isolated causal experiment.",
            }
        )
    trial = (
        read(root / "constrained-trial.json")
        if (root / "constrained-trial.json").exists()
        else None
    )
    old = read(repo / "docs/evidence/phase_4j66/constrained-epigraph.json")
    oldunc = read(repo / "docs/evidence/phase_4j66/block-uncertainty.json")
    with np.load(root / "physical-arrays.npz", allow_pickle=False) as saved:
        data = {
            k: (saved[k].tolist() if k == "names" else saved[k]) for k in saved.files
        }

    design = read(root / "actual-weighted-design.json")
    worst = []
    for row in outcomes:
        if "postfit_predictions" not in row:
            continue
        full = row["postfit_predictions"]["full26"]
        worst.append(
            {
                "replicate": row["replicate"],
                "eligible": row["eligible"],
                "energy_RMSE": full["energy_RMSE"],
                "force_RMS": full["force_RMS"],
                "maximum_local_force_error": full["maximum_local_force_error"],
            }
        )
    numerical = {
        "qualified_solver_fits": sum(
            r.get("solver_accuracy_qualified", False) for r in outcomes
        ),
        "fit_status_counts": dict(
            Counter(r.get("solver_status", "not_executed") for r in outcomes)
        ),
        "at_most_eligible_ignoring_status": sum(
            not (set(r["failure_categories"]) - {"solver_accuracy"}) for r in outcomes
        ),
    }
    result = {
        "J66": {
            "eligible": oldunc["eligible_replicates"],
            "eligible_fraction": oldunc["eligible_fraction"],
            "nominal": {
                k: old["original_unit_verification"][k]
                for k in (
                    "energy_RMSE",
                    "force_RMS",
                    "maximum_absolute_energy_error",
                    "maximum_local_force_error",
                )
            },
            "condition": old["actual_weighted_design"]["scaled_condition"],
            "failure_counts": oldunc["failure_counts"],
        },
        "J66_frozen_prediction_on_expanded_data": report_predictions(
            data, old["solution"]["values"], plan["configuration_order"]
        ),
        "J67": {
            "eligible": summary["eligible_replicates"],
            "eligible_fraction": summary["eligible_fraction"],
            "condition": design["design"]["scaled_condition"],
            "failure_counts": summary["failure_counts"],
            "nominal_by_subset": trial["predictions"] if trial else None,
        },
        "recurrent_original_bound_parameters": dict(failures.most_common()),
        "block_omission_associations": omission,
        "solver_status_separation": numerical,
        "postfit_full26_extremes": worst,
        "weak_combination": design["weak_direction_effects"],
        "historical_stress_check_is_qualification": False,
        "new_observations_are_independent_conformational_basins": False,
    }
    # Incompatible samples have no fitted vector. Report their retained minimax
    # predictions separately, without including them in fitted uncertainty or
    # pretending that their omitted observations constrained the optimization.
    result["draw_prediction_roles"] = [
        {
            "replicate": row["replicate"],
            "role": "constrained_fit"
            if row["fit_executed"]
            else "incompatible_sample_minimax_only",
        }
        for row in outcomes
    ]
    with tempfile.TemporaryDirectory(prefix="j67-diagnostic-draws-") as temporary:
        directory = Path(temporary) / "draws"
        unpack(root / "bootstrap-solves.tar.xz", directory)
        result["incompatible_draw_minimax_predictions"] = [
            {
                "replicate": row["replicate"],
                "scope": "Incompatible sampled minimax diagnostic; no constrained fit, no qualifying uncertainty contribution.",
                "predictions": report_predictions(
                    data,
                    read(directory / f"replicate-{row['replicate']:02d}-rho.json")[
                        "selected_values"
                    ],
                    plan["configuration_order"],
                ),
            }
            for row in outcomes
            if not row["fit_executed"]
        ]
    if trial:
        values = trial["selected_values"]
        result["stretching_corrections"] = {
            names[-3]: values[-3],
            names[-2]: values[-2],
        }
        model = load_model(repo)
        oldplan = read(repo / "docs/evidence/phase_4j65/plan/scientific-plan.json")
        calibration = read(repo / "docs/evidence/phase_4j64/baseline-calibration.json")
        params = dict(zip(names[:-1], values[:-1], strict=True))
        curvature = {}
        for name in ("SiC", "CH"):
            spec = next(r for r in oldplan["paths"] if r["name"] == name)
            geometries = [
                path(model, oldplan["center_coordinates"], spec, h)
                for h in (0.003, -0.003)
            ]
            derivatives = [
                -np.sum(
                    evaluate(model, g["coordinates"], params)["forces"] * g["tangent"]
                )
                for g in geometries
            ]
            value = float((derivatives[0] - derivatives[1]) / 0.006)
            reference = calibration["paths"][name]["protocol_responses"][
                "TZVP_dense_tight"
            ]["curvature_from_path_gradient"]
            curvature[name] = {
                "delta_K2": params["delta_K2:" + name],
                "total_MM_curvature": value,
                "QM_curvature": reference,
                "MM_minus_QM": value - reference,
                "units": "kcal/mol/angstrom^2",
                "scope": "Correlated training diagnostic, not an independent observation or isolated source constant.",
            }
        result["stretching_curvature_diagnostics"] = curvature
    if trial is None:
        values = read(root / "nominal-rho.json")["selected_values"]
    minimax = read(root / "nominal-rho.json")
    result["minimax_diagnostic_by_subset"] = report_predictions(
        data, minimax["selected_values"], plan["configuration_order"]
    )
    result["minimax_interpretation"] = (
        "Uniform performance-limit feasibility diagnostic; not a constrained least-squares trial or an accepted parameter vector."
    )
    local = np.array(minimax["checks"]["full"]["local_vector_force_errors"])
    order = np.argsort(local.ravel())[-20:][::-1]
    result["limiting_local_residuals"] = [
        {
            "configuration": plan["configuration_order"][i],
            "atom_id": int(j + 1),
            "source_type": load_model(repo)["typing"][str(j + 1)],
            "vector_force_error": float(local[i, j]),
        }
        for i, j in (np.unravel_index(index, local.shape) for index in order)
    ]
    result["minimax_bound_activity"] = [
        {**specs[i], "value": minimax["selected_values"][i]}
        for i in minimax["checks"]["full"]["active_bound_indices"]
    ]
    result["weak_combination"]["top_scaled_parameters"] = [
        {
            "name": names[i],
            "scaled_component": design["weak_direction_effects"]["scaled_combination"][
                i
            ],
        }
        for i in np.argsort(
            np.abs(design["weak_direction_effects"]["scaled_combination"])
        )[-10:][::-1]
    ]
    fitted = [r["values"] for r in outcomes if r["fit_executed"]]
    if len(fitted) >= 2:
        descriptive = intervals(fitted, values, names)
        descriptive["scope"] = (
            "All executed fits, including failed/inaccurate/bound-active replicates. Descriptive sensitivity only, NOT qualifying uncertainty or survivor intervals."
        )
        with np.load(root / "physical-arrays.npz", allow_pickle=False) as arrays:
            data = {
                k: (arrays[k].tolist() if k == "names" else arrays[k])
                for k in arrays.files
            }
        descriptive["prediction_spread"] = prediction_uncertainty(
            data, fitted, values, plan["configuration_order"]
        )
        _e, nominal_f = predictions(data, values)
        changes = np.array([predictions(data, v)[1] - nominal_f for v in fitted])
        descriptive["atom_vector_prediction95"] = np.quantile(
            np.linalg.norm(changes, axis=3), 0.95, axis=0
        ).tolist()
        result["all_fitted_descriptive_sensitivity"] = descriptive
    write(root / "block-identifiability-diagnosis.json", signed(result))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j67"))
    a = p.parse_args()
    run(a.repo.resolve(), a.evidence.resolve())
