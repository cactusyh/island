"""Conditional J66 fits on the exact retained block draws; no native execution."""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed, write
from island.pcff_qm_pilot.hmds_constrained import (
    canonical,
    physical_problem,
    save_system,
    solve,
    verify,
)
from island.pcff_qm_pilot.hmds_design import metrics
from island.pcff_qm_pilot.hmds_fit import performance, prediction_uncertainty
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j66 import full_evaluator_check, reconstruct


def attempt(problem, lo, hi, root, name, kind="rho", interior=0, witness=None):
    require(
        not (root / (name + "-solver.json")).exists(),
        "additive_solve_attempt_identity_required",
    )
    system = canonical(problem, lo, hi, kind=kind, interior=interior)
    save_system(root / (name + "-canonical.npz"), system)
    result = solve(system)
    write(root / (name + "-solver.json"), signed(result))
    checked = verify(problem, result["values"], lo, hi, result["rho"])
    selection = {"rule": "raw solver values", "blend_fraction": 0.0}
    if (
        kind == "least_squares"
        and not checked["original_limits_satisfied"]
        and witness is not None
    ):
        reference = np.array(witness, float)
        reference[:-1] = np.clip(reference[:-1], lo, hi)
        require(
            verify(problem, reference, lo, hi)["original_limits_satisfied"],
            "repair_witness_not_feasible",
        )
        raw = np.array(result["values"])
        for alpha in (1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3):
            candidate = (1 - alpha) * raw + alpha * reference
            repaired = verify(problem, candidate, lo, hi)
            if repaired["original_limits_satisfied"]:
                selection = {
                    "rule": "first original-unit feasible convex blend in frozen sequence",
                    "blend_fraction": alpha,
                    "raw_values": raw.tolist(),
                    "objective_change": repaired["weighted_objective"]
                    - checked["weighted_objective"],
                }
                result = {**result, "values": candidate.tolist()}
                checked = repaired
                break
    return {
        "solution": result,
        "original_unit_verification": checked,
        "numerical_selection": selection,
    }


def intervals(values, nominal, names):
    array = np.array(values)
    limits = np.quantile(array, [0.025, 0.975], axis=0)
    covariance = np.cov(array, rowvar=False)
    sd = np.sqrt(np.maximum(0, np.diag(covariance)))
    corr = np.divide(
        covariance,
        np.outer(sd, sd),
        out=np.zeros_like(covariance),
        where=np.outer(sd, sd) > 0,
    )
    return {
        "parameter_intervals": [
            {
                "name": n,
                "nominal": v,
                "absolute_standard_deviation": s,
                "lower": l,
                "upper": u,
            }
            for n, v, s, l, u in zip(names, nominal, sd, *limits, strict=True)
        ],
        "covariance": covariance.tolist(),
        "correlations": corr.tolist(),
        "scope": "Unqualified survivor-conditional small-block intervals when overall eligibility fails; no model-discrepancy/basis/physical uncertainty.",
    }


def run(repo, root):
    plan, model, observations, data, old = reconstruct(repo, root)
    lo = np.array([s["lower"] for s in plan["parameter_specifications"]])
    hi = np.array([s["upper"] for s in plan["parameter_specifications"]])
    require(
        read(root / "nominal-rho.json")["original_unit_verification"][
            "original_limits_satisfied"
        ],
        "J66_feasible_nominal_required",
    )
    problem = physical_problem(data)
    trial = attempt(
        problem,
        lo,
        hi,
        root,
        "constrained-epigraph",
        kind="least_squares",
        witness=read(root / "nominal-rho.json")["solution"]["values"],
    )
    values = trial["solution"]["values"]
    trial["complete_evaluator_check"] = full_evaluator_check(
        model, observations, data, values
    )
    trial["performance"] = performance(data, values, plan["configuration_order"])
    trial["actual_weighted_design"] = metrics(problem["w"])
    nominal = read(repo / "docs/evidence/phase_4j65/trial-fit.json")["fit"]["values"]
    trial["change_from_J65"] = (np.array(values) - nominal).tolist()
    trial["parameter_names"] = data["names"] + ["energy_offset"]
    write(root / "constrained-epigraph.json", signed(trial))
    interior = attempt(problem, lo, hi, root, "interior-rho", interior=1e-6)
    interior["restriction"] = (
        "Original bounds reduced by original bound_activity_fraction=1e-6 of their widths; no bounds widened."
    )
    write(root / "interior-rho.json", signed(interior))
    require(
        trial["original_unit_verification"]["original_limits_satisfied"],
        "J66_constrained_trial_not_verified",
    )
    attempts = []
    grouped = {
        b: [i for i, r in enumerate(observations) if r["block"] == b]
        for b in old["blocks"]
    }
    directory = root / "bootstrap"
    directory.mkdir(exist_ok=True)
    for original in old["attempts"]:
        number = original["replicate"]
        name = f"replicate-{number:02d}"
        indices = [i for b in original["sampled_blocks"] for i in grouped[b]]
        sample = physical_problem(data, indices)
        design = metrics(sample["w"])
        feasibility = attempt(sample, lo, hi, directory, name + "-rho")
        result = {
            "replicate": number,
            "sampled_blocks": original["sampled_blocks"],
            "indices": indices,
            "feasibility": feasibility,
            "actual_weighted_design": design,
            "eligible": False,
            "values": None,
            "failure_categories": [],
            "J65_eligible": original["eligible"],
        }
        if design["rank"] != 59:
            result["failure_categories"].append("rank")
        if design["scaled_condition"] > 10000:
            result["failure_categories"].append("conditioning")
        if feasibility["original_unit_verification"]["original_limits_satisfied"]:
            fit = attempt(
                sample,
                lo,
                hi,
                directory,
                name + "-fit",
                kind="least_squares",
                witness=feasibility["solution"]["values"],
            )
            result["fit"] = fit
            result["values"] = fit["solution"]["values"]
            check = fit["original_unit_verification"]
            if fit["solution"]["status"] != "Solved":
                result["failure_categories"].append("solver_status")
            if not check["original_limits_satisfied"]:
                result["failure_categories"].append("performance_or_bounds")
            if check["active_bound_indices"]:
                result["failure_categories"].append("active_bounds")
            result["original22_prediction_errors"] = performance(
                data, result["values"], plan["configuration_order"]
            )
        else:
            result["failure_categories"].append("no_verified_feasible_witness")
        result["eligible"] = not result["failure_categories"]
        attempts.append(result)
        write(directory / (name + ".json"), signed(result))
        print(number, result["failure_categories"], flush=True)
    survivors = [r["values"] for r in attempts if r["eligible"]]
    summary = {
        "attempts": attempts,
        "replacement_resamples": 0,
        "all64_original_draws_retained": True,
        "eligible_replicates": len(survivors),
        "eligible_fraction": len(survivors) / 64,
        "failure_counts": dict(
            Counter(c for r in attempts for c in r["failure_categories"])
        ),
        "qualified": False,
        "parameter_intervals": None,
        "prediction_spread": None,
    }
    if len(survivors) >= 2:
        summary.update(intervals(survivors, values, data["names"] + ["energy_offset"]))
        spread = prediction_uncertainty(
            data, survivors, values, plan["configuration_order"]
        )
        summary["prediction_spread"] = spread
        summary["qualified"] = bool(
            len(survivors) / 64 >= 0.8
            and spread["energy_RMS95"] <= 0.1
            and spread["force_RMS95"] <= 0.5
        )
    write(root / "block-uncertainty.json", signed(summary))
    write(
        root / "decision.json",
        signed(
            {
                "training_feasible": True,
                "constrained_trial_verified": True,
                "nominal_bound_active": bool(
                    trial["original_unit_verification"]["active_bound_indices"]
                ),
                "uncertainty_qualified": summary["qualified"],
                "accepted_candidates": [],
                "validation_QM": 0,
                "new_QM": 0,
                "MD_steps": 0,
                "production_validated": False,
                "simulation_readiness": "not_established",
                "automatic_assignment": False,
            }
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j66")
    )
    args = parser.parse_args()
    run(args.repo.resolve(), args.evidence.resolve())
