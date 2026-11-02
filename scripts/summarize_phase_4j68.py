"""Bind the executed J68 comparisons and bounded scientific stopping decision."""

import importlib.metadata
import platform
import sys
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, signed, write
from island.pcff_qm_pilot.hmds_constrained import physical_problem
from island.pcff_qm_pilot.hmds_local_baseline import extend
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j68 import ROOT, reconstruct


def summarize(repo):
    root = repo / ROOT
    _plan, model, rows, data, _lo, _hi = reconstruct(repo)
    specs = read(repo / "docs/evidence/phase_4j67/plan/scientific-plan.json")[
        "parameter_specifications"
    ]
    summary = {}
    violations = {}
    for attempt in read(root / "solver-attempts.json")["attempts"]:
        name = attempt["name"]
        r = read(root / "solves" / (name + ".json"))
        vals = np.array(r["solution"]["values"])
        d = (
            extend(
                model,
                rows,
                data,
                "bridge_angle" if name.startswith("bridge_angle") else "SiO_stretch",
            )
            if name.startswith(("bridge_angle", "SiO_stretch"))
            else data
        )
        p = physical_problem(d)
        m = p["w"] / np.linalg.norm(p["w"], axis=0)
        singular = np.linalg.svd(m, compute_uv=False)
        checks = r["checks"]["subsets"]
        limiting = sorted(
            [
                {"configuration": row["id"], "atom": j + 1, "vector_error": v}
                for row, errors in zip(
                    rows, checks["full26"]["local_vector_force_errors"], strict=True
                )
                for j, v in enumerate(errors)
            ],
            key=lambda x: -x["vector_error"],
        )[:20]
        df = d["force"] @ vals[:-1]
        absolute = np.sum(np.linalg.norm(d["force"] * vals[:-1], axis=2), axis=-1)
        summary[name] = {
            "rho": r["solution"]["rho"],
            "status": r["solution"]["status"],
            "lower": r["certificate"]["lower_bound"],
            "upper": r["certificate"]["primal_upper_bound"],
            "certificate_arithmetic_allowance": r["certificate"][
                "arithmetic_allowance"
            ],
            "classification": r["classification"],
            "performance_limits_satisfied": r["checks"]["performance_limits_satisfied"],
            "original_bounds_satisfied": r["checks"]["original_bounds_satisfied"],
            "subsets": {
                k: {
                    metric: c[metric]
                    for metric in (
                        "energy_RMSE",
                        "force_RMS",
                        "maximum_absolute_energy_error",
                        "maximum_local_force_error",
                        "required_rho",
                    )
                }
                for k, c in checks.items()
            },
            "new4": {
                k: r["predictions"]["new4"][k]
                for k in (
                    "energy_RMSE",
                    "force_RMS",
                    "maximum_absolute_energy_error",
                    "maximum_local_force_error",
                )
            },
            "limiting_local_errors": limiting,
            "parameter_values": vals.tolist(),
            "parameter_names": d["names"] + ["energy_offset"],
            "scaled_condition": float(singular[0] / singular[-1]),
            "rank": int(np.linalg.matrix_rank(m)),
            "scaled_singular_values": singular.tolist(),
            "maximum_sum_adjustable_column_force_norms": float(np.max(absolute)),
            "maximum_norm_of_adjustable_sum": float(np.max(np.linalg.norm(df, axis=2))),
            "corrections": dict(zip(d["names"][-2:], vals[-3:-1], strict=True))
            if len(d["names"]) == 60
            else {},
            "energy_only_intercept_kcal_mol": float(vals[-1]),
            "complete_evaluator": r["complete_evaluator_check"],
            "runtime_seconds": r["solution"]["runtime_seconds"],
            "resources": r["resources"],
        }
        violations[name] = [
            {
                **v,
                "name": specs[v["index"]]["name"],
                "units": specs[v["index"]]["units"],
            }
            for v in r["checks"]["outside_original_bounds"]
        ]
    write(
        root / "executed-comparison.json",
        signed(
            {
                "models": summary,
                "scope": "All added coefficients and outside-original-bound vectors are diagnostic only. Two separate variants, never combined.",
                "solver_call_count": 8,
                "new_QM": 0,
                "MD_steps": 0,
            }
        ),
    )
    write(root / "original-bound-violations.json", signed(violations))
    environment = {
        "python": platform.python_version(),
        "executable": sys.executable,
        "versions": {
            p: importlib.metadata.version(p)
            for p in ("numpy", "scipy", "clarabel", "scs", "cffi", "pycparser")
        },
        "worker_count": 1,
        "BLAS_threads": 1,
        "limits": "4GiB RLIMIT_AS, solver180s, subprocess210s;12-call cap; stopped at8 calls",
        "actual_solver_successes": 8,
        "failed_processes": 1,
        "raw_state_recoveries_without_new_solve": 1,
        "historical_observations_reconstructed": 26,
        "historical_J67_new_gradients_reconstructed": 4,
        "new_native_QM_attempts": 0,
        "new_MD_steps": 0,
        "new_bootstrap_draws": 0,
        "reserved_validation_QM": 0,
    }
    write(root / "execution-environment.json", signed(environment))
    decision = {
        "schema": "island.j68-decision.v1",
        "training_feasible_original_bounded": False,
        "training_performance_feasible_original_bound_free": False,
        "original_span_with_fixed_baseline": "checked_insufficient_for_four_limits",
        "bounds_only_explanation": "rejected; rho drops only0.000143858 when removing all original bounds",
        "hypothesis_results": {
            "bridge_angle": "checked_incompatible",
            "SiO_stretch": "checked_incompatible; strong partial residual reduction",
        },
        "next_step": "Stop this allocation. Audit the source Si-O baseline/equilibrium dependencies and limiting residuals in the existing26 observations before a separately frozen model revision. The tested shared Si-O polynomial is quantitatively relevant but insufficient alone. No automatic combination, further QM, bound revision, validation or parameter release follows.",
        "next_evidence_requirement": "A separately frozen conservative scalar-potential model with inventory-supported derivatives spanning the remaining joint-040/C3, joint-028 and joint-002 local residual directions, and independent source-equilibrium/Class-II dependency evidence. Re-test unchanged protected22/full26 limits on existing data first; more observations cannot remove their established incompatibility.",
        "unique_physical_cause": "not_demonstrated",
        "minimum_bound_relaxation": "not_activated; no bound relaxation can fix original bound-free span incompatibility",
        "candidate_qualified": False,
        "accepted_coefficients": 0,
        "automatic_assignment": False,
        "parameter_completeness": "incomplete",
        "backend_correctness": "research_derivatives_and_complete_evaluator_independently_checked",
        "protocol_specific_target_reproduction": "fails_simultaneous_training_limits",
        "basis_sensitivity": "historical_failed_or_unestablished_qualifications_unchanged",
        "chemical_transfer": "not_established",
        "scientific_validation": "not_established",
        "cross_engine_equivalence": "not_claimed",
        "validation_activated": False,
        "validation_consumed": False,
        "new_QM": 0,
        "validation_QM": 0,
        "MD_steps": 0,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "original_author_artifact_status": "not_yet_supplied",
        "author_confirmations": "28 user-attested",
        "J66_original22_feasibility": "preserved",
        "J67_expanded_incompatibility": "preserved",
        "J67_timing_deviation": "preserved unchanged in parent native evidence and docs",
        "scope": "Neutral-singlet explicit-H27 HMDS existing26 deformations of retained center; no independent basins, native PCFF equivalence, PDMS/capped-DP2 transfer or condensed-phase validation.",
        "comparison_sha256": pin(root / "executed-comparison.json"),
    }
    write(root / "decision.json", signed(decision))
    print(
        {
            k: (s["rho"], s["lower"], s["limiting_local_errors"][:3])
            for k, s in summary.items()
        }
    )


if __name__ == "__main__":
    summarize(Path.cwd())


def remaining_localization(
    model, rows, data, values, comparison, fixed_rows, projection
):
    """Decompose limiting amended-model atoms without testing another potential."""
    from island.pcff_qm_pilot.hmds_local_baseline import jet_correction

    records = []
    for limiting in comparison["limiting_local_errors"]:
        index = next(
            i for i, r in enumerate(rows) if r["id"] == limiting["configuration"]
        )
        atom = limiting["atom"] - 1
        row = rows[index]
        components = {
            k: np.array(v["forces"])[atom]
            for k, v in fixed_rows[index]["families"].items()
        }
        for j, name in enumerate(data["names"]):
            family = "adjustable:" + (
                name.split(":")[0] if not name.startswith("delta_K2") else name
            )
            components.setdefault(family, np.zeros(3))
            components[family] += data["force"][index, atom, :, j] * values[j]
        _, added = jet_correction(model, row["coordinates"], "SiO_stretch")
        components["research_SiO_polynomial"] = added[atom] @ values[58:60]
        total = sum(components.values())
        residual = total - data["QM_force"][index, atom]
        direction = residual / np.linalg.norm(residual)
        records.append(
            {
                **limiting,
                "source_type": model["typing"][str(atom + 1)],
                "residual_vector": residual.tolist(),
                "residual_norm_check": float(np.linalg.norm(residual)),
                "original_weighted_LS_force_residual": projection["checks"]["subsets"][
                    "full26"
                ]["signed_Cartesian_force_errors"][index][atom],
                "sum_component_norms": float(
                    sum(np.linalg.norm(v) for v in components.values())
                ),
                "norm_of_sum": float(np.linalg.norm(total)),
                "families": [
                    {
                        "family": k,
                        "vector": v.tolist(),
                        "projection_on_residual_unit": float(v @ direction),
                    }
                    for k, v in components.items()
                ],
            }
        )
    return {
        "rows": records,
        "scope": "Post-diagnostic localization only. Original fixed/adjustable families plus the same frozen SiO scalar potential; no third hypothesis, fit or new observations.",
    }
