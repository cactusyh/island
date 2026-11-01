"""Block-identifiability, verified physical constraints and targeted geometry design."""

import argparse
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from island.pcff_completion import signed, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_baseline import evaluate, extended_parts
from island.pcff_qm_pilot.hmds_constrained import (
    canonical,
    dual_lower_bound,
    physical_problem,
    verify,
)
from island.pcff_qm_pilot.hmds_design import metrics, pool
from island.pcff_qm_pilot.hmds_paths import equivalent, leakage, path
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j66 import full_evaluator_check, reconstruct


def weak(problem, names, count=3):
    norms = np.maximum(np.linalg.norm(problem["w"], axis=0), 1e-30)
    _, sv, vt = np.linalg.svd(problem["w"] / norms, full_matrices=False)
    rows = []
    for vector, value in zip(vt[-count:], sv[-count:], strict=True):
        direction = vector / norms
        top = np.argsort(abs(vector))[-8:][::-1]
        e = problem["e"] @ direction
        f = (problem["f"] @ direction).reshape(problem["n"], problem["atoms"], 3)
        rows.append(
            {
                "singular_value": float(value),
                "scaled_direction": vector.tolist(),
                "physical_direction": direction.tolist(),
                "top_scaled_parameters": [
                    {
                        "name": names[i],
                        "scaled_component": float(vector[i]),
                        "physical_component": float(direction[i]),
                    }
                    for i in top
                ],
                "energy_sensitivity": e.tolist(),
                "local_force_sensitivity": np.linalg.norm(f, axis=2).tolist(),
                "normalization": "Unit Euclidean direction in column-scaled coefficients. Physical coefficient step=scaled_direction/column_norm; not a confidence interval.",
            }
        )
    return rows


def audit_block_draws(old, current, names):
    old_bound = {r["replicate"] for r in old["attempts"] if r["bound_activity"]}
    old_cond = {
        r["replicate"] for r in old["attempts"] if r["scaled_condition"] > 10000
    }
    per_block = []
    for b in old["blocks"]:
        absent = [r for r in old["attempts"] if b not in r["sampled_blocks"]]
        present = [r for r in old["attempts"] if b in r["sampled_blocks"]]
        per_block.append(
            {
                "block": b,
                "omitted_draws": len(absent),
                "J65_omitted_bound_fraction": sum(
                    bool(r["bound_activity"]) for r in absent
                )
                / max(1, len(absent)),
                "J65_present_bound_fraction": sum(
                    bool(r["bound_activity"]) for r in present
                )
                / max(1, len(present)),
                "J66_omitted_failures": sum(
                    not r["eligible"]
                    for r in current["attempts"]
                    if b not in r["sampled_blocks"]
                ),
                "J66_present_failures": sum(
                    not r["eligible"]
                    for r in current["attempts"]
                    if b in r["sampled_blocks"]
                ),
                "J65_omitted_conditions": [r["scaled_condition"] for r in absent],
                "interpretation": "Correlated block-inclusion associations, not independent causal effects.",
            }
        )
    counts = Counter(n for r in old["attempts"] for n in r["bound_activity"])
    newcounts = Counter(
        names[j]
        for r in current["attempts"]
        if "fit" in r
        for j in r["fit"]["original_unit_verification"]["active_bound_indices"]
    )
    overlaps = Counter(
        "+".join(sorted(r["failure_categories"])) or "eligible"
        for r in current["attempts"]
    )
    return {
        "J65": {
            "total": 64,
            "bound_active": len(old_bound),
            "conditioning_failures": len(old_cond),
            "overlap": sorted(old_bound & old_cond),
            "bound_only": sorted(old_bound - old_cond),
            "conditioning_only": sorted(old_cond - old_bound),
            "eligible": sum(r["eligible"] for r in old["attempts"]),
            "bound_parameter_recurrence": dict(counts),
        },
        "J66_failure_combinations": dict(overlaps),
        "J66_bound_parameter_recurrence": dict(newcounts),
        "per_block": per_block,
    }


def run(repo, root):
    plan, model, observations, data, old = reconstruct(repo, root)
    current = read(root / "block-uncertainty.json")
    if (root / "bootstrap-solves.tar.xz").exists():
        with tempfile.TemporaryDirectory(prefix="j66-diagnosis-records-") as temp:
            directory = Path(temp) / "records"
            unpack(root / "bootstrap-solves.tar.xz", directory)
            current["attempts"] = [
                read(directory / f"replicate-{i:02d}.json") for i in range(64)
            ]
    trial = read(root / "constrained-epigraph.json")
    lo = np.array([s["lower"] for s in plan["parameter_specifications"]])
    hi = np.array([s["upper"] for s in plan["parameter_specifications"]])
    names = data["names"] + ["energy_offset"]
    problem = physical_problem(data)
    audit = audit_block_draws(old, current, names)
    weak_rows = []
    for original in sorted(
        old["attempts"], key=lambda r: r["scaled_condition"], reverse=True
    )[:8]:
        indices = [
            i
            for b in original["sampled_blocks"]
            for i, r in enumerate(observations)
            if r["block"] == b
        ]
        p = physical_problem(data, indices)
        weak_rows.append(
            {
                "replicate": original["replicate"],
                "omitted_blocks": sorted(
                    set(old["blocks"]) - set(original["sampled_blocks"])
                ),
                "condition": original["scaled_condition"],
                "weak": weak(p, names),
                "multiplicities": dict(Counter(original["sampled_blocks"])),
            }
        )
    audit["weak_draws"] = weak_rows
    audit["nominal_weak"] = weak(problem, names)
    certificate = []
    for r in current["attempts"]:
        if "no_verified_feasible_witness" in r["failure_categories"]:
            p = physical_problem(data, r["indices"])
            system = canonical(p, lo, hi)
            bound = dual_lower_bound(system, r["feasibility"]["solution"], p, lo, hi)
            certificate.append(
                {
                    "replicate": r["replicate"],
                    "certificate": bound,
                    "interpretation": "Infeasible sampled-weight problem only"
                    if bound["sufficient_margin"]
                    else "Inconclusive; no independently feasible witness retained. Not physical-model failure.",
                }
            )
    audit["failed_draw_numerical_bounds"] = certificate
    # Verify SCS after an explicitly reported bound projection, without altering limits.
    scs = read(root / "nominal-rho-SCS.json")["solution"]
    v = np.array(scs["values"])
    v[:-1] = np.clip(v[:-1], lo, hi)
    cross = {
        "raw_status": scs["status"],
        "raw_bound_violation": True,
        "bound_projection_maximum": float(np.max(abs(v - np.array(scs["values"])))),
        "projected_original_unit_verification": verify(problem, v, lo, hi),
        "complete_evaluator_check": full_evaluator_check(model, observations, data, v),
    }
    write(root / "independent-solver-cross-check.json", signed(cross))
    physical_checks = {}
    for label in ("nominal-rho", "constrained-epigraph", "interior-rho"):
        values = read(root / (label + ".json"))["solution"]["values"]
        direct = full_evaluator_check(model, observations, data, values)
        physical_checks[label] = direct
    selected = [r for r in current["attempts"] if r["eligible"]]
    for r in selected:
        physical_checks[f"bootstrap-{r['replicate']}"] = full_evaluator_check(
            model, observations, data, r["values"]
        )
    write(root / "complete-evaluator-verification.json", signed(physical_checks))
    groups = {}
    for label, isjoint in (("primary", False), ("joint", True)):
        rows = [
            r
            for r in trial["performance"]["per_configuration"]
            if r["id"].startswith("joint") == isjoint
        ]
        forces = np.array([r["Cartesian_force_errors"] for r in rows])
        groups[label] = {
            "configurations": len(rows),
            "energy_RMSE": float(
                np.sqrt(np.mean([r["energy_error_kcal_mol"] ** 2 for r in rows]))
            ),
            "force_RMS": float(np.sqrt(np.mean(forces**2))),
            "maximum_local_vector_force": float(np.linalg.norm(forces, axis=2).max()),
        }
    oldplan = read(repo / "docs/evidence/phase_4j65/plan/scientific-plan.json")
    calibration = read(repo / "docs/evidence/phase_4j64/baseline-calibration.json")
    params = dict(zip(data["names"], trial["solution"]["values"][:-1], strict=True))
    curvature = {}
    for name in ("SiC", "CH"):
        spec = next(s for s in oldplan["paths"] if s["name"] == name)
        geometries = [
            path(model, oldplan["center_coordinates"], spec, h) for h in (0.003, -0.003)
        ]
        derivatives = [
            -np.sum(evaluate(model, g["coordinates"], params)["forces"] * g["tangent"])
            for g in geometries
        ]
        value = float((derivatives[0] - derivatives[1]) / 0.006)
        reference = calibration["paths"][name]["protocol_responses"][
            "TZVP_dense_tight"
        ]["curvature_from_path_gradient"]
        curvature[name] = {
            "delta_K2": params["delta_K2:" + name],
            "full_model_curvature": value,
            "fixed_plus_twice_correction": calibration["paths"][name]["fixed_curvature"]
            + 2 * params["delta_K2:" + name],
            "QM_curvature": reference,
            "model_minus_QM": value - reference,
            "units": "kcal/mol/angstrom^2",
            "use": "correlated training diagnostic only",
        }
    write(root / "block-identifiability-diagnosis.json", signed(audit))
    write(
        root / "training-comparison.json",
        signed(
            {
                "groups": groups,
                "stretch_curvature": curvature,
                "training_feasibility": "verified original-unit and complete-evaluator witnesses; not model/data/bounds incompatibility",
                "solver_precision": "Raw AlmostSolved status and strict-gate failure retained. Frozen convex numerical selection produces verified witness; not proof of exact LS optimality.",
                "source_or_bounds_changed": False,
            }
        ),
    )
    # Score only geometry features against the measured weak subspace. No QM.
    directions = np.array(
        [r["physical_direction"] for r in audit["nominal_weak"]]
        + [r["weak"][-1]["physical_direction"] for r in weak_rows]
    )
    candidates = []
    training = [r["coordinates"] for r in observations]
    reserved = oldplan["validation_coordinates"]
    for row in pool(model, oldplan["center_coordinates"], 48):
        if any(equivalent(model, row["coordinates"], x) for x in training + reserved):
            continue
        p = extended_parts(model, row["coordinates"])
        e = np.r_[p["energy_jacobian"], 1.0] / (0.1 * np.sqrt(22))
        f = np.c_[p["force_jacobian"].reshape(81, 58), np.zeros(81)] / (
            0.5 * np.sqrt(81 * 22)
        )
        information = np.vstack([e, f]) @ directions.T
        candidates.append((row, information))
    selected = []
    accum = np.zeros((0, len(directions)))
    for _ in range(4):
        scores = [
            (
                float(
                    np.linalg.slogdet(
                        np.eye(len(directions)) * 1e-12
                        + np.vstack([accum, m]).T @ np.vstack([accum, m])
                    )[1]
                ),
                i,
            )
            for i, (_, m) in enumerate(candidates)
        ]
        _, i = max(scores)
        row, m = candidates.pop(i)
        accum = np.vstack([accum, m])
        selected.append(
            {**row, "target_weak_direction_norms": np.linalg.norm(m, axis=0).tolist()}
        )
    leak = leakage(model, training + [r["coordinates"] for r in selected], reserved)
    additions = []
    for row in selected:
        p = extended_parts(model, row["coordinates"])
        additions.extend(
            [
                np.r_[p["energy_jacobian"], 1.0] / (0.1 * np.sqrt(22)),
                *np.c_[p["force_jacobian"].reshape(81, 58), np.zeros(81)]
                / (0.5 * np.sqrt(81 * 22)),
            ]
        )
    write(
        root / "geometry-only-next-design.json",
        signed(
            {
                "selected": selected,
                "maximum_proposed": 4,
                "new_QM": 0,
                "score": "Greedy logdet in measured nominal3 and worst8 draw weak directions; diagnostic ridge1e-12 in design score only, no parameter regularization.",
                "weak_directions": directions.tolist(),
                "before": metrics(problem["w"]),
                "after": metrics(np.vstack([problem["w"], additions])),
                "leakage": leak,
                "validation_still_untouched": True,
                "scope": "Bounded targeted redundancy proposal, not guaranteed cure for uncertainty or model discrepancy; freeze actual future protocol/objective before QM. No basis ladder, new parameter family or random sampling campaign.",
            }
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j66"))
    a = p.parse_args()
    run(a.repo.resolve(), a.evidence.resolve())
