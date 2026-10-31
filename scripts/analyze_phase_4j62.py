"""Reconstruct executed controls and compare the preregistered nested models."""

import argparse
import copy
import json
import re
from decimal import Decimal
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed, write
from island.pcff_qm_pilot.angular_fit import Problem
from island.pcff_qm_pilot.molecular_control import response
from island.pcff_qm_pilot.molecular_native import normalized
from scripts.audit_phase_4j62 import build


def read(p):
    return json.loads(Path(p).read_text())


def sample(work, key):
    r = read(work / "requests" / (key + ".json"))
    return {"key": key, **normalized(work / "outputs" / key, r)}


def independent(directory, n):
    text = (directory / "psi4.out").read_text()
    e = Decimal(re.findall(r"@(?:DF-)?RKS Final Energy:\s*([-+\d.]+)", text)[-1])
    lines = (directory / "gradient.dat").read_text().splitlines()
    require(
        lines[1] == "symmetry 0" and int(lines[2]) == n * 3, "independent_native_shape"
    )
    g = np.zeros((n, 3))
    seen = set()
    eh = Decimal("4.3597447222060e-18") * Decimal("6.02214076e23") / Decimal(4184)
    bohr = Decimal("5.29177210544e-11") * Decimal("1e10")
    for line in lines[3:]:
        if not line.strip():
            continue
        i, j, v = line.split()
        i, j = int(i), int(j)
        require((i, j) not in seen, "independent_duplicate_axis")
        seen.add((i, j))
        g[i, j] = float(-Decimal(v) * eh / bohr)
    require(len(seen) == n * 3, "independent_native_coverage")
    return float(e * eh), g


def controls(repo, root, work):
    old = repo / "docs/evidence/phase_4j61"
    plan = read(root / "plan/scientific-plan.json")
    criteria = plan["control"]["precision"]
    keys = read(root / "plan/control-jobs.json")["jobs"]
    samples = {k: sample(work, k) for k in keys}
    qn = np.array(read(root / "nonstretch-projection.json")["Q"])
    qs = np.array(read(old / "stretching-subspace.json")["Q"])
    q = np.column_stack([qs, qn])
    center = plan["control"]["center"]
    oldh = read(old / "projected-Hessian.json")
    steps = {}
    gates = {}
    for h in (0.002, 0.004):
        a = [samples[f"control/h{h:g}-d{d}-s+1"]["forces"] for d in range(4)]
        b = [samples[f"control/h{h:g}-d{d}-s-1"]["forces"] for d in range(4)]
        r = response(a, b, h)
        rs = np.array(
            oldh["steps"][f"{h:g}"]["Cartesian_gradient_response_kcal_mol_A2"]
        )
        full = np.column_stack([rs, r])
        raw = q.T @ full
        sym = (raw + raw.T) / 2
        cross = qs.T @ r - (qn.T @ rs).T
        diag = np.array(
            [
                (
                    samples[f"control/h{h:g}-d{d}-s+1"]["energy"]
                    - 2 * center["energy"]
                    + samples[f"control/h{h:g}-d{d}-s-1"]["energy"]
                )
                / h**2
                for d in range(4)
            ]
        )
        errors = diag - np.diag(qn.T @ r)
        steps[f"{h:g}"] = {
            "Cartesian_response": r.tolist(),
            "combined_raw_projected": raw.tolist(),
            "combined_symmetric_projected": sym.tolist(),
            "nonstretch_diagonal_energy_second_differences": diag.tolist(),
            "energy_gradient_curvature_diagonal_errors": errors.tolist(),
            "stretch_nonstretch_reciprocity_error": cross.tolist(),
            "cross_block_RMS_error": float(np.sqrt(np.mean(cross**2))),
            "combined_eigenvalues": np.linalg.eigvalsh(sym).tolist(),
        }
        gates[f"cross_block_h{h:g}"] = np.sqrt(np.mean(cross**2)) <= criteria[
            "cross_block_absolute_RMS_kcal_mol_A2"
        ] + criteria["cross_block_relative"] * np.sqrt(np.mean((qs.T @ r) ** 2))
        gates[f"nonstretch_symmetry_h{h:g}"] = np.sqrt(
            np.mean(((qn.T @ r) - (qn.T @ r).T) ** 2)
        ) <= criteria["symmetry_absolute_RMS_kcal_mol_A2"] + criteria[
            "symmetry_relative"
        ] * np.sqrt(np.mean((qn.T @ r) ** 2))
        gates[f"energy_curvature_h{h:g}"] = bool(
            np.all(
                abs(errors)
                <= criteria["energy_curvature_absolute_kcal_mol_A2"]
                + criteria["energy_curvature_relative"] * abs(np.diag(qn.T @ r))
            )
        )
    fine = np.array(steps["0.002"]["Cartesian_response"])
    coarse = np.array(steps["0.004"]["Cartesian_response"])
    diff = np.sqrt(np.mean((fine - coarse) ** 2, axis=0))
    scale = np.sqrt(np.mean(fine**2, axis=0))
    gates["two_step_columns"] = bool(
        np.all(
            diff
            <= criteria["two_step_absolute_column_RMS_kcal_mol_A2"]
            + criteria["two_step_relative"] * scale
        )
    )
    sensitivity = {}
    for name, spec in plan["control"]["sensitivity"].items():
        directions = spec["directions"]
        a = [samples[f"sensitivity/{name}-d{d}-s+1"]["forces"] for d in directions]
        b = [samples[f"sensitivity/{name}-d{d}-s-1"]["forces"] for d in directions]
        r = response(a, b, 0.002)
        err = np.sqrt(np.mean((r - fine[:, directions]) ** 2, axis=0))
        ref = scale[directions]
        limit = (
            criteria[name + "_absolute_column_RMS_kcal_mol_A2"]
            + criteria[name + "_relative"] * ref
        )
        gates[name + "_columns"] = bool(np.all(err <= limit))
        sensitivity[name] = {
            "directions": directions,
            "Cartesian_response": r.tolist(),
            "column_RMS_errors": err.tolist(),
            "allowed_column_RMS_errors": limit.tolist(),
            "coverage": "Representative subset only; not a full alternate-protocol Hessian.",
        }
    checks = []
    for key, s in samples.items():
        e, f = independent(work / "outputs" / key, 9)
        de = abs(e - s["energy"])
        df = float(np.max(abs(f - s["forces"])))
        require(de < 2e-7 and df < 1e-10, "independent_native_conversion")
        checks.append(
            {"key": key, "energy_error_kcal_mol": de, "force_max_error_kcal_mol_A": df}
        )
    report = signed(
        {
            "schema": "island.j62-nonstretch-control.v1",
            "steps": steps,
            "two_step_column_RMS_errors": diff.tolist(),
            "fine_column_RMS_scales": scale.tolist(),
            "sensitivity": sensitivity,
            "gates": {k: bool(v) for k, v in gates.items()},
            "reference_precision_passed": bool(all(gates.values())),
            "independent_native_checks": checks,
            "Q_combined": q.tolist(),
            "units": "kcal/mol/angstrom^2; unweighted normalized Cartesian directions",
            "minimum_claim": "Twelve sampled directions only; no verified full-molecular minimum or isolated bond/angle/torsion coefficient claim.",
        }
    )
    write(root / "angular-torsional-control.json", report)
    return report


def nested(repo, root, report):
    base, initial = build(repo)
    fine = report["steps"]["0.002"]
    problem = Problem(
        base,
        report["Q_combined"],
        fine["combined_symmetric_projected"],
        fine["Cartesian_response"],
    )
    values = np.r_[initial, np.zeros(4)]
    rows = {}
    for name, count in [("M0", 24), ("M1", 25), ("M2", 27), ("M3", 28)]:
        result = problem.fit(values, count)
        result["reference_precision_passed"] = report["reference_precision_passed"]
        write(root / "nested" / f"{name}.json", signed(result))
        rows[name] = result
        values = np.array([result["parameters"][k] for k in problem.names])
        print(
            json.dumps(
                {
                    "model": name,
                    "training": {
                        k: result["training"][k] for k in ["energy_RMSE", "force_RMS"]
                    },
                    "stationary_force_RMS": result["stationary_force_RMS"],
                    "curvature": result["projected_curvature_relative_error"],
                    "gates": result["training_gates"],
                }
            ),
            flush=True,
        )
    qualifying = [
        k
        for k, v in rows.items()
        if v["training_qualifies"] and report["reference_precision_passed"]
    ]
    chosen = qualifying[0] if qualifying else None
    uncertainty = {
        "status": "not_run_no_training_qualifier",
        "scope": "Every model reports absolute residual-covariance sigma proxies, including near-zero coefficients; these are not scientific confidence intervals.",
    }
    if chosen:
        result = rows[chosen]
        values = np.array([result["parameters"][k] for k in problem.names])
        rng = np.random.default_rng(4062)
        samples = []
        failures = []
        for iteration in range(100):
            indices = [
                0,
                *rng.integers(
                    1, len(problem.training), len(problem.training) - 1
                ).tolist(),
            ]
            p = copy.copy(problem)
            p.training = [problem.training[i] for i in indices]
            p.configurations = [problem.configurations[i] for i in indices]
            p._last = None
            trial = p.fit(values, result["parameter_count"], maximum_evaluations=100)
            if (
                trial["optimizer"]["success"]
                and trial["rank"] == result["parameter_count"]
                and trial["column_scaled_condition"] < 1000
            ):
                samples.append(
                    [
                        trial["parameters"][k]
                        for k in problem.names[: result["parameter_count"]]
                    ]
                )
            else:
                failures.append(iteration)
        sigma = np.std(samples, axis=0, ddof=1)
        v = values[: result["parameter_count"]]
        weak = [k for k, s, x in zip(problem.names, sigma, v) if s > 0.5 * abs(x)]
        uncertainty = {
            "status": "passed" if len(samples) >= 90 and not weak else "failed",
            "samples": samples,
            "sigma_absolute": sigma.tolist(),
            "names": problem.names[: result["parameter_count"]],
            "valid_replicates": len(samples),
            "failed_replicates": failures,
            "weak_relative_sigma_parameters": weak,
            "correlated_stationary_and_curvature_block_fixed": True,
            "validation_used": False,
            "method_source_and_model_systematics_included": False,
        }
    write(root / "nested-uncertainty.json", signed(uncertainty))
    selection = signed(
        {
            "models": list(rows),
            "chosen_for_conditional_fresh_validation": chosen
            if uncertainty["status"] == "passed"
            else None,
            "training_qualifiers": qualifying,
            "validation_used_for_model_selection": False,
            "decision": "Si-H stops at this bounded milestone; no automatic further model expansion.",
            "accepted_coefficients": {},
            "scientific_validation": "not_established",
            "parameter_completeness": "incomplete",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "nested-model-selection.json", selection)
    return selection


def hmds_report(root, work):
    center = read(root / "HMDS-center.json")
    last = center["optimization_result"]["calls"][-1]
    keys = [
        k
        for k, v in read(work / "ledger.json")["calls"].items()
        if k.startswith("hmds") and v["status"] == "completed"
    ]
    samples = {k: sample(work, k) for k in keys}
    checks = []
    for k, s in samples.items():
        e, f = independent(work / "outputs" / k, 27)
        require(
            abs(e - s["energy"]) < 2e-7 and np.max(abs(f - s["forces"])) < 1e-10,
            "HMDS_independent_conversion",
        )
        checks.append(
            {
                "key": k,
                "energy_error": abs(e - s["energy"]),
                "force_max_error": float(np.max(abs(f - s["forces"]))),
            }
        )
    q = np.array(center["directions"])
    h = 0.003
    rows = []
    for d in range(3):
        a, b = samples[f"hmds-response/d{d}-s+1"], samples[f"hmds-response/d{d}-s-1"]
        r = response([a["forces"]], [b["forces"]], h)[:, 0]
        rows.append(
            {
                "direction": d,
                "Cartesian_response": r.tolist(),
                "projected_diagonal_kcal_mol_A2": float(q[:, d] @ r),
                "energy_second_difference_kcal_mol_A2": (
                    a["energy"] - 2 * samples[last]["energy"] + b["energy"]
                )
                / h**2,
            }
        )
    a, b = samples["hmds-reference/tzvp-d0-s+1"], samples["hmds-reference/tzvp-d0-s-1"]
    r = response([a["forces"]], [b["forces"]], h)[:, 0]
    svp = np.array(rows[0]["Cartesian_response"])
    sensitivity = {
        "center_gradient_RMS_difference_kcal_mol_A": float(
            np.sqrt(
                np.mean(
                    (
                        np.array(samples["hmds-reference/tzvp-center"]["forces"])
                        - samples[last]["forces"]
                    )
                    ** 2
                )
            )
        ),
        "relative_SiC_deformation_energy_differences_kcal_mol": [
            (s["energy"] - samples["hmds-reference/tzvp-center"]["energy"])
            - (
                samples[f"hmds-response/d0-s{sign:+d}"]["energy"]
                - samples[last]["energy"]
            )
            for s, sign in [(a, 1), (b, -1)]
        ],
        "SiC_response_RMS_difference_kcal_mol_A2": float(
            np.sqrt(np.mean((r - svp) ** 2))
        ),
        "SiC_response_relative_difference": float(
            np.linalg.norm(r - svp) / np.linalg.norm(r)
        ),
        "TZVP_Cartesian_response": r.tolist(),
    }
    write(
        root / "HMDS-reference-results.json",
        signed(
            {
                "center_stationary_candidate": center["stationary_candidate"],
                "responses": rows,
                "SVP_TZVP_sensitivity": sensitivity,
                "independent_native_checks": checks,
                "accuracy_scope": "Exploratory SVP/TZVP contrast at the same geometry. One-step responses, no grid/step convergence, final QM accuracy or classical coefficient inference established.",
                "MM_fitting_performed": False,
                "minimum_verified": False,
                "methyl_typing_and_charge_source": "Pinned J54 PCFF-only candidate audit; no hsi values transferred.",
                "next_dependency": "Verify remaining methyl Class-II source coefficients and coupled dependency model. Then obtain converged Si-C/C-H and methyl/Si-O torsional reference responses at one justified method, with independent deformations and held-out geometries. Interior capped-DP2 is a separate environment.",
            }
        ),
    )
    return samples


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j62"))
    p.add_argument("--control-only", action="store_true")
    a = p.parse_args()
    r = controls(a.repo, a.evidence, a.work)
    print(json.dumps({"reference_gates": r["gates"]}), flush=True)
    nested(a.repo, a.evidence, r)
    if not a.control_only:
        hmds_report(a.evidence, a.work)
