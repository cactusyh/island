"""Numerical J61 curvature gates and separately constrained research refits."""

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from island.pcff_completion import CompletionError, digest, pin, require, signed, write
from island.pcff_qm_pilot.hamiltonian import load_model
from island.pcff_qm_pilot.stretch_baseline import (
    directional_curvature,
    source_baseline,
)
from island.pcff_qm_pilot.stretch_fit import (
    Problem,
    leakage,
    metrics,
    remaining_discrepancy,
)
from scripts.independent_phase_4j59 import native_values
from scripts.phase_4j61 import imported, read, stationary


def curvature_report(work, evidence):
    work, evidence = Path(work), Path(evidence)
    plan = read(evidence / "plan/scientific-plan.json")
    verification = read(evidence / "plan/verification-criteria.json")
    center = stationary(work, evidence)["center"]
    q = np.asarray(read(evidence / "stretching-subspace.json")["Q"])
    rows = {}
    native_checks = []
    keys = read(evidence / "plan/gradient-jobs.json")["jobs"]
    samples = {key: imported(work, key) for key in keys}
    # Independent reconstruction reads native text/ASCII matrices, not importer arrays.
    for key in [*keys, *read(work / "optimization/calls.json")]:
        sample = samples.get(key) or imported(work, key)
        energy, force = native_values(work / "outputs" / key)
        de = abs(energy - sample["energy"])
        df = float(np.max(abs(force - np.asarray(sample["forces"]))))
        require(
            de < verification["independent_native_energy_error_kcal_mol"]
            and df < verification["independent_native_force_error_kcal_mol_A"],
            "J61_independent_native_units_failed",
        )
        native_checks.append(
            {"key": key, "energy_error_kcal_mol": de, "force_max_error_kcal_mol_A": df}
        )
    for h in plan["projection"]["steps_angstrom"]:
        plus = [samples[f"curvature/h{h:g}-d{d}-s+1"]["forces"] for d in range(8)]
        minus = [samples[f"curvature/h{h:g}-d{d}-s-1"]["forces"] for d in range(8)]
        response, raw, sym = directional_curvature(plus, minus, q, h)
        diagonal = np.array(
            [
                (
                    samples[f"curvature/h{h:g}-d{d}-s+1"]["energy"]
                    - 2 * center["energy"]
                    + samples[f"curvature/h{h:g}-d{d}-s-1"]["energy"]
                )
                / h**2
                for d in range(8)
            ]
        )
        # Reconstruct the response from independently read native ASCII gradients.
        independent_response = np.column_stack(
            [
                -(
                    native_values(work / "outputs" / f"curvature/h{h:g}-d{d}-s+1")[1]
                    - native_values(work / "outputs" / f"curvature/h{h:g}-d{d}-s-1")[1]
                ).ravel()
                / (2 * h)
                for d in range(8)
            ]
        )
        rows[f"{h:g}"] = {
            "step_angstrom": h,
            "Cartesian_gradient_response_kcal_mol_A2": response.tolist(),
            "projected_raw_kcal_mol_A2": raw.tolist(),
            "projected_symmetric_kcal_mol_A2": sym.tolist(),
            "projected_eigenvalues_kcal_mol_A2": np.linalg.eigvalsh(sym).tolist(),
            "antisymmetric_relative_Frobenius": float(
                np.linalg.norm(raw - raw.T) / np.linalg.norm(sym)
            ),
            "independent_response_max_error": float(
                np.max(abs(independent_response - response))
            ),
            "energy_second_difference_diagonal": diagonal.tolist(),
            "energy_gradient_curvature_diagonal_errors": (
                diagonal - np.diag(sym)
            ).tolist(),
        }
    fine = np.asarray(rows["0.002"]["projected_symmetric_kcal_mol_A2"])
    coarse = np.asarray(rows["0.004"]["projected_symmetric_kcal_mol_A2"])
    relative_steps = float(np.linalg.norm(fine - coarse) / np.linalg.norm(fine))
    sensitivities = {}
    for name in ("basis", "grid", "scf"):
        directions = plan["projection"]["sensitivity_directions"]
        a = np.column_stack(
            [
                -(
                    np.asarray(samples[f"sensitivity/{name}-d{d}-s+1"]["forces"])
                    - np.asarray(samples[f"sensitivity/{name}-d{d}-s-1"]["forces"])
                ).ravel()
                / 0.004
                for d in directions
            ]
        )
        b = np.asarray(rows["0.002"]["Cartesian_gradient_response_kcal_mol_A2"])[
            :, directions
        ]
        delta = a - b
        sensitivities[name] = {
            "directions": directions,
            "Cartesian_response": a.tolist(),
            "difference": delta.tolist(),
            "Cartesian_relative_Frobenius": float(
                np.linalg.norm(delta) / np.linalg.norm(b)
            ),
            "projected_relative_Frobenius": float(
                np.linalg.norm(q.T @ delta) / np.linalg.norm(q.T @ b)
            ),
            "center_force_RMS_difference": float(
                np.sqrt(
                    np.mean(
                        (
                            np.asarray(samples[f"sensitivity/{name}-center"]["forces"])
                            - center["forces"]
                        )
                        ** 2
                    )
                )
            ),
            "coverage": "two predeclared representative stretch directions, not a full alternate-protocol Hessian",
        }
    force0 = np.asarray(center["forces"]).ravel()
    directional_gradient = float(-q[:, 0] @ force0)
    checks = []
    for h in (0.001, 0.0005):
        secant = (
            samples[f"finite-difference/h{h:g}-s+1"]["energy"]
            - samples[f"finite-difference/h{h:g}-s-1"]["energy"]
        ) / (2 * h)
        checks.append(
            {
                "step_angstrom": h,
                "energy_secant_kcal_mol_A": secant,
                "analytic_directional_gradient": directional_gradient,
                "absolute_error": abs(secant - directional_gradient),
            }
        )
    basis_checks = []
    for i in range(4):
        a, b = samples[f"fresh/validation-{i}"], samples[f"fresh/basis-validation-{i}"]
        de = (b["energy"] - samples["sensitivity/basis-center"]["energy"]) - (
            a["energy"] - center["energy"]
        )
        df = float(np.sqrt(np.mean((np.asarray(a["forces"]) - b["forces"]) ** 2)))
        basis_checks.append(
            {
                "key": a["key"],
                "relative_energy_difference_kcal_mol": de,
                "full_gradient_RMS_difference_kcal_mol_A": df,
            }
        )
    criteria = plan["projection"]["noise_checks"]
    gates = {
        "two_step_curvature": relative_steps
        <= criteria["steps_relative_Frobenius_max"],
        "curvature_symmetry": max(
            r["antisymmetric_relative_Frobenius"] for r in rows.values()
        )
        <= criteria["symmetry_relative_Frobenius_max"],
        "energy_curvature_diagonal": all(
            abs(v)
            <= verification["energy_curvature_diagonal_absolute_allowance_kcal_mol_A2"]
            + verification["energy_curvature_diagonal_relative_max"] * abs(h)
            for r in rows.values()
            for v, h in zip(
                r["energy_gradient_curvature_diagonal_errors"],
                np.diag(r["projected_symmetric_kcal_mol_A2"]),
                strict=True,
            )
        ),
        "energy_gradient_sign_and_units": max(c["absolute_error"] for c in checks)
        <= verification["central_energy_gradient_max_error_kcal_mol_A"],
        **{
            name + "_representative_curvature": max(
                s["Cartesian_relative_Frobenius"], s["projected_relative_Frobenius"]
            )
            <= criteria[
                ("SCF" if name == "scf" else name) + "_representative_relative_max"
            ]
            for name, s in sensitivities.items()
        },
        "fresh_validation_basis_energy": max(
            abs(c["relative_energy_difference_kcal_mol"]) for c in basis_checks
        )
        <= verification["fresh_validation_basis_relative_energy_kcal_mol"],
        "fresh_validation_basis_gradient": max(
            c["full_gradient_RMS_difference_kcal_mol_A"] for c in basis_checks
        )
        <= verification["fresh_validation_basis_full_gradient_RMS_kcal_mol_A"],
    }
    result = signed(
        {
            "schema": "island.j61-projected-curvature.v1",
            "center_identity": digest(center["coordinates"]),
            "subspace_sha256": pin(evidence / "stretching-subspace.json"),
            "steps": rows,
            "step_relative_Frobenius": relative_steps,
            "Richardson_estimated_fine_truncation_relative": relative_steps / 3,
            "Richardson_limitation": "Two-step O(h^2) estimate; SCF/grid noise need not scale as h^2.",
            "sensitivity": sensitivities,
            "energy_gradient_checks": checks,
            "fresh_validation_basis_accuracy": basis_checks,
            "gates": gates,
            "reference_curvature_gate_passed": all(gates.values()),
            "positive_stretch_curvature": bool(np.min(np.linalg.eigvalsh(fine)) > 0),
            "full_molecular_minimum_verified": False,
            "coverage": "Eight-dimensional unweighted stretch subspace only. Bending/torsional stability and transferable force-field adequacy are not established.",
        }
    )
    write(
        evidence / "independent-native-checks.json",
        signed(
            {
                "production_importer_used_for_reconstruction": False,
                "checks": native_checks,
                "passed": True,
            }
        ),
    )
    write(evidence / "projected-Hessian.json", result)
    return result, samples


def fit_training(repo, evidence, reference, samples):
    repo, evidence = Path(repo), Path(evidence)
    require(
        reference["reference_curvature_gate_passed"], "reference_curvature_gate_failed"
    )
    old = repo / "docs/evidence/phase_4j60"
    plan = read(evidence / "plan/scientific-plan.json")
    model = load_model(repo)
    design = read(old / "joint-target-design.json")
    data = read(old / "stages/tight-d-joint-fit/normalized-reference.json")
    training = [
        {**data["samples"]["cc5/" + s["id"]], "key": "J60/" + s["id"]}
        for s in design["selected_training"]
    ]
    training.extend(samples[f"fresh/training-{i}"] for i in range(4))
    validation = [samples[f"fresh/validation-{i}"] for i in range(4)]
    historical = [
        {**data["samples"]["cc5/" + s["id"]], "key": "J60/" + s["id"]}
        for s in design["held_out"]
    ]
    separation = leakage(training, validation, historical)
    write(evidence / "fresh-validation-separation.json", signed(separation))
    center = read(evidence / "stationary-center.json")["center"]
    q = np.asarray(read(evidence / "stretching-subspace.json")["Q"])
    target = reference["steps"]["0.002"]["projected_symmetric_kcal_mol_A2"]
    problem = Problem(model, training, center, q, target)
    old_fit = read(old / "experimental-joint-fit.json")
    initial = np.r_[
        source_baseline(model), [old_fit["coefficients"][k] for k in problem.names[4:]]
    ]
    criteria = plan["fit"]["acceptance"]
    statistical = read(evidence / "plan/fit-statistical-criteria.json")
    _, initial_jacobian = problem.system(initial)
    norms = np.linalg.norm(initial_jacobian, axis=0)
    require(np.all(norms > 1e-12), "pre_fit_zero_parameter_column")
    singular = np.linalg.svd(initial_jacobian / norms, compute_uv=False)
    initial_rank = int(sum(singular > singular[0] * 1e-10))
    initial_condition = float(singular[0] / singular[-1])
    write(
        evidence / "pre-fit-identifiability.json",
        signed(
            {
                "rank": initial_rank,
                "parameters": 24,
                "scaled_condition": initial_condition,
                "scaled_singular_values": singular.tolist(),
                "source_baseline_and_J60_start_only": True,
                "fresh_validation_used": False,
            }
        ),
    )
    require(
        initial_rank == 24 and initial_condition < 1000, "pre_fit_joint_unidentifiable"
    )
    # No validation data enter either objective, initialization, bounds or weights.
    authority = signed(
        {
            "plan_sha256": pin(evidence / "plan/scientific-plan.json"),
            "curvature_identity": reference["identity"],
            "training_geometry_identities": [
                digest(s["coordinates"]) for s in training
            ],
            "validation_geometry_identities": [
                digest(s["coordinates"]) for s in validation
            ],
            "validation_values_used_for_fit_selection": False,
            "source_and_J60_records_modified": False,
        }
    )
    write(evidence / "plan/fit-authority.json", authority)
    coupling_only = problem.fit(initial, fixed_baseline=True)
    write(evidence / "trials/coupling-only.json", signed(coupling_only))
    joint = problem.fit(initial)
    write(evidence / "trials/baseline-joint.json", signed(joint))
    # A source-centered algorithmic start is not an assigned set of zero rows.
    alternative_initial = np.r_[source_baseline(model), np.full(20, 0.1)]
    restarted = problem.fit(alternative_initial)
    write(evidence / "trials/source-centered-restart.json", signed(restarted))
    parameter_vector = np.array([joint["parameters"][k] for k in problem.names])
    write(
        evidence / "remaining-training-discrepancy.json",
        signed(remaining_discrepancy(problem, parameter_vector)),
    )
    restart_vector = np.array([restarted["parameters"][k] for k in problem.names])
    coarse_problem = copy.copy(problem)
    coarse_problem.curvature = np.asarray(
        reference["steps"]["0.004"]["projected_symmetric_kcal_mol_A2"]
    )
    coarse_problem._last = None
    coarse = coarse_problem.fit(parameter_vector)
    write(evidence / "trials/coarse-curvature.json", signed(coarse))
    coarse_values = np.array([coarse["parameters"][k] for k in problem.names])
    # Whole-geometry bootstrap, no re-use of displaced QM gradients as independent observations.
    rng = np.random.default_rng(4061)
    bootstrap, failures = [], []
    for iteration in range(100):
        indices = [0, *rng.integers(1, len(training), len(training) - 1).tolist()]
        resampled = copy.copy(problem)
        resampled.training = [problem.training[i] for i in indices]
        resampled.configurations = [problem.configurations[i] for i in indices]
        resampled._last = None
        try:
            fitted = resampled.fit(parameter_vector, maximum_evaluations=100)
            bootstrap.append([fitted["parameters"][k] for k in problem.names])
        except CompletionError as exc:
            failures.append({"replicate": iteration, "reason": str(exc)})
        if (iteration + 1) % 10 == 0:
            write(
                evidence / "trials/bootstrap-progress.json",
                signed(
                    {
                        "completed_replicates": iteration + 1,
                        "coefficient_samples": bootstrap,
                        "failed_replicates": failures,
                    }
                ),
            )
            print(
                json.dumps(
                    {"bootstrap_completed": iteration + 1, "valid": len(bootstrap)}
                ),
                flush=True,
            )
    values = np.asarray(bootstrap)
    require(len(values) >= 2, "insufficient_bootstrap_replicates")
    sigma = np.std(values, axis=0, ddof=1)
    covariance = np.cov(values, rowvar=False)
    d = np.sqrt(np.maximum(np.diag(covariance), 0))
    correlation = np.divide(
        covariance,
        np.outer(d, d),
        out=np.zeros_like(covariance),
        where=np.outer(d, d) > 0,
    )
    weak = [
        k
        for k, s, v in zip(problem.names, sigma, parameter_vector, strict=True)
        if s > 0.5 * abs(v)
    ]
    write(
        evidence / "training-bootstrap.json",
        signed(
            {
                "names": problem.names,
                "replicates": 100,
                "valid_replicates": len(values),
                "failed_replicates": failures,
                "seed": 4061,
                "coefficient_samples": values.tolist(),
                "sigma": sigma.tolist(),
                "correlation": correlation.tolist(),
                "weak_relative_sigma_parameters": weak,
                "validation_used": False,
                "scope": "Whole-training-configuration resampling with stationary/curvature information fixed as a correlated external block. No independent finite-difference measurement, source/basis/method systematic or transfer uncertainty claim.",
            }
        ),
    )
    joint["fresh_validation"] = metrics(
        model, validation, parameter_vector, joint["training_reference_offset"]
    )
    joint["J60_training_only"] = metrics(
        model, training[:18], parameter_vector, joint["training_reference_offset"]
    )
    joint["historical_inspected_J60_validation"] = metrics(
        model, historical, parameter_vector, joint["training_reference_offset"]
    )
    p_old = initial
    old_offset = old_fit["energy_offset"]
    original = {
        "training": metrics(model, training, p_old, old_offset),
        "fresh_validation": metrics(model, validation, p_old, old_offset),
    }
    p_coupling = np.array([coupling_only["parameters"][k] for k in problem.names])
    coupling_only["fresh_validation"] = metrics(
        model, validation, p_coupling, coupling_only["training_reference_offset"]
    )
    reduction = (
        1 - joint["weighted_residual_SS"] / coupling_only["weighted_residual_SS"]
    )
    natural_scales = np.r_[[1, 1, 100, 100], np.ones(20)]
    restart_error = float(
        np.max(abs(restart_vector - parameter_vector) / natural_scales)
    )
    _, full_jacobian = problem.system(parameter_vector)
    protocol_responses = {}
    for name, response in reference["sensitivity"].items():
        directions = response["directions"]
        reduced_difference = q.T @ np.asarray(response["difference"])
        delta_h = np.zeros((8, 8))
        mask = np.zeros((8, 8), bool)
        for i in range(8):
            for j in range(8):
                estimates = []
                if j in directions:
                    estimates.append(reduced_difference[i, directions.index(j)])
                if i in directions:
                    estimates.append(reduced_difference[j, directions.index(i)])
                if estimates:
                    mask[i, j] = True
                    delta_h[i, j] = np.mean(estimates)
        target_shift = np.zeros(full_jacobian.shape[0])
        center_change = np.asarray(
            samples[f"sensitivity/{name}-center"]["forces"]
        ) - np.asarray(center["forces"])
        target_shift[-63:-36] = -center_change.ravel() / 0.05 / np.sqrt(27)
        target_shift[-36:] = -delta_h[problem.upper] / 30
        dp = np.linalg.lstsq(full_jacobian, -target_shift, rcond=1e-10)[0]
        protocol_responses[name] = {
            "measured_curvature_mask": mask.tolist(),
            "first_order_parameter_changes": dict(
                zip(problem.names, dp.tolist(), strict=True)
            ),
            "interpretation": "Local response to only measured representative curvature entries and center forces; unmeasured protocol changes are not known. Not a full alternate-protocol fit, a bound, or a transferable uncertainty estimate.",
        }
    write(
        evidence / "joint-parameter-jacobian.json",
        signed(
            {
                "names": problem.names,
                "matrix": full_jacobian.tolist(),
                "shape": list(full_jacobian.shape),
                "row_order": "For each of22 training configurations: relative-energy row,27 Cartesian-force rows; then27 stationary-force rows and36 upper-triangle curvature rows, all predeclared weighted.",
                "parameters_checksum": digest(joint["parameters"]),
            }
        ),
    )
    gates = {
        "reference_curvature": reference["reference_curvature_gate_passed"],
        "identifiable_full24": joint["rank"] == 24 and joint["scaled_condition"] < 1000,
        "no_active_bounds": not joint["active_bounds"],
        "fresh_energy": joint["fresh_validation"]["energy_RMSE"]
        <= criteria["validation_energy_RMSE_kcal_mol"],
        "fresh_forces": joint["fresh_validation"]["force_RMS"]
        <= criteria["validation_force_RMS_kcal_mol_A"],
        "stationarity": joint["stationary_force_RMS"]
        <= criteria["stationary_force_RMS_kcal_mol_A"],
        "curvature": joint["projected_curvature_relative_Frobenius"]
        <= criteria["projected_Hessian_relative_Frobenius"],
        "material_training_reduction": reduction
        >= criteria["material_training_weighted_SS_reduction_min_fraction"],
        "bootstrap": len(values) >= statistical["bootstrap_minimum_valid"] and not weak,
        "source_centered_restart": restart_error
        <= statistical["restart_agreement_max_relative_to_natural_scale"],
    }
    result = signed(
        {
            "schema": "island.j61-baseline-joint-refit.v1",
            "fit_authority_identity": authority["identity"],
            "original_J60": original,
            "coupling_only_refit_same_constraints": coupling_only,
            "baseline_adjustable_joint": joint,
            "training_weighted_SS_reduction_fraction": reduction,
            "source_centered_restart_parameter_differences": dict(
                zip(
                    problem.names,
                    (restart_vector - parameter_vector).tolist(),
                    strict=True,
                )
            ),
            "coarse_step_parameter_differences": dict(
                zip(
                    problem.names,
                    (coarse_values - parameter_vector).tolist(),
                    strict=True,
                )
            ),
            "source_centered_restart_scaled_max_difference": restart_error,
            "partial_protocol_parameter_sensitivities": protocol_responses,
            "gates": gates,
            "accepted_coefficients": joint["parameters"] if all(gates.values()) else {},
            "scope": "Exact nine-atom Si-H disiloxane, PBE0/cc-pV(5+d)Z bounded research extension only; no Si-C or methyl-PDMS transfer",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(evidence / "baseline-versus-joint-fit.json", result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=("curvature", "fit"))
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j61"))
    a = p.parse_args()
    reference, samples = curvature_report(a.work, a.evidence)
    print(json.dumps({"reference_gates": reference["gates"]}), flush=True)
    if a.operation == "fit":
        result = fit_training(a.repo, a.evidence, reference, samples)
        print(
            json.dumps(
                {
                    "fit_gates": result["gates"],
                    "accepted_coefficients": len(result["accepted_coefficients"]),
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
