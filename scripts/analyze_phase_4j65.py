"""Native controls, exact-protocol intake and the preregistered coupled trial."""

import argparse
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm.reference import CONVERSION
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_fit import (
    accepts,
    features,
    infer,
    performance,
    prediction_uncertainty,
    uncertainty,
    weighted_problem,
)
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_paths import leakage, path
from island.pcff_qm_pilot.hmds_precision import (
    decimal_energy,
    path_response,
    protocol_contrast,
    read,
)
from scripts.analyze_phase_4j63 import gate, independent_native, sample
from scripts.analyze_phase_4j64 import force_difference
from scripts.independent_phase_4j64 import (
    ascii_gradient,
    independent_path,
    native_decimal_total,
)


def independent_value(source, key, value):
    directory = source / "outputs" / key
    energy = native_decimal_total(directory)
    factor = 4.3597447222060e-18 * 6.02214076e23 / 4184
    bohr = 5.29177210544e-11 * 1e10
    force = -ascii_gradient(directory) * factor / bohr
    require(
        np.max(abs(force - np.asarray(value["forces"]))) < 1e-9,
        "independent_native_force_conversion_or_sign",
    )
    require(
        abs(float(energy) * factor - value["energy"]) < 1e-8,
        "independent_native_energy_conversion",
    )
    return {
        "native_decimal_total_hartree": str(energy),
        "Hartree_to_kcal_mol": factor,
        "bohr_to_A": bohr,
        "maximum_force_reconstruction_error": float(
            np.max(abs(force - value["forces"]))
        ),
    }


def parent(repo, work):
    destination = work / "parents/J64"
    if not destination.exists():
        unpack(repo / "docs/evidence/phase_4j64/native-outputs.tar.xz", destination)
    return destination


def relative(energy, center):
    with localcontext() as ctx:
        ctx.prec = 50
        return float(
            (Decimal(energy) - Decimal(center))
            * Decimal(str(CONVERSION["energy_factor"]))
        )


def observation(repo, work, plan, row):
    wanted = request(plan, row["coordinates"])
    if row["id"] in plan["reuse"]:
        source = parent(repo, work)
        binding = plan["reuse"][row["id"]]
        key = binding["key"]
        producer = {
            **binding,
            "kind": "reused_original_execution",
            "new_native_attempt": False,
        }
        require(
            pin(source / "outputs" / key / "result.json")
            == binding["producer_result_sha256"],
            "reused_producer_mutation",
        )
    else:
        source, key = work, "training/" + row["id"]
        if not (source / "outputs" / key / "result.json").exists():
            return None
        producer = {
            "kind": "new_J65_execution",
            "key": key,
            "producer_result_sha256": pin(source / "outputs" / key / "result.json"),
        }
    original = read(source / "requests" / (key + ".json"))
    require(
        original == wanted, "exact_geometry_and_protocol_reuse_required", id=row["id"]
    )
    value = sample(source, key)
    check = {
        **independent_native(source, key, value),
        "independent_decimal_ASCII": independent_value(source, key, value),
    }
    energy = str(decimal_energy(source / "outputs" / key)[0])
    center = str(decimal_energy(parent(repo, work) / "outputs/exact/center")[0])
    return {
        "id": row["id"],
        "block": row["block"],
        "coordinates": value["coordinates"],
        "forces": value["forces"],
        "relative_energy_kcal_mol": relative(energy, center),
        "native_energy_hartree_decimal": energy,
        "native_center_energy_hartree_decimal": center,
        "original_normalized_sample": value,
        "original_request": original,
        "producer": producer,
        "independent_native_check": check,
    }


def assemble(repo, work, plan, include_joints=True):
    rows = plan["primary"] + (plan["joints"] if include_joints else [])
    observations = [observation(repo, work, plan, row) for row in rows]
    return [r for r in observations if r is not None], [
        row["id"] for row, r in zip(rows, observations, strict=True) if r is None
    ]


def primary_gates(repo, root, work):
    plan = read(root / "plan/scientific-plan.json")
    model = load_model(repo)
    data, missing = assemble(repo, work, plan, False)
    report = {
        "scientific_plan_identity": plan["identity"],
        "observations": len(data),
        "missing": missing,
        "joints_authorized": False,
        "checks": [],
    }
    if not missing:
        indexed = {r["id"]: r for r in data}
        x = np.array(plan["center_coordinates"])
        rules = plan["precision"]["inherited"]["energy_gradient_chain_rule"]
        for spec in plan["paths"]:
            for h in spec["steps"]:
                rows = [indexed["center"]] + [
                    indexed[f"paths/{spec['name']}-h{h:g}-s{s:+d}"] for s in (1, -1)
                ]
                geos = [path(model, x, spec, v) for v in (0, h, -h)]
                native = [
                    (
                        Decimal(r["native_energy_hartree_decimal"]),
                        float(r["native_energy_hartree_decimal"]),
                    )
                    for r in rows
                ]
                response = path_response(
                    [r["original_normalized_sample"] for r in rows], geos, native, h
                )
                criterion = gate(
                    response["curvature_from_energies"],
                    response["curvature_from_path_gradient"],
                    rules[spec["kind"] + "_curvature_floor"],
                    rules["relative"],
                )
                report["checks"].append(
                    {
                        "path": spec["name"],
                        "step": h,
                        "response": response,
                        "energy_gradient_consistency": criterion,
                    }
                )
        report["joints_authorized"] = all(
            r["energy_gradient_consistency"]["passed"] for r in report["checks"]
        )
    result = signed(report)
    write(root / "primary-numerical-gates.json", result)
    return result


def controls(repo, root, work):
    plan = read(root / "plan/scientific-plan.json")
    model = load_model(repo)
    old = parent(repo, work)
    x = np.array(plan["center_coordinates"])
    rows = []
    for prediction in plan["cross_site_predictions"]:
        spec = prediction["spec"]
        name = spec["name"]
        locations = [
            (old, "qz/center"),
            (work, "controls/" + name + "-plus"),
            (work, "controls/" + name + "-minus"),
        ]
        if not all((w / "outputs" / k / "result.json").exists() for w, k in locations):
            return None
        values = [sample(w, k) for w, k in locations]
        native_checks = [
            independent_value(w, k, v)
            for (w, k), v in zip(locations, values, strict=True)
        ]
        for w, k in locations[1:]:
            protocol_contrast(
                read(old / "requests/qz/center.json"),
                read(w / "requests" / (k + ".json")),
                ["coordinates"],
            )
        response = path_response(
            values,
            [path(model, x, spec, v) for v in (0, 0.003, -0.003)],
            [decimal_energy(w / "outputs" / k) for w, k in locations],
            0.003,
        )
        predicted = prediction["predicted_QZVP_total_curvature"]
        actual = response["curvature_from_path_gradient"]
        independent_derivatives = []
        for (w, k), s in zip(locations[1:], (1, -1), strict=True):
            coordinates = read(w / "requests" / (k + ".json"))["coordinates"]
            tangent, _ = independent_path(model, x, coordinates, spec)
            independent_derivatives.append(
                float(
                    np.sum(ascii_gradient(w / "outputs" / k) * tangent)
                    * native_checks[0]["Hartree_to_kcal_mol"]
                    / native_checks[0]["bohr_to_A"]
                )
            )
        reconstructed = (
            independent_derivatives[0] - independent_derivatives[1]
        ) / 0.006
        require(abs(reconstructed - actual) < 1e-9, "independent_cross_site_curvature")
        rows.append(
            {
                "frozen_prediction": prediction,
                "response": response,
                "observed_minus_predicted_curvature": actual - predicted,
                "relative_error": (actual - predicted) / predicted,
                "corrections_refitted": False,
                "independent_native_checks": native_checks,
                "independent_ASCII_path_curvature": reconstructed,
                "scope": "Nearly symmetry-equivalent other sites in same conformer; path-alignment RMS equivalence passes1e-3A and fails1e-6A. Not independent basins or transferable validation.",
            }
        )
    key = "controls/QZVP-dense-center"
    if not (work / "outputs" / key / "result.json").exists():
        return None
    original = read(old / "requests/qz/center.json")
    new = read(work / "requests" / (key + ".json"))
    contrast = protocol_contrast(original, new, ["grid"])
    left, right = sample(old, "qz/center"), sample(work, key)
    answer = signed(
        {
            "scientific_plan_identity": plan["identity"],
            "cross_site": rows,
            "QZVP_grid_only": {
                "contrast": contrast,
                "force_difference": force_difference(left, right),
                "energy_change_kcal_mol": relative(
                    str(decimal_energy(work / "outputs" / key)[0]),
                    str(decimal_energy(old / "outputs/qz/center")[0]),
                ),
                "basis_changed": False,
                "stationary_point_established": False,
            },
            "historical_basis_failures_unchanged": True,
            "J64_complete_force_fit_authorization": False,
            "historical_decisions_changed": False,
            "method_specific_target_is_separate": True,
        }
    )
    write(root / "targeted-controls.json", answer)
    return answer


def numerical(repo, root, work, data):
    plan = read(root / "plan/scientific-plan.json")
    p = plan["precision"]
    expected = [
        k
        for k, j in read(work / "jobs.json")["jobs"].items()
        if j["stage"] == "precision"
    ]
    missing = [
        k for k in expected if not (work / "outputs" / k / "result.json").exists()
    ]
    result = {
        "scientific_plan_identity": plan["identity"],
        "missing": missing,
        "eligible": False,
    }
    if missing:
        return signed(result)
    indexed = {r["id"]: r for r in data}
    center = indexed["center"]
    joint = indexed[p["FD_center_id"]]
    grid = []
    for label, row in (("center", center), ("joint040", joint)):
        key = "precision/grid-" + label
        r = sample(work, key)
        protocol_contrast(
            row["original_request"], read(work / "requests" / (key + ".json")), ["grid"]
        )
        delta = force_difference(row["original_normalized_sample"], r)
        rmsref = float(np.sqrt(np.mean(np.asarray(r["forces"]) ** 2)))
        allowance = p["force_local_floor"] + p["relative"] * np.linalg.norm(
            r["forces"], axis=1
        )
        grid.append(
            {
                "id": row["id"],
                "force_difference": delta,
                "force_RMS_passed": delta["RMS_kcal_mol_A"]
                <= p["force_RMS_floor"] + p["relative"] * rmsref,
                "all_local_forces_passed": bool(
                    np.all(np.array(delta["per_atom_vector_norm"]) <= allowance)
                ),
            }
        )
    native = lambda key: str(decimal_energy(work / "outputs" / key)[0])
    dense_delta = joint["relative_energy_kcal_mol"]
    finer_delta = relative(
        native("precision/grid-joint040"), native("precision/grid-center")
    )
    energy_grid = gate(
        dense_delta, finer_delta, p["grid_energy_floor_kcal_mol"], p["relative"]
    )
    scf = sample(work, "precision/scf-center")
    protocol_contrast(
        center["original_request"],
        read(work / "requests/precision/scf-center.json"),
        ["thresholds"],
    )
    sf = force_difference(center["original_normalized_sample"], scf)
    scf_energy = abs(
        relative(
            native("precision/scf-center"), center["native_energy_hartree_decimal"]
        )
    )
    scf_pass = (
        scf_energy <= p["SCF_energy_floor_kcal_mol"]
        and sf["RMS_kcal_mol_A"]
        <= p["force_RMS_floor"]
        + p["relative"] * float(np.sqrt(np.mean(np.asarray(scf["forces"]) ** 2)))
        and np.all(
            np.array(sf["per_atom_vector_norm"])
            <= p["force_local_floor"]
            + p["relative"] * np.linalg.norm(scf["forces"], axis=1)
        )
    )
    direction = np.array(p["FD_direction"])
    analytic = float(-np.sum(np.array(joint["forces"]) * direction))
    fd = []
    for h in p["FD_steps_A"]:
        derivative = relative(
            native(f"precision/FD-{h}-plus"), native(f"precision/FD-{h}-minus")
        ) / (2 * h)
        fd.append(
            {
                "step_A": h,
                "energy_derived_gradient": derivative,
                "Cartesian_gradient_contraction": analytic,
                "check": gate(
                    derivative,
                    analytic,
                    p["FD_absolute_floor_kcal_mol_A"],
                    p["FD_relative"],
                ),
            }
        )
    result.update(
        grid=grid,
        relative_energy_grid=energy_grid,
        SCF={
            "energy_change": scf_energy,
            "force_difference": sf,
            "passed": bool(scf_pass),
        },
        finite_difference=fd,
        eligible=bool(
            all(r["force_RMS_passed"] and r["all_local_forces_passed"] for r in grid)
            and energy_grid["passed"]
            and scf_pass
            and all(r["check"]["passed"] for r in fd)
        ),
    )
    return signed(result)


def trial(repo, root, work):
    plan = read(root / "plan/scientific-plan.json")
    model = load_model(repo)
    observations, missing = assemble(repo, work, plan)
    gates = numerical(repo, root, work, observations)
    write(root / "fitting-numerical-checks.json", gates)
    dataset = signed(
        {
            "scientific_plan_identity": plan["identity"],
            "protocol": plan["protocol"],
            "observations": observations,
            "missing": missing,
            "mixed_protocols": False,
            "reused": sum(
                r["producer"]["kind"] == "reused_original_execution"
                for r in observations
            ),
        }
    )
    write(root / "fitting-dataset.json", dataset)
    primary = primary_gates(repo, root, work)
    decision = {
        "scientific_plan_identity": plan["identity"],
        "fit_executed": False,
        "validation_authorized": False,
        "missing": missing,
        "numerical_eligible": gates["eligible"],
        "physical_validation": "not_established",
        "basis_qualification": "not_established_historical_failures_unchanged",
        "parameter_completeness": "incomplete",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    if missing or not gates["eligible"] or not primary["joints_authorized"]:
        write(root / "trial-fit.json", signed(decision))
        return decision
    data = features(model, observations)
    matrix, residual = weighted_problem(data, plan["fit"])
    np.savez_compressed(
        root / "actual-fitting-system.npz",
        weighted_Jacobian=matrix,
        QM_minus_verified_fixed=residual,
        names=np.array(data["names"] + ["energy_offset"]),
    )
    fit = infer(data, plan["parameter_specifications"], plan["fit"])
    scores = performance(data, fit["values"], [r["id"] for r in observations])
    spread = uncertainty(
        data, observations, fit["values"], plan["parameter_specifications"], plan["fit"]
    )
    write(root / "block-uncertainty.json", signed(spread))
    eligible = (
        fit["solver_success"]
        and fit["identifiable"]
        and not fit["active_bound_parameters"]
        and accepts(scores, plan["fit"])
        and spread["qualified"]
    )
    decision.update(
        fit_executed=True,
        fit=fit,
        training=scores,
        training_performance_passed=accepts(scores, plan["fit"]),
        uncertainty_qualified=spread["qualified"],
        validation_authorized=bool(eligible),
        coefficient_origin="fitted_external method-specific research trial; not native PCFF source rows",
        accepted_candidates=[],
        scope="Exact27atom HMDS, frozen denseTZVP target; scientific/chemical transfer validation separate",
        baseline_diagnostics={
            "selected_protocol_parent": "phase_4j64/TZVP_dense_tight",
            "conditional_curvature_estimates": read(
                repo / "docs/evidence/phase_4j64/baseline-calibration.json"
            )["conditional_candidates"]["TZVP_dense_tight"],
            "use_in_objective": "diagnostic only; no double-counted gradient constraints and no QZVP imposed correction",
            "joint_fitted_corrections": dict(
                zip(data["names"][-2:], fit["values"][-3:-1], strict=True)
            ),
        },
    )
    write(root / "trial-fit.json", signed(decision))
    return decision


def validation_preflight(repo, root, work):
    plan = read(root / "plan/scientific-plan.json")
    decision = read(root / "trial-fit.json")
    require(decision["validation_authorized"], "untouched_validation_not_activated")
    model = load_model(repo)
    seen = []
    sources = []
    for file in (work / "requests").rglob("*.json"):
        key = str(file.relative_to(work / "requests"))[:-5]
        if (work / "outputs" / key / "result.json").exists():
            require(not key.startswith("validation/"), "validation_already_observed")
            sample(
                work, key
            )  # Verify completed native hashes/mapping before screening.
            seen.append(read(file)["coordinates"])
            sources.append({"key": key, "request_sha256": pin(file)})
    seen.extend(row["coordinates"] for row in plan["primary"])
    result = signed(
        {
            "scientific_plan_identity": plan["identity"],
            "trial_fit_sha256": pin(root / "trial-fit.json"),
            "model_weights_bounds_frozen_before_validation": True,
            "observed_sources": sources,
            "leakage": leakage(model, seen, plan["validation_coordinates"]),
        }
    )
    write(root / "validation-activation.json", result)
    return result


def validation(repo, root, work):
    plan = read(root / "plan/scientific-plan.json")
    decision = read(root / "trial-fit.json")
    require(decision["validation_authorized"], "untouched_validation_not_activated")
    model = load_model(repo)
    observed = [
        read(f)["coordinates"]
        for f in (work / "requests").rglob("*.json")
        if not str(f.relative_to(work / "requests")).startswith("validation/")
        and (
            work
            / "outputs"
            / str(f.relative_to(work / "requests"))[:-5]
            / "result.json"
        ).exists()
    ]
    observed += [r["coordinates"] for r in plan["primary"]]
    check = leakage(model, observed, plan["validation_coordinates"])
    rows = []
    center = str(decimal_energy(parent(repo, work) / "outputs/exact/center")[0])
    for i, geo in enumerate(plan["validation_coordinates"]):
        key = f"validation/config-{i}"
        value = sample(work, key)
        require(
            read(work / "requests" / (key + ".json")) == request(plan, geo),
            "validation_protocol_mismatch",
        )
        rows.append(
            {
                "id": key,
                "coordinates": geo,
                "forces": value["forces"],
                "relative_energy_kcal_mol": relative(
                    str(decimal_energy(work / "outputs" / key)[0]), center
                ),
            }
        )
    data = features(model, rows)
    score = performance(data, decision["fit"]["values"], [r["id"] for r in rows])
    bootstrap = read(root / "block-uncertainty.json")
    spread = prediction_uncertainty(
        data,
        [r["values"] for r in bootstrap["attempts"] if r["eligible"]],
        decision["fit"]["values"],
        [r["id"] for r in rows],
    )
    limits = plan["fit"]["uncertainty"]
    prediction_passed = (
        spread["energy_RMS95"] <= limits["prediction_energy_RMS95_max"]
        and spread["force_RMS95"] <= limits["prediction_force_RMS95_max"]
    )
    result = signed(
        {
            "scientific_plan_identity": plan["identity"],
            "untouched_values_inspected_only_after_training_eligibility": True,
            "leakage": check,
            "performance": score,
            "method_specific_target_passed": bool(
                accepts(score, plan["fit"]) and prediction_passed
            ),
            "prediction_uncertainty": spread,
            "physical_validation": "not_established",
            "automatic_assignment": False,
        }
    )
    write(root / "untouched-validation.json", result)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "stage",
        choices=["controls", "primary", "trial", "validation", "validation-preflight"],
    )
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j65"))
    p.add_argument("--work", type=Path, required=True)
    a = p.parse_args()
    function = {
        "controls": controls,
        "primary": primary_gates,
        "trial": trial,
        "validation": validation,
        "validation-preflight": validation_preflight,
    }[a.stage]
    result = function(a.repo.resolve(), a.evidence.resolve(), a.work.resolve())
    print(
        {
            k: v
            for k, v in (result or {}).items()
            if k
            in (
                "fit_executed",
                "validation_authorized",
                "joints_authorized",
                "observations",
                "missing",
                "numerical_eligible",
            )
        }
    )
