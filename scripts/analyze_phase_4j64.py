"""Replay objective-specific HMDS precision and conditional curvature calibration."""

import argparse
import json
import tempfile
from itertools import combinations
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_baseline import (
    IDENTITY,
    NAMES,
    calibrate,
    corrections,
    extended_parts,
    objective_status,
)
from island.pcff_qm_pilot.hmds_design import metrics
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_paths import leakage, path
from island.pcff_qm_pilot.hmds_precision import (
    decimal_energy,
    path_response,
    protocol_contrast,
    read,
)
from scripts.analyze_phase_4j63 import gate, sample
from scripts.curvature_phase_4j63 import unknown_null_check


def force_difference(left, right):
    delta = np.array(left["forces"]) - np.array(right["forces"])
    local = np.linalg.norm(delta, axis=1)
    return {
        "Cartesian_force_difference": delta.tolist(),
        "RMS_kcal_mol_A": float(np.sqrt(np.mean(delta**2))),
        "maximum_component": float(abs(delta).max()),
        "maximum_local_vector": float(local.max()),
        "worst_stable_atom_id": int(local.argmax()) + 1,
        "per_atom_vector_norm": local.tolist(),
    }


def triple(work, keys, model, center, spec, h):
    samples = [sample(work, key) for key in keys]
    geometries = [path(model, center, spec, v) for v in (0, h, -h)]
    energies = [decimal_energy(work / "outputs" / key) for key in keys]
    return samples, geometries, path_response(samples, geometries, energies, h)


def compare(
    left_work, left_keys, right_work, right_keys, model, x, spec, h, allowed, criteria
):
    left, geo, lr = triple(left_work, left_keys, model, x, spec, h)
    right, _, rr = triple(right_work, right_keys, model, x, spec, h)
    contrasts = [
        protocol_contrast(
            read(left_work / "requests" / (lk + ".json")),
            read(right_work / "requests" / (rk + ".json")),
            allowed,
        )
        for lk, rk in zip(left_keys, right_keys, strict=True)
    ]
    points = []
    for i in (1, 2):
        l, r = left[i], right[i]
        f = force_difference(l, r)
        local = np.array(f["per_atom_vector_norm"])
        allowance = criteria["Cartesian_force_floor_kcal_mol_A"] + criteria[
            "relative"
        ] * np.linalg.norm(r["forces"], axis=1)
        le = (
            decimal_energy(left_work / "outputs" / left_keys[i])[0]
            - decimal_energy(left_work / "outputs" / left_keys[0])[0]
        )
        re = (
            decimal_energy(right_work / "outputs" / right_keys[i])[0]
            - decimal_energy(right_work / "outputs" / right_keys[0])[0]
        )
        ef = 627.5094740630558
        points.append(
            {
                "force": f,
                "maximum_local_vector_gate": bool(np.all(local <= allowance)),
                "Cartesian_components": gate(
                    l["forces"],
                    r["forces"],
                    criteria["Cartesian_force_floor_kcal_mol_A"],
                    criteria["relative"],
                ),
                "relative_deformation_energy": gate(
                    float(le) * ef,
                    float(re) * ef,
                    criteria["energy_floor_kcal_mol"],
                    criteria["relative"],
                ),
                "path_force_or_torque": gate(
                    -np.sum(np.asarray(l["forces"]) * geo[i]["tangent"]),
                    -np.sum(np.asarray(r["forces"]) * geo[i]["tangent"]),
                    criteria["path_force_floor"],
                    criteria["relative"],
                ),
            }
        )
    curvature = gate(
        lr["curvature_from_path_gradient"],
        rr["curvature_from_path_gradient"],
        criteria[spec["kind"] + "_curvature_floor"],
        criteria["relative"],
    )
    return {
        "protocol_contrasts": contrasts,
        "left_keys": left_keys,
        "right_keys": right_keys,
        "left_response": lr,
        "right_response": rr,
        "points": points,
        "curvature": curvature,
        "full_gate_passed": curvature["passed"]
        and all(
            p["maximum_local_vector_gate"]
            and all(
                p[k]["passed"]
                for k in (
                    "Cartesian_components",
                    "relative_deformation_energy",
                    "path_force_or_torque",
                )
            )
            for p in points
        ),
    }


def extended_matrix(model, geometries, anchor):
    a = extended_parts(model, anchor)
    rows = []
    for x in geometries:
        p = extended_parts(model, x)
        rows.extend(
            [
                (p["energy_jacobian"] - a["energy_jacobian"]) / 0.05,
                *p["force_jacobian"].reshape(81, -1) / 0.5,
            ]
        )
    return np.array(rows), a["columns"]


def design_update(model, x, old):
    reference = read(old / "normalized-reference.json")["samples"]
    current = [
        x,
        *[
            s["coordinates"]
            for k, s in sorted(reference.items())
            if k.startswith("paths/")
        ],
    ]
    old_design = read(old / "joint-design.json")
    selected = old_design["selected_additional_configurations"]
    a, names = extended_matrix(model, current, x)
    blocks = {
        r["id"]: extended_matrix(model, [r["coordinates"]], x)[0]
        for r in old_design["candidate_pool"]
    }
    chosen = [r["id"] for r in selected]
    full = metrics(np.vstack([a, *[blocks[n] for n in chosen]]))
    attempts = []
    if not full["identifiable_for_screen"]:
        for _ in range(3):
            options = [
                (metrics(np.vstack([a, *[blocks[n] for n in chosen], b])), name)
                for name, b in blocks.items()
                if name not in chosen
            ]
            m, name = min(
                options, key=lambda v: (-v[0]["rank"], v[0]["scaled_condition"])
            )
            chosen.append(name)
            attempts.append({"added": name, "metrics": m})
            if m["identifiable_for_screen"]:
                break
    smallest, final = chosen, metrics(np.vstack([a, *[blocks[n] for n in chosen]]))
    for n in range(len(chosen) + 1):
        options = []
        for subset in combinations(chosen, n):
            m = metrics(np.vstack([a, *[blocks[key] for key in subset]]))
            if m["identifiable_for_screen"]:
                options.append((m["scaled_condition"], subset, m))
        if options:
            _, smallest, final = min(options, key=lambda row: row[0])
            break
    pool = {r["id"]: r for r in old_design["candidate_pool"]}
    reserved = read(old / "plan/reserved-validation.json")["coordinates"]
    return signed(
        {
            "Hamiltonian": IDENTITY,
            "names": names,
            "current": metrics(a),
            "J63_five": {"ids": [r["id"] for r in selected], "metrics": full},
            "bounded_additional_screen": attempts,
            "selected_geometry_only": [pool[n] for n in smallest],
            "after_selected": final,
            "minimality_scope": "Minimum subset within retained five plus at most3 greedy additions from original48-point pool; not global minimality.",
            "diagnostic_scales": {
                "relative_energy_kcal_mol": 0.05,
                "force_kcal_mol_A": 0.5,
            },
            "QM_executed_on_joint_or_validation": False,
            "fitting_weights_selected": False,
            "validation": leakage(
                model, [*current, *[pool[n]["coordinates"] for n in smallest]], reserved
            ),
            "validation_data_supplied": False,
        }
    )


def analyze(repo, root, work):
    old = repo / "docs/evidence/phase_4j63"
    plan = read(root / "plan/scientific-plan.json")
    model = load_model(repo)
    x = np.array(plan["center_coordinates"])
    sil = next(s for s in plan["paths"] if s["name"] == "siloxane")
    with tempfile.TemporaryDirectory(prefix="j64-old-native-") as tmp:
        historical = Path(tmp) / "native"
        unpack(old / "native-outputs.tar.xz", historical)
        # The additive correction preserves the successful zero-path roundoff
        # control and uses the exact stored center for strict protocol contrasts.
        dense_center = (
            "exact/center"
            if (root / "plan/exact-center-correction.json").exists()
            else "dense/center"
        )
        dense_keys = [dense_center, "dense/plus", "dense/minus"]
        if not all((work / "outputs" / k / "result.json").exists() for k in dense_keys):
            return {"status": "dense_data_missing"}
        samples, _, r = triple(work, dense_keys, model, x, sil, 0.02)
        precision = plan["precision"]
        chain = precision["energy_gradient_chain_rule"]
        consistency = gate(
            r["curvature_from_energies"],
            r["curvature_from_path_gradient"],
            chain["rotation_curvature_floor"],
            chain["relative"],
        )
        contrast = compare(
            historical,
            [
                "controls/scf-center",
                "controls/scf-siloxane-s+1",
                "controls/scf-siloxane-s-1",
            ],
            work,
            dense_keys,
            model,
            x,
            sil,
            0.02,
            ["grid"],
            precision["grid"],
        )
        dense = signed(
            {
                "response": r,
                "energy_gradient_consistency": consistency,
                "grid_only_contrast": contrast,
                "J63_failed_gates_unchanged": True,
                "rounding_sufficient_to_explain_mismatch": abs(
                    r["curvature_from_energies"] - r["curvature_from_path_gradient"]
                )
                <= r["rounding"]["printed_decimal_rounding_bound"]
                + r["rounding"]["binary_serialization_arithmetic_bound"],
            }
        )
        write(root / "dense-grid-decision.json", dense)
        if (work / "outputs/dense/coarse-plus/result.json").exists() and (
            work / "outputs/dense/coarse-minus/result.json"
        ).exists():
            _, _, coarse = triple(
                work,
                [dense_center, "dense/coarse-plus", "dense/coarse-minus"],
                model,
                x,
                sil,
                0.04,
            )
            write(
                root / "torsion-step-control.json",
                signed(
                    {
                        "coarse": coarse,
                        "fine": r,
                        "two_step_gradient_curvature": gate(
                            coarse["curvature_from_path_gradient"],
                            r["curvature_from_path_gradient"],
                            precision["two_step"]["rotation_curvature_floor"],
                            precision["two_step"]["relative"],
                        ),
                        "coarse_energy_gradient": gate(
                            coarse["curvature_from_energies"],
                            coarse["curvature_from_path_gradient"],
                            chain["rotation_curvature_floor"],
                            chain["relative"],
                        ),
                    }
                ),
            )
        if all(
            (work / "outputs/rigid" / k / "result.json").exists()
            for k in ("center", "plus", "minus")
        ):
            rotation = np.array(plan["proper_rigid_rotation"])
            geos = [path(model, x, sil, s) for s in (0, 0.02, -0.02)]
            rotated_geos = [
                {
                    **g,
                    **{
                        n: g[n] @ rotation.T
                        for n in ("coordinates", "tangent", "acceleration")
                    },
                }
                for g in geos
            ]
            keys = ["rigid/center", "rigid/plus", "rigid/minus"]
            rotated = [sample(work, k) for k in keys]
            rr = path_response(
                rotated,
                rotated_geos,
                [decimal_energy(work / "outputs" / k) for k in keys],
                0.02,
            )
            checks = []
            for s, t, dk, rk in zip(samples, rotated, dense_keys, keys, strict=True):
                energy_difference = (
                    float(
                        decimal_energy(work / "outputs" / rk)[0]
                        - decimal_energy(work / "outputs" / dk)[0]
                    )
                    * 627.5094740630558
                )
                checks.append(
                    {
                        "energy_difference_kcal_mol": energy_difference,
                        "forces": force_difference(
                            {"forces": np.array(s["forces"]) @ rotation.T}, t
                        ),
                    }
                )
            write(
                root / "rigid-rotation-control.json",
                signed(
                    {
                        "proper_rotation": rotation.tolist(),
                        "orthogonality_error": float(
                            abs(rotation.T @ rotation - np.eye(3)).max()
                        ),
                        "determinant": float(np.linalg.det(rotation)),
                        "points": checks,
                        "rotated_response": rr,
                        "original_response": r,
                        "energy_gradient_consistency": gate(
                            rr["curvature_from_energies"],
                            rr["curvature_from_path_gradient"],
                            chain["rotation_curvature_floor"],
                            chain["relative"],
                        ),
                        "interpretation": "Whole-molecule proper rotation preserves the graph and path; finite atom-centered quadrature may show orientation sensitivity. This alone is not an engine-defect diagnosis.",
                    }
                ),
            )
        if not (work / "outputs/qz/center/result.json").exists():
            return {"dense_consistency": consistency["passed"]}
        qz = sample(work, "qz/center")
        center_key = "optimization/step-031"
        basis = signed(
            {
                "comparisons": {
                    name: {
                        "only_orbital_basis": protocol_contrast(
                            read(historical / "requests" / (key + ".json")),
                            read(work / "requests/qz/center.json"),
                            ["basis"],
                        ),
                        "forces": force_difference(sample(historical, key), qz),
                    }
                    for name, key in (
                        ("TZVP_to_QZVP", center_key),
                        ("TZVPP_to_QZVP", "controls/basis-center"),
                    )
                },
                "QZVP": {
                    k: read(work / "outputs/qz/center/result.json")[k]
                    for k in ("basis_functions", "basis_shells", "elapsed_seconds")
                },
                "auxiliary": plan["basis_composition"]["def2-universal-jkfit"],
                "scope": "Common retained center force sensitivity only; neither QZVP stationary point nor all-path or complete-basis convergence.",
            }
        )
        write(root / "QZVP-center-decision.json", basis)
        completed = set(read(work / "ledger.json")["calls"])
        needed = {
            f"{prefix}/{name}-{sign}"
            for prefix in ("tight", "dense", "qz")
            for name in ("SiC", "CH")
            for sign in ("plus", "minus")
        }
        if not needed <= completed:
            return {
                "dense_consistency": consistency["passed"],
                "QZVP_center": "completed",
            }
        rows, conditions = {}, {}
        for spec in plan["paths"]:
            if spec["kind"] != "stretch":
                continue
            name = spec["name"]
            primary_keys = [
                center_key,
                f"paths/{name}-h0.003-s+1",
                f"paths/{name}-h0.003-s-1",
            ]
            tight_keys = ["controls/scf-center"]
            # Old center belongs to old work; new matching endpoints belong to new work.
            # Copy exact verified historical center into the analysis temp only.
            combined = Path(tmp) / ("combined-" + name)
            import shutil

            for base, keys in (
                (historical, tight_keys),
                (work, [f"tight/{name}-plus", f"tight/{name}-minus"]),
            ):
                for key in keys:
                    shutil.copytree(base / "outputs" / key, combined / "outputs" / key)
                    write(
                        combined / "requests" / (key + ".json"),
                        read(base / "requests" / (key + ".json")),
                    )
            tight_keys += [f"tight/{name}-plus", f"tight/{name}-minus"]
            dkeys = [dense_center, f"dense/{name}-plus", f"dense/{name}-minus"]
            qkeys = ["qz/center", f"qz/{name}-plus", f"qz/{name}-minus"]
            controls = {
                "SCF": compare(
                    historical,
                    primary_keys,
                    combined,
                    tight_keys,
                    model,
                    x,
                    spec,
                    0.003,
                    ["thresholds"],
                    precision["scf"],
                ),
                "grid": compare(
                    combined,
                    tight_keys,
                    work,
                    dkeys,
                    model,
                    x,
                    spec,
                    0.003,
                    ["grid"],
                    precision["grid"],
                ),
                "basis_QZVP": compare(
                    historical,
                    primary_keys,
                    work,
                    qkeys,
                    model,
                    x,
                    spec,
                    0.003,
                    ["basis"],
                    precision["basis"],
                ),
            }
            null = unknown_null_check(model, x, spec, 0.003)
            fixed = null["verified_fixed_response"]["curvature_from_path_gradient"]
            pc = corrections(model, x)
            gp, gm = [path(model, x, spec, s * 0.003) for s in (1, -1)]
            cp, cm = [corrections(model, g["coordinates"]) for g in (gp, gm)]
            sensitivity = (
                -np.einsum("iap,ia->p", cp["force_jacobian"], gp["tangent"])
                + np.einsum("iap,ia->p", cm["force_jacobian"], gm["tangent"])
            ) / 0.006
            scalar_sensitivity = (
                cp["energy_jacobian"]
                - 2 * pc["energy_jacobian"]
                + cm["energy_jacobian"]
            ) / 0.003**2
            expected = np.array([2.0, 0.0] if name == "SiC" else [0.0, 2.0])
            require(
                np.max(abs(sensitivity - expected)) < 1e-8
                and np.max(abs(scalar_sensitivity - expected)) < 1e-8,
                "quadratic_curvature_sensitivity_error",
            )
            protocols = {
                "TZVP_standard": controls["SCF"]["left_response"],
                "TZVP_tight": controls["SCF"]["right_response"],
                "TZVP_dense_tight": controls["grid"]["right_response"],
                "QZVP_standard": controls["basis_QZVP"]["right_response"],
                "TZVPP_standard": read(old / "deformation-reference-results.json")[
                    "controls"
                ]["basis"]["paths"][name]["control_response"],
            }
            chain_gates = {
                n: gate(
                    rp["curvature_from_energies"],
                    rp["curvature_from_path_gradient"],
                    chain["stretch_curvature_floor"],
                    chain["relative"],
                )
                for n, rp in protocols.items()
            }
            conditions[name] = all(
                c["curvature"]["passed"] for c in controls.values()
            ) and all(c["passed"] for c in chain_gates.values())
            rows[name] = {
                "controls": controls,
                "energy_gradient_checks": chain_gates,
                "protocol_responses": protocols,
                "unknown_null": null,
                "fixed_curvature": fixed,
                "baseline_sensitivity": sensitivity.tolist(),
                "independent_scalar_sensitivity": scalar_sensitivity.tolist(),
                "curvature_objective_precision_passed": conditions[name],
                "conditional_delta_K2": {
                    n: (rp["curvature_from_path_gradient"] - fixed) / 2
                    for n, rp in protocols.items()
                },
            }
        estimates = {}
        for protocol in rows["SiC"]["protocol_responses"]:
            values = calibrate(
                [
                    rows[n]["protocol_responses"][protocol][
                        "curvature_from_path_gradient"
                    ]
                    for n in ("SiC", "CH")
                ],
                [rows[n]["fixed_curvature"] for n in ("SiC", "CH")],
                [rows[n]["baseline_sensitivity"] for n in ("SiC", "CH")],
            )
            estimates[protocol] = dict(zip(NAMES, values.tolist(), strict=True))
        coarse_sensitivity = {
            n: {
                h: (
                    r["QM_total_path_curvature"]
                    - r["model_total_path_curvature_independent_of56_unknowns"]
                )
                / 2
                for h, r in record["steps"].items()
            }
            for n, record in read(old / "source-baseline-path-curvature.json")[
                "paths"
            ].items()
        }
        calibration = signed(
            {
                "Hamiltonian": IDENTITY,
                "paths": rows,
                "conditional_candidates": estimates,
                "retained_two_step_candidates": coarse_sensitivity,
                "units": "kcal/mol/angstrom^2",
                "prefactor": "delta_E=delta_K2*dr^2, no1/2; delta_curvature=2*delta_K2",
                "equilibrium_and_higher_order_rows_changed": False,
                "all_fixed_contributions_included": True,
                "unknown_coefficients_assigned": False,
                "uncertainty": "Step, SCF, grid and basis sensitivities separately reported; no independent statistical uncertainty or complete-basis bound.",
                "scope": "Conditional total-curvature calibration for exact HMDS path/protocol only; not force-qualified, transferable, native source values or validated complete parameters.",
                "objective_status": objective_status(
                    True,
                    all(conditions.values()),
                    consistency["passed"] and contrast["full_gate_passed"],
                    True,
                ),
                "accepted_for_automatic_assignment": [],
                "full56_coefficient_fit_executed": False,
                "source_baseline_reference": model["identity"],
            }
        )
        write(root / "baseline-calibration.json", calibration)
        if not (root / "updated-joint-design.json").exists():
            write(root / "updated-joint-design.json", design_update(model, x, old))
        return calibration["objective_status"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j64")
    )
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            analyze(args.repo.resolve(), args.evidence.resolve(), args.work.resolve()),
            sort_keys=True,
        )
    )
