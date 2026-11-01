"""Reconstruct and compare graph-controlled molecular HMDS responses; never fit them."""

import argparse
import json
import re
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed, write
from island.pcff_qm.reference import CONVERSION
from island.pcff_qm_pilot.hmds_paths import path, response
from island.pcff_qm_pilot.molecular_native import normalized


def read(p):
    return json.loads(Path(p).read_text())


def sample(work, key):
    return normalized(work / "outputs" / key, read(work / "requests" / (key + ".json")))


def independent_native(work, key, expected):
    """Text/ASCII reconstruction, independent of native(), matrix_ascii(), and npy."""
    directory = work / "outputs" / key
    out = (directory / "psi4.out").read_text()
    energy = float(re.findall(r"@(?:DF-)?RKS Final Energy:\s*([-+\d.]+)", out)[-1])
    lines = (directory / "gradient.dat").read_text().splitlines()
    require(int(lines[2]) == 81, "independent_gradient_dimension")
    g = np.empty((27, 3))
    seen = set()
    for row in lines[3:]:
        i, j, value = row.split()
        i, j = int(i), int(j)
        require(
            0 <= i < 27 and 0 <= j < 3 and (i, j) not in seen,
            "independent_gradient_axes",
        )
        seen.add((i, j))
        g[i, j] = float(value)
    require(len(seen) == 81 and np.isfinite(g).all(), "independent_gradient_coverage")
    ef = 4.3597447222060e-18 * 6.02214076e23 / 4184
    lf = 5.29177210544e-11 * 1e10
    require(
        abs(energy * ef - expected["energy"]) < 2e-7
        and np.allclose(-g * ef / lf, expected["forces"], atol=1e-10, rtol=0),
        "independent_native_conversion_or_sign",
    )
    return {
        "energy_difference_kcal_mol": abs(energy * ef - expected["energy"]),
        "maximum_force_difference_kcal_mol_A": float(
            np.max(abs(-g * ef / lf - expected["forces"]))
        ),
    }


def gate(actual, reference, floor, relative):
    a, r = np.asarray(actual), np.asarray(reference)
    difference = np.abs(a - r)
    allowed = floor + relative * abs(r)
    index = int(np.argmax(difference / allowed))
    return {
        "passed": bool(np.all(difference <= allowed)),
        "maximum_absolute_error": float(difference.max()),
        "maximum_error_to_allowance": float((difference / allowed).max()),
        "worst_flat_index": index,
        "absolute_floor": floor,
        "relative": relative,
    }


def analyze(model, plan, center, samples):
    x = np.array(center["coordinates"])
    primary, precision = {}, {}
    for spec in plan["paths"]:
        name = spec["name"]
        p0 = path(model, x, spec, 0)
        steps = {}
        for h in spec["steps"]:
            pp, pm = [path(model, x, spec, s * h) for s in (1, -1)]
            plus, minus = [samples[f"paths/{name}-h{h:g}-s{s:+d}"] for s in (1, -1)]
            steps[str(h)] = response(center, plus, minus, p0, pp, pm, h)
        coarse, fine = [steps[str(h)] for h in spec["steps"]]
        criterion = plan["precision"]["two_step"]
        floor = criterion[spec["kind"] + "_curvature_floor"]
        two = gate(
            coarse["curvature_from_path_gradient"],
            fine["curvature_from_path_gradient"],
            floor,
            criterion["relative"],
        )
        chain = plan["precision"]["energy_gradient_chain_rule"]
        chain_checks = {}
        for h, row in steps.items():
            chain_checks[h] = {
                "energy_gradient": gate(
                    row["curvature_from_energies"],
                    row["curvature_from_path_gradient"],
                    chain[spec["kind"] + "_curvature_floor"],
                    chain["relative"],
                ),
                "curved_path_chain_rule": gate(
                    row["chain_rule_sum"],
                    row["curvature_from_path_gradient"],
                    chain[spec["kind"] + "_curvature_floor"],
                    chain["relative"],
                ),
            }
        primary[name] = {
            "kind": spec["kind"],
            "parameter_units": spec["units"],
            "steps": steps,
            "two_step_curvature": two,
            "chain_checks": chain_checks,
            "two_step_full_Cartesian_response_RMS_difference": float(
                np.sqrt(
                    np.mean(
                        (
                            np.array(coarse["Cartesian_gradient_response"])
                            - fine["Cartesian_gradient_response"]
                        )
                        ** 2
                    )
                )
            ),
            "precision_passed": two["passed"]
            and all(
                g["passed"] for checks in chain_checks.values() for g in checks.values()
            ),
        }
    for name, control in plan["controls"].items():
        criteria = plan["precision"][name]
        other_center = samples[f"controls/{name}-center"]
        paths = {}
        for target in control["paths"]:
            spec = next(s for s in plan["paths"] if s["name"] == target)
            h = min(spec["steps"])
            points = {}
            for sign in (1, -1):
                baseline = samples[f"paths/{target}-h{h:g}-s{sign:+d}"]
                other = samples[f"controls/{name}-{target}-s{sign:+d}"]
                p = path(model, x, spec, sign * h)
                bforce, oforce = np.array(baseline["forces"]), np.array(other["forces"])
                local = np.linalg.norm(bforce - oforce, axis=1)
                points[str(sign)] = {
                    "relative_deformation_energy": gate(
                        baseline["energy"] - center["energy"],
                        other["energy"] - other_center["energy"],
                        criteria["energy_floor_kcal_mol"],
                        criteria["relative"],
                    ),
                    "Cartesian_force": gate(
                        bforce,
                        oforce,
                        criteria["Cartesian_force_floor_kcal_mol_A"],
                        criteria["relative"],
                    ),
                    "local_atom_force_norm": gate(
                        np.linalg.norm(bforce, axis=1),
                        np.linalg.norm(oforce, axis=1),
                        criteria["Cartesian_force_floor_kcal_mol_A"],
                        criteria["relative"],
                    ),
                    "path_force_or_torque": gate(
                        float(np.sum(bforce * p["tangent"])),
                        float(np.sum(oforce * p["tangent"])),
                        criteria["path_force_floor"],
                        criteria["relative"],
                    ),
                    "Cartesian_force_RMS_error": float(
                        np.sqrt(np.mean((bforce - oforce) ** 2))
                    ),
                    "maximum_local_vector_force_error": float(local.max()),
                    "maximum_local_vector_force_site": int(local.argmax()) + 1,
                }
            pp, pm = [path(model, x, spec, s * h) for s in (1, -1)]
            rc = response(
                other_center,
                samples[f"controls/{name}-{target}-s+1"],
                samples[f"controls/{name}-{target}-s-1"],
                path(model, x, spec, 0),
                pp,
                pm,
                h,
            )
            rb = primary[target]["steps"][str(h)]
            curvature = gate(
                rb["curvature_from_path_gradient"],
                rc["curvature_from_path_gradient"],
                criteria[spec["kind"] + "_curvature_floor"],
                criteria["relative"],
            )
            # Vector errors, rather than a stretch-dominated aggregate, enter the gate.
            vector_gate = all(
                np.max(
                    np.linalg.norm(
                        np.array(samples[f"paths/{target}-h{h:g}-s{s:+d}"]["forces"])
                        - samples[f"controls/{name}-{target}-s{s:+d}"]["forces"],
                        axis=1,
                    )
                    / (
                        criteria["Cartesian_force_floor_kcal_mol_A"]
                        + criteria["relative"]
                        * np.linalg.norm(
                            samples[f"controls/{name}-{target}-s{s:+d}"]["forces"],
                            axis=1,
                        )
                    )
                )
                <= 1
                for s in (1, -1)
            )
            paths[target] = {
                "points": points,
                "control_response": rc,
                "curvature": curvature,
                "full_Cartesian_response_RMS_difference": float(
                    np.sqrt(
                        np.mean(
                            (
                                np.array(rb["Cartesian_gradient_response"])
                                - rc["Cartesian_gradient_response"]
                            )
                            ** 2
                        )
                    )
                ),
                "maximum_local_vector_force_gate": bool(vector_gate),
                "passed": curvature["passed"]
                and vector_gate
                and all(
                    v["passed"]
                    for p in points.values()
                    for v in p.values()
                    if isinstance(v, dict) and "passed" in v
                ),
            }
        precision[name] = {
            "paths": paths,
            "passed": all(r["passed"] for r in paths.values()),
            "protocol": other_center["protocol"],
            "same_geometries": True,
        }
    return signed(
        {
            "schema": "island.j63-HMDS-path-precision.v1",
            "primary": primary,
            "controls": precision,
            "bounded_precision_passed": all(
                v["precision_passed"] for v in primary.values()
            )
            and all(v["passed"] for v in precision.values()),
            "basis_scope": "TZVP/TZVPP additional polarization; not a higher-cardinal or complete-basis test",
            "method_accuracy_established": False,
            "minimum_verified": False,
            "curvature_units": "kcal/mol/angstrom^2 (stretch), kcal/mol/radian^2 (rotation)",
            "response_units": "kcal/mol/angstrom per path-unit",
            "conversion": CONVERSION,
            "interpretation": "Total molecular path curvature includes gradient dot path acceleration; no isolated physical coefficient inferred.",
        }
    )


if __name__ == "__main__":
    from island.pcff_qm_pilot.hmds_hamiltonian import load_model

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    p.add_argument("--work", type=Path, required=True)
    a = p.parse_args()
    ledger = read(a.work / "ledger.json")
    samples = {
        k: sample(a.work, k)
        for k, v in ledger["calls"].items()
        if v["status"] == "completed"
    }
    center = read(a.evidence / "stationary-center.json")["center"]
    result = analyze(
        load_model(a.repo),
        read(a.evidence / "plan/scientific-plan.json"),
        center,
        samples,
    )
    write(a.evidence / "deformation-reference-results.json", result)
    print(json.dumps({"bounded_precision_passed": result["bounded_precision_passed"]}))
