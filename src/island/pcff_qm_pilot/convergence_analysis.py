"""Same-geometry relative energies, gradients and mixed scan derivatives."""

import json
import re
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, require, signed
from island.pcff_qm.geometry import measure
from island.pcff_qm.reference import CONVERSION

from .campaign import ENGINE_BOHR_ANGSTROM, verify_completed
from .importer import UNITS, matrix_ascii
from .importer import native as j58_native


def mixed(samples, repo, protocol):
    repo = Path(repo)
    out = {}
    for name, h in (("angle-angle", 0.08), ("bond-bond", 0.04)):
        keys = [f"{protocol}/{name}-{i:04}" for i in (0, 4, 20, 24)]
        if not all(k in samples for k in keys):
            continue
        e = [samples[k]["energy"] for k in keys]
        energy = (e[3] - e[2] - e[1] + e[0]) / (4 * h * h)
        scan = json.loads(
            (repo / f"docs/evidence/phase_4j58/plan/{name}.json").read_text()
        )
        # Independent geometry path tangent from three original neighboring scan
        # points. These derivatives include actual J58 constraint correction.
        grid = np.array([p["coordinates"] for p in scan["points"]]).reshape(5, 5, 9, 3)
        derivatives = []
        for i, j in ((0, 0), (0, 4), (4, 0), (4, 4)):
            tangent = (
                (-3 * grid[0, j] + 4 * grid[1, j] - grid[2, j]) / (h)
                if i == 0
                else (3 * grid[4, j] - 4 * grid[3, j] + grid[2, j]) / (h)
            )
            # The five-point interval is h/2: denominator 2*(h/2)=h.
            forces = np.array(samples[f"{protocol}/{name}-{i * 5 + j:04}"]["forces"])
            derivatives.append(float(np.sum(-forces * tangent)))
        gradient = (
            (derivatives[3] - derivatives[2]) + (derivatives[1] - derivatives[0])
        ) / (4 * h)
        out[name] = {
            "energy_mixed_secant": energy,
            "gradient_mixed_secant": gradient,
            "units": "kcal/mol/radian^2"
            if name == "angle-angle"
            else "kcal/mol/angstrom^2",
            "range": [-h, h],
            "interpretation": "finite-range total-QM path mixed derivative, includes all interactions and geometry-path terms; not an MM cross coefficient",
            "gradient_path_tangent": "three-neighbor one-sided finite geometry path tangent at original grid corners; secant/truncation difference retained",
        }
    return out


def analyze(repo, evidence, normalized, j58work):
    repo, evidence, j58work = Path(repo), Path(evidence), Path(j58work)
    plan = json.loads((evidence / "plan/convergence-plan.json").read_text())
    s = dict(normalized["samples"])
    # Reuse immutable historical results, never rerun or blend their protocols.
    for label in plan["geometries"]:
        if label == "center":
            key = "angle-angle/point-0012"
        else:
            name, number = label.rsplit("-", 1)
            key = name + "/point-" + number
        req = json.loads((j58work / "requests" / (key + ".json")).read_text())
        r, a = j58_native(j58work / "outputs" / key, req)
        s["svp/" + label] = {
            "energy": r["energy_hartree"] * CONVERSION["energy_factor"],
            "forces": (
                -a["gradient"]
                * CONVERSION["energy_factor"]
                / CONVERSION["distance_factor"]
            ).tolist(),
        }
    comparisons = []
    mixed_results = {p: mixed(s, repo, p) for p in ("svp", "tz", "qz", "grid")}
    pairs = [
        ("svp", "tz", "basis_SVP_TZVP_with_DF_algorithm_control"),
        ("tz", "qz", "basis_TZVP_QZVP"),
        ("qz", "grid", "grid_99x590_150x974"),
        ("grid", "scf", "SCF_thresholds"),
        ("tz", "pk-control", "DF_PK_algorithm"),
        ("grid", "method", "method_PBE0_B3LYP"),
    ]
    for a, b, label in pairs:
        center_a = s[a + "/center"]
        center_b = s[b + "/center"]
        records = []
        common = [
            k.split("/", 1)[1]
            for k in s
            if k.startswith(a + "/") and b + "/" + k.split("/", 1)[1] in s
        ]
        for point in sorted(common):
            u, v = s[a + "/" + point], s[b + "/" + point]
            ea = u["energy"] - center_a["energy"]
            eb = v["energy"] - center_b["energy"]
            gf = np.array(v["forces"]) - np.array(u["forces"])
            relative = gf - (
                np.array(center_b["forces"]) - np.array(center_a["forces"])
            )
            records.append(
                {
                    "geometry": point,
                    "relative_energy_a": ea,
                    "relative_energy_b": eb,
                    "deformation_energy_difference": eb - ea,
                    "gradient_RMS_difference": float(np.sqrt(np.mean(gf**2))),
                    "deformation_gradient_RMS_difference": float(
                        np.sqrt(np.mean(relative**2))
                    ),
                    "gradient_max_difference": float(np.max(np.abs(gf))),
                }
            )
        differences = {}
        if a in mixed_results and b in mixed_results:
            for scan in mixed_results[a].keys() & mixed_results[b].keys():
                u = mixed_results[a][scan]["energy_mixed_secant"]
                v = mixed_results[b][scan]["energy_mixed_secant"]
                floor = 1.0 if scan == "angle-angle" else 5.0
                differences[scan] = {
                    "absolute_difference": v - u,
                    "relative_difference_with_predeclared_floor": abs(v - u)
                    / max(abs(v), floor),
                }
        max_e = max(abs(r["deformation_energy_difference"]) for r in records)
        max_g = max(r["gradient_RMS_difference"] for r in records)
        if label == "basis_TZVP_QZVP":
            passed = (
                max_e <= 0.05
                and max_g <= 0.5
                and all(
                    v["relative_difference_with_predeclared_floor"] <= 0.1
                    for v in differences.values()
                )
            )
        elif label == "grid_99x590_150x974":
            passed = (
                max_e <= 0.005
                and max_g <= 0.02
                and all(
                    v["relative_difference_with_predeclared_floor"] <= 0.02
                    for v in differences.values()
                )
            )
        elif label == "SCF_thresholds":
            passed = max_e <= 0.0001 and max_g <= 0.005
        elif label == "DF_PK_algorithm":
            passed = max_e <= 0.005 and max_g <= 0.02
        else:
            passed = None
        comparisons.append(
            {
                "effect": label,
                "protocols": [a, b],
                "geometries": records,
                "max_deformation_energy_difference_kcal_mol": max_e,
                "max_gradient_RMS_difference_kcal_mol_A": max_g,
                "mixed_derivative_differences": differences,
                "predeclared_criteria_passed": passed,
            }
        )
    first = s["grid/angle-angle-0000"]
    repeat = s["thread-control/angle-angle-0000"]
    repeat_de = repeat["energy"] - first["energy"]
    repeat_dg = float(
        np.sqrt(np.mean((np.array(repeat["forces"]) - np.array(first["forces"])) ** 2))
    )
    return signed(
        {
            "schema": "island.j59-QM-convergence.v1",
            "thread_repeat_control": {
                "energy_difference_kcal_mol": repeat_de,
                "gradient_RMS_difference_kcal_mol_A": repeat_dg,
                "within_existing_grid_error_budget": abs(repeat_de) <= 0.005
                and repeat_dg <= 0.02,
                "purpose": "identical geometry/method/basis/grid technical repeat at 2 versus 8 threads; not a QM adequacy criterion",
            },
            "plan_identity": plan["identity"],
            "comparisons": comparisons,
            "total_QM_mixed_derivatives": mixed_results,
            "protocols_separate": True,
            "absolute_electronic_energy_offsets_not_accuracy_criteria": True,
            "basis_convergence_established_for_objectives": next(
                r["predeclared_criteria_passed"]
                for r in comparisons
                if r["effect"] == "basis_TZVP_QZVP"
            ),
            "final_QM_adequacy": "not_established",
            "fit_authorized": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def center_report(work, repo):
    work, repo = Path(work), Path(repo)
    result_path = work / "outputs/optimize/qz-dense/result.json"
    native_result = (
        json.loads(result_path.read_text()) if result_path.exists() else None
    )
    r = native_result or {
        "status": "interrupted_without_native_metadata",
        "operation": "optimize",
    }
    out = (work / "outputs/optimize/qz-dense/psi4.out").read_text()
    request = json.loads((work / "requests/optimize/qz-dense.json").read_text())
    if r["status"] == "completed":
        verify_completed(result_path.parent, digest(request))
    # This installed analytic-gradient route prints actual SCF geometries, not
    # Input QCSchema blocks. Retain its printed precision without claiming a
    # full-precision trajectory that was never exported.
    geometries = re.findall(
        r"Geometry \(in (Angstrom|Bohr)\), charge = (\d+), multiplicity = (\d+):((?:(?!Geometry \(in ).)*?)Running in",
        out,
        re.DOTALL,
    )
    energies = [
        float(v) for v in re.findall(r"@(?:DF-)?RKS Final Energy:\s*([-+\d.]+)", out)
    ]
    bases = re.findall(r"==> Primary Basis <==\s*Basis Set:\s*(\S+)", out)
    functionals = re.findall(r"Composite Functional:\s*(\S+)", out)
    require(
        bases
        and functionals
        and all(v.lower() == request["basis"] for v in bases)
        and all(v.lower() == request["method"] for v in functionals),
        "optimization_native_protocol_mutation",
    )
    trajectory = []
    require(len(geometries) >= len(energies), "optimization_trajectory_coverage")
    for (unit, charge, spin, body), energy in zip(geometries, energies):
        rows = re.findall(
            r"^\s+(SI|O|H)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s+[-+\d.]+\s*$",
            body,
            re.MULTILINE,
        )
        require(
            [row[0] for row in rows] == [s.upper() for s in request["elements"]]
            and int(charge) == request["charge"]
            and int(spin) == request["multiplicity"],
            "optimization_native_mapping_mutation",
        )
        xyz = np.asarray([[float(v) for v in row[1:]] for row in rows])
        if unit == "Bohr":
            xyz *= ENGINE_BOHR_ANGSTROM
        trajectory.append(
            {
                "coordinates": xyz.tolist(),
                "energy_hartree": energy,
                "driver": "gradient",
                "source": "native SCF Geometry table (12 decimal places) and following converged Final Energy; no optimization success inferred",
            }
        )
    old = json.loads(
        (repo / "docs/evidence/phase_4j58/plan/angle-angle.json").read_text()
    )
    result = {
        "optimization_status": r["status"],
        "native_result": native_result,
        "native_completed_SCF_trajectory": trajectory,
        "SCF_convergence_markers": out.count("Energy and wave function converged."),
        "unfinished_geometry_count": len(geometries) - len(energies),
        "parent_J58_center_sha256": old["center_native_result_sha256"],
        "historical_scan_modified": False,
    }
    if trajectory:
        best = min(trajectory, key=lambda row: row["energy_hartree"])
        old_energy = json.loads((work / "outputs/grid/center/result.json").read_text())[
            "energy_hartree"
        ]
        result["lowest_completed_SCF_center_energy_change_kcal_mol"] = (
            best["energy_hartree"] - old_energy
        ) * CONVERSION["energy_factor"]
        result["lowest_completed_SCF_geometry"] = best["coordinates"]
        result["lowest_completed_SCF_is_an_optimized_minimum_claim"] = False
    if r["status"] == "completed":
        x = np.array(r["geometry_bohr"]) * ENGINE_BOHR_ANGSTROM
        steps = re.findall(
            r"^\s+(\d+)\s+(-\d+\.\d+)\s+(.+?)\s+~\s*$", out, re.MULTILINE
        )
        require(
            r["scf_converged"] is True
            and r["units"] == UNITS
            and r["axis_order"] == "atom-major Cartesian xyz"
            and "Final optimized geometry and variables:" in out
            and steps
            and int(steps[-1][0]) <= 20
            and steps[-1][2].count("*") == 4
            and abs(energies[-1] - r["energy_hartree"]) < 2e-10
            and np.allclose(
                trajectory[0]["coordinates"], request["coordinates"], atol=1e-9, rtol=0
            )
            and np.allclose(trajectory[-1]["coordinates"], x, atol=1e-9, rtol=0),
            "optimization_native_completion_mutation",
        )
        result["final_native_convergence_row"] = {
            "step": int(steps[-1][0]),
            "raw_measures": steps[-1][2],
            "four_active_force_displacement_criteria_met": True,
            "units": "atomic units; Delta E criterion inactive (o)",
        }
        before = np.array(old["center_coordinates"])
        a = x - x.mean(axis=0)
        b = before - before.mean(axis=0)
        u, _sv, vt = np.linalg.svd(a.T @ b)
        rotation = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
        result["RMS_aligned_cartesian_change_A"] = float(
            np.sqrt(np.mean((a @ rotation - b) ** 2))
        )
        result["minimum_stability_claimed"] = False
        result["stability_note"] = (
            "optimization stationarity only; no new Hessian/local-minimum certification"
        )
        coords = [
            {"kind": "bond", "sites": [1, 2]},
            {"kind": "bond", "sites": [2, 3]},
            {"kind": "angle", "sites": [1, 2, 3]},
        ]
        old_result = json.loads((work / "outputs/grid/center/result.json").read_text())
        old_force = (
            -np.load(work / "outputs/grid/center/gradient.npy", allow_pickle=False)
            * CONVERSION["energy_factor"]
            / CONVERSION["distance_factor"]
        )
        new_gradient = np.load(
            work / "outputs/optimize/qz-dense/gradient.npy", allow_pickle=False
        )
        require(
            new_gradient.shape == (9, 3)
            and np.isfinite(new_gradient).all()
            and np.allclose(
                new_gradient,
                matrix_ascii(work / "outputs/optimize/qz-dense/gradient.dat", (9, 3)),
                atol=1e-15,
                rtol=0,
            ),
            "optimization_native_gradient_mutation",
        )
        new_force = (
            -new_gradient
            * CONVERSION["energy_factor"]
            / CONVERSION["distance_factor"]
        )
        result["same_protocol_center_energy_lowering_kcal_mol"] = (
            r["energy_hartree"] - old_result["energy_hartree"]
        ) * CONVERSION["energy_factor"]
        result["old_center_force_RMS_kcal_mol_A"] = float(
            np.sqrt(np.mean(old_force**2))
        )
        result["new_center_force_RMS_kcal_mol_A"] = float(
            np.sqrt(np.mean(new_force**2))
        )
        result["new_center_force_max_kcal_mol_A"] = float(np.max(np.abs(new_force)))
        result.update(
            {
                "new_center_identity": signed(
                    {
                        "parent": old["target"]["identity"],
                        "coordinates": x.tolist(),
                        "method": "PBE0",
                        "basis": "def2-QZVP",
                        "grid": [150, 974],
                        "SCF": "DF/def2-universal-jkfit",
                        "native_job_identity": r["job_identity"],
                    }
                )["identity"],
                "coordinates": x.tolist(),
                "internal_coordinate_changes": [
                    {
                        "coordinate": c,
                        "J58": measure(before, list(range(1, 10)), c),
                        "QZ": measure(x, list(range(1, 10)), c),
                    }
                    for c in coords
                ],
                "RMS_cartesian_change_A": float(np.sqrt(np.mean((x - before) ** 2))),
                "new_25_point_scan_executed": False,
                "policy": "center-only optimization; not mixed with old-center scan; revised scan requires separately bound target/geometry identity",
            }
        )
    return signed(result)


def revised_scans(repo, distribution, work):
    """Additive QZ-center targets/scans; generated coordinates, not new QM data."""
    from island.pcff_completion import digest, require
    from island.pcff_qm.geometry import geometries

    from .campaign import NAMES, selected

    center = center_report(work, repo)
    require(
        center["optimization_status"] == "completed",
        "optimized_QZ_center_not_available",
    )
    system, assessment, targets = selected(repo, distribution, center["coordinates"])
    results = {}
    for name, target in zip(NAMES, targets, strict=True):
        p = {k: v for k, v in target.payload.items() if k != "identity"}
        p["schema"] = "island.j59-QZ-center-target.v1"
        p["previous_target_identity"] = target.payload["identity"]
        p["center_record_identity"] = center["identity"]
        p["qm_protocol"] = {
            **p["qm_protocol"],
            "authority": "J59 proposed convergence study; not author/source authority",
            "basis_set": "def2-QZVP",
            "grid": [150, 974],
            "SCF_algorithm": "DF",
            "SCF_auxiliary_basis": "def2-universal-jkfit",
        }
        p["identity"] = digest(p)
        results[name] = signed(
            {
                "schema": "island.j59-QZ-center-scan.v1",
                "target": p,
                "center_record_identity": center["identity"],
                "points": geometries(system, assessment, p["scan"]),
                "QM_scan_executed": False,
                "fit_authorized": False,
                "historical_data_modified": False,
            }
        )
    return results
