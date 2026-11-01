"""Independent Cartesian/ASCII/decimal and unknown-null baseline reconstructions."""

import argparse
import json
import re
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed, write
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_independent import scalar
from island.pcff_qm_pilot.hmds_paths import fragment, path
from island.pcff_qm_pilot.hmds_precision import read


def independent_path(model, center, geometry, spec):
    """Cross products at actual endpoints, independent of path's sin/cos derivatives."""
    x, y = np.asarray(center), np.asarray(geometry)
    i, j = [n - 1 for n in spec["bond"]]
    moving = np.array(fragment(model, spec["bond"], spec["moving_root"])) - 1
    axis = x[j] - x[i]
    axis /= np.linalg.norm(axis)
    tangent, acceleration = np.zeros((27, 3)), np.zeros((27, 3))
    if spec["kind"] == "stretch":
        tangent[moving] = axis
    else:
        tangent[moving] = np.cross(axis, y[moving] - x[i])
        acceleration[moving] = np.cross(axis, tangent[moving])
    return tangent, acceleration


def ascii_gradient(directory):
    lines = (directory / "gradient.dat").read_text().splitlines()
    require(int(lines[2]) == 81, "independent81_gradient_required")
    entries = {}
    for line in lines[3:]:
        i, j, value = line.split()
        key = int(i), int(j)
        require(key not in entries, "independent_duplicate_gradient_component")
        entries[key] = float(value)
    require(set(entries) == set(np.ndindex(27, 3)), "independent_gradient_coverage")
    result = np.array([[entries[i, j] for j in range(3)] for i in range(27)])
    require(np.isfinite(result).all(), "independent_gradient_nonfinite")
    return result


def native_decimal_total(directory):
    """Separate line parser; do not reuse production decimal_energy()."""
    values = [
        line.split("=", 1)[1].strip()
        for line in (directory / "psi4.out").read_text().splitlines()
        if line.strip().startswith("Total Energy =")
    ]
    require(len(values) == 1, "independent_unique_total_energy_required")
    value = Decimal(values[0])
    require(value.is_finite(), "independent_native_energy_nonfinite")
    return value


def check(repo, root, work):
    model = load_model(repo)
    plan = read(root / "plan/scientific-plan.json")
    x = np.array(plan["center_coordinates"])
    ef = 4.3597447222060e-18 * 6.02214076e23 / 4184
    lf = 5.29177210544e-11 * 1e10
    r = read(root / "baseline-calibration.json")
    result = {}
    for name, spec in ((p["name"], p) for p in plan["paths"]):
        keys = (
            ["dense/plus", "dense/minus"]
            if name == "siloxane"
            else [f"dense/{name}-plus", f"dense/{name}-minus"]
            if spec["kind"] == "stretch"
            else []
        )
        if not keys:
            continue
        h = 0.02 if spec["kind"] == "rotation" else 0.003
        points, derivs = [], []
        for sign, key in zip((1, -1), keys, strict=True):
            geometry = read(work / "requests" / (key + ".json"))["coordinates"]
            tangent, acceleration = independent_path(model, x, geometry, spec)
            actual = path(model, x, spec, sign * h)
            require(
                np.max(abs(actual["tangent"] - tangent)) < 1e-12
                and np.max(abs(actual["acceleration"] - acceleration)) < 1e-12,
                "independent_path_derivative_mismatch",
            )
            # Direct ASCII native gradient contraction, then independent unit conversion.
            g = ascii_gradient(work / "outputs" / key)
            derivative = float(np.sum(g * tangent) * ef / lf)
            derivs.append(derivative)
            points.append(
                {
                    "key": key,
                    "tangent": tangent.tolist(),
                    "acceleration": acceleration.tolist(),
                    "path_energy_derivative": derivative,
                    "path_force_or_torque": -derivative,
                }
            )
        curvature = (derivs[0] - derivs[1]) / (2 * h)
        expected = (
            read(root / "dense-grid-decision.json")["response"]
            if name == "siloxane"
            else r["paths"][name]["protocol_responses"]["TZVP_dense_tight"]
        )
        require(
            abs(curvature - expected["curvature_from_path_gradient"]) < 1e-9,
            "independent_ASCII_path_response_mismatch",
        )
        totals = [
            native_decimal_total(work / "outputs" / key)
            for key in ["exact/center", *keys]
        ]
        with localcontext() as context:
            context.prec = 60
            # Independently derive Hartree -> kcal/mol, then convert the small
            # native second difference, without converting electronic offsets.
            conversion = Decimal("4.3597447222060e-18") * Decimal("6.02214076e23")
            conversion /= Decimal(4184)
            energy_curvature = (
                (totals[1] + totals[2] - 2 * totals[0])
                / Decimal(str(h)) ** 2
                * conversion
            )
        energy_error = float(energy_curvature) - expected["curvature_from_energies"]
        require(abs(energy_error) < 1e-9, "independent_decimal_curvature_mismatch")
        result[name] = {
            "points": points,
            "independent_curvature": curvature,
            "difference_from_production": curvature
            - expected["curvature_from_path_gradient"],
            "native_decimal_energies_hartree": list(map(str, totals)),
            "independent_energy_curvature": float(energy_curvature),
            "energy_curvature_difference_from_production": energy_error,
            "independent_energy_conversion_factor": str(conversion),
        }
    fixed_checks = {}
    for spec in plan["paths"]:
        if spec["kind"] != "stretch":
            continue
        name = spec["name"]
        rows = {}
        for h in spec["steps"]:
            positions = [path(model, x, spec, s * h)["coordinates"] for s in (0, 1, -1)]
            scalars = [scalar(model, y) for y in positions]
            fixed_curvature = (
                scalars[1]["fixed_energy"]
                - 2 * scalars[0]["fixed_energy"]
                + scalars[2]["fixed_energy"]
            ) / h**2
            nulls = {
                n: (
                    scalars[1]["columns"][n]
                    - 2 * scalars[0]["columns"][n]
                    + scalars[2]["columns"][n]
                )
                / h**2
                for n in scalars[0]["columns"]
            }
            require(
                max(abs(v) for v in nulls.values()) < 1e-7,
                "independent_unknown_null_failure",
            )
            # Independent explicit correction sum using raw Euclidean bond lengths.
            corrections = []
            for y in positions:
                value = 0.0
                for row in model["interaction_inventory"]:
                    if row["family"] != "quartic_bond":
                        continue
                    i, j = row["sites"]
                    types = {model["typing"][str(i)], model["typing"][str(j)]}
                    if types == ({"sio", "c3"} if name == "SiC" else {"c3", "hc"}):
                        value += (
                            np.linalg.norm(y[i - 1] - y[j - 1])
                            - row["assignment"]["values"][0]
                        ) ** 2
                corrections.append(value)
            sensitivity = (corrections[1] - 2 * corrections[0] + corrections[2]) / h**2
            require(
                abs(sensitivity - 2) < 1e-8, "independent_quadratic_prefactor_failure"
            )
            qm = (
                read(
                    repo / "docs/evidence/phase_4j63/deformation-reference-results.json"
                )["primary"][name]["steps"][str(h)]["curvature_from_energies"]
                if h == 0.006
                else r["paths"][name]["protocol_responses"]["TZVP_dense_tight"][
                    "curvature_from_energies"
                ]
            )
            rows[str(h)] = {
                "all_fixed_scalar_curvature": fixed_curvature,
                "all56_scalar_unknown_nulls": nulls,
                "independent_quadratic_sensitivity": sensitivity,
                "conditional_scalar_curvature_estimate": (qm - fixed_curvature)
                / sensitivity,
                "QM_scalar_energy_curvature": qm,
                "QM_protocol": "TZVP150x974_standard_SCF"
                if h == 0.006
                else "TZVP250x1454_tight_SCF",
                "scalar_vs_force_curvature_scope": "Finite-step energy and gradient curvatures differ by expected O(h^2) quartic/path truncation; no bare-source-bond substitution.",
            }
        fixed_checks[name] = rows
    # Native calculation success is independently asserted by the importer elsewhere;
    # this records no synthetic energy, force or physical coefficient.
    native = {}
    for key in read(work / "ledger.json")["calls"]:
        text = (work / "outputs" / key / "psi4.out").read_text()
        values = re.findall(r"@(?:DF-)?RKS Final Energy:\s*([-+\d.]+)", text)
        require(values, "independent_native_energy_missing")
        native[key] = {
            "printed_final_energy_hartree": values[-1],
            "independent_native_total_energy_hartree": str(
                native_decimal_total(work / "outputs" / key)
            ),
        }
    out = signed(
        {
            "paths": result,
            "unknown_null_fixed_baseline_checks": fixed_checks,
            "native_energy_strings": native,
            "physical56_coefficients_selected": False,
        }
    )
    write(root / "independent-path-and-baseline-checks.json", out)
    return out


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
            {
                "paths": list(
                    check(
                        args.repo.resolve(),
                        args.evidence.resolve(),
                        args.work.resolve(),
                    )["paths"]
                )
            }
        )
    )
