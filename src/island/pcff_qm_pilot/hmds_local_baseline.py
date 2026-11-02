"""Localized scalar-potential hypotheses and independent fixed-family audit.

These additive polynomials do not propagate a changed source equilibrium into
Class-II dependencies. They never alter charges, rows, or automatic assignment.
"""

from collections import defaultdict
from itertools import combinations

import numpy as np

from island.pcff_completion import require

from .hamiltonian import COULOMB, Jet, bend, cos_phi, harmonic, length
from .hmds_hamiltonian import parts
from .hmds_paths import adjacency
from .molecular_control import cartesian
from .residual import angle, bond

SITES = {"bridge_angle": [(2, 5, 6)], "SiO_stretch": [(2, 5), (5, 6)]}


def references(model, variant):
    require(variant in SITES, "J68_unknown_correction")
    family = "quartic_angle" if variant == "bridge_angle" else "quartic_bond"
    rows = {r["id"]: r for r in model["interaction_inventory"]}
    return [
        {
            "sites": list(s),
            "source_row": family + ":" + ",".join(map(str, min(s, s[::-1]))),
            "equilibrium": rows[family + ":" + ",".join(map(str, min(s, s[::-1])))][
                "assignment"
            ]["values"][0],
            "types": [model["typing"][str(i)] for i in s],
        }
        for s in SITES[variant]
    ]


def correction(model, xyz, variant):
    """E=a sum(q-q0)+b sum(q-q0)^2; F=-dE/dCartesian, no 1/2."""
    x = np.asarray(xyz, float)
    require(x.shape == (27, 3) and np.isfinite(x).all(), "J68_correction_coordinates")
    e, f = np.zeros(2), np.zeros((27, 3, 2))
    for row in references(model, variant):
        value, gradient = (angle if variant == "bridge_angle" else bond)(
            x, row["sites"]
        )
        delta = value - row["equilibrium"]
        e += [delta, delta**2]
        f[:, :, 0] -= gradient
        f[:, :, 1] -= 2 * delta * gradient
    return e, f


def scalar(model, xyz, variant):
    """Separate scalar equation for central finite-difference derivative checks."""
    x = np.asarray(xyz)
    q = []
    for row in references(model, variant):
        sites = [i - 1 for i in row["sites"]]
        if variant == "bridge_angle":
            u, v = x[sites[0]] - x[sites[1]], x[sites[2]] - x[sites[1]]
            value = np.arccos(
                np.clip(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1)
            )
        else:
            value = np.linalg.norm(x[sites[0]] - x[sites[1]])
        q.append(value - row["equilibrium"])
    q = np.array(q)
    return np.array([sum(q), q @ q])


def derivative_check(model, xyz, variant, step=1e-5):
    e, f = correction(model, xyz, variant)
    numerical = np.zeros_like(f)
    x = np.array(xyz, float)
    for i in range(27):
        for j in range(3):
            plus, minus = x.copy(), x.copy()
            plus[i, j] += step
            minus[i, j] -= step
            numerical[i, j] = -(
                scalar(model, plus, variant) - scalar(model, minus, variant)
            ) / (2 * step)
    jets = cartesian(xyz)
    features = [Jet(0, np.zeros(81)), Jet(0, np.zeros(81))]
    for row in references(model, variant):
        d = (bend if variant == "bridge_angle" else length)(jets, row["sites"]) - row[
            "equilibrium"
        ]
        features[0] += d
        features[1] += d * d
    return {
        "scalar_energy_error": float(np.max(abs(e - scalar(model, xyz, variant)))),
        "force_finite_difference_error": float(np.max(abs(f - numerical))),
        "independent_Jet_energy_error": float(np.max(abs(e - [j.v for j in features]))),
        "independent_Jet_force_error": float(
            np.max(abs(f - np.stack([-j.d.reshape(27, 3) for j in features], axis=-1)))
        ),
        "finite_difference_step_angstrom": step,
    }


def extend(model, rows, data, variant):
    kernels = [correction(model, r["coordinates"], variant) for r in rows]
    return {
        **data,
        "names": data["names"] + [variant + ":linear", variant + ":quadratic"],
        "energy": np.c_[data["energy"], np.array([k[0] for k in kernels])],
        "force": np.concatenate(
            [data["force"], np.array([k[1] for k in kernels])], axis=-1
        ),
    }


def span_sensitivity(problem, amended):
    scale = np.maximum(np.linalg.norm(problem["w"], axis=0), 1e-30)
    base = problem["w"] / scale
    # Added physical columns precede the same energy intercept.
    c = amended["w"][:, -3:-1]
    cn = np.linalg.norm(c, axis=0)
    if np.any(cn < 1e-15):
        return {
            "independent": False,
            "reason": "zero added column",
            "added_column_norms": cn.tolist(),
        }
    normalized = c / cn
    off = normalized - base @ np.linalg.lstsq(base, normalized, rcond=None)[0]
    singular = np.linalg.svd(off, compute_uv=False)
    residual = problem["y"] - base @ np.linalg.lstsq(base, problem["y"], rcond=None)[0]
    return {
        "added_column_norms": cn.tolist(),
        "projected_columns": off.tolist(),
        "outside_old_span_norms": np.linalg.norm(off, axis=0).tolist(),
        "projected_singular_values": singular.tolist(),
        "independent": bool(singular[-1] > 1e-7),
        "target_projection_dot_products": (off.T @ residual).tolist(),
        "target_projected_capture_norm": float(
            np.linalg.norm(off @ np.linalg.lstsq(off, residual, rcond=None)[0])
        ),
    }


def fixed_families(model, xyz):
    """Independent Jet decomposition using actual source rows and pair conventions."""
    x = cartesian(xyz)
    rows = model["interaction_inventory"]
    lookup = {r["id"]: r for r in rows}
    zero = lambda: Jet(0, np.zeros(81))
    families = defaultdict(zero)

    def delta(f, s):
        key = f + ":" + ",".join(map(str, min(tuple(s), tuple(s[::-1]))))
        return (length if f == "quartic_bond" else bend)(x, s) - lookup[key][
            "assignment"
        ]["values"][0]

    def label(r):
        f, s = r["family"], r["sites"]
        if f == "quartic_bond" and set(s) in ({2, 5}, {5, 6}):
            return "quartic_bond:SiO"
        if f == "quartic_angle" and s == [2, 5, 6]:
            return "quartic_angle:SiOSi"
        return f

    for r in rows:
        f, s, a = r["family"], r["sites"], r["assignment"]
        if a is None:
            continue
        v = a.get("values")
        if f in ("quartic_bond", "quartic_angle"):
            q = delta(f, s)
            families[label(r)] += sum(v[k] / 4.184 * q ** (k + 1) for k in (1, 2, 3))
        elif f in ("bond-bond", "bond-angle"):
            i, j, k = s
            dl, dr = delta("quartic_bond", [i, j]), delta("quartic_bond", [j, k])
            families[f] += (
                v[0] * dl * dr
                if f == "bond-bond"
                else (v[0] * dl + v[1] * dr) * delta("quartic_angle", s)
            ) / 4.184
        elif f == "angle-angle":
            i, j, k, l = s
            families[f] += (
                v[0]
                / 4.184
                * delta("quartic_angle", [i, j, k])
                * delta("quartic_angle", [k, j, l])
            )
        elif f == "torsion_3":
            c = cos_phi(x, s)
            families[f] += sum(
                v[2 * (n - 1)] / 4.184 * (1 + harmonic(c, n)) for n in (1, 2, 3)
            )
    adj = adjacency(model)
    for i, j in combinations(model["atom_order"], 2):
        excluded = j in adj[i] or any(j in adj[k] for k in adj[i])
        if excluded:
            continue
        r = length(x, [i, j])
        ri, ei = lookup[f"nonbond(9-6):{i}"]["assignment"]["values"]
        rj, ej = lookup[f"nonbond(9-6):{j}"]["assignment"]["values"]
        rm = ((ri**6 + rj**6) / 2) ** (1 / 6)
        ep = 2 * np.sqrt(ei * ej) / 4.184 * ri**3 * rj**3 / (ri**6 + rj**6)
        families["LJ96"] += ep * (2 * (rm / r) ** 9 - 3 * (rm / r) ** 6)
        families["Coulomb"] += (
            COULOMB * model["charges"][str(i)] * model["charges"][str(j)] / r
        )
    out = {
        k: {"energy": j.v, "forces": (-j.d.reshape(27, 3)).tolist()}
        for k, j in families.items()
    }
    reference = parts(model, xyz)
    esum = sum(v["energy"] for v in out.values())
    fsum = sum(np.array(v["forces"]) for v in out.values())
    require(
        abs(esum - reference["fixed_energy"]) < 1e-9
        and np.max(abs(fsum - reference["fixed_forces"])) < 1e-9,
        "J68_fixed_family_closure",
    )
    return {
        "families": out,
        "energy_sum_error": float(abs(esum - reference["fixed_energy"])),
        "force_sum_error": float(np.max(abs(fsum - reference["fixed_forces"]))),
    }


def jet_correction(model, xyz, variant):
    """Independent automatic derivatives for complete amended-MM verification."""
    x = cartesian(xyz)
    out = [Jet(0, np.zeros(81)), Jet(0, np.zeros(81))]
    for row in references(model, variant):
        d = (bend if variant == "bridge_angle" else length)(x, row["sites"]) - row[
            "equilibrium"
        ]
        out[0] += d
        out[1] += d * d
    return np.array([j.v for j in out]), np.stack(
        [-j.d.reshape(27, 3) for j in out], axis=-1
    )
