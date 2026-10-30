"""Independent scalar Hamiltonian, normals/atan2 geometry and FD controls.

This module imports neither the J60 Hamiltonian nor its Cartesian kernels.
Synthetic coefficients used in controls are not source/fitted coefficients.
"""

from collections import defaultdict

import numpy as np

from scripts.independent_phase_4j59 import source_energy


def scalar(model, xyz):
    x = np.asarray(xyz, float)
    rows = model["interaction_inventory"]
    lookup = {r["id"]: r for r in rows}

    def r(s):
        return float(np.linalg.norm(x[s[0] - 1] - x[s[1] - 1]))

    def t(s):
        a, b, c = [x[i - 1] for i in s]
        u = a - b
        v = c - b
        return float(np.arctan2(np.linalg.norm(np.cross(u, v)), u @ v))

    def phi(s):
        a, b, c, d = [x[i - 1] for i in s]
        axis = c - b
        axis /= np.linalg.norm(axis)
        u = a - b
        v = d - c
        u = u - (u @ axis) * axis
        v = v - (v @ axis) * axis
        return float(np.arctan2(np.cross(axis, u) @ v, u @ v))

    def delta(f, s):
        k = f + ":" + ",".join(map(str, min(tuple(s), tuple(s[::-1]))))
        return (r(s) if f == "quartic_bond" else t(s)) - lookup[k]["assignment"][
            "values"
        ][0]

    def dr(s):
        return delta("quartic_bond", s)

    def dt(s):
        return delta("quartic_angle", s)

    energy = sum(source_energy(rows, x).values())
    columns = defaultdict(float)
    for row in rows:
        f = row["family"]
        sites = list(row["sites"])
        a = row["assignment"]
        if f == "torsion_3":
            v = a["values"]
            p = phi(sites)
            energy += sum(
                v[2 * n - 2] * (1 + np.cos(n * p - v[2 * n - 1])) / 4.184
                for n in (1, 2, 3)
            )
            continue
        if a is not None:
            continue
        labels = tuple(model["typing"][str(i)] for i in sites)
        if f == "angle-angle":
            key = (
                f
                + ":"
                + ",".join(min(labels, (labels[3], labels[1], labels[2], labels[0])))
            )
            i, j, k, l = sites
            columns[key] += dt([i, j, k]) * dt([k, j, l])
            continue
        if labels[::-1] < labels:
            sites = sites[::-1]
            labels = labels[::-1]
        i, j, k, l = sites
        key = f + ":" + ",".join(labels)
        p = phi(sites)
        if f == "bond-bond_1_3":
            columns[key] += dr([i, j]) * dr([k, l])
        elif f == "angle-angle-torsion_1":
            columns[key] += dt([i, j, k]) * dt([j, k, l]) * np.cos(p)
        else:
            ds = (
                [dr([i, j]), dr([k, l])]
                if f == "end_bond-torsion_3"
                else [dr([j, k])]
                if f == "middle_bond-torsion_3"
                else [dt([i, j, k]), dt([j, k, l])]
            )
            for block, d in enumerate(ds):
                for n in (1, 2, 3):
                    columns[key + f":block{block}:n{n}"] += d * np.cos(n * p)
    for pair in model["table"]:
        if pair["family"] != "LJ9-6+Coulomb" or pair["topological_distance"] < 3:
            continue
        c = pair["coefficient"]
        distance = r(pair["sites"])
        ri, ei = c["ri_ei"]
        rj, ej = c["rj_ej"]
        sig = ((ri**6 + rj**6) / 2) ** (1 / 6)
        eps = 2 * np.sqrt(ei * ej) / 4.184 * ri**3 * rj**3 / (ri**6 + rj**6)
        energy += (
            eps * (2 * (sig / distance) ** 9 - 3 * (sig / distance) ** 6)
            + 332.06371 * float(c["qi"]) * float(c["qj"]) / distance
        )
    return energy, dict(columns)


def check(model, xyz):
    from island.pcff_qm_pilot.hamiltonian import parts

    p = parts(model, xyz)
    e, c = scalar(model, xyz)
    names = p["columns"]
    require = lambda ok: (
        None if ok else (_ for _ in ()).throw(ValueError("independent_J60_mismatch"))
    )
    energy_error = abs(e - p["fixed_energy"])
    column_error = max(abs(c[k] - p["energy_jacobian"][i]) for i, k in enumerate(names))
    require(energy_error < 1e-10 and column_error < 1e-12)
    errors = []
    for h in (1e-4, 5e-5):
        x = np.array(xyz, float)
        fixed = np.zeros_like(x)
        force = np.zeros((9, 3, len(names)))
        for i in range(9):
            for j in range(3):
                a = x.copy()
                b = x.copy()
                a[i, j] += h
                b[i, j] -= h
                ea, ca = scalar(model, a)
                eb, cb = scalar(model, b)
                fixed[i, j] = -(ea - eb) / (2 * h)
                force[i, j] = [-(ca[k] - cb[k]) / (2 * h) for k in names]
        ef = float(np.max(abs(fixed - p["fixed_forces"])))
        ec = float(np.max(abs(force - p["force_jacobian"])))
        require(ef < 2e-5 and ec < 2e-6)
        errors.append(
            {
                "h_angstrom": h,
                "fixed_force_max_error": ef,
                "all_column_force_max_error": ec,
            }
        )
    return {
        "independent_fixed_energy_error": energy_error,
        "independent_columns_max_error": column_error,
        "finite_differences": errors,
    }
