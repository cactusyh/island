"""Independent scalar reconstruction of the J63 research equations (no Jets)."""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations

import numpy as np


def scalar(model, xyz):
    x = np.asarray(xyz, float)
    rows = model["interaction_inventory"]
    lookup = {r["id"]: r for r in rows}
    components = defaultdict(float)
    columns = defaultdict(float)

    @lru_cache(None)
    def bond(s):
        return np.linalg.norm(x[s[0] - 1] - x[s[1] - 1])

    @lru_cache(None)
    def angle(s):
        a, b, c = x[np.array(s) - 1]
        u = a - b
        v = c - b
        return np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))

    @lru_cache(None)
    def cosine(s):
        i, j, k, l = x[np.array(s) - 1]
        axis = k - j
        axis /= np.linalg.norm(axis)
        u = i - j
        v = l - k
        u -= u @ axis * axis
        v -= v @ axis * axis
        return u @ v / np.linalg.norm(u) / np.linalg.norm(v)

    def dr(s):
        s = tuple(s)
        key = "quartic_bond:" + ",".join(map(str, min(s, s[::-1])))
        return bond(s) - lookup[key]["assignment"]["values"][0]

    def dt(s):
        s = tuple(s)
        key = "quartic_angle:" + ",".join(map(str, min(s, s[::-1])))
        return angle(s) - lookup[key]["assignment"]["values"][0]

    for row in rows:
        f, s, a = row["family"], tuple(row["sites"]), row["assignment"]
        if f in ("atom_type", "mass", "charge", "nonbond(9-6)"):
            continue
        if a is not None:
            v = np.array(a["values"]) / 4.184
            if f == "quartic_bond":
                d = dr(s)
                e = sum(v[n] * d ** (n + 1) for n in (1, 2, 3))
            elif f == "quartic_angle":
                d = dt(s)
                e = sum(v[n] * d ** (n + 1) for n in (1, 2, 3))
            elif f == "bond-bond":
                e = v[0] * dr(s[:2]) * dr(s[1:])
            elif f == "bond-angle":
                e = (v[0] * dr(s[:2]) + v[1] * dr(s[1:])) * dt(s)
            elif f == "angle-angle":
                i, j, k, l = s
                e = v[0] * dt((i, j, k)) * dt((k, j, l))
            elif f == "torsion_3":
                c = cosine(s)
                phi = np.arccos(np.clip(c, -1, 1))
                e = sum(v[2 * (n - 1)] * (1 + np.cos(n * phi)) for n in (1, 2, 3))
            else:
                raise ValueError(f)
            components[f] += e
            continue
        labels = tuple(model["typing"][str(i)] for i in s)
        reverse = (
            (labels[3], labels[1], labels[2], labels[0])
            if f == "angle-angle"
            else labels[::-1]
        )
        if reverse < labels:
            labels = reverse
            s = (s[3], s[1], s[2], s[0]) if f == "angle-angle" else s[::-1]
        name = f + ":" + ",".join(labels)
        i, j, k, l = s
        if f == "angle-angle":
            columns[name] += dt((i, j, k)) * dt((k, j, l))
        elif f == "bond-bond_1_3":
            columns[name] += dr((i, j)) * dr((k, l))
        elif f == "angle-angle-torsion_1":
            columns[name] += dt((i, j, k)) * dt((j, k, l)) * cosine(s)
        else:
            ds = (
                [dr((i, j)), dr((k, l))]
                if f == "end_bond-torsion_3"
                else [dr((j, k))]
                if f == "middle_bond-torsion_3"
                else [dt((i, j, k)), dt((j, k, l))]
            )
            phi = np.arccos(np.clip(cosine(s), -1, 1))
            for b, d in enumerate(ds):
                for n in (1, 2, 3):
                    columns[name + f":block{b}:n{n}"] += d * np.cos(n * phi)
    adj = {i: set() for i in range(1, 28)}
    for i, j, _ in model["bonds"]:
        adj[i].add(j)
        adj[j].add(i)
    for i, j in combinations(range(1, 28), 2):
        if j in adj[i] or any(j in adj[k] for k in adj[i]):
            continue
        r = bond((i, j))
        ri, ei = lookup[f"nonbond(9-6):{i}"]["assignment"]["values"]
        rj, ej = lookup[f"nonbond(9-6):{j}"]["assignment"]["values"]
        rm6 = (ri**6 + rj**6) / 2
        eps = 2 * np.sqrt(ei * ej) / 4.184 * (ri * rj) ** 3 / (ri**6 + rj**6)
        components["LJ96_Coulomb"] += (
            eps * (2 * rm6**1.5 / r**9 - 3 * rm6 / r**6)
            + 332.06371 * model["charges"][str(i)] * model["charges"][str(j)] / r
        )
    return {
        "fixed_energy": float(sum(components.values())),
        "components": dict(components),
        "columns": dict(columns),
    }


def energy(model, xyz, coefficients):
    p = scalar(model, xyz)
    return p["fixed_energy"] + sum(v * coefficients[n] for n, v in p["columns"].items())
