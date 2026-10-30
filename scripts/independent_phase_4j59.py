"""Independent scalar source equations and native numerical reconstruction.

Does not use the production residual evaluator or native importer. No new
physical coefficients: source inputs only; outputs are verification records.
"""

import argparse
import json
import re
from decimal import Decimal
from pathlib import Path

import numpy as np


def source_energy(rows, xyz):
    x = np.asarray(xyz, float)
    lookup = {r["id"]: r for r in rows}
    out = {}

    def length(ids):
        return float(np.sqrt(np.sum((x[ids[0] - 1] - x[ids[1] - 1]) ** 2)))

    def bend(ids):
        a, b, c = [x[i - 1] for i in ids]
        u = a - b
        v = c - b
        return float(np.arctan2(np.linalg.norm(np.cross(u, v)), u @ v))

    for r in rows:
        f = r["family"]
        a = r["assignment"]
        ids = r["sites"]
        if f not in ("quartic_bond", "quartic_angle", "bond-bond", "bond-angle"):
            continue
        v = a["values"]
        if f in ("quartic_bond", "quartic_angle"):
            d = (length(ids) if f == "quartic_bond" else bend(ids)) - v[0]
            e = d * d * (v[1] + d * (v[2] + d * v[3]))
        else:
            i, j, k = ids
            deps = r["dependencies"]
            dl = length([i, j]) - lookup[deps[0]]["assignment"]["values"][0]
            dr = length([j, k]) - lookup[deps[1]]["assignment"]["values"][0]
            e = (
                v[0] * dl * dr
                if f == "bond-bond"
                else (v[0] * dl + v[1] * dr)
                * (bend(ids) - lookup[deps[2]]["assignment"]["values"][0])
            )
        out[r["id"]] = e / 4.184
    return out


def numerical_force(rows, xyz, step):
    x = np.array(xyz, float)
    f = np.zeros_like(x)
    for i in range(len(x)):
        for j in range(3):
            a = x.copy()
            b = x.copy()
            a[i, j] += step
            b[i, j] -= step
            f[i, j] = -(
                sum(source_energy(rows, a).values())
                - sum(source_energy(rows, b).values())
            ) / (2 * step)
    return f


def native_values(directory):
    path = Path(directory)
    out = (path / "psi4.out").read_text()
    energy = Decimal(re.findall(r"@(?:DF-)?RKS Final Energy:\s*([-+\d.]+)", out)[-1])
    lines = (path / "gradient.dat").read_text().splitlines()
    if lines[1] != "symmetry 0" or int(lines[2]) != 27:
        raise ValueError("independent_matrix_shape")
    g = np.zeros((9, 3))
    seen = set()
    eh = Decimal("4.3597447222060e-18") * Decimal("6.02214076e23") / Decimal(4184)
    bohr = Decimal("5.29177210544e-11") * Decimal("1e10")
    for line in lines[3:]:
        if not line.strip():
            continue
        i, j, v = line.split()
        i, j = int(i), int(j)
        if (i, j) in seen:
            raise ValueError("duplicate_native_axis")
        seen.add((i, j))
        g[i, j] = float(-Decimal(v) * eh / bohr)
    if len(seen) != 27:
        raise ValueError("independent_matrix_coverage")
    return float(energy * eh), g


def check(work, normalized):
    work = Path(work)
    data = json.loads(Path(normalized).read_text())
    records = []
    for key, sample in data["samples"].items():
        energy, forces = native_values(work / "outputs" / key)
        error = float(np.max(np.abs(forces - np.array(sample["forces"]))))
        de = abs(energy - sample["energy"])
        if de > 2e-7 or error > 1e-10:
            raise ValueError("independent_unit_disagreement")
        records.append(
            {
                "job": key,
                "energy_error_kcal_mol": de,
                "force_max_error_kcal_mol_A": error,
            }
        )
    return {
        "production_importer_used": False,
        "formula": "Decimal(SI hartree*Avogadro/4184); force=-gradient*(energy_factor/Bohr_A)",
        "samples": records,
        "passed": True,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--normalized", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.write_text(
        json.dumps(check(a.work, a.normalized), indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
