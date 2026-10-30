"""Declared J60 isolated Class-II research Hamiltonian, not a PCFF profile.

All fixed values come from J59's pinned source rows. Every missing coefficient
is an explicit linear variable; omission is an error in a complete evaluation.
Cartesian derivatives use forward differentiation, not numerical differencing.
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed

from .residual import FIXED, fixed_rows, validate_model

IDENTITY = "island_research_class2_sih_disiloxane_j60_v1"
SOURCE_MODEL_IDENTITY = (
    "daea1878c79887938aa0016e7187d91aa7106b38b5c17467ea5bc4ea954c4ed2"
)
COULOMB = 332.06371  # pinned LAMMPS real unit implementation, dielectric=1


class Jet:
    """Scalar with exact first-order Cartesian chain rule."""

    def __init__(self, value, derivative):
        self.v = float(value)
        self.d = np.asarray(derivative, float)

    def __add__(self, other):
        if isinstance(other, Jet):
            return Jet(self.v + other.v, self.d + other.d)
        return Jet(self.v + other, self.d)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, Jet):
            return Jet(self.v * other.v, self.d * other.v + other.d * self.v)
        return Jet(self.v * other, self.d * other)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        return Jet(self.v**exponent, exponent * self.v ** (exponent - 1) * self.d)

    def __truediv__(self, other):
        return self * (other**-1 if isinstance(other, Jet) else 1 / other)

    def __rtruediv__(self, other):
        return self**-1 * other


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def subtract(a, b):
    return [x - y for x, y in zip(a, b, strict=True)]


def norm(a):
    v = dot(a, a)
    require(v.v > 1e-20, "singular_internal_coordinate")
    return v**0.5


def geometry(xyz):
    x = np.asarray(xyz, float)
    require(x.shape == (9, 3) and np.isfinite(x).all(), "geometry_mapping_mutation")
    eye = np.eye(x.size)
    return [[Jet(x[i, j], eye[3 * i + j]) for j in range(3)] for i in range(9)]


def length(x, sites):
    a, b = sites
    return norm(subtract(x[a - 1], x[b - 1]))


def bend(x, sites):
    a, b, c = sites
    u, v = subtract(x[a - 1], x[b - 1]), subtract(x[c - 1], x[b - 1])
    z = dot(u, v) / norm(u) / norm(v)
    require(abs(z.v) < 1 - 1e-12, "singular_angle")
    return Jet(np.arccos(z.v), -z.d / np.sqrt(1 - z.v * z.v))


def cos_phi(x, sites):
    """LAMMPS class2 cis=0 cosine; defined through projected end bonds.

    Signed phi is atan2(axis dot (end1 cross end2), end1 dot end2).
    LAMMPS source flips the sign when (end1 cross axis) dot end2 >0.
    Our energy kernels are even in phi (source phases here are zero).
    """
    i, j, k, l = sites
    axis = subtract(x[k - 1], x[j - 1])
    axis = [v / norm(axis) for v in axis]
    a, b = subtract(x[i - 1], x[j - 1]), subtract(x[l - 1], x[k - 1])
    p = subtract(a, [dot(a, axis) * v for v in axis])
    q = subtract(b, [dot(b, axis) * v for v in axis])
    return dot(p, q) / norm(p) / norm(q)


def harmonic(c, n):
    return {1: lambda: c, 2: lambda: 2 * c * c - 1, 3: lambda: 4 * c * c * c - 3 * c}[
        n
    ]()


def canonical(row, types):
    sites = list(row["sites"])
    labels = tuple(types[str(i)] for i in sites)
    if row["family"] == "angle-angle":
        reverse = (labels[3], labels[1], labels[2], labels[0])
        return sites, min(labels, reverse)
    if labels[::-1] < labels:
        return sites[::-1], labels[::-1]
    return sites, labels


def parts(model, xyz):
    """Fixed E/F and full unknown energy/force Jacobian; no missing value filled."""
    validate_model(model, SOURCE_MODEL_IDENTITY)
    x = geometry(xyz)
    rows = model["interaction_inventory"]
    lookup = {r["id"]: r for r in rows}
    types = model["typing"]

    # Source equilibrium values are distances/radians, not energy converted.
    def delta(f, sites):
        key = f + ":" + ",".join(map(str, min(tuple(sites), tuple(sites[::-1]))))
        eq = lookup[key]["assignment"]["values"][0]
        return (length(x, sites) if f == "quartic_bond" else bend(x, sites)) - eq

    def dr(s):
        return delta("quartic_bond", s)

    def dt(s):
        return delta("quartic_angle", s)

    fixed = fixed_rows(rows, xyz)
    fixed_energy = fixed["energy"]
    fixed_force = np.array(fixed["forces"])
    components = dict(fixed["components"])
    columns = defaultdict(lambda: Jet(0, np.zeros(27)))
    torsion = Jet(0, np.zeros(27))
    for row in rows:
        f = row["family"]
        s = row["sites"]
        a = row["assignment"]
        if f == "torsion_3":
            c = cos_phi(x, s)
            values = a["values"]
            require(
                all(abs(values[i]) < 1e-14 for i in (1, 3, 5)),
                "unsupported_source_phase",
            )
            torsion += sum(
                values[2 * (n - 1)] / 4.184 * (1 + harmonic(c, n)) for n in (1, 2, 3)
            )
        elif f == "angle-angle":
            require(a is None, "unexpected_new_source_row")
            s, labels = canonical(row, types)
            i, j, k, l = s
            columns[f + ":" + ",".join(labels)] += dt([i, j, k]) * dt([k, j, l])
        elif f in (
            "bond-bond_1_3",
            "end_bond-torsion_3",
            "middle_bond-torsion_3",
            "angle-torsion_3",
            "angle-angle-torsion_1",
        ):
            require(a is None, "unexpected_new_source_row")
            s, labels = canonical(row, types)
            i, j, k, l = s
            key = f + ":" + ",".join(labels)
            if f == "bond-bond_1_3":
                columns[key] += dr([i, j]) * dr([k, l])
            elif f == "angle-angle-torsion_1":
                columns[key] += dt([i, j, k]) * dt([j, k, l]) * cos_phi(x, s)
            else:
                ds = (
                    [dr([i, j]), dr([k, l])]
                    if f == "end_bond-torsion_3"
                    else [dr([j, k])]
                    if f == "middle_bond-torsion_3"
                    else [dt([i, j, k]), dt([j, k, l])]
                )
                c = cos_phi(x, s)
                for block, d in enumerate(ds):
                    for n in (1, 2, 3):
                        columns[key + f":block{block}:n{n}"] += d * harmonic(c, n)
    components["source_phase_transformed_main_torsion"] = torsion.v
    fixed_energy += torsion.v
    fixed_force -= torsion.d.reshape(9, 3)
    pair = Jet(0, np.zeros(27))
    pair_rows = []
    for row in model["table"]:
        if row["family"] != "LJ9-6+Coulomb":
            continue
        distance = row["topological_distance"]
        weight = 0 if distance in (1, 2) else 1
        c = row["coefficient"]
        r = length(x, row["sites"])
        ri, ei = c["ri_ei"]
        rj, ej = c["rj_ej"]
        ei /= 4.184
        ej /= 4.184
        rm = ((ri**6 + rj**6) / 2) ** (1 / 6)
        ep = 2 * np.sqrt(ei * ej) * ri**3 * rj**3 / (ri**6 + rj**6)
        require(r.v < 100, "isolated_cutoff_exceeded")
        energy = weight * (
            ep * (2 * (rm / r) ** 9 - 3 * (rm / r) ** 6)
            + COULOMB * float(c["qi"]) * float(c["qj"]) / r
        )
        pair += energy
        pair_rows.append(
            {
                "sites": row["sites"],
                "graph_distance": distance,
                "weight": weight,
                "energy": energy.v,
            }
        )
    components["isolated_LJ96_Coulomb"] = pair.v
    fixed_energy += pair.v
    fixed_force -= pair.d.reshape(9, 3)
    names = sorted(columns)
    return {
        "fixed_energy": fixed_energy,
        "fixed_forces": fixed_force,
        "components": components,
        "columns": names,
        "energy_jacobian": np.array([columns[k].v for k in names]),
        "force_jacobian": np.stack(
            [-columns[k].d.reshape(9, 3) for k in names], axis=-1
        ),
        "pairs": pair_rows,
        "accounting_complete": True,
        "coefficients_supplied": False,
    }


def evaluate(model, xyz, coefficients=None):
    p = parts(model, xyz)
    require(
        coefficients is not None and set(coefficients) == set(p["columns"]),
        "explicit_all_unknown_coefficients_required",
        required=p["columns"],
    )
    v = np.array([coefficients[k] for k in p["columns"]], float)
    require(np.isfinite(v).all(), "nonfinite_parameter")
    return {
        "energy": float(p["fixed_energy"] + p["energy_jacobian"] @ v),
        "forces": p["fixed_forces"] + p["force_jacobian"] @ v,
        "hamiltonian": IDENTITY,
        "parameter_origin": "research_extension_candidate",
    }


def accounting(model):
    """Every original J59 record receives an explicit selected-model treatment."""
    table = []
    for r in model["table"]:
        f = r["family"]
        if f == "torsion-torsion_1":
            status = "model_inapplicable"
            reason = "Not in declared LAMMPS class2 sum; research scope differs from literal FRC section. Not an assigned zero."
        elif (
            f == "wilson_out_of_plane" or r.get("id") == "improper_domain_applicability"
        ):
            status = "structurally_inapplicable"
            reason = "Pinned msi2lmp generates Wilson only at degree-three centers; this graph has degrees 1,2,4. All degree-four AA terms retained."
        elif f in FIXED or f in (
            "atom_type",
            "mass",
            "charge",
            "nonbond(9-6)",
            "LJ9-6+Coulomb",
            "torsion_3",
        ):
            status = "verified_fixed"
            reason = "Source row plus explicit J60 executable convention."
        else:
            status = "explicit_fitted_unknown"
            reason = (
                "Missing numerical coefficient is a joint-model variable, not zero."
            )
        table.append({**r, "J60_status": status, "J60_reason": reason})
    return signed(
        {
            "schema": "island.j60-hamiltonian-accounting.v1",
            "internal_identity": IDENTITY,
            "parent_J59_model_identity": model["identity"],
            "table": table,
            "complete_functional_accounting": True,
            "numerical_parameter_completeness": "incomplete",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def load_model(repo):
    return json.loads(
        (Path(repo) / "docs/evidence/phase_4j59/interaction-table.json").read_text()
    )
