"""Validation-only Class II energy kernels; no production force evaluator.

Small immutable terms never retain an FRC catalog. Forces below are central
finite differences for verification, not a dynamics implementation.
"""

from dataclasses import dataclass
from math import acos, asin, atan2, cos, isfinite, sqrt

import numpy as np

from .source import boundary, require

# sites, coefficients, equilibrium dependencies, in H3 physical role order.
SHAPES = {
    "quartic_bond": (2, 4, 0),
    "quartic_angle": (3, 4, 0),
    "bond-bond": (3, 1, 2),
    "bond-angle": (3, 2, 3),
    "torsion_3": (4, 6, 0),
    "middle_bond-torsion_3": (4, 3, 1),
    "end_bond-torsion_3": (4, 6, 2),
    "angle-torsion_3": (4, 6, 2),
    "angle-angle-torsion_1": (4, 1, 2),
    "bond-bond_1_3": (4, 1, 2),
    "angle-angle": (4, 1, 2),
}


def numeric(values):
    require(
        all(type(v) in (int, float) and isfinite(v) for v in values),
        "Finite nonboolean numbers required",
    )


def norm(v):
    r = float(np.linalg.norm(v))
    require(isfinite(r) and r > 1e-12, "Singular/overflow geometry")
    return r


def angle(a, b, c):
    u = a - b
    v = c - b
    cosine = float(np.dot(u, v) / (norm(u) * norm(v)))
    require(isfinite(cosine) and abs(cosine) < 1 - 1e-14, "Singular angle")
    return acos(cosine)


def torsion(x):
    t = x[2] - x[1]
    t = t / norm(t)
    v1 = x[0] - x[1]
    v3 = x[3] - x[2]
    a = v1 - np.dot(v1, t) * t
    b = v3 - np.dot(v3, t) * t
    norm(a)
    norm(b)
    return atan2(float(np.dot(np.cross(t, a), b)), float(np.dot(a, b)))


@dataclass(frozen=True)
class Class2Term:
    family: str
    sites: tuple
    coefficients: tuple
    equilibria: tuple = ()

    def _shapes(self):
        return SHAPES

    @boundary
    def __post_init__(self):
        require(
            self.family in self._shapes(), "Unsupported Class II term/applicability"
        )
        require(
            all(
                type(x) is tuple
                for x in (self.sites, self.coefficients, self.equilibria)
            ),
            "Owned tuple term data required",
        )
        require(
            tuple(map(len, (self.sites, self.coefficients, self.equilibria)))
            == self._shapes()[self.family],
            "Incorrect term dimensions",
        )
        require(
            all(type(i) is int for i in self.sites)
            and len(set(self.sites)) == len(self.sites),
            "Unique integer sites required",
        )
        numeric(self.coefficients + self.equilibria)

    @boundary
    def energy(self, coordinates, *, coordinate_unit="angstrom"):
        require(coordinate_unit == "angstrom", "Explicit angstrom coordinates required")
        require(
            type(coordinates) is dict
            and all(type(i) is int for i in coordinates)
            and set(coordinates) == set(self.sites),
            "Exact term coordinate coverage required",
        )
        raw = [coordinates[i] for i in self.sites]
        require(
            all(
                len(v) == 3 and all(not isinstance(k, (bool, np.bool_)) for k in v)
                for v in raw
            ),
            "Invalid coordinates",
        )
        x = np.array(raw, dtype=float)
        require(
            x.shape == (len(self.sites), 3) and bool(np.isfinite(x).all()),
            "Invalid coordinates",
        )
        with np.errstate(all="raise"):
            e = self._energy(x)
        require(isfinite(e), "Nonfinite term energy")
        return float(e)

    def _energy(self, x):
        p = self.coefficients
        q = self.equilibria
        f = self.family
        if f == "wilson_out_of_plane":
            arms = [x[k] - x[1] for k in (0, 2, 3)]
            unit = [v / norm(v) for v in arms]
            signed = float(np.dot(unit[0], np.cross(unit[1], unit[2])))
            chi = (
                sum(
                    asin(signed / norm(np.cross(unit[i], unit[j])))
                    for i, j in ((1, 2), (2, 0), (0, 1))
                )
                / 3
            )
            return p[0] * (chi - p[1]) ** 2
        if f == "quartic_bond":
            d = norm(x[0] - x[1]) - p[0]
            return sum(p[n - 1] * d**n for n in (2, 3, 4))
        if f == "quartic_angle":
            d = angle(*x) - p[0]
            return sum(p[n - 1] * d**n for n in (2, 3, 4))
        if f == "angle-angle":
            return (
                p[0]
                * (angle(x[0], x[1], x[2]) - q[0])
                * (angle(x[2], x[1], x[3]) - q[1])
            )
        if f in ("bond-bond", "bond-angle"):
            dl = norm(x[0] - x[1]) - q[0]
            dr = norm(x[1] - x[2]) - q[1]
            if f == "bond-bond":
                return p[0] * dl * dr
            return (p[0] * dl + p[1] * dr) * (angle(*x) - q[2])
        if f == "bond-bond_1_3":
            return p[0] * (norm(x[0] - x[1]) - q[0]) * (norm(x[2] - x[3]) - q[1])
        phi = torsion(x)
        if f == "torsion_3":
            return sum(
                p[2 * n - 2] * (1 - cos(n * phi - p[2 * n - 1])) for n in (1, 2, 3)
            )
        basis = tuple(cos(n * phi) for n in (1, 2, 3))
        left = sum(a * b for a, b in zip(p[:3], basis))
        if f == "middle_bond-torsion_3":
            return (norm(x[1] - x[2]) - q[0]) * left
        if f == "end_bond-torsion_3":
            return (norm(x[0] - x[1]) - q[0]) * left + (norm(x[2] - x[3]) - q[1]) * sum(
                a * b for a, b in zip(p[3:], basis)
            )
        dl = angle(*x[:3]) - q[0]
        dr = angle(*x[1:]) - q[1]
        if f == "angle-angle-torsion_1":
            return p[0] * dl * dr * cos(phi)
        return dl * left + dr * sum(a * b for a, b in zip(p[3:], basis))


class SourceClass2Term(Class2Term):
    """Expanded kernel contract; the historical Class2Term domain is unchanged."""

    def _shapes(self):
        return {**SHAPES, "wilson_out_of_plane": (4, 2, 0)}


@boundary
def finite_difference_forces(term, coordinates, *, displacement):
    numeric((displacement,))
    require(displacement > 0, "Positive FD displacement required")
    term.energy(coordinates)
    x = {sid: np.array(coordinates[sid], dtype=float, copy=True) for sid in term.sites}
    result = {sid: np.zeros(3) for sid in term.sites}
    for sid in term.sites:
        for axis in range(3):
            original = x[sid][axis]
            x[sid][axis] = original + displacement
            plus = term.energy(x)
            x[sid][axis] = original - displacement
            minus = term.energy(x)
            x[sid][axis] = original
            result[sid][axis] = -(plus - minus) / (2 * displacement)
    return result


@boundary
def pair_energy(
    distance,
    *,
    rmin_i,
    epsilon_i,
    rmin_j,
    epsilon_j,
    charge_i,
    charge_j,
    lj_weight,
    coulomb_weight,
):
    """Isolated all-pairs 9-6/Coulomb kernel with explicit weights; no cutoff."""
    values = (
        distance,
        rmin_i,
        epsilon_i,
        rmin_j,
        epsilon_j,
        charge_i,
        charge_j,
        lj_weight,
        coulomb_weight,
    )
    numeric(values)
    require(
        distance > 0
        and rmin_i > 0
        and rmin_j > 0
        and epsilon_i >= 0
        and epsilon_j >= 0,
        "Invalid pair data",
    )
    require(0 <= lj_weight <= 1 and 0 <= coulomb_weight <= 1, "Invalid pair weights")
    r6 = (rmin_i**6 + rmin_j**6) / 2
    rmin = r6 ** (1 / 6)
    eps = sqrt(epsilon_i * epsilon_j) * rmin_i**3 * rmin_j**3 / r6
    lj = lj_weight * eps * (2 * (rmin / distance) ** 9 - 3 * (rmin / distance) ** 6)
    coul = coulomb_weight * (332.06371 * 4.184) * charge_i * charge_j / distance
    require(isfinite(lj) and isfinite(coul), "Nonfinite pair energy")
    return {"lj_9_6": lj, "coulomb": coul, "total": lj + coul}
