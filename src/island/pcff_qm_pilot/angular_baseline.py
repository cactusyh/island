"""Small preregistered J62 research changes; historical source rows immutable.

Exact mixed Cartesian/parameter chain rules propagate local angle-equilibrium
shifts through every angle-dependent fixed and fitted Class-II term.
"""

import numpy as np

from island.pcff_completion import require

from .hamiltonian import canonical, harmonic
from .stretch_baseline import Configuration as StretchConfiguration
from .stretch_baseline import bond_kind

EXTRA_NAMES = (
    "main_n3_delta_kcal_mol",
    "HSiH_theta0_delta_rad",
    "HSiO_theta0_delta_rad",
    "SiOSi_K2_delta_kcal_mol_rad2",
)
IDENTITY = "island_research_class2_sih_disiloxane_j62_nested_v1"


class Mixed:
    """Value, Cartesian gradient, parameter gradient and mixed derivative."""

    def __init__(self, v, x=None, p=None, xp=None):
        self.v = float(v)
        self.x = np.zeros(27) if x is None else np.asarray(x)
        self.p = np.zeros(28) if p is None else np.asarray(p)
        self.xp = np.zeros((27, 28)) if xp is None else np.asarray(xp)

    @classmethod
    def geometric(cls, j):
        return cls(j.v, j.d)

    @classmethod
    def variable(cls, v, index):
        p = np.zeros(28)
        p[index] = 1
        return cls(v, p=p)

    def __add__(self, other):
        o = other if isinstance(other, Mixed) else Mixed(other)
        return Mixed(self.v + o.v, self.x + o.x, self.p + o.p, self.xp + o.xp)

    __radd__ = __add__

    def __neg__(self):
        return Mixed(-self.v, -self.x, -self.p, -self.xp)

    def __sub__(self, o):
        return self + -o

    def __rsub__(self, o):
        return -self + o

    def __mul__(self, other):
        o = other if isinstance(other, Mixed) else Mixed(other)
        return Mixed(
            self.v * o.v,
            self.x * o.v + o.x * self.v,
            self.p * o.v + o.p * self.v,
            self.xp * o.v
            + o.xp * self.v
            + np.outer(self.x, o.p)
            + np.outer(o.x, self.p),
        )

    __rmul__ = __mul__

    def __pow__(self, n):
        if n == 0:
            return Mixed(1)
        a = self.v**n
        b = n * self.v ** (n - 1)
        c = n * (n - 1) * self.v ** (n - 2) if n >= 2 else 0
        return Mixed(
            a, b * self.x, b * self.p, b * self.xp + c * np.outer(self.x, self.p)
        )


class Configuration:
    def __init__(self, model, xyz):
        self.base = StretchConfiguration(model, xyz)
        self.model = model
        self.xyz = np.array(xyz)
        self.names = [*self.base.names, *EXTRA_NAMES]

    def evaluate(self, parameters):
        p = np.asarray(parameters, float)
        require(
            p.shape == (28,) and np.isfinite(p).all(),
            "explicit_all28_parameters_required",
        )
        b = self.base.evaluate(p[:24])
        dp = np.r_[b["energy_parameter_jacobian"], np.zeros(4)]
        xp = np.pad(-b["force_parameter_jacobian"].reshape(27, 24), ((0, 0), (0, 4)))
        energy = Mixed(b["energy"], -b["forces"].ravel(), dp, xp)
        variables = [Mixed.variable(v, i) for i, v in enumerate(p)]

        def key(f, s):
            return f + ":" + ",".join(map(str, min(tuple(s), tuple(s[::-1]))))

        def dr(s):
            return (
                Mixed.geometric(self.base.bond_geometry[key("quartic_bond", s)])
                - variables[bond_kind(self.model, s)]
            )

        def dt(s, shifted):
            row = self.base.lookup[key("quartic_angle", s)]
            t = (
                Mixed.geometric(self.base.angle_geometry[row["id"]])
                - row["assignment"]["values"][0]
            )
            labels = tuple(self.model["typing"][str(i)] for i in s)
            labels = min(labels, labels[::-1])
            if shifted:
                if labels == ("hsi", "sio", "hsi"):
                    t -= variables[25]
                elif labels == ("hsi", "sio", "osi"):
                    t -= variables[26]
            return t

        def coefficient(name):
            return variables[self.names.index(name)]

        for row in self.model["interaction_inventory"]:
            f, s, a = row["family"], row["sites"], row["assignment"]
            if f == "torsion_3":
                energy += variables[24] * Mixed.geometric(
                    1 + harmonic(self.base.cosines[tuple(s)], 3)
                )
                continue
            if f == "quartic_angle":
                old, new = dt(s, False), dt(s, True)
                v = a["values"]
                energy += sum(
                    v[n] / 4.184 * (new ** (n + 1) - old ** (n + 1)) for n in (1, 2, 3)
                )
                if tuple(self.model["typing"][str(i)] for i in s) == (
                    "sio",
                    "osi",
                    "sio",
                ):
                    energy += variables[27] * new**2
            elif f == "bond-angle":
                i, j, k = s
                v = a["values"]
                energy += (v[0] / 4.184 * dr([i, j]) + v[1] / 4.184 * dr([j, k])) * (
                    dt(s, True) - dt(s, False)
                )
            elif (
                f in ("angle-angle", "angle-torsion_3", "angle-angle-torsion_1")
                and a is None
            ):
                s, labels = canonical(row, self.model["typing"])
                i, j, k, l = s
                name = f + ":" + ",".join(labels)
                if f == "angle-angle":
                    energy += coefficient(name) * (
                        dt([i, j, k], True) * dt([k, j, l], True)
                        - dt([i, j, k], False) * dt([k, j, l], False)
                    )
                elif f == "angle-angle-torsion_1":
                    energy += (
                        coefficient(name)
                        * (
                            dt([i, j, k], True) * dt([j, k, l], True)
                            - dt([i, j, k], False) * dt([j, k, l], False)
                        )
                        * Mixed.geometric(self.base.cosines[tuple(s)])
                    )
                else:
                    for block, sites in enumerate(([i, j, k], [j, k, l])):
                        for n in (1, 2, 3):
                            energy += (
                                coefficient(name + f":block{block}:n{n}")
                                * (dt(sites, True) - dt(sites, False))
                                * Mixed.geometric(
                                    harmonic(self.base.cosines[tuple(s)], n)
                                )
                            )
        return {
            "energy": energy.v,
            "forces": -energy.x.reshape(9, 3),
            "energy_parameter_jacobian": energy.p,
            "force_parameter_jacobian": -energy.xp.reshape(9, 3, 28),
            "model": IDENTITY,
        }


def source_extensions(model):
    amps = [
        r["assignment"]["values"][4] / 4.184
        for r in model["interaction_inventory"]
        if r["family"] == "torsion_3"
    ]
    require(np.ptp(amps) < 1e-12, "main_torsion_type_symmetry_unbound")
    angles = [
        r["assignment"]["values"][1] / 4.184
        for r in model["interaction_inventory"]
        if r["family"] == "quartic_angle"
        and tuple(model["typing"][str(i)] for i in r["sites"]) == ("sio", "osi", "sio")
    ]
    require(len(angles) == 1, "central_angle_mapping")
    return amps[0], angles[0]
