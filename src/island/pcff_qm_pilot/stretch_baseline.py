"""Separate J61 baseline-adjustable Class-II Si-H research Hamiltonian.

J60/source records are immutable inputs. Equilibrium changes enter *all* bond
deviations, including fixed and fitted couplings. Analytic Cartesian/parameter
and mixed force/parameter derivatives use the existing Cartesian Jet kernels.
"""

from collections import defaultdict

import numpy as np

from island.pcff_completion import require

from .hamiltonian import (
    SOURCE_MODEL_IDENTITY,
    Jet,
    bend,
    canonical,
    cos_phi,
    geometry,
    harmonic,
    length,
    parts,
)
from .residual import FIXED, fixed_rows, validate_model

IDENTITY = "island_research_class2_sih_disiloxane_j61_stretch_v1"
BASELINE_NAMES = ("r0_SiO_A", "r0_SiH_A", "K2_SiO_kcal_mol_A2", "K2_SiH_kcal_mol_A2")


def zero():
    return Jet(0, np.zeros(27))


def bond_kind(model, sites):
    types = sorted(model["typing"][str(i)] for i in sites)
    require(types in (["osi", "sio"], ["hsi", "sio"]), "unsupported_stretch_type")
    return 0 if types[0] == "osi" else 1


def source_baseline(model):
    validate_model(model, SOURCE_MODEL_IDENTITY)
    p = np.zeros(4)
    for row in model["interaction_inventory"]:
        if row["family"] == "quartic_bond":
            k = bond_kind(model, row["sites"])
            v = row["assignment"]["values"]
            p[k], p[k + 2] = v[0], v[1] / 4.184
    return p


class Configuration:
    """Cache geometric invariants; evaluate an explicitly supplied full24 vector."""

    def __init__(self, model, xyz):
        validate_model(model, SOURCE_MODEL_IDENTITY)
        self.model = model
        self.xyz = np.asarray(xyz, float)
        self.x = geometry(xyz)
        self.parent = parts(model, xyz)
        original = fixed_rows(model["interaction_inventory"], xyz)
        self.constant = Jet(
            self.parent["fixed_energy"] - original["energy"],
            (-self.parent["fixed_forces"] + np.asarray(original["forces"])).ravel(),
        )
        self.columns = self.parent["columns"]
        self.names = [*BASELINE_NAMES, *self.columns]
        self.lookup = {r["id"]: r for r in model["interaction_inventory"]}
        self.bond_geometry = {}
        self.angle_geometry = {}
        self.cosines = {}
        for r in model["interaction_inventory"]:
            f, sites = r["family"], r["sites"]
            if f == "quartic_bond":
                self.bond_geometry[r["id"]] = length(self.x, sites)
            elif f == "quartic_angle":
                self.angle_geometry[r["id"]] = bend(self.x, sites)
            elif f in (
                "torsion_3",
                "angle-angle-torsion_1",
                "angle-torsion_3",
                "end_bond-torsion_3",
                "middle_bond-torsion_3",
            ):
                self.cosines[tuple(sites)] = cos_phi(self.x, sites)
                self.cosines[tuple(sites[::-1])] = self.cosines[tuple(sites)]

    def evaluate(self, parameters):
        validate_model(self.model, SOURCE_MODEL_IDENTITY)
        p = np.asarray(parameters, float)
        require(
            p.shape == (24,) and np.isfinite(p).all(),
            "explicit_all24_parameters_required",
        )

        def key(f, sites):
            return f + ":" + ",".join(map(str, min(tuple(sites), tuple(sites[::-1]))))

        def dr(sites):
            k = bond_kind(self.model, sites)
            d = self.bond_geometry[key("quartic_bond", sites)] - p[k]
            dp = [zero() for _ in range(4)]
            dp[k] = Jet(-1, np.zeros(27))
            return d, dp

        def dt(sites):
            name = key("quartic_angle", sites)
            return (
                self.angle_geometry[name] - self.lookup[name]["assignment"]["values"][0]
            )

        energy = self.constant
        derivatives = [zero() for _ in range(4)]
        features = defaultdict(zero)
        feature_dp = defaultdict(lambda: [zero() for _ in range(4)])
        for row in self.model["interaction_inventory"]:
            f, sites = row["family"], row["sites"]
            if f in FIXED:
                values = row["assignment"]["values"]
                if f == "quartic_bond":
                    q, dq = dr(sites)
                    kind = bond_kind(self.model, sites)
                    k2, k3, k4 = p[kind + 2], values[2] / 4.184, values[3] / 4.184
                    energy += k2 * q**2 + k3 * q**3 + k4 * q**4
                    slope = 2 * k2 * q + 3 * k3 * q**2 + 4 * k4 * q**3
                    for j in range(4):
                        derivatives[j] += dq[j] * slope
                    derivatives[kind + 2] += q**2
                elif f == "quartic_angle":
                    q = dt(sites)
                    energy += sum(values[j] / 4.184 * q ** (j + 1) for j in (1, 2, 3))
                else:
                    i, j, k = sites
                    left, dl = dr([i, j])
                    right, dright = dr([j, k])
                    v = np.asarray(values) / 4.184
                    if f == "bond-bond":
                        energy += v[0] * left * right
                        for n in range(4):
                            derivatives[n] += v[0] * (dl[n] * right + left * dright[n])
                    else:
                        t = dt(sites)
                        energy += (v[0] * left + v[1] * right) * t
                        for n in range(4):
                            derivatives[n] += (v[0] * dl[n] + v[1] * dright[n]) * t
                continue
            if row["assignment"] is not None:
                continue
            s, labels = canonical(row, self.model["typing"])
            i, j, k, l = s
            name = f + ":" + ",".join(labels)
            if f == "angle-angle":
                features[name] += dt([i, j, k]) * dt([k, j, l])
            elif f == "angle-angle-torsion_1":
                features[name] += dt([i, j, k]) * dt([j, k, l]) * self.cosines[tuple(s)]
            elif f == "bond-bond_1_3":
                a, da = dr([i, j])
                b, db = dr([k, l])
                features[name] += a * b
                for n in range(4):
                    feature_dp[name][n] += da[n] * b + a * db[n]
            else:
                ds = (
                    [dr([i, j]), dr([k, l])]
                    if f == "end_bond-torsion_3"
                    else [dr([j, k])]
                    if f == "middle_bond-torsion_3"
                    else [
                        (dt([i, j, k]), [zero() for _ in range(4)]),
                        (dt([j, k, l]), [zero() for _ in range(4)]),
                    ]
                )
                for block, (d, dd) in enumerate(ds):
                    for n in (1, 2, 3):
                        h = harmonic(self.cosines[tuple(s)], n)
                        column = name + f":block{block}:n{n}"
                        features[column] += d * h
                        for v in range(4):
                            feature_dp[column][v] += dd[v] * h
        require(set(features) == set(self.columns), "partial_coupling_accounting")
        for name, value in zip(self.columns, p[4:], strict=True):
            energy += value * features[name]
            for j in range(4):
                derivatives[j] += value * feature_dp[name][j]
        derivatives.extend(features[k] for k in self.columns)
        return {
            "energy": energy.v,
            "forces": -energy.d.reshape(9, 3),
            "energy_parameter_jacobian": np.array([v.v for v in derivatives]),
            "force_parameter_jacobian": np.stack(
                [-v.d.reshape(9, 3) for v in derivatives], axis=-1
            ),
            "model": IDENTITY,
        }


def stretching_projection(xyz, bonds):
    """Unweighted Cartesian metric; rigid SVD projection then positive-diagonal QR."""
    x = np.asarray(xyz, float)
    require(x.shape == (9, 3) and np.isfinite(x).all(), "projection_geometry")
    centered = x - x.mean(axis=0)
    rigid = np.column_stack(
        [np.tile(e, (9, 1)).ravel() for e in np.eye(3)]
        + [np.cross(e, centered).ravel() for e in np.eye(3)]
    )
    u, sv, _ = np.linalg.svd(rigid, full_matrices=False)
    require(sum(sv > sv[0] * 1e-12) == 6, "rigid_mode_rank")
    projector = np.eye(27) - u @ u.T
    raw = []
    for i, j in bonds:
        vector = x[i - 1] - x[j - 1]
        length_ = np.linalg.norm(vector)
        require(length_ > 1e-8, "projection_singular_bond")
        d = np.zeros_like(x)
        d[i - 1], d[j - 1] = vector / length_, -vector / length_
        raw.append(d.ravel())
    b = projector @ np.column_stack(raw)
    sv_b = np.linalg.svd(b, compute_uv=False)
    require(sum(sv_b > sv_b[0] * 1e-12) == 8, "stretch_subspace_rank")
    q, r = np.linalg.qr(b, mode="reduced")
    q *= np.where(np.diag(r) < 0, -1, 1)
    return {"Q": q, "rigid": u, "projector": projector, "bond_singular_values": sv_b}


def directional_curvature(plus, minus, q, step):
    require(step > 0 and np.isfinite(step), "invalid_curvature_step")
    plus, minus = np.asarray(plus), np.asarray(minus)
    require(plus.shape == minus.shape == (8, 9, 3), "directional_response_mapping")
    response = -(plus - minus).reshape(8, 27).T / (2 * step)
    reduced = q.T @ response
    return response, reduced, (reduced + reduced.T) / 2
