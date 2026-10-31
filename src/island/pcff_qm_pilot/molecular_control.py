"""Graph-bound, dimension-independent coordinates for bounded research controls.

Historical nine-atom kernels and identities are unchanged.
"""

import numpy as np

from island.pcff_completion import require

from .hamiltonian import Jet, bend, cos_phi, dot, norm, subtract


def cartesian(xyz):
    x = np.asarray(xyz, float)
    require(
        x.ndim == 2 and x.shape[1] == 3 and len(x) >= 3 and np.isfinite(x).all(),
        "invalid_molecular_geometry",
    )
    eye = np.eye(x.size)
    return [[Jet(x[i, j], eye[3 * i + j]) for j in range(3)] for i in range(len(x))]


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def torsion(x, sites):
    i, j, k, l = sites
    axis = subtract(x[k - 1], x[j - 1])
    axis = [v / norm(axis) for v in axis]
    a, b = subtract(x[i - 1], x[j - 1]), subtract(x[l - 1], x[k - 1])
    p = subtract(a, [dot(a, axis) * v for v in axis])
    q = subtract(b, [dot(b, axis) * v for v in axis])
    c = cos_phi(x, sites)
    s = dot(axis, cross(p, q)) / norm(p) / norm(q)
    return Jet(np.arctan2(s.v, c.v), (c.v * s.d - s.v * c.d) / (c.v * c.v + s.v * s.v))


def rigid_projection(xyz):
    x = np.asarray(xyz, float)
    x = x - x.mean(axis=0)
    a = np.column_stack(
        [np.tile(e, (len(x), 1)).ravel() for e in np.eye(3)]
        + [np.cross(e, x).ravel() for e in np.eye(3)]
    )
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    require(sum(s > s[0] * 1e-12) == 6, "rigid_rank_deficiency")
    return u, np.eye(x.size) - u @ u.T


def nonstretch_projection(xyz, stretch, seeds):
    x = cartesian(xyz)
    u, p = rigid_projection(xyz)
    stretch = np.asarray(stretch, float)
    require(
        stretch.shape[0] == len(x) * 3
        and np.allclose(stretch.T @ stretch, np.eye(stretch.shape[1]), atol=1e-10),
        "invalid_stretch_mapping",
    )
    raw = np.column_stack(
        [
            (bend(x, s["sites"]) if s["kind"] == "angle" else torsion(x, s["sites"])).d
            for s in seeds
        ]
    )
    projected = (p - stretch @ stretch.T) @ raw
    sv = np.linalg.svd(projected, compute_uv=False)
    require(sum(sv > sv[0] * 1e-10) == len(seeds), "nonstretch_rank_deficiency")
    q, r = np.linalg.qr(projected, mode="reduced")
    q *= np.where(np.diag(r) < 0, -1, 1)
    return {
        "Q": q,
        "raw_coordinate_gradients": raw,
        "projected_coordinate_gradients": projected,
        "singular_values": sv,
        "rigid_basis": u,
    }


def response(plus, minus, step):
    a, b = np.asarray(plus, float), np.asarray(minus, float)
    require(
        a.shape == b.shape
        and a.ndim == 3
        and a.shape[2] == 3
        and step > 0
        and np.isfinite(a).all()
        and np.isfinite(b).all(),
        "response_mapping_or_units",
    )
    return -(a - b).reshape(len(a), -1).T / (2 * step)
