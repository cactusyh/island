"""Tree-fragment internal-coordinate scans, with independent measurement checks."""

from itertools import pairwise, product

import numpy as np

from island.pcff_completion import _finite, _verified, digest, require


def graph_record(system, p):
    return {
        "graph_identity": p["graph_identity"],
        "atoms": [
            {
                "id": i,
                "element": system.topology.sites[i].element,
                "formal_charge": system.topology.sites[i].formal_charge,
            }
            for i in p["atom_order"]
        ],
        "bonds": sorted(
            [
                min(b.site1, b.site2),
                max(b.site1, b.site2),
                float(b.order),
                bool(b.aromatic),
            ]
            for b in system.topology.bonds.values()
        ),
    }


def measure(xyz, order, coordinate):
    v = np.asarray([xyz[order.index(i)] for i in coordinate["sites"]])
    kind = coordinate["kind"]
    if kind == "bond":
        return float(np.linalg.norm(v[1] - v[0]))
    if kind == "angle":
        a, b = v[0] - v[1], v[2] - v[1]
        return float(
            np.arccos(
                np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)
            )
        )
    require(kind == "torsion", "unsupported_scan_coordinate")
    axis = v[2] - v[1]
    axis /= np.linalg.norm(axis)
    a, b = v[0] - v[1], v[3] - v[2]
    a -= np.dot(a, axis) * axis
    b -= np.dot(b, axis) * axis
    require(np.linalg.norm(a) * np.linalg.norm(b) > 1e-12, "degenerate_scan_axis")
    return float(np.arctan2(np.dot(np.cross(axis, a), b), np.dot(a, b)))


def _branch(system, order, left, right):
    adj = {i: set() for i in order}
    for b in system.topology.bonds.values():
        adj[b.site1].add(b.site2)
        adj[b.site2].add(b.site1)
    require(right in adj[left], "scan_path_not_in_graph")
    adj[left].remove(right)
    adj[right].remove(left)
    group, todo = set(), [right]
    while todo:
        i = todo.pop()
        if i not in group:
            group.add(i)
            todo.extend(adj[i] - group)
    require(left not in group, "cyclic_scan_domain_forbidden")
    return [order.index(i) for i in sorted(group)]


def _adjust(system, xyz, order, coord, delta):
    sites, kind = coord["sites"], coord["kind"]
    adjacency = {frozenset((b.site1, b.site2)) for b in system.topology.bonds.values()}
    require(
        all(frozenset(pair) in adjacency for pair in pairwise(sites)),
        "scan_path_not_in_graph",
    )
    if kind == "bond":
        left, right = sites
    else:
        left, right = sites[1:3]
    group = _branch(system, order, left, right)
    origin = xyz[order.index(left)].copy()
    axis = xyz[order.index(right)] - origin
    if kind == "angle":
        axis = np.cross(xyz[order.index(sites[0])] - origin, axis)
    require(np.linalg.norm(axis) > 1e-12, "degenerate_scan_axis")
    axis /= np.linalg.norm(axis)
    if kind == "bond":
        xyz[group] += delta * axis
    else:
        if kind == "torsion":
            origin = xyz[order.index(right)].copy()
        v = xyz[group] - origin
        xyz[group] = origin + (
            v * np.cos(delta)
            + np.cross(axis, v) * np.sin(delta)
            + np.outer(v @ axis, axis) * (1 - np.cos(delta))
        )


def geometries(system, assessment, scan):
    """Cartesian product of explicit displacements; no optimization or QM data.

    Sequential rigid subtrees are corrected iteratively when coordinates couple.
    The recorded coordinate residual is independently remeasured, never assumed.
    """
    p = _verified(system, assessment)
    order = p["atom_order"]
    xyz = np.array([system.coordinates.get(i).tolist() for i in order])
    coords = scan["coordinates"]
    require(0 < len(coords) <= 3, "scan_dimension_not_supported")
    baseline, grids = [], []
    for c in coords:
        require(
            c["units"] == ("angstrom" if c["kind"] == "bond" else "radian"),
            "scan_units_mutation",
        )
        expected = {"bond": 2, "angle": 3, "torsion": 4}.get(c["kind"])
        require(
            len(c["sites"]) == expected and len(set(c["sites"])) == expected,
            "scan_sites_required",
        )
        values = _finite(c["displacements"])
        require(
            values.ndim == 1
            and values.size >= 3
            and len(set(values.tolist())) == len(values),
            "scan_grid_incomplete",
        )
        baseline.append(measure(xyz, order, c))
        grids.append(values.tolist())
    samples, seen = [], set()
    for k, offsets in enumerate(product(*grids)):
        x = xyz.copy()
        goals = np.asarray(baseline) + offsets
        for c, goal in zip(coords, goals, strict=True):
            require(c["kind"] != "bond" or goal > 0, "nonpositive_deformed_bond")
            require(c["kind"] != "angle" or 0 < goal < np.pi, "invalid_deformed_angle")
        for _ in range(40):
            for c, goal in zip(coords, goals, strict=True):
                delta = goal - measure(x, order, c)
                if c["kind"] == "torsion":
                    delta = (delta + np.pi) % (2 * np.pi) - np.pi
                _adjust(system, x, order, c, delta)
            actual = np.array([measure(x, order, c) for c in coords])
            residuals = actual - goals
            for j, c in enumerate(coords):
                if c["kind"] == "torsion":
                    residuals[j] = (residuals[j] + np.pi) % (2 * np.pi) - np.pi
            if np.max(np.abs(residuals)) < 1e-10:
                break
        require(np.max(np.abs(residuals)) < 1e-10, "coupled_scan_not_reconstructed")
        require(np.isfinite(x).all(), "nonfinite_scan_geometry")
        geometry_hash = digest(x.tolist())
        require(geometry_hash not in seen, "duplicate_scan_geometry")
        seen.add(geometry_hash)
        samples.append(
            {
                "id": f"point-{k:04d}",
                "coordinates": x.tolist(),
                "atom_order": order,
                "geometry_checksum": geometry_hash,
                "offsets": list(offsets),
                "baseline_coordinates": baseline,
                "measured_coordinates": actual.tolist(),
                "coordinate_residuals": residuals.tolist(),
                "split": "held_out" if k % 5 == 1 else "train",
                "expected_output_filename": f"point-{k:04d}.out",
                "calculation_status": "input_only_not_executed",
            }
        )
    return samples


def coordinate(kind, sites, step=None):
    step = (
        step
        if step is not None
        else {"bond": 0.02, "angle": 0.04, "torsion": 0.2}[kind]
    )
    return {
        "kind": kind,
        "sites": list(sites),
        "units": "angstrom" if kind == "bond" else "radian",
        "displacements": [j * step for j in (-2, -1, 0, 1, 2)],
        "scan_range": [-2 * step, 2 * step],
        "scan_spacing": step,
        "reference": "displacement from exact J54 input geometry; not an optimized equilibrium",
    }


def target_scan(row):
    a, b, c, d = row["sites"]
    if row["family"] == "angle-angle":
        # J54 dependencies are (a,b,c) and (c,b,d). Move independent outer branches.
        coords = [coordinate("angle", [c, b, a]), coordinate("angle", [c, b, d])]
    else:
        require(row["family"] == "bond-bond_1_3", "unsupported_missing_term_family")
        coords = [coordinate("bond", [b, a]), coordinate("bond", [c, d])]
    return {
        "coordinates": coords,
        "kind": "independent_paired_deformation",
        "joint_fit": False,
        "target_coefficients": 1,
        "policy": "isolate this directed coefficient or supply an explicit full-rank source-bound residual; no combined fit inferred",
    }
