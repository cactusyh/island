"""Geometry-only coupled design screening; no physical coefficient inference."""

from itertools import combinations

import numpy as np

from island.pcff_completion import digest, require, signed

from .hmds_hamiltonian import parts
from .hmds_paths import fragment, path


def matrix(model, geometries, anchor=None):
    anchor = np.asarray(geometries[0] if anchor is None else anchor)
    a = parts(model, anchor)
    rows = []
    for x in geometries:
        p = parts(model, x)
        require(p["columns"] == a["columns"], "coupled_design_mapping")
        rows.extend(
            [
                (p["energy_jacobian"] - a["energy_jacobian"]) / 0.05,
                *p["force_jacobian"].reshape(81, -1) / 0.5,
            ]
        )
    return np.array(rows), a["columns"]


def metrics(a):
    a = np.array(a)
    norms = np.linalg.norm(a, axis=0)
    scaled = a / np.maximum(norms, 1e-30)
    _, sv, vt = np.linalg.svd(scaled, full_matrices=False)
    rank = int(np.sum(sv > sv[0] * 1e-10)) if sv[0] > 0 else 0
    corr = scaled.T @ scaled
    np.fill_diagonal(corr, 0)
    i, j = np.unravel_index(np.argmax(abs(corr)), corr.shape)
    return {
        "columns": a.shape[1],
        "rows": a.shape[0],
        "rank": rank,
        "column_norms": norms.tolist(),
        "scaled_condition": float(sv[0] / max(sv[-1], 1e-300)),
        "singular_values": sv.tolist(),
        "weakest_scaled_combination": vt[-1].tolist(),
        "strongest_column_correlation": {
            "indices": [int(i), int(j)],
            "cosine": float(corr[i, j]),
        },
        "identifiable_for_screen": bool(
            rank == a.shape[1]
            and sv[0] / max(sv[-1], 1e-300) < 1000
            and np.all(norms > 1e-12)
        ),
    }


def angular(model, xyz, sites, delta):
    i, j, k = sites
    x = np.asarray(xyz, float)
    moving = fragment(model, [j, i], i)
    indices = np.array(moving) - 1
    u = x[i - 1] - x[j - 1]
    v = x[k - 1] - x[j - 1]
    axis = np.cross(u, v)
    axis /= np.linalg.norm(axis)
    r = x[indices] - x[j - 1]
    c = np.cos(delta)
    s = np.sin(delta)
    y = x.copy()
    y[indices] = (
        x[j - 1] + c * r + s * np.cross(axis, r) + (1 - c) * np.outer(r @ axis, axis)
    )
    return y


def pool(model, center, count=48):
    rng = np.random.default_rng(14063)
    rows = []
    stretches = [[2, 1], [1, 10], [5, 6]]
    angles = [[1, 2, 3], [1, 2, 5], [2, 5, 6], [2, 1, 10], [10, 1, 11]]
    for i in range(count):
        x = np.array(center)
        actions = []
        for bond in stretches:
            d = float(rng.uniform(-0.035, 0.035))
            spec = {"kind": "stretch", "bond": bond, "moving_root": bond[1]}
            x = path(model, x, spec, d)["coordinates"]
            actions.append(
                {"kind": "stretch", "sites": bond, "value": d, "units": "angstrom"}
            )
        for sites in angles:
            d = float(rng.uniform(-0.06, 0.06))
            x = angular(model, x, sites, d)
            actions.append(
                {"kind": "bend", "sites": sites, "value": d, "units": "radian"}
            )
        for bond in ([2, 1], [5, 6]):
            d = float(rng.uniform(-0.5, 0.5))
            spec = {"kind": "rotation", "bond": bond, "moving_root": bond[1]}
            x = path(model, x, spec, d)["coordinates"]
            actions.append(
                {"kind": "rotation", "sites": bond, "value": d, "units": "radian"}
            )
        rows.append(
            {
                "id": f"joint-{i:03d}",
                "coordinates": x.tolist(),
                "geometry_identity": digest(x.tolist()),
                "actions": actions,
                "QM_data_supplied": False,
            }
        )
    return rows


def screen(model, center, current, maximum_additions=6):
    a, names = matrix(model, current, anchor=center)
    before = metrics(a)
    candidates = pool(model, center)
    blocks = [matrix(model, [r["coordinates"]], anchor=center)[0] for r in candidates]
    chosen = []
    remaining = list(range(len(blocks)))
    combined = a.copy()
    history = []
    # Rank/conditioning only, no reference energies, forces or validation values.
    for _ in range(maximum_additions):
        if metrics(combined)["identifiable_for_screen"]:
            break
        options = []
        for i in remaining:
            m = metrics(np.vstack([combined, blocks[i]]))
            options.append(
                ((m["rank"], -np.log10(max(m["scaled_condition"], 1))), i, m)
            )
        _, i, m = max(options, key=lambda v: v[0])
        chosen.append(i)
        remaining.remove(i)
        combined = np.vstack([combined, blocks[i]])
        history.append({"id": candidates[i]["id"], "metrics": m})
    # Exhaust all subsets of the bounded selected pool to minimize cardinality
    # within this pool. This is not a global experimental-design minimality proof.
    smallest = []
    for n in range(1, len(chosen) + 1):
        passing = []
        for subset in combinations(chosen, n):
            m = metrics(np.vstack([a, *[blocks[i] for i in subset]]))
            if m["identifiable_for_screen"]:
                passing.append((m["scaled_condition"], subset, m))
        if passing:
            _, smallest, after = min(passing, key=lambda v: v[0])
            break
    else:
        after = metrics(combined)
    if before["identifiable_for_screen"]:
        smallest = []
        after = before
    return signed(
        {
            "names": names,
            "current": before,
            "current_geometry_count": len(current),
            "candidate_pool": candidates,
            "selected_additional_configurations": [candidates[i] for i in smallest],
            "after_selected": after,
            "greedy_history": history,
            "best_attempted_configurations": [candidates[i] for i in chosen],
            "passing_additional_design_found": bool(smallest)
            or before["identifiable_for_screen"],
            "smallest_design_scope": "Minimum-cardinality passing subset of at most6 greedily selected graph-controlled candidates from48, not a global minimum proof.",
            "reference_values_used_for_selection": False,
            "physical_parameters_fitted": False,
            "design_weights": "diagnostic E/0.05 kcal/mol, F/0.5 kcal/mol/A; not an adopted fitting objective",
            "unknown_nuisance_terms_omitted": False,
            "known_unknown_count": len(names),
        }
    )
