"""Strict unit conversions and statistics for independently executed single points."""

import hashlib
import json
import math
from pathlib import Path

import numpy as np

KCAL_TO_KJ = 4.184
ANGSTROM_PER_NM = 10.0
ATM_TO_BAR = 1.01325
KJ_MOL_NM3_TO_BAR = 16.605390671738468
COULOMB_OPENMM = 138.93545764438198


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def signed(value):
    body = {k: v for k, v in value.items() if k != "identity"}
    digest = hashlib.sha256(
        json.dumps(
            body, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()
    return dict(body, identity=digest)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parse_dump(path):
    lines = Path(path).read_text().splitlines()
    require(lines[:1] == ["ITEM: TIMESTEP"], "Missing dump header")
    require(lines.count("ITEM: TIMESTEP") == 1, "Expected exactly one retained frame")
    n = int(lines[3])
    header = next(i for i, v in enumerate(lines) if v.startswith("ITEM: ATOMS "))
    names = lines[header].split()[2:]
    require(
        names == ["id", "type", "q", "x", "y", "z", "fx", "fy", "fz"],
        "Changed dump mapping",
    )
    rows = np.asarray([list(map(float, v.split())) for v in lines[header + 1 :]])
    require(
        rows.shape == (n, 9) and np.isfinite(rows).all(),
        "Truncated/nonfinite force dump",
    )
    ids = rows[:, 0].astype(int)
    require(
        np.array_equal(ids, np.arange(1, n + 1)),
        "Duplicate, missing or reordered atom IDs",
    )
    box = np.asarray([list(map(float, v.split())) for v in lines[5:8]])
    return {
        "atom_ids": ids.tolist(),
        "charges_e": rows[:, 2].tolist(),
        "positions_nm": (rows[:, 3:6] / ANGSTROM_PER_NM).tolist(),
        "box_lengths_nm": ((box[:, 1] - box[:, 0]) / ANGSTROM_PER_NM).tolist(),
        "forces_kj_mol_nm": (rows[:, 6:9] * KCAL_TO_KJ * ANGSTROM_PER_NM).tolist(),
    }


def force_metrics(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    require(
        a.shape == b.shape and a.ndim == 2 and a.shape[1] == 3 and len(a),
        "Wrong force-vector mapping",
    )
    require(np.isfinite(a).all() and np.isfinite(b).all(), "Nonfinite force vector")
    difference = a - b
    norms = np.linalg.norm(difference, axis=1)
    magnitude = np.linalg.norm(a, axis=1)
    scale = float(np.sqrt(np.mean(a * a)))
    return {
        "units": "kJ/(mol nm)",
        "per_atom_difference": difference.tolist(),
        "rms_component_error": float(np.sqrt(np.mean(difference * difference))),
        "rms_atom_norm_error": float(np.sqrt(np.mean(norms * norms))),
        "maximum_component_error": float(np.max(np.abs(difference))),
        "maximum_atom_norm_error": float(max(norms)),
        "atom_norm_error_percentiles": {
            str(p): float(np.percentile(norms, p)) for p in (0, 50, 90, 95, 99, 100)
        },
        "maximum_reference_atom_force": float(max(magnitude)),
        "rms_reference_force": scale,
        "relative_rms_component_error": float(
            np.sqrt(np.mean(difference * difference)) / scale
        )
        if scale
        else None,
        "total_force_difference": np.sum(difference, axis=0).tolist(),
    }


def compare(a, b):
    require(a["frame_identity"] == b["frame_identity"], "Different physical frames")
    metrics = force_metrics(a["forces_kj_mol_nm"], b["forces_kj_mol_nm"])
    energy = a["energy_kj_mol"] - b["energy_kj_mol"]
    components = {
        k: a["components_kj_mol"][k] - b["components_kj_mol"][k]
        for k in ("bond", "angle", "proper", "improper", "lj", "coulomb")
    }
    f1, f2 = (
        np.asarray(a["forces_kj_mol_nm"]) / 10,
        np.asarray(b["forces_kj_mol_nm"]) / 10,
    )
    # Preserve J22/J40 bounds, including their original force units and rtol.
    passed = (
        abs(energy) <= 1e-5
        and max(map(abs, components.values())) <= 1e-5
        and np.allclose(f1, f2, atol=1e-5, rtol=2e-10)
    )
    return {
        "signed_energy_error_kj_mol": energy,
        "absolute_energy_error_kj_mol": abs(energy),
        "relative_energy_error": abs(energy / a["energy_kj_mol"])
        if a["energy_kj_mol"]
        else None,
        "relative_coulomb_energy_error": abs(
            components["coulomb"] / a["components_kj_mol"]["coulomb"]
        )
        if a["components_kj_mol"]["coulomb"]
        else None,
        "signed_component_errors_kj_mol": components,
        "forces": metrics,
        "energy_bound_kj_mol": 1e-5,
        "force_bound_kj_mol_angstrom": 1e-5,
        "force_relative_bound": 2e-10,
        "passed": bool(passed),
    }


def analytic_ewald_terms(charges, alpha_nm, volume_nm3, coulomb_constant):
    require(alpha_nm > 0 and volume_nm3 > 0, "Invalid Ewald geometry")
    return {
        "self_kj_mol": -coulomb_constant
        * alpha_nm
        * math.fsum(q * q for q in charges)
        / math.sqrt(math.pi),
        "background_kj_mol": -coulomb_constant
        * math.pi
        * math.fsum(charges) ** 2
        / (2 * alpha_nm**2 * volume_nm3),
        "conducting_surface_kj_mol": 0.0,
        "surface_convention": "3D periodic conducting (tin-foil), no slab correction",
    }
