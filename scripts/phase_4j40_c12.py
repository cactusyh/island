"""Explicit special-pair completion for a single-frame LAMMPS comparison.

This is a comparison adapter, not a change to a production force field or a
claim that unmodified LAMMPS implements untruncated OpenMM exceptions.
"""

import copy
import math

import numpy as np

from scripts.phase_4j40_c11_numeric import (
    KJ_MOL_NM3_TO_BAR,
    compare,
    require,
    signed,
)

# Pinned LAMMPS real-unit Coulomb prefactor, unchanged from C11.
COULOMB = 138.935456264


def binding(backend, frame, setting):
    """Validate identities before computing any completion."""
    require(frame == signed(frame), "Changed frame identity")
    for key in (
        "final_graph_identity",
        "assignment_identity",
        "charge_identity",
        "typed_graph_identity",
        "numerical_data",
    ):
        require(frame[key] == backend[key], "Changed source/graph/charge: " + key)
    data = backend["numerical_data"]
    ids = [s["id"] for s in data["sites"]]
    require(ids == sorted(set(ids)) == frame["atom_ids"], "Changed atom ordering")
    xyz = np.asarray(frame["positions_nm"], dtype=float)
    box = np.asarray(frame["box_lengths_nm"], dtype=float)
    require(
        xyz.shape == (len(ids), 3) and np.isfinite(xyz).all(), "Invalid coordinates"
    )
    require(
        box.shape == (3,) and np.isfinite(box).all() and (box > 0).all(), "Invalid box"
    )
    require(
        all(
            type(site["charge"]) in (int, float) and math.isfinite(site["charge"])
            for site in data["sites"]
        ),
        "Invalid source charge",
    )
    expected_mixing = "arithmetic" if data["family"] == "GAFF2" else "geometric"
    require(
        setting.get("mixing", expected_mixing) == expected_mixing, "Changed LJ mixing"
    )
    lj14, coul14 = backend["config"]["one_four_scaling"]
    require(
        setting.get("lj14", lj14) == lj14 and setting.get("coul14", coul14) == coul14,
        "Changed native special scaling",
    )
    cutoff = setting["cutoff_nm"]
    require(
        type(cutoff) in (int, float)
        and math.isfinite(cutoff)
        and 0 < cutoff < min(box) / 2,
        "Unsafe cutoff",
    )
    require(backend["config"]["switching_policy"] == "none", "Unsupported switching")
    require(data["family"] in ("GAFF2", "OPLS-AA"), "Unsupported interaction model")
    # Explicit membership is independently checked against shortest graph paths.
    adjacency = {i: set() for i in ids}
    for bond in data["bonds"]:
        a, b = bond["sites"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    paths = {}
    for start in ids:
        reached, frontier = {start}, {start}
        for distance in (1, 2, 3):
            frontier = set.union(set(), *(adjacency[i] for i in frontier)) - reached
            reached |= frontier
            for end in frontier:
                if start < end:
                    paths[(start, end)] = distance
    expected = {}
    lj14, coul14 = backend["config"]["one_four_scaling"]
    for pair, distance in paths.items():
        expected[pair] = (lj14, coul14) if distance == 3 else (0, 0)
    actual = {tuple(p["sites"]): (p["lj"], p["coulomb"]) for p in data["pairs"]}
    require(
        len(actual) == len(data["pairs"]) and actual == expected,
        "Changed special-pair membership/scaling",
    )
    return data, xyz, box, paths


def completion(backend, frame, setting, native):
    """Missing pair energy, force and pair virial, in kJ/mol and nm.

    For r>rc, LAMMPS omits both the special LJ term and its direct-space
    Ewald expression: C*q_i*q_j*(erfc(alpha*r) - (1-f_c))/r.
    The erfc tail is included, not approximated as zero. The reciprocal-space
    output is untouched. The log's rounded alpha is recorded as a limitation.
    """
    data, xyz, box, paths = binding(backend, frame, setting)
    require(native["frame_identity"] == frame["identity"], "Changed native frame")
    require(native["setting"] == setting, "Changed native cutoff/tolerance")
    require(native["atom_ids"] == frame["atom_ids"], "Changed native atom ordering")
    require(
        native["charges_e"] == [s["charge"] for s in data["sites"]],
        "Changed native charges",
    )
    alpha = native["alpha_nm_from_8_digit_log"]
    require(math.isfinite(alpha) and alpha > 0, "Invalid Ewald alpha")
    index = {sid: i for i, sid in enumerate(frame["atom_ids"])}
    sites = {s["id"]: s for s in data["sites"]}
    forces = np.zeros_like(xyz)
    virial = np.zeros((3, 3))
    terms = []
    for pair in data["pairs"]:
        a, b = pair["sites"]
        delta = xyz[index[b]] - xyz[index[a]]
        delta -= box * np.rint(delta / box)
        r = np.linalg.norm(delta)
        require(
            r > 0 and not np.any(np.isclose(abs(delta), box / 2, atol=1e-12, rtol=0)),
            "Coincident pair or ambiguous periodic image",
        )
        if r <= setting["cutoff_nm"]:
            continue
        x, y = sites[a], sites[b]
        sigma = (
            (x["sigma"] + y["sigma"]) / 2
            if data["family"] == "GAFF2"
            else math.sqrt(x["sigma"] * y["sigma"])
        ) / 10
        eps = math.sqrt(x["epsilon"] * y["epsilon"]) * pair["lj"]
        sr6 = (sigma / r) ** 6
        lj = 4 * eps * (sr6**2 - sr6)
        dlj = -24 * eps * (2 * sr6**2 - sr6) / r
        q = COULOMB * x["charge"] * y["charge"]
        erfctail = math.erfc(alpha * r)
        c = erfctail - (1 - pair["coulomb"])
        coul = q * c / r
        dcoul = q * (
            -2 * alpha / math.sqrt(math.pi) * math.exp(-((alpha * r) ** 2)) / r
            - c / r**2
        )
        force_lj = dlj * delta / r
        force_coul = dcoul * delta / r
        force = force_lj + force_coul  # Force on a = dU/dr * (r_b-r_a)/r.
        forces[index[a]] += force
        forces[index[b]] -= force
        tensor = -np.outer(delta, force)
        virial += tensor
        terms.append(
            {
                "atom_ids": [a, b],
                "path_length": paths[(a, b)],
                "lj_scaling": pair["lj"],
                "coulomb_scaling": pair["coulomb"],
                "distance_nm": float(r),
                "lj_kj_mol": float(lj),
                "coulomb_kj_mol": float(coul),
                "erfc_tail": erfctail,
                "force_on_first_kj_mol_nm": force.tolist(),
                "force_lj_on_first_kj_mol_nm": force_lj.tolist(),
                "force_coulomb_on_first_kj_mol_nm": force_coul.tolist(),
                "virial_lj_kj_mol": (-np.outer(delta, force_lj)).tolist(),
                "virial_coulomb_kj_mol": (-np.outer(delta, force_coul)).tolist(),
                "virial_kj_mol": tensor.tolist(),
            }
        )
    return signed(
        {
            "schema": "island.j40-c12-special-pair-completion.v1",
            "frame_identity": frame["identity"],
            "setting": setting,
            "numerical_data_identity": signed(data)["identity"],
            "native_identity": signed(native)["identity"],
            "alpha_nm_from_8_digit_log": alpha,
            "pairs": terms,
            "energy_lj_kj_mol": math.fsum(p["lj_kj_mol"] for p in terms),
            "energy_coulomb_kj_mol": math.fsum(p["coulomb_kj_mol"] for p in terms),
            "forces_kj_mol_nm": forces.tolist(),
            "virial_kj_mol": virial.tolist(),
            "pressure_increment_bar": float(
                np.trace(virial) / (3 * np.prod(box)) * KJ_MOL_NM3_TO_BAR
            ),
            "method": "explicit analytical single-frame comparison completion; native output retained",
        }
    )


def corrected(native, record):
    require(record == signed(record), "Changed completion identity")
    require(
        record["native_identity"] == signed(native)["identity"], "Changed native output"
    )
    original_forces = np.asarray(native["forces_kj_mol_nm"])
    added_forces = np.asarray(record["forces_kj_mol_nm"])
    require(
        original_forces.shape == added_forces.shape
        and np.isfinite(original_forces).all()
        and np.isfinite(added_forces).all(),
        "Changed force-vector mapping",
    )
    result = copy.deepcopy(native)
    result["engine"] = "LAMMPS plus explicit special-pair completion (composite)"
    result["energy_kj_mol"] += (
        record["energy_lj_kj_mol"] + record["energy_coulomb_kj_mol"]
    )
    result["components_kj_mol"]["lj"] += record["energy_lj_kj_mol"]
    result["components_kj_mol"]["coulomb"] += record["energy_coulomb_kj_mol"]
    result["forces_kj_mol_nm"] = (
        np.asarray(native["forces_kj_mol_nm"]) + record["forces_kj_mol_nm"]
    ).tolist()
    result["configurational_pressure_bar"] += record["pressure_increment_bar"]
    tensor = np.asarray(record["virial_kj_mol"])
    # No fabricated full tensor: expose the independently calculated increment.
    result["special_pair_virial_increment_kj_mol"] = tensor.tolist()
    result["special_pair_completion_identity"] = record["identity"]
    return result


def comparison(openmm, native, record):
    return {
        "native": compare(openmm, native),
        "composite": compare(openmm, corrected(native, record)),
    }
