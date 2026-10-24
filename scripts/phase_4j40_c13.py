"""Conservative numerical characterization; no new acceptance bounds are fitted."""

import numpy as np

from scripts.phase_4j40_c11_numeric import compare, require, signed
from scripts.phase_4j40_c12 import corrected

POLICY = signed(
    {
        "schema": "island.j40-c13-policy.v1",
        "status": "candidate_not_qualified",
        "energy_absolute_kj_mol": 1e-5,
        "force_absolute_kj_mol_nm": 1e-4,
        "force_relative": 2e-10,
        "pressure_bound_bar": None,
        "justification": "Preserve J22/J40 bounds; no evidence-derived relaxation. A pressure bound is not established.",
        "classes": [
            "bitwise_equality",
            "roundoff_level_equality",
            "bounded_numerical_equivalence",
            "unresolved_backend_discrepancy",
        ],
        "roundoff_policy": "not established for heterogeneous engine algorithms",
        "bounded_equivalence_policy": "not established without independent convergence and repeated-frame qualification",
    }
)


def characterize(row, frame):
    require(row["frame_identity"] == frame["identity"], "Frame identity mismatch")
    require(signed(frame) == frame, "Mutated frame")
    a = row["openmm"]
    b = corrected(row["lammps"], row["completion"])
    m = compare(a, b)
    raw_a, raw_b = a["raw_groups_kj_mol"], b["raw_groups_kj_mol"]
    identical = (
        np.float64(a["energy_kj_mol"]).tobytes()
        == np.float64(b["energy_kj_mol"]).tobytes()
        and a["components_kj_mol"] == b["components_kj_mol"]
        and np.asarray(a["forces_kj_mol_nm"], dtype="<f8").tobytes()
        == np.asarray(b["forces_kj_mol_nm"], dtype="<f8").tobytes()
    )
    return signed(
        {
            "name": row["name"],
            "frame_identity": frame["identity"],
            "atoms": len(frame["atom_ids"]),
            "box_lengths_nm": frame["box_lengths_nm"],
            "metrics": m,
            "classification": "bitwise_equality"
            if identical
            else "unresolved_backend_discrepancy",
            "existing_energy_force_bounds_passed": m["passed"],
            "real_space_signed_difference_kj_mol": raw_a["coulomb_real_plus_exception"]
            - raw_b["real"],
            "reciprocal_signed_difference_kj_mol": raw_a[
                "coulomb_reciprocal_including_self_background"
            ]
            - raw_b["reciprocal"],
            "decomposition_caveat": "Native alpha and exception decompositions differ; separate pieces are descriptive, not equal observables. Native raw groups exclude analytical completion.",
            "virial_tensor_error": None,
            "scalar_pressure_error_bar": None,
            "tensor_pressure_error_bar": None,
            "pressure_diagnostic": "pressure_unavailable: OpenMM full periodic reciprocal virial not independently measured",
            "settings": {"openmm": a["setting"], "lammps": b["setting"]},
            "mesh": a["mesh"],
            "alpha_nm": a["alpha_nm"],
            "comparison_semantics": b["engine"],
        }
    )


def validate_summary(value):
    require(value["schema"] == "island.j40-c13-study.v1", "Wrong summary schema")
    require(signed(value) == value, "Mutated summary")
    require(value["policy"] == POLICY, "Mutated policy")
    require(bool(value["rows"]), "Empty matrix")
    for row in value["rows"]:
        require(signed(row) == row, "Mutated characterization")
        require(
            row["classification"]
            in ("bitwise_equality", "unresolved_backend_discrepancy"),
            "Unsupported qualification",
        )
    for key, expected in {
        "ewald_status": "ewald_unvalidated",
        "pressure_status": "pressure_unavailable",
        "cutoff_recovery": "not_executed",
        "long_npt": "not_executed",
        "production_validated": False,
        "simulation_readiness": "not_established",
        "scientific_acceptance": False,
    }.items():
        require(value[key] == expected, "Unsupported readiness: " + key)
    return True
