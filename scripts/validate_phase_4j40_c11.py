"""Reconstruct C11 statistics from fresh raw, unit-checked engine outputs."""

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.phase_4j40_c11_numeric import (
    compare,
    parse_dump,
    require,
    sha,
    signed,
    write,
)


def characterize(root):
    study = json.loads((root / "study.json").read_text())
    baseline = study["baseline"]
    a, b = baseline["openmm"], baseline["lammps"]
    metrics = compare(a, b)
    backend = json.loads((root / "inputs/backend.json").read_text())
    sites = sorted(backend["numerical_data"]["sites"], key=lambda s: s["id"])
    dump = parse_dump(root / "baseline/lammps/forces.dump")
    require(
        dump["charges_e"] == [s["charge"] for s in sites],
        "Changed native charge vector",
    )
    require(dump["forces_kj_mol_nm"] == b["forces_kj_mol_nm"], "Changed raw forces")
    old = json.loads((root / "inputs/historical-comparison.json").read_text())
    require(
        abs(a["energy_kj_mol"] - old["openmm_energy_kj_mol"]) <= 1e-12,
        "Historical Ewald baseline not reproduced",
    )
    charges = dump["charges_e"]
    rows = [
        {
            "name": r["name"],
            "setting": r["setting"],
            "comparison": compare(r["openmm"], r["lammps"]),
        }
        for r in study["rows"]
    ]
    for row, original in zip(rows, study["rows"], strict=True):
        require(
            row["comparison"] == original["comparison"], "Changed numerical comparison"
        )
    baseline_series = [
        r for r in rows if r["name"] == "baseline" or r["name"].startswith("tolerance-")
    ]
    return signed(
        {
            "schema": "island.j40-c11-characterization.v2",
            "system": "GAFF2 linear PE; two neutral C8H18 chains",
            "graph_identity": backend["final_graph_identity"],
            "backend_identity": backend["identity"],
            "charge_identity": backend["charge_identity"],
            "charge_origin": "RESP/Psi4 exact graph; not AM1-BCC or native OPLS-AA",
            "atom_count": len(sites),
            "box_lengths_nm": dump["box_lengths_nm"],
            "cutoff_nm": 0.5,
            "switching": "disabled",
            "neighbor_skin_nm_lammps": 0.2,
            "baseline_electrostatics": "OpenMM Ewald versus LAMMPS Ewald (not PME)",
            "total_energy_kj_mol": a["energy_kj_mol"],
            "coulomb_energy_kj_mol": a["components_kj_mol"]["coulomb"],
            "absolute_energy_error_kj_mol": metrics["absolute_energy_error_kj_mol"],
            "relative_energy_error": metrics["relative_energy_error"],
            "force_errors": metrics["forces"],
            "charge_range_e": [min(charges), max(charges)],
            "charge_absolute_range_e": [min(map(abs, charges)), max(map(abs, charges))],
            "maximum_lammps_atom_force_kj_mol_nm": float(
                max(np.linalg.norm(dump["forces_kj_mol_nm"], axis=1))
            ),
            "tolerance_series": baseline_series,
            "all_comparisons": rows,
            "pressure_strain": study["pressure_strain"],
            "raw_baseline_dump_sha256": sha(root / "baseline/lammps/forces.dump"),
            "local_timestep_error_proxy": {
                "timestep_ps": 0.0005,
                "max_single_step_velocity_difference_nm_ps": float(
                    np.max(
                        np.abs(
                            np.asarray(a["forces_kj_mol_nm"])
                            - np.asarray(b["forces_kj_mol_nm"])
                        )
                        / np.asarray([s["mass"] for s in sites])[:, None]
                    )
                    * 0.0005
                ),
                "max_single_step_position_difference_nm": float(
                    0.5
                    * np.max(
                        np.abs(
                            np.asarray(a["forces_kj_mol_nm"])
                            - np.asarray(b["forces_kj_mol_nm"])
                        )
                        / np.asarray([s["mass"] for s in sites])[:, None]
                    )
                    * 0.0005**2
                ),
                "limitation": "instantaneous force-difference acceleration proxy; not an integrated trajectory error bound or MD validation",
            },
            "units_correction": "historical kJ/(mol Angstrom) force error multiplied by 10, not divided; raw real-unit forces multiplied by 4.184*10",
            "scope": "fixed system size; box changes are strained single-point controls, not additional polymer sizes; no system-size bound established",
            "classification": "tolerance_sensitive_discrepancy_with_unresolved_residual_and_geometry_dependence",
            "bounded_numerical_difference_established": False,
            "scientific_acceptance": False,
            "ewald_status": "ewald_unvalidated",
            "pressure_status": "pressure_unavailable",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    write(a.output, characterize(a.study))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
