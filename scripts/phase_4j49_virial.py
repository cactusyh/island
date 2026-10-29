"""Affine virial and precise cause ledger; pressure remains a diagnostic proxy."""

import argparse
import json
from pathlib import Path

import numpy as np

from island.charge_references.records import unpack
from island.forcefields.opls_registry import (
    OPLSAssets,
    prepare_polymer_forcefield,
    use_opls_assets,
)
from island.forcefields.selection.inventory import digest
from island.pdms import restore_pdms_system
from scripts.phase_4j47 import Reference, write
from scripts.qualify_pdms_public import LMP_COULOMB, MM_COULOMB


def audit(root, assets):
    root = Path(root)
    records = {}
    for dp in (0, 2, 4, 6):
        system = restore_pdms_system(unpack((root / f"dp{dp}/system.json").read_text()))
        with use_opls_assets(assets):
            prepared = prepare_polymer_forcefield(system, forcefield="opls")
        backend = prepared.validate_integrity(system).validate_integrity(system)
        xyz = np.array(
            [system.coordinates.get(a["site_id"]) for a in backend.payload["mapping"]]
        )
        ref = Reference(backend.payload, xyz)
        try:
            virial = ref.virial()
        finally:
            ref.close()
        modes = json.loads((root / f"components-dp{dp}.json").read_text())
        records[str(dp)] = {
            "selected_profile": prepared.payload["selected_profile"],
            "OpenMM_affine_virial_at_1e-4_and_5e-5": virial,
            "modes": {
                name: {
                    "LAMMPS_virial_kj_mol": v["LAMMPS"][
                        "virial_kj_mol_xx_yy_zz_xy_xz_yz"
                    ],
                    "absolute_errors_kj_mol": [
                        abs(a["values_at_1e-4_and_5e-5"][-1] - b)
                        for a, b in zip(
                            virial, v["LAMMPS"]["virial_kj_mol_xx_yy_zz_xy_xz_yz"]
                        )
                    ],
                    "configurational_pressure_proxy_atm": v["LAMMPS"][
                        "configurational_pressure_proxy_atm"
                    ],
                    "pressure_validated": False,
                    "NPT_permitted": False,
                }
                for name, v in modes.items()
            },
            "virial_acceptance_bound": None,
            "reason": "No inherited virial/pressure scientific gate; raw measured residuals do not authorize pressure/NPT",
        }
    result = {
        "schema": "island.j49-discrepancy-diagnosis.v1",
        "cases": records,
        "generic_charged_controls_also_fail_author_numerics": True,
        "PDMS_charge_mapping_not_root_cause": True,
        "classification": {
            "bonded_angle_torsion_LJ": "precision",
            "author_PPPM_vs_OpenMM_Ewald": "Ewald/PPPM convention",
            "mesh_convergence": "reciprocal tolerance",
            "pair_table_energy_force_interpolation": "precision",
            "real_erfc_polynomial_energy_derivative_inconsistency": "precision",
            "Coulomb_constant_residual": "unit conversion",
            "real_vs_reciprocal_partition": "Ewald/PPPM convention",
            "source_author_replacement_numerical_contract": "source ambiguity",
            "pressure_scientific_gate": "unresolved",
        },
        "verified_equal_inputs": [
            "final graph",
            "atom order",
            "coordinates",
            "box",
            "masses",
            "charge literals",
            "1-2 exclusions",
            "1-3 exclusions",
            "1-4 LJ and Coulomb 0.5",
            "geometric sigma/epsilon",
            "LJ r-squared switching 9-11 A",
            "Coulomb cutoff 11 A",
            "dielectric 1",
            "periodic conducting boundary",
            "energy/force output units",
        ],
        "measured_unequal_numerics": {
            "author": {
                "solver": "PPPM",
                "tolerance": 1e-4,
                "table": "default",
                "erfc": "Abramowitz-Stegun polynomial",
            },
            "OpenMM": {
                "solver": "Ewald",
                "tolerance": 1e-7,
                "erfc": "libm; exact analytic equation",
            },
            "exact_kernel_control": {
                "solver": "Ewald",
                "tolerance": 1e-7,
                "alpha_per_angstrom": float(np.sqrt(-np.log(2e-7)) / 11),
                "reciprocal_max": [20, 20, 20],
                "table": 0,
                "erfc": "libm; explicitly modified executable",
            },
        },
        "Coulomb_constants_kj_mol_angstrom_e2": {
            "OpenMM": MM_COULOMB,
            "LAMMPS": LMP_COULOMB,
        },
        "no_energy_or_force_output_corrections": True,
        "no_charge_or_parameter_changes": True,
        "raw_diagnostic_static_and_workflow_passed": True,
        "original_author_settings_still_blocked": True,
        "Ewald_PPPM_scientific_closure": False,
        "backend_pressure_available": False,
    }
    return {**result, "identity": digest(result)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    a = p.parse_args(argv)
    write(
        a.directory / "diagnosis.json",
        audit(a.directory, OPLSAssets(a.archive, a.archive_sha256)),
    )


if __name__ == "__main__":
    main()
