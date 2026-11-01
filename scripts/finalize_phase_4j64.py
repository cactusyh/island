"""Summarize executed HMDS observations without promoting conditional candidates."""

import argparse
import re
import tempfile
from decimal import Decimal
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, signed, write
from island.pcff_qm.reference import CONVERSION
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_precision import read


def summarize(root):
    baseline = read(root / "baseline-calibration.json")
    dense = read(root / "dense-grid-decision.json")
    qz = read(root / "QZVP-center-decision.json")
    design = read(root / "updated-joint-design.json")
    compute = read(root / "compute-accounting.json")
    # Keep this repeat separate from the exact matched grid-only contrast.
    refs = read(root / "normalized-reference.json")["samples"]
    independent = read(root / "independent-native-checks.json")["samples"]
    pair = ["dense/center", "exact/center"]
    energies = [
        Decimal(independent[k]["native_total_energy_hartree_decimal"]) for k in pair
    ]
    shift = float(energies[1] - energies[0]) * CONVERSION["energy_factor"]
    force_delta = np.array(refs[pair[1]]["forces"]) - refs[pair[0]]["forces"]
    with tempfile.TemporaryDirectory(prefix="j64-grid-inventory-") as tmp:
        native = Path(tmp) / "native"
        unpack(root / "native-outputs.tar.xz", native)
        grids = {}
        for key in refs:
            file = native / "outputs" / key / "psi4.out"
            text = file.read_text()
            grids[key] = {
                "native_output_sha256": pin(file),
                **{
                    label: re.findall(re.escape(label) + r"\s*=\s*([^\n]+)", text)
                    for label in (
                        "Total Points",
                        "Total Blocks",
                        "Pruning Scheme",
                        "Radial Points",
                        "Spherical Points",
                    )
                },
            }
    write(root / "native-grid-inventory.json", signed({"jobs": grids}))
    write(
        root / "center-repeatability-control.json",
        signed(
            {
                "jobs": pair,
                "native_decimal_energies_hartree": list(map(str, energies)),
                "exact_minus_zero_path_energy_kcal_mol": shift,
                "energy_curvature_change_from_center_replacement_at0.02rad": -2
                * shift
                / 0.02**2,
                "maximum_input_coordinate_difference_A": float(
                    abs(
                        np.array(refs[pair[0]]["coordinates"])
                        - refs[pair[1]]["coordinates"]
                    ).max()
                ),
                "force_RMS_difference_kcal_mol_A": float(
                    np.sqrt(np.mean(force_delta**2))
                ),
                "maximum_local_force_difference": float(
                    np.linalg.norm(force_delta, axis=1).max()
                ),
                "grid_inventories": {k: grids[k] for k in pair},
                "scope": "One-input-ulp repeat/control with observed grid-screening inventory differences; not a pure arithmetic experiment or proof of the sole cause of a torsional discrepancy. Exact-center data alone enter matched main contrasts.",
            }
        ),
    )
    sensitivities = {}
    for name, row in baseline["paths"].items():
        estimates = row["conditional_delta_K2"]
        primary = estimates["TZVP_standard"]
        values = list(estimates.values())
        sensitivities[name] = {
            "protocol_estimates": estimates,
            "absolute_shift_from_TZVP_standard": {
                k: v - primary for k, v in estimates.items()
            },
            "protocol_range": [float(min(values)), float(max(values))],
            "maximum_absolute_protocol_shift": float(
                max(abs(v - primary) for v in values)
            ),
            "retained_two_step_estimates": baseline["retained_two_step_candidates"][
                name
            ],
            "precision_gates": {k: v["curvature"] for k, v in row["controls"].items()},
            "full_force_gates": {
                k: v["full_gate_passed"] for k, v in row["controls"].items()
            },
            "statistical_standard_error": None,
            "uncertainty_scope": "Finite-step and protocol sensitivities, correlated observations; not independent errors to add in quadrature, statistical confidence intervals or a complete-basis bound.",
        }
    precision = baseline["objective_status"]
    selected = [r["id"] for r in design["selected_geometry_only"]]
    next_data = signed(
        {
            "HMDS_reference": {
                "dense_energy_gradient_check_passed": dense[
                    "energy_gradient_consistency"
                ]["passed"],
                "dense_grid_only_full_check_passed": dense["grid_only_contrast"][
                    "full_gate_passed"
                ],
                "QZVP_center_force_differences": {
                    k: v["forces"] for k, v in qz["comparisons"].items()
                },
                "remaining_scope": "No all-path QZVP or complete-basis claim; force adequacy and soft-torsion precision remain separate from stretching-curvature calibration.",
                "smallest_discriminating_QZVP_numerical_control_before_full_force_fitting": {
                    "new_gradients": 1,
                    "geometry": "Exact retained J63 center; no optimization",
                    "method_basis": "PBE0/def2-QZVP",
                    "auxiliary_basis": "def2-universal-jkfit",
                    "grid": [250, 1454],
                    "SCF_thresholds": [1e-10, 1e-8],
                    "contrast": "The retained QZVP150x974 center with identical SCF, auxiliary basis and coordinates; isolates QZVP integration-grid force sensitivity before attributing the whole center-force shift to the orbital basis.",
                    "execution": "Not executed or added to the J64 campaign. Requires its own larger-grid resource preflight and frozen future plan; one such control would not establish orbital-basis convergence.",
                },
            },
            "baseline": {
                "conditional_calibration_completed": precision[
                    "baseline_calibration_completed"
                ],
                "physical_interpretation": "Two total-curvature corrections under the fixed research model; do not establish bare source-row errors or force qualification.",
                "next_step": "Keep equilibrium distances/charges/higher-order rows fixed. A later force-fit protocol must declare how these two independently constrained, protocol-specific curvature variables enter its objective and uncertainty; no additional parameter family is released here.",
                "smallest_independent_site_check": {
                    "specification": "plan/proposed-independent-site-check.json",
                    "new_gradients": 4,
                    "paths": ["SiC(6-7) +/-0.003angstrom", "CH(7-19) +/-0.003angstrom"],
                    "reuse_center": True,
                    "activation": "Only after the applicable stretching-curvature precision gates pass and a new validation-only protocol is frozen. These unknown-null observations test shared corrections without fitting56 couplings or releasing equilibria.",
                    "scope": "Other graph-screened sites in the same conformer, not a new conformational basin or complete force validation.",
                },
            },
            "joint_training": {
                "geometry_only_ids": selected,
                "additional_QM_count": len(selected),
                "rank": design["after_selected"]["rank"],
                "scaled_condition": design["after_selected"]["scaled_condition"],
                "execution": "Not activated. Conditioning is not reference or baseline validation. Choose and freeze one adequately checked protocol and energy/force/curvature objective before this campaign.",
                "existing_primary_geometries": 17,
            },
            "validation": {
                "separate_QM_count": 3,
                "reserved_file": "../phase_4j63/plan/reserved-validation.json",
                "still_unexecuted": True,
                "scope": "Fresh graph-screened deformations; not independently sampled conformational basins. No reference values used to select model, weights or bounds.",
            },
            "capped_DP2": "Separate exact graph/type/charge/dependency inventory and data for sio(osi,osi,c3,c3), O-Si-O and Si-O-Si-O; HMDS-only correction values do not supply those environments.",
        }
    )
    write(root / "next-data-requirements.json", next_data)
    decision = signed(
        {
            "Hamiltonian": baseline["Hamiltonian"],
            "counted_process_attempts": compute["counted_process_attempts"],
            "successful_native_gradients": compute["successful_gradients"],
            "objective_status": precision,
            "precision_scope": {
                "soft_torsion_reference_precision_flag": "Matched TZVP grid/SCF and energy-gradient consistency at the retained +/-0.02rad path only; not all-path orbital-basis or full-force adequacy.",
                "soft_torsion_basis_precision": "not_established; historical TZVP/TZVPP full-force failure unchanged, no new QZVP torsion endpoints",
                "stretching_reference_precision_flag": "The declared step, matched grid/SCF/basis curvature and energy-gradient criteria for two unknown-null observables; full Cartesian-force tests remain separately reported.",
                "reference_method": "PBE0 without dispersion, a declared research protocol; not author-confirmed final parameterization or complete-basis accuracy",
            },
            "conditional_candidates": baseline["conditional_candidates"],
            "candidate_origin": "fitted_external curvature-calibrated research-extension candidates",
            "accepted_scope": "Conditional total-curvature diagnostic for exact27-atom neutral-singlet HMDS at the explicitly listed protocol/path, only where its objective-specific gates pass.",
            "candidate_scientific_validation": "not_established",
            "force_qualified_candidates": [],
            "transferable_candidates": [],
            "automatic_assignment_candidates": [],
            "full56_coupling_fit_executed": False,
            "source_equilibrium_or_higher_rows_changed": False,
            "source_charges_or_pair_policy_changed": False,
            "parameter_completeness": "incomplete",
            "typing_status": "complete",
            "charge_status": "complete",
            "lammps_export_status": "not_run",
            "lammps_runtime_status": "not_run",
            "openmm_build_status": "not_run",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "PCFF_MD_executed": False,
            "new_registered_profiles": [],
            "SiH_investigation_reopened": False,
            "user_attested_confirmation_count": 28,
            "original_author_artifact_status": "not_yet_supplied",
            "physical56_coefficients_still_missing": 56,
            "conditional_candidate_checksum": baseline["identity"],
            "next_data_identity": next_data["identity"],
        }
    )
    write(root / "candidate-and-scope-decision.json", decision)
    write(
        root / "uncertainty-and-provenance.json",
        signed(
            {
                "sensitivity_records": sensitivities,
                "candidate_checksum": baseline["identity"],
                "normalized_native_reference_sha256": pin(
                    root / "normalized-reference.json"
                ),
                "raw_archive_sha256": pin(root / "native-outputs.tar.xz"),
                "data_license_scope": "Native calculations generated for this user-authorized research task; engine/source licenses and historical authority remain separately recorded.",
                "equilibrium_and_force_validation_supplied": False,
            }
        ),
    )
    return decision


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j64"))
    a = p.parse_args()
    result = summarize(a.evidence.resolve())
    print(result["objective_status"])
