"""Summarize actual bounded calculations and the explicit Si-H scope decision."""

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.angular_baseline import source_extensions
from island.pcff_qm_pilot.angular_fit import Problem
from scripts.audit_phase_4j62 import build


def read(p):
    return json.loads(Path(p).read_text())


def geometry_audit(root):
    source = read(root / "HMDS-domain-audit.json")
    center = read(root / "HMDS-center.json")
    x = np.array(center["center"]["coordinates"])
    order = source["atom_order"]
    require(
        center["center"]["graph_identity"] == source["graph_identity"]
        and center["center"]["atom_order"] == order,
        "HMDS_final_graph_mapping_changed",
    )
    bonded = {tuple(sorted((i, j))) for i, j, _ in source["bonds"]}
    rows, categories = [], defaultdict(list)
    for i, j in sorted(bonded):
        pair = "-".join(sorted((source["elements"][i - 1], source["elements"][j - 1])))
        distance = float(np.linalg.norm(x[i - 1] - x[j - 1]))
        categories[pair].append(distance)
        rows.append({"sites": [i, j], "length_A": distance})
    other = sorted(
        (float(np.linalg.norm(x[i - 1] - x[j - 1])), i, j)
        for i in order
        for j in order
        if i < j and (i, j) not in bonded
    )
    require(
        x.shape == (27, 3) and np.isfinite(x).all() and other[0][0] > 0.8,
        "HMDS_final_geometry_clash",
    )
    result = signed(
        {
            "graph_identity": source["graph_identity"],
            "source_graph_sha256": pin(root / "HMDS-domain-audit.json"),
            "geometry_identity": digest(x.tolist()),
            "atom_order": order,
            "elements": source["elements"],
            "formal_charge": 0,
            "multiplicity": 1,
            "source_valence_and_environments_preserved": True,
            "bond_lengths": rows,
            "category_ranges_A": {k: [min(v), max(v)] for k, v in categories.items()},
            "closest_nonbonded_pair": {
                "distance_A": other[0][0],
                "sites": other[0][1:],
            },
            "optimization_converged": center["optimization_result"]["status"]
            == "completed",
            "stationary_candidate": center["stationary_candidate"],
            "interpretation": "Distance/valence compatibility screen of the unchanged explicit graph; not a quantum bond-order determination or minimum proof.",
        }
    )
    write(root / "HMDS-geometry-validation.json", result)
    return result


def finalize(repo, root, work):
    geometry_audit(root)
    ledger = read(work / "ledger.json")["calls"]
    require(len(ledger) <= 48, "J62_budget_exceeded")
    groups = Counter(k.split("/")[0] for k in ledger)
    statuses = Counter(v["status"] for v in ledger.values())
    native = [
        read(work / "outputs" / k / "result.json")
        for k, v in ledger.items()
        if v["status"] == "completed"
    ]
    write(
        root / "compute-accounting.json",
        signed(
            {
                "schema": "island.j62-compute-accounting.v1",
                "attempted_energy_gradient_evaluations": len(ledger),
                "completed_gradients": statuses["completed"],
                "statuses": dict(statuses),
                "groups": dict(groups),
                "primary_control_gradients": groups["control"],
                "sensitivity_gradients": groups["sensitivity"],
                "HMDS_optimization_gradients": groups["hmds-opt"],
                "HMDS_optimizer_status": read(work / "hmds-optimization/result.json")[
                    "status"
                ],
                "gradient_elapsed_seconds": {
                    "sum": sum(r["elapsed_seconds"] for r in native),
                    "minimum": min(r["elapsed_seconds"] for r in native),
                    "maximum": max(r["elapsed_seconds"] for r in native),
                },
                "limit": 48,
                "numerical_Hessian_hidden_gradient_calls": 0,
                "full_J57_inventory_executed": False,
                "native_file_inventory_sha256": pin(
                    root / "native-file-inventory.json"
                ),
                "failed_attempts_preserved": True,
                "all_native_processes_reaped": True,
            }
        ),
    )
    control = read(root / "angular-torsional-control.json")
    base, initial = build(repo)
    fine = control["steps"]["0.002"]
    problem = Problem(
        base,
        control["Q_combined"],
        fine["combined_symmetric_projected"],
        fine["Cartesian_response"],
    )
    x = np.r_[initial, np.zeros(4)]
    h, _, r = problem.curvature(x)
    target = np.array(fine["combined_symmetric_projected"])
    ref = np.array(fine["Cartesian_response"])
    compare = {
        "unchanged_J61_soft_optimum": {
            "combined_curvature_relative_error": float(
                np.linalg.norm(h - target) / np.linalg.norm(target)
            ),
            "nonstretch_column_RMS_errors": np.sqrt(
                np.mean((r[:, 8:] - ref) ** 2, axis=0)
            ).tolist(),
            "MM_projected_curvature": h.tolist(),
        }
    }
    for name in ("M0", "M1", "M2", "M3"):
        t = read(root / "nested" / f"{name}.json")
        compare[name] = {
            k: t[k]
            for k in (
                "parameter_count",
                "rank",
                "column_scaled_condition",
                "active_bounds",
                "weighted_SS",
                "training",
                "stationary_force_RMS",
                "projected_curvature_relative_error",
                "nonstretch_response_column_RMS_errors",
                "residual_covariance_sigma_proxy",
                "training_gates",
            )
        }
    write(
        root / "research-model-provenance.json",
        signed(
            {
                "model": "island_research_class2_sih_disiloxane_j62_nested_v1",
                "scope": "exact nine-atom Si-H disiloxane neutral singlet only",
                "origin": "fitted_research_extension_trials",
                "original_native_PCFF_source_rows": False,
                "source_rows_overwritten": False,
                "parameter_names": problem.names,
                "lower_bounds": problem.lower.tolist(),
                "upper_bounds": problem.upper_bounds.tolist(),
                "J61_units_record": "../phase_4j61/parameter-specifications.json",
                "extension_units": [
                    "kcal/mol",
                    "radian",
                    "radian",
                    "kcal/mol/radian^2",
                ],
                "scientific_plan_sha256": pin(root / "plan/scientific-plan.json"),
                "frozen_model_implementation_sha256": pin(
                    root / "plan/model-implementation-pins.json"
                ),
                "actual_native_reference_sha256": pin(
                    root / "normalized-reference.json"
                ),
                "known_source_Hamiltonian_sha256": pin(
                    repo / "docs/evidence/phase_4j60/hamiltonian-definition.json"
                ),
                "trials": {
                    name: pin(root / "nested" / f"{name}.json")
                    for name in ("M0", "M1", "M2", "M3")
                },
                "uncertainty_record": "nested-uncertainty.json",
                "accepted_coefficients": {},
                "validation_used_for_parameter_selection": False,
            }
        ),
    )
    amplitude, k2 = source_extensions(base.model)
    write(
        root / "control-diagnosis.json",
        signed(
            {
                "schema": "island.j62-control-diagnosis.v1",
                "total_QM_and_total_MM_compared": True,
                "comparison": compare,
                "incremental_training_only_evidence": {
                    current: {
                        "compared_to": prior,
                        "weighted_SS_reduction_fraction": 1
                        - compare[current]["weighted_SS"]
                        / compare[prior]["weighted_SS"],
                        "force_RMS_change_kcal_mol_A": compare[current]["training"][
                            "force_RMS"
                        ]
                        - compare[prior]["training"]["force_RMS"],
                        "energy_RMSE_change_kcal_mol": compare[current]["training"][
                            "energy_RMSE"
                        ]
                        - compare[prior]["training"]["energy_RMSE"],
                        "stationarity_change_kcal_mol_A": compare[current][
                            "stationary_force_RMS"
                        ]
                        - compare[prior]["stationary_force_RMS"],
                    }
                    for prior, current in (("M0", "M1"), ("M1", "M2"), ("M2", "M3"))
                },
                "infeasibility_proved": False,
                "verified_source_main_n3_kcal_mol": amplitude,
                "verified_source_central_K2_kcal_mol_rad2": k2,
                "control_eigenvalues_kcal_mol_A2": np.linalg.eigvalsh(target).tolist(),
                "interpretation": "Nested improvements discriminate the tested parameter directions only. A molecular response is not an isolated source coefficient; projection fractions overlap, and remaining error does not uniquely identify its physical cause. Fixed charges, exclusions, pair settings and phase conventions were not adjusted.",
                "coverage": "12 of21 internal directions at the J61 geometry; four projected nonstretch directions plus8 stretches. Not a verified full minimum.",
                "J61_conclusions_rewritten": False,
            }
        ),
    )
    selection = read(root / "nested-model-selection.json")
    require(
        selection["chosen_for_conditional_fresh_validation"] is None,
        "finalize_conditional_validation_before_scope_decision",
    )
    hmds = read(root / "HMDS-reference-results.json")
    write(
        root / "HMDS-next-experiment-decision.json",
        signed(
            {
                "graph_identity": read(root / "HMDS-domain-audit.json")[
                    "graph_identity"
                ],
                "observed_basis_contrast": hmds["SVP_TZVP_sensitivity"],
                "center_status": "unfinished SVP optimization; not stationary",
                "data_integrity": "21 native gradients reconstructed independently",
                "precision_not_established": [
                    "two-step differentiation",
                    "full response grid/SCF precision",
                    "basis convergence beyond SVP/TZVP",
                    "QM-method adequacy",
                    "full molecular stability",
                ],
                "smallest_next_numerical_dependency": "At the exact retained HMDS graph, use the last SVP geometry as an initial guess for a separately budgeted PBE0/def2-TZVP stationary-candidate optimization. The measured 3.84 kcal/mol/A gradient-RMS shift makes completion of the SVP optimization alone insufficient. Retain all callbacks; do not presume convergence or final TZVP adequacy.",
                "next_response_data": "At a genuinely method-matched candidate, take a second predeclared step for full Si-C, C-H and methyl/Si-O torsional Cartesian responses, and representative grid/SCF and larger-basis checks. The present one-step Cartesian torsion-seed response includes angular and other molecular contributions and is not an isolated torsion constant.",
                "independent_MM_dependency": "Bind all fixed/nuisance methyl contributions and reduce the 26 directed missing source queries only by verified source/type symmetries; then select joint bond-angle-angle-torsional deformations from the complete Jacobian and reserve independent deformed conformers. Missing MM coefficients do not prohibit the stationary-candidate QM calculation.",
                "next_interior_environment": "Separately scoped capped-DP2: sio with two osi and two c3 neighbors; O-Si-O and Si-O-Si-O interior responses not supplied by HMDS.",
                "SiH_transfer": False,
                "additional_J62_calculations_activated": False,
                "basis_label_claim": "def2-TZVP is a cost-informed next reference candidate, not a converged final methyl protocol or complete-basis result",
            }
        ),
    )
    write(
        root / "fitting-decision.json",
        signed(
            {
                "schema": "island.j62-fitting-decision.v1",
                "SiH_scope_decision": "Stop at the bounded angular/torsional milestone; no automatic additional Si-H campaign.",
                "training_only_feasibility_witness": read(
                    root / "training-weighting-tradeoff.json"
                )["stationarity_and_curvature_feasible_witness_found"],
                "reference_precision_passed": control["reference_precision_passed"],
                "training_qualifiers": selection["training_qualifiers"],
                "fresh_validation_status": "not_run_no_qualifying_model",
                "fresh_validation_reserved": 3,
                "trial_parameter_counts": [24, 25, 27, 28],
                "accepted_coefficients": {},
                "fit_completion": "bounded nested trial fits completed; no scientific acceptance",
                "scientific_validation": "not_established",
                "parameter_completeness": "incomplete",
                "typing_status": "complete",
                "charge_status": "complete",
                "QM_execution": "passed",
                "QM_data_integrity": "passed",
                "fit_completed": True,
                "fit_validated": False,
                "production_readiness": "not_established",
                "lammps_export_status": "not_run",
                "lammps_runtime_status": "not_run",
                "openmm_build_status": "not_run",
                "cross_engine_equivalence": "not_claimed",
                "HMDS_QM_execution": "completed exploratory gradients",
                "HMDS_stationary_candidate": hmds["center_stationary_candidate"],
                "HMDS_MM_fit": "not_run_missing_couplings_and_reference_precision",
                "methyl_transfer_from_SiH": False,
                "next_dependency_record": "HMDS-next-experiment-decision.json",
                "chemical_scope_record": "chemical-scope-and-next-data.json",
                "user_attested_confirmations": 28,
                "original_author_artifact_status": "not_yet_supplied",
                "PCFF_MD_executed": False,
                "new_registered_profiles": [],
                "production_validated": False,
                "simulation_readiness": "not_established",
                "ordinary_API": 'forcefield="pcff"',
            }
        ),
    )
    return {
        "new_gradients": len(ledger),
        "accepted_coefficients": 0,
        "HMDS_stationary_candidate": hmds["center_stationary_candidate"],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j62"))
    a = p.parse_args()
    print(json.dumps(finalize(a.repo, a.evidence, a.work), sort_keys=True))
