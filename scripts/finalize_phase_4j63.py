"""Summarize executed HMDS work and the numerical prerequisite for joint fitting."""

import argparse
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, signed, write
from scripts.analyze_phase_4j63 import gate, read


def finalize(root, work):
    ledger = read(work / "ledger.json")["calls"]
    statuses = Counter(r["status"] for r in ledger.values())
    native = [
        read(work / "outputs" / k / "result.json")
        for k, v in ledger.items()
        if v["status"] == "completed"
    ]
    center = read(root / "stationary-center.json")
    paths = (
        read(root / "deformation-reference-results.json")
        if (root / "deformation-reference-results.json").exists()
        else None
    )
    design = (
        read(root / "joint-design.json")
        if (root / "joint-design.json").exists()
        else None
    )
    fixed = (
        read(root / "source-baseline-compatibility.json")
        if (root / "source-baseline-compatibility.json").exists()
        else None
    )
    trajectory_path = root / "optimizer-model-diagnosis.json"
    trajectory = read(trajectory_path) if trajectory_path.exists() else None
    curvature_path = root / "source-baseline-path-curvature.json"
    curvature = read(curvature_path) if curvature_path.exists() else None
    reservation = root / "plan/reserved-validation.json"
    reserved_count = (
        len(read(reservation)["coordinates"]) if reservation.exists() else 0
    )
    failed = []
    if paths:
        failed += [
            f"two-step/chain-rule/{k}"
            for k, v in paths["primary"].items()
            if not v["precision_passed"]
        ]
        for control, row in paths["controls"].items():
            failed += [
                f"{control}/{k}" for k, v in row["paths"].items() if not v["passed"]
            ]
    obstacles = []
    if failed:
        obstacles.append("Frozen reference objectives fail: " + ", ".join(failed))
    if design and not design["current"]["identifiable_for_screen"]:
        obstacles.append(
            f"Complete primary Jacobian rank {design['current']['rank']} has scaled "
            f"condition {design['current']['scaled_condition']:.6g}, above 1000."
        )
    if curvature:
        for name, row in curvature["paths"].items():
            fine = row["steps"]["0.003"]
            obstacles.append(
                f"{name}: total QM minus fixed-model curvature "
                f"{fine['QM_minus_model_curvature']:.6g} kcal/mol/A2, while all 56 "
                "unknown coupling second path derivatives are independently zero. "
                f"TZVPP minus TZVP curvature {row['basis_polarization_curvature_shift']:.6g}."
            )
    precision_entry = {
        "question": "Resolve the soft siloxane energy/gradient discrepancy at identical geometries, without changing the stationary center or acceptance threshold.",
        "exact_geometries": [
            "stationary-center",
            "paths/siloxane-h0.02-s+1",
            "paths/siloxane-h0.02-s-1",
        ],
        "next_comparison": "A separately budgeted PBE0/def2-TZVP unpruned250x1454-grid center and +/-0.02rad pair, using the retained tighter-SCF setting1e-12/1e-10. Compare with the retained150x974 tighter-SCF triple to isolate the grid effect. Preserve both original failed gates.",
        "gradient_calls": 3,
        "execution_status": "not_activated",
        "scope": "This tests one local soft response; not full molecular stability, method adequacy or all-objective reference closure.",
    }
    if paths:
        scf_response = paths["controls"]["scf"]["paths"]["siloxane"]["control_response"]
        rule = read(root / "plan/scientific-plan.json")["precision"][
            "energy_gradient_chain_rule"
        ]
        scf_chain = gate(
            scf_response["curvature_from_energies"],
            scf_response["curvature_from_path_gradient"],
            rule["rotation_curvature_floor"],
            rule["relative"],
        )
        precision_entry["retained_tighter_SCF_energy_gradient_check"] = scf_chain
        if scf_chain["passed"]:
            precision_entry.update(
                exact_geometries=[
                    "paths/siloxane-h0.04-s+1",
                    "paths/siloxane-h0.04-s-1",
                ],
                next_comparison="The retained tighter-SCF fine triple passes the original energy/gradient allowance. Obtain only the two coarse0.04rad gradients at the same150x974 grid and SCF1e-12/1e-10, then reconstruct a separately identified two-step protocol. Do not mix settings into the original dataset or erase its failed gates.",
                gradient_calls=2,
            )
    next_data = signed(
        {
            "schema": "island.j63-specific-HMDS-next-experiment.v1",
            "stationary_candidate_obtained": center["stationary_candidate"],
            "current_failed_reference_objectives": failed,
            "precision_entry_experiment": precision_entry,
            "basis_entry_experiment": {
                "question": "Does the local Cartesian-force polarization sensitivity persist with the inventoried1079-function def2-QZVP orbital basis?",
                "exact_geometry": "stationary-center",
                "settings": "PBE0, same150x974 grid, same SCF1e-10/1e-8, same def2-universal-jkfit, no dispersion; only orbital basis changes.",
                "gradient_calls": 1,
                "execution_status": "not_activated; requires a new cost/memory plan, no automatic basis ladder",
                "scope": "A center-force discriminant, not a completed deformation-energy/force basis gate. Any later paired paths require their own bounded plan.",
            },
            "baseline_requirement_before_more_coupling_QM": "Resolve the demonstrated total-model Si-C/C-H curvature gaps through independently justified source/baseline accounting or a separately identified and preregistered HMDS baseline-adjustable research variant. Missing couplings cannot change these path curvatures. Do not automatically release rows or interpret a molecular curvature as an isolated bond coefficient.",
            "conditional_joint_design": []
            if design is None
            else [r["id"] for r in design["selected_additional_configurations"]],
            "conditional_joint_coordinates": "joint-design.json; geometry-only selection, no returned QM",
            "fresh_validation_reservation": "plan/reserved-validation.json"
            if reserved_count
            else None,
            "fresh_validation_QM_supplied": False,
            "coupling_data_alone_resolve_baseline_gap": False,
            "capped_DP2_separate_requirement": "New exact capped-DP2 graph, typing/charge/source closure and responses for sio(osi,osi,c3,c3), O-Si-O and Si-O-Si-O; terminal HMDS cannot identify those environments.",
            "new_QM_or_baseline_parameters_activated": False,
        }
    )
    write(root / "next-data-requirements.json", next_data)
    write(
        root / "compute-accounting.json",
        signed(
            {
                "schema": "island.j63-compute-accounting.v1",
                "attempted_energy_gradient_evaluations": len(ledger),
                "completed_gradients": statuses["completed"],
                "statuses": dict(statuses),
                "groups": dict(Counter(k.split("/")[0] for k in ledger)),
                "limit": 80,
                "optimizer_limit": 40,
                "optimizer_status": center["optimizer_status"],
                "numerical_process_returncode": center["numerical_process_returncode"],
                "numerical_Hessian_calls": 0,
                "new_SiH_QM_calls": 0,
                "J57_full_inventory_executed": False,
                "sum_gradient_seconds": sum(v["elapsed_seconds"] for v in native),
                "gradient_seconds_range": [
                    min(v["elapsed_seconds"] for v in native),
                    max(v["elapsed_seconds"] for v in native),
                ],
                "native_inventory_sha256": pin(root / "native-file-inventory.json"),
                "failed_interrupted_attempts_preserved": True,
            }
        ),
    )
    decision = signed(
        {
            "schema": "island.j63-HMDS-fitting-decision.v1",
            "graph_identity": center["center"]["graph_identity"],
            "scope": "Exact 27-atom neutral singlet HMDS, terminal sio(osi,c3,c3,c3) and c3(sio,hc,hc,hc) only. No methyl transfer from Si-H.",
            "stationary_candidate": center["stationary_candidate"],
            "minimum_verified": False,
            "functional_accounting": "complete research Hamiltonian with 56 explicit unknown coefficients; no defaults for missing physical values",
            "independent_coefficient_count": 56,
            "directed_missing_queries": 26,
            "QM_execution": "completed"
            if statuses["completed"] == len(ledger)
            else "contains_unsuccessful_attempts",
            "QM_data_integrity": "independently_reconstructed",
            "optimizer_status": center["optimizer_status"],
            "optimizer_process_returncode": center["numerical_process_returncode"],
            "bounded_reference_precision_passed": False
            if paths is None
            else paths["bounded_precision_passed"],
            "failed_reference_objectives": failed,
            "actual_design_rank": None if design is None else design["current"]["rank"],
            "actual_design_scaled_condition": None
            if design is None
            else design["current"]["scaled_condition"],
            "actual_design_identifiable": False
            if design is None
            else design["current"]["identifiable_for_screen"],
            "proposed_additional_design_count": None
            if design is None
            else len(design["selected_additional_configurations"]),
            "source_baseline_optimistic_joint_lower_bound": None
            if fixed is None
            else fixed["optimistic_joint_column_space_lower_bound"],
            "optimizer_development_lower_bound": None
            if trajectory is None
            else trajectory["optimistic_unconstrained_lower_bound"],
            "unknown_null_stretch_curvature": None
            if curvature is None
            else {
                name: {
                    k: row[k]
                    for k in (
                        "QM_total_path_curvature",
                        "model_total_path_curvature_independent_of56_unknowns",
                        "QM_minus_model_curvature",
                        "relative_curvature_difference",
                    )
                }
                | {"TZVPP_minus_TZVP": entry["basis_polarization_curvature_shift"]}
                for name, entry in curvature["paths"].items()
                for row in [entry["steps"]["0.003"]]
            },
            "fit_readiness": "not_established",
            "fit_completed": False,
            "fit_validated": False,
            "accepted_coefficients": {},
            "rejected_physical_candidates": [],
            "reason": " ".join(obstacles)
            + " No physical fit or validation is inferred from a column-space projection or optimizer convergence.",
            "next_data_requirements_identity": next_data["identity"],
            "fresh_validation": "graph-screened new deformed configurations reserved, QM not executed, not used for weights/bounds"
            if reserved_count
            else "not_activated_stationary_stage_prerequisite_unsatisfied",
            "fresh_validation_reserved_count": reserved_count,
            "baseline_rows_released": False,
            "source_charges_modified": False,
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
            "user_attested_confirmations": 28,
            "original_author_artifact_status": "not_yet_supplied",
            "ordinary_API": 'forcefield="pcff"',
            "next_interior_domain": "Capped-DP2 sio(osi,osi,c3,c3), interior O-Si-O and Si-O-Si-O require separate graph-bound sources and reference data; HMDS cannot constrain them.",
        }
    )
    write(root / "fitting-decision.json", decision)
    return decision


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    p.add_argument("--work", type=Path, required=True)
    a = p.parse_args()
    print(finalize(a.evidence, a.work)["failed_reference_objectives"])
