"""Conservative derived-profile qualification, with no native registry overlay."""

from copy import deepcopy
from math import isfinite

from . import FAMILY, FLAGS, PROFILE, digest, require, signed
from .fitting import LIMITS
from .inventory import scope_check

# No physical extension is registered by the offline J78 phase.
QUALIFIED_BUNDLES = ()
GATES = (
    "complete_scope_inventory",
    "all_required_terms_resolved",
    "training_energy",
    "heldout_energy",
    "training_force",
    "heldout_force",
    "maximum_energy_error",
    "maximum_local_force_error",
    "equilibrium_geometry",
    "stationarity",
    "torsion_scan",
    "curvature",
    "conditioning",
    "parameter_uncertainty",
    "no_unstable_geometry",
    "no_unphysical_curvature",
    "charge_conservation",
    "exclusions_and_one_four_complete",
    "Class_II_coverage",
    "native_PCFF_immutable",
    "provenance_and_hashes",
    "independent_replay",
    "QM_protocol_adequacy",
    "source_and_model_uncertainty",
)


def evaluate_gates(evidence):
    """Numeric gates require actual values, never truthy 'passed' placeholders."""
    require(evidence["limits"] == LIMITS, "acceptance_threshold_mutation")
    gates = {k: False for k in GATES}
    for k in (
        "complete_scope_inventory",
        "all_required_terms_resolved",
        "no_unstable_geometry",
        "no_unphysical_curvature",
        "charge_conservation",
        "exclusions_and_one_four_complete",
        "Class_II_coverage",
        "native_PCFF_immutable",
        "provenance_and_hashes",
        "independent_replay",
        "QM_protocol_adequacy",
        "source_and_model_uncertainty",
    ):
        gates[k] = evidence.get(k) is True
    metric_names = {
        "training_energy": "energy_RMSE_kcal_mol",
        "heldout_energy": "energy_RMSE_kcal_mol",
        "training_force": "force_RMS_kcal_mol_A",
        "heldout_force": "force_RMS_kcal_mol_A",
        "maximum_energy_error": "maximum_energy_error_kcal_mol",
        "maximum_local_force_error": "maximum_local_force_error_kcal_mol_A",
        "stationarity": "stationarity_RMS_kcal_mol_A",
        "torsion_scan": "torsion_RMSE_kcal_mol",
        "curvature": "relative_curvature_error",
        "conditioning": "scaled_condition_max",
        "parameter_uncertainty": "uncertainty_fraction_of_bound_width_max",
    }
    for k, threshold in metric_names.items():
        value = evidence.get("metrics", {}).get(k)
        gates[k] = (
            type(value) in (int, float)
            and isfinite(value)
            and 0 <= value <= LIMITS[threshold]
        )
    geometry = evidence.get("geometry_errors", {})
    gates["equilibrium_geometry"] = all(
        type(geometry.get(k)) in (float, int)
        and isfinite(geometry[k])
        and 0 <= geometry[k] <= LIMITS[k]
        for k in ("equilibrium_bond_error_A", "equilibrium_angle_error_rad")
    )
    # Full-rank information must be provided independently of condition numbers.
    gates["conditioning"] &= (
        type(evidence.get("metrics", {}).get("conditioning")) in (int, float)
        and evidence["metrics"]["conditioning"] >= 1
        and type(evidence.get("rank")) is int
        and type(evidence.get("columns")) is int
        and evidence["rank"] == evidence["columns"]
        and evidence["columns"] > 0
    )
    passed = all(gates.values())
    return signed(
        {
            "gates": gates,
            "failed_gates": [k for k in GATES if not gates[k]],
            "software_contract_passed": passed,
            "scientific_qualification": passed
            and evidence.get("evidence_kind")
            == "independently_verified_QM_and_complete_model",
            "synthetic_data_cannot_qualify_physics": True,
        }
    )


def candidate(inv, targets, split, fit_contract, sources):
    parameters = [
        {
            "parameter_id": t["id"],
            "origin": "qm_fitted_extension",
            "values": None,
            "status": "unresolved",
            "family": FAMILY,
        }
        for t in inv["terms"]
    ]
    return signed(
        {
            "schema": "island.pdms-extension-candidate.v1",
            "profile": PROFILE,
            "family": FAMILY,
            "source_identities": sources,
            "chemical_scope": inv["graphs"],
            "inventory_identity": inv["identity"],
            "QM_protocol": targets["protocol"],
            "QM_target_manifest_identity": targets["identity"],
            "training_dataset_identity": digest(split["physical_samples"]["training"]),
            "validation_dataset_identity": digest(
                split["physical_samples"]["validation"]
            ),
            "dataset_status": "not_supplied; empty_role_hashes_are_not_observations",
            "parameters": parameters,
            "parameter_checksum": digest(parameters),
            "fitting_contract_identity": fit_contract["identity"],
            "fitting_solver_identity": fit_contract["solver"],
            "fitted_coefficients": 0,
            "uncertainty": "not_identifiable_no_data",
            "conditioning": "not_computable_no_design_matrix",
            "readiness_state": "derived_profile_blocked",
            **FLAGS,
        }
    )


def reports(inv, targets, bundle, *, native_immutable):
    ev = {
        "limits": deepcopy(LIMITS),
        "evidence_kind": "no_physical_dataset",
        "complete_scope_inventory": True,
        "native_PCFF_immutable": native_immutable,
        "provenance_and_hashes": True,
        "metrics": {},
        "geometry_errors": {},
    }
    gates = evaluate_gates(ev)
    blockers = [
        "execution_profile_absent",
        "completed_QM_dataset_absent",
        "finite_physical_parameter_bounds_not_supplied",
        "all_extension_parameters_unresolved",
        "reference_equilibria_and_Class_II_dependency_closure_absent",
        "charge_LJ_design_and_charge_conservation_unverified",
        "native_exclusions_and_model_policy_not_inferred",
        "no_qualified_runtime_adapter_or_registered_extension",
    ]
    return {
        "held-out-validation.json": signed(
            {
                "status": "not_evaluated_no_data",
                "training_error": None,
                "heldout_error": None,
                "physical_samples": 0,
                "reserved_validation_executed": False,
                "limits": LIMITS,
            }
        ),
        "conditioning-and-uncertainty.json": signed(
            {
                "status": "not_identifiable_no_data",
                "rank": None,
                "singular_values": None,
                "condition": None,
                "correlation": None,
                "parameter_uncertainty": None,
                "design_matrix_identity": None,
                "bootstrap_or_perturbation_runs_on_real_data": 0,
                "synthetic_controls_are_separate": True,
            }
        ),
        "qualification-report.json": signed(
            {
                "schema": "island.pdms-extension-qualification.v1",
                "profile": PROFILE,
                "candidate_identity": bundle["identity"],
                "decision": "derived_profile_blocked",
                "qualified": False,
                "planned_QM_targets": targets["planned_target_requests"],
                "completed_QM_targets": 0,
                "fitted_coefficients": 0,
                "unresolved_unique_terms": inv["unique_terms"],
                "unresolved_scalar_parameters": sum(
                    len(t["fields"]) for t in inv["terms"]
                ),
                "unresolved_policies": [p["id"] for p in inv["required_policies"]],
                "gate_evidence": ev,
                "gate_report": gates,
                **FLAGS,
            }
        ),
        "failure-report.json": signed(
            {
                "decision": "derived_profile_blocked",
                "blockers": blockers,
                "failed_acceptance_gates": gates["failed_gates"],
                "missing_term_ids": [t["id"] for t in inv["terms"]],
                "required_next_evidence": "authorized, protocol-pinned independent QM/reference campaign; explicit finite bounds; complete derived convention and dependency contract",
                "native_PCFF_conclusion_preserved": {
                    "exact_native_row_resolved": 0,
                    "unresolved_instances": 378,
                },
                **FLAGS,
            }
        ),
    }


def prepare_polymer_forcefield(
    system_payload, *, forcefield="pcff", bundle, qualification
):
    """Family-only audit selector; never silently use an incomplete extension."""
    require(forcefield == "pcff", "wrong_forcefield_family")
    scope_check(system_payload)
    require(
        bundle["family"] == FAMILY
        and bundle["profile_origin"] == "qm_fitted_extension"
        and bundle["native_pcff_complete"] is False,
        "profile_family_mixing",
    )
    require(
        bundle["parameter_checksum"] == digest(bundle["parameters"]),
        "parameter_checksum_mutation",
    )
    require(
        qualification["candidate_identity"] == bundle["identity"],
        "qualification_binding_mutation",
    )
    # An outer hash or manufactured passing gates cannot grant owned registry authority.
    allowed = bundle["identity"] in QUALIFIED_BUNDLES
    require(
        allowed and qualification.get("qualified") is True,
        "derived_profile_not_qualified",
        profile=PROFILE,
        decision="derived_profile_blocked",
        missing_terms=[
            p["parameter_id"] for p in bundle["parameters"] if p["values"] is None
        ],
    )
    raise ValueError("No J78 runtime adapter or execution profile is registered")
