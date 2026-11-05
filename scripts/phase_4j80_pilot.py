"""J80 fenced source intake and blocked minimal QM-pilot contracts.

Gate A failure stops scientific work. The remaining outputs describe inactive
requirements, not completed source adjudication, runnable QM inputs or fits.
"""

import argparse
import json
import subprocess
from pathlib import Path

from island.pcff_pdms.pilot_intake import (
    CONVENTIONS,
    ROW_CLASSES,
    classify,
    convention_gate,
    pilot_gate,
    source_gate,
)
from scripts.phase_4j73_source_intake import (
    digest,
    exact,
    pin,
    read,
    require,
    signed,
    verify_signed,
    write,
)
from scripts.phase_4j76_source_resolution import scope

BASE = "1676a780d1f5aa1ce46e707a13aadeff14ccf39d"
BRANCH = "codex/phase-4j80-authorized-pdms-source-and-qm-pilot"
ROOT = Path("docs/evidence/phase_4j80")
FREEZE_SHA = "9b0bb8b0b68665f367d12f157a2f714fe7ce18133f05b98c43550722bcdd5385"
IMPLEMENTATION = (
    "src/island/pcff_pdms/pilot_intake/__init__.py",
    "src/island/pcff_pdms/pilot_intake/gates.py",
    "scripts/phase_4j80_pilot.py",
    "scripts/validate_phase_4j80_receipt.py",
    "tests/test_phase_4j80.py",
    "docs/phase_4j80.md",
)
INPUTS = (
    "docs/evidence/phase_4j79/receipt.json",
    "docs/evidence/phase_4j79/manifest.json",
    "docs/evidence/phase_4j79/literature-evidence.json",
    "docs/evidence/phase_4j79/j78-literature-crosswalk.json",
    "docs/evidence/phase_4j79/convention-compatibility.json",
    "docs/evidence/phase_4j79/reconciliation-report.json",
    "docs/evidence/phase_4j79/j78-input.json",
    "docs/evidence/phase_4j75/source-file-inventory.json",
    "docs/evidence/phase_4j78/qualification-report.json",
    "docs/evidence/phase_4j78/fitting-contract.json",
    "docs/evidence/phase_4j72/receipt.json",
    "docs/evidence/phase_4j72_c1/quarantine-manifest.json",
)
SCOPE = {
    "J80_QM_attempts": 0,
    "completed_QM_requests": 0,
    "J80_MD_steps": 0,
    "converter_executions": 0,
    "scientific_evaluator_executions": 0,
    "physical_fitting_executions": 0,
    "fitted_parameters": 0,
    "accepted_parameters": 0,
    "automatic_atom_type_assignment": False,
    "automatic_assignment_enabled": False,
    "runtime_registration": False,
    "production_forcefield_export": False,
    "production_MD_execution": False,
    "silent_analog_fallback": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "native_exact_rows_resolved": 0,
    "native_unresolved_instances": 378,
    "historical_reserved_validation_activated": False,
    "synthetic_controls_count_as_scientific_data": False,
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
}
BLOCKS = ("bonded", "angles", "torsions", "nonbonded", "charges", "cross_terms")


def capture(repo, root):
    scope()
    require(not root.exists(), "do_not_overwrite_capture")
    parent = repo / "docs/evidence/phase_4j79"
    evidence = read(parent / "literature-evidence.json")
    crosswalk = read(parent / "j78-literature-crosswalk.json")
    inv = read(parent / "j78-input.json")
    for value in (evidence, crosswalk, read(parent / "receipt.json")):
        verify_signed(value)
    write(
        root / "parent-inputs.json",
        {
            "literature": evidence,
            "crosswalk": crosswalk,
            "comparison": read(parent / "convention-compatibility.json"),
            "J79_result": read(parent / "reconciliation-report.json"),
            "graphs": inv["graphs"],
            "fitting_contract": read(
                repo / "docs/evidence/phase_4j78/fitting-contract.json"
            ),
        },
    )
    write(
        root / "source-catalog.json",
        read(repo / "docs/evidence/phase_4j75/source-file-inventory.json"),
    )
    write(root / "local-source-check.json", read(Path("/tmp/j80-source-recovery.json")))
    write(root / "source-access.json", read(Path("/tmp/j80-source-access.json")))
    # No synthetic audit or original asset is inserted when the primary fails.
    write(root / "trusted-authority-audits.json", {})
    history = {}
    tree = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    for row in tree.splitlines():
        meta, name = row.split("\t", 1)
        p = repo / name
        history[name] = (
            digest(str(p.readlink())) if meta.startswith("120000") else pin(p)
        )
    write(root / "historical-files.json", history)
    files = (
        "parent-inputs.json",
        "source-catalog.json",
        "local-source-check.json",
        "source-access.json",
        "trusted-authority-audits.json",
        "historical-files.json",
    )
    write(
        root / "frozen-inputs.json",
        signed(
            {
                "schema": "island.j80-input-freeze.v1",
                "base": BASE,
                "fetched_origin_main": subprocess.check_output(
                    ["git", "rev-parse", "origin/main"], cwd=repo, text=True
                ).strip(),
                "scope": SCOPE,
                "original_inputs": {n: pin(repo / n) for n in INPUTS},
                "capture_files": {n: pin(root / n) for n in files},
            }
        ),
    )
    return {
        "freeze_sha256": pin(root / "frozen-inputs.json"),
        "historical_files": len(history),
    }


def inputs(root):
    scope()
    require(
        pin(root / "frozen-inputs.json") == FREEZE_SHA,
        "fabricated_or_resigned_input_freeze",
    )
    frozen = read(root / "frozen-inputs.json")
    verify_signed(frozen)
    exact(frozen["scope"], SCOPE, "scientific_scope_mutation")
    require(frozen["base"] == BASE, "base_mutation")
    for n, sha in frozen["capture_files"].items():
        require(pin(root / n) == sha, "captured_source_hash_mutation:" + n)
    return read(root / "parent-inputs.json")


def reconstruct(root):
    parent = inputs(root)
    literature = parent["literature"]
    reported = literature["reported_rows"]
    ids = [r["id"] for r in reported]
    local = read(root / "local-source-check.json")
    access = read(root / "source-access.json")
    require(
        not access["full_text_recovered"],
        "new_primary_requires_independent_source_intake_not_historical_claim_promotion",
    )
    require(
        not local["new_or_modified_files"] and not local["missing_files"],
        "local_catalog_changed_requires_source_review",
    )
    catalog = read(root / "source-catalog.json")
    provenance = signed(
        {
            "schema": "island.j80-source-provenance.v1",
            "citation": literature["bibliography"],
            "primary_source_status": "unavailable",
            "artifacts": [],
            "required_artifact_fields": [
                "immutable_filename",
                "sha256",
                "kind",
                "authority",
                "authority_audit_id",
                "citation",
                "source_version",
                "publication_year",
                "license_or_access_status",
            ],
            "required_row_fields": [
                "candidate_id",
                "source_sha256",
                "artifact_id",
                "location",
                "reporting_mode",
                "source_family",
                "raw_values_and_units",
            ],
            "row_evidence": [
                {
                    "candidate_id": r["id"],
                    "J79_record_identity": r["identity"],
                    "artifact_id": None,
                    "source_sha256": None,
                    "source_version": None,
                    "publication_year": 1997,
                    "citation": literature["bibliography"],
                    "location": None,
                    "reporting_mode": "unverified_user_reported",
                    "source_family": "published_PDMS_PCFF_candidate_not_native_row_authority",
                    "raw_values_and_units": {
                        "reported_values": r["raw_reported_values"],
                        "reported_units": r["reported_units"],
                        "verified_source_units": None,
                    },
                }
                for r in reported
            ],
            "context_files": [
                {
                    "relative_path": r["relative_path"],
                    "sha256": r["sha256"],
                    "bytes": r["bytes"],
                    "citation": r["citation"],
                    "source_family": r["source_family"],
                    "source_version": r["version"],
                    "license_status": r["license_status"],
                    "role": "context_only_not_authoritative_for_the_20_reported_candidates",
                }
                for r in catalog["files"]
            ],
            "current_local_hash_check": local,
            "publisher_access": access,
            "metadata_is_not_parameter_evidence": True,
            "authority_audits_available": 0,
        }
    )
    trusted = read(root / "trusted-authority-audits.json")
    a = source_gate(provenance, {}, trusted, ids)
    require(not a["passed"], "unexpected_source_qualification")
    contract = signed(
        {
            "schema": "island.j80-convention-contract.v1",
            "fields": {
                k: {
                    "value": None,
                    "evidence_refs": [],
                    "required_value_fields": list(v),
                    "status": "not_evaluated_Gate_A_blocked",
                }
                for k, v in CONVENTIONS.items()
            },
            "J79_comparison_reference": {
                "generic_12_6_incompatible": True,
                "isolated_9_6_exists_but_not_called": True,
                "paper_mixing_rule": None,
                "capped_charge_if_repeat_claims_transferred_e": "+0.240",
                "source_verified": False,
            },
            "no_cross_family_or_9_6_to_12_6_substitution": True,
            "internal_types_are_not_element_sequence_aliases": True,
            "runtime_use_policy": "only independently exact_native_row may be considered in a future authorized phase; J80 registers nothing",
        }
    )
    b = convention_gate(contract, a)
    c = pilot_gate(a, b)
    candidates = [
        {
            **r,
            "classification": classify(r, a, b),
            "source_verified": False,
            "future_runtime_consideration": False,
        }
        for r in provenance["row_evidence"]
    ]
    target_rows = []
    for old in parent["crosswalk"]["rows"]:
        temp = {
            "candidate_id": old["J78_parameter_id"],
            "source_family": "unqualified_extension_target",
            "reporting_mode": "unverified_user_reported"
            if old["reported_candidate_ids"]
            else "missing_value",
            "end_group": old["end_group_involved"],
        }
        target_rows.append(
            {
                "J78_parameter_id": old["J78_parameter_id"],
                "J78_parameter_identity": old["J78_parameter_identity"],
                "J79_row_identity": old["identity"],
                "J79_status": old["status"],
                "ordered_types": old["ordered_types"],
                "family": old["family"],
                "reported_candidate_ids": old["reported_candidate_ids"],
                "cross_term": old["cross_term"],
                "end_group": old["end_group_involved"],
                "classification": classify(temp, a, b),
                "resolved": False,
                "classification_scope": "missing target value; no analog or inferred coefficient supplied",
            }
        )
    row_report = signed(
        {
            "schema": "island.j80-row-classification.v1",
            "allowed_classes": list(ROW_CLASSES),
            "candidate_rows": candidates,
            "candidate_class_counts": {
                s: sum(r["classification"] == s for r in candidates)
                for s in ROW_CLASSES
            },
            "extension_targets": target_rows,
            "target_class_counts": {
                s: sum(r["classification"] == s for r in target_rows)
                for s in ROW_CLASSES
            },
            "verified_literature_rows": 0,
            "blocked_candidate_rows": len(candidates),
            "unresolved_cross_terms": sum(r["cross_term"] for r in target_rows),
            "unresolved_end_group_terms": sum(r["end_group"] for r in target_rows),
            "unresolved_extension_terms": len(target_rows),
        }
    )
    require(
        len(candidates) == 20 and len(target_rows) == 129, "incomplete_row_coverage"
    )
    graphs = [g for g in parent["graphs"] if g["fixture"] in ("hmds", "dp2", "dp4")]
    request = signed(
        {
            "schema": "island.j80-QM-pilot-requests.v1",
            "status": "not_generated_Gate_A_blocked",
            "execution_enabled": False,
            "execution_profile": None,
            "requests": [],
            "completed_requests": [],
            "QM_input_decks_generated": 0,
            "prospective_contract_only": {
                "max_requests": 36,
                "full_423_campaign_allowed": False,
                "molecules": [
                    {
                        "name": g["fixture"],
                        "graph_identity": g["graph_identity"],
                        "explicit_atoms": g["sites"],
                        "support_reference": "immutable J78 graph/scope validation, not QM data",
                    }
                    for g in graphs
                ],
                "charge": 0,
                "multiplicity": 1,
                "geometry_scope": "finite linear trimethylsilyl-capped explicit-H; no periodic calculation",
                "training_topologies": ["hmds", "dp2"],
                "exclusive_validation_family": "independent_coupled_deformations_never_used_in_training",
                "held_out_topology": "dp4",
                "geometry_families_required": [
                    "equilibrium",
                    "Si-O_Si-C_C-H_bond_distortions",
                    "O-Si-O_Si-O-Si_C-Si-O_H-C-Si_angle_distortions",
                    "representative_torsion_scans",
                    "separate_intermolecular_nonbonded_cases",
                    "separate_electrostatic_validation",
                    "coupled_cross_term_deformations",
                ],
                "protocol_required_before_requests": {
                    k: None
                    for k in (
                        "software",
                        "version",
                        "method",
                        "basis",
                        "dispersion",
                        "convergence",
                        "constraints",
                        "independent_coordinates",
                        "atom_mapping",
                        "input_hash",
                        "output_hash",
                    )
                },
                "data_leakage_gate": "disjoint sample IDs, geometry hashes and conformer groups; exclusive validation family and held-out graph",
                "historical_reserved_validation_activated": False,
            },
            "parameter_blocks": list(BLOCKS),
            "coupled_or_joint_fit_allowed": False,
            "fitting_enabled": False,
        }
    )
    conditioning = signed(
        {
            "schema": "island.j80-identifiability.v1",
            "computed": False,
            "reason": "Gate A blocked; no independent QM observations or supplied design matrix; no physical fitting",
            "blocks": [
                {
                    "block": name,
                    "status": "unresolved_not_computed",
                    "rank": None,
                    "condition_number": None,
                    "null_space_dimension": None,
                    "parameter_correlation": None,
                    "design_matrix_sha256": None,
                    "bounds_identity": None,
                    "uncertainty": None,
                    "identifiable": False,
                    "fitted_parameters": 0,
                }
                for name in BLOCKS
            ],
            "reuse": [
                "island.pcff_pdms.extension.fitting.check_bounds",
                "island.pcff_pdms.extension.fitting.split_check",
                "island.pcff_qm_pilot.hmds_design.metrics",
            ],
            "pre_fit_contract": {
                "parameter_order_fixed": True,
                "explicit_finite_physical_bounds_required": True,
                "rank_nullspace_condition_correlation_before_fit": True,
                "rank_deficient_or_ill_conditioned_blocks_remain_unresolved": True,
                "source_of_thresholds": "immutable J78 fitting-contract.json; no widened thresholds",
                "J78_limits": parent["fitting_contract"]["limits"],
                "synthetic_diagnostics_are_scientific_evidence": False,
            },
        }
    )
    return {
        "source-provenance-manifest.json": provenance,
        "gate-report.json": signed(
            {"A": a, "B": b, "C": c, "scientific_work_stopped_at": "Gate_A"}
        ),
        "convention-contract.json": contract,
        "candidate-row-classification.json": row_report,
        "qm-pilot-request-manifest.json": request,
        "identifiability-and-conditioning.json": conditioning,
        "result.json": signed(
            {
                "decision": "qm_pilot_blocked",
                "source_recovery_status": "primary_source_unavailable",
                "convention_status": "not_evaluated_Gate_A_blocked",
                "verified_literature_rows": 0,
                "blocked_rows": 20,
                "unresolved_cross_terms": 87,
                "unresolved_end_group_terms": 69,
                "unresolved_extension_terms": 129,
                "unresolved_scalar_parameters": 365,
                "unresolved_policies": 11,
                "completed_QM_requests": 0,
                "fitted_parameters": 0,
                "accepted_parameters": 0,
                "pilot_requests_generated": 0,
                "scientific_work_stopped_at": "Gate_A",
                "scope": SCOPE,
            }
        ),
    }


def replay(root, build=False):
    documents = reconstruct(root)
    for n, value in documents.items():
        if build:
            require(not (root / n).exists(), "do_not_overwrite_attempt")
            write(root / n, value)
        exact(read(root / n), value, "resigned_gate_semantic_mutation:" + n)
    return documents["result.json"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("capture", "build", "replay"))
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    a = p.parse_args()
    try:
        scope()
        value = (
            capture(a.repo, a.directory)
            if a.mode == "capture"
            else replay(a.directory, a.mode == "build")
        )
        print(json.dumps(value, sort_keys=True))
        return 0
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as e:
        print(json.dumps({"status": "blocked", "reason": str(e), "scope": SCOPE}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
