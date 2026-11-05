"""Evidence gates, with no QM, evaluator, fitter, exporter or registry runner.

Authority audits are independent, hash-pinned inputs, never a boolean supplied
alongside a candidate. Production J80 has no primary artifact or authority
audit. Synthetic tests of complete contracts cannot create scientific evidence.
"""

import hashlib
import math
from copy import deepcopy
from decimal import Decimal, InvalidOperation

from island.pcff_pdms.extension import digest, require

ROW_CLASSES = (
    "exact_native_row",
    "exact_native_row_convention_blocked",
    "literature_candidate_verified",
    "literature_candidate_convention_blocked",
    "inferred_or_analog_value_blocked",
    "end_group_scope_blocked",
)
CONVENTIONS = {
    "vdw_function": (
        "source_family",
        "source_equation",
        "target_family",
        "target_equation",
    ),
    "energy_and_length_units": ("energy", "length", "angle_deviation"),
    "sigma_epsilon_minimum": (
        "sigma_meaning",
        "epsilon_meaning",
        "minimum_distance_definition",
    ),
    "combination_rules": ("distance_equation", "energy_equation"),
    "pair_exclusions_and_1_4": (
        "one_two",
        "one_three",
        "one_four_LJ",
        "one_four_Coulomb",
    ),
    "charges_and_neutrality": (
        "origin",
        "atom_type_charge_mapping",
        "type_counts_by_scope",
        "repeat_charge_e",
        "capped_chain_charge_e",
        "end_group_policy",
    ),
    "bond_and_angle_polynomials": (
        "bond_equation",
        "angle_equation",
        "coefficient_order",
        "prefactor",
        "equilibrium_dependencies",
        "angular_coordinate",
    ),
    "torsion": (
        "equation",
        "sign",
        "phase",
        "periodicities",
        "normalization",
        "angle_origin",
        "orientation",
    ),
    "improper_out_of_plane": ("equation_or_documented_nonapplicability", "ordering"),
    "cross_terms": (
        "equations",
        "ordering",
        "equilibrium_dependencies",
        "complete_required_coverage",
    ),
    "internal_and_end_group_typing": (
        "source_type_ids",
        "chemical_environments",
        "exact_graph_mapping",
        "cap_transfer_supported",
    ),
    "finite_and_periodic_scope": (
        "finite_capped",
        "periodic_support",
        "allowed_graph_identities",
        "electrostatics",
        "cutoff",
        "tail",
    ),
}
VDW_9_6 = "epsilon_ij*[2*(sigma_ij/r_ij)^9-3*(sigma_ij/r_ij)^6]"


def source_gate(manifest, contents, trusted_audits, expected_ids):
    """Validate primary bytes and audited row locations, without opening files.

    ``trusted_audits`` must come from the separately pinned intake freeze, not
    from this manifest. Hashes alone cannot supply authorship or authority.
    """
    errors = []
    artifacts = manifest.get("artifacts", [])
    rows = manifest.get("row_evidence", [])
    artifact_ids = [a.get("id") for a in artifacts]
    row_ids = [r.get("candidate_id") for r in rows]
    shas = [a.get("sha256") for a in artifacts]
    if not artifacts:
        errors.append("primary_artifact_unavailable")
    if not all(artifact_ids) or len(set(artifact_ids)) != len(artifact_ids):
        errors.append("duplicate_or_missing_source_identity")
    if len(set(shas)) != len(shas):
        errors.append("duplicate_source_bytes_identity")
    if not expected_ids or len(set(expected_ids)) != len(expected_ids):
        errors.append("expected_candidate_inventory_required")
    if len(row_ids) != len(set(row_ids)) or set(row_ids) != set(expected_ids):
        errors.append("incomplete_or_duplicate_candidate_coverage")
    accepted, assertions = {}, {}
    for artifact in artifacts:
        aid = artifact.get("id")
        sha = artifact.get("sha256")
        raw = contents.get(aid)
        if not (
            isinstance(sha, str)
            and len(sha) == 64
            and isinstance(raw, bytes)
            and hashlib.sha256(raw).hexdigest() == sha
        ):
            errors.append("source_hash_missing_or_mutated:" + str(aid))
            continue
        if artifact.get("kind") == "primary_PDF" and not raw.startswith(b"%PDF-"):
            errors.append("primary_PDF_format_mismatch:" + str(aid))
            continue
        audit = trusted_audits.get(artifact.get("authority_audit_id"))
        if not audit or audit.get("artifact_metadata") != artifact:
            errors.append("independent_primary_authority_audit_required:" + str(aid))
            continue
        if (
            artifact.get("kind")
            not in (
                "primary_PDF",
                "original_forcefield_file",
                "equivalent_primary_artifact",
            )
            or artifact.get("authority")
            not in ("original_publisher", "original_author", "audited_native_release")
            or not artifact.get("immutable_filename")
            or not artifact.get("source_version")
            or type(artifact.get("publication_year")) is not int
            or not artifact.get("license_or_access_status")
            or not all(
                artifact.get("citation", {}).get(k)
                for k in ("title", "authors", "identifier")
            )
        ):
            errors.append("incomplete_primary_provenance:" + str(aid))
            continue
        accepted[aid] = audit
        for key, value in audit.get("convention_assertions", {}).items():
            if key in assertions and assertions[key] != value:
                errors.append("conflicting_source_conventions:" + key)
            else:
                assertions[key] = value
    verified = []
    for row in rows:
        cid, aid = row.get("candidate_id"), row.get("artifact_id")
        audit = accepted.get(aid)
        if not audit:
            errors.append("candidate_primary_evidence_missing:" + str(cid))
            continue
        if (
            row.get("reporting_mode") != "directly_reported"
            or not row.get("location")
            or row.get("source_sha256") != audit["artifact_metadata"]["sha256"]
            or audit.get("verified_rows", {}).get(cid) != row
        ):
            errors.append("row_location_or_transcription_not_audited:" + str(cid))
            continue
        verified.append(cid)
    return {
        "gate": "A",
        "passed": not errors,
        "status": "passed" if not errors else "blocked",
        "blocking_reasons": sorted(set(errors)),
        "verified_candidate_ids": sorted(verified),
        "verified_rows": {
            r["candidate_id"]: deepcopy(r)
            for r in rows
            if r["candidate_id"] in verified
        },
        "verified_row_count": len(verified),
        "audited_convention_assertions": deepcopy(assertions),
        "artifact_hashes_verified": len(accepted),
    }


def known(value):
    if value is None or value == "" or value == [] or value == {}:
        return False
    if isinstance(value, str) and value in ("unknown", "unresolved", "not_established"):
        return False
    if isinstance(value, float) and not math.isfinite(value):
        return False
    if isinstance(value, dict):
        return all(known(v) for v in value.values())
    if isinstance(value, list):
        return all(known(v) for v in value)
    return True  # explicit zero exclusions or false applicability are meaningful


def convention_gate(contract, source):
    if not source["passed"]:
        return {
            "gate": "B",
            "passed": False,
            "status": "not_evaluated",
            "blocking_reasons": ["Gate_A_blocked"],
            "resolved_fields": [],
        }
    fields = contract.get("fields", {})
    errors = []
    resolved = []
    if set(fields) != set(CONVENTIONS):
        errors.append("incomplete_convention_coverage")
    for name, required in CONVENTIONS.items():
        field = fields.get(name, {})
        value = field.get("value")
        if (
            not isinstance(value, dict)
            or set(value) != set(required)
            or not known(value)
            or not field.get("evidence_refs")
        ):
            errors.append("unresolved_convention:" + name)
            continue
        if digest(field) != digest(source["audited_convention_assertions"].get(name)):
            errors.append("unaudited_or_resigned_convention:" + name)
            continue
        resolved.append(name)
    if not errors:
        units = fields["energy_and_length_units"]["value"]
        if (
            units["energy"] not in ("kJ/mol", "kcal/mol")
            or units["length"] not in ("angstrom", "nm")
            or units["angle_deviation"] not in ("radian", "degree")
        ):
            errors.append("explicit_supported_energy_length_angle_units_required")
        vdw = fields["vdw_function"]["value"]
        if (
            vdw["source_family"] != "lj_9_6"
            or vdw["target_family"] != "lj_9_6"
            or vdw["source_equation"] != VDW_9_6
            or vdw["target_equation"] != VDW_9_6
        ):
            errors.append("no_9_6_to_12_6_coercion_or_analog_form")
        if (
            fields["sigma_epsilon_minimum"]["value"]["sigma_meaning"]
            != "minimum_energy_distance"
        ):
            errors.append("9_6_minimum_distance_required")
        charges = fields["charges_and_neutrality"]["value"]
        if any(
            type(charges[k]) not in (int, float) or charges[k] != 0
            for k in ("repeat_charge_e", "capped_chain_charge_e")
        ):
            errors.append("neutral_repeat_and_capped_chain_required")
        mapping = charges["atom_type_charge_mapping"]
        counts = charges["type_counts_by_scope"]
        try:
            require(
                isinstance(mapping, dict) and isinstance(counts, dict),
                "charge_mapping_required",
            )
            require(
                set(counts) == {"repeat", "hmds", "dp2", "dp4"},
                "charge_scope_counts_required",
            )
            q = {t: Decimal(str(v)) for t, v in mapping.items()}
            require(all(v.is_finite() for v in q.values()), "finite_charges_required")
            for composition in counts.values():
                require(
                    isinstance(composition, dict) and bool(composition),
                    "type_counts_required",
                )
                require(
                    all(
                        t in q and type(n) is int and n > 0
                        for t, n in composition.items()
                    ),
                    "exact_integer_type_counts_required",
                )
                require(
                    sum(q[t] * n for t, n in composition.items()) == 0,
                    "charge_neutrality_arithmetic_failed",
                )
        except (ValueError, TypeError, InvalidOperation):
            errors.append("charge_mapping_or_neutrality_unresolved")
        if not fields["cross_terms"]["value"]["complete_required_coverage"]:
            errors.append("required_Class_II_dependencies_incomplete")
        typing = fields["internal_and_end_group_typing"]["value"]
        if not typing["cap_transfer_supported"] or any(
            t in ("Si", "O", "C", "H") for t in typing["source_type_ids"]
        ):
            errors.append("exact_source_types_and_cap_environments_required")
        if fields["finite_and_periodic_scope"]["value"]["finite_capped"] is not True:
            errors.append("finite_capped_scope_required")
    return {
        "gate": "B",
        "passed": not errors,
        "status": "passed" if not errors else "blocked",
        "blocking_reasons": sorted(set(errors)),
        "resolved_fields": resolved,
    }


def classify(row, source, convention):
    reference = source.get("verified_rows", {}).get(row.get("candidate_id"), {})
    if row.get("end_group") and not (
        source["passed"] and reference.get("audited_end_group_scope") is True
    ):
        return "end_group_scope_blocked"
    if row.get("reporting_mode") not in (
        "directly_reported",
        "unverified_user_reported",
    ) or row.get("source_family") in (
        "analog",
        "GAFF",
        "OPLS",
        "CVFF",
        "COMPASS",
        "CFF91",
    ):
        return "inferred_or_analog_value_blocked"
    verified = (
        source["passed"]
        and row.get("candidate_id") in source["verified_candidate_ids"]
        and all(
            digest(row.get(k)) == digest(reference.get(k))
            for k in ("source_family", "reporting_mode", "native_identity")
        )
    )
    native = row.get("source_family") == "native_PCFF"
    proof = row.get("native_identity", {})
    types = proof.get("ordered_source_types", [])
    complete = (
        all(
            isinstance(proof.get(k), str) and known(proof[k])
            for k in ("upstream_release", "source_row_id", "equation", "ordering")
        )
        and all(
            isinstance(proof.get(k), dict) and known(proof[k])
            for k in (
                "exact_environment",
                "units_and_phase",
                "exclusion_and_mixing_policy",
                "provenance",
            )
        )
        and isinstance(types, list)
        and bool(types)
        and all(
            isinstance(t, str) and t not in ("Si", "O", "C", "H") and "*" not in t
            for t in types
        )
    )
    exact_native = native and verified and complete
    if exact_native:
        return (
            "exact_native_row"
            if convention["passed"]
            else "exact_native_row_convention_blocked"
        )
    if native:
        return "inferred_or_analog_value_blocked"
    return (
        "literature_candidate_verified"
        if verified and convention["passed"]
        else "literature_candidate_convention_blocked"
    )


def pilot_gate(source, convention, authorization=None):
    errors = []
    if not source["passed"]:
        errors.append("Gate_A_blocked")
    if not convention["passed"]:
        errors.append("Gate_B_not_passed")
    if errors:
        return {
            "gate": "C",
            "passed": False,
            "status": "not_evaluated",
            "blocking_reasons": errors,
            "launch_allowed": False,
            "decision": "qm_pilot_blocked",
        }
    # No launcher exists here. A source/convention pass is not an execution
    # profile, a complete dataset, an identifiable fit or an acceptance grant.
    if not authorization:
        errors.append("explicit_QM_execution_profile_missing")
    return {
        "gate": "C",
        "passed": False,
        "status": "blocked",
        "blocking_reasons": errors
        or ["dataset_and_pre_fit_identifiability_not_established"],
        "launch_allowed": False,
        "decision": "qm_pilot_blocked",
    }


def validate_split(samples):
    from island.pcff_pdms.extension.fitting import split_check

    roles = {
        r: [s for s in samples if s["role"] == r]
        for r in ("training", "validation", "test")
    }
    require(sum(map(len, roles.values())) == len(samples), "invalid_dataset_role")
    split_check(roles["training"], roles["validation"], roles["test"])
    train_families = {s["geometry_family"] for s in roles["training"]}
    valid_families = {s["geometry_family"] for s in roles["validation"]}
    require(
        valid_families - train_families, "exclusive_validation_geometry_family_required"
    )
    train_topologies = {s["graph_identity"] for s in roles["training"]}
    require(
        any(s["graph_identity"] not in train_topologies for s in roles["test"]),
        "held_out_topology_required",
    )
    return True


def synthetic_diagnostics(a, parameters, bounds):
    """Manufactured matrix screening only; never fits or returns coefficients."""
    import numpy as np

    from island.pcff_pdms.extension.fitting import LIMITS, check_bounds
    from island.pcff_qm_pilot.hmds_design import metrics

    check_bounds(parameters, bounds)
    a = np.asarray(a, dtype=float)
    require(
        a.ndim == 2
        and min(a.shape) > 0
        and a.shape[1] == len(parameters)
        and np.isfinite(a).all(),
        "finite_nonempty_matrix_required",
    )
    diag = metrics(a)
    nullity = len(parameters) - diag["rank"]
    condition = diag["scaled_condition"] if nullity == 0 else None
    blockers = []
    if nullity:
        blockers.append("rank_deficient_parameter_block")
    if any(norm <= 1e-12 for norm in diag["column_norms"]):
        blockers.append("insufficient_column_norm")
    if condition is None or condition >= LIMITS["scaled_condition_max"]:
        blockers.append("conditioning_gate_failed")
    if (
        abs(diag["strongest_column_correlation"]["cosine"])
        > LIMITS["parameter_correlation_max"]
    ):
        blockers.append("parameter_correlation_gate_failed")
    return {
        "evidence_kind": "synthetic_contract_control",
        "scientific_QM_data": False,
        "design_matrix_sha256": digest(a.tolist()),
        "parameter_order": parameters,
        "bounds_identity": digest(bounds),
        "diagnostics": diag,
        "null_space_dimension": nullity,
        "condition_number": condition,
        "condition_status": "singular_infinite" if condition is None else "finite",
        "identifiable": not blockers,
        "blocking_reasons": blockers,
        "fitted_parameters": [],
        "accepted_parameters": [],
    }
