"""Deterministic evidence normalization, with independent source/convention gates."""

import tarfile
import zipfile
from copy import deepcopy
from decimal import Decimal

from island.pcff_pdms.extension import digest, require
from island.pcff_pdms.pilot_intake import gates

from . import parsing
from .schema import PARSER_VERSION, STATUSES, VERSION, validate

CROSS_TERMS = (
    "bond-bond",
    "bond-bond_1_3",
    "bond-angle",
    "angle-angle",
    "end-bond-torsion",
    "middle-bond-torsion",
    "angle-torsion",
    "angle-angle-torsion",
    "torsion-torsion",
)
SCOPE = {
    "J81_QM_attempts": 0,
    "completed_QM_requests": 0,
    "J81_MD_steps": 0,
    "converter_executions": 0,
    "fitting_executions": 0,
    "evaluator_executions": 0,
    "reserved_validation_executions": 0,
    "parameter_generation_executions": 0,
    "fitted_parameters": 0,
    "accepted_parameters": 0,
    "automatic_atom_type_assignment": False,
    "automatic_assignment": False,
    "runtime_registration": False,
    "production_parameter_export": False,
    "physical_MD_execution": False,
    "silent_analog_fallback": False,
    "synthetic_fixtures_are_scientific_data": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "native_PCFF_complete": False,
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
}


def capture(manifest, root):
    """Read only explicitly supplied files after installing the inherited fence."""
    from scripts.phase_4j76_source_resolution import scope

    scope()
    validate(manifest)
    require(
        len(manifest["artifacts"]) <= 256 and len(manifest["candidates"]) <= 4096,
        "bundle_size_limit",
    )
    ids = [a["id"] for a in manifest["artifacts"]]
    shas = [a["sha256"] for a in manifest["artifacts"]]
    candidates = manifest["candidates"]
    require(
        len(ids) == len(set(ids)) and len(shas) == len(set(shas)),
        "duplicate_source_identity",
    )
    cids = [c["id"] for c in candidates]
    require(len(cids) == len(set(cids)), "duplicate_candidate_identity")
    require(
        all(c["source_id"] in ids for c in candidates),
        "candidate_source_identity_missing",
    )
    locations = [(c["source_id"], digest(c["source_location"])) for c in candidates]
    require(
        len(locations) == len(set(locations)), "duplicate_candidate_location_conflict"
    )
    artifacts = []
    for a in manifest["artifacts"]:
        item = {
            "metadata": deepcopy(a),
            "status": "source_present_unparsed",
            "errors": [],
            "extraction": None,
        }
        try:
            raw = parsing.local_bytes(root, a["source_path"])
            require(
                parsing.sha(raw) == a["sha256"] and len(raw) == a["file_size"],
                "source_checksum_or_size_mutation",
            )
            selected = [c for c in candidates if c["source_id"] == a["id"]]
            item["extraction"] = parsing.extract(raw, a, selected)
            item["status"] = (
                "source_parsed"
                if item["extraction"]["documents"]
                and not item["extraction"]["failures"]
                else "source_present_unparsed"
            )
        except (
            ValueError,
            OSError,
            KeyError,
            TypeError,
            tarfile.TarError,
            zipfile.BadZipFile,
        ) as e:
            item["errors"] = [str(e)]
            if getattr(e, "reason", str(e)) == "source_missing":
                item["status"] = "source_missing"
        artifacts.append(item)
    return {
        "schema": VERSION,
        "parser_version": PARSER_VERSION,
        "manifest": deepcopy(manifest),
        "artifacts": artifacts,
    }


def bridge(snapshot):
    """Bind all interpreted metadata and exact excerpts into the J80 audit format."""
    artifacts, rows = [], []
    for item in snapshot["artifacts"]:
        a = item["metadata"]
        kind = (
            "primary_PDF"
            if a["source_type"] == "primary_pdf"
            else "original_forcefield_file"
            if a["source_type"] == "original_pcff_frc"
            else "equivalent_primary_artifact"
        )
        artifacts.append(
            {
                "id": a["id"],
                "sha256": a["sha256"],
                "kind": kind,
                "authority": a["authority"],
                "authority_audit_id": a["id"],
                "immutable_filename": a["source_path"],
                "source_version": a["source_version"],
                "publication_year": a["bibliography"]["year"],
                "citation": {
                    "title": a["bibliography"]["title"],
                    "authors": a["bibliography"]["authors"],
                    "identifier": a["bibliography"]["identifier"],
                },
                "license_or_access_status": a["license_or_redistribution_status"],
                "J81_metadata": deepcopy(a),
            }
        )
    for c in snapshot["manifest"]["candidates"]:
        item = next(
            i for i in snapshot["artifacts"] if i["metadata"]["id"] == c["source_id"]
        )
        fact = (item["extraction"] or {}).get("facts", {}).get(c["id"])
        if not fact:
            continue
        rows.append(
            {
                "candidate_id": c["id"],
                "artifact_id": c["source_id"],
                "source_sha256": item["metadata"]["sha256"],
                "location": c["source_location"],
                "reporting_mode": c["reporting_mode"],
                "source_family": item["metadata"]["source_family"],
                "raw_values_and_units": deepcopy(c),
                "native_identity": c["native_identity"],
                "audited_end_group_scope": c["chain_scope"]
                in ("both", "capped_end_group"),
            }
        )
    return {"artifacts": artifacts, "row_evidence": rows}


def field_checks(row, contract):
    """Compare row claims with the actual contract; never infer missing conventions."""
    errors, fields = [], contract.get("fields", {})

    def value(name):
        return fields.get(name, {}).get("value") or {}

    for k in (
        "units",
        "functional_form",
        "mixing_rule",
        "exclusion_and_1_4_policy",
        "charge_convention",
        "torsion_convention",
        "chemical_environment",
    ):
        if not gates.known(row[k]):
            errors.append("unresolved_row_field:" + k)
    expected = {
        "units": value("energy_and_length_units"),
        "mixing_rule": value("combination_rules"),
        "exclusion_and_1_4_policy": value("pair_exclusions_and_1_4"),
        "charge_convention": value("charges_and_neutrality"),
        "torsion_convention": value("torsion"),
    }
    units = row["units"]
    expected_units = expected.pop("units")
    if not isinstance(units, dict) or any(
        units.get(k) != v for k, v in expected_units.items()
    ):
        errors.append("row_convention_mismatch:units")
    per_value = units.get("per_value") if isinstance(units, dict) else None
    if (
        not isinstance(per_value, list)
        or len(per_value) != len(row["numerical_values"])
        or not all(isinstance(v, str) and gates.known(v) for v in per_value)
    ):
        errors.append("explicit_per_value_units_required")
    for k, v in expected.items():
        if not v or row[k] != v:
            errors.append("row_convention_mismatch:" + k)
    family = row["parameter_family"]
    equations = {
        "bond": value("bond_and_angle_polynomials").get("bond_equation"),
        "angle": value("bond_and_angle_polynomials").get("angle_equation"),
        "torsion": value("torsion").get("equation"),
        "nonbonded": value("vdw_function").get("source_equation"),
        "improper": value("improper_out_of_plane").get(
            "equation_or_documented_nonapplicability"
        ),
        "charge": "source_type_charge",
    }
    cross = value("cross_terms").get("equations")
    if isinstance(cross, dict):
        equations.update(cross)
    if not equations.get(family) or row["functional_form"] != equations.get(family):
        errors.append("row_equation_mismatch:" + family)
    types = row["atom_type_sequence"]
    source_types = value("internal_and_end_group_typing").get("source_type_ids", [])
    if not types or any(
        t in ("Si", "C", "H", "O", "X") or "*" in t or t not in source_types
        for t in types
    ):
        errors.append("exact_source_atom_types_unresolved")
    if row["native_identity"]:
        proof = row["native_identity"]
        if (
            proof.get("ordered_source_types") != types
            or proof.get("equation") != row["functional_form"]
            or proof.get("exact_environment") != row["chemical_environment"]
        ):
            errors.append("native_identity_does_not_match_actual_row")
    return sorted(errors)


def source_evidence(snapshot, audits):
    """Audit assertions are independently hash-pinned, never manifest attestations.

    Replay authenticates the captured extraction instead of pretending to rehash
    an unavailable original. Live source byte checks remain a separate operation.
    """
    mapped = bridge(snapshot)
    errors, verified, assertions = [], [], {}
    if not mapped["artifacts"]:
        errors.append("primary_artifact_unavailable")
    if not snapshot["manifest"]["candidates"]:
        errors.append("candidate_inventory_missing")
    for x, item in zip(mapped["artifacts"], snapshot["artifacts"], strict=True):
        audit = audits.get(x["id"], {})
        if (
            audit.get("capture_identity") != digest(snapshot)
            or audit.get("artifact_metadata") != x
            or audit.get("evidence_kind") != "independently_reviewed_primary_artifact"
        ):
            errors.append("independent_primary_authority_audit_required:" + x["id"])
            continue
        if item["status"] != "source_parsed" or item["errors"]:
            errors.append("source_not_successfully_parsed:" + x["id"])
        if x["authority"] not in (
            "original_publisher",
            "original_author",
            "audited_native_release",
        ):
            errors.append("independent_primary_authority_required")
        if x["J81_metadata"]["authorization_status"] not in (
            "authorized_local_use",
            "public_access",
        ):
            errors.append("source_use_authorization_missing")
        if not x["citation"]["identifier"]:
            errors.append("persistent_source_identity_required_for_verification")
        for k, v in audit.get("convention_assertions", {}).items():
            if k in assertions and assertions[k] != v:
                errors.append("conflicting_source_conventions:" + k)
            assertions[k] = deepcopy(v)
        for r in (r for r in mapped["row_evidence"] if r["artifact_id"] == x["id"]):
            if (
                r["reporting_mode"] != "directly_reported"
                or audit.get("verified_rows", {}).get(r["candidate_id"]) != r
            ):
                errors.append(
                    "row_location_or_transcription_not_audited:" + r["candidate_id"]
                )
            else:
                verified.append(r["candidate_id"])
    if len(mapped["row_evidence"]) != len(snapshot["manifest"]["candidates"]):
        errors.append("incomplete_parsed_candidate_coverage")
    return {
        "gate": "A",
        "passed": not errors,
        "status": "passed" if not errors else "blocked",
        "blocking_reasons": sorted(set(errors)),
        "verified_candidate_ids": verified,
        "verified_rows": {
            r["candidate_id"]: r
            for r in mapped["row_evidence"]
            if r["candidate_id"] in verified
        },
        "verified_row_count": len(verified),
        "audited_convention_assertions": assertions,
        "artifact_hashes_verified": sum(
            i["status"] == "source_parsed" for i in snapshot["artifacts"]
        ),
    }


def normalize(snapshot, audits, expected_ids=()):
    """No files or scientific backends are accessed by this deterministic step."""
    manifest = snapshot["manifest"]
    validate(manifest)
    require(
        isinstance(audits, dict) and all(isinstance(v, dict) for v in audits.values()),
        "independent_audit_schema_required",
    )
    mapped = bridge(snapshot)
    a = source_evidence(snapshot, audits)
    contract = manifest["convention_contract"]
    b = gates.convention_gate(contract, a)
    normalized = []
    for c in manifest["candidates"]:
        item = next(
            i for i in snapshot["artifacts"] if i["metadata"]["id"] == c["source_id"]
        )
        fact = (item["extraction"] or {}).get("facts", {}).get(c["id"])
        problems = (
            field_checks(c, contract)
            if a["passed"]
            else ["source_authority_not_verified; conventions_not_evaluated"]
        )
        status = "source_parsed" if fact else item["status"]
        if fact and not a["passed"]:
            status = "literature_candidate_convention_blocked"
        if fact and a["passed"]:
            row = next(
                r for r in mapped["row_evidence"] if r["candidate_id"] == c["id"]
            )
            row["end_group"] = c["chain_scope"] in ("both", "capped_end_group")
            checked = deepcopy(b)
            checked["passed"] = b["passed"] and not problems
            status = gates.classify(row, a, checked)
            if status == "exact_native_row_convention_blocked":
                status = "literature_candidate_convention_blocked"
        if c["reporting_mode"] in ("reconstructed", "inferred") or item["metadata"][
            "source_family"
        ] in ("OPLS", "GAFF", "CVFF", "COMPASS", "CFF91", "analog"):
            status = "inferred_or_analog_value_blocked"
        if c["chain_scope"] in ("capped_end_group", "both") and not (
            a["passed"] and b["passed"]
        ):
            status = "end_group_scope_blocked"
        require(status in STATUSES, "unknown_row_status")
        normalized.append(
            {
                **deepcopy(c),
                "source_sha256": item["metadata"]["sha256"],
                "source_family": item["metadata"]["source_family"],
                "extracted_fact": fact,
                "declared_field_diagnostics": field_checks(c, contract) if fact else [],
                "declared_diagnostics_are_authority": False,
                "raw_values": list(c["numerical_values"]),
                "interpreted_values": [str(Decimal(v)) for v in c["numerical_values"]]
                if fact
                else None,
                "original_units": deepcopy(c["units"]),
                "converted_values": None,
                "unit_conversion": "not performed; source units preserved; no 9-6 to 12-6 conversion",
                "interpretation_status": "decimal_transcribed_no_parameter_acceptance"
                if fact
                else "unparsed",
                "verification_status": status,
                "blocking_reasons": sorted(
                    set(
                        problems
                        + item["errors"]
                        + list(
                            (item["extraction"] or {})
                            .get("failures", {})
                            .get(c["id"], "")
                            .splitlines()
                        )
                    )
                ),
            }
        )
    verified = [
        r
        for r in normalized
        if r["verification_status"]
        in ("literature_candidate_verified", "exact_native_row")
    ]
    ready = (
        a["passed"]
        and b["passed"]
        and len(verified) == len(normalized)
        and bool(expected_ids)
        and set(expected_ids) <= {r["id"] for r in verified}
    )
    parsed = any(
        i["extraction"] and i["extraction"]["facts"] for i in snapshot["artifacts"]
    )
    decision = (
        "source_parsed_ready_for_qm_authorization"
        if ready
        else "source_parsed_convention_blocked"
        if parsed
        else "source_ingestion_blocked"
    )
    return {
        "decision": decision,
        "source_gate": a,
        "convention_comparison": b,
        "candidate_rows": normalized,
        "candidate_status_counts": {
            s: sum(r["verification_status"] == s for r in normalized) for s in STATUSES
        },
        "verified_literature_rows": len(verified),
        "native_exact_rows": sum(
            r["verification_status"] == "exact_native_row" for r in normalized
        ),
        "blocked_supplied_rows": len(normalized) - len(verified),
        "blocked_rows": len(set(expected_ids) - {r["id"] for r in verified})
        + sum(r["id"] not in expected_ids and r not in verified for r in normalized),
        "source_ingestion_status": "source_missing"
        if not snapshot["artifacts"]
        else "source_parsed"
        if parsed
        else "source_present_unparsed",
        "historical_candidates_still_blocked": sorted(
            set(expected_ids) - {r["id"] for r in verified}
        ),
        "missing_cross_term_classes": [
            f
            for f in CROSS_TERMS
            if not any(r["parameter_family"] == f for r in verified)
        ],
        "J80_preserved": {
            "decision": "qm_pilot_blocked",
            "verified_literature_rows": 0,
            "blocked_candidate_rows": 20,
            "unresolved_cross_terms": 87,
            "unresolved_end_group_terms": 69,
            "unresolved_extension_terms": 129,
        },
        "unresolved_cross_terms": 87,
        "unresolved_end_group_terms": 69,
        "scope": deepcopy(SCOPE),
        "qm_authorization_granted": False,
        "qualification_note": "Parsing/readiness for a future authorization review never resolves J78 targets or enables execution.",
    }
