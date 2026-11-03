"""Offline, additive J72 source adjudication. No QM, fitting, MD or conversion."""

import argparse
from collections import Counter
from itertools import combinations
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.literature_intake import candidate_sites
from island.pcff_qm_pilot.source_adjudication import (
    AUTHORITY,
    CATEGORIES,
    compare,
    execution_deviation,
    normalized_cell,
    verify_signed,
)
from scripts.replay_phase_4j71 import replay as replay_parent

ROOT = Path("docs/evidence/phase_4j72")
REQUIRED = BASE = "e0c82269164781e0aa12881ebb6b308e21735e84"
BRANCH = "codex/phase-4j72-pcff-source-conflict-adjudication"
FETCHED_MAIN = "4c9066e4501591211049c44401e82b3110118f58"
FROZEN_DECLARATIONS = {
    "decision": "source_conflict_unresolved",
    "parameter_completeness": "incomplete",
    "accepted_coefficients": 0,
    "automatic_assignment": False,
    "new_QM": 0,
    "QM_fits": 0,
    "MD_steps": 0,
    "reserved_validation_activated": False,
    "reserved_validation_consumed": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "scientific_validation": "not_established",
    "chemical_transfer": "not_established",
    "original_author_artifact_status": "not_yet_supplied",
    "author_confirmations": "28 user-attested",
}
DECLARATIONS = {
    **FROZEN_DECLARATIONS,
    "new_QM": 4,
    "HMDS_new_QM": 0,
    "MD_steps": None,
    "MD_steps_verified_minimum": 46,
    "MD_steps_possible_maximum": 48,
    "authorized_new_QM": 0,
    "authorized_MD_steps": 0,
    "zero_execution_scope_met": False,
}
PARENT = Path("docs/evidence/phase_4j71")
# Frozen review contract, pinned after writing and before adjudication execution.
CONTRACT_SHA = "f2c52ba6eb3f0bfae459bb538170cc0019ddfc713fe2f3990d5b770e6bc6b28e"


def inputs(repo, root):
    require(pin(root / "contract.json") == CONTRACT_SHA, "J72_frozen_contract_mutation")
    plan = read(root / "contract.json")
    verify_signed(plan)
    require(
        plan["base_commit"] == BASE
        and plan["required_J71"] == REQUIRED
        and plan["branch"] == BRANCH
        and plan["fetched_origin_main"] == FETCHED_MAIN,
        "J72_ancestry",
    )
    require(plan["scientific_scope"] == FROZEN_DECLARATIONS, "J72_contract_scope")
    for name, expected in plan["parent_pins"].items():
        require(pin(repo / name) == expected, "J72_parent_mutation", path=name)
    replay_parent(repo, repo / PARENT)
    records = read(repo / PARENT / "parameter-records.json")["records"]
    cells = read(repo / PARENT / "extraction-ledger.json")["cells"]
    reading = read(repo / PARENT / "reading-status.json")
    sources = {r["sha256"]: r for r in reading["files"]}
    for r in reading["files"]:
        if r["scope"] == "primary_directory":
            sources[r["sha256"]] = r
    matrix = read(repo / "docs/evidence/phase_4j70/completeness-matrix.json")
    require(
        len(records) == 467 and len(matrix["instances"]) == 378, "J72_parent_inventory"
    )
    return plan, records, cells, sources, matrix


def source_context(record):
    """Bounded scalar-comparison scope, not a source-to-native compatibility claim."""
    src, model = record["source_file_sha256"][:12], record["source_force_field_family"]
    if model == "PolCA":
        return {
            "contract_scope": "PolCA2021 final model; main and associated SI",
            "charge_nonbonded": "PolCA2021 Eqs12-15 LJ12-6/Lorentz-Berthelot, final fixed charges; 1-4 excluded (PDF p5)",
            "source_target": "PolCA2021 final organosilicon model; molecular torsion profiles and liquid properties, not HMDS J68 residual selection",
            "location": "main PDF p5 Eqs9-15; p8 Table5; SI p65 TableS72",
        }
    if src in ("6e86aa640854", "8b0e8b74303b") and model == "COMPASS":
        return {
            "contract_scope": "benchmark2025 nominal COMPASS labels versus topology-bound associated archive",
            "charge_nonbonded": "benchmark COMPASS fixed charges and LJ9-6; native PCFF contract not established",
            "source_target": "same associated benchmark COMPASS representation, source version/correction unavailable",
            "location": "benchmark SI pp4-5 TableS3; archive param/COMPASS/system.in.settings and topology",
        }
    return {
        "contract_scope": f"{src}:{model}:source-specific environment; no cross-source equivalence declared",
        "charge_nonbonded": None,
        "source_target": None,
        "location": "original J71 source location and complete reading notes; unbound assumptions remain unknown",
    }


def comparison_record(record, cell, source):
    context = source_context(record)
    family = record["interaction_family"]
    elements = candidate_sites(cell)
    ordered = list(record["ordered_sites"])
    # Only symmetric scalar bond/angle values are canonicalized. Cross terms
    # retain their coefficient/reference order; no guessed reversal is applied.
    canonical_order = ordered
    if (
        context["contract_scope"].startswith("benchmark2025")
        and family in ("quartic_bond", "quartic_angle")
        and None not in elements
    ):
        canonical_order = min(elements, elements[::-1])
    energy = record["scalar_energy_equation"]
    issues = list(record["missing_convention_fields"])
    if record["extraction_confidence"] == "extraction_unresolved":
        issues.append("unresolved_extraction")
    # The source's UA hybrid-style spelling is inconsistent with a blanket
    # LJ12-6 interpretation. Keep J71 unchanged; block scalar comparison here.
    if (
        record["source_force_field_family"] == "UA"
        and family == "LJ12_6"
        and "lj/class2" in record["location"]["raw_snippet"]
    ):
        energy = None
        issues.append("UA_printed_lj_class2_style_vs_LJ12_6_contract_unresolved")
    dependencies = (
        None
        if "exact_source_equilibrium_dependencies" in issues
        else record["cross_term_dependencies"]
    )
    prefactor = {
        "harmonic_angle": "1/2",
        "opls_torsion": "1/2",
        "smith_torsion": "-1/2 with source sine-squared modulation",
    }.get(family, "as printed in scalar equation; no additional prefactor")
    contract = {
        "interaction_family": family,
        "ordered_sites": canonical_order,
        "chemical_environment": context["contract_scope"],
        "resolution": record["source_model_resolution"],
        "forcefield_family": record["source_force_field_family"],
        "energy_equation": energy,
        "force_equation": record["analytic_force_equation"],
        "phase": None
        if any("phase" in i for i in issues)
        else record["phase_convention"],
        "prefactor": prefactor,
        "multiplicity": record["multiplicity"],
        "orientation": record["orientation_reversal_rule"],
        "equilibrium_dependencies": dependencies,
        "charge_nonbonded_contract": context["charge_nonbonded"],
        "source_target": context["source_target"],
    }
    return signed(
        {
            "schema": "island.j72-comparison-record.v1",
            "record_id": record["record_id"],
            "J71_identity": record["identity"],
            "source_file_sha256": record["source_file_sha256"],
            "citation": record["citation"],
            "location": record["location"],
            "raw_value": record["raw_value"],
            "raw_unit": record["raw_unit"],
            "raw_ordered_sites": ordered,
            "raw_energy_equation": record["scalar_energy_equation"],
            "fields": record["fields"],
            "contract": contract,
            "equilibrium_values": record["equilibrium_values"],
            "comparison_values": [
                normalized_cell(v, u)
                for v, u in zip(record["raw_value"], record["raw_unit"], strict=True)
            ],
            "context_evidence_location": context["location"],
            "source_chemical_domain": source.get(
                "chemical_domain", "not_established_in_archive_metadata"
            ),
            "source_fitting_context": {
                "9c6d6315b91e": {
                    "location": "PDF p3 SectionIII; p5 preceding SectionIV; Appendix p11",
                    "context": "QM energies, first/second derivatives and ESP, generic scaling and refined nonbonded terms; methodology described as CFF93-type; resulting Appendix explicitly named PCFF. Polysilane environments are not automatically oxygenated sio types.",
                },
                "65bda3a5b14c": {
                    "location": "PDF p1 abstract; p3 Section2 and Eq1b; Appendix pp21-22",
                    "context": "CFF-type siloxane valence fitted to ab initio energies/derivatives; van der Waals fitted to liquid density/cohesive energy density at atmospheric pressure and one temperature. No single row-specific temperature or omitted coupling value inferred.",
                },
                "a1415925e044": {
                    "location": "PDF p3 SectionIIIA-B, Eqs1-5 and Table3; p5 Table4",
                    "context": "ESP-fitted molecular charges; silicon exp6 from SiH4 pair binding; equilibrium geometries and DMTS conformational energetics. Half-scaled 1-4 nonbonded terms and source molecule-dependent charges differ from PCFF.",
                },
            }.get(
                record["source_file_sha256"][:12],
                {
                    "location": context["location"],
                    "context": context["source_target"],
                },
            ),
            "source_temperature": None,
            "temperature_status": "no single parameter-row fitting temperature established; no temperature equivalence inferred",
            "missing_conventions": sorted(set(issues)),
            "authority_level_for_printed_own_contract": 3
            if "pdf_page" in record["location"]
            else None,
            "own_contract_authority_qualification": "publication is evidence for its printed cells, not proof of complete native PCFF contract; the archive has no ranked compatible-PCFF contract or established native author/version lineage",
            "authority_level_for_J70_target": 5,
            "exact_native_environment": False,
            "native_pcff_authority": False,
            "selection_basis": "source_evidence",
            "execution_eligible": False,
        }
    )


def screening_key(record, cell):
    elements = candidate_sites(cell)
    if not elements or None in elements:
        return None
    family = record["interaction_family"]
    if family == "angle-angle":
        return family, tuple(
            min(elements, [elements[3], elements[1], elements[2], elements[0]])
        )
    return family, tuple(min(elements, elements[::-1]))


def requirements(target):
    evidence = [
        "hash-bound exact native FRC row or exact same-contract primary artifact/publication",
        "explicit c3/sio/osi/hc typing and complete ordered stable-site mapping; element resemblance is insufficient",
        "complete scalar equation, analytic derivative, units, prefactor, phase, multiplicity and reversal contract",
        "bind all source equilibrium dependencies to their exact native row identities",
    ]
    if target["family"] == "angle-angle-torsion_1":
        evidence.append(
            "authoritative native literal Phi versus cosine AAT transformation, or explicit source correction; no assumed equivalence"
        )
    if target["family"] == "angle-angle":
        evidence.append(
            "AA center j/shared k and each M component/reference pair; no unrestricted four-site reversal"
        )
    return evidence


def convention_findings(records, original_findings):
    """Separate an observed mismatch from an unproved suspected source error."""
    si = lambda r: r["source_file_sha256"].startswith("6e86aa640854")
    definitions = [
        (
            "Si_C_H_reference_order",
            ["equilibrium_dependency_conflict", "insufficient_evidence"],
            lambda r: si(r) and r["record_id"] in ("cell-0241", "cell-0242"),
            "Obtain the ordered bond-equilibrium and left/right coefficient contract or an explicit source revision. Do not swap references by chemical intuition.",
        ),
        (
            "SI_improper_dihedral_command_mapping",
            ["equation_conflict", "source_authority_conflict", "insufficient_evidence"],
            lambda r: (
                si(r)
                and r["source_force_field_family"] == "COMPASS"
                and r["location"].get("row_label") in ("si_c_c_h", "h_c_h_h")
            ),
            "Obtain an authoritative improper/AA mapping and M/reference ordering; dihedral commands after an improper style are not that mapping.",
        ),
        (
            "UA_nonbonded_contract",
            ["equation_conflict", "insufficient_evidence"],
            lambda r: (
                r["source_force_field_family"] == "UA"
                and r["interaction_family"] == "LJ12_6"
            ),
            "Resolve printed hybrid/lj-class2 commands versus archive LJ12-6 and all pair-specific conventions. The J72 blanket-LJ comparison guard does not rewrite J71.",
        ),
        (
            "PolCA_OH_H_printed_length",
            ["insufficient_evidence"],
            lambda r: r["record_id"] == "cell-0133",
            "Confirm the printed .945 nm against an author correction or exact parameter artifact. No alternative decimal is supplied, so a unit/transcription error is suspected, not established.",
        ),
        (
            "native_Phi_versus_cosine_AAT",
            ["equation_conflict", "insufficient_evidence", "source_authority_conflict"],
            lambda r: r["interaction_family"] == "angle-angle-torsion_1",
            "Establish native Phi/AAT equation, units, phase and dependency authority. Literal Phi and cos(Phi) are different functions; no equivalent transformation or source correction is supplied.",
        ),
    ]
    require(len(original_findings) == len(definitions), "J72_convention_inventory")
    return [
        signed(
            {
                "id": name,
                "J71_finding": original,
                "categories": categories,
                "candidate_record_ids": [
                    r["record_id"] for r in records if predicate(r)
                ],
                "status": "unresolved",
                "selected_record": None,
                "required_evidence": evidence,
                "native_PCFF_rows_resolved": 0,
                "numeric_comparison_performed": False,
            }
        )
        for (name, categories, predicate, evidence), original in zip(
            definitions, original_findings, strict=True
        )
    ]


def reports(repo, root):
    plan, records, cells, sources, matrix = inputs(repo, root)
    normalized = [
        comparison_record(r, c, sources[r["source_file_sha256"]])
        for r, c in zip(records, cells, strict=True)
    ]
    comparisons = []
    for i, j in combinations(range(len(records)), 2):
        key = screening_key(records[i], cells[i])
        if key is not None and key == screening_key(records[j], cells[j]):
            comparisons.append(compare(normalized[i], normalized[j]))
    by_id = {r["record_id"]: r for r in normalized}
    old_conflicts = read(repo / PARENT / "source-conflicts.json")
    conventions = convention_findings(
        records, old_conflicts["other_convention_obstructions"]
    )
    duplicate_replays = []
    for duplicate in old_conflicts["duplicate_sources"]:
        a, b = (by_id[k] for k in duplicate["records"])
        checked = compare(a, b)
        require(
            checked["numerically_compared"]
            and checked["equal_with_declared_arithmetic_allowance"],
            "J72_duplicate_replay_failure",
        )
        # Independent verification: equality of source decimal cells in the
        # same printed unit is sufficient, without normalized-cell arithmetic.
        require(
            a["raw_unit"] == b["raw_unit"] and a["raw_value"] == b["raw_value"],
            "J72_independent_duplicate_replay",
        )
        duplicate_replays.append(
            signed(
                {
                    "schema": "island.j72-equivalence-resolution.v1",
                    "alternatives": [a["record_id"], b["record_id"]],
                    "alternative_J71_identities": [
                        a["J71_identity"],
                        b["J71_identity"],
                    ],
                    "citations": [
                        {
                            "citation": x["citation"],
                            "location": x["location"],
                            "source_file_sha256": x["source_file_sha256"],
                        }
                        for x in (a, b)
                    ],
                    "transformation": "identity on C0..C5 in the same PolCA2021 RB Eq11, ordered sites and kJ/mol; comparison E_kcal=E_kJ/4.184",
                    "comparison": checked,
                    "independent_replay": "source raw decimal strings and unit arrays agree exactly before normalization",
                    "resolution": "duplicate_equivalence_confirmed; no unequal numeric conflict repaired",
                    "native_PCFF_rows_resolved": 0,
                    "target_instances_resolved": 0,
                }
            )
        )
    findings = []
    for conflict in old_conflicts["unresolved_conflicts"]:
        a, b = (by_id[r["record_id"]] for r in conflict["records"])
        checked = compare(a, b)
        require(
            checked["numerically_compared"]
            and "same_row_numeric_conflict" in checked["categories"],
            "J72_nominal_angle_conflict_not_reproduced",
        )
        findings.append(
            signed(
                {
                    "id": conflict["id"],
                    "alternatives": [a["record_id"], b["record_id"]],
                    "comparison": checked,
                    "categories": [
                        "same_row_numeric_conflict",
                        "source_authority_conflict",
                    ],
                    "status": "source_conflict_unresolved",
                    "scope": "nominal benchmark COMPASS scalar angle assignment, not native PCFF row authority",
                    "highest_authority": "published SI and explicitly associated archive disagree; no authoritative revision/correction supplied",
                    "selected_record": None,
                    "required_evidence": conflict["required_resolution"],
                }
            )
        )
    mappings = read(repo / PARENT / "mapping-report.json")
    parent_by_id = {t["id"]: t for t in matrix["instances"]}
    rows = []
    for old in mappings["instances"]:
        target = parent_by_id[old["id"]]
        candidates = []
        for m in old["candidate_mappings"]:
            r = by_id[m["record_id"]]
            categories = {"chemical_environment_conflict", "source_authority_conflict"}
            if r["contract"]["forcefield_family"] != "PCFF":
                categories.add("forcefield_family_conflict")
            if r["contract"]["resolution"] != "all_atom":
                categories.add("resolution_conflict")
            if r["missing_conventions"]:
                categories.update(
                    ("equilibrium_dependency_conflict", "insufficient_evidence")
                )
            if target["family"] == "angle-angle-torsion_1":
                categories.update(
                    ("phase_conflict", "equation_conflict", "unit_conflict")
                )
            if r["contract"]["interaction_family"] == "class2_AA_triplet":
                categories.add("equation_conflict")
            candidates.append(
                {
                    "record_id": r["record_id"],
                    "J71_identity": r["J71_identity"],
                    "comparison_identity": r["identity"],
                    "J71_mapping_status": m["status"],
                    "mapping_basis": "inherited element/order screening only, not an exact native source lookup",
                    "authority_level_for_target": 5,
                    "compatibility_classification": "analog_unqualified",
                    "conflict_categories": sorted(categories),
                    "exact_native": False,
                    "external_complete_row": False,
                    "missing_conventions": r["missing_conventions"],
                    "required_evidence": requirements(target),
                }
            )
        if not candidates:
            status = "source_missing"
        elif target["family"] == "angle-angle-torsion_1":
            status = "convention_unresolved"
        else:
            status = "analog_unqualified"
        best = min((r["authority_level_for_target"] for r in candidates), default=None)
        rows.append(
            signed(
                {
                    "id": target["id"],
                    "J70_identity": target["identity"],
                    "J71_instance_identity": old["source_identity"],
                    "J71_status": "unresolved",
                    "J71_resolved": old["resolved"],
                    "family": target["family"],
                    "stable_sites": target["stable_sites"],
                    "ordered_types": target["ordered_types"],
                    "symmetry_equivalent_types": target["symmetry_equivalent_types"],
                    "source_lookup": target["source_lookup"],
                    "equilibrium_dependencies": target["dependencies"],
                    "candidate_records": candidates,
                    "conflict_categories": sorted(
                        {c for m in candidates for c in m["conflict_categories"]}
                    )
                    or ["insufficient_evidence"],
                    "highest_authority_level_for_target": best,
                    "highest_authority_records": [
                        r["record_id"]
                        for r in candidates
                        if r["authority_level_for_target"] == best
                    ],
                    "selected_source": None,
                    "exact_enough_for_native_PCFF": False,
                    "external_research_complete_row": False,
                    "remaining_missing_conventions": requirements(target),
                    "required_additional_evidence": requirements(target),
                    "QM_fit_required": False,
                    "QM_fit_status": "not_established; source absence does not authorize fitting; inherited exact J70 targets remain unexecuted",
                    "final_J72_status": status,
                    "resolved": False,
                }
            )
        )
    require(
        len(rows) == len({r["id"] for r in rows}) == 378, "J72_exact_matrix_coverage"
    )
    categories = Counter(c for p in comparisons for c in p["categories"])
    summary = {
        "J71_records_preserved": len(records),
        "J71_numerical_cells_preserved": sum(
            v is not None for r in records for v in r["raw_value"]
        ),
        "comparison_pairs": len(comparisons),
        "numerically_compatible_pairs": sum(
            p["numerically_compared"] for p in comparisons
        ),
        "equal_scalar_pairs": sum(
            p["equal_with_declared_arithmetic_allowance"] is True for p in comparisons
        ),
        "pair_conflict_category_counts_overlapping": {
            c: categories[c] for c in CATEGORIES
        },
        "documented_unequal_source_assignment_conflicts": len(findings),
        "convention_obstruction_counts_overlapping": dict(
            Counter(c for f in conventions for c in f["categories"])
        ),
        "unequal_source_conflicts_mathematically_resolved": 0,
        "real_duplicate_equivalences_replayed": len(duplicate_replays),
        "synthetic_controls_count_as_real_resolutions": False,
        "exact_target_instances_resolved": 0,
        "native_PCFF_instances_resolved": 0,
        "external_complete_research_instances": 0,
        "non_native_reference_records": sum(
            r["source_force_field_family"] != "PCFF" for r in records
        ),
        "PCFF_publication_reference_records": sum(
            r["source_force_field_family"] == "PCFF" for r in records
        ),
        "analog_only_instances": sum(
            r["final_J72_status"] == "analog_unqualified" for r in rows
        ),
        "instances_without_candidate_records": sum(
            not r["candidate_records"] for r in rows
        ),
        "instances_unresolved": len(rows),
        "row_status_counts": dict(Counter(r["final_J72_status"] for r in rows)),
        "instance_conflict_category_counts_overlapping": dict(
            Counter(c for r in rows for c in r["conflict_categories"])
        ),
        "rows_requiring_exact_primary_or_author_evidence": 378,
        "rows_definitely_requiring_QM_fitting": 0,
        "rows_explicitly_blocked_by_native_AAT_equation_phase_units": 60,
        "additional_order_and_dependency_blockers": "all candidate mappings require complete native equilibrium/source typing; classifications overlap",
    }
    bundles = []
    for model in sorted({r["contract"]["forcefield_family"] for r in normalized}):
        catalog = [r for r in normalized if r["contract"]["forcefield_family"] == model]
        bundles.append(
            {
                "source_family": model,
                "record_ids": [r["record_id"] for r in catalog],
                "original_J71_identities": [r["J71_identity"] for r in catalog],
                "source_hashes": sorted({r["source_file_sha256"] for r in catalog}),
                "resolutions": sorted({r["contract"]["resolution"] for r in catalog}),
                "role": "separate non-executable source evidence catalog; PCFF publication rows are not a native FRC replacement",
                "native_rows": [],
                "external_complete_rows": [],
                "parameter_completeness": "incomplete",
                "coherent_complete_contract": False,
                "execution_profile": None,
            }
        )
    return {
        "parent-reconstruction.json": signed(
            {
                "schema": "island.j72-parent-reconstruction.v1",
                "required_J71": REQUIRED,
                "J71_reading_evidence": {
                    "main_papers": 12,
                    "supporting_documents": 3,
                    "pages": 287,
                    "extracted_rows": 467,
                    "nonblank_numeric_cells": 1470,
                    "new_J72_full_document_read_claims": 0,
                    "original_hash_recheck_command": "scripts/read_phase_4j71_sources.py",
                    "original_hash_recheck_scope": "25 files and 124 exact archive rows; optional originals check, distinct from offline replay",
                },
                "source_files": [
                    {
                        k: s.get(k)
                        for k in (
                            "relative_path",
                            "sha256",
                            "title",
                            "authors",
                            "journal_year",
                            "doi",
                            "actually_read_successfully",
                        )
                    }
                    for s in read(repo / PARENT / "reading-status.json")["files"]
                ],
                "graph_identity": matrix["graph_identity"],
                "atom_order": matrix["atom_order"],
                "parent_pins": plan["parent_pins"],
                "observation_policy": "historical 26 observations and original producers unchanged; J68 incompatibility independently retained, not refitted",
            }
        ),
        "execution-scope-report.json": execution_deviation(root),
        "comparison-normalized-records.json": signed(
            {
                "schema": "island.j72-normalized-records.v1",
                "records": normalized,
                "raw_parent_sha256": plan["parent_pins"][
                    str(PARENT / "parameter-records.json")
                ],
                "normalization": "separate comparison fields only; 48-digit Decimal arithmetic with fixed pi; raw source strings unchanged",
            }
        ),
        "authority-hierarchy.json": signed(
            {
                "schema": "island.j72-authority.v1",
                "levels": {str(k): v for k, v in AUTHORITY.items()},
                "native_eligible_levels": [1, 2, 3],
                "external_research_only_levels": [4],
                "never_native_levels": [5, 6],
                "tie_policy": "retain all alternatives; no averaging, residual/MD selection or unproved preference",
                "exact_primary_override": None,
            }
        ),
        "source-conflict-matrix.json": signed(
            {
                "schema": "island.j72-conflicts.v1",
                "comparisons": comparisons,
                "documented_source_disagreements": findings,
                "inherited_convention_obstructions": old_conflicts[
                    "other_convention_obstructions"
                ],
                "summary": summary,
                "comparison_scope": "all unordered pairs sharing an inherited element/order interaction signature; structural categories are screening exclusions, not proof of erroneous source values; all other pairs not_applicable; all 467 records normalized",
                "values_averaged": False,
                "preferred_values_selected": False,
            }
        ),
        "mathematical-resolutions.json": signed(
            {
                "schema": "island.j72-resolutions.v1",
                "resolutions": duplicate_replays,
                "unequal_conflicts_resolved": [],
                "scope": "real duplicate identity transformations only; manufactured unit/radian controls are in tests and are not source resolutions",
            }
        ),
        "convention-adjudications.json": signed(
            {
                "schema": "island.j72-convention-adjudications.v1",
                "findings": conventions,
                "classification_policy": "categories describe demonstrated representation differences and explicitly missing evidence; suspected printed-value error alone is insufficient_evidence, not a proved unit correction",
            }
        ),
        "row-adjudication.json": signed(
            {
                "schema": "island.j72-row-adjudication.v1",
                "graph_identity": matrix["graph_identity"],
                "atom_order": matrix["atom_order"],
                "parent_matrix_identity": matrix["identity"],
                "instances": rows,
                "family_counts": dict(Counter(r["family"] for r in rows)),
                "summary": summary,
            }
        ),
        "unresolved-evidence-requirements.json": signed(
            {
                "schema": "island.j72-evidence-requirements.v1",
                "instances": [
                    {
                        "id": r["id"],
                        "status": r["final_J72_status"],
                        "required_evidence": r["required_additional_evidence"],
                        "candidate_record_ids": [
                            m["record_id"] for m in r["candidate_records"]
                        ],
                        "QM_fit_required": False,
                    }
                    for r in rows
                ],
                "author_conflicts": findings,
                "QM_policy": "No new QM/fitting allocation. Exact inherited J70 target templates remain unexecuted; first resolve source/convention authority.",
            }
        ),
        "candidate-source-bundles.json": signed(
            {
                "schema": "island.j72-candidate-bundles.v1",
                "catalogs": bundles,
                "parameter_completeness": "incomplete",
                "remaining_row_ids": [r["id"] for r in rows],
                "resolved_instances": 0,
                "zero_fill": False,
                "wildcard_fallback": False,
                "cross_family_overlay": False,
                "execution_profile": None,
            }
        ),
        "strict-converter-report.json": signed(
            {
                "schema": "island.j72-converter-gate.v1",
                "eligible": False,
                "executed": False,
                "attempts": 0,
                "ignore_used": False,
                "complete": False,
                "reason": "378 unresolved instances; strict conversion not eligible; no permissive conversion",
                "historical_failure": "J70 returncode17: Unable to find endbonddihedral data for hc c3 sio c3; not a new J72 attempt",
            }
        ),
        "static-closure-report.json": signed(
            {
                "schema": "island.j72-static-closure-gate.v1",
                "eligible": False,
                "executed": False,
                "complete": False,
                "retained_fixture_count": 26,
                "reason": "No complete coherent source bundle; available-subset historical closure is distinct from complete static closure",
            }
        ),
        "canonical-capability.json": signed(
            {
                "schema": "island.j72-canonical-capability.v1",
                "ordinary_api": 'prepare_polymer(..., forcefield="pcff")',
                "ordinary_behavior_modified": False,
                "registered_profiles": [],
                "automatic_assignment": False,
                "manual_profile_required": False,
                "cross_family_fallback": False,
                "missing_row_ids": [r["id"] for r in rows],
                "source_status": "source_conflict_unresolved",
                "execution_profile_created": False,
            }
        ),
        "decision.json": signed(
            {
                "schema": "island.j72-decision.v1",
                **DECLARATIONS,
                "summary": summary,
                "converter_complete": False,
                "static_complete_closure": False,
                "execution_profile_created": False,
                "J68_incompatibility": "preserved independently; literature adjudication does not change its observations, model span or diagnostic outcomes",
                "bounded_next_step": "Obtain the exact native Si/O Class-II rows, typing, equilibrium dependencies and Phi/AAT authority; ask the benchmark authors for corrected SI/archive angle assignments and reference order. No new QM, fitting or MD is justified by this authority audit alone.",
            }
        ),
    }


def human_report(result):
    s = result["decision.json"]["summary"]
    out = [
        "# J72 source conflict adjudication",
        "",
        "J71's 12 papers and three supporting documents (287 pages) are retained reading evidence, not 15 newly read documents in J72. All 467 extraction identities and 1,470 numerical cells are preserved. Original file hashes and archive lines are checked separately from offline reconstruction. No new parameter extraction or HMDS fitting/validation. The initial historical suite violated the zero-execution scope: four SQM AM1 minimizations and 46 recorded MD steps, with up to two uncommitted steps possible. These J31 CHO outputs are quarantined and unused in adjudication; execution-scope-report.json reconstructs the actual counts.",
        "",
        "## Counts and scope",
        "",
    ]
    out += [f"- {k}: {v}" for k, v in s.items()]
    out += [
        "",
        "Pair category counts overlap. They describe structural compatibility screening, not a count of incorrect publications. Numerical comparison occurs only after compatible scalar source contracts are established. The benchmark comparison scope is its nominal COMPASS labels and archive topology, never a native PCFF equivalence contract.",
        "",
        "## Unresolved printed assignments",
        "",
        "| Evidence | Alternatives | Result | Required evidence |",
        "|---|---|---|---|",
    ]
    for f in result["source-conflict-matrix.json"]["documented_source_disagreements"]:
        out.append(
            f"| {f['id']} | {', '.join(f['alternatives'])} | {f['status']}; neither selected | {f['required_evidence']} |"
        )
    out += [
        "",
        "The C–Si–O SI row is 113.1855/36.2069/−20.3939/20.0172 while topology-bound archive type102 (O–Si–C) is 114.9060/23.0218/−31.3993/24.9814. C–Si–C has the opposite assignments (archive type29). The order reversal of a scalar angle is exact; it does not prove that swapping SI labels is an authoritative correction. No explicit revision is supplied.",
        "",
        "## Mathematical equivalence and remaining conventions",
        "",
        "Sixteen PolCA RB main/SI pairs agree exactly in their original kJ/mol cells and replay after E_kcal=E_kJ/4.184. These are duplicate confirmations, not repaired unequal conflicts or PCFF rows. Original alternatives and exact citation locations are retained in mathematical-resolutions.json. Degree/radian, coefficient angular-power and energy conversions have independent manufactured controls; those controls are not real-source corrections.",
        "",
    ]
    out += [
        f"- {x}"
        for x in result["source-conflict-matrix.json"][
            "inherited_convention_obstructions"
        ]
    ]
    out += [
        "",
        "The UA printed lj/class2 command versus blanket LJ12-6 description is additionally blocked in comparison metadata, without changing J71 records. Unknown temperature, source target, charge/exclusion contract and equilibrium dependencies remain unknown. No equality or numeric disagreement is inferred across those gaps.",
        "",
        "## Authority and missing instances",
        "",
        "Native FRC > same-contract author-primary > same-contract published table > explicitly compatible external research > analog/transfer > unexecuted QM target. A publication's authority over its own printed cells does not establish exact c3/sio/osi/hc mapping. Native PCFF/CFF91, CFF93, COMPASS, Smith PDMS, HUANG, OPLS, GAFF, PolCA/UA, coarse-grained and QM corrections remain separate contracts. No cross-family overlay is authorized.",
        "",
        "All 378 instances retain exact graph sites, original identities, source lookups, symmetry and equilibrium dependencies in row-adjudication.json. Element-only candidates are analogs. AAT additionally lacks native literal-Phi/cosine equation, phase and unit authority. Missing rows are not zeros. Source absence is not proof that a row requires QM fitting; exact J70 target templates remain unexecuted.",
        "",
        "## Decision and capability",
        "",
        'Decision: **source_conflict_unresolved**. Zero exact native or complete external instances resolved. All 378 remain unresolved. No execution profile, strict conversion, complete static closure or parameter assignment capability. Ordinary prepare_polymer(..., forcefield="pcff") remains unchanged and fails closed through its existing capability checks. All coefficients remain unaccepted. Chemical transfer, scientific validation and production readiness remain unestablished. Obtain exact primary rows and an authoritative benchmark correction before considering a separately scoped next phase.',
        "",
    ]
    return "\n".join(out)


def replay(repo, root, deep_parent=False):
    expected = reports(repo, root)
    for name, value in expected.items():
        require(read(root / name) == value, "J72_semantic_reconstruction", path=name)
    require(
        (root / "source-conflict-report.md").read_text() == human_report(expected),
        "J72_human_report_mutation",
    )
    if deep_parent:
        from scripts.validate_phase_4j71_receipt import validate

        validate(repo, repo / PARENT)
    return {
        "offline": True,
        **DECLARATIONS,
        "records": 467,
        "numerical_cells": 1470,
        "instances": 378,
        "resolved_instances": 0,
        "real_duplicate_equivalences": 16,
        "unequal_conflicts_resolved": 0,
        "converter_executed": False,
        "static_complete_closure": False,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--replay", action="store_true")
    p.add_argument("--deep-parent", action="store_true")
    a = p.parse_args()
    repo, root = a.repo.resolve(), a.root.resolve()
    if a.replay:
        print(replay(repo, root, a.deep_parent))
    else:
        result = reports(repo, root)
        for name, value in result.items():
            write(root / name, value)
        (root / "source-conflict-report.md").write_text(human_report(result))
        print(result["decision.json"]["summary"])
