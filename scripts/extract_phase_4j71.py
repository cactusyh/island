"""Offline reconstruction of J71 reviewed literature cells; no scientific execution."""

import argparse
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.literature_intake import (
    candidate_sites,
    coherent_bundle_gate,
    compatibility,
    make_record,
    semantic_declarations,
)

ROOT = Path("docs/evidence/phase_4j71")
BASE = "4c9066e4501591211049c44401e82b3110118f58"
REQUIRED = "9040f4bf6b9948c68dabcd85d7dd25a40e5a6536"
BRANCH = "codex/phase-4j71-pcff-literature-parameter-extraction"
PINS = {
    "source-inventory.json": "2960f685489b8edd284d6520f8c7636f00ac73eaddbb06f33d964d0de9353e78",
    "extraction-ledger.json": "1b5d395a48027d523119cd22b57b3e2e7c2b32090e1cdadd7e5f2ba58a0e5076",
    "reading-status.json": "4bb39ab0b7c645a49b6e4640a1e89d2630f358ea45fe0899def8cdaacb934d7b",
    "archive-members.json": "6d6fba418bdc9485ceae85f816e0ca63d00865dfb2685fdbc37317cd24767715",
}
DECLARATIONS = {
    "decision": "source_conflict",
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


def inputs(repo, root):
    require(
        pin(root / "contract.json")
        == "7ddc92287514a14c8d846e336948de5fa82f84200372a86a07b2d5a30471748a",
        "J71_contract_mutation",
    )
    for name, expected in PINS.items():
        require(pin(root / name) == expected, "J71_frozen_review_mutation", path=name)
    plan = read(root / "contract.json")
    require(
        plan["base_commit"] == BASE
        and plan["required_J70"] == REQUIRED
        and plan["source_pins"] == PINS
        and plan["new_QM"] == plan["MD_steps"] == plan["QM_fits"] == 0,
        "J71_contract_semantics",
    )
    for name, expected in plan["parent_pins"].items():
        require(pin(repo / name) == expected, "J71_parent_mutation", path=name)
    ledger = read(root / "extraction-ledger.json")
    reading = read(root / "reading-status.json")
    sources = {r["sha256"][:12]: r for r in reading["files"]}
    # The complete primary copy owns the metadata for byte-identical checkpoints.
    for r in reading["files"]:
        if r["scope"] == "primary_directory":
            sources[r["sha256"][:12]] = r
    require(
        len({c["cell_id"] for c in ledger["cells"]}) == len(ledger["cells"]),
        "J71_duplicate_cells",
    )
    return plan, ledger["cells"], reading, sources


def source_conflicts(cells, records):
    by_id = {c["cell_id"]: r for c, r in zip(cells, records, strict=True)}

    def select(source, command, label, sub=None):
        return next(
            c
            for c in cells
            if c["source"] == source
            and c.get("command") == command
            and c["location"]["row_label"] == label
            and c.get("substyle") == sub
        )

    findings = []
    for label, archive_id in (("c_si_o", "102"), ("c_si_c", "29")):
        si = select("6e86aa640854", "angle_coeff", label)
        arc = select("8b0e8b74303b", "angle_coeff", archive_id)
        findings.append(
            {
                "id": "benchmark_SI_archive_angle_" + label,
                "status": "source_conflict",
                "records": [by_id[c["cell_id"]] for c in (si, arc)],
                "interpretation": "SI label and topology-bound archive assign different angle coefficients; source identity preserved, no relabeling or selection",
                "required_resolution": "author correction of TablesS3 labels and topology/source-type mapping; full source/converter contract",
                "strongest_PCFF_authority": "Sun/Rigby1997 primary siloxane Appendix for its printed scalar rows; no authority for omitted couplings",
            }
        )
    duplicates = []
    for c, r in zip(cells, records, strict=True):
        if c["source"] != "5372747c08f6" or c["family"] != "RB_torsion":
            continue
        other = next(
            o
            for o in records
            if o["source_file_sha256"].startswith("1e95027a1202")
            and o["ordered_sites"] == r["ordered_sites"]
        )
        require(other["raw_value"] == r["raw_value"], "PolCA_duplicate_disagreement")
        duplicates.append(
            {
                "records": [r["record_id"], other["record_id"]],
                "status": "duplicate_equal_values_preserved",
            }
        )
    return signed(
        {
            "schema": "island.j71-source-conflicts.v1",
            "unresolved_conflicts": findings,
            "duplicate_sources": duplicates,
            "other_convention_obstructions": [
                "SI COMPASS Si-C-H reference lengths are reversed relative to printed ordered sites; archive topology cannot by itself resolve converter orientation",
                "SI after improper_style class2 prints dihedral_coeff si_c_c_h and h_c_h_h; these are not an authoritative improper/AA mapping",
                "SI UA hybrid 9-6/12-6 description versus archive lj/cut/coul/cut 12-6",
                "PolCA Table3 OH-H 0.945 nm retained as printed; decimal/unit issue unresolved, no inferred correction",
                "native FRC literal Phi/AAT authority versus Sun1995 cosine AAT is an inherited unresolved convention, not silently reconciled",
            ],
            "different_models_not_actual_conflicts": [
                "Sun1995 polysilane Si-C r0=1.899 versus more precise archive1.8995 is different reporting/model context",
                "Sun1997 Si4 epsilon0.198/rmin4.405 versus Si4c epsilon0.131/rmin4.290 are different explicitly typed environments",
                "Smith2004 Si-O r0=1.651 and Si-O-Si theta0=137.63 versus Sun1997 1.640/159.0 use different targets and exp6/1-4 conventions",
                "PolCA Si-OB-Si149.5 and RB torsions belong to united-atom model, not a competing native PCFF row",
            ],
            "values_averaged": False,
            "preferred_values_selected": False,
        }
    )


# Locations established during the full reading, including sources that cannot
# supply an HMDS Si/O row. Property/geometry tables are not force constants.
TABLE_REVIEW = {
    "00b279decf05": "PDF pp5,7–15 TablesI–VIII; atom types/increments, ester/aromatic valence and conformational/vibrational tests. No Si-containing parameter table; carbonate/ester O is not osi.",
    "ac85e79667e1": "PDF p3 CFF93 energy equations; p9 Table10 force-field parameters, Tables1–9 structure/energy/frequency tests. Carbonate/aromatic sites, no siloxane rows.",
    "9a62dafe2236": "PDF pp2–3 functional forms; p4 Table1 typed equivalences; p7 Table2 increments; p16 Table8 LJ comparisons; AppendixB pp24–25 alkane/benzene numeric tables. Not a native PCFF siloxane coupling source; Si4c/O2z typing is not an element wildcard.",
    "ee8614dcf423": "PDF pp3–4 method and typing; Tables1–26 model compounds and validation; Appendix pp16–18 contains P/N/organic valence/nonbonded and couplings. No Si/O siloxane row authority.",
    "f9277e1c608f": "PDF pp3–5 method; Appendix pp17–20 atom equivalences, increments, quartic bonds/angles, torsions, nonbond and cross terms for alcohols/ethers. Oxygen o2h/o2e and alpha-carbon c4o must not be equated to siloxane osi; original text font uses special decimal glyphs, no OCR-based PCFF transfer.",
    "a8a4961f8c1b": "PDF pp3–5 equations and ESP/refinement; p5 Table4 charges; AppendixA/Table12 pp11–12 nitrate types, increments, valence and coupling parameters. Nitrate/ester oxygen is chemically distinct from bridging siloxane O.",
    "c518da36d84d": "All7 pages. Tables1–4 crystal geometry and RIS/polymer outputs, not a new numerical force-constant set. pcff2 is explicitly COMPASS1.0; fitted RIS energies are not Class-II interaction coefficients.",
    "04ecee672431": "All10 pages. Force-field comparisons and aggregate property tests; Fig7 charge plots are not digitized. Ref14 Huang2017 primary source is cited but not supplied; Ref22 binds the local associated benchmark archive.",
    "451cf754daab": "All13 rendered pages manually read, no OCR. PDF p2 energy functions; p3 Table1 nonbonded organic atoms; p9 Table6 compares AMBER/CHARMM22 and is not OPLS source authority; p10 Table7 alkane torsions. No Si, all-atom organic OPLS is not PCFF.",
    "65bda3a5b14c": "All23 pages. Section2/PDF p3 Eq1a–b scalar energy, Eq2 increments; Appendix pp21–22 has28 extracted rows. General cross terms in Eq1b do not supply missing numeric AA/EBT/MBT/AT/AAT coefficients; absence is not zero. PDF p16 trans=180 degrees.",
    "9c6d6315b91e": "All12 pages. SectionIII/PDF p3 scalar functions; p4 Table5 refined nonbond values; Appendix p11 has77 extracted rows. PDMS here means poly(dimethylsilane), a Si–Si backbone; shared Si–C/H terms do not establish siloxane-equilibrium dependencies.",
    "a1415925e044": "All11 pages. PDF p3 Eqs1–5 and Table3 molecular charges (HMDS column separately retained); p5 Table4 nonbond/bond/bend/torsion constants and footnotes. exp6, scaled 1–4 and Eq6 sine-squared modulation define another model, not missing Class-II rows.",
    "5372747c08f6": "All16 pages. Section2.3/PDF p5 Eqs9–15 rigid bonds, half-prefactor harmonic angles, RB torsions, LJ12-6/Lorentz–Berthelot; p6 Tables1–2 intermediate DDEC versus final model charges; p7 Tables3–4, p8 Table5, p9 Table6. Si charge-by-neutrality remains specific to PolCA and is not imported.",
    "1e95027a1202": "All76 pages, including rotation-aware extraction of pp41–76. Property-data TablesS1–S71 are fitting/reference data rather than PCFF coefficients; sectionsS3–S4 discuss bond geometry and torsional fits; p65 TableS72 duplicates all16 final RB rows. No energy curve digitization or new fitting.",
    "6e86aa640854": "All14 pages. TableS1 p2 UA; S2 p3 HUANG; S3 pp4–5 COMPASS; S4 p6 OPLS; S5 p7 Dreiding. All numeric model commands were extracted. Conflicting labels, duplicated commands and the improper/dihedral command inconsistency are preserved, not repaired.",
}


def reading_notes(reading, records):
    notes = []
    for source in reading["files"]:
        prefix = source["sha256"][:12]
        if source["scope"] != "primary_directory" or prefix not in TABLE_REVIEW:
            continue
        source_records = [
            r for r in records if r["source_file_sha256"] == source["sha256"]
        ]
        notes.append(
            {
                "source_file_sha256": source["sha256"],
                "citation": {
                    k: source.get(k)
                    for k in ("title", "authors", "journal_year", "doi")
                },
                "page_coverage": source["page_coverage"],
                "force_field_identity": source["force_field"],
                "model_resolution": source["resolution"],
                "chemical_domain": source["chemical_domain"],
                "table_and_equation_review": TABLE_REVIEW[prefix],
                "parameter_candidate_records": [r["record_id"] for r in source_records],
                "scalar_functions_in_extracted_rows": sorted(
                    {r["scalar_energy_equation"] for r in source_records}
                ),
                "exclusion_policy": "different chemical environments, model resolution and force-field families remain separate; non-target tables are indexed rather than transferred",
            }
        )
    return signed(
        {
            "schema": "island.j71-complete-reading-notes.v1",
            "documents": notes,
            "unique_primary_PDFs_read": 15,
            "primary_PDF_pages_read": 287,
            "main_papers": 12,
            "supporting_documents": 3,
            "full_papers_committed": False,
            "OCR_used": False,
            "archive_review_scope": "all member bytes hashed; topology and applicable coefficient rows parsed; README/small inputs read; large generic unused tables and production coordinates not claimed as manually read papers",
        }
    )


def reports(repo, root):
    plan, cells, reading, sources = inputs(repo, root)
    records = [make_record(c, sources[c["source"]]) for c in cells]
    parent = read(repo / "docs/evidence/phase_4j70/completeness-matrix.json")
    matrix = parent["instances"]
    require(
        len(matrix) == 378 and len({r["id"] for r in matrix}) == 378,
        "J71_parent_378_instances",
    )
    candidate_rows, evidence = [], []
    for cell, record in zip(cells, records, strict=True):
        outcomes = []
        for target in matrix:
            status = compatibility(record, target, candidate_sites(cell))
            if status != "not_applicable":
                row = {
                    "record_id": record["record_id"],
                    "target_id": target["id"],
                    "graph_sites": target["stable_sites"],
                    "target_ordered_types": target["ordered_types"],
                    "symmetry_class": target["family"]
                    + ":"
                    + ",".join(min(target["symmetry_equivalent_types"])),
                    "status": status,
                    "coherent_contract": False,
                    "reason": "element/order comparison only; source family, exact environment, native phase and equilibrium dependency contract not interchangeable",
                }
                candidate_rows.append(row)
                outcomes.append(row)
        evidence.append(
            {
                **record,
                "target_J70_row_ids": [m["target_id"] for m in outcomes],
                "mapping_statuses": sorted({m["status"] for m in outcomes})
                or ["not_applicable"],
                "required_author_confirmation": "exact siloxane typing; ordered equilibrium references; native versus COMPASS provenance; Phi/AAT and multiplicity contract",
                "QM_fit_requirement": "unresolved source authority does not authorize fitting; exact retained J70 targets remain unexecuted",
            }
        )
    complete = coherent_bundle_gate(matrix, candidate_rows, records)
    require(not complete, "J71_unexpected_source_complete")
    instances = [
        {
            "id": t["id"],
            "family": t["family"],
            "stable_sites": t["stable_sites"],
            "ordered_types": t["ordered_types"],
            "source_identity": t["identity"],
            "symmetry_equivalent_types": t["symmetry_equivalent_types"],
            "dependencies": t["dependencies"],
            "candidate_mappings": [
                m for m in candidate_rows if m["target_id"] == t["id"]
            ],
            "exact_matches": [],
            "resolved": False,
            "remaining_evidence": "exact primary numerical source and complete convention/dependency contract, not another element-only analogue",
        }
        for t in matrix
    ]
    conflicts = source_conflicts(cells, records)
    bundle = signed(
        {
            "schema": "island.j71-partial-literature-bundle.v1",
            "parameter_completeness": "incomplete",
            "coherent_contract": False,
            "native_pcff_rows": [],
            "author_primary_rows": [],
            "external_executable_rows": [],
            "reference_catalog": [r["identity"] for r in records],
            "reference_catalog_role": "non-executable extracted evidence; not supplied complete interaction rows",
            "remaining_row_ids": [t["id"] for t in matrix],
            "zero_fill": False,
            "wildcard_fallback": False,
            "execution_profile": None,
            "automatic_assignment": False,
        }
    )
    decision = signed(
        {
            "schema": "island.j71-decision.v1",
            **DECLARATIONS,
            "extracted_candidate_rows": len(records),
            "extracted_nonblank_numeric_cells": sum(
                v is not None for r in records for v in r["raw_value"]
            ),
            "resolved_instances": 0,
            "remaining_instances": 378,
            "converter_complete": False,
            "static_complete_closure": False,
            "execution_profile_created": False,
            "bounded_next_step": "Resolve benchmark TablesS3/archive label and equilibrium-order contradictions with its authors; obtain exact native Si/O coupling rows and Phi/AAT authority. No new QM, fitting or generic sampling allocation. J68 incompatibility remains independently preserved.",
        }
    )
    semantic_declarations(decision)
    capability = signed(
        {
            "schema": "island.j71-canonical-capability.v1",
            "ordinary_api": 'prepare_polymer(..., forcefield="pcff")',
            "manual_profile_required": False,
            "manual_charge_method_required": False,
            "registered_profiles": [],
            "execution_profile_created": False,
            "automatic_assignment": False,
            "cross_family_fallback": False,
            "bounded_domain_complete": False,
            "missing_row_ids": [t["id"] for t in matrix],
            "source_status": "source_conflict_and_parameter_incomplete",
            "historical_live_API_check": {
                "path": "docs/evidence/phase_4j70/canonical-capability.json",
                "sha256": plan["parent_pins"][
                    "docs/evidence/phase_4j70/canonical-capability.json"
                ],
            },
            "ordinary_behavior_modified": False,
        }
    )
    return {
        "document-reading-notes.json": reading_notes(reading, records),
        "literature-file-inventory.json": signed(
            {
                "schema": "island.j71-enriched-file-inventory.v1",
                "files": reading["files"],
                "inventory_before_reading": PINS["source-inventory.json"],
            }
        ),
        "parameter-records.json": signed(
            {
                "schema": "island.j71-parameter-records.v1",
                "records": records,
                "accepted_coefficients": 0,
            }
        ),
        "source-evidence-table.json": signed(
            {"schema": "island.j71-source-evidence.v1", "rows": evidence}
        ),
        "mapping-report.json": signed(
            {
                "schema": "island.j71-378-mapping.v1",
                "parent_matrix_identity": parent["identity"],
                "graph_identity": parent["graph_identity"],
                "atom_order": parent["atom_order"],
                "instances": instances,
                "family_counts": dict(Counter(t["family"] for t in matrix)),
                "candidate_status_counts": dict(
                    Counter(m["status"] for m in candidate_rows)
                ),
                "resolved_instances": 0,
                "absence_policy": "All candidate/target pairs not present in the sparse lists have status not_applicable; no implicit exact mappings",
            }
        ),
        "source-conflicts.json": conflicts,
        "partial-literature-bundle.json": bundle,
        "strict-converter-report.json": signed(
            {
                "schema": "island.j71-converter-gate.v1",
                "eligible": False,
                "executed": False,
                "attempts": 0,
                "complete": False,
                "ignore_used": False,
                "reason": "no complete coherent bundle; J71 converter must not run",
                "historical_failure": "J70 returncode17: Unable to find endbonddihedral data for hc c3 sio c3; historical evidence unchanged",
            }
        ),
        "static-closure-report.json": signed(
            {
                "schema": "island.j71-static-closure-gate.v1",
                "eligible": False,
                "executed": False,
                "complete": False,
                "retained_fixture_count": 26,
                "reason": "missing rows and conventions; historical available-subset closure is not complete literature bundle closure",
            }
        ),
        "canonical-capability.json": capability,
        "qm-targets.json": signed(
            {
                "schema": "island.j71-unexecuted-targets.v1",
                "inherited_targets_file": "docs/evidence/phase_4j70/qm-fit-targets.json",
                "sha256": plan["parent_pins"][
                    "docs/evidence/phase_4j70/qm-fit-targets.json"
                ],
                "exact_target_ids": [
                    r["id"]
                    for r in read(
                        repo / "docs/evidence/phase_4j70/qm-fit-targets.json"
                    )["targets"]
                ],
                "status": "unexecuted_not_authorized; missing frozen fitting protocol",
                "source_missing_is_not_necessarily_QM_fit_required": True,
                "new_QM": 0,
                "QM_fits": 0,
            }
        ),
        "decision.json": decision,
    }


def human_table(reports):
    groups = [
        "Si–O bond terms",
        "Si–O–Si angle terms",
        "O–Si–O terms",
        "Si–C and C–H terms",
        "Bond-bond terms",
        "Bond-angle terms",
        "Angle-angle terms",
        "End-bond-torsion terms",
        "Middle-bond-torsion terms",
        "Angle-torsion terms",
        "Angle-angle-torsion terms",
        "Torsion phases and AAT conventions",
        "Charges and nonbonded terms",
    ]
    rows = {g: [] for g in groups}
    for r in reports["source-evidence-table.json"]["rows"]:
        f = r["interaction_family"]
        sites = "-".join(r["ordered_sites"])
        if (
            r["source_file_sha256"].startswith("8b0e8b74303b")
            and r["source_force_field_family"] == "COMPASS"
        ):
            lookup = {"6": "C", "13": "H", "39": "O", "44": "Si"}
            sites = "-".join(lookup.get(t, t) for t in r["ordered_sites"])
        group = {
            "bond-bond": groups[4],
            "bond-bond_1_3": groups[4],
            "bond-angle": groups[5],
            "angle-angle": groups[6],
            "class2_AA_triplet": groups[6],
            "end_bond-torsion_3": groups[7],
            "middle_bond-torsion_3": groups[8],
            "angle-torsion_3": groups[9],
            "angle-angle-torsion_1": groups[10],
        }.get(f)
        if group is None:
            if "torsion" in f or "improper" in f or "out_of_plane" in f:
                group = groups[11]
            elif "bond" in f and f != "bond_increment" and f != "nonbond(9-6)":
                group = groups[0] if "Si" in sites and "O" in sites else groups[3]
            elif "angle" in f:
                group = (
                    groups[1]
                    if sites.lower() in ("si-o-si", "si4-o2-si4", "si-ob-si")
                    else groups[2]
                    if sites.lower() in ("o-si-o", "o2-si4-o2")
                    else groups[3]
                )
            else:
                group = groups[12]
        loc = r["location"]
        where = (
            f"PDF p{loc['pdf_page']}, {loc['table']}, row {loc['row_label']}"
            if "pdf_page" in loc
            else f"{loc['archive_member']}:{loc['line']}"
        )
        values = "; ".join(
            f"{k}={v if v is not None else 'blank (unresolved)'} {u}"
            for k, v, u in zip(r["fields"], r["raw_value"], r["raw_unit"], strict=True)
        )
        rows[group].append(
            f"| {r['record_id']} | {r['source_force_field_family']}; {r['source_model_resolution']} | {where} | {sites} | {values} | {', '.join(r['mapping_statuses'])} |"
        )
    out = [
        "# J71 extracted source evidence",
        "",
        "All values below are raw printed cells. Normalized fields, exact SHA identities, equations, analytic chain-rule forces, dependency blockers and target instance IDs are in the JSON records. A blank is not zero. Archive numeric IDs retain their topology-bound source meaning. No row is accepted for automatic assignment.",
        "",
    ]
    for g in groups:
        out += (
            [
                "## " + g,
                "",
                "| Record | Source family/resolution | Location | Raw sites | Raw fields | J70 compatibility |",
                "|---|---|---|---|---|---|",
            ]
            + rows[g]
            + [""]
        )
    return "\n".join(out).rstrip() + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    result = reports(args.repo.resolve(), args.root.resolve())
    for name, data in result.items():
        write(args.root / name, data)
    (args.root / "source-evidence-table.md").write_text(human_table(result))
    print(
        "J71 source_conflict; incomplete; no converter/QM/fitting/MD/validation execution"
    )


if __name__ == "__main__":
    main()
