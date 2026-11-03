"""Exact PCFF source-row adjudication, without scientific execution.

Capture reads only hash-verified FRC text at the exact J75 source root. Replay
needs only frozen metadata. Source absence and convention absence are separate.
"""

import argparse
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from scripts import phase_4j75_source_path as j75
from scripts.phase_4j73_source_intake import (
    IntakeError,
    digest,
    exact,
    pin,
    read,
    require,
    signed,
    verify_signed,
    write,
)

BASE = "31735a5f3236b6231cf507e57bcab0db7f87606e"
MAIN = "c118ddc06a7bd2532bea93c5e8d6ce17727726ef"
BRANCH = "codex/phase-4j76-pcff-exact-row-resolution"
ROOT = Path("docs/evidence/phase_4j76")
FREEZE_SHA = "bc7f1bbe5cc4f77e58a1dd43255c456a2e2fd3e0530db2f7c9f6d85aea381c69"
PARENTS = (
    "docs/evidence/phase_4j73/source-intake.json",
    *[
        str(j75.ROOT / n)
        for n in (
            "receipt.json",
            "manifest.json",
            "frozen-inputs.json",
            "source-root-manifest.json",
            "source-file-inventory.json",
            "extracted-parameter-records.json",
            "source-parsing-facts.json",
            "row-resolution-matrix.json",
            "source-conflicts.json",
        )
    ],
    "docs/evidence/phase_4j72/source-conflict-matrix.json",
    "docs/evidence/phase_4j72/row-adjudication.json",
    "docs/evidence/phase_4j72/receipt.json",
    "docs/evidence/phase_4j72_c1/quarantine-manifest.json",
    "docs/evidence/phase_4j72_c1/receipt.json",
    *j75.REPLAY_FILES,
)
REPLAY_FILES = ("scripts/phase_4j76_source_resolution.py", *j75.REPLAY_FILES)
STATES = (
    "exact_native_row_resolved",
    "exact_native_row_convention_blocked",
    "compass_only",
    "external_only",
    "analog_only",
    "source_missing",
    "conflicting_source",
    "unresolved",
)
FIELDS = (
    "functional_form",
    "units_and_conversion",
    "parameter_ordering",
    "sign_and_phase_multiplicity",
    "orientation",
    "equilibrium_references",
    "Class_II_dependencies",
    "atom_typing_and_environment",
    "charges",
    "LJ_form_and_mixing",
    "one_two_rules",
    "one_three_rules",
    "one_four_rules",
    "exclusions_and_scaling",
    "long_range_electrostatics",
    "cutoff_and_tail",
)
TARGET_FAMILIES = {
    "angle-angle",
    "bond-bond_1_3",
    "end_bond-torsion_3",
    "middle_bond-torsion_3",
    "angle-torsion_3",
    "angle-angle-torsion_1",
}
SCOPE = {
    "J76_new_QM_attempts": 0,
    "J76_MD_steps": 0,
    "J76_LAMMPS_executions": 0,
    "J76_OpenMM_executions": 0,
    "J76_msi2lmp_executions": 0,
    "converter_attempts": 0,
    "fitting_executions": 0,
    "parameter_generation_executions": 0,
    "evaluator_executions": 0,
    "reserved_validation_executions": 0,
    "accepted_coefficients": 0,
    "executable_profile_created": False,
    "automatic_assignment_enabled": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "chemical_transfer": "not_established",
    "scientific_validation": "not_established",
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
    "author_confirmed_reported": 28,
    "author_confirmed_audited": 0,
    "original_author_artifact_status": "not_yet_supplied",
    "source_adjudication_excluded_quarantine": True,
}


def scope():
    j75.scope()


def graph_from_requests(requests):
    sites = {}
    for t in requests:
        for site in t["chemical_environment"]:
            if site["id"] in sites:
                exact(site, sites[site["id"]], "inconsistent_stable_graph_site")
            sites[site["id"]] = site
    require(len(sites) == 27, "HMDS_27_site_coverage")
    for site in sites.values():
        for n in site["neighbors"]:
            require(
                n["id"] in sites and n["type"] == sites[n["id"]]["type"],
                "graph_neighbor_type",
            )
            require(
                any(x["id"] == site["id"] for x in sites[n["id"]]["neighbors"]),
                "graph_edge_not_symmetric",
            )
    require(sum(len(s["neighbors"]) for s in sites.values()) == 52, "HMDS_26_bonds")
    exact(
        dict(Counter(s["element"] for s in sites.values())),
        {"C": 6, "H": 18, "O": 1, "Si": 2},
        "HMDS_elements",
    )
    return [sites[i] for i in sorted(sites)]


def parse_frc(text, source, facts):
    """Independent literal row index, including rows outside J75's Si filter.

    All target-family numeric rows are indexed. Unrelated coefficient values
    are hash-bound without copying the force-field library into evidence.
    J75 retained values and identities are linked by exact line and hash.
    """
    section, namespace = "", ""
    headers, types, equivalences, rows, notes = {}, {}, {}, [], []
    retained = defaultdict(list)
    for f in facts:
        if f["source_path"] == source["relative_path"]:
            retained[f["location"].get("line")].append(f["record_id"])
    for n, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        w = s.split()
        if not w:
            continue
        if s.startswith("#"):
            section, namespace = w[0][1:], w[1] if len(w) > 1 else "not_stated"
            headers.setdefault(
                section + ":" + namespace, {"section_line": n, "statements": []}
            )
        if s.startswith(("#", "!", ">")):
            if s.startswith(">") or (
                s.startswith("!") and any(k in s for k in ("Ver Ref", "LEFT", "RIGHT"))
            ):
                headers[section + ":" + namespace]["statements"].append(
                    {"line": n, "text": s}
                )
            if re.search(
                r"unit|kcal|angstrom|radian|degree|coulomb|exclu|1[-–]2|1[-–]3|1[-–]4|scaling|long.range|cutoff|cut off|tail|ewald",
                s,
                re.IGNORECASE,
            ):
                notes.append({"line": n, "text": s[:320]})
            continue
        if not re.fullmatch(r"\d+(?:\.\d+)?", w[0]) or len(w) < 3:
            continue
        if section == "atom_types" and len(w) >= 6:
            require(w[2] not in types, "duplicate_or_conflicting_source_typing")
            types[w[2]] = {
                "line": n,
                "version": w[0],
                "reference": w[1],
                "mass_raw": w[3],
                "element": w[4],
                "valence_raw": w[5],
                "description": " ".join(w[6:]),
            }
        elif section == "equivalence" and "auto" not in namespace and len(w) >= 8:
            mapping = dict(
                zip(
                    ("nonbonded", "bond", "angle", "torsion", "improper"),
                    w[3:8],
                    strict=True,
                )
            )
            declaration = {
                "line": n,
                "version": w[0],
                "reference": w[1],
                "namespace": namespace,
                "mapping": mapping,
            }
            if w[2] in equivalences:
                prior = equivalences[w[2]]
                require(
                    all(prior[k] == v for k, v in mapping.items())
                    and prior["namespace"] == namespace,
                    "duplicate_or_conflicting_source_equivalence",
                )
                prior["source_declarations"].append(declaration)
            else:
                equivalences[w[2]] = {
                    "line": n,
                    "namespace": namespace,
                    **mapping,
                    "source_declarations": [declaration],
                }
        elif section in TARGET_FAMILIES:
            require(len(w) >= 7, "malformed_Class_II_source_row")
            ordered, values = w[2:6], w[6:]
            require(
                all(
                    re.fullmatch(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", v)
                    for v in values
                ),
                "nonliteral_source_coefficient",
            )
            rows.append(
                signed(
                    {
                        "source_path": source["relative_path"],
                        "source_sha256": source["sha256"],
                        "source_family": source["source_family"],
                        "section": section,
                        "namespace": namespace,
                        "line": n,
                        "row_version": w[0],
                        "source_reference": w[1],
                        "ordered_types": ordered,
                        "coefficient_count": len(values),
                        "coefficient_lexemes_sha256": digest(values),
                        "source_line_sha256": j75.old.sha_bytes(line.encode()),
                        "retained_J75_fact_ids": retained.get(n, []),
                        "wildcard": any("*" in t for t in ordered),
                        "automatic_namespace": "auto" in namespace,
                    }
                )
            )
    return signed(
        {
            "source": source,
            "headers": headers,
            "atom_types": types,
            "ordinary_equivalences": equivalences,
            "target_family_row_index": rows,
            "convention_keyword_statements": notes,
            "all_source_lines_read": len(text.splitlines()),
            "source_not_executed": True,
        }
    )


def capture(repo, root, source):
    scope()
    j75.check_declared(source)
    # Verifies all 27 inventory entries, including source-root identity.
    j75.originals(repo / j75.ROOT)
    require(not root.exists(), "do_not_overwrite_J76_capture")
    require(
        not root.resolve().is_relative_to(Path(source)), "output_inside_source_root"
    )
    parents = {n: pin(repo / n) for n in PARENTS}
    request = read(repo / PARENTS[0])
    facts = read(repo / j75.ROOT / "extracted-parameter-records.json")["records"]
    sources = read(repo / j75.ROOT / "source-file-inventory.json")["files"]
    parsed = []
    for source_file in sources:
        if not source_file["relative_path"].endswith(".frc"):
            continue
        p = Path(source) / source_file["relative_path"]
        require(
            not p.is_symlink() and pin(p) == source_file["sha256"],
            "authoritative_FRC_hash_or_alias_mutation",
        )
        parsed.append(parse_frc(p.read_text(), source_file, facts))
    j75.originals(repo / j75.ROOT)
    root.mkdir(parents=True)
    write(
        root / "source-capture.json",
        signed(
            {
                "schema": "island.j76-source-capture.v1",
                "exact_source_root": source,
                "J75_root_manifest_sha256": pin(
                    repo / j75.ROOT / "source-root-manifest.json"
                ),
                "J75_inventory_files": 27,
                "J75_retained_source_facts": len(facts),
                "new_PDF_extractions": 0,
                "FRC_sources": parsed,
                "graph_identity": request["graph_identity"],
                "stable_graph": graph_from_requests(request["requests"]),
                "scope": SCOPE,
            }
        ),
    )
    tree = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    history = {}
    for line in tree.splitlines():
        meta, n = line.split("\t", 1)
        p = repo / n
        history[n] = (
            j75.old.sha_bytes(str(p.readlink()).encode())
            if meta.startswith("120000")
            else pin(p)
        )
    write(root / "historical-files.json", history)
    write(
        root / "frozen-inputs.json",
        signed(
            {
                "schema": "island.j76-input-freeze.v1",
                "base_commit": BASE,
                "fetched_origin_main": MAIN,
                "parents": parents,
                "captured_files": {
                    n: pin(root / n)
                    for n in ("source-capture.json", "historical-files.json")
                },
                "scope": SCOPE,
            }
        ),
    )
    print(
        json.dumps(
            {
                "freeze_sha256": pin(root / "frozen-inputs.json"),
                "history": len(history),
                "FRC_sources": len(parsed),
            }
        )
    )


def inputs(repo, root):
    scope()
    require(
        pin(root / "frozen-inputs.json") == FREEZE_SHA,
        "fabricated_or_resigned_input_freeze",
    )
    freeze = read(root / "frozen-inputs.json")
    verify_signed(freeze)
    require(
        freeze["base_commit"] == BASE and freeze["fetched_origin_main"] == MAIN,
        "base_identity",
    )
    exact(freeze["scope"], SCOPE, "resigned_scientific_scope")
    require(
        set(freeze["parents"]) == set(PARENTS)
        and set(freeze["captured_files"])
        == {"source-capture.json", "historical-files.json"},
        "explicit_input_inventory",
    )
    for n, sha in freeze["parents"].items():
        require(pin(repo / n) == sha, "historical_input_mutation:" + n)
    for n, sha in freeze["captured_files"].items():
        require(pin(root / n) == sha, "source_capture_mutation:" + n)
    capture = read(root / "source-capture.json")
    verify_signed(capture)
    j75.check_declared(capture["exact_source_root"])
    return capture


def query_source(source, target):
    """Only source-declared ordinary equivalence; never wildcard or analogy."""
    role = "angle" if target["family"] == "angle-angle" else "torsion"
    equiv = source["ordinary_equivalences"]
    mapped = []
    equivalence_rows = []
    complete = True
    for types in target["symmetry_equivalent_types"]:
        values = []
        for t in types:
            if t not in equiv:
                complete = False
                break
            values.append(equiv[t][role])
            equivalence_rows.append(
                {
                    "type": t,
                    "role": role,
                    "line": equiv[t]["line"],
                    "maps_to": equiv[t][role],
                    "source_declarations": equiv[t]["source_declarations"],
                }
            )
        if len(values) == len(types):
            mapped.append(values)
    exact_types = [*target["symmetry_equivalent_types"], *mapped]
    matches = [
        r
        for r in source["target_family_row_index"]
        if r["section"] == target["family"]
        and r["ordered_types"] in exact_types
        and not r["wildcard"]
        and not r["automatic_namespace"]
    ]
    return {
        "source_path": source["source"]["relative_path"],
        "source_sha256": source["source"]["sha256"],
        "source_family": source["source"]["source_family"],
        "PCFF_upstream_identity_complete": source["source"][
            "native_source_identity_complete"
        ],
        "role": role,
        "source_equivalence_rows": sorted(
            {digest(e): e for e in equivalence_rows}.values(),
            key=lambda e: (e["type"], e["line"]),
        ),
        "equivalence_mapping_complete": complete,
        "source_order_queries": exact_types,
        "indexed_family_rows": sum(
            r["section"] == target["family"] for r in source["target_family_row_index"]
        ),
        "exact_numeric_row_matches": matches,
        "wildcard_or_analog_used": False,
    }


def validate_candidates(candidates):
    require(
        len(candidates) == len({r["id"] for r in candidates}),
        "duplicate_candidate_identity",
    )
    for c in candidates:
        require(
            c.get("selection_basis")
            not in {
                "residual_improvement",
                "average",
                "fitted_value",
                "generic_default",
            },
            "unqualified_parameter_selection",
        )
        require(
            c.get("value_origin")
            not in {
                "zero_fill",
                "wildcard_fallback",
                "guessed_mixing_rule",
                "parameter_generation",
            },
            "missing_value_substitution",
        )


def adjudicate(candidates, historical_status, required_contract=None):
    """Synthetic-testable gate; real candidates retain explicit failed checks."""
    validate_candidates(candidates)
    exact_rows = [c for c in candidates if c["exact_types"] and c["exact_environment"]]
    native = [
        c
        for c in exact_rows
        if c["source_family"] == "native_PCFF" and c["authority_verified"]
    ]
    if native:
        # Different contracts are not compared numerically or averaged.
        contracts = {digest(c["contract"]) for c in native}
        if len(contracts) != 1 or len({digest(c["raw_values"]) for c in native}) != 1:
            return "conflicting_source"
        closed = [
            c
            for c in native
            if required_contract is not None
            and set(required_contract) == set(FIELDS)
            and c["contract"] == required_contract
            and set(c["closure"]) == set(FIELDS)
            and all(
                v.get("status")
                in {"explicit_source_supported", "explicit_source_not_applicable"}
                and v.get("source_references")
                and v.get("value") == required_contract[k]
                for k, v in c["closure"].items()
            )
        ]
        return (
            "exact_native_row_resolved"
            if closed
            else "exact_native_row_convention_blocked"
        )
    if any(
        c["source_family"] == "COMPASS" and c["authority_verified"] for c in exact_rows
    ):
        return "compass_only"
    if any(
        c["source_family"] == "external_PDMS" and c["authority_verified"]
        for c in exact_rows
    ):
        return "external_only"
    return {
        "analog_only": "analog_only",
        "convention_blocked": "unresolved",
        "source_missing": "source_missing",
    }[historical_status]


def fact_summary(f):
    return {
        "id": f["record_id"],
        "J75_fact_identity": f["identity"],
        "source_file": f["source_path"],
        "source_sha256": f["source_sha256"],
        "container_sha256": f["container_sha256"],
        "source_family": f["source_family"],
        "source_model_identity": f["source_model_identity"],
        "source_authority": f["source_authority"],
        "document_identity": f["citation"],
        "source_release": f["version"],
        "location": f["location"],
        "ordered_atom_types": f["ordered_atom_types"],
        "interaction_family": f["interaction_family"],
        "raw_values": f["coefficients_raw"],
        "raw_units": f["units"],
        "unit_conversion": {"status": "not_applied_unclosed_contract", "factor": None},
        "functional_form": f["analytic_equation"],
        "force_equation_metadata": f["analytic_force_equation"],
        "equilibrium_value": f["equilibrium_value"],
        "phase_multiplicity": f["phase_or_multiplicity"],
        "parameter_ordering": f.get(
            "source_field_labels", "see original ordered J75 record"
        ),
        "sign_convention": f["sign_convention"],
        "orientation": f["orientation_order_convention"],
        "dependencies": f["equilibrium_dependencies"],
        "exclusions": f["applicable_exclusions"],
        "charge_origin": f["charge_origin"],
        "LJ_form_and_mixing": f["LJ_equation_and_mixing"],
        "exact_graph_applicability": f["exact_graph_applicability"],
        "original_record_reference": str(j75.ROOT / "extracted-parameter-records.json")
        + "#/records/"
        + f["record_id"],
    }


def closure_for(target, queries):
    native = next(q for q in queries if q["source_path"] == "pcff.frc")
    refs = [
        {
            "source_file": native["source_path"],
            "sha256": native["source_sha256"],
            "section": target["family"],
            "exact_numeric_rows": len(native["exact_numeric_row_matches"]),
        }
    ]
    return {
        field: {
            "status": "unresolved_for_exact_native_row",
            "source_references": refs,
            "value": None,
            "reason": "no exact numeric native PCFF row; generic section headers and historical model conventions cannot supply missing row authority",
            "historical_context_reference": "required_evidence/equation_and_units"
            if field
            in {
                "functional_form",
                "units_and_conversion",
                "parameter_ordering",
                "sign_and_phase_multiplicity",
                "orientation",
            }
            else "required_evidence/equilibria"
            if field in {"equilibrium_references", "Class_II_dependencies"}
            else "global_conventions",
        }
        for field in FIELDS
    }


def global_conventions(native, graph):
    src = native["source"]
    headers = native["headers"]
    support = lambda sections: [
        {"source_file": "pcff.frc", "sha256": src["sha256"], "section": s, **h}
        for s, h in headers.items()
        if s.split(":")[0] in sections
    ]
    typing = []
    for site in graph:
        row = native["atom_types"].get(site["type"])
        typing.append(
            {
                "site_id": site["id"],
                "type": site["type"],
                "source_typing_row": row,
                "element_and_valence_match": bool(
                    row
                    and row["element"] == site["element"]
                    and row["valence_raw"] == str(len(site["neighbors"]))
                ),
                "does_not_supply_missing_interaction": True,
            }
        )
    return {
        "units": {
            "status": "native_file_unit_authority_not_declared; historical units retained only as context",
            "keyword_statements": native["convention_keyword_statements"],
        },
        "bond_function": {
            "status": "literal_source_form; not complete missing-row contract",
            "evidence": support({"quartic_bond"}),
        },
        "angle_function": {
            "status": "literal_source_form; angular-unit authority separate",
            "evidence": support({"quartic_angle"}),
        },
        "torsion_function": {
            "status": "literal_source_phase_fields; no cross-term phase inferred",
            "evidence": support({"torsion_3"}),
        },
        "Class_II_cross_terms": {
            "status": "family equations present; exact requested numeric rows absent",
            "evidence": support(TARGET_FAMILIES | {"bond-bond", "bond-angle"}),
            "AAT": "Literal Phi-Phi1(0) header is not an authorized cosine replacement; numeric missing rows and Phi/AAT authority remain separate blockers",
        },
        "atom_typing": {
            "status": "source declarations and stable-site element/valence comparison; no automatic assignment",
            "sites": typing,
        },
        "charges": {
            "status": "unresolved_native_exact_graph_charge_origin",
            "evidence": support({"bond_increments"}),
            "restriction": "historical fixed charges and external observed charges are not author-primary charge rules",
        },
        "Lennard_Jones_form": {
            "status": "literal 9-6 source equation",
            "evidence": support({"nonbond(9-6)"}),
        },
        "mixing_rule": {
            "status": "literal source sixth-power mixing statements; no guessed rules",
            "evidence": support({"nonbond(9-6)"}),
        },
        **{
            name: {
                "status": "unresolved_native_source_contract",
                "source_evidence": [],
                "restriction": "not supplied by native FRC; external implementation settings cannot close PCFF convention authority",
            }
            for name in (
                "one_two_rules",
                "one_three_rules",
                "one_four_rules",
                "exclusions_and_scaling",
                "long_range_electrostatics",
                "cutoff",
                "tail",
            )
        },
    }


def reconstruct(repo, root):
    capture = inputs(repo, root)
    request = read(repo / PARENTS[0])
    previous = read(repo / j75.ROOT / "row-resolution-matrix.json")
    facts = read(repo / j75.ROOT / "extracted-parameter-records.json")["records"]
    by_fact = {f["record_id"]: f for f in facts}
    require(len(by_fact) == len(facts) == 1012, "duplicate_or_incomplete_J75_facts")
    by_previous = {t["id"]: t for t in previous["instances"]}
    require(
        len(by_previous)
        == len(request["requests"])
        == len({t["id"] for t in request["requests"]})
        == 378,
        "incomplete_378_requests",
    )
    native = next(
        s for s in capture["FRC_sources"] if s["source"]["relative_path"] == "pcff.frc"
    )
    globals_ = global_conventions(native, capture["stable_graph"])
    rows = []
    for t in request["requests"]:
        prior = by_previous[t["id"]]
        queries = [query_source(s, t) for s in capture["FRC_sources"]]
        # An unexpected exact discovery must receive a new evidence audit, not
        # be promoted by an existing generic equation or a local file name.
        require(
            not any(q["exact_numeric_row_matches"] for q in queries),
            "unexpected_exact_source_row_requires_explicit_contract_audit",
        )
        candidates = []
        for old_candidate in prior["candidate_source_records"]:
            f = by_fact[old_candidate["record_id"]]
            family = (
                "native_PCFF"
                if f["source_path"] == "pcff.frc"
                else "COMPASS"
                if f["source_family"] == j75.FAMILIES[1]
                else "external_PDMS"
                if f["source_family"]
                in {j75.FAMILIES[2], j75.FAMILIES[3], j75.FAMILIES[4]}
                else "analog"
            )
            candidates.append(
                {
                    "id": f["record_id"],
                    "source_family": family,
                    "exact_types": f["ordered_atom_types"]
                    in t["symmetry_equivalent_types"],
                    "exact_environment": f["exact_graph_applicability"],
                    "authority_verified": False,
                    "contract": {
                        "family": f["interaction_family"],
                        "equation": f["analytic_equation"],
                        "units": f["units"],
                    },
                    "raw_values": f["coefficients_raw"],
                    "closure": {},
                    "checks": {
                        k: old_candidate[k]
                        for k in (
                            "exact_type_match",
                            "exact_graph_match",
                            "function_form_match",
                            "unit_match",
                            "phase_multiplicity_match",
                            "orientation_order_match",
                            "equilibrium_dependency_match",
                            "exclusion_match",
                            "charge_convention_match",
                            "LJ_and_mixing_match",
                        )
                    },
                    "catalog_reference": "candidate-source-catalog.json#/records/"
                    + f["record_id"],
                    "eligibility": False,
                }
            )
        state = adjudicate(candidates, prior["resolution_status"])
        rows.append(
            signed(
                {
                    "id": t["id"],
                    "J70_identity": t["J70_identity"],
                    "J72_identity": t["J72_identity"],
                    "J73_identity": t["identity"],
                    "J75_identity": prior["identity"],
                    "J75_status": prior["resolution_status"],
                    "status": state,
                    "ordered_atom_types": t["ordered_types"],
                    "stable_site_ids": t["stable_sites"],
                    "chemical_environment": t["chemical_environment"],
                    "graph_identity": request["graph_identity"],
                    "interaction_family": t["family"],
                    "symmetry_equivalent_types": t["symmetry_equivalent_types"],
                    "source_queries": queries,
                    "candidate_records": candidates,
                    "exact_native_row_present": False,
                    "closed_convention_fields": 0,
                    "convention_closure": closure_for(t, queries),
                    "historical_equation_and_units_context": t["required_evidence"][
                        "equation_and_units"
                    ],
                    "historical_equilibrium_context": t["required_evidence"][
                        "equilibria"
                    ],
                    "remaining_source_evidence": t["required_evidence"],
                    "blockers": [
                        "exact native numeric row absent",
                        "whole scalar/unit/phase/orientation/dependency authority absent for this row",
                        "native graph-wide charge/nonbonded/exclusion/runtime conventions not closed",
                    ],
                    "historical_conflict_categories": t[
                        "inherited_conflict_categories"
                    ],
                    "primary_artifact_supplied": False,
                    "author_confirmed_audited": False,
                }
            )
        )
    counts = {state: sum(r["status"] == state for r in rows) for state in STATES}
    exact(
        counts,
        {
            **dict.fromkeys(STATES, 0),
            "analog_only": 294,
            "source_missing": 24,
            "unresolved": 60,
        },
        "real_source_statuses_changed",
    )
    family_counts = dict(Counter(r["interaction_family"] for r in rows))
    reports = {}
    for name, selector in (
        ("native-pcff", lambda f: f["source_path"] == "pcff.frc"),
        ("compass", lambda f: f["source_family"] == j75.FAMILIES[1]),
        (
            "external-pdms",
            lambda f: (
                f["source_family"]
                in {j75.FAMILIES[2], j75.FAMILIES[3], j75.FAMILIES[4]}
            ),
        ),
    ):
        selected = [f for f in facts if selector(f)]
        reports[name + "-coverage.json"] = signed(
            {
                "schema": "island.j76-coverage.v1",
                "source_family": name,
                "retained_facts": len(selected),
                "fact_ids": [f["record_id"] for f in selected],
                "exact_native_requests_resolved": 0,
                "exact_same_graph_source_rows_supplied": 0,
                "unresolved_native_requests": 378,
                "unresolved_row_ids": [r["id"] for r in rows],
                "missing_Class_II_instances": family_counts,
                "scope_of_counts": "exact J73/J75 native PCFF requests; no statement of completeness of a different model",
                "COMPASS_or_external_values_promoted_to_PCFF": False,
                "CFF91_provenance_policy": "cff91.frc indexed separately, never authorized as native PCFF solely by name",
                "scope": SCOPE,
            }
        )
    requirements = [
        {
            "id": r["id"],
            "request_identity": r["identity"],
            "status": r["status"],
            "ordered_types": r["ordered_atom_types"],
            "stable_sites": r["stable_site_ids"],
            "native_source_row_needed": {
                "family": r["interaction_family"],
                "numeric_values_supplied": False,
                "complete_source_equation_and_units_needed": True,
            },
            "required_artifact": "original author/publisher same-contract row/table/source file, immutable SHA-256, release, license, exact location and complete site/convention authority",
            "required_fields": list(FIELDS),
            "historical_required_evidence": r["remaining_source_evidence"],
            "prohibited_substitutions": [
                "analog",
                "COMPASS overlay",
                "external PDMS overlay",
                "zero fill",
                "default mixing",
                "fitting",
                "averaging",
            ],
        }
        for r in sorted(rows, key=lambda r: r["id"])
    ]
    original_conflicts = read(
        repo / "docs/evidence/phase_4j72/source-conflict-matrix.json"
    )
    return {
        "row-resolution-matrix.json": signed(
            {
                "schema": "island.j76-row-resolution.v1",
                "instances": rows,
                "counts": counts,
                "family_counts": family_counts,
                "scope": SCOPE,
            }
        ),
        **reports,
        "candidate-source-catalog.json": signed(
            {
                "schema": "island.j76-candidate-catalog.v1",
                "records": {f["record_id"]: fact_summary(f) for f in facts},
                "source_assets_redistributed": False,
                "source_value_generation": False,
                "scope": SCOPE,
            }
        ),
        "convention-closure-matrix.json": signed(
            {
                "schema": "island.j76-convention-closure.v1",
                "global_conventions": globals_,
                "requests": [
                    {
                        "id": r["id"],
                        "status": r["status"],
                        "native_row_present": False,
                        "closure": r["convention_closure"],
                    }
                    for r in rows
                ],
                "complete_native_rows": 0,
                "complete_bundle": False,
                "scope": SCOPE,
            }
        ),
        "source-conflicts.json": signed(
            {
                "schema": "island.j76-conflicts.v1",
                "original_J72_matrix": original_conflicts,
                "original_J72_sha256": pin(
                    repo / "docs/evidence/phase_4j72/source-conflict-matrix.json"
                ),
                "equal_duplicate_equivalences": 16,
                "unequal_unresolved_COMPASS_assignments": 2,
                "new_resolutions": 0,
                "unequal_source_values_preserved": True,
                "native_row_conflicts": 0,
                "interpretation": "Historical unequal COMPASS alternatives remain analog conflicts; no exact native candidate exists to adjudicate. J76 never compares unlike equations numerically.",
                "scope": SCOPE,
            }
        ),
        "remaining-evidence-requests.json": signed(
            {
                "schema": "island.j76-evidence-requests.v1",
                "requests": requirements,
                "shared_native_convention_blockers": globals_,
                "author_attestations_are_not_artifacts": True,
                "scope": SCOPE,
            }
        ),
        "capability-report.json": signed(
            {
                "schema": "island.j76-capability.v1",
                "decision": "parameter_incomplete",
                "row_counts": counts,
                "converter_status": "not_run_prohibited",
                "static_closure_status": "not_run_prohibited",
                "ordinary_pcff_behavior_changed": False,
                "scope": SCOPE,
            }
        ),
    }


def replay(repo, root, build=False):
    documents = reconstruct(repo, root)
    for name, expected in documents.items():
        if build:
            require(not (root / name).exists(), "do_not_overwrite_adjudication_attempt")
            write(root / name, expected)
        exact(read(root / name), expected, "resigned_scientific_evidence:" + name)
    return {
        "decision": "parameter_incomplete",
        "status_counts": documents["row-resolution-matrix.json"]["counts"],
        "scope": SCOPE,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("capture", "build", "replay"))
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--directory", type=Path, default=ROOT)
    parser.add_argument("--source", default=j75.AUTHORITATIVE)
    args = parser.parse_args()
    try:
        scope()
        j75.check_declared(args.source)
        result = (
            capture(args.repo, args.directory, args.source)
            if args.mode == "capture"
            else replay(args.repo, args.directory, args.mode == "build")
        )
        print(json.dumps(result, sort_keys=True))
        return 0
    except (IntakeError, OSError, KeyError, TypeError, RuntimeError) as error:
        print(json.dumps({"status": "blocked", "reason": str(error), "scope": SCOPE}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
