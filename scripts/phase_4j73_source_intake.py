"""Deterministic requests for primary evidence; never parameters or execution."""

import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

BASE = "73eb64e02ab5d970e02f444e1b41c8183c49034d"
BRANCH = "codex/phase-4j73-authoritative-pcff-source-intake"
ROOT = Path("docs/evidence/phase_4j73")
FREEZE_SHA = "e0ee9909fa9d333fed4f14a75015ba3eccf8f149f4a85cb82dc71585b33b3801"
INPUTS = (
    "docs/evidence/phase_4j69/source-dependency-graph.json",
    "docs/evidence/phase_4j69/source-row-inventory.json",
    "docs/evidence/phase_4j69/pair-source-dependencies.json",
    "docs/evidence/phase_4j70/completeness-matrix.json",
    "docs/evidence/phase_4j71/parameter-records.json",
    "docs/evidence/phase_4j72/authority-hierarchy.json",
    "docs/evidence/phase_4j72/row-adjudication.json",
    "docs/evidence/phase_4j72/source-conflict-matrix.json",
    "docs/evidence/phase_4j72/unresolved-evidence-requirements.json",
    "docs/evidence/phase_4j72/receipt.json",
    "docs/evidence/phase_4j72_c1/quarantine-manifest.json",
    "docs/evidence/phase_4j72_c1/contamination-exclusion.json",
    "docs/evidence/phase_4j72_c1/receipt.json",
)
COUNTS = {"analog_unqualified": 294, "convention_unresolved": 60, "source_missing": 24}
ROW_STATES = {
    "author_confirmed_reported": "historical_28_queries_not_mapped_to_this_instance",
    "author_confirmed_audited": False,
    "primary_artifact_supplied": False,
    "exact_row_resolved": False,
    "convention_resolved": False,
    "executable_bundle_qualified": False,
}
SCOPE = {
    "authorized_new_QM": 0,
    "authorized_MD_steps": 0,
    "J73_new_QM_attempts": 0,
    "J73_MD_steps": 0,
    "J73_OpenMM_executions": 0,
    "J73_LAMMPS_executions": 0,
    "J73_AmberTools_executions": 0,
    "J73_converter_attempts": 0,
    "J73_fitting_executions": 0,
    "J73_parameter_generation_executions": 0,
    "J73_reserved_validation_executions": 0,
    "HMDS_new_QM": 0,
    "reserved_validation_QM": 0,
    "reserved_validation_activated": False,
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
    "scientific_evidence_contaminated": False,
    "source_adjudication_excluded_quarantine": True,
    "accepted_coefficients": 0,
    "automatic_assignment": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "parameter_completeness": "incomplete",
    "converter_complete": False,
    "chemical_transfer": "not_established",
    "scientific_validation": "not_established",
    "decision": "source_conflict_unresolved",
    "original_author_artifact_status": "not_yet_supplied",
}
GATE = {
    "phase": "J73_evidence_and_intake_requests_only",
    "profile": None,
    "parameter_bundle": None,
    "converter_result": None,
    "automatic_assignment_path": None,
    "MD_workflow": None,
    "executable_bundle_qualified": False,
    "all_378_exact_rows_required_for_future_qualification": True,
    "resolved_rows_alone_authorize_execution": False,
    "future_requirements": [
        "all exact rows and dependencies independently source-audited under one coherent contract",
        "a separately authorized phase for strict conversion and independent static closure",
        "separate bounded-domain registration review; no J73 execution authority",
    ],
}


class IntakeError(ValueError):
    pass


def require(ok, reason):
    if not ok:
        raise IntakeError(reason)


def scope():
    from scripts.phase_4j73_scope import install

    install()


def pin(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def signed(value):
    return {**value, "identity": digest(value)}


def verify_signed(value):
    require(
        value.get("identity")
        == digest({k: v for k, v in value.items() if k != "identity"}),
        "record_identity_mutation",
    )


def read(path):
    def unique(pairs):
        result = {}
        for k, v in pairs:
            require(k not in result, "duplicate_JSON_key:" + k)
            result[k] = v
        return result

    return json.loads(Path(path).read_text(), object_pairs_hook=unique)


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )


def exact(actual, expected, reason):
    # JSON booleans cannot substitute for integer execution/coverage counts.
    require(type(actual) is type(expected), reason)
    if isinstance(expected, dict):
        require(set(actual) == set(expected), reason)
        for k in expected:
            exact(actual[k], expected[k], reason + ":" + k)
    elif isinstance(expected, list):
        require(len(actual) == len(expected), reason)
        for a, b in zip(actual, expected, strict=True):
            exact(a, b, reason)
    else:
        require(actual == expected, reason)


def parents(repo, root):
    scope()
    require(
        pin(root / "frozen-inputs.json") == FREEZE_SHA,
        "frozen_parent_contract_mutation",
    )
    frozen = read(root / "frozen-inputs.json")
    verify_signed(frozen)
    require(
        frozen["base_commit"] == BASE and frozen["branch"] == BRANCH, "base_identity"
    )
    require(set(frozen["inputs"]) == set(INPUTS), "explicit_parent_inventory")
    result = {}
    for name, sha in frozen["inputs"].items():
        require(pin(repo / name) == sha, "historical_source_hash_mutation:" + name)
        result[name] = read(repo / name)
        verify_signed(result[name])
    return result


def artifact_contract():
    return signed(
        {
            "schema": "island.j73-primary-artifact-contract.v1",
            "purpose": "future source intake requirements; no accepted artifact or numerical row in J73",
            "accepted_origins": [
                "original_author_file",
                "original_publisher_file",
                "exact_native_PCFF_FRC_release",
            ],
            "required_artifact_fields": [
                "original_filename",
                "sha256",
                "citation",
                "license_status",
                "redistribution_status",
                "origin",
                "provenance_chain",
                "source_version",
                "locations",
                "target_ids",
            ],
            "required_citation_fields": [
                "title",
                "authors",
                "publisher_or_journal",
                "year",
                "DOI_or_persistent_identifier",
            ],
            "required_location_fields": [
                "page",
                "section",
                "table",
                "equation",
                "source_row_or_code_lines",
            ],
            "location_rule": "every field must be a nonempty exact reference or documented not_applicable reason; unknown references remain unresolved",
            "required_conventions": [
                "units",
                "scalar_energy_equation",
                "analytic_force_equation",
                "prefactor",
                "phase",
                "multiplicity",
                "orientation_reversal",
                "equilibrium_dependencies",
                "exclusions_12_13_14",
                "charge_origin",
                "LJ_equation_and_mixing",
                "exact_graph_applicability",
            ],
            "raw_values_separate_from_normalized": True,
            "no_inferred_analog_values": True,
            "no_averaging": True,
            "no_zero_filling": True,
            "no_wildcard_fallback": True,
            "no_cross_family_substitution": True,
            "hash_is_not_author_authentication": True,
            "provenance_review_required_before_trust_pin": True,
            "trusted_primary_artifact_pins": {},
            "J73_new_source_record_acceptance": False,
            "new_trust_pins_require_separately_reviewed_contract_revision": True,
        }
    )


def artifact_preflight(record, directory):
    """Reject incomplete/forged submissions; no trust-on-first-use authority.

    This phase has no independently authenticated new artifact pins. Even a
    byte-hashed fabricated file cannot bootstrap its own author authority.
    """
    scope()
    contract = artifact_contract()
    require(set(record) == set(contract["required_artifact_fields"]), "artifact_schema")
    name, sha = record["original_filename"], record["sha256"]
    require(
        type(sha) is str and re.fullmatch("[0-9a-f]{64}", sha), "artifact_hash_required"
    )
    require(
        type(name) is str and Path(name).name == name and name not in ("", ".", ".."),
        "immutable_original_filename_required",
    )
    file = Path(directory) / name
    require(
        file.is_file() and not file.is_symlink() and pin(file) == sha,
        "artifact_bytes_or_filename_mutation",
    )
    require(record["origin"] in contract["accepted_origins"], "analog_not_primary")
    require(
        set(record["citation"]) == set(contract["required_citation_fields"])
        and all(record["citation"].values()),
        "citation_metadata_required",
    )
    require(
        record["locations"]
        and all(
            set(x) == set(contract["required_location_fields"]) and all(x.values())
            for x in record["locations"]
        ),
        "exact_source_locations_required",
    )
    require(
        all(
            record[k]
            for k in (
                "license_status",
                "redistribution_status",
                "provenance_chain",
                "source_version",
                "target_ids",
            )
        ),
        "artifact_provenance_required",
    )
    require(
        sha in contract["trusted_primary_artifact_pins"],
        "untrusted_primary_artifact_no_self_signed_authority",
    )


def reconstruct(repo, root):
    p = parents(repo, root)
    graph = p[INPUTS[0]]
    inventory = p[INPUTS[1]]
    pair = p[INPUTS[2]]
    matrix = p[INPUTS[3]]
    literature = {r["record_id"]: r for r in p[INPUTS[4]]["records"]}
    adjudication = p[INPUTS[6]]
    targets = {r["id"]: r for r in matrix["instances"]}
    old = {r["id"]: r for r in adjudication["instances"]}
    require(
        len(targets) == len(old) == 378
        and set(old) == set(targets) == set(graph["missing_interactions"]),
        "exact_378_instance_inventory",
    )
    require(
        Counter(r["final_J72_status"] for r in old.values()) == COUNTS,
        "unchanged_J72_classifications",
    )
    require(
        adjudication["summary"]["real_duplicate_equivalences_replayed"] == 16
        and adjudication["summary"]["documented_unequal_source_assignment_conflicts"]
        == 2,
        "historical_conflicts_preserved",
    )
    c1 = p[INPUTS[-1]]
    for key in (
        "accepted_coefficients",
        "production_validated",
        "automatic_assignment",
        "simulation_readiness",
        "scientific_evidence_contaminated",
        "source_adjudication_excluded_quarantine",
    ):
        exact(c1[key], SCOPE[key], "C1_inherited_scope:" + key)
    for key, value in (
        ("unintended_external_QM_attempts", 4),
        ("unintended_MD_steps_verified_minimum", 46),
        ("unintended_MD_steps_possible_maximum", 48),
        ("zero_execution_scope_met", False),
    ):
        exact(c1[key], value, "C1_historical_violation:" + key)
    nodes = graph["nodes"]
    typing = inventory["typing"]
    bonds = [n["stable_sites"] for n in nodes if n["family"] == "quartic_bond"]
    require(
        len(typing) == 27
        and len(bonds) == 26
        and graph["graph_identity"] == matrix["graph_identity"],
        "exact_graph_identity",
    )
    adjacent = defaultdict(list)
    for i, j in bonds:
        adjacent[i].append(j)
        adjacent[j].append(i)
    elements = {"c3": "C", "sio": "Si", "osi": "O", "hc": "H"}
    rows, groups = [], defaultdict(list)
    shared = {
        "exclusions_12_13_14": {
            "status": "source_authority_required",
            "required": "original Coulomb and LJ scaling/exclusion rules for each of 1-2, 1-3 and 1-4 pairs; separate special-pair rules and exception handling",
            "historical_source_authority": pair["source_pair_exclusion_authority"],
        },
        "charge_origin": {
            "status": "primary_confirmation_required",
            "historical_context_only": pair["charge_origin"],
            "charge_identity": pair["charge_identity"],
            "required": "original charge or directed bond-increment rows, type/orientation rules, net-charge provenance, no normalization",
        },
        "LJ_equation_and_mixing": {
            "status": "primary_confirmation_required",
            "historical_context_only": {
                "functional_form": pair["pairs"][0]["functional_form"],
                "mixing": pair["pairs"][0]["mixing"],
            },
            "required": "native 9-6 scalar LJ equation, epsilon/rmin versus sigma definition, units, sixth-power or other exact source mixing rule and exclusions; no 12-6 overlay",
        },
        "provenance": "original author/publisher artifact with immutable filename, SHA-256, citation, release/revision, license/redistribution status and exact source location",
    }
    for key in sorted(targets):
        target, inherited = targets[key], old[key]
        require(
            inherited["J70_identity"] == target["identity"]
            and inherited["stable_sites"] == target["stable_sites"]
            and inherited["ordered_types"] == target["ordered_types"],
            "parent_exact_mapping",
        )
        verify_signed(target)
        verify_signed(inherited)
        sites = set(target["stable_sites"])
        context = {}
        for label, families in {
            "bond": {"quartic_bond"},
            "angle": {"quartic_angle"},
            "torsion": {"torsion_3"},
            "Class_II": {
                "angle-angle",
                "bond-bond",
                "bond-angle",
                "bond-bond_1_3",
                "end_bond-torsion_3",
                "middle_bond-torsion_3",
                "angle-torsion_3",
                "angle-angle-torsion_1",
            },
        }.items():
            context[label] = [
                {"id": n["id"], "source_row_ids": n["source_rows"]}
                for n in nodes
                if n["family"] in families and set(n["stable_sites"]) <= sites
            ]
        context["improper"] = {
            "retained_J69_applicability": graph["Wilson"],
            "required_evidence": "authoritative improper applicability and ordering; absence in the J69 declared degree-three inventory is not permission to zero-fill any source family",
        }
        candidates = []
        for candidate in inherited["candidate_records"]:
            original = literature[candidate["record_id"]]
            require(
                candidate["J71_identity"] == original["identity"],
                "candidate_original_identity",
            )
            candidates.append(
                {
                    "record_id": original["record_id"],
                    "J71_identity": original["identity"],
                    "source_file_sha256": original["source_file_sha256"],
                    "citation": original["citation"],
                    "location": original["location"],
                    "source_force_field_family": original["source_force_field_family"],
                    "source_model_resolution": original["source_model_resolution"],
                    "inherited_J72_candidate": candidate,
                    "role": "unresolved_alternative_context_only_not_primary_PCFF_authority",
                }
            )
        group_key = (
            inherited["final_J72_status"],
            target["family"],
            tuple(target["ordered_types"]),
        )
        group_id = "request-group-" + digest(group_key)[:16]
        row = signed(
            {
                "schema": "island.j73-source-intake-request.v1",
                "id": key,
                "J70_identity": target["identity"],
                "J72_identity": inherited["identity"],
                "existing_native_source_lookup": target["source_lookup"],
                "J72_status": inherited["final_J72_status"],
                "request_group": group_id,
                "family": target["family"],
                "stable_sites": target["stable_sites"],
                "ordered_types": target["ordered_types"],
                "symmetry_equivalent_types": target["symmetry_equivalent_types"],
                "chemical_environment": [
                    {
                        "id": i,
                        "element": elements[typing[str(i)]],
                        "type": typing[str(i)],
                        "neighbors": [
                            {"id": j, "type": typing[str(j)]}
                            for j in sorted(adjacent[i])
                        ],
                    }
                    for i in target["stable_sites"]
                ],
                "interaction_context": context,
                "required_evidence": {
                    "force_constants": {
                        "required_fields": target["contract"]["fields"],
                        "value_status": "not_supplied; exact primary numeric cells required",
                        "value": None,
                    },
                    "equilibria": {
                        "historical_dependency_context_only": target["dependencies"],
                        "required": "original equilibrium values, units and source row IDs for every dependency; full Class-II propagation contract; historical references are not new authoritative values",
                    },
                    "equation_and_units": {
                        "historical_unresolved_contract": target["contract"],
                        "required": "complete exact scalar energy and analytic force equations, units, polynomial prefactors, native Phi/AAT, phase, multiplicity, reference order and reversal/orientation; documented not_applicable reasons where appropriate",
                    },
                    "source_location": {
                        "required": "exact original table/page/section/equation/row or source-code function and line range",
                        "authoritative_location": None,
                        "status": "not_supplied",
                    },
                    "applicability": {
                        "graph_identity": matrix["graph_identity"],
                        "resolution": "all_atom_explicit_H",
                        "charge": 0,
                        "multiplicity": 1,
                        "atom_order": matrix["atom_order"],
                        "required": "explicit native PCFF/cff91 typing/equivalence and full neighborhood applicability for these exact stable sites; no element-sequence analogy",
                    },
                    "common_contract_reference": "source-artifact-requests.json#/shared_requirements",
                    "J72_remaining_missing_conventions": inherited[
                        "remaining_missing_conventions"
                    ],
                },
                "unresolved_alternatives": candidates,
                "inherited_conflict_categories": inherited["conflict_categories"],
                "statuses": ROW_STATES.copy(),
                "supplied_artifacts": [],
                "source_records": [],
                "QM_fit_requirement": "not_established; source absence does not authorize fitting",
            }
        )
        rows.append(row)
        groups[group_key].append(row)
    requests = []
    for (status, family, types), members in sorted(groups.items()):
        required = [
            f"Complete original PCFF/FRC release or original author/publisher numerical table for {family}, ordered types {' '.join(types)}, including every listed instance and exact equivalence/reversal rules.",
            "Original scalar-equation and dependency documentation, including bond/angle equilibria and coefficients; source-code record with immutable revision and exact lines if the convention is implemented in code.",
            "Original typing dictionary and graph applicability, charge/bond-increment origin, 1-2/1-3/1-4 LJ/Coulomb rules, LJ equation and mixing contract.",
        ]
        if status == "convention_unresolved":
            required.append(
                "Author-primary native literal Phi/AAT geometry, phase, units and prefactor definition, or exact published transformation/correction with its original supplement/code; cosine COMPASS/research AAT is not native authority."
            )
        if status == "analog_unqualified":
            required.append(
                "An exact same-contract chemical-environment row replacing analogy; no transfer of polysilane/COMPASS/other model values. Preserve all conflicting source alternatives."
            )
        if status == "source_missing":
            required.append(
                "Original missing Si-centered angle-angle coefficient table/source rows; no matching numerical source was available in J72."
            )
        requests.append(
            signed(
                {
                    "id": members[0]["request_group"],
                    "J72_status": status,
                    "family": family,
                    "ordered_types": list(types),
                    "instance_ids": [r["id"] for r in members],
                    "required_original_files_tables_supplements_code": required,
                    "original_primary_filename": None,
                    "original_primary_location": None,
                    "unknown_location_policy": "request exact filename/table/page/revision from author or publisher; do not invent a citation location",
                    "existing_alternative_record_ids_not_native_authority": sorted(
                        {
                            c["record_id"]
                            for r in members
                            for c in r["unresolved_alternatives"]
                        }
                    ),
                }
            )
        )
    statuses = signed(
        {
            "schema": "island.j73-status-separation.v1",
            "author_confirmed_reported": {
                "count": 28,
                "authority": "user_attestation_only",
                "per_instance_primary_mapping": None,
                "historical_C1_receipt_identity": c1["identity"],
            },
            "author_confirmed_audited": {"count": 0, "status": "not_established"},
            "primary_artifact_supplied": {"count": 0, "status": "not_yet_supplied"},
            "exact_row_resolved": {"instances": 0},
            "convention_resolved": {"complete_instance_contracts": 0},
            "executable_bundle_qualified": False,
            "no_status_implies_another_without_its_own_evidence": True,
        }
    )
    package = signed(
        {
            "schema": "island.j73-source-intake.v1",
            "base_commit": BASE,
            "graph_identity": matrix["graph_identity"],
            "atom_order": matrix["atom_order"],
            "family_counts": matrix["family_counts"],
            "unresolved_counts": COUNTS,
            "requests": rows,
            "scope": SCOPE,
            "no_execution_gate": GATE,
        }
    )
    artifact_requests = signed(
        {
            "schema": "island.j73-artifact-requests.v1",
            "incomplete_native_source_context": {
                "filename": "pcff.frc",
                "sha256": matrix["instances"][0]["source_lookup"]["source_sha256"],
                "citation": next(
                    r["citation"]
                    for r in inventory["rows"].values()
                    if r.get("source_file") == "pcff.frc"
                ),
                "role": "existing incomplete distribution candidate; not a supplied new primary artifact",
            },
            "shared_requirements": shared,
            "groups": requests,
            "required_instances": 378,
            "new_primary_artifacts_read": 0,
        }
    )
    return {
        "source-intake.json": package,
        "artifact-contract.json": artifact_contract(),
        "source-artifact-requests.json": artifact_requests,
        "status-separation.json": statuses,
        "no-execution-gate.json": signed(
            {
                "schema": "island.j73-no-execution-gate.v1",
                **GATE,
                "scope": SCOPE,
                "quarantine_manifest_sha256": pin(repo / INPUTS[-3]),
                "J72_C1_receipt_sha256": pin(repo / INPUTS[-1]),
                "historical_equal_duplicate_equivalences": 16,
                "historical_unequal_COMPASS_assignments_unresolved": 2,
            }
        ),
    }


def validate_requests(package, expected):
    verify_signed(package)
    rows = package["requests"]
    require(
        len(rows) == 378
        and len({r["id"] for r in rows}) == 378
        and len({r["identity"] for r in rows}) == 378,
        "incomplete_or_duplicate_row_coverage",
    )
    for row in rows:
        verify_signed(row)
        require(
            row["source_records"] == [] and row["supplied_artifacts"] == [],
            "new_or_fabricated_source_record_forbidden_in_J73",
        )
        exact(row["statuses"], ROW_STATES, "no_status_promotion")
    exact(package["scope"], SCOPE, "scientific_scope_mutation")
    exact(package["no_execution_gate"], GATE, "premature_profile_activation")
    # Regenerate requests from hash-pinned historical inputs, without invoking
    # adjudication, row parsers, parameter builders or a scientific evaluator.
    exact(package, expected, "resigned_request_or_hidden_conflict_mutation")


def render(documents):
    p = documents["source-intake.json"]
    native = documents["source-artifact-requests.json"][
        "incomplete_native_source_context"
    ]
    lines = [
        "# J73 authoritative source-evidence request",
        "",
        "This requests original evidence; it supplies no parameter values and resolves no row.",
        "",
        "294 analog-only, 60 convention-blocked and 24 source-missing instances remain unresolved.",
        "The 28 user-attested confirmations are provenance metadata; zero have audited primary artifacts.",
        "",
        f"Existing incomplete native candidate: `{native['citation']['path']}`, SHA-256 `{native['sha256']}`, LAMMPS revision `{native['citation']['revision']}`. Its retained Biosym/cff91 provenance does not supply the missing rows. Each JSON request binds its original exact native lookup and per-instance search location. Request the complete authoritative original release or an explicitly documented revision with the missing rows; do not substitute a similarly named force-field file.",
        "",
        "## Required for every instance",
        "",
    ]
    for key, value in documents["source-artifact-requests.json"][
        "shared_requirements"
    ].items():
        lines += [
            f"- **{key}**: {value if isinstance(value, str) else value['required']}"
        ]
    lines += [
        "",
        "Each original file needs an immutable filename/SHA-256, author/publisher provenance, citation, version, license/redistribution status and exact page/table/section/equation or code lines. Unknown original filenames/locations are explicitly not supplied; none is invented. Existing alternative citations below identify historical context, not native row authority.",
        "",
        "No analogy, averaging, zero-fill or cross-family substitution is accepted. The native Phi/AAT and equilibrium-dependency definitions, pair exclusions, charge origin and LJ/mixing authority must be supplied together. Passing intake cannot authorize execution in J73.",
        "",
    ]
    by_id = {r["id"]: r for r in p["requests"]}
    for group in documents["source-artifact-requests.json"]["groups"]:
        lines += [
            f"## {group['id']}: {group['family']} / {' '.join(group['ordered_types'])}",
            "",
            f"Inherited status: `{group['J72_status']}`; {len(group['instance_ids'])} exact instances.",
            "",
        ]
        lines += [
            "- " + item
            for item in group["required_original_files_tables_supplements_code"]
        ]
        lines += [
            "",
            "Exact instance IDs (the JSON requests bind ordered sites, full neighborhoods, bond/angle/torsion/improper/Class-II context and source equilibrium dependencies):",
            "",
        ]
        lines += [f"- `{key}`" for key in group["instance_ids"]]
        alternatives = {
            c["record_id"]: c
            for key in group["instance_ids"]
            for c in by_id[key]["unresolved_alternatives"]
        }
        lines += ["", "Historical alternatives remain unqualified:", ""]
        if not alternatives:
            lines += [
                "- No candidate numerical source record exists; original missing source rows are required."
            ]
        for key, candidate in sorted(alternatives.items()):
            citation, loc = candidate["citation"], candidate["location"]
            lines += [
                f"- `{key}`: {citation['title']} ({citation.get('authors', '')}); {citation.get('doi')}; PDF page {loc.get('pdf_page')}, table {loc.get('table')}, equation {loc.get('equation')}, row {loc.get('row_label')}. Source family/resolution: {candidate['source_force_field_family']}/{candidate['source_model_resolution']}. These locations do not resolve the exact PCFF request."
            ]
        lines += [""]
    lines += [
        "## Unchanged decision and scope",
        "",
        "`source_conflict_unresolved`; accepted_coefficients=0; automatic_assignment=False; production_validated=False; simulation_readiness=not_established. No profile, parameter bundle, converter result or MD workflow is created. Chemical transfer and scientific validation remain unestablished.",
        "",
        "J73 actual scientific executions: zero QM, MD, OpenMM, LAMMPS, AmberTools, converters, fitting, parameter generation and reserved validation. J72's four unintended SQM attempts and 46–48 MD steps remain preserved as a historical violation, excluded from scientific evidence by the unchanged C1 quarantine audit.",
        "",
    ]
    return "\n".join(lines)


def replay(repo, root):
    documents = reconstruct(repo, root)
    validate_requests(
        read(root / "source-intake.json"), documents["source-intake.json"]
    )
    for name, expected in documents.items():
        exact(read(root / name), expected, "offline_reconstruction:" + name)
    require(
        (root / "evidence-request-report.md").read_text() == render(documents),
        "request_report_mutation",
    )
    return {
        "offline": True,
        "instances": 378,
        "unresolved_counts": COUNTS,
        "scope": SCOPE,
    }


def build(repo, root):
    require(
        not (root / "source-intake.json").exists(), "do_not_overwrite_an_intake_attempt"
    )
    documents = reconstruct(repo, root)
    for name, value in documents.items():
        write(root / name, value)
    (root / "evidence-request-report.md").write_text(render(documents))
    return replay(repo, root)
