"""Offline J79 literature claims, topology crosswalk and convention audit.

No physical evaluator, fitter, converter or runtime registry is called. The
primary paper is unavailable: reported numbers are claims, not extracted rows.
"""

import argparse
import json
import subprocess
from collections import Counter
from decimal import Decimal, localcontext
from pathlib import Path

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

BASE = "3c9f9fb5066ac2132d8a9aa46f228a3e53a77ddd"
BRANCH = "codex/phase-4j79-pdms-literature-reconciliation"
ROOT = Path("docs/evidence/phase_4j79")
FREEZE_SHA = "ba846da89eacb520c97cb99399a1d89b035afa05294327a32eca3ac98c913207"
INVENTORY = "docs/evidence/phase_4j78/unique-term-inventory.json"
SOURCE_INVENTORY = "docs/evidence/phase_4j75/source-file-inventory.json"
IMPLEMENTATION = (
    "scripts/phase_4j79_literature.py",
    "scripts/validate_phase_4j79_receipt.py",
    "tests/test_phase_4j79.py",
    "docs/phase_4j79.md",
)
MODEL_FILES = (
    "src/island/forcefields/parameters/models.py",
    "src/island/forcefields/pcff/model.py",
    "src/island/forcefields/pcff/terms.py",
    "src/island/forcefields/pcff/registry.py",
    "src/island/pcff_pdms/extension/inventory.py",
    "src/island/pcff_pdms/extension/qualification.py",
)
STATUSES = (
    "literature_candidate_exact_form",
    "literature_candidate_convention_blocked",
    "end_group_unresolved",
    "cross_term_unresolved",
    "not_covered_by_source",
)
PRIMARY = {
    "authors": ["Lydia Fritz", "Dieter Hofmann"],
    "title": "Molecular dynamics simulations of the transport of water-ethanol mixtures through polydimethylsiloxane membranes",
    "journal": "Polymer",
    "volume": "38",
    "year": 1997,
    "pages": "1035–1045",
    "PII": "S0032-3861(96)00600-3",
    "DOI": "10.1016/S0032-3861(96)00600-3",
}
VDW = "epsilon_ij*[2*(sigma_ij/r_ij)^9-3*(sigma_ij/r_ij)^6]"
SCOPE = {
    "J79_new_QM_attempts": 0,
    "J79_MD_steps": 0,
    "converter_executions": 0,
    "fitting_attempted": False,
    "fitting_executions": 0,
    "parameter_generation_executions": 0,
    "scientific_evaluator_executions": 0,
    "reserved_validation_executions": 0,
    "accepted_coefficients": 0,
    "native_exact_rows_resolved": 0,
    "native_unresolved_instances": 378,
    "runtime_registered": False,
    "automatic_assignment_enabled": False,
    "derived_profile_executable": False,
    "native_pcff_complete": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "chemical_transfer": "not_established",
    "scientific_validation": "not_established",
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
    "source_adjudication_excluded_quarantine": True,
}


def reported_rows():
    """Exact decimal strings from the request; no paper provenance is implied."""
    data = []

    def add(family, atoms, fields, values):
        data.append(
            {
                "id": family + ":" + atoms,
                "family": family,
                "ordered_elements": atoms.split("-"),
                "raw_reported_values": dict(zip(fields, values.split(), strict=True)),
            }
        )

    for atom, q in (
        ("Si", "+0.640"),
        ("O", "-0.440"),
        ("C", "-0.259"),
        ("H", "+0.053"),
    ):
        add("partial_charge", atom, ["q"], q)
    for atoms, values in (
        ("Si-O", "1.6562 306.1232 -517.3424 673.7067"),
        ("Si-C", "1.9073 157.0049 -237.7023 356.0328"),
        ("C-H", "1.1010 345.0000 -691.8900 844.6000"),
    ):
        add("quartic_bond", atoms, ["b0", "K2", "K3", "K4"], values)
    for atoms, values in (
        ("Si-O-Si", "157.0260 9.0740 -19.5576 8.5000"),
        ("O-Si-O", "110.6930 70.3069 -6.9375 0"),
        ("C-Si-O", "114.9060 23.0218 -31.3993 -924.9814"),
        ("H-C-Si", "111.5360 30.2481 -15.5255 0"),
        ("C-Si-C", "114.9060 23.0218 -31.3993 24.9814"),
    ):
        add("quartic_angle", atoms, ["theta0", "H2", "H3", "H4"], values)
    for atoms, values in (
        ("Si-O-Si-O", "0.3000 0.3658 0"),
        ("Si-O-Si-C", "0 0 -0.1300"),
        ("H-C-Si-O", "-1.3513 0 -0.0580"),
        ("H-C-Si-C", "-1.3513 0 -0.0580"),
    ):
        add("torsion_3", atoms, ["V1", "V2", "V3"], values)
    for atom, values in (
        ("Si", "4.2840 0.07000"),
        ("O", "3.3500 0.24000"),
        ("C", "4.0100 0.05400"),
        ("H", "2.9950 0.02000"),
    ):
        add("lj_9_6", atom, ["sigma", "epsilon"], values)
    return data


def convert(value, unit):
    """Unit arithmetic only; never changes the force function or infers a unit."""
    require(
        unit in ("angstrom_to_nm", "kJ_to_kcal"), "explicit_unit_conversion_required"
    )
    with localcontext() as context:
        context.prec = 28
        return str(
            Decimal(value) * Decimal("0.1")
            if unit == "angstrom_to_nm"
            else Decimal(value) / Decimal("4.184")
        )


def torsion_contract(row):
    require(
        row.get("source_verified") is True
        and row.get("source_equation")
        and row.get("prefactor") is not None
        and row.get("phase_convention")
        and row.get("angle_origin")
        and row.get("orientation_rule")
        and row.get("multiplicities") == [1, 2, 3],
        "torsion_sign_phase_prefactor_unresolved",
    )
    return row


def literature():
    rows = []
    for raw in reported_rows():
        values = raw["raw_reported_values"]
        derived = {}
        for k, v in values.items():
            if k in ("b0", "sigma"):
                derived[k + "_nm"] = {
                    "value": convert(v, "angstrom_to_nm"),
                    "factor": "0.1",
                    "condition": "reported angstrom unit; primary unit unverified",
                }
            elif k not in ("q", "theta0"):
                derived[k + "_energy_kcal"] = {
                    "value": convert(v, "kJ_to_kcal"),
                    "factor": "1/4.184",
                    "condition": "IF source energy is kJ/mol; coordinate powers unchanged; not assignable",
                }
        rows.append(
            signed(
                {
                    **raw,
                    "citation": PRIMARY,
                    "evidence_origin": "user_reported_transcription_not_primary_extraction",
                    "source_family_reported": "PDMS PCFF as described by Fritz/Hofmann; not native row authority",
                    "source_hash": None,
                    "source_page": None,
                    "source_table": None,
                    "source_equation_location": None,
                    "source_text": None,
                    "source_units": None,
                    "source_verified": False,
                    "reported_units": {
                        "energy": "kJ/mol reported; unverified for each coefficient",
                        "length": "angstrom",
                        "theta0": "degree"
                        if raw["family"] == "quartic_angle"
                        else None,
                        "charge": "elementary_charge" if "q" in values else None,
                    },
                    "angular_force_constant_coordinate": "unresolved_degree_vs_radian",
                    "source_equation": None,
                    "reported_equation": VDW if raw["family"] == "lj_9_6" else None,
                    "prefactor": None,
                    "phase_convention": None,
                    "multiplicities": None,
                    "angle_origin": None,
                    "orientation_rule": None,
                    "atom_order": raw["ordered_elements"],
                    "atom_typing": "element_labels_only; source row/type IDs absent",
                    "derived_conditional_conversions": derived,
                    "warning": (
                        "C-Si-O H4=-924.9814 retained exactly; verify original table, do not repair by analogy"
                        if raw["id"] == "quartic_angle:C-Si-O"
                        else "original table, units and conventions required"
                    ),
                    "usable_for_runtime_assignment": False,
                }
            )
        )
    return signed(
        {
            "schema": "island.j79-literature-evidence.v1",
            "bibliography": PRIMARY,
            "reading_status": "metadata_only; primary_full_text_not_available_not_read",
            "primary_artifact_sha256": None,
            "parameter_rows_source_extracted": 0,
            "reported_rows": rows,
            "paper_mixing_rule": None,
            "paper_dielectric": None,
            "paper_cutoff": None,
            "paper_exclusions": None,
            "paper_angle_coefficient_radians_verified": False,
            "paper_torsion_sign_phase_verified": False,
            "license": "publisher access restricted; redistribution permission not established",
        }
    )


def elements(types):
    return [t.split("_")[1] for t in types]


def candidates(term, rows):
    seq = elements(term["ordered_types"])
    return [
        r["id"]
        for r in rows
        if r["family"] == term["family"]
        and seq in (r["ordered_elements"], r["ordered_elements"][::-1])
    ]


def crosswalk(j78, evidence):
    require(len(j78["terms"]) == 129, "incomplete_J78_terms")
    ids = [t["id"] for t in j78["terms"]]
    require(len(set(ids)) == len(ids), "duplicate_term_identity")
    rows = []
    for t in j78["terms"]:
        family = t["family"]
        cand = candidates(t, evidence["reported_rows"])
        cap = any(x.endswith("_T") for x in t["ordered_types"])
        cross = family not in (
            "partial_charge",
            "lj_9_6",
            "quartic_bond",
            "quartic_angle",
            "torsion_3",
        )
        status = (
            "cross_term_unresolved"
            if cross
            else "end_group_unresolved"
            if cap
            else "literature_candidate_convention_blocked"
            if cand
            else "not_covered_by_source"
        )
        instances = [i for i in j78["instances"] if i["parameter_id"] == t["id"]]
        rows.append(
            signed(
                {
                    "J78_parameter_id": t["id"],
                    "J78_parameter_identity": t["identity"],
                    "ordered_types": t["ordered_types"],
                    "ordered_elements": elements(t["ordered_types"]),
                    "family": family,
                    "J78_equation": t["equation"],
                    "J78_units": t["units"],
                    "reported_candidate_ids": cand,
                    "element_sequence_match_only": bool(cand),
                    "exact_chemical_environment_verified": False,
                    "end_group_involved": cap,
                    "cross_term": cross,
                    "status": status,
                    "blocking_reasons": [
                        "primary tables/equations and source atom types unavailable",
                        *(["trimethylsilyl cap transfer unproven"] if cap else []),
                        *(
                            [
                                "no cross-term values supplied; absence in paper not yet verified"
                            ]
                            if cross
                            else []
                        ),
                        *(
                            ["no reported matching element tuple"]
                            if not cand and not cross
                            else []
                        ),
                    ],
                    "instance_count": len(instances),
                    "J78_instance_ids": [i["id"] for i in instances],
                    "J78_instance_identities_sha256": digest(
                        [i["identity"] for i in instances]
                    ),
                    "graph_identities": sorted(
                        {i["graph_identity"] for i in instances}
                    ),
                    "native_row_resolved": False,
                    "extension_term_resolved": False,
                }
            )
        )
    return signed(
        {
            "schema": "island.j79-crosswalk.v1",
            "source_J78_identity": j78["identity"],
            "unique_terms": len(rows),
            "site_instances": sum(r["instance_count"] for r in rows),
            "status_counts": {s: sum(r["status"] == s for r in rows) for s in STATUSES},
            "status_precedence": "cross_term, end_group, convention_blocked, not_covered; orthogonal flags retain overlap",
            "rows": rows,
        }
    )


def neutrality(j78):
    q = {
        r["ordered_elements"][0]: Decimal(r["raw_reported_values"]["q"])
        for r in reported_rows()
        if r["family"] == "partial_charge"
    }
    repeat = {"Si": 1, "O": 1, "C": 2, "H": 6}
    result = {"repeat_SiOC2H6_charge_e": str(sum(q[k] * n for k, n in repeat.items()))}
    for graph in j78["graphs"]:
        counts = Counter(elements([t for _, t in graph["types"]]))
        result[graph["fixture"]] = {
            "element_counts": dict(counts),
            "charge_if_repeat_values_transferred_e": str(
                sum(q[k] * n for k, n in counts.items())
            ),
            "transfer_allowed": False,
            "reason": "neutral repeat does not establish capped oligomer charge policy; no renormalization",
        }
    return result


def conventions(j78, models):
    require(
        'functional_form: str = field(default="lj_12_6"' in models[MODEL_FILES[0]],
        "generic_LJ_model_changed",
    )
    require(
        "(rmin / distance) ** 9" in models[MODEL_FILES[2]],
        "isolated_9_6_contract_changed",
    )
    require('"lj_9_6"' in models[MODEL_FILES[4]], "J78_9_6_contract_changed")
    return signed(
        {
            "generic_parameter_model": {
                "file": MODEL_FILES[0],
                "class": "LennardJonesParameter",
                "form": "4*epsilon*((sigma/r)^12-(sigma/r)^6)",
                "compatible_with_reported_paper_9_6": False,
                "reason": "different powers and distance meaning; unit conversion cannot repair this",
            },
            "isolated_PCFF_model": {
                "files": list(MODEL_FILES[1:3]),
                "form": "epsilon*(2*(rmin/r)^9-3*(rmin/r)^6)",
                "mixing": {
                    "rmin_ij": "((ri^6+rj^6)/2)^(1/6)",
                    "epsilon_ij": "2*sqrt(eps_i*eps_j)*ri^3*rj^3/(ri^6+rj^6)",
                },
                "called": False,
                "paper_mixing_compatibility": "unresolved_primary_equation_required",
                "execution_or_assignment_qualification": False,
            },
            "J78_proposal": {
                "file": MODEL_FILES[4],
                "form": "lj_9_6 minimum-distance",
                "symbolic_form_match_to_reported_paper": True,
                "paper_sigma": "minimum-distance IF reported 9-6 equation is verified; not 12-6 sigma",
                "convention_compatibility": "blocked; same powers alone insufficient",
            },
            "paper_conventions": {
                k: {"value": None, "status": "primary_source_required"}
                for k in (
                    "mixing_rule",
                    "dielectric",
                    "cutoff",
                    "tail",
                    "exclusions_1_2",
                    "exclusions_1_3",
                    "scaling_1_4",
                    "charges_origin",
                    "long_range_electrostatics",
                    "angle_coefficient_degree_or_radian",
                    "torsion_phase_sign_prefactor",
                    "orientation",
                )
            },
            "J78_policy_audit": [
                {
                    "J78_policy_id": p["id"],
                    "J78_status": p["status"],
                    "J79_status": "unresolved",
                    "source_value": None,
                }
                for p in j78["required_policies"]
            ],
            "torsion": {
                "reported_coefficients_preserve_sign": True,
                "paper_three_term_equation": None,
                "J78_expression": "sum(n=1..3,Vn*(1-cos(n*phi-phase_n)))",
                "translation": None,
                "reason": "prefactor, dihedral origin, sign, phases and reversal not verified; do not infer phase from signed Vn",
            },
            "angle": {
                "reported_equilibrium_unit": "degree",
                "coefficient_angular_coordinate": None,
                "conversion_rule_if_verified": "theta_rad=theta_deg*pi/180; K_rad,n=K_deg,n*(180/pi)^n",
                "applied_to_coefficients": False,
            },
            "charge_arithmetic_only": neutrality(j78),
            "runtime_registration": False,
            "automatic_assignment_enabled": False,
        }
    )


def reconstruct(root):
    scope()
    require(
        pin(root / "frozen-inputs.json") == FREEZE_SHA,
        "fabricated_or_resigned_input_freeze",
    )
    frozen = read(root / "frozen-inputs.json")
    verify_signed(frozen)
    exact(frozen["scope"], SCOPE, "execution_scope_mutation")
    for n, sha in frozen["capture_files"].items():
        require(pin(root / n) == sha, "captured_source_hash_mutation:" + n)
    j78 = read(root / "j78-input.json")
    evidence = literature()
    walk = crosswalk(j78, evidence)
    used = {c for r in walk["rows"] for c in r["reported_candidate_ids"]}
    cross = [r for r in walk["rows"] if r["cross_term"]]
    end = [r for r in walk["rows"] if r["end_group_involved"]]
    counts = {
        "literature_rows_inspected": 20,
        "primary_source_rows_inspected": 0,
        "reported_rows_element_sequence_matching": len(used),
        "literature_rows_chemically_matching": 0,
        "literature_rows_convention_compatible": 0,
        "literature_rows_blocked_by_convention": 20,
        "native_exact_rows_resolved": 0,
        "unresolved_cross_terms": len(cross),
        "unresolved_end_group_terms": len(end),
        "unresolved_end_group_non_cross_terms": sum(not r["cross_term"] for r in end),
        "unresolved_extension_terms": 129,
        "unresolved_scalar_parameters": 365,
        "unresolved_policies": 11,
        "source_verified_numerical_cells": 0,
        "reported_numerical_cells": sum(
            len(r["raw_reported_values"]) for r in evidence["reported_rows"]
        ),
    }
    missing_classes = [
        "bond-bond",
        "bond-angle",
        "angle-angle",
        "torsion-bond",
        "torsion-angle",
        "improper/out-of-plane",
        "torsion-torsion",
    ]
    unresolved = signed(
        {
            "required_classes": [
                {
                    "family": f,
                    "status": "unresolved",
                    "reason": "not supplied in reported claims; paper coverage not verifiable without full text",
                }
                for f in missing_classes
            ],
            "unresolved_term_ids": [r["J78_parameter_id"] for r in walk["rows"]],
            "end_group_term_ids": [r["J78_parameter_id"] for r in end],
            "cross_term_ids": [r["J78_parameter_id"] for r in cross],
            "end_groups": "HMDS/DP2/DP4/DP6 cap coverage not established; internal-repeat transfer forbidden",
            "required_evidence": [
                "hash-bound original Fritz/Hofmann paper, pages 1035–1045 and parameter tables",
                "exact quartic equations, coefficient units and degree/radian definition",
                "exact three-term Fourier equation, phases, prefactor and dihedral origin",
                "9-6 combining rule, exclusions/scaling, dielectric, electrostatic and cutoff/tail settings",
                "source atom-type IDs and internal versus trimethylsilyl end-group applicability",
                "all required Class-II cross terms and equilibrium dependencies; absent does not mean zero",
                "verify C-Si-O H4=-924.9814 against the original cell",
            ],
            "cross_term_count_includes_end_group_overlap": True,
            "counts": counts,
        }
    )
    return {
        "literature-evidence.json": evidence,
        "j78-literature-crosswalk.json": walk,
        "convention-compatibility.json": conventions(
            j78, read(root / "model-inputs.json")
        ),
        "unresolved-terms.json": unresolved,
        "reconciliation-report.json": signed(
            {
                "decision": "literature_reconciliation_partial",
                "reason": "reported claims reconciled offline; primary paper inaccessible, all source/convention gates blocked",
                "counts": counts,
                "crosswalk_status_counts": walk["status_counts"],
                "scope": SCOPE,
                "claims_are_not_extracted_source_values": True,
                "paper_read": False,
                "QM_targets_completed": 0,
                "QM_targets_planned_unchanged": 423,
            }
        ),
    }


def capture(repo, root):
    scope()
    require(not root.exists(), "do_not_overwrite_capture")
    inv = read(repo / INVENTORY)
    verify_signed(inv)
    # Retain exact original unique/instance identities, order, neighborhoods and
    # dependencies. Coordinates are neither copied nor used.
    compact = {k: inv[k] for k in ("identity", "terms", "graphs", "required_policies")}
    compact["instances"] = [
        {
            k: i[k]
            for k in (
                "id",
                "identity",
                "parameter_id",
                "fixture",
                "sites",
                "source_order_types",
                "graph_identity",
                "equilibrium_parameter_dependencies",
            )
        }
        for i in inv["instances"]
    ]
    write(root / "j78-input.json", compact)
    write(root / "model-inputs.json", {n: (repo / n).read_text() for n in MODEL_FILES})
    source_inv = read(repo / SOURCE_INVENTORY)
    write(
        root / "local-source-search.json",
        {
            "method": "inspect immutable J75 inventory identities; no unrelated papers parsed",
            "source_inventory_sha256": pin(repo / SOURCE_INVENTORY),
            "files": [
                {
                    "relative_path": r["relative_path"],
                    "sha256": r["sha256"],
                    "citation": r["citation"],
                }
                for r in source_inv["files"]
            ],
            "requested_primary_present_in_inventory": False,
            "primary_file_opened": False,
        },
    )
    # Retain only bibliographic facts/hashes, not publisher full-text assets.
    metadata = read(Path("/tmp/j79-source-access-0"))["message"]
    write(
        root / "source-access.json",
        {
            "crossref_response_sha256": pin(Path("/tmp/j79-source-access-0")),
            "crossref_identity": {
                k: metadata[k]
                for k in (
                    "DOI",
                    "title",
                    "author",
                    "publisher",
                    "page",
                    "volume",
                    "published",
                    "link",
                    "license",
                )
            },
            "access_attempts": read(Path("/tmp/j79-source-access.json"))
            + read(Path("/tmp/j79-source-elsevier.json")),
            "publisher_XML_response_sha256": pin(Path("/tmp/j79-source-elsevier-0")),
            "publisher_XML_contains_full_text": False,
            "primary_artifact_sha256": None,
            "full_text_read": False,
        },
    )
    history = {}
    tree = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    for line in tree.splitlines():
        meta, name = line.split("\t", 1)
        p = repo / name
        history[name] = (
            digest(str(p.readlink())) if meta.startswith("120000") else pin(p)
        )
    write(root / "historical-files.json", history)
    names = (
        "j78-input.json",
        "model-inputs.json",
        "local-source-search.json",
        "source-access.json",
        "historical-files.json",
    )
    write(
        root / "frozen-inputs.json",
        signed(
            {
                "schema": "island.j79-input-freeze.v1",
                "base": BASE,
                "fetched_origin_main": "594593ee656aace09133e0eba99648225c303ca3",
                "scope": SCOPE,
                "original_inputs": {
                    n: pin(repo / n)
                    for n in (INVENTORY, SOURCE_INVENTORY, *MODEL_FILES)
                },
                "capture_files": {n: pin(root / n) for n in names},
            }
        ),
    )
    return {
        "freeze_sha256": pin(root / "frozen-inputs.json"),
        "historical_files": len(history),
    }


def replay(root, build=False):
    docs = reconstruct(root)
    for name, value in docs.items():
        if build:
            require(not (root / name).exists(), "do_not_overwrite_attempt")
            write(root / name, value)
        exact(read(root / name), value, "resigned_literature_semantic_mutation:" + name)
    return docs["reconciliation-report.json"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("capture", "build", "replay"))
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    a = p.parse_args()
    try:
        scope()
        print(
            json.dumps(
                capture(a.repo, a.directory)
                if a.mode == "capture"
                else replay(a.directory, a.mode == "build"),
                sort_keys=True,
            )
        )
        return 0
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as e:
        print(json.dumps({"status": "blocked", "reason": str(e), "scope": SCOPE}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
