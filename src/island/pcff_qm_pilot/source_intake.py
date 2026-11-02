"""Immutable, source-separated Class-II intake; never grants runtime assignment.

An equation is not a coefficient. Missing records contain nulls, not zeros.
The retained research cosine AAT convention is not native FRC authority.
"""

import json
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256

import numpy as np

from island.exceptions import PCFFError
from island.forcefields.pcff.class2 import FAMILIES, FIELDS, resolve
from island.pcff_completion import canonical, digest, require, signed

from .hmds_hamiltonian import FORMS, validate

SCHEMA = "island.pcff-source-row.v1"
CATEGORIES = (
    "native_pcff",
    "author_primary",
    "literature_external",
    "research_derived",
)
COUNTS = {
    "angle-angle": 78,
    "bond-bond_1_3": 60,
    "end_bond-torsion_3": 60,
    "middle_bond-torsion_3": 60,
    "angle-torsion_3": 60,
    "angle-angle-torsion_1": 60,
}
POLICY = "exact_then_ordinary_equivalence_no_wildcard_no_lower_family_no_zero_fill"
ROW_KEYS = frozenset(
    (
        "coefficients",
        "contract",
        "dependencies",
        "eligible_static",
        "family",
        "graph_identity",
        "id",
        "identity",
        "normalized",
        "ordered_types",
        "required_evidence",
        "schema",
        "selection_policy",
        "source",
        "source_location",
        "source_lookup",
        "source_namespace",
        "stable_sites",
        "status",
        "strict_msi2lmp_usable",
        "symmetry_equivalent_types",
        "zero_filled",
    )
)


def orientations(family, types):
    other = (
        [types[3], types[1], types[2], types[0]]
        if family == "angle-angle"
        else types[::-1]
    )
    return [types, other]


def contract(family):
    units = list(FAMILIES[family][2])
    aat = family == "angle-angle-torsion_1"
    unresolved = "torsion" in family
    return {
        "fields": FIELDS[family],
        "source_units": units,
        "unit_symbols": {"E": "kcal/mol", "L": "angstrom", "A": "radian"},
        "scalar_energy": "K dt_IJK dt_JKL Phi (literal FRC); cosine research form is separate"
        if aat
        else FORMS[family],
        "analytic_force": "F_a = -sum_q (dE/dq)(dq/dx_a); all coefficients share the same scalar E",
        "prefactor": 1.0,
        "orientation": "J center K shared, exchange I/L only"
        if family == "angle-angle"
        else "reverse consecutive I-J-K-L; exchange left/right 3-coefficient blocks for EBT/AT",
        "multiplicity": "one per enumerated site; AA three angle pairs per neighbor triple; torsion n=1,2,3",
        "phase": "native Phi/AAT authority unresolved"
        if unresolved
        else "cos(n Phi), zero cross-term phase; no main-torsion phase shift applied to cross terms",
        "equation_status": "convention_unresolved"
        if unresolved
        else "declared_FRC_equation",
        "units_status": "convention_unresolved" if aat else "complete",
        "equilibrium_policy": "source bond/angle references, no independent fitting or global equilibrium revision",
    }


@dataclass(frozen=True)
class SourceRow:
    """Serialized immutable row; payload access returns a fresh mutable copy."""

    json_text: str

    def __post_init__(self):
        p = self.payload
        require(set(p) == ROW_KEYS, "source_row_schema_fields")
        require(
            type(p["eligible_static"]) is bool
            and type(p["strict_msi2lmp_usable"]) is bool,
            "source_row_boolean_schema",
        )
        require(p.get("schema") == SCHEMA, "source_row_schema")
        require(
            p["id"] == p["family"] + ":" + ",".join(map(str, p["stable_sites"])),
            "stable_source_row_id",
        )
        require(
            p["identity"] == digest({k: v for k, v in p.items() if k != "identity"}),
            "source_row_identity",
        )
        require(
            p["contract"] == contract(p["family"]),
            "source_row_units_equation_phase_multiplicity",
        )
        require(
            p["symmetry_equivalent_types"]
            == orientations(p["family"], p["ordered_types"]),
            "source_row_orientation",
        )
        require(not any("*" in t for t in p["ordered_types"]), "wildcard_forbidden")
        require(
            p["zero_filled"] is False and p["selection_policy"] == POLICY,
            "zero_fill_or_fallback_forbidden",
        )
        if p["source"] is None:
            require(
                p["coefficients"] is None
                and p["normalized"] is None
                and not p["eligible_static"],
                "missing_row_not_zero",
            )
            require(
                p["status"]
                == ("qm_fit_required" if p["family"] == "angle-angle" else "missing")
                and p["strict_msi2lmp_usable"] is False
                and p["source_location"] is None,
                "missing_source_declaration",
            )
        else:
            s = p["source"]
            require(
                s["category"] in CATEGORIES and s["source_family"] == "pcff",
                "cross_family_source_forbidden",
            )
            require(
                sha256(s["raw_file"].encode()).hexdigest() == s["sha256"],
                "source_file_hash_mutation",
            )
            # Re-parse the original file; supplied normalized cells are never trusted.
            from island.forcefields.pcff.class2 import inspect_pcff_class2
            from island.forcefields.pcff.source import PCFFSource, records

            src = PCFFSource(s["raw_file"].encode(), s["sha256"])
            cat = inspect_pcff_class2(src)["records"]
            decision = resolve(
                p["family"],
                p["ordered_types"],
                cat,
                records(src.inventory, "equivalence"),
            )
            require(decision["status"] == "assigned", "exact_source_row_required")
            ids = sorted({r["record_id"] for r in decision["selected"]})
            require(len(ids) == 1, "duplicate_source_rows")
            require(ids == s["row_ids"], "source_row_mapping_mutation")
            raw = [cat[r] for r in ids]
            require(
                p["coefficients"] == [r["original_values"] for r in raw],
                "raw_coefficients_mutation",
            )
            require(
                p["normalized"]
                == {
                    "values": decision["normalized_values"],
                    "units": raw[0]["normalized_units"],
                    "conversion_factors": raw[0]["conversion_factors"],
                },
                "source_normalization_mutation",
            )
            require(
                p["source_location"]
                == [
                    {"line": r["line"], "raw": r["raw"], "namespace": r["namespace"]}
                    for r in raw
                ],
                "source_location_mutation",
            )
            require(
                p["eligible_static"]
                == (
                    p["contract"]["equation_status"] != "convention_unresolved"
                    and s["category"] != "research_derived"
                ),
                "row_convention_gate",
            )
            require(
                p["status"]
                == (
                    "convention_unresolved"
                    if not p["eligible_static"]
                    else {
                        "native_pcff": "native_source",
                        "author_primary": "author_primary",
                        "literature_external": "literature_external",
                        "research_derived": "convention_unresolved",
                    }[s["category"]]
                ),
                "source_category_mutation",
            )
            require(
                p["strict_msi2lmp_usable"] == p["eligible_static"],
                "converter_row_compatibility",
            )

    @property
    def payload(self):
        return json.loads(self.json_text)

    def validate_against(self, expected):
        require(self.payload == expected, "source_row_graph_dependency_mutation")
        return self.payload


def intake_row(template, raw, *, expected_sha256, category, source_family="pcff"):
    """Parse one exact FRC row, preserving originals apart from normalization.

    A parsed row does not grant converter or preparation authority. Separate
    complete-file conversion and native-convention gates remain mandatory.
    """
    from island.forcefields.pcff.class2 import inspect_pcff_class2
    from island.forcefields.pcff.source import PCFFSource, records

    require(sha256(raw).hexdigest() == expected_sha256, "source_file_hash_mutation")
    require(
        category in CATEGORIES and source_family == "pcff",
        "cross_family_source_forbidden",
    )
    src = PCFFSource(raw, expected_sha256)
    try:
        cat = inspect_pcff_class2(src)["records"]
    except PCFFError as error:
        require(False, "source_format_invalid", detail=str(error))
    # Even identical duplicate highest-version source rows are ambiguous intake.
    decision = resolve(
        template["family"],
        template["ordered_types"],
        cat,
        records(src.inventory, "equivalence"),
    )
    require(decision["status"] == "assigned", "exact_source_row_required")
    ids = sorted({r["record_id"] for r in decision["selected"]})
    require(len(ids) == 1, "duplicate_source_rows")
    rows = [cat[r] for r in ids]
    p = {k: v for k, v in template.items() if k != "identity"}
    p.update(
        source={
            "category": category,
            "source_family": source_family,
            "sha256": expected_sha256,
            "raw_file": raw.decode(),
            "row_ids": ids,
        },
        source_location=[
            {"line": r["line"], "raw": r["raw"], "namespace": r["namespace"]}
            for r in rows
        ],
        coefficients=[r["original_values"] for r in rows],
        normalized={
            "values": decision["normalized_values"],
            "units": rows[0]["normalized_units"],
            "conversion_factors": rows[0]["conversion_factors"],
        },
        eligible_static=template["contract"]["equation_status"]
        != "convention_unresolved"
        and category != "research_derived",
        status={
            "native_pcff": "native_source",
            "author_primary": "author_primary",
            "literature_external": "literature_external",
            "research_derived": "convention_unresolved",
        }[category],
    )
    p["strict_msi2lmp_usable"] = p["eligible_static"]
    if not p["eligible_static"]:
        p["status"] = "convention_unresolved"
    return SourceRow(canonical(signed(p)))


def missing_rows(model):
    validate(model)
    inventory = {r["id"]: r for r in model["interaction_inventory"]}
    out = []
    for r in inventory.values():
        if r["assignment"] is not None:
            continue
        f = r["family"]
        require(f in COUNTS, "unexpected_missing_family")
        types = [model["typing"][str(i)] for i in r["sites"]]
        deps = [
            {
                "interaction": d,
                "source_rows": inventory[d]["assignment"]["row_ids"],
                "equilibrium": inventory[d]["assignment"]["values"][0],
                "units": "angstrom" if d.startswith("quartic_bond:") else "radian",
            }
            for d in r["dependencies"]
        ]
        p = signed(
            {
                "schema": SCHEMA,
                "id": r["id"],
                "family": f,
                "graph_identity": model["graph_identity"],
                "stable_sites": r["sites"],
                "ordered_types": types,
                "symmetry_equivalent_types": orientations(f, types),
                "source_namespace": "cff91",
                "source_lookup": {
                    "family": f,
                    "namespace": "cff91",
                    "source_sha256": model["source_hashes"]["pcff.frc"],
                    "row_id": None,
                    "status": "missing",
                    "per_instance_search_evidence": "source-review.json: parsed_PCFF_files/result/queries/"
                    + r["id"],
                },
                "source": None,
                "source_location": None,
                "coefficients": None,
                "normalized": None,
                "dependencies": deps,
                "contract": contract(f),
                "selection_policy": POLICY,
                "zero_filled": False,
                "status": "qm_fit_required" if f == "angle-angle" else "missing",
                "required_evidence": [
                    "exact authoritative numerical row",
                    "raw source hash and location",
                    "complete equation, units, phase, orientation and equilibrium interpretation",
                ],
                "strict_msi2lmp_usable": False,
                "eligible_static": False,
            }
        )
        out.append(SourceRow(canonical(p)).payload)
    require(Counter(r["family"] for r in out) == COUNTS, "missing_instance_enumeration")
    require(len({r["id"] for r in out}) == 378, "duplicate_interaction")
    return out


def completeness(model, rows):
    expected = missing_rows(model)
    require(
        len(rows) == len(expected) and len({r["id"] for r in rows}) == len(rows),
        "partial_or_duplicate_bundle",
    )
    require(
        {r["id"] for r in rows} == {r["id"] for r in expected},
        "missing_instance_mapping",
    )
    byid = {r["id"]: r for r in expected}
    for row in rows:
        parsed = SourceRow(canonical(row))
        # Graph/reference metadata is invariant even for a future supplied source.
        p = parsed.payload
        for key in (
            "family",
            "graph_identity",
            "stable_sites",
            "ordered_types",
            "symmetry_equivalent_types",
            "dependencies",
            "contract",
            "source_lookup",
            "source_namespace",
        ):
            require(
                p[key] == byid[p["id"]][key],
                "source_row_graph_dependency_mutation",
                field=key,
            )
    sources = {
        (r["source"]["category"], r["source"]["sha256"]) for r in rows if r["source"]
    }
    require(len(sources) <= 1, "native_external_source_mixing_forbidden")
    return all(r["eligible_static"] and r["strict_msi2lmp_usable"] for r in rows)


def require_complete(model, rows, converter, *, full_inventory_closed=False):
    require(
        completeness(model, rows),
        "parameter_incomplete",
        missing=[r["id"] for r in rows if not r["eligible_static"]],
    )
    require(full_inventory_closed is True, "full_ClassII_inventory_not_closed")
    require(
        converter["complete"] is True and converter["deterministic"] is True,
        "strict_converter_incomplete",
    )
    return True


def converter_gate(returncode, stdout, emitted_sections, expected_sections):
    """A success string cannot compensate for missing tables or missing rows."""
    missing = any(
        t in stdout.lower() for t in ("unable to find", "missing", "cannot find")
    )
    return (
        returncode == 0
        and not missing
        and set(emitted_sections) == set(expected_sections)
    )


def finite_static(energy, forces):
    require(
        np.isfinite(energy) and np.isfinite(forces).all(), "nonfinite_static_evaluator"
    )
    return True
