"""OPLS-only, source-bound exchange contracts. Numeric payloads are ephemeral."""

import json
import re
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from math import isfinite

from island.forcefields.assignment_core.schema import identity, require
from island.forcefields.bundle_adapters.schema import keys, loads, text

PROVIDERS = ("polypargen", "ligpargen", "oplsaa_db", "audited_bundle")
RESULTS = (
    "provenance_only",
    "provisional_parameterization",
    "audited_exact_graph",
    "parameter_bundle_incomplete",
    "source_bundle_required",
)
FLAGS = {
    "accepted_coefficients": 0,
    "runtime_registration": False,
    "automatic_assignment_enabled": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "chemical_family_transfer": "not_established",
}
FAMILIES = ("bond", "angle", "torsion", "improper", "charge", "nonbonded")
DIMENSIONS = ("atom_typing", *FAMILIES, "interaction_14", "topology")
ARITY = dict(zip(FAMILIES, (2, 3, 4, 4, 1, 1), strict=True))
# No numerical force-field parameters occur in these software form contracts.
FORMS = {
    "bond_harmonic": {
        "family": "bond",
        "fields": ["r0", "k"],
        "units": ["nm", "kJ/mol/nm^2"],
        "equation": "k*(r-r0)^2/2",
        "ordering": "reversal_equivalent",
    },
    "angle_harmonic": {
        "family": "angle",
        "fields": ["theta0", "k"],
        "units": ["degree", "kJ/mol/radian^2"],
        "equation": "k*(theta-theta0)^2/2; displacements in radians",
        "ordering": "reversal_equivalent",
    },
    "torsion_periodic": {
        "family": "torsion",
        "fields": ["phase", "k", "periodicity"],
        "units": ["degree", "kJ/mol", "integer"],
        "equation": "sum k*(1+cos(n*phi-phase))",
        "ordering": "ordered",
    },
    "torsion_rb": {
        "family": "torsion",
        "fields": ["c" + str(i) for i in range(6)],
        "units": ["kJ/mol"] * 6,
        "equation": "sum Cn*cos(phi-pi)^n; n=0..5",
        "ordering": "ordered",
    },
    "torsion_opls": {
        "family": "torsion",
        "fields": ["v" + str(i) for i in range(1, 5)],
        "units": ["kJ/mol"] * 4,
        "equation": "(V1*(1+cos(phi))+V2*(1-cos(2*phi))+V3*(1+cos(3*phi))+V4*(1-cos(4*phi)))/2",
        "ordering": "ordered",
    },
    "improper_harmonic": {
        "family": "improper",
        "fields": ["phase", "k"],
        "units": ["degree", "kJ/mol/radian^2"],
        "equation": "k*wrapped(phi-phase)^2/2",
        "ordering": "ordered_source",
    },
    "charge": {
        "family": "charge",
        "fields": ["q"],
        "units": ["e"],
        "equation": "source_atomic_charge",
        "ordering": "site",
    },
    "lj12_6": {
        "family": "nonbonded",
        "fields": ["sigma", "epsilon"],
        "units": ["nm", "kJ/mol"],
        "equation": "4*epsilon*((sigma/r)^12-(sigma/r)^6)",
        "ordering": "site",
    },
}
CONVENTIONS = (
    "energy_unit",
    "length_unit",
    "charge_unit",
    "angle_equilibrium_unit",
    "angle_force_unit",
    "nonbonded",
    "mixing_rule",
    "exclusions",
    "scale14_lj",
    "scale14_coulomb",
    "long_range",
    "cutoff",
    "tail",
    "torsion_phase",
    "torsion_normalization",
    "improper_ordering",
    "atom_typing",
    "charge_origin",
)


def number(value):
    require(isinstance(value, str), "numeric_literals_must_preserve_source_strings")
    require(
        re.fullmatch(
            r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", value
        )
        is not None,
        "unsupported_engine_numeric_literal",
    )
    try:
        result = Decimal(value)
    except InvalidOperation:
        raise ValueError("invalid_numeric_literal") from None
    require(result.is_finite(), "nonfinite_numeric_literal")
    native = float(result)
    require(
        isfinite(native) and (native != 0 or result == 0),
        "native_numeric_overflow_or_underflow",
    )
    return result


def convention_blockers(c):
    keys(c, CONVENTIONS, "complete_convention_contract_required")
    for value in c.values():
        text(value, "missing_convention")
    expected = {
        "energy_unit": "kJ/mol",
        "length_unit": "nm",
        "charge_unit": "e",
        "angle_equilibrium_unit": "degree",
        "angle_force_unit": "radian",
        "nonbonded": "lj12_6_sigma_zero",
        "mixing_rule": "geometric_sigma_epsilon",
        "exclusions": "graph_distance_1_2",
        "torsion_phase": "cos(n*phi-phase)",
        "torsion_normalization": "form_explicit",
        "improper_ordering": "source_ordered_center_second",
    }
    blockers = [
        "unsupported_convention:" + k
        for k, v in c.items()
        if v.strip().lower() in ("unknown", "unresolved", "default", "missing", "n/a")
        or (k in expected and expected[k] != v)
    ]
    for name in ("scale14_lj", "scale14_coulomb"):
        value = number(c[name])
        require(0 <= value <= 1, "invalid_14_scaling")
    return sorted(blockers)


def row_key(family, sites):
    sites = tuple(sites)
    if family in ("bond", "angle", "torsion"):
        sites = min(sites, sites[::-1])
    return family + ":" + ",".join(map(str, sites))


def validate_row(row, atom_order):
    keys(
        row,
        (
            "id",
            "family",
            "sites",
            "form",
            "values",
            "applicability",
            "reason",
            "source",
        ),
        "canonical_row_schema",
    )
    family = row["family"]
    require(family in FAMILIES, "OPLS_Class_I_only")
    require(
        isinstance(row["sites"], list)
        and len(row["sites"]) == ARITY[family]
        and all(type(i) is int and i in atom_order for i in row["sites"])
        and len(set(row["sites"])) == len(row["sites"]),
        "invalid_ordered_sites",
    )
    text(row["id"], "record_identity_required")
    keys(row["source"], ("path", "sha256", "pointer"), "source_row_pointer_required")
    text(row["source"]["pointer"], "source_row_pointer_required")
    require(
        row["applicability"] in ("required", "not_applicable"), "invalid_applicability"
    )
    require(isinstance(row["values"], list), "numeric_array_required")
    if row["applicability"] == "not_applicable":
        require(
            family == "improper" and not row["values"] and row["form"] is None,
            "only_explicit_improper_omission_supported",
        )
        text(row["reason"], "evidenced_omission_reason_required")
        return
    require(
        row["form"] in FORMS and FORMS[row["form"]]["family"] == family,
        "unsupported_functional_form",
    )
    require(row["values"], "missing_values_no_zero_substitution")
    require(
        len(row["values"]) == 1 or row["form"] == "torsion_periodic",
        "duplicate_parameter_alias",
    )
    for term in row["values"]:
        require(
            isinstance(term, list) and len(term) == len(FORMS[row["form"]]["fields"]),
            "incomplete_numeric_fields",
        )
        values = [number(v) for v in term]
        if family == "nonbonded":
            require(values[0] > 0 and values[1] >= 0, "invalid_nonbonded_bounds")
        if family in ("bond", "angle"):
            require(values[0] > 0 and values[1] >= 0, "invalid_primary_bounds")
        if row["form"] == "torsion_periodic":
            require(
                re.fullmatch(r"[+]?[0-9]+", term[2]) is not None and values[2] > 0,
                "invalid_periodicity",
            )
    require(
        len({tuple(t) for t in row["values"]}) == len(row["values"]),
        "duplicate_fourier_component",
    )


@dataclass(frozen=True)
class CanonicalOPLS:
    """Immutable JSON snapshot. Never persist this coefficient-bearing object in Git."""

    json_text: str = field(repr=False)

    @classmethod
    def create(cls, value):
        return cls(
            json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        )

    @property
    def payload(self):
        return loads(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


def sealed(value):
    return {**value, "identity": identity(value)}


def verify_seal(value):
    require(
        value.get("identity")
        == identity({k: v for k, v in value.items() if k != "identity"}),
        "receipt_identity_mutation",
    )
