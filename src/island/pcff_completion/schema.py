"""Portable request schema. Semantic source/graph/unit checks run at intake."""

from . import VERSION

HASH = {"type": "string", "pattern": "^[0-9a-f]{64}$"}
TEXT = {"type": "string", "minLength": 1}
VECTOR = {"type": "array", "minItems": 1, "items": {"type": "number"}}
SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://island.invalid/schema/pcff-class2-completion-v1",
    "title": "External Class-II data-bound fitting request",
    "type": "object",
    "required": [
        "schema",
        "family",
        "candidate_profile",
        "source_hashes",
        "graph_identity",
        "typing_identity",
        "charge_identity",
        "atom_order",
        "charge_vector",
        "domain",
        "target_id",
        "query",
        "equation_identity",
        "files",
        "samples_file",
        "provenance",
        "units",
        "fit",
        "energy_scope",
        "isolation_evidence_file",
    ],
    "properties": {
        "schema": {"const": VERSION},
        "family": {"const": "pcff"},
        "candidate_profile": {"const": "island_external_pcff_pdms_class2_v1"},
        "source_hashes": {
            "type": "object",
            "minProperties": 3,
            "additionalProperties": HASH,
        },
        "graph_identity": HASH,
        "typing_identity": HASH,
        "charge_identity": HASH,
        "equation_identity": HASH,
        "atom_order": {
            "type": "array",
            "minItems": 1,
            "uniqueItems": True,
            "items": {"type": "integer"},
        },
        "charge_vector": {
            "type": "object",
            "minProperties": 1,
            "additionalProperties": {"type": "number"},
        },
        "domain": {"type": "object"},
        "target_id": TEXT,
        "query": TEXT,
        "files": {"type": "object", "minProperties": 1, "additionalProperties": HASH},
        "samples_file": TEXT,
        "energy_scope": {
            "enum": ["total_reference", "isolated_target_component", None]
        },
        "isolation_evidence_file": {"type": ["string", "null"]},
        "units": {
            "const": {
                "coordinates": "angstrom",
                "energy": "kcal/mol",
                "forces": "kcal/mol/angstrom",
                "hessian": "kcal/mol/angstrom^2",
                "angles": "radian",
                "charge": "e",
            }
        },
        "provenance": {"$ref": "#/$defs/provenance"},
        "fit": {"$ref": "#/$defs/fit"},
    },
    "additionalProperties": False,
    "$defs": {
        "provenance": {
            "type": "object",
            "required": [
                "kind",
                "citation",
                "license",
                "license_file",
                "permission_status",
                "data_origin",
                "reference_version",
            ],
            "properties": {
                "kind": {"enum": ["qm", "experimental", "reference_engine"]},
                "citation": TEXT,
                "license": TEXT,
                "license_file": TEXT,
                "permission_status": {"enum": ["licensed", "permission_granted"]},
                "data_origin": TEXT,
                "reference_version": TEXT,
                "method": {"type": ["string", "null"]},
                "basis": {"type": ["string", "null"]},
                "software": {"type": ["string", "null"]},
            },
            "allOf": [
                {
                    "if": {"properties": {"kind": {"const": "qm"}}},
                    "then": {
                        "required": ["method", "basis", "software"],
                        "properties": {"method": TEXT, "basis": TEXT, "software": TEXT},
                    },
                }
            ],
        },
        "fit": {
            "type": "object",
            "required": [
                "objective",
                "fields",
                "initial_values",
                "lower_bounds",
                "upper_bounds",
                "constraints",
                "optimizer",
                "weights",
            ],
            "properties": {
                "objective": {
                    "const": "weighted_isolated_component_energy_least_squares"
                },
                "fields": {"type": "array", "minItems": 1, "items": TEXT},
                "initial_values": VECTOR,
                "lower_bounds": VECTOR,
                "upper_bounds": VECTOR,
                "constraints": {
                    "const": [
                        "source_equation",
                        "fixed_verified_equilibria",
                        "exact_graph_domain",
                    ]
                },
                "optimizer": {"const": "numpy.linalg.lstsq_SVD_v1"},
                "weights": {
                    "type": "object",
                    "minProperties": 1,
                    "additionalProperties": {"type": "number", "exclusiveMinimum": 0},
                },
            },
        },
        "sample": {
            "type": "object",
            "required": ["id", "coordinates", "atom_order", "split"],
            "properties": {
                "id": TEXT,
                "coordinates": {
                    "type": "array",
                    "items": {
                        "type": "array",
                        "minItems": 3,
                        "maxItems": 3,
                        "items": {"type": "number"},
                    },
                },
                "atom_order": {"type": "array", "items": {"type": "integer"}},
                "split": {"enum": ["train", "held_out"]},
                "energy": {"type": ["number", "null"]},
                "forces": {"type": ["array", "null"]},
                "hessian": {"type": ["array", "null"]},
                "vibrations": {"type": ["object", "null"]},
                "observables": {"type": ["array", "null"]},
            },
            "description": "At least one reference channel; exact graph dimensions and units validated by intake.",
        },
        "derived_record": {
            "type": "object",
            "required": [
                "origin",
                "formula",
                "source_equation",
                "input_rows",
                "conversion",
                "implementation",
                "uncertainty",
                "independent_reconstruction",
                "identity",
            ],
            "properties": {"origin": {"const": "derived_external"}, "identity": HASH},
            "description": "No coefficient derivation adapter is currently authorized. This schema is not authority to create a term.",
        },
        "fitted_record": {
            "type": "object",
            "required": [
                "origin",
                "contract",
                "dataset_identity",
                "parameters",
                "objective",
                "weights",
                "constraints",
                "bounds",
                "initial_values",
                "optimizer",
                "convergence",
                "uncertainty",
                "held_out_validation",
                "implementation",
                "final_parameter_checksum",
                "identity",
            ],
            "properties": {
                "origin": {"const": "fitted_external"},
                "identity": HASH,
                "final_parameter_checksum": HASH,
            },
        },
    },
    "description": "An unfilled template is deliberately invalid. Hash identity is canonical JSON SHA-256 without trailing newline. Intake alone cannot grant fit/profile/MD authority.",
}
