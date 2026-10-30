"""JSON schemas for immutable, deliberately unfilled targets and raw intake."""

from .targets import SETTINGS

TARGET_SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "island.pcff-qm-target.v1",
    "type": "object",
    "required": [
        "schema",
        "identity",
        "query",
        "interaction_type",
        "atom_type_environment",
        "molecular_fixture",
        "graph_identity",
        "typing_identity",
        "charge_identity",
        "target_site",
        "functional_form",
        "charge_model",
        "units",
        "required_deformation",
        "identifiability_requirements",
        "held_out_validation_requirements",
        *SETTINGS,
    ],
    "properties": {
        "schema": {"const": "island.pcff-qm-target.v1"},
        "status": {"const": "qm_fit_required"},
        **{k: {"type": ["null", "string", "object"]} for k in SETTINGS},
    },
    "description": "Null QM settings are explicit planning gaps, never computational defaults.",
}
REFERENCE_SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "island.pcff-qm-reference.v1",
    "type": "object",
    "required": [
        "schema",
        "family",
        "target_identity",
        "graph_identity",
        "typing_identity",
        "charge_identity",
        "charge_vector",
        "source_hashes",
        "domain",
        "atom_order",
        "units",
        "files",
        "methods",
        "samples_file",
        "license_file",
        "permission",
        "output_verification_file",
    ],
    "properties": {
        "schema": {"const": "island.pcff-qm-reference.v1"},
        "family": {"const": "pcff"},
        "files": {
            "type": "object",
            "additionalProperties": {"type": "string", "pattern": "^[0-9a-f]{64}$"},
        },
        "methods": {"type": "object", "minProperties": 1},
    },
    "$defs": {
        "sample": {
            "type": "object",
            "required": [
                "id",
                "graph_identity",
                "bond_order_identity",
                "atom_order",
                "coordinates",
                "charge_state",
                "method_id",
                "split",
            ],
            "anyOf": [
                {"required": [k]}
                for k in ("energy", "forces", "hessian", "frequencies")
            ],
        },
        "output_verification": {
            "type": "object",
            "required": [
                "importer_file",
                "importer_sha256",
                "implementation_version",
                "records",
            ],
        },
    },
}
