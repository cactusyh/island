"""Versioned source-bundle vocabulary and strict offline metadata validation."""

import re
from datetime import datetime

from island.pcff_pdms.extension import require
from island.pcff_pdms.pilot_intake.gates import CONVENTIONS

VERSION = "island.j81-source-bundle.v1"
PARSER_VERSION = "island.j81-parser.v1"
SOURCE_TYPES = (
    "primary_pdf",
    "original_pcff_frc",
    "materials_studio_export",
    "parameter_table",
    "authorized_archive",
)
STATUSES = (
    "source_missing",
    "source_present_unparsed",
    "source_parsed",
    "literature_candidate_verified",
    "literature_candidate_convention_blocked",
    "exact_native_row",
    "inferred_or_analog_value_blocked",
    "end_group_scope_blocked",
)
ARTIFACT_FIELDS = (
    "id",
    "source_path",
    "source_type",
    "sha256",
    "file_size",
    "bibliography",
    "source_version",
    "acquisition_method",
    "authorization_status",
    "license_or_redistribution_status",
    "extraction_tool",
    "extraction_version",
    "extraction_timestamp",
    "parser_version",
    "authority",
    "source_family",
)
ROW_FIELDS = (
    "id",
    "source_id",
    "source_location",
    "raw_text",
    "atom_type_sequence",
    "element_sequence",
    "parameter_family",
    "numerical_values",
    "units",
    "functional_form",
    "mixing_rule",
    "exclusion_and_1_4_policy",
    "charge_convention",
    "torsion_convention",
    "chain_scope",
    "reporting_mode",
    "chemical_environment",
    "native_identity",
)


def schema():
    artifact = {k: {"type": "string", "minLength": 1} for k in ARTIFACT_FIELDS}
    artifact.update(
        source_type={"enum": list(SOURCE_TYPES)},
        sha256={"type": "string", "pattern": "^[0-9a-f]{64}$"},
        file_size={"type": "integer", "minimum": 0},
        source_path={"type": "string", "minLength": 1},
        bibliography={
            "type": "object",
            "required": ["title", "authors", "year", "identifier"],
            "properties": {
                "title": {"type": "string", "minLength": 1},
                "authors": {
                    "type": "array",
                    "minItems": 1,
                    "items": {"type": "string", "minLength": 1},
                },
                "year": {"type": "integer", "minimum": 1800, "maximum": 2200},
                "identifier": {"type": ["string", "null"]},
            },
            "additionalProperties": False,
        },
        parser_version={"const": PARSER_VERSION},
        source_family={
            "enum": [
                "native_PCFF",
                "published_PDMS",
                "external_PDMS",
                "COMPASS",
                "OPLS",
                "GAFF",
                "CVFF",
                "CFF91",
                "analog",
                "unknown",
            ]
        },
        authority={
            "enum": [
                "original_publisher",
                "original_author",
                "audited_native_release",
                "user_attested",
                "unknown",
            ]
        },
        authorization_status={
            "enum": [
                "authorized_local_use",
                "public_access",
                "unknown",
                "not_authorized",
            ]
        },
        extraction_timestamp={"type": "string", "format": "date-time"},
    )
    row = {k: {} for k in ROW_FIELDS}
    row.update(
        id={"type": "string", "minLength": 1},
        source_id={"type": "string", "minLength": 1},
        source_location={"type": "object"},
        native_identity={"type": "object"},
        chemical_environment={"type": "object"},
        parameter_family={"type": "string", "minLength": 1},
        raw_text={"type": "string", "minLength": 1},
        atom_type_sequence={"type": "array", "items": {"type": "string"}},
        element_sequence={"type": "array", "items": {"type": "string"}},
        numerical_values={"type": "array", "minItems": 1, "items": {"type": "string"}},
        reporting_mode={
            "enum": ["directly_reported", "reconstructed", "inferred", "unknown"]
        },
        chain_scope={"enum": ["internal", "capped_end_group", "both", "unknown"]},
    )
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": VERSION,
        "type": "object",
        "additionalProperties": False,
        "required": ["schema", "artifacts", "candidates", "convention_contract"],
        "properties": {
            "schema": {"const": VERSION},
            "artifacts": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": list(ARTIFACT_FIELDS),
                    "properties": artifact,
                    "additionalProperties": False,
                },
            },
            "candidates": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": list(ROW_FIELDS),
                    "properties": row,
                    "additionalProperties": False,
                },
            },
            "convention_contract": {
                "type": "object",
                "required": ["fields"],
                "properties": {
                    "fields": {
                        "type": "object",
                        "additionalProperties": False,
                        "properties": {
                            k: {
                                "type": "object",
                                "required": ["value", "evidence_refs"],
                                "properties": {
                                    "value": {"type": ["object", "null"]},
                                    "evidence_refs": {
                                        "type": "array",
                                        "items": {"type": "string"},
                                    },
                                },
                            }
                            for k in CONVENTIONS
                        },
                    }
                },
            },
        },
    }


def validate(value, rule=None, location="bundle"):
    """Small deterministic validator for the published schema; no remote refs."""
    rule = schema() if rule is None else rule
    types = {
        "object": dict,
        "array": list,
        "string": str,
        "integer": int,
        "null": type(None),
    }
    wanted = rule.get("type")
    if wanted:
        options = wanted if isinstance(wanted, list) else [wanted]
        require(
            any(type(value) is types[t] for t in options), "schema_type:" + location
        )
    if "const" in rule:
        require(value == rule["const"], "schema_constant:" + location)
    if "enum" in rule:
        require(value in rule["enum"], "schema_enum:" + location)
    if type(value) is dict:
        require(
            set(rule.get("required", [])) <= value.keys(), "schema_missing:" + location
        )
        properties = rule.get("properties", {})
        if rule.get("additionalProperties") is False:
            require(set(value) <= properties.keys(), "schema_unknown:" + location)
        for key in value:
            validate(value[key], properties.get(key, {}), location + "." + key)
    if type(value) is list:
        require(len(value) >= rule.get("minItems", 0), "schema_empty:" + location)
        for i, item in enumerate(value):
            validate(item, rule.get("items", {}), location + "." + str(i))
    if type(value) is str:
        require(len(value) >= rule.get("minLength", 0), "schema_empty:" + location)
        require(
            not rule.get("pattern") or re.fullmatch(rule["pattern"], value),
            "schema_pattern:" + location,
        )
        if rule.get("format") == "date-time":
            require(
                datetime.fromisoformat(value).tzinfo is not None,
                "timezone_required",
            )
    if type(value) is int:
        require(
            rule.get("minimum", value) <= value <= rule.get("maximum", value),
            "schema_range:" + location,
        )
