"""Portable JSON Schema for a reviewable, externally pinned source package."""

import json

from island.opls_conventions import REQUIRED_SETTINGS

VERSION = "island.opls-source-intake.v1"
STATUSES = (
    "not_supplied",
    "supplied_unreadable",
    "supplied_incomplete",
    "supplied_conflicting",
    "source_complete_not_executed",
    "source_complete_executed",
    "source_complete_independently_verified",
)
REQUIRED_CONVENTIONS = tuple(
    sorted(
        set(REQUIRED_SETTINGS)
        | {
            "parameter_functional_form",
            "charge_method",
            "exclusions",
        }
    )
)


def obj(properties):
    return {
        "type": "object",
        "properties": properties,
        "required": list(properties),
        "additionalProperties": False,
    }


TEXT = {"type": "string", "minLength": 1}
HASH = {"type": "string", "pattern": "^[0-9a-f]{64}$"}
DOI = {"type": "string", "pattern": r"^10\.\d{4,9}/\S+$"}
EVIDENCE = obj(
    {
        "file": TEXT,
        "sha256": HASH,
        "location": TEXT,
        "first_line": {"type": "integer", "minimum": 1},
        "raw_text": TEXT,
        "doi": DOI,
        "citation": TEXT,
    }
)
CONVERSION = obj({"formula": TEXT, "implementation": TEXT, "version": TEXT})
RECORD = obj(
    {
        "key": TEXT,
        "family": {"const": "opls"},
        "labels": {"type": "array", "items": TEXT},
        "functional_form": TEXT,
        "resolution": {"enum": ["source_resolved", "converter_resolved"]},
        "value": TEXT,
        "original_value": TEXT,
        "units": TEXT,
        "conversion": CONVERSION,
        "evidence": {"type": "array", "items": EVIDENCE, "minItems": 1},
        "reconstruction": obj(
            {
                "code": EVIDENCE,
                "result": EVIDENCE,
                "implementation": TEXT,
                "version": TEXT,
                "value": TEXT,
                "units": TEXT,
            }
        ),
        "evidence_checksum": HASH,
    }
)
DOMAIN = obj(
    {
        "name": TEXT,
        "system_file": TEXT,
        "graph_identity": HASH,
        "typed_graph_identity": HASH,
        "atom_types": {"type": "array", "items": TEXT, "minItems": 1},
        "charges": {"type": "array", "items": {"type": "number"}, "minItems": 1},
        "charge_origin": TEXT,
        "charge_method": TEXT,
        "charge_units": {"const": "elementary_charge"},
        "rows": {"type": "array", "items": RECORD, "minItems": 1},
    }
)
SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "urn:island:opls-source-intake:v1",
    "x-required-conventions": list(REQUIRED_CONVENTIONS),
    "x-source-statuses": list(STATUSES),
    "x-qualification-policy": {
        "source_complete_not_executed": "complete contract; no execution or registration",
        "source_complete_executed": "requires a pinned executable qualification adapter",
        "source_complete_independently_verified": "also requires independent numerical and workflow reconstruction",
        "package_status_assertions": "forbidden; status is computed by ISLAND",
    },
    "description": (
        "Structural source review only; completeness does not register a profile or "
        "authorize execution. Manifest bytes require an independent SHA-256 pin. "
        "All other physical files, including extracted evidence and code, are inventoried."
    ),
    **obj(
        {
            "schema": {"const": VERSION},
            "family": {"const": "opls"},
            "source_kind": {"const": "external_literature"},
            "proposed_profile": TEXT,
            "source_version": TEXT,
            "doi": DOI,
            "citation": TEXT,
            "files": {
                "type": "array",
                "minItems": 1,
                "items": obj(
                    {
                        "path": TEXT,
                        "size": {"type": "integer", "minimum": 0},
                        "sha256": HASH,
                        "role": {
                            "enum": [
                                "paper",
                                "supplement",
                                "author_parameters",
                                "author_input",
                                "author_data",
                                "converter",
                                "version",
                                "permission",
                                "extracted_evidence",
                                "independent_reconstruction",
                                "system",
                            ]
                        },
                    }
                ),
            },
            "authorization": obj(
                {
                    "status": {"enum": ["permitted", "unknown", "prohibited"]},
                    "scope": {"const": "local_research"},
                    "rights_holder": TEXT,
                    "redistribution": {
                        "enum": ["permitted", "not_established", "prohibited"]
                    },
                    "evidence": EVIDENCE,
                    "covered_files": {
                        "type": "array",
                        "minItems": 1,
                        "items": obj({"file": TEXT, "sha256": HASH}),
                    },
                }
            ),
            "version_evidence": EVIDENCE,
            "conventions": {"type": "array", "minItems": 1, "items": RECORD},
            "domains": {"type": "array", "minItems": 1, "items": DOMAIN},
        }
    ),
}

_SCHEMA_JSON = json.dumps(SCHEMA, sort_keys=True, allow_nan=False)


def source_intake_schema():
    """An owned copy; consumers cannot mutate the runtime intake contract."""
    return json.loads(_SCHEMA_JSON)
