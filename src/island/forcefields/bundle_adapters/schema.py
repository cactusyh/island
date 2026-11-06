"""Strict local bundle contracts. Hashes establish identity, never authorship."""

import json
from datetime import datetime
from hashlib import sha256
from pathlib import PurePosixPath

from island.forcefields.assignment_core.schema import (
    FAMILIES,
    identity,
    is_hash,
    require,
)

BUNDLE = "island.audited-parameter-bundle.v1"
SOURCE = "island.parameter-source-records.v1"
AUDIT = "island.parameter-bundle-audit.v1"
CAPTURE = "island.j84.captured-bundle.v1"
ADAPTER = "island.j84.canonical-json.v1"
ORIGINS = (
    "vendor_stock",
    "author_primary",
    "authorized_export",
    "authorized_table",
    "synthetic_control",
)
MANIFEST_FIELDS = (
    "schema",
    "bundle_id",
    "profile_id",
    "source_uri",
    "source_sha256",
    "family",
    "version",
    "files",
    "unit_convention",
    "functional_form_convention",
    "mixing_rule",
    "policy_14",
    "charge_convention",
    "provenance",
    "ingestion_timestamp",
    "adapter_version",
    "synthetic",
    "evidence_status",
    "topology_support",
    "identity",
)
PROVENANCE_FIELDS = (
    "citation",
    "authors",
    "year",
    "persistent_id",
    "license",
    "redistribution",
    "acquisition_method",
    "authorization",
    "origin",
    "upstream_identity",
)
FILE_FIELDS = ("path", "sha256", "size", "role", "family", "version")


def keys(value, expected, reason):
    require(isinstance(value, dict) and set(value) == set(expected), reason)


def text(value, reason):
    require(isinstance(value, str) and bool(value.strip()), reason)


def timestamp(value):
    text(value, "timestamp_missing")
    require(value.endswith("Z"), "timestamp_requires_UTC_Z")
    datetime.fromisoformat(value)


def loads(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate_JSON_key:" + key)
            result[key] = value
        return result

    return json.loads(
        raw,
        object_pairs_hook=pairs,
        parse_constant=lambda _: require(False, "nonfinite_JSON"),
    )


def sha(raw):
    return sha256(raw).hexdigest()


def source_digest(files):
    return identity(sorted((f["path"], f["sha256"], f["size"]) for f in files))


def manifest_identity(manifest):
    return identity({k: v for k, v in manifest.items() if k != "identity"})


def relative_path(value):
    text(value, "source_path_missing")
    path = PurePosixPath(value)
    require(
        not path.is_absolute()
        and ".." not in path.parts
        and str(path) == value
        and "\\" not in value,
        "unsafe_source_path",
    )
    require(value not in (".", ""), "unsafe_source_path")
    return path


def validate_manifest(m):
    keys(m, MANIFEST_FIELDS, "bundle_schema_fields")
    require(
        m["schema"] == BUNDLE and m["adapter_version"] == ADAPTER,
        "unsupported_bundle_adapter",
    )
    for name in (
        "bundle_id",
        "profile_id",
        "source_uri",
        "version",
        "mixing_rule",
        "charge_convention",
    ):
        text(m[name], "missing_bundle_metadata:" + name)
    require(m["family"] in FAMILIES, "unknown_forcefield_family")
    require(type(m["synthetic"]) is bool, "synthetic_flag_requires_boolean")
    timestamp(m["ingestion_timestamp"])
    keys(
        m["unit_convention"],
        ("energy_unit", "length_unit", "angular_unit"),
        "unit_convention_schema",
    )
    keys(m["policy_14"], ("lj", "coulomb"), "policy_14_schema")
    for values in (m["unit_convention"], m["policy_14"]):
        for value in values.values():
            text(value, "missing_convention_value")
    keys(m["provenance"], PROVENANCE_FIELDS, "license_or_provenance_metadata_missing")
    for name in PROVENANCE_FIELDS:
        if name != "year":
            text(m["provenance"][name], "missing_provenance:" + name)
    require(
        type(m["provenance"]["year"]) is int
        and 1800 <= m["provenance"]["year"] <= 9999,
        "publication_year_missing",
    )
    require(isinstance(m["files"], list) and m["files"], "source_file_list_empty")
    for f in m["files"]:
        keys(f, FILE_FIELDS, "source_file_schema")
        relative_path(f["path"])
        require(
            is_hash(f["sha256"]) and type(f["size"]) is int and f["size"] > 0,
            "source_file_hash_or_size_missing",
        )
        require(
            f["role"] in ("record_source", "primary_artifact"),
            "unknown_source_file_role",
        )
        require(
            f["family"] == m["family"] and f["version"] == m["version"],
            "cross_family_or_version_overlay",
        )
    require(
        len({f["path"] for f in m["files"]}) == len(m["files"]), "duplicate_source_path"
    )
    require(
        len({f["sha256"] for f in m["files"]}) == len(m["files"]),
        "duplicate_source_identity",
    )
    require(
        any(f["role"] == "record_source" for f in m["files"]),
        "no_parseable_record_source",
    )
    require(
        source_digest(m["files"]) == m["source_sha256"],
        "aggregate_source_hash_mismatch",
    )
    require(manifest_identity(m) == m["identity"], "manifest_identity_mismatch")
    require(
        isinstance(m["functional_form_convention"], list), "functional_forms_schema"
    )
    require(
        isinstance(m["topology_support"], list)
        and m["topology_support"]
        and set(m["topology_support"])
        <= {"finite_linear", "finite_branched", "finite_crosslinked"},
        "topology_scope_invalid",
    )
    require(
        len(set(m["topology_support"])) == len(m["topology_support"]),
        "duplicate_topology_scope",
    )


def schemas():
    """JSON Schemas for exchange structure; runtime also verifies semantic bindings."""
    string = {"type": "string", "minLength": 1}
    checksum = {"type": "string", "pattern": "^[0-9a-f]{64}$"}
    when = {"type": "string", "format": "date-time", "pattern": "Z$"}

    def obj(properties):
        return {
            "type": "object",
            "properties": properties,
            "required": list(properties),
            "additionalProperties": False,
        }

    def array(item):
        return {"type": "array", "items": item}

    def document(value):
        return {"$schema": "https://json-schema.org/draft/2020-12/schema", **value}

    form = obj(
        {
            "family": string,
            "equation": string,
            "fields": array(string),
            "units": array(string),
            "ordering": string,
            "equilibrium_fields": array(string),
        }
    )
    file = obj(
        {
            "path": string,
            "sha256": checksum,
            "size": {"type": "integer", "minimum": 1},
            "role": {"enum": ["record_source", "primary_artifact"]},
            "family": {"enum": list(FAMILIES)},
            "version": string,
        }
    )
    provenance = obj(
        {
            k: {"type": "integer", "minimum": 1800, "maximum": 9999}
            if k == "year"
            else string
            for k in PROVENANCE_FIELDS
        }
    )
    bundle = document(
        obj(
            {
                "schema": {"const": BUNDLE},
                "bundle_id": string,
                "profile_id": string,
                "source_uri": string,
                "source_sha256": checksum,
                "family": {"enum": list(FAMILIES)},
                "version": string,
                "files": {**array(file), "minItems": 1},
                "unit_convention": obj(
                    {k: string for k in ("energy_unit", "length_unit", "angular_unit")}
                ),
                "functional_form_convention": array(form),
                "mixing_rule": string,
                "policy_14": obj({"lj": string, "coulomb": string}),
                "charge_convention": string,
                "provenance": provenance,
                "ingestion_timestamp": when,
                "adapter_version": {"const": ADAPTER},
                "synthetic": {"type": "boolean"},
                "evidence_status": string,
                "topology_support": {
                    "type": "array",
                    "items": {
                        "enum": [
                            "finite_linear",
                            "finite_branched",
                            "finite_crosslinked",
                        ]
                    },
                    "minItems": 1,
                    "uniqueItems": True,
                },
                "identity": checksum,
            }
        )
    )
    values = {
        "atom_type": obj(
            {
                "id": string,
                "atom_type": string,
                "smarts": string,
                "priority": {"type": "integer"},
                "overrides": array(string),
            }
        ),
        "form": form,
        "convention": obj({"name": string, "value": string}),
        "graph_scope": obj({"graph_identity": checksum}),
        "parameter": obj(
            {
                "id": string,
                "family": string,
                "atom_types": array(string),
                "environment": checksum,
                "form": form,
                "values": array({"type": "string"}),
                "dependency_keys": array(checksum),
                "applicability": {"enum": ["required", "not_applicable"]},
                "absence_reason": {"type": "string"},
                "equilibrium_references": {"type": "array"},
            }
        ),
    }
    row = {
        "oneOf": [
            obj(
                {
                    "record_id": string,
                    "kind": {"const": k},
                    "location": string,
                    "origin": {"const": "directly_reported"},
                    "value": v,
                }
            )
            for k, v in values.items()
        ]
    }
    source = document(
        obj(
            {
                "schema": {"const": SOURCE},
                "family": {"enum": list(FAMILIES)},
                "version": string,
                "records": array(row),
            }
        )
    )
    claim = obj({"record_id": string, "claim_sha256": checksum})
    review = obj(
        {
            "path": string,
            "sha256": checksum,
            "provenance_verified": {"const": True},
            "license_verified": {"const": True},
            "claims": array(claim),
        }
    )
    audit = document(
        obj(
            {
                "schema": {"const": AUDIT},
                "manifest_identity": checksum,
                "reviewer": string,
                "review_timestamp": when,
                "authority_basis": string,
                "policy": {"const": "independent_source_review"},
                "synthetic": {"type": "boolean"},
                "source_reviews": array(review),
            }
        )
    )
    return {"bundle": bundle, "source_records": source, "audit": audit}
