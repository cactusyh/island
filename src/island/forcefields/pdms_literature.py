"""Explicit local literature ingestion; data integrity never grants MD authority.

No discovery, PDF inference, parameter matching, conversion or numerical execution
occurs here. A caller supplies a pinned source manifest and, separately, any
graph-bound C5 bundle. The C5 contract and native force fields are unchanged.
"""

import re
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.dynamics._checkpoint_data import strict_load
from island.exceptions import IslandError
from island.workflows import storage

from .external_pdms import (
    LIMITS,
    PROFILE,
    ROWS_SCHEMA,
    UNITS,
    load_external_pdms_bundle,
    number,
    require,
    source_row,
)
from .external_pdms_adapter import FAMILY, ExternalPDMSOptions
from .preparation import ForceFieldRequest, prepare_forcefield

SOURCE_SCHEMA = "island_pdms_literature_source_v1"
RECORD_SCHEMA = "island_pdms_literature_record_v1"
RESULT_SCHEMA = "island_pdms_literature_ingestion_v1"
CLASSIFICATIONS = (
    "pdms_specific_external",
    "compass_class2_external",
    "pcff_extension_external",
    "opls_compatible_external",
    "incomplete_literature_source",
    "unusable_or_unlicensed",
)
SOURCE_FIELDS = {
    "schema",
    "bibliography",
    "classification",
    "source_family",
    "model_name",
    "license",
    "files",
    "primary_source",
    "parameterization",
    "functional_form",
    "units",
    "charge_method",
    "atom_type_definitions",
    "reported_validation_properties",
}
INPUT_ERRORS = (
    ValueError,
    TypeError,
    KeyError,
    AttributeError,
    IndexError,
    OSError,
    UnicodeError,
    ArithmeticError,
    IslandError,
)


@dataclass(frozen=True)
class LiteratureSourceRecord:
    """Owned source inventory with exact bytes, provenance and extraction results."""

    json_text: str

    def __post_init__(self):
        p = self.payload
        require(p["schema"] == RECORD_SCHEMA, "wrong_literature_record_schema")
        require(p["classification"] in CLASSIFICATIONS, "unknown_source_classification")
        require(p["execution_authorized"] is False, "source_cannot_grant_execution")

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return sha256(self.json_text.encode()).hexdigest()


@dataclass(frozen=True)
class LiteratureIngestion:
    """An immutable report, never an execution profile or parameter assignment."""

    json_text: str

    def __post_init__(self):
        p = self.payload
        require(p["schema"] == RESULT_SCHEMA, "wrong_ingestion_schema")
        require(
            p["status"]
            in (
                "asset_unavailable",
                "parameter_incomplete",
                "scientifically_unqualified",
            )
            and p["acceptance_gate"] == 1
            and p["execution_profile"] is None
            and p["backend_executions"] == []
            and p["limitations"] == LIMITS,
            "literature_cannot_promote_readiness",
        )
        if p["source_record"] is not None:
            LiteratureSourceRecord(pack(p["source_record"]))

    @property
    def payload(self):
        return unpack(self.json_text)


def pinned_bytes(path, expected):
    """Check an independent pin before parsing or interpreting any asset."""
    require(
        type(expected) is str and re.fullmatch(r"[0-9a-f]{64}", expected),
        "independent_sha256_required",
    )
    raw = Path(path).read_bytes()
    require(sha256(raw).hexdigest() == expected, "asset_sha256_mismatch")
    return raw


def reconstruct_rows(document, assets):
    """Independently extract decimal tokens, then compare with the C5 reader.

    PDF bytes are not text tables. Suppliers must provide all numerical rows as
    independently cited text; no figures or neighboring rows are interpreted.
    """
    require(
        type(document["rows"]) is dict and document["rows"], "numerical_rows_missing"
    )
    extracted, errors = {}, []
    for rid, row in document["rows"].items():
        try:
            citations, values = [], []
            for c in row["citations"]:
                raw = assets[c["file"]]
                require(not raw.startswith(b"%PDF"), "pdf_not_reconstructable_table")
                lines = raw.decode("utf-8").splitlines()
                require(
                    type(c["line"]) is int and 1 <= c["line"] <= len(lines),
                    "one_based_line_required",
                )
                text = lines[c["line"] - 1]
                require(text == c["raw"], "raw_citation_changed")
                require(
                    type(c["columns"]) is list
                    and bool(c["columns"])
                    and all(type(i) is int and i >= 0 for i in c["columns"]),
                    "numeric_columns_required",
                )
                tokens = text.split("!", 1)[0].split()
                decimals = [tokens[i] for i in c["columns"]]
                values.extend(number(v) for v in decimals)
                citations.append(
                    {
                        "source_file": c["file"],
                        "line": c["line"],
                        "raw_text": text,
                        "numeric_columns": c["columns"],
                        "decimal_values": decimals,
                    }
                )
            require(
                bool(citations) and values == [number(v) for v in row["values"]],
                "extracted_decimal_mismatch",
            )
            require(
                values == source_row(row, document, assets),
                "independent_reader_mismatch",
            )
            extracted[rid] = {"row": row, "citations": citations}
        except INPUT_ERRORS as error:
            errors.append(
                {"stage": "source_row_incomplete", "row": rid, "context": str(error)}
            )
    return extracted, errors


def inspect_source(path, pin):
    """Read only the explicit manifest and its hash-listed relative dependencies."""
    diagnostics, assets, extracted, document, manifest = [], {}, {}, None, None
    classification = "incomplete_literature_source"
    manifest_hash = None
    try:
        raw = pinned_bytes(path, pin)
        manifest_hash = sha256(raw).hexdigest()
        manifest = strict_load(raw.decode("utf-8"))
        require(
            type(manifest) is dict and set(manifest) == SOURCE_FIELDS,
            "source_manifest_fields_invalid",
        )
        require(manifest["schema"] == SOURCE_SCHEMA, "source_schema_invalid")
        require(manifest["classification"] in CLASSIFICATIONS, "classification_invalid")
        b = manifest["bibliography"]
        require(
            set(b) == {"doi", "title", "authors", "year"}
            and all(type(b[k]) is str and b[k].strip() for k in ("doi", "title"))
            and type(b["authors"]) is list
            and bool(b["authors"])
            and all(type(a) is str and a.strip() for a in b["authors"])
            and type(b["year"]) is int
            and b["year"] > 0,
            "bibliography_incomplete",
        )
        lic = manifest["license"]
        licensed = (
            type(lic) is dict
            and set(lic)
            == {"identifier", "authorized_use", "redistribution", "evidence_file"}
            and type(lic["identifier"]) is str
            and bool(lic["identifier"].strip())
            and lic["authorized_use"] is True
            and lic["redistribution"] in ("allowed", "local_only")
        )
        if not licensed or manifest["classification"] == "unusable_or_unlicensed":
            classification = "unusable_or_unlicensed"
            raise ValueError("license_or_authorization_missing")
        family = manifest["source_family"]
        require(
            type(family) is str and family.startswith("external:") and len(family) > 9,
            "external_family_required",
        )
        require(
            type(manifest["model_name"]) is str
            and bool(manifest["model_name"].strip()),
            "model_name_missing",
        )
        require(
            manifest["functional_form"] == "class2" and manifest["units"] == UNITS,
            "functional_form_or_units_unsupported",
        )
        method = manifest["parameterization"]
        require(
            set(method) == {"method", "quantum_chemistry", "condensed_phase_targets"}
            and type(method["method"]) is str
            and bool(method["method"].strip())
            and type(method["condensed_phase_targets"]) is list
            and (
                method["quantum_chemistry"] is None
                or type(method["quantum_chemistry"]) is str
            ),
            "parameterization_metadata_missing",
        )
        require(
            type(manifest["charge_method"]) is str
            and bool(manifest["charge_method"].strip())
            and type(manifest["atom_type_definitions"]) is dict
            and bool(manifest["atom_type_definitions"])
            and all(
                type(v) is str and bool(v.strip())
                for v in manifest["atom_type_definitions"].values()
            )
            and type(manifest["reported_validation_properties"]) is list,
            "chemical_provenance_missing",
        )
        files = manifest["files"]
        require(type(files) is dict and bool(files), "source_files_missing")
        for name, f in files.items():
            try:
                require(
                    set(f) == {"sha256", "family", "role", "provenance"}
                    and f["family"] == family,
                    "mixed_source_or_file_fields_invalid",
                )
                require(
                    type(f["provenance"]) is str
                    and bool(f["provenance"].strip())
                    and f["role"]
                    in (
                        "paper",
                        "supplement",
                        "raw_table",
                        "row_table",
                        "code",
                        "license",
                        "benchmark",
                    ),
                    "file_provenance_invalid",
                )
                assets[name] = pinned_bytes(
                    storage.child(Path(path).parent, name), f["sha256"]
                )
            except INPUT_ERRORS as error:
                diagnostics.append(
                    {
                        "stage": "source_dependency_invalid",
                        "file": name,
                        "context": str(error),
                    }
                )
        require(
            lic["evidence_file"] in assets
            and files[lic["evidence_file"]]["role"] == "license"
            and bool(assets[lic["evidence_file"]].strip()),
            "license_evidence_missing",
        )
        primary = manifest["primary_source"]
        require(
            primary in assets and files[primary]["role"] == "row_table",
            "reconstructable_parameter_table_missing",
        )
        document = strict_load(assets[primary].decode("utf-8"))
        require(
            set(document) == {"schema", "source_family", "rows", "provenance"}
            and document["schema"] == ROWS_SCHEMA
            and document["source_family"] == family
            and type(document["provenance"]) is str
            and bool(document["provenance"].strip()),
            "source_row_document_invalid",
        )
        extracted, errors = reconstruct_rows(document, assets)
        diagnostics.extend(errors)
        if not diagnostics:
            classification = manifest["classification"]
    except INPUT_ERRORS as error:
        diagnostics.append(
            {"stage": "literature_source_incomplete", "context": str(error)}
        )
    return LiteratureSourceRecord(
        pack(
            {
                "schema": RECORD_SCHEMA,
                "manifest_sha256": manifest_hash,
                "manifest": manifest,
                "classification": classification,
                "source_rows": document.get("rows") if type(document) is dict else None,
                "extracted_rows": extracted,
                "diagnostics": diagnostics,
                "execution_authorized": False,
            }
        )
    )


def ingest_literature(
    system, *, source=None, source_sha256=None, bundle=None, bundle_sha256=None
):
    """Fail closed on absent/incomplete literature, provenance and graph coverage.

    A data-valid C5 import still uses J42's blocked preparation adapter. There is
    no qualified numerical path without independently executed acceptance.
    """
    preparation = prepare_forcefield(
        system, ForceFieldRequest(FAMILY, ExternalPDMSOptions())
    )
    p = preparation.payload
    status, data_status = "asset_unavailable", "incomplete"
    diagnostics, source_record, validation = [], None, None
    if source is None:
        diagnostics.append({"stage": "explicit_literature_source_required"})
    else:
        source_record = inspect_source(source, source_sha256).payload
        diagnostics.extend(source_record["diagnostics"])
        status = "parameter_incomplete"
        if not diagnostics and bundle is not None:
            try:
                pinned_bytes(bundle, bundle_sha256)
                candidate, outcome = load_external_pdms_bundle(
                    bundle, system, expected_sha256=bundle_sha256
                )
                validation = outcome.payload
                diagnostics.extend(validation["diagnostics"])
                if outcome.complete:
                    c, m = candidate.payload, source_record["manifest"]
                    require(
                        source_record["classification"] in CLASSIFICATIONS[:4],
                        "incomplete_source_cannot_qualify_bundle",
                    )
                    require(
                        c["source_family"] == m["source_family"]
                        and c["files"] == m["files"]
                        and c["primary_source"] == m["primary_source"]
                        and c["source_rows"]
                        == {
                            k: v["row"]
                            for k, v in source_record["extracted_rows"].items()
                        },
                        "mixed_or_uncited_bundle_source",
                    )
                    require(
                        all(
                            a["type"] in m["atom_type_definitions"] for a in c["atoms"]
                        ),
                        "atom_type_definition_missing",
                    )
                    preparation = prepare_forcefield(
                        system,
                        ForceFieldRequest(
                            FAMILY, ExternalPDMSOptions(bundle, bundle_sha256)
                        ),
                    )
                    status, data_status = "scientifically_unqualified", "data_valid"
            except INPUT_ERRORS as error:
                diagnostics.append(
                    {"stage": "literature_bundle_incomplete", "context": str(error)}
                )
        elif not diagnostics:
            diagnostics.append({"stage": "explicit_graph_bound_bundle_required"})
    if p["status"] == "typing_blocked":
        status = "parameter_incomplete"
        diagnostics.extend(p["diagnostics"])
    return LiteratureIngestion(
        pack(
            {
                "schema": RESULT_SCHEMA,
                "profile_candidate": PROFILE,
                "status": status,
                "data_status": data_status,
                "source_record": source_record,
                "bundle_validation": validation,
                "preparation": preparation.payload,
                "missing_interactions": preparation.payload["missing_interactions"],
                "unresolved_dependencies": preparation.payload[
                    "unresolved_dependencies"
                ],
                "diagnostics": diagnostics,
                "execution_profile": None,
                "backend_executions": [],
                "native_pcff_acceptance_gate": 1,
                "acceptance_gate": 1,
                "limitations": LIMITS,
            }
        )
    )
