"""Audited Xiang Table S4 data, isolated from native Foyer OPLS-AA.

This is a queryable literature supplement, never an executable parameter source.
J44/J45 evidence stays authoritative and unchanged. Units, terminal assignment
and pair policies are intentionally unresolved; no conversion is performed.
"""

from dataclasses import dataclass
from decimal import Decimal

from island.charge_references.records import pack, unpack
from island.exceptions import IslandError
from island.forcefields.external_pdms import LIMITS, digest
from island.forcefields.literature_closure import (
    CATALOG_SHA,
    SELECTION_SHA,
    DomainSelection,
    close_domain,
)
from island.forcefields.pdms_literature import INPUT_ERRORS, pinned_bytes
from island.workflows.storage import child

SCHEMA = "island_external_opls_pdms_literature_supplement_v1"
DATA_ID = "island_external_opls_compatible_pdms_xiang2025_supplement_v1"
SOURCE_FILES = ("pcff/jp4c08471.pdf", "pcff/jp4c08471_si_001.pdf")
KINDS = {
    "pair_coeff": "nonbonded",
    "bond_coeff": "bond",
    "angle_coeff": "angle",
    "dihedral_coeff": "proper_torsion",
}


class OPLSPDMSLiteratureError(IslandError):
    """Classified asset/integrity error with the complete source diagnostics."""

    def __init__(self, status, stage, diagnostics=()):
        self.status = status
        self.diagnostics = list(diagnostics) or [{"stage": stage}]
        super().__init__(status + ": " + stage)


def _row(kind, labels, values, columns, evidence):
    # Numeric tokens are source decimals, not unit-converted coefficients.
    row = {
        "kind": kind,
        "literal_labels": labels,
        "original_values": values,
        "numeric_columns_zero_based": columns,
        "source_evidence": evidence,
        "source_units": None,
        "converted_values": None,
        "conversion_code_version": None,
        "explicit_zero_value_indices": [
            i for i, value in enumerate(values) if Decimal(value) == 0
        ],
        "accepted_assignment": None,
    }
    row["evidence_checksum"] = digest(row)
    return row


def _document(selection):
    if type(selection) is not DomainSelection:
        raise OPLSPDMSLiteratureError(
            "source_incompatible", "pinned_j45_selection_required"
        )
    # Recheck independent J44/J45 pins, including re-signed evidence mutations.
    selection = DomainSelection(selection.json_text, selection.catalog_text)
    p = selection.payload
    rows = []
    for evidence in p["candidate_rows"]:
        command, *tokens = evidence["raw_text"].split()
        labels = tokens[:2] if command == "pair_coeff" else tokens[0].split("_")
        rows.append(
            _row(
                KINDS[command],
                labels,
                evidence["original_values"],
                evidence["numeric_columns"],
                evidence,
            )
        )
    for evidence in p["citations"]:
        tokens = evidence["raw_text"].split()
        if (
            evidence["source_file"] == SOURCE_FILES[1]
            and evidence["page"] == 6
            and len(tokens) == 3
            and tokens[0] in ("Si", "O", "C", "H")
        ):
            row = _row("partial_charge", tokens[:1], tokens[1:], [1, 2], evidence)
            row["column_names"] = ["inside monomer", "end monomer"]
            row["charge_method_assignment"] = None
            row["evidence_checksum"] = digest(
                {k: v for k, v in row.items() if k != "evidence_checksum"}
            )
            rows.append(row)
    if len(rows) != 21:
        raise OPLSPDMSLiteratureError(
            "source_incompatible", "table_s4_evidence_coverage_changed"
        )
    files = {f["filename"]: f for f in selection.catalog["files"]}
    result = {
        "schema": SCHEMA,
        "data_identity": DATA_ID,
        "source_family": "external:opls_compatible_pdms",
        "j44_catalog_sha256": CATALOG_SHA,
        "j45_selection_sha256": SELECTION_SHA,
        "sources": {name: files[name] for name in SOURCE_FILES},
        "rows": rows,
        "unresolved_policies": p["required_unresolved_policies"],
        "status": "parameter_incomplete",
        "profile_identity": None,
        "parameter_assignment": None,
        "exact_charge_vector": None,
        "execution_authorized": False,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    result["identity"] = digest(result)
    return result


@dataclass(frozen=True)
class OPLSPDMSLiteratureSupplement:
    """Immutable data-only view; payload copies cannot mutate pinned evidence."""

    json_text: str
    selection: DomainSelection

    def __post_init__(self):
        try:
            valid = pack(self.payload) == pack(_document(self.selection))
        except INPUT_ERRORS as error:
            raise OPLSPDMSLiteratureError(
                "source_incompatible", "malformed_or_mutated_opls_pdms_evidence"
            ) from error
        if not valid:
            raise OPLSPDMSLiteratureError(
                "source_incompatible", "opls_pdms_rows_or_authority_changed"
            )

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    def lookup(self, kind, labels):
        """Literal source-label lookup; no typing, reversal, wildcard or fallback."""
        return tuple(
            row
            for row in self.payload["rows"]
            if row["kind"] == kind and row["literal_labels"] == list(labels)
        )

    def source_diagnostics(self, references):
        errors = []
        for name, source in self.payload["sources"].items():
            try:
                pinned_bytes(child(references, name), source["sha256"])
            except INPUT_ERRORS as error:
                errors.append(
                    {
                        "source_file": name,
                        "expected_sha256": source["sha256"],
                        "status": "asset_unavailable"
                        if isinstance(error, FileNotFoundError)
                        else "source_incompatible",
                        "reason": str(error),
                    }
                )
        return errors

    def verify_sources(self, references):
        errors = self.source_diagnostics(references)
        if errors:
            status = (
                "source_incompatible"
                if any(e["status"] == "source_incompatible" for e in errors)
                else "asset_unavailable"
            )
            raise OPLSPDMSLiteratureError(status, "source_verification_failed", errors)

    def assess(self, system, *, references):
        """Reuse J45's graph/order/charge/dependency audit without new authority."""
        return close_domain(
            system, self.selection, source_errors=self.source_diagnostics(references)
        )


def build_opls_pdms_literature(selection):
    """Expose the historical decimal records; this performs no source discovery."""
    try:
        return OPLSPDMSLiteratureSupplement(pack(_document(selection)), selection)
    except OPLSPDMSLiteratureError:
        raise
    except INPUT_ERRORS as error:
        raise OPLSPDMSLiteratureError(
            "source_incompatible", "invalid_opls_pdms_selection"
        ) from error


def load_opls_pdms_literature(path, expected_sha256, *, selection, references):
    """Load an explicit export with an independent pin and verified local PDFs."""
    try:
        raw = pinned_bytes(path, expected_sha256)
        result = OPLSPDMSLiteratureSupplement(raw.decode("utf-8"), selection)
        result.verify_sources(references)
        return result
    except OPLSPDMSLiteratureError:
        raise
    except FileNotFoundError as error:
        raise OPLSPDMSLiteratureError(
            "asset_unavailable", "explicit_supplement_missing"
        ) from error
    except INPUT_ERRORS as error:
        raise OPLSPDMSLiteratureError(
            "source_incompatible", "supplement_pin_or_record_invalid"
        ) from error
