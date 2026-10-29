"""Evidence-only PDMS supplementation from one pinned, unmodified FRC.

This is neither the Sun paper model nor an executable external bundle. Exact
source decimals may close transcription gaps; absent couplings stay absent.
"""

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.exceptions import IslandError, PCFFError

from .external_pdms import LIMITS, digest, requirements
from .pcff.automatic import assign_automatic_pcff_charges, type_pcff_atoms
from .pcff.class2 import FAMILIES, inspect_pcff_class2, resolve
from .pcff.fallbacks import MSI_POLICY
from .pcff.mixed_sip import DISTRIBUTION_PINS
from .pcff.source import records
from .pdms_literature import INPUT_ERRORS, pinned_bytes

SCHEMA = "island_pdms_source_completion_data_v2"
DATA_ID = "island_external_pdms_pcff_source_supplement_v2"
POLICY = "literal_then_ordinary_family_equivalence_no_wildcards_no_lower_forms_v1"


class PDMSCompletionError(IslandError):
    def __init__(self, status, stage):
        self.status = status
        self.diagnostics = [{"stage": stage}]
        super().__init__(status + ": " + stage)


class PDMSCompletionBlockedError(IslandError):
    def __init__(self, result):
        self.result = result
        self.status = result.payload["status"]
        self.diagnostics = result.payload["missing_dependencies"]
        super().__init__(self.status + ": complete PDMS source is still required")


def source_row(row, source_hash):
    """Retain whole rows, including explicit zeros, with identity conversions."""
    arity = FAMILIES[row["section"]][0]
    words = row["raw"].split("!", 1)[0].split()
    values = words[2 + arity :]
    units = row["source_units"]
    if len(values) != len(units):
        units = units[: len(values)]
    result = {
        "id": row["id"],
        "source_file": "pcff.frc",
        "source_sha256": source_hash,
        "line": row["line"],
        "raw": row["raw"],
        "family": row["section"],
        "namespace": row["namespace"],
        "version": row["version"],
        "reference": row["reference"],
        "types": row["types"],
        "fields": row["fields"][: len(values)],
        "original_values": values,
        "units": units,
        "converted_values": values,
        "conversion": "identity_source_decimal_v1; no mixing with paper values",
        "abbreviated_symmetric_row": row["expansion"],
        "explicit_zero_columns": [i for i, v in enumerate(values) if Decimal(v) == 0],
    }
    return {**result, "evidence_checksum": digest(result)}


def _assess(system, distribution):
    distribution.validate_integrity()
    try:
        required = requirements(system)
    except (PCFFError, *INPUT_ERRORS) as error:
        raise PDMSCompletionError(
            "typing_blocked", "final_graph_invalid:" + str(error)
        ) from error
    source = distribution.source
    typing = type_pcff_atoms(system, source, profile="island_pcff_source_graph_v1")
    if not typing.complete:
        raise PDMSCompletionError("typing_blocked", "source_graph_typing_incomplete")
    tp = typing.payload
    labels = tp["assignments"]
    catalog = inspect_pcff_class2(source)["records"]
    eq = records(source.inventory, "equivalence")
    queries, interactions, used, gaps = {}, {}, set(), []
    for key, request in required["interactions"].items():
        family, sites = request["family"], request["sites"]
        supplied = [labels[i] for i in sites]
        query_key = family + ":" + ",".join(supplied)
        if query_key not in queries:
            decision = resolve(family, supplied, catalog, eq)
            queries[query_key] = decision
            used.update(r["record_id"] for r in decision.get("selected", []))
        decision = queries[query_key]
        interactions[key] = {
            **request,
            "query": query_key,
            "status": decision["status"],
        }
        if decision["status"] != "assigned":
            gaps.append(
                {
                    "interaction": key,
                    "query": query_key,
                    "sites": sites,
                    "types": supplied,
                    "ordinary_types": decision["resolved_types"],
                    "family": family,
                    "status": decision["status"],
                    "reason": decision["reason"],
                    "dependencies": request["dependencies"],
                }
            )
    for request in interactions.values():
        request["dependencies_complete"] = all(
            interactions[d]["status"] == "assigned" for d in request["dependencies"]
        )
    charges = assign_automatic_pcff_charges(
        system, typing, resolution_policy=MSI_POLICY
    )
    # Native increments are retained as a diagnostic, never passed to an engine.
    charge_record = charges.payload["native_charge_record"]
    if not charges.complete:
        gaps.append({"family": "bond_increments", "reason": "native_charge_incomplete"})
    counts = {}
    for request in interactions.values():
        family = request["family"]
        counts.setdefault(family, Counter())[request["status"]] += 1
    rows = {
        rid: source_row(catalog[rid], distribution.hashes["pcff.frc"])
        for rid in sorted(used)
    }
    reference_ids = {r["reference"] for r in rows.values()}
    source_references = {
        s["namespace"]: {"line": s["line"], "raw_lines": s["lines"]}
        for s in source.inventory["sections"]
        if s["name"] == "reference" and s["namespace"] in reference_ids
    }
    # Conventions below describe a converter oracle, not recovered paper-model
    # authority. Missing source constants cannot be replaced by initialization.
    unresolved = [
        "no complete independently qualified PDMS execution bundle",
        "paper 1997 and pinned FRC differ; no numerical overlay is authorized",
        "torsion/angle-angle-torsion equation and unit mapping needs numerical qualification",
        "1-4 and mixing converter conventions are comparisons, not evidence for the 1997 paper model",
        "empty torsion-torsion section does not establish physical inapplicability",
    ]
    p = {
        "schema": SCHEMA,
        "data_identity": DATA_ID,
        "status": "parameter_incomplete",
        "source_family": "pcff_single_pinned_source_data_only",
        "source_hashes": distribution.hashes,
        "selection_policy": POLICY,
        "requirements": required,
        "typing_identity": typing.identity,
        "typing": tp,
        "source_rows": rows,
        "source_references": source_references,
        "queries": queries,
        "interactions": interactions,
        "inapplicable": required["inapplicable"],
        "native_charge_identity": charges.identity,
        "native_charge_record": charge_record,
        "exact_native_charge_vector": charges.charges if charges.complete else None,
        "charge_origin": "same pinned FRC bond increments; zero atomic base, no normalization",
        "coverage": {f: dict(c) for f, c in sorted(counts.items())},
        "missing_dependencies": gaps,
        "unresolved_model_dependencies": unresolved,
        "profile_identity": None,
        "parameter_assignment": None,
        "execution_authorized": False,
        "native_pcff_completion": False,
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    return {**p, "identity": digest(p)}


@dataclass(frozen=True)
class PDMSCompletionRecord:
    """Immutable data assessment. Replay, not a self-signed checksum, is authority."""

    json_text: str

    def __post_init__(self):
        try:
            p = self.payload
            valid = (
                p["schema"] == SCHEMA
                and p["data_identity"] == DATA_ID
                and p["source_hashes"] == DISTRIBUTION_PINS
                and p["selection_policy"] == POLICY
                and p["status"] == "parameter_incomplete"
                and p["identity"]
                == digest({k: v for k, v in p.items() if k != "identity"})
                and p["execution_authorized"] is False
                and p["acceptance_gate"] == 1
                and p["profile_identity"] is None
                and p["parameter_assignment"] is None
                and p["native_pcff_completion"] is False
                and p["limitations"] == LIMITS
            )
        except INPUT_ERRORS as error:
            raise PDMSCompletionError(
                "source_incompatible", "malformed_source_record"
            ) from error
        if not valid:
            raise PDMSCompletionError(
                "source_incompatible", "changed_source_record_or_authority"
            )

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    def validate_integrity(self, system, distribution):
        if self.payload != assess_pdms_source_completion(system, distribution).payload:
            raise PDMSCompletionError(
                "source_incompatible", "graph_rows_charges_or_evidence_changed"
            )

    def require_executable(self):
        raise PDMSCompletionBlockedError(self)


def assess_pdms_source_completion(system, distribution):
    """Validate final graph and audit one complete source, without preparing MD."""
    try:
        return PDMSCompletionRecord(pack(_assess(system, distribution)))
    except PDMSCompletionError:
        raise
    except (PCFFError, *INPUT_ERRORS) as error:
        raise PDMSCompletionError(
            "source_incompatible", "source_assessment_failed:" + str(error)
        ) from error


def load_pdms_source_completion(path, expected_sha256, *, system, distribution):
    """Explicit local path and independent pin; relocation reruns source selection."""
    try:
        raw = pinned_bytes(Path(path), expected_sha256)
        record = PDMSCompletionRecord(raw.decode())
        record.validate_integrity(system, distribution)
        return record
    except PDMSCompletionError:
        raise
    except FileNotFoundError as error:
        raise PDMSCompletionError(
            "asset_unavailable", "source_record_missing"
        ) from error
    except (PCFFError, *INPUT_ERRORS) as error:
        raise PDMSCompletionError(
            "source_incompatible", "source_record_pin_or_replay_failed"
        ) from error


def file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()
