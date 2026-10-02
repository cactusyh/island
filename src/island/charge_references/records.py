"""Owned data-only whole-oligomer AM1-BCC evidence; no charge-transfer engine."""

from dataclasses import dataclass
from math import isfinite
from pathlib import Path

from island.core import MolecularSystem
from island.exceptions import ChargeReferenceError
from island.forcefields.ambertools.models import (
    AmberToolsPreparationResult,
    digest,
    validate_preparation_record,
)
from island.forcefields.charges.models import (
    ChargeAssignment,
    ChargeAssignmentDiagnostic,
    ChargeAssignmentResult,
    ChargeCoverage,
    ComponentChargeDiagnostic,
)
from island.forcefields.typing.signatures import graph_signature
from island.workflows import storage
from island.workflows.bundle import record, system_data, system_from

from .correspondence import repeat_correspondence, require

REFERENCE_SCHEMA = "island_oligomer_charge_reference_v1"


def charge_from(p):
    p = dict(p)
    p["assignments"] = {k: ChargeAssignment(**v) for k, v in p["assignments"].items()}
    p["coverage"] = ChargeCoverage(**p["coverage"])
    p["diagnostics"] = tuple(ChargeAssignmentDiagnostic(**v) for v in p["diagnostics"])
    p["component_diagnostics"] = tuple(
        ComponentChargeDiagnostic(**v) for v in p["component_diagnostics"]
    )
    return ChargeAssignmentResult(**p)


def pack(payload):
    encoded = storage.encode(payload)
    return storage.json_bytes(
        {"payload": encoded, "sha256": storage.checksum(storage.json_bytes(encoded))}
    ).decode()


def unpack(raw):
    from island.dynamics._checkpoint_data import strict_load

    envelope = strict_load(raw)
    require(set(envelope) == {"payload", "sha256"}, "Malformed record envelope")
    require(
        envelope["sha256"] == storage.checksum(storage.json_bytes(envelope["payload"])),
        "Record checksum changed",
    )
    return storage.decode(envelope["payload"])


@dataclass(frozen=True)
class ChargeReference:
    """JSON-owned evidence. Checksums prove consistency, not calculation authenticity."""

    json_text: str

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return storage.checksum(
            storage.json_bytes(storage.encode(unpack(self.json_text)))
        )

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p)
                == {
                    "schema",
                    "system",
                    "system_sha256",
                    "correspondence",
                    "charge_result",
                    "import_content",
                    "preparation",
                    "seeds",
                    "unit",
                    "production_validated",
                    "simulation_readiness",
                },
                "Unexpected reference fields",
            )
            require(
                p["schema"] == REFERENCE_SCHEMA and p["unit"] == "elementary_charge",
                "Unsupported schema/unit",
            )
            require(
                p["production_validated"] is False
                and p["simulation_readiness"] == "not_established",
                "Readiness promotion forbidden",
            )
            require(
                p["system_sha256"]
                == storage.checksum(storage.json_bytes(storage.encode(p["system"]))),
                "Source system changed",
            )
            system = system_from(p["system"])
            corr = repeat_correspondence(system)
            require(
                corr["compatible"] and corr == p["correspondence"],
                "Source correspondence invalid or changed",
            )
            require(
                system.number_of_sites <= 100, "AM1-BCC reference exceeds 100 sites"
            )
            seeds = p["seeds"]
            generation = system.metadata["polymer"]["coordinate_generation"]
            require(
                generation["method"] == "local_templates_self_avoiding_random_walk"
                and generation["success"] is True,
                "Require successful local-template construction",
            )
            require(
                set(seeds) == {"template_seed", "assembly_seed"}, "Unexpected seeds"
            )
            for k, v in seeds.items():
                require(
                    type(v) is int and 0 <= v < 2**31 and generation[k] == v,
                    "Construction seed mismatch",
                )
            charge = charge_from(p["charge_result"])
            charge.validate_integrity(system)
            require(
                charge.complete
                and charge.total_formal_charge == 0
                and charge.total_within_tolerance,
                "Invalid charge outcome",
            )
            require(
                all(isfinite(a.charge) for a in charge.assignments.values()),
                "Nonfinite charge",
            )
            content = p["import_content"]
            prep = p["preparation"]
            r = prep["record"]
            require(
                digest(content) == r["imported_result_signature"],
                "Original imported signature changed",
            )
            require(
                content["charge_signature"] == charge.result_signature
                and content["graph"] == graph_signature(system.topology),
                "Imported charge/graph identity mismatch",
            )
            require(
                r["requested_force_field"] == "gaff2"
                and r["charge_method"] == "am1bcc",
                "Require GAFF2 whole-oligomer AM1-BCC",
            )
            require(
                len(content["mapping"]) == system.number_of_sites,
                "Imported mapping count differs from source",
            )
            validate_preparation_record(
                system,
                r,
                prep["record_signature"],
                imported_signature=r["imported_result_signature"],
                imported_provenance=content["provenance"],
                imported_mapping=dict(content["mapping"]),
                source_sha256=content["source_sha256"],
                charge_result=charge,
            )

        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(f"Invalid charge reference: {error}") from error


def create_charge_reference(system, preparation):
    """Bind a real validated preparation; preserve raw charge and original signatures."""
    try:
        require(
            isinstance(system, MolecularSystem)
            and isinstance(preparation, AmberToolsPreparationResult),
            "Expected MolecularSystem and AmberToolsPreparationResult",
        )
        preparation.validate_integrity(system)
        generation = system.metadata["polymer"]["coordinate_generation"]
        source = system_data(system)
        imported = preparation.imported_result
        payload = {
            "schema": REFERENCE_SCHEMA,
            "system": source,
            "system_sha256": storage.checksum(
                storage.json_bytes(storage.encode(source))
            ),
            "correspondence": repeat_correspondence(system),
            "charge_result": record(imported.charge_result),
            "import_content": imported._content(),
            "preparation": {
                "record": dict(preparation.record),
                "record_signature": preparation.record_signature,
            },
            "seeds": {k: generation[k] for k in ("template_seed", "assembly_seed")},
            "unit": "elementary_charge",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        result = ChargeReference(pack(payload))
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot create reference: {error}") from error


def save_record(record, path):
    """Validate first, then exclusive durable publication; no overwrite option."""
    try:
        from .audit import ChargeAudit
        from .conservation_audit import ConservationAudit
        from .observations import RawChargeObservation
        from .pe_template import (
            PEChargePrediction,
            PETemplateModel,
            PETemplateValidation,
        )
        from .projection import ChargeProjection

        require(
            isinstance(
                record,
                (
                    ChargeReference,
                    ChargeAudit,
                    RawChargeObservation,
                    ChargeProjection,
                    ConservationAudit,
                    PETemplateModel,
                    PEChargePrediction,
                    PETemplateValidation,
                ),
            ),
            "Expected charge evidence record",
        )
        record.validate_integrity()
        storage.publish(Path(path), record.json_text.encode())
    except Exception as error:
        raise ChargeReferenceError(f"Cannot publish charge record: {error}") from error


def load_charge_reference(path):
    try:
        result = ChargeReference(Path(path).read_text())
        result.validate_integrity()
        return result
    except Exception as error:
        raise ChargeReferenceError(f"Cannot load charge reference: {error}") from error
