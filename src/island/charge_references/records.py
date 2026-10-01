"""Owned data-only whole-oligomer AM1-BCC evidence; no charge-transfer engine."""

from dataclasses import dataclass
from math import isfinite
from pathlib import Path

from island.core import MolecularSystem
from island.exceptions import ChargeReferenceError
from island.forcefields.ambertools.models import AmberToolsPreparationResult, digest
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
                r["schema"]
                in {
                    "island_ambertools_preparation_v2",
                    "island_ambertools_preparation_v3",
                },
                "Unsupported preparation schema",
            )
            require(
                r["engine_version"] == ("2" if r["schema"].endswith("v2") else "3"),
                "Unsupported preparation engine",
            )
            require(
                [stage["stage"] for stage in r["stages"]]
                == ["antechamber", "parmchk2", "tleap"]
                and all(
                    type(stage["returncode"]) is int and stage["returncode"] == 0
                    for stage in r["stages"]
                ),
                "Failed or incomplete preparation stages",
            )
            if r["schema"].endswith("v3"):
                from island.forcefields.ambertools.policy import policy_record

                policy = r["size_policy"]
                require(
                    policy
                    == policy_record(
                        policy["max_atoms"], "am1bcc", system.number_of_sites
                    )
                    and r["provided_charge_input"] is None,
                    "Preparation size/charge policy disagrees",
                )
            sqm_hash = r["charge_outcome"]["sqm_out_sha256"]
            require(
                type(sqm_hash) is str
                and len(sqm_hash) == 64
                and set(sqm_hash) <= set("0123456789abcdef"),
                "Missing SQM output checksum",
            )
            require(
                content["provenance"]["force_field"] == "gaff2"
                and content["provenance"]["charge_method"] == "AM1-BCC",
                "Imported method declarations disagree",
            )
            require(
                charge.tolerance == r["charge_validation_tolerance_e"],
                "Charge tolerance disagrees",
            )
            require(
                r["input_coordinate_signature"]
                == digest(r["input_coordinates_angstrom"]),
                "Preparation coordinate digest changed",
            )
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
                digest(r) == prep["record_signature"],
                "Original preparation signature changed",
            )
            require(
                digest(content["provenance"]["ambertools_preparation"])
                == digest(
                    {k: v for k, v in r.items() if k != "imported_result_signature"}
                ),
                "Preparation copies disagree",
            )
            require(
                r["requested_force_field"] == "gaff2"
                and r["charge_method"] == "am1bcc",
                "Require GAFF2 whole-oligomer AM1-BCC",
            )
            require(
                r["charge_outcome"]["mode"] == "am1bcc"
                and r["charge_outcome"]["qm_run"] is True
                and r["charge_outcome"]["convergence_marker"].lower()
                == "calculation completed",
                "Failed/missing QM outcome",
            )
            require(
                r["input_coordinates_angstrom"]
                == {
                    str(s): list(system.coordinates.get(s))
                    for s in system.topology.sites
                },
                "Original input coordinates changed",
            )
            require(
                len(content["mapping"]) == system.number_of_sites
                and set(dict(content["mapping"])) == set(range(system.number_of_sites)),
                "Source index mapping is not a bijection",
            )
            require(
                dict(content["mapping"]) == r["lineage"]["prmtop_index_to_site_id"],
                "Mapping disagrees",
            )
            require(
                set(dict(content["mapping"]).values()) == set(system.topology.sites),
                "Mapping coverage changed",
            )
            require(
                content["source_sha256"] == r["artifact_sha256"]["result.prmtop"],
                "Source checksum mismatch",
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

        require(
            isinstance(record, (ChargeReference, ChargeAudit)),
            "Expected reference or audit record",
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
