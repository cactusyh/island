"""Explicit external-PDMS preparation boundary, without implicit source discovery.

C5's data schema remains authoritative. No complete authorized bundle has been
supplied to qualify a numerical adapter. These results are blocked preparation
records, never PreparedForceField objects or executable parameter assignments.
"""

import re
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from island.exceptions import ForceFieldError, IslandError

FAMILY = "external_pdms_class2"
ADAPTER = "island_external_pdms_preparation_v1"
SCHEMA = "island_external_pdms_preparation_result_v1"


@dataclass(frozen=True)
class ExternalPDMSOptions:
    """Caller-selected path and independent pin; neither is inferred or searched.

    None/None explicitly requests the unavailable-asset diagnostic. A path
    without a valid independent pin is refused before opening the candidate.
    Source/evidence assets are resolved relative to the supplied bundle path.
    """

    bundle_path: Path | None = None
    expected_sha256: str | None = None

    def __post_init__(self):
        if self.bundle_path is not None:
            if (
                not isinstance(self.bundle_path, (str, Path))
                or not str(self.bundle_path).strip()
            ):
                raise ValueError("explicit_bundle_path_required")
            object.__setattr__(self, "bundle_path", Path(self.bundle_path))
        if self.expected_sha256 is not None and type(self.expected_sha256) is not str:
            raise ValueError("independent_bundle_sha256_must_be_string")


class ExternalPDMSBlockedError(ForceFieldError):
    """A classified preparation result cannot authorize a numerical backend."""

    def __init__(self, result, *, validation_error=None):
        self.result = result
        # Retain the original result even if its envelope was deliberately
        # tampered with after construction. Corruption must not reach an engine.
        try:
            payload = result.payload
        except (ValueError, TypeError, KeyError, AttributeError, IslandError):
            payload = {}
        if type(payload) is not dict:
            payload = {}
        self.status = payload.get("status", "parameter_incomplete")
        if self.status not in (
            "asset_unavailable",
            "parameter_incomplete",
            "typing_blocked",
        ):
            self.status = "parameter_incomplete"
        self.diagnostics = deepcopy(payload.get("diagnostics", []))
        self.validation_diagnostic = (
            None
            if validation_error is None
            else {
                "stage": "blocked_record_invalid",
                "exception_type": type(validation_error).__name__,
                "context": str(validation_error),
            }
        )
        super().__init__(
            f"{self.status}: external PDMS execution is blocked; "
            "inspect result.payload for interactions, dependencies and diagnostics"
        )


@dataclass(frozen=True)
class ExternalPDMSPreparation:
    """JSON-owned, checksummed blocked report; payload access returns a copy."""

    json_text: str

    def __post_init__(self):
        self.validate_integrity()

    @property
    def payload(self):
        from island.charge_references.records import unpack

        return unpack(self.json_text)

    @property
    def identity(self):
        from hashlib import sha256

        self.validate_integrity()
        return sha256(self.json_text.encode()).hexdigest()

    def validate_integrity(self, system=None):
        from .external_pdms import LIMITS, PROFILE, require, requirements

        p = self.payload
        require(p["schema"] == SCHEMA and p["adapter"] == ADAPTER, "adapter_changed")
        require(p["family"] == FAMILY, "external_family_changed")
        require(p["requested_profile"] == PROFILE, "requested_profile_changed")
        require(
            p["status"]
            in ("asset_unavailable", "parameter_incomplete", "typing_blocked")
            and bool(p["diagnostics"]),
            "unclassified_or_promoted_preparation",
        )
        require(
            p["execution_authorized"] is False
            and p["execution_profile"] is None
            and p["parameter_assignment"] is None
            and p["native_pcff_acceptance_gate"] == 1
            and p["acceptance_gate"] == 1
            and p["backend_executions"] == []
            and p["limitations"] == LIMITS,
            "blocked_record_cannot_grant_execution_authority",
        )
        diagnostics = p["diagnostics"]
        require(
            type(diagnostics) is list
            and all(
                type(d) is dict and type(d.get("stage")) is str for d in diagnostics
            ),
            "malformed_blocked_diagnostics",
        )
        typing_stage = any(d["stage"] == "final_graph_invalid" for d in diagnostics)
        if p["status"] == "typing_blocked":
            require(
                p["requirements"] is None
                and len(diagnostics) == 1
                and typing_stage
                and set(diagnostics[0]) == {"stage", "context"}
                and type(diagnostics[0]["context"]) is str
                and bool(diagnostics[0]["context"])
                and p["missing_interactions"] == {}
                and p["unresolved_dependencies"] == []
                and p["bundle_identity"] is None
                and p["import_validation"] is None,
                "typing_blocked_record_changed",
            )
            # No valid graph inventory exists. Repeating the rejected graph
            # validation would replace its recorded outcome with a raw error.
            return
        require(
            type(p["requirements"]) is dict and not typing_stage,
            "blocked_status_requirements_mismatch",
        )
        if system is not None:
            require(p["requirements"] == requirements(system), "stale_final_graph")

    def require_executable(self, system=None):
        try:
            self.validate_integrity(system)
        except (
            ValueError,
            TypeError,
            KeyError,
            AttributeError,
            IndexError,
            IslandError,
        ) as error:
            raise ExternalPDMSBlockedError(self, validation_error=error) from error
        raise ExternalPDMSBlockedError(self)


class ExternalPDMSAdapter:
    """Named public adapter. Only the unavailable/incomplete boundary is qualified.

    A complete schema import alone cannot enable export, Context creation or MD.
    A future numerical adapter needs an explicitly supplied authorized asset and
    executed independent source/energy/force/continuation qualification first.
    """

    def prepare(self, system, options):
        from island.charge_references.records import pack

        from .external_pdms import (
            LIMITS,
            PROFILE,
            load_external_pdms_bundle,
            requirements,
        )

        if type(options) is not ExternalPDMSOptions:
            raise TypeError("ExternalPDMSOptions required")
        # Revalidate even a deserialized or deliberately tampered options object.
        options = ExternalPDMSOptions(options.bundle_path, options.expected_sha256)
        diagnostics, inventory, candidate, imported = [], None, None, None
        status = "asset_unavailable"
        try:
            inventory = requirements(system)
        except (
            ValueError,
            TypeError,
            KeyError,
            AttributeError,
            IndexError,
            IslandError,
        ) as error:
            status = "typing_blocked"
            diagnostics.append({"stage": "final_graph_invalid", "context": str(error)})
        if inventory is not None:
            if options.bundle_path is None:
                diagnostics.append({"stage": "explicit_bundle_path_required"})
            elif not re.fullmatch(r"[0-9a-f]{64}", options.expected_sha256 or ""):
                status = "parameter_incomplete"
                diagnostics.append({"stage": "independent_bundle_sha256_required"})
            else:
                try:
                    bundle, report = load_external_pdms_bundle(
                        options.bundle_path,
                        system,
                        expected_sha256=options.expected_sha256,
                    )
                    imported = report.payload
                    status = imported["status"]
                    diagnostics.extend(imported["diagnostics"])
                    if bundle is not None:
                        candidate = bundle.payload
                    if report.complete:
                        # Data integrity is not an executed backend qualification.
                        status = "asset_unavailable"
                        diagnostics.append(
                            {"stage": "backend_qualification_unavailable"}
                        )
                except (
                    ValueError,
                    TypeError,
                    KeyError,
                    OSError,
                    UnicodeError,
                    IslandError,
                ) as error:
                    status = "parameter_incomplete"
                    diagnostics.append(
                        {"stage": "malformed_bundle", "context": str(error)}
                    )
        coverage = {} if inventory is None else inventory["interactions"]
        supplied = candidate.get("interactions", {}) if type(candidate) is dict else {}
        if type(supplied) is not dict:
            supplied = {}
        missing = {k: v for k, v in coverage.items() if supplied.get(k) is None}
        dependencies = [
            {"consumer": k, "dependency": d}
            for k, v in coverage.items()
            for d in v["dependencies"]
            if k in missing or d in missing
        ]
        return ExternalPDMSPreparation(
            pack(
                {
                    "schema": SCHEMA,
                    "adapter": ADAPTER,
                    "family": FAMILY,
                    "requested_profile": PROFILE,
                    "status": status,
                    "requirements": inventory,
                    "bundle_identity": None
                    if imported is None
                    else imported["bundle_identity"],
                    "import_validation": imported,
                    "missing_interactions": missing,
                    "unresolved_dependencies": dependencies,
                    "diagnostics": diagnostics,
                    "execution_authorized": False,
                    "execution_profile": None,
                    "parameter_assignment": None,
                    "native_pcff_acceptance_gate": 1,
                    "acceptance_gate": 1,
                    "backend_executions": [],
                    "limitations": LIMITS,
                }
            )
        )
