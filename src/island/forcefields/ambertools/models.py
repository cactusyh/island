"""Options and auditable wrapper for one externally prepared Amber import."""

import hashlib
import json
from copy import deepcopy
from dataclasses import dataclass, field
from math import isfinite
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal

from island.core import MolecularSystem
from island.exceptions import AmberToolsInputError, InvalidAmberImportResultError
from island.forcefields.amber import ImportedAmberResult
from island.forcefields.parameterized import ParameterizedSystem


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


@dataclass(frozen=True)
class AmberToolsOptions:
    """Explicit GAFF family, charge method, and bounded external execution."""

    force_field: Literal["gaff", "gaff2"]
    charge_method: Literal["provided", "am1bcc"]
    provided_charges: dict[int, float] | None = None
    timeout_seconds: float = 600.0
    charge_tolerance: float = 1e-4
    work_root: Path | None = None
    amberhome: Path | None = None
    retain_success_artifacts: bool = False

    def __post_init__(self) -> None:
        if self.force_field not in ("gaff", "gaff2"):
            raise AmberToolsInputError("force_field must be 'gaff' or 'gaff2'")
        if self.charge_method not in ("provided", "am1bcc"):
            raise AmberToolsInputError("charge_method must be 'provided' or 'am1bcc'")
        if (self.charge_method == "provided") != (self.provided_charges is not None):
            raise AmberToolsInputError(
                "provided charge mode requires an exact stable-site charge mapping; "
                "AM1-BCC must not receive provided charges"
            )
        if self.provided_charges is not None:
            object.__setattr__(self, "provided_charges", MappingProxyType(
                deepcopy(dict(self.provided_charges))
            ))
        for name in ("timeout_seconds", "charge_tolerance"):
            value = getattr(self, name)
            if (not isinstance(value, (int, float)) or isinstance(value, bool)
                or not isfinite(value) or value <= 0):
                raise AmberToolsInputError(f"{name} must be positive and finite")
        if self.work_root is not None:
            object.__setattr__(self, "work_root", Path(self.work_root))
        if self.amberhome is not None:
            object.__setattr__(self, "amberhome", Path(self.amberhome))

    def __deepcopy__(self, memo: dict[int, object]) -> "AmberToolsOptions":
        copied = type(self)(
            self.force_field, self.charge_method,
            None if self.provided_charges is None else deepcopy(dict(self.provided_charges), memo),
            self.timeout_seconds, self.charge_tolerance,
            self.work_root, self.amberhome, self.retain_success_artifacts,
        )
        memo[id(self)] = copied
        return copied


@dataclass(frozen=True)
class AmberToolsPreparationResult:
    """Signed preparation record around, not instead of, ImportedAmberResult."""

    imported_result: ImportedAmberResult
    record: dict[str, Any]
    record_signature: str
    schema: str = field(default="island_ambertools_preparation_v1", init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "record", MappingProxyType(deepcopy(dict(self.record))))

    def __deepcopy__(self, memo: dict[int, object]) -> "AmberToolsPreparationResult":
        copied = type(self)(
            deepcopy(self.imported_result, memo),
            deepcopy(dict(self.record), memo), self.record_signature,
        )
        memo[id(self)] = copied
        return copied

    def validate_integrity(self, system: MolecularSystem) -> None:
        self.imported_result.validate_integrity(system)
        try:
            consistent = (
                self.record.get("schema") == self.schema
                and self.record.get("imported_result_signature")
                == self.imported_result.result_signature
                and self.record_signature == digest(dict(self.record))
            )
        except (TypeError, ValueError) as error:
            raise InvalidAmberImportResultError(
                f"Malformed AmberTools preparation record: {error}"
            ) from error
        if not consistent:
            raise InvalidAmberImportResultError(
                "AmberTools preparation record or imported-result signature changed"
            )
        coords = {
            str(site_id): system.coordinates.get(site_id).tolist()
            for site_id in sorted(system.topology.sites)
        }
        if self.record.get("input_coordinate_signature") != digest(coords):
            raise InvalidAmberImportResultError(
                "Preparation coordinates changed; imported numerical parameters "
                "remain graph-compatible, but this preparation record is stale"
            )

    def to_parameterized_system(self, system: MolecularSystem) -> ParameterizedSystem:
        """Create a validated owned snapshot with independent preparation provenance."""
        self.validate_integrity(system)
        snapshot = self.imported_result.to_parameterized_system(system)
        snapshot.metadata["ambertools_preparation"] = deepcopy(dict(self.record))
        snapshot.metadata["ambertools_preparation"]["record_signature"] = (
            self.record_signature
        )
        return snapshot
