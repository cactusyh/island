"""Options and auditable wrapper for one externally prepared Amber import."""

import hashlib
import json
import re
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from math import isfinite
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal

from island.core import MolecularSystem
from island.exceptions import (
    AmberToolsInputError,
    ChargeAssignmentError,
    InvalidAmberImportResultError,
)
from island.forcefields.amber import ImportedAmberResult
from island.forcefields.parameterized import ParameterizedSystem

from .policy import DEFAULT_MAX_ATOMS, policy_record, validate_size_policy

LEGACY_PREPARATION_SCHEMA = "island_ambertools_preparation_v2"
PREPARATION_SCHEMA = "island_ambertools_preparation_v3"

SERIALIZATION_TOLERANCE = 1e-5  # elementary charges, per site
SQM_SUCCESS = re.compile(r"calculation\s+completed", re.IGNORECASE)
HEX = set("0123456789abcdef")


def _invalid(message: str) -> None:
    raise InvalidAmberImportResultError(f"AmberTools preparation integrity: {message}")


def _sha256(value: object, label: str) -> None:
    if (not isinstance(value, str) or len(value) != 64
        or any(character not in HEX for character in value)):
        _invalid(f"{label} must be a lowercase SHA-256 digest")


def _mapping(value: object, label: str) -> Mapping:
    if not isinstance(value, Mapping):
        _invalid(f"{label} must be a mapping")
    return value


def _site_bijection(
    value: object, site_ids: set[int], label: str, *, zero_based: bool | None = None,
) -> None:
    records = _mapping(value, label)
    if (any(not isinstance(v, int) or isinstance(v, bool)
            for v in records.values())
        or set(records.values()) != site_ids or len(records) != len(site_ids)):
        _invalid(f"{label} must cover each stable site ID exactly once")
    if any(type(index) is not int for index in records):
        _invalid(f"{label} must use integer source indices")
    if zero_based is True and set(records) != set(range(len(site_ids))):
        _invalid(f"{label} must use contiguous zero-based source indices")
    if zero_based is False and set(records) != set(range(1, len(site_ids) + 1)):
        _invalid(f"{label} must use contiguous one-based MOL2 indices")


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
    max_atoms: int = DEFAULT_MAX_ATOMS
    charge_source: str | None = None

    def __post_init__(self) -> None:
        validate_size_policy(self.max_atoms, self.charge_method)
        if self.charge_source is not None and (self.charge_method != "provided" or type(self.charge_source) is not str or not self.charge_source.strip()):
            raise AmberToolsInputError("charge_source must be nonempty text for provided charges only")
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
            self.max_atoms, self.charge_source,
        )
        memo[id(self)] = copied
        return copied


@dataclass(frozen=True)
class AmberToolsPreparationResult:
    """Signed preparation record around, not instead of, ImportedAmberResult."""

    imported_result: ImportedAmberResult
    record: dict[str, Any]
    record_signature: str
    schema: str = field(default=PREPARATION_SCHEMA, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.record, Mapping):
            _invalid("record must be a mapping")
        object.__setattr__(self, "record", MappingProxyType(deepcopy(dict(self.record))))
        object.__setattr__(self, "schema", self.record.get("schema"))

    def __deepcopy__(self, memo: dict[int, object]) -> "AmberToolsPreparationResult":
        copied = type(self)(
            deepcopy(self.imported_result, memo),
            deepcopy(dict(self.record), memo), self.record_signature,
        )
        memo[id(self)] = copied
        return copied

    def validate_integrity(self, system: MolecularSystem) -> None:
        """Reject malformed reconstructed data through the import error boundary."""
        try:
            self._validate_integrity(system)
        except InvalidAmberImportResultError:
            raise
        except (AmberToolsInputError, ChargeAssignmentError, TypeError, ValueError, KeyError, AttributeError, OverflowError) as error:
            raise InvalidAmberImportResultError(
                f"AmberTools preparation integrity: malformed record: {error}"
            ) from error

    def _validate_integrity(self, system: MolecularSystem) -> None:
        if not isinstance(self.imported_result, ImportedAmberResult):
            _invalid("imported_result must be an ImportedAmberResult")
        self.imported_result.validate_integrity(system)
        record = dict(self.record)
        if record.get("schema") not in (PREPARATION_SCHEMA, LEGACY_PREPARATION_SCHEMA):
            _invalid("unsupported or missing preparation schema")
        if record.get("engine_version") != ("3" if record["schema"] == PREPARATION_SCHEMA else "2"):
            _invalid("unsupported or missing AmberTools engine version")
        _sha256(self.record_signature, "record_signature")
        if record.get("imported_result_signature") != self.imported_result.result_signature:
            _invalid("outer imported-result signature does not match imported_result")
        inner = self.imported_result.provenance.get("ambertools_preparation")
        if not isinstance(inner, Mapping):
            _invalid("signed imported provenance lacks an AmberTools preparation payload")
        shared = {key: value for key, value in record.items()
                  if key != "imported_result_signature"}
        if digest(shared) != digest(dict(inner)):
            _invalid("outer preparation record contradicts signed imported provenance")
        try:
            consistent = self.record_signature == digest(record)
        except (TypeError, ValueError) as error:
            _invalid(f"malformed record content: {error}")
        if not consistent:
            _invalid("record content signature changed")

        requested = record.get("requested_force_field")
        method = record.get("charge_method")
        if (not isinstance(requested, str) or requested not in {"gaff", "gaff2"}
            or not isinstance(method, str) or method not in {"provided", "am1bcc"}):
            _invalid("force-field or charge-method declaration is missing/unsupported")
        if (self.imported_result.provenance.get("force_field") != requested
            or self.imported_result.provenance.get("charge_method")
            != ("provided" if method == "provided" else "AM1-BCC")):
            _invalid("force-field/charge-method declarations contradict imported source")
        if record["schema"] == LEGACY_PREPARATION_SCHEMA:
            validate_size_policy(DEFAULT_MAX_ATOMS, method, len(system.topology.sites), error_type=InvalidAmberImportResultError)
            if "size_policy" in record or "provided_charge_input" in record:
                _invalid("v2 records cannot declare a v3 execution policy")
        else:
            policy = _mapping(record.get("size_policy"), "size_policy")
            if type(policy.get("actual_atoms")) is not int or policy != policy_record(policy.get("max_atoms"), method, len(system.topology.sites)):
                _invalid("size policy contradicts count or charge method")
            supplied = record.get("provided_charge_input")
            if method == "am1bcc":
                if supplied is not None:
                    _invalid("AM1-BCC cannot claim supplied charge provenance")
            else:
                supplied = _mapping(supplied, "provided_charge_input")
                if set(supplied) != {"charges", "sha256", "source", "method_provenance", "charges_file_sha256"}:
                    _invalid("invalid supplied charge record structure")
                charges = _mapping(supplied["charges"], "supplied charges")
                ids = set(system.topology.sites)
                if set(charges) != {str(s) for s in ids}:
                    _invalid("supplied charge coverage differs")
                if any(type(q) not in (int, float) or not isfinite(q) for q in charges.values()):
                    _invalid("supplied charges must be finite")
                if digest(dict(charges)) != supplied["sha256"]:
                    _invalid("supplied charge checksum mismatch")
                from island.forcefields.charges.engines import ProvidedChargeEngine
                ProvidedChargeEngine().assign(system, {int(s): q for s, q in charges.items()},
                                             source="Recorded supplied charges", tolerance=record["charge_validation_tolerance_e"])
                for site in ids:
                    if abs(charges[str(site)] - self.imported_result.charge_result.assignments[site].charge) > SERIALIZATION_TOLERANCE:
                        _invalid("supplied charge differs from imported charge")
                if not isinstance(supplied["source"], str) or not supplied["source"].strip() or supplied["method_provenance"] != "unverified_user_supplied":
                    _invalid("supplied charge source must explicitly retain unverified method provenance")
                serialized = "\n".join(f"{charges[str(site)]:.10f}" for site in sorted(ids)) + "\n"
                expected_sha = hashlib.sha256(serialized.encode()).hexdigest()
                if supplied["charges_file_sha256"] != expected_sha or record["artifact_sha256"].get("charges.txt") != expected_sha:
                    _invalid("supplied charge file checksum mismatch")
        outcome = _mapping(record.get("charge_outcome"), "charge_outcome")
        if outcome.get("mode") != method or outcome.get("qm_run") is not (method == "am1bcc"):
            _invalid("charge outcome contradicts selected charge method")
        if method == "am1bcc":
            _sha256(outcome.get("sqm_out_sha256"), "charge_outcome.sqm_out_sha256")
            if (not isinstance(outcome.get("convergence_marker"), str)
                or SQM_SUCCESS.fullmatch(outcome["convergence_marker"]) is None):
                _invalid("AM1-BCC completion marker is absent or invalid")
            if "serialization_tolerance_e" in outcome:
                _invalid("AM1-BCC outcome must not claim provided-charge serialization")
        elif any(key in outcome for key in (
            "sqm_out_sha256", "convergence_marker", "sqm_warning_lines",
        )):
            _invalid("provided-charge record must not claim SQM execution")
        for label in ("charge_validation_tolerance_e", "provided_charge_tolerance_e"):
            value = record.get(label)
            if (not isinstance(value, (int, float)) or isinstance(value, bool)
                or not isfinite(value) or value <= 0):
                _invalid(f"{label} is missing or invalid")
        if (record["provided_charge_tolerance_e"] != SERIALIZATION_TOLERANCE
            or (method == "provided" and (
                isinstance(outcome.get("serialization_tolerance_e"), bool)
                or outcome.get("serialization_tolerance_e") != SERIALIZATION_TOLERANCE
            ))):
            _invalid("provided-charge serialization tolerance is missing or inconsistent")
        if record["charge_validation_tolerance_e"] != self.imported_result.charge_result.tolerance:
            _invalid("charge-validation tolerance contradicts imported charge result")

        coordinates = _mapping(record.get("input_coordinates_angstrom"),
                               "input_coordinates_angstrom")
        ids = set(system.topology.sites)
        if set(coordinates) != {str(site_id) for site_id in ids}:
            _invalid("input coordinate coverage differs from authoritative sites")
        for site_id, xyz in coordinates.items():
            if (not isinstance(xyz, (list, tuple)) or len(xyz) != 3
                or any(not isinstance(value, (int, float)) or isinstance(value, bool)
                       or not isfinite(value) for value in xyz)):
                _invalid(f"input coordinate {site_id} is not a finite 3-vector")
        if record.get("input_coordinate_signature") != digest(dict(coordinates)):
            _invalid("recorded coordinate content and digest disagree")
        coords = {
            str(site_id): system.coordinates.get(site_id).tolist()
            for site_id in sorted(system.topology.sites)
        }
        if record["input_coordinate_signature"] != digest(coords):
            _invalid("preparation coordinates changed; record is stale")
        prepared_coordinates = _mapping(
            record.get("antechamber_coordinates_angstrom"),
            "antechamber_coordinates_angstrom",
        )
        if set(prepared_coordinates) != {str(site_id) for site_id in ids}:
            _invalid("antechamber coordinate coverage differs from authoritative sites")
        for site_id, xyz in prepared_coordinates.items():
            if (not isinstance(xyz, (list, tuple)) or len(xyz) != 3
                or any(not isinstance(value, (int, float)) or isinstance(value, bool)
                       or not isfinite(value) for value in xyz)):
                _invalid(f"antechamber coordinate {site_id} is not a finite 3-vector")

        assigned_cip = _mapping(record.get("expected_cip_by_site"),
                                "expected_cip_by_site")
        if any(not isinstance(site_id, int) or site_id not in ids
               or not isinstance(label, str) or label not in {"R", "S"}
               or (system.topology.sites[site_id].metadata.get("cip_label")
                   not in (None, label))
               for site_id, label in assigned_cip.items()):
            _invalid("expected CIP assignments contradict authoritative sites")
        for site_id, site in system.topology.sites.items():
            label = site.metadata.get("cip_label")
            if label in {"R", "S"} and assigned_cip.get(site_id) != label:
                _invalid(f"expected CIP assignment missing for site {site_id}")

        lineage = _mapping(record.get("lineage"), "lineage")
        from island.forcefields.ambertools.lineage import generated_atom_name

        expected_names = {
            generated_atom_name(system.topology.sites[site_id].element, index): site_id
            for index, site_id in enumerate(sorted(ids))
        }
        name_map = lineage.get("input_name_to_site_id")
        if name_map != expected_names:
            _invalid("input generated atom names or stable-site mapping changed")
        _site_bijection(lineage.get("typed_mol2_index_to_site_id"), ids,
                        "typed_mol2_index_to_site_id", zero_based=False)
        _site_bijection(lineage.get("prmtop_index_to_site_id"), ids,
                        "prmtop_index_to_site_id", zero_based=True)
        if dict(lineage["prmtop_index_to_site_id"]) != dict(self.imported_result.mapping):
            _invalid("prmtop lineage contradicts signed imported atom mapping")

        artifacts = _mapping(record.get("artifact_sha256"), "artifact_sha256")
        for name in ("typed.mol2", "typed.frcmod", "result.prmtop",
                     "result.rst7", "leap.in", "leap.log"):
            _sha256(artifacts.get(name), f"artifact_sha256[{name}]")
        if artifacts["result.prmtop"] != self.imported_result.source_sha256:
            _invalid("prmtop artifact checksum contradicts imported source checksum")
        for name in ("input_mol2_sha256", "input_lineage_sha256"):
            _sha256(record.get(name), name)
        for name in ("force_field_data", "leaprc"):
            entry = _mapping(record.get(name), name)
            _sha256(entry.get("sha256"), f"{name}.sha256")
            if not isinstance(entry.get("path"), str) or not entry["path"]:
                _invalid(f"{name}.path is missing")
        if not isinstance(record["force_field_data"].get("header"), str) or not record["force_field_data"]["header"]:
            _invalid("force-field data header is missing")
        home = record.get("amberhome")
        if not isinstance(home, str) or not home.startswith("/"):
            _invalid("selected AMBERHOME is missing or not absolute")
        if (not record["force_field_data"]["path"].endswith(f"/{requested}.dat")
            or not record["leaprc"]["path"].endswith(f"/leaprc.{requested}")):
            _invalid("force-field data/leaprc paths contradict requested family")
        if (record["force_field_data"]["path"] != f"{home}/dat/leap/parm/{requested}.dat"
            or record["leaprc"]["path"] != f"{home}/dat/leap/cmd/leaprc.{requested}"):
            _invalid("force-field data/leaprc are outside selected AMBERHOME")
        executables = _mapping(record.get("executables"), "executables")
        executable_hashes = _mapping(record.get("executable_sha256"),
                                     "executable_sha256")
        versions = _mapping(record.get("tool_versions"), "tool_versions")
        for name in ("antechamber", "parmchk2", "tleap"):
            if (not isinstance(executables.get(name), str) or not executables[name]
                or not isinstance(versions.get(name), str) or not versions[name]):
                _invalid(f"{name} executable/version provenance is incomplete")
            _sha256(executable_hashes.get(name), f"executable_sha256[{name}]")
            if executables[name] != f"{home}/bin/{name}":
                _invalid(f"{name} executable contradicts selected AMBERHOME")
        stages = record.get("stages")
        if (not isinstance(stages, (list, tuple)) or len(stages) != 3
            or [stage.get("stage") if isinstance(stage, Mapping) else None
                for stage in stages] != ["antechamber", "parmchk2", "tleap"]):
            _invalid("required antechamber/parmchk2/tleap stage records are incomplete")
        required_options = (
            {"-i": "input.mol2", "-fi": "mol2", "-o": "typed.mol2", "-fo": "mol2",
             "-at": requested, "-c": "rc" if method == "provided" else "bcc",
             "-nc": str(sum(site.formal_charge for site in system.topology.sites.values())),
             "-m": "1", "-s": "2", "-j": "4", "-du": "yes", "-pf": "no",
             **({"-cf": "charges.txt"} if method == "provided" else {})},
            {"-i": "typed.mol2", "-f": "mol2", "-o": "typed.frcmod",
             "-s": "1" if requested == "gaff" else "2"},
            {"-f": "leap.in"},
        )
        for stage, expected_options in zip(stages, required_options, strict=True):
            command = stage.get("command")
            if not isinstance(command, (list, tuple)) or not command:
                _invalid(f"{stage['stage']} command is missing or empty")
            if any(not isinstance(item, str) for item in command):
                _invalid(f"{stage['stage']} command contains nontext arguments")
            if (type(stage.get("returncode")) is not int or stage["returncode"] != 0
                or command[0] != executables[stage["stage"]]):
                _invalid(f"{stage['stage']} did not record successful selected executable")
            # These three tools use option/value pairs. Parse once so a second,
            # contradictory flag cannot hide behind one matching adjacent pair.
            if len(command) % 2 != 1:
                _invalid(f"{stage['stage']} command has an option without a value")
            options = {}
            for flag, value in zip(command[1::2], command[2::2], strict=True):
                if not flag.startswith("-") or not value:
                    _invalid(f"{stage['stage']} command has a malformed option/value")
                if flag in options:
                    _invalid(f"{stage['stage']} command has duplicate option {flag}")
                options[flag] = value
            if any(options.get(flag) != value for flag, value in expected_options.items()):
                _invalid(f"{stage['stage']} command contradicts requested settings")
            if stage["stage"] == "antechamber" and method == "am1bcc" and "-cf" in options:
                _invalid("AM1-BCC command must not supply a provided-charge file")

    def to_parameterized_system(self, system: MolecularSystem) -> ParameterizedSystem:
        """Create a validated owned snapshot with independent preparation provenance."""
        self.validate_integrity(system)
        snapshot = self.imported_result.to_parameterized_system(system)
        snapshot.metadata["ambertools_preparation"] = deepcopy(dict(self.record))
        snapshot.metadata["ambertools_preparation"]["record_signature"] = (
            self.record_signature
        )
        return snapshot
