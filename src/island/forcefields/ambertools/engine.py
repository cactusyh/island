"""Bounded, auditable AmberTools execution feeding the existing Amber importer."""

import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np

from island.core import MolecularSystem
from island.exceptions import (
    AmberImportError,
    AmberToolsInputError,
    AmberToolsStageError,
    AmberToolsUnavailableError,
    IncompleteChargeAssignmentError,
    InvalidChargeDefinitionError,
)
from island.forcefields.amber import import_amber_prmtop
from island.forcefields.ambertools.lineage import (
    _coordinate_cip,
    parse_and_validate_mol2,
    prepare_input,
)
from island.forcefields.ambertools.models import (
    PREPARATION_SCHEMA,
    SERIALIZATION_TOLERANCE,
    SQM_SUCCESS,
    AmberToolsOptions,
    AmberToolsPreparationResult,
    digest,
)
from island.forcefields.charges import ProvidedChargeEngine

UNRESOLVED = re.compile(r"ATTN\s*,?\s*need\s+revision|missing\s+parameter|"
                        r"parameter\s+not\s+found", re.IGNORECASE)
FATAL = re.compile(r"\b(?:fatal|error:|failed|failure)\b", re.IGNORECASE)
NONZERO_ERRORS = re.compile(r"\bErrors\s*[=:]\s*[1-9][0-9]*\b", re.IGNORECASE)
SQM_FAILURE = re.compile(
    r"not\s+converged|failed\s+to\s+converge|scf\s+convergence\s+failure",
    re.IGNORECASE,
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _header(path: Path) -> str:
    return next(
        (line.strip()[:300] for line in path.read_text(errors="replace").splitlines()
         if line.strip()),
        "unavailable (no nonblank header)",
    )


@dataclass(frozen=True)
class _Toolchain:
    executables: dict[str, str]
    amberhome: Path
    leaprc: Path
    data_file: Path
    versions: dict[str, str]
    package: dict[str, str] = field(default_factory=dict)


def _probe_version(executable: str) -> str:
    try:
        result = subprocess.run(
            [executable, "-h"], capture_output=True, text=True,
            timeout=5, check=False, shell=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return f"unavailable ({type(error).__name__})"
    text = (result.stdout + "\n" + result.stderr).strip()
    for line in text.splitlines():
        cleaned = re.sub(r"\x1b\[[0-9;]*m", "", line)
        if re.search(
            r"(?:antechamber|parmchk2|tleap|AmberTools)\s+"
            r"(?:version\s*)?v?\d+(?:\.\d+)+", cleaned, re.IGNORECASE,
        ):
            return cleaned[:300]
    return "unavailable (help output has no executable version banner)"


def _discover(options: AmberToolsOptions) -> _Toolchain:
    found = {name: shutil.which(name) for name in ("antechamber", "parmchk2", "tleap")}
    missing = [name for name, path in found.items() if path is None]
    if missing:
        raise AmberToolsUnavailableError(
            "AmberTools executable(s) missing: " + ", ".join(missing)
            + "; install AmberTools and set AMBERHOME/PATH"
        )
    home = options.amberhome or (Path(os.environ["AMBERHOME"])
                                 if os.environ.get("AMBERHOME") else
                                 Path(found["tleap"]).resolve().parent.parent)
    home = home.resolve()
    bin_directory = home / "bin"
    resolved_executables = {
        name: str(Path(path).resolve()) for name, path in found.items()
    }
    inconsistent = [name for name, path in resolved_executables.items()
                    if Path(path).parent != bin_directory]
    if inconsistent:
        raise AmberToolsUnavailableError(
            "AmberTools executables do not all belong to selected AMBERHOME/bin "
            f"({bin_directory}): {inconsistent}; mixed installations are unsupported"
        )
    leaprc = home / "dat" / "leap" / "cmd" / f"leaprc.{options.force_field}"
    data_file = home / "dat" / "leap" / "parm" / f"{options.force_field}.dat"
    for path in (leaprc, data_file):
        if not path.is_file() or path.stat().st_size == 0:
            raise AmberToolsUnavailableError(
                f"Selected {options.force_field} data file missing/empty: {path}"
            )
    if data_file.name not in leaprc.read_text(errors="replace"):
        raise AmberToolsUnavailableError(
            f"{leaprc} does not explicitly load {data_file.name}; selected "
            "force-field data identity cannot be verified"
        )
    package: dict[str, str] = {}
    package_records = sorted((home / "conda-meta").glob("ambertools-*.json"))
    if len(package_records) == 1:
        try:
            payload = json.loads(package_records[0].read_text())
            if payload.get("name") == "ambertools":
                package = {
                    key: str(payload[key]) for key in ("version", "build", "channel")
                    if key in payload
                }
                package["metadata_path"] = str(package_records[0])
                package["metadata_sha256"] = _sha(package_records[0])
        except (OSError, ValueError, TypeError):
            package = {}
    return _Toolchain(
        resolved_executables, home, leaprc.resolve(), data_file.resolve(),
        {name: _probe_version(path) for name, path in found.items()},
        package,
    )


def _run_stage(
    stage: str, command: list[str], directory: Path, timeout: float,
    *, environment: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Execute one process group without a shell; kill descendants on timeout."""
    started = time.monotonic()
    try:
        process = subprocess.Popen(
            command, cwd=directory, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, shell=False,
            start_new_session=True,
            env=environment,
        )
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = process.communicate()
            (directory / f"{stage}.stdout.log").write_text(stdout)
            (directory / f"{stage}.stderr.log").write_text(stderr)
            raise AmberToolsStageError(
                stage, f"timed out after {timeout} seconds; process group killed",
                artifact_dir=str(directory), command=tuple(command),
                stdout=stdout, stderr=stderr,
            ) from None
    except OSError as error:
        raise AmberToolsStageError(
            stage, f"could not start executable: {error}",
            artifact_dir=str(directory), command=tuple(command),
        ) from error
    (directory / f"{stage}.stdout.log").write_text(stdout)
    (directory / f"{stage}.stderr.log").write_text(stderr)
    if process.returncode != 0:
        raise AmberToolsStageError(
            stage, f"exit code {process.returncode}: {stderr[-1200:] or stdout[-1200:]}",
            artifact_dir=str(directory), command=tuple(command),
            stdout=stdout, stderr=stderr, returncode=process.returncode,
        )
    if (UNRESOLVED.search(stdout + "\n" + stderr)
        or NONZERO_ERRORS.search(stdout + "\n" + stderr)
        or any(
        FATAL.search(line) for line in (stdout + "\n" + stderr).splitlines()
        if "0 errors" not in line.lower()
    )):
        raise AmberToolsStageError(
            stage, "tool reported unresolved parameters or a fatal diagnostic",
            artifact_dir=str(directory), command=tuple(command),
            stdout=stdout, stderr=stderr, returncode=process.returncode,
        )
    return {
        "stage": stage, "command": command, "returncode": process.returncode,
        "elapsed_seconds": round(time.monotonic() - started, 6),
        "stdout_sha256": _sha(directory / f"{stage}.stdout.log"),
        "stderr_sha256": _sha(directory / f"{stage}.stderr.log"),
        "stdout_tail": stdout[-4000:], "stderr_tail": stderr[-4000:],
    }


def _require_artifact(stage: str, path: Path, directory: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise AmberToolsStageError(
            stage, f"required artifact missing or empty: {path.name}",
            artifact_dir=str(directory),
        )
    return _sha(path)


def _check_charge_values(
    expected: dict[int, float], actual: dict[int, float], stage: str,
    directory: Path, tolerance: float = SERIALIZATION_TOLERANCE,
) -> None:
    if set(expected) != set(actual):
        raise AmberToolsStageError(
            stage, "charge site coverage changed", artifact_dir=str(directory)
        )
    mismatches = [
        (site_id, expected[site_id], actual[site_id])
        for site_id in sorted(expected)
        if abs(expected[site_id] - actual[site_id]) > tolerance
    ]
    if mismatches:
        raise AmberToolsStageError(
            stage,
            f"per-site charges changed beyond {tolerance} e: {mismatches[:8]}",
            artifact_dir=str(directory),
        )


def _prmtop_lineage(
    prmtop: Path, restart: Path, system: MolecularSystem,
    names: dict[str, int], expected_cip: dict[int, str], directory: Path,
) -> tuple[dict[int, int], dict[int, float]]:
    try:
        import parmed as pmd
    except ImportError as error:
        raise AmberToolsUnavailableError("ParmEd is required; install island[amber]") from error
    try:
        parm = pmd.load_file(str(prmtop), str(restart))
    except Exception as error:
        raise AmberToolsStageError(
            "tleap", f"cannot parse prmtop/restart: {error}",
            artifact_dir=str(directory),
        ) from error
    if len(parm.atoms) != len(names) or parm.coordinates is None:
        raise AmberToolsStageError(
            "tleap", "prmtop/restart atom inventory or coordinates missing",
            artifact_dir=str(directory),
        )
    try:
        coordinates = np.asarray(parm.coordinates, dtype=float)
        finite = bool(np.all(np.isfinite(coordinates)))
    except (TypeError, ValueError) as error:
        raise AmberToolsStageError(
            "tleap", f"restart coordinates are not numeric: {error}",
            artifact_dir=str(directory),
        ) from error
    if coordinates.shape != (len(parm.atoms), 3) or not finite:
        bad_rows = (
            np.flatnonzero(~np.all(np.isfinite(coordinates), axis=1)).tolist()
            if coordinates.shape == (len(parm.atoms), 3) else []
        )
        raise AmberToolsStageError(
            "tleap", "restart coordinates must be a finite N x 3 array; "
            f"shape={coordinates.shape}, nonfinite_source_indices={bad_rows[:8]}",
            artifact_dir=str(directory),
        )
    mapping: dict[int, int] = {}
    charges: dict[int, float] = {}
    positions = {}
    for atom in parm.atoms:
        name = atom.name
        if name not in names or name in {parm.atoms[i].name for i in mapping}:
            raise AmberToolsStageError(
                "tleap", f"lost/duplicate/truncated generated atom name {name}",
                artifact_dir=str(directory),
            )
        site_id = names[name]
        site = system.topology.sites[site_id]
        if atom.atomic_number != site.atomic_number:
            raise AmberToolsStageError(
                "tleap", f"element changed for atom {name}/site {site_id}",
                artifact_dir=str(directory),
            )
        mapping[atom.idx] = site_id
        charge = float(atom.charge)
        if not np.isfinite(charge):
            raise AmberToolsStageError(
                "tleap", f"nonfinite charge at source atom {atom.idx}",
                artifact_dir=str(directory),
            )
        charges[site_id] = charge
        positions[site_id] = tuple(float(x) for x in coordinates[atom.idx])
    if set(mapping.values()) != set(system.topology.sites):
        raise AmberToolsStageError(
            "tleap", "final atom-name lineage is not bijective",
            artifact_dir=str(directory),
        )
    try:
        _coordinate_cip(system, positions, expected_cip, stage="tleap restart")
    except AmberToolsInputError as error:
        raise AmberToolsStageError(
            "tleap", str(error), artifact_dir=str(directory)
        ) from error
    return mapping, charges


class AmberToolsParameterizationEngine:
    """Run GAFF/GAFF2 preparation, then delegate parameter parsing to Phase 4D1."""

    engine_version = "2"

    def parameterize(
        self, system: MolecularSystem, options: AmberToolsOptions,
    ) -> AmberToolsPreparationResult:
        if not isinstance(options, AmberToolsOptions):
            raise TypeError("options must be AmberToolsOptions")
        if options.charge_method == "provided":
            try:
                provided_result = ProvidedChargeEngine().assign(
                    system, dict(options.provided_charges),
                    source="AmberTools provided-charge input",
                    tolerance=options.charge_tolerance,
                )
            except (InvalidChargeDefinitionError, IncompleteChargeAssignmentError) as error:
                raise AmberToolsInputError(f"Invalid provided charges: {error}") from error
            charges = {site_id: record.charge
                       for site_id, record in provided_result.assignments.items()}
        else:
            charges = None
        prepared = prepare_input(system, charges)
        chain = _discover(options)
        if options.work_root is not None:
            options.work_root.mkdir(parents=True, exist_ok=True)
        directory = Path(tempfile.mkdtemp(
            prefix="island-ambertools-", dir=options.work_root,
        ))
        try:
            return self._prepare_in_directory(
                system, options, prepared, charges, chain, directory,
            )
        except AmberToolsStageError:
            raise  # Failed job artifacts stay in this unique directory.
        except (AmberToolsInputError, AmberImportError, ValueError, KeyError, OSError) as error:
            raise AmberToolsStageError(
                "validation/import", str(error), artifact_dir=str(directory),
            ) from error
        finally:
            # Only clean a successful job. Failed jobs are always retained.
            if (not options.retain_success_artifacts
                and (directory / "ISLAND_SUCCESS").is_file()):
                shutil.rmtree(directory)

    def _prepare_in_directory(
        self, system: MolecularSystem, options: AmberToolsOptions,
        prepared: object, charges: dict[int, float] | None,
        chain: _Toolchain, directory: Path,
    ) -> AmberToolsPreparationResult:
        names = prepared.names
        environment = {**os.environ, "AMBERHOME": str(chain.amberhome)}
        (directory / "input.mol2").write_text(prepared.mol2_text)
        from json import dumps

        (directory / "lineage.json").write_text(dumps({
            "schema": "island_ambertools_lineage_v1",
            "input_name_to_site_id": names,
            "input_index_to_site_id": {
                str(i): site_id for i, site_id in enumerate(sorted(system.topology.sites), 1)
            },
        }, sort_keys=True, indent=2))
        command = [
            chain.executables["antechamber"], "-i", "input.mol2", "-fi", "mol2",
            "-o", "typed.mol2", "-fo", "mol2", "-at", options.force_field,
            "-c", "rc" if charges is not None else "bcc",
            "-nc", str(prepared.formal_charge), "-m", "1", "-s", "2",
            "-j", "4", "-du", "yes", "-pf", "no",
        ]
        if charges is not None:
            (directory / "charges.txt").write_text("\n".join(
                f"{charges[names[name]]:.10f}" for name in names
            ) + "\n")
            command.extend(("-cf", "charges.txt"))
        stages = [_run_stage(
            "antechamber", command, directory, options.timeout_seconds,
            environment=environment,
        )]
        typed = directory / "typed.mol2"
        _require_artifact("antechamber", typed, directory)
        try:
            typed_index, typed_charges, typed_positions = parse_and_validate_mol2(
                typed, system, names, prepared.expected_cip, stage="antechamber",
            )
        except (ValueError, KeyError, AmberToolsInputError) as error:
            raise AmberToolsStageError(
                "antechamber", str(error), artifact_dir=str(directory),
            ) from error
        if charges is not None:
            _check_charge_values(charges, typed_charges, "antechamber", directory)
            if (directory / "sqm.out").exists():
                raise AmberToolsStageError(
                    "antechamber", "provided-charge mode unexpectedly ran SQM",
                    artifact_dir=str(directory),
                )
            charge_outcome = {"mode": "provided", "qm_run": False,
                              "serialization_tolerance_e": SERIALIZATION_TOLERANCE}
        else:
            sqm = directory / "sqm.out"
            _require_artifact("antechamber", sqm, directory)
            sqm_text = sqm.read_text(errors="replace")
            if (not SQM_SUCCESS.search(sqm_text) or FATAL.search(sqm_text)
                or SQM_FAILURE.search(sqm_text)):
                raise AmberToolsStageError(
                    "antechamber", "AM1-BCC SQM convergence/completion is unverified",
                    artifact_dir=str(directory),
                    stdout=sqm_text[-4000:],
                )
            charge_outcome = {"mode": "am1bcc", "qm_run": True,
                              "sqm_out_sha256": _sha(sqm),
                              "convergence_marker": SQM_SUCCESS.search(sqm_text).group(0),
                              "sqm_warning_lines": [line for line in sqm_text.splitlines()
                                                    if "warn" in line.lower()]}
        parmchk_command = [
            chain.executables["parmchk2"], "-i", "typed.mol2", "-f", "mol2",
            "-o", "typed.frcmod", "-s",
            "1" if options.force_field == "gaff" else "2",
        ]
        stages.append(_run_stage(
            "parmchk2", parmchk_command, directory, options.timeout_seconds,
            environment=environment,
        ))
        frcmod = directory / "typed.frcmod"
        _require_artifact("parmchk2", frcmod, directory)
        frcmod_text = frcmod.read_text(errors="replace")
        if UNRESOLVED.search(frcmod_text):
            raise AmberToolsStageError(
                "parmchk2", "unresolved ATTN, need revision or missing parameter",
                artifact_dir=str(directory), stdout=frcmod_text[-4000:],
            )
        estimates = [line.strip() for line in frcmod_text.splitlines()
                     if any(token in line.lower() for token in (
                         "analogy", "estimated", "empirical",
                     ))]
        leaprc = str(chain.leaprc).replace('"', '\\"')
        script = (
            f'source "{leaprc}"\n'
            "loadamberparams typed.frcmod\n"
            "mol = loadmol2 typed.mol2\n"
            "check mol\n"
            "saveamberparm mol result.prmtop result.rst7\n"
            "quit\n"
        )
        (directory / "leap.in").write_text(script)
        stages.append(_run_stage(
            "tleap", [chain.executables["tleap"], "-f", "leap.in"],
            directory, options.timeout_seconds,
            environment=environment,
        ))
        prmtop = directory / "result.prmtop"
        restart = directory / "result.rst7"
        _require_artifact("tleap", prmtop, directory)
        _require_artifact("tleap", restart, directory)
        leap_log = directory / "leap.log"
        _require_artifact("tleap", leap_log, directory)
        log = leap_log.read_text(errors="replace")
        if UNRESOLVED.search(log) or NONZERO_ERRORS.search(log) or any(
            FATAL.search(line) for line in log.splitlines()
            if "0 errors" not in line.lower()
        ):
            raise AmberToolsStageError(
                "tleap", "LEaP reported unresolved parameters or errors",
                artifact_dir=str(directory), stdout=log[-4000:],
            )
        mapping, final_charges = _prmtop_lineage(
            prmtop, restart, system, names, prepared.expected_cip, directory,
        )
        _check_charge_values(typed_charges, final_charges, "tleap", directory)
        if charges is not None:
            _check_charge_values(charges, final_charges, "tleap", directory)
        try:
            imported = import_amber_prmtop(
                system, prmtop, mapping,
                source=f"AmberTools {options.force_field} generated prmtop",
                force_field=options.force_field,
                charge_method="provided" if charges is not None else "AM1-BCC",
                charge_tolerance=options.charge_tolerance,
            )
        except (AmberImportError, IncompleteChargeAssignmentError) as error:
            raise AmberToolsStageError(
                "import", str(error), artifact_dir=str(directory),
            ) from error
        input_coordinates = {
            str(site_id): system.coordinates.get(site_id).tolist()
            for site_id in sorted(system.topology.sites)
        }
        record = {
            "schema": PREPARATION_SCHEMA,
            "engine_version": self.engine_version,
            "requested_force_field": options.force_field,
            "charge_method": options.charge_method,
            "charge_outcome": charge_outcome,
            "input_coordinate_signature": digest(input_coordinates),
            "input_coordinates_angstrom": input_coordinates,
            "expected_cip_by_site": prepared.expected_cip,
            "antechamber_coordinates_angstrom": {
                str(site_id): list(position)
                for site_id, position in sorted(typed_positions.items())
            },
            "input_mol2_sha256": _sha(directory / "input.mol2"),
            "input_lineage_sha256": _sha(directory / "lineage.json"),
            "lineage": {
                "input_name_to_site_id": names,
                "typed_mol2_index_to_site_id": typed_index,
                "prmtop_index_to_site_id": mapping,
            },
            "tool_versions": chain.versions,
            "tool_package": chain.package,
            "executables": chain.executables,
            "executable_sha256": {
                name: _sha(Path(path)) for name, path in chain.executables.items()
            },
            "amberhome": str(chain.amberhome),
            "force_field_data": {
                "path": str(chain.data_file), "header": _header(chain.data_file),
                "sha256": _sha(chain.data_file),
            },
            "leaprc": {"path": str(chain.leaprc), "sha256": _sha(chain.leaprc)},
            "stages": stages,
            "parmchk2_estimation_lines": estimates,
            "warnings": [line for line in log.splitlines()
                         if "warn" in line.lower()],
            "artifact_sha256": {name: _sha(directory / name) for name in (
                "typed.mol2", "typed.frcmod", "result.prmtop", "result.rst7",
                "leap.in", "leap.log",
            )},
            "artifact_dir": str(directory) if options.retain_success_artifacts else None,
            "charge_validation_tolerance_e": options.charge_tolerance,
            "provided_charge_tolerance_e": SERIALIZATION_TOLERANCE,
        }
        imported = replace(imported, provenance={
            **dict(imported.provenance), "ambertools_preparation": record,
        }, result_signature="")
        imported = replace(imported, result_signature=imported.content_signature())
        imported.validate_integrity(system)
        record = {**record, "imported_result_signature": imported.result_signature}
        result = AmberToolsPreparationResult(imported, record, digest(record))
        result.validate_integrity(system)
        (directory / "ISLAND_SUCCESS").write_text(result.record_signature + "\n")
        return result
