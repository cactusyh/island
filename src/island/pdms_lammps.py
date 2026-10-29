"""Checksummed static ISLAND → pinned LAMMPS execution, never implicit MD.

Raw engine output is preserved. A successful process is not numerical parity.
Relocation revalidates the supplied archive and regenerates every input byte.
"""

import os
import platform
import re
import subprocess
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from island.charge_references.records import pack, unpack
from island.exceptions import IslandError
from island.forcefields.external_pdms import LIMITS, digest
from island.forcefields.pdms_lammps import ForceFieldRequest, prepare_forcefield
from island.forcefields.pdms_opls import PDMSOPLSError, check, export_lammps
from island.pdms import restore_pdms_system
from island.workflows.bundle import system_data

COMPONENTS = (
    "bond",
    "angle",
    "proper_torsion",
    "improper",
    "lj",
    "coulomb_real_corrected",
    "coulomb_reciprocal_self_background",
    "total",
)
# Existing qualification limits, in kJ/mol and kJ/(mol angstrom).
ENERGY_BOUND = 3e-4
FORCE_BOUND = 3e-4
STATIC = """compute virial all pressure NULL virial
thermo 1
thermo_style custom step ebond eangle edihed eimp evdwl ecoul elong pe c_virial[1] c_virial[2] c_virial[3] c_virial[4] c_virial[5] c_virial[6]
thermo_modify format float %.17g
run 0
write_dump all custom state.dump id type mol q x y z fx fy fz modify sort id format float %.17g
"""


def pin(path):
    return sha256(Path(path).read_bytes()).hexdigest()


@dataclass(frozen=True)
class LAMMPSOptions:
    executable: object
    expected_sha256: str

    def validate(self):
        try:
            path = Path(self.executable).resolve(strict=True)
            check(
                re.fullmatch("[0-9a-f]{64}", self.expected_sha256 or "") is not None,
                "independent_lammps_hash_required",
                "asset_unavailable",
            )
            check(
                pin(path) == self.expected_sha256,
                "lammps_executable_hash_mismatch",
                "asset_unavailable",
            )
            return path
        except (OSError, TypeError, ValueError) as error:
            raise PDMSOPLSError(
                "asset_unavailable", "lammps_unavailable", reason=str(error)
            ) from error


def _execute(engine, root, input_name):
    executable = engine.validate()
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    try:
        proc = subprocess.run(
            [
                str(executable),
                "-in",
                input_name,
                "-log",
                "engine.log",
                "-screen",
                "engine.screen",
            ],
            cwd=root,
            env=env,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PDMSOPLSError(
            "converter_unavailable", "lammps_execution_failed", reason=str(error)
        ) from error
    (root / "process.stdout").write_text(proc.stdout)
    (root / "process.stderr").write_text(proc.stderr)
    return proc.returncode


def engine_record(engine):
    path = engine.validate()
    try:
        proc = subprocess.run(
            [str(path), "-h"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PDMSOPLSError(
            "asset_unavailable", "lammps_version_failed", reason=str(error)
        ) from error
    check(proc.returncode == 0, "lammps_version_failed", "asset_unavailable")
    help_text = proc.stdout.replace(str(path), "<pinned-lammps>")
    return {
        "executable_sha256": engine.expected_sha256,
        "help_text": help_text,
        "help_sha256": sha256(help_text.encode()).hexdigest(),
        "version_revision_build_packages": "see complete checksummed help_text",
        "platform": platform.platform(),
        "machine": platform.machine(),
        "units": "real",
        "energy_conversion": "kcal/mol * 4.184 = kJ/mol",
        "force_conversion": "kcal/(mol angstrom) * 4.184",
        "virial_conversion": "atm * angstrom^3 * 4.184/68568.415",
        "precision": "see help_text FFT precision and compile flags; output %.17g",
        "processes": 1,
        "OMP_NUM_THREADS": 1,
    }


def parse_static(root, preparation):
    root = Path(root)
    lines = (root / "engine.screen").read_text().splitlines()
    header = next(
        (
            n
            for n, line in enumerate(lines)
            if line.strip().startswith("Step ") and "E_bond" in line
        ),
        None,
    )
    check(header is not None, "lammps_thermo_header_missing")
    values = [float(v) for v in lines[header + 1].split()]
    check(
        len(values) == 15 and np.isfinite(values).all(),
        "lammps_nonfinite_or_incomplete_thermo",
    )
    energy = dict(zip(COMPONENTS, (v * 4.184 for v in values[1:9])))
    p = preparation.payload
    dump = (root / "state.dump").read_text().splitlines()
    start = next(n for n, s in enumerate(dump) if s.startswith("ITEM: ATOMS"))
    check(
        dump[start].split()[2:]
        == ["id", "type", "mol", "q", "x", "y", "z", "fx", "fy", "fz"],
        "dump_columns_changed",
    )
    rows = [[float(v) for v in line.split()] for line in dump[start + 1 :]]
    check(len(rows) == len(p["mapping"]), "dump_atom_count_changed")
    for n, (row, atom) in enumerate(zip(rows, p["mapping"]), 1):
        check(
            row[:4]
            == [
                n,
                {"H": 1, "O": 2, "Si": 3, "C": 4}[atom["element"]],
                atom["molecule_id"],
                float(atom["charge_e"]),
            ],
            "dump_atom_type_order_or_charge_mismatch",
        )
    forces = np.array([row[7:] for row in rows]) * 4.184
    volume = np.prod(
        [
            float(dump[5 + i].split()[1]) - float(dump[5 + i].split()[0])
            for i in range(3)
        ]
    )
    stress = np.array(values[9:])
    virial = stress * volume * 4.184 / 68568.415
    return {
        "energy_kj_mol": energy,
        "forces_kj_mol_angstrom": forces.tolist(),
        "coordinates_angstrom": [row[4:7] for row in rows],
        "virial_kj_mol_xx_yy_zz_xy_xz_yz": virial.tolist(),
        "configurational_pressure_proxy_atm": stress.tolist(),
        "pressure_status": "diagnostic_configurational_proxy_only; backend_pressure_unavailable",
        "force_field_pressure_validated": False,
    }


def seal(root):
    root = Path(root)
    files = {
        str(p.relative_to(root)): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file() and p != root / "manifest.json"
    }
    value = {"schema": "island_pdms_lammps_workflow_manifest_v1", "files": files}
    (root / "manifest.json").write_text(pack({**value, "identity": digest(value)}))
    return pin(root / "manifest.json")


def run_static(system, prepared, directory, *, engine):
    """Source-authored PPPM 1e-4. No reciprocal, charge or pair overrides."""
    backend = prepared.validate_integrity(system)
    info = engine_record(engine)  # explicit engine pin checked before creating files
    root = Path(directory)
    check(not root.exists(), "output_directory_exists")
    root.mkdir(parents=True)
    export_lammps(system, backend, root / "lammps")
    (root / "system.json").write_text(pack(system_data(system)))
    (root / "assignment.json").write_text(prepared.json_text)
    (root / "source-preparation.json").write_text(backend.json_text)
    (root / "engine.json").write_text(pack(info))
    work = root / "lammps"
    (work / "static.in").write_text((work / "system.in").read_text() + STATIC)
    code = _execute(engine, work, "static.in")
    if code == 0:
        try:
            observed = parse_static(work, backend)
            result = {
                "status": "static_executed_unqualified",
                "observed": observed,
                "diagnostics": [],
            }
        except (ValueError, IndexError, StopIteration, OSError, PDMSOPLSError) as error:
            result = {
                "status": "parameter_incomplete",
                "diagnostics": [
                    {"stage": "lammps_output_invalid", "reason": str(error)}
                ],
            }
    else:
        result = {
            "status": "converter_unavailable",
            "diagnostics": [{"stage": "lammps_process_failed", "returncode": code}],
        }
    result.update(
        {
            "profile_identity": backend.payload["profile_identity"],
            "assignment_identity": prepared.payload["identity"],
            "executable_sha256": engine.expected_sha256,
            "MD_executed": False,
            "NPT_executed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "limitations": LIMITS,
        }
    )
    (root / "result.json").write_text(pack({**result, "identity": digest(result)}))
    return {**result, "manifest_sha256": seal(root)}


def reconstruct(directory, *, manifest_sha256, request, engine):
    """Never trust JSON as authority: reprepare and regenerate source-bound inputs."""
    root = Path(directory).resolve()
    try:
        check(
            pin(root / "manifest.json") == manifest_sha256,
            "workflow_manifest_hash_mismatch",
        )
        manifest = unpack((root / "manifest.json").read_text())
        check(
            manifest["schema"] == "island_pdms_lammps_workflow_manifest_v1"
            and manifest["identity"]
            == digest({k: v for k, v in manifest.items() if k != "identity"}),
            "workflow_manifest_identity_mismatch",
        )
        names = set(manifest["files"])
        actual = {
            str(p.relative_to(root))
            for p in root.rglob("*")
            if p.is_file() and p != root / "manifest.json"
        }
        check(actual == names, "workflow_membership_changed")
        for name, expected in manifest["files"].items():
            path = root / name
            check(
                not Path(name).is_absolute()
                and ".." not in Path(name).parts
                and path.resolve().is_relative_to(root)
                and not path.is_symlink(),
                "invalid_workflow_path",
            )
            check(pin(path) == expected, "workflow_file_mutation", file=name)
        system = restore_pdms_system(unpack((root / "system.json").read_text()))
        check(type(request) is ForceFieldRequest, "explicit_request_required")
        prepared = prepare_forcefield(system, request)
        backend = prepared.validate_integrity(system)
        check(
            prepared.json_text == (root / "assignment.json").read_text(),
            "stale_assignment",
        )
        check(
            backend.json_text == (root / "source-preparation.json").read_text(),
            "stale_source_preparation",
        )
        engine.validate()
        recorded = unpack((root / "engine.json").read_text())
        check(recorded == engine_record(engine), "engine_build_or_platform_changed")
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "lammps"
            export_lammps(system, backend, target)
            for name in ("system.data", "system.in", "coefficients.in"):
                check(
                    (target / name).read_bytes()
                    == (root / "lammps" / name).read_bytes(),
                    "regenerated_input_mismatch",
                    file=name,
                )
            check(
                (root / "lammps/static.in").read_text()
                == (target / "system.in").read_text() + STATIC,
                "static_protocol_mutation",
            )
        report = unpack((root / "result.json").read_text())
        check(
            report["identity"]
            == digest({k: v for k, v in report.items() if k != "identity"}),
            "result_identity_mismatch",
        )
        check(
            report["MD_executed"] is False
            and report["NPT_executed"] is False
            and report["production_validated"] is False
            and report["simulation_readiness"] == "not_established"
            and report["limitations"] == LIMITS,
            "forged_workflow_authority",
        )
        check(
            report["profile_identity"] == backend.payload["profile_identity"]
            and report["assignment_identity"] == prepared.payload["identity"]
            and report["executable_sha256"] == engine.expected_sha256,
            "result_source_profile_or_engine_mismatch",
        )
        check(
            report["status"]
            in (
                "static_executed_unqualified",
                "parameter_incomplete",
                "converter_unavailable",
            ),
            "forged_workflow_status",
        )
        if report["status"] == "static_executed_unqualified":
            check(report["diagnostics"] == [], "forged_static_diagnostics")
            check(
                report["observed"] == parse_static(root / "lammps", backend),
                "result_output_mismatch",
            )
        return system, prepared, report
    except PDMSOPLSError:
        raise
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        IslandError,
    ) as error:
        raise PDMSOPLSError(
            "parameter_incomplete", "workflow_reconstruction_failed", reason=str(error)
        ) from error


def run_md(system, prepared, directory, *, engine, **kwargs):
    """Fail before engine invocation/files/checkpoints while J47 gates are open."""
    prepared.require_md(system)
