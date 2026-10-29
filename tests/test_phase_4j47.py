"""Real pinned external PDMS rows, static execution and fail-closed MD gates."""

import os
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import numpy as np
import pytest

from island.charge_references.records import pack, unpack
from island.core import SimulationBox
from island.forcefields.pdms_lammps import (
    ForceFieldRequest,
    OPLSOptions,
    Preparation,
    prepare_forcefield,
)
from island.forcefields.pdms_opls import PROFILE, PDMSOPLSError, export_lammps
from island.pdms_lammps import LAMMPSOptions, pin, reconstruct, run_md, run_static, seal
from scripts.qualify_pdms_public import molecule

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "island_ffref/pcff/PDMS_benchmark-main.zip"
ARCHIVE_PIN = "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4"
LMP = ROOT.parent / "island-validation/phase4j22/lmp-kspace"
LMP_PIN = "6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191"


def fixture(dp):
    system, anchors = molecule(dp)
    xyz = np.array([system.coordinates.get(i) for i in system.topology.sites])
    system.coordinates.translate(20 - xyz.mean(axis=0))
    return system, ForceFieldRequest(
        "opls", OPLSOptions(PROFILE, ARCHIVE, ARCHIVE_PIN, anchors)
    )


@pytest.fixture(scope="module")
def assets():
    if not ARCHIVE.is_file() or not LMP.is_file():
        pytest.skip("Explicit local research ZIP and pinned LAMMPS required")
    assert pin(ARCHIVE) == ARCHIVE_PIN and pin(LMP) == LMP_PIN
    return LAMMPSOptions(LMP, LMP_PIN)


@pytest.fixture(scope="module")
def workflow(assets, tmp_path_factory):
    system, request = fixture(2)
    prepared = prepare_forcefield(system, request)
    root = tmp_path_factory.mktemp("j47") / "workflow"
    result = run_static(system, prepared, root, engine=assets)
    assert result["status"] == "static_executed_unqualified", result
    return system, request, prepared, root, result


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_domain_identity_and_report(assets, dp):
    system, request = fixture(dp)
    prepared = prepare_forcefield(system, request)
    p = prepared.payload
    backend = prepared.validate_integrity(system)
    assert p["status"] == "prepared_static_only"
    assert p["profile_identity"] == PROFILE
    assert p["ordered_charge_literals"] == [
        r["charge_e"] for r in backend.payload["mapping"]
    ]
    assert p["assigned_interactions"] == backend.payload["interactions"]
    assert p["charge_origin"]["normalization"] is False
    assert p["charge_origin"]["redistribution"] is False
    assert not p["stages"]["executable_MD"]
    assert p["stages"]["parameter_complete"]
    assert p["production_validated"] is False
    with pytest.raises(FrozenInstanceError):
        request.options.profile = "native_opls"


@pytest.mark.parametrize(
    "field,value",
    [
        ("profile", "native_oplsaa"),
        ("lj_mixing", "arithmetic"),
        ("special_lj", (0, 0, 1)),
        ("special_coulomb", (0, 0, 1)),
        ("lj_inner_angstrom", 8),
        ("cutoff_angstrom", 12),
        ("reciprocal_method", "ewald"),
        ("reciprocal_tolerance", 1e-5),
        ("expected_sha256", "0" * 64),
        ("expected_sha256", None),
        ("archive_path", "/missing/pdms.zip"),
        ("use_scope", "production"),
    ],
)
def test_convention_and_asset_negative_controls(assets, field, value):
    system, request = fixture(2)
    changed = replace(request, options=replace(request.options, **{field: value}))
    prepared = prepare_forcefield(system, changed)
    assert prepared.payload["status"] in (
        "typing_blocked",
        "parameter_incomplete",
        "asset_unavailable",
    )
    with pytest.raises(PDMSOPLSError) as error:
        prepared.validate_integrity(system)
    assert error.value.result is prepared
    assert error.value.diagnostics == prepared.payload["diagnostics"]


def test_pcff_family_rejected(assets):
    system, request = fixture(0)
    p = prepare_forcefield(system, replace(request, family="pcff"))
    assert p.payload["status"] == "typing_blocked"
    assert (
        p.payload["diagnostics"][0]["stage"] == "external_opls_profile_family_mismatch"
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "box",
        "coordinate",
        "order",
        "charge",
        "mass",
        "status",
        "rows",
        "source",
        "profile",
    ],
)
def test_stale_and_forged_assignment(assets, mutation):
    system, request = fixture(2)
    p = prepare_forcefield(system, request)
    value = p.payload
    if mutation == "box":
        system.box = SimulationBox(41, 40, 40)
    elif mutation == "coordinate":
        system.coordinates.translate((0.001, 0, 0))
    elif mutation == "order":
        sites = system.topology.sites
        atom = sites.pop(next(iter(sites)))
        sites[atom.id] = atom
    elif mutation == "charge":
        value["ordered_charge_literals"][0] = "0"
    elif mutation == "mass":
        value["source_identity"] = "different-mass"
    elif mutation == "status":
        value["stages"]["executable_MD"] = True
    elif mutation == "rows":
        value["assigned_interactions"]["angle"].pop()
    elif mutation == "source":
        value["archive_sha256"] = "0" * 64
    else:
        value["profile_identity"] = "native_oplsaa"
    if mutation not in ("box", "coordinate", "order"):
        p = Preparation(pack(value), request)
    with pytest.raises(PDMSOPLSError):
        p.validate_integrity(system)


def test_source_preparation_charge_and_row_forgery(assets):
    from island.forcefields.pdms_opls import PDMSOPLSPreparation

    system, request = fixture(2)
    backend = prepare_forcefield(system, request).validate_integrity(system)
    for field in ("charge", "rows"):
        p = backend.payload
        if field == "charge":
            p["mapping"][0]["charge_e"] = "0"
        else:
            p["interactions"]["proper_torsion"].pop()
        with pytest.raises(PDMSOPLSError):
            PDMSOPLSPreparation(pack(p), backend.options).require_executable(system)


def test_missing_source_rows_and_malformed_zip(assets, tmp_path):
    from zipfile import ZipFile

    system, request = fixture(2)
    malformed = tmp_path / "malformed.zip"
    malformed.write_bytes(b"not a zip")
    incomplete = tmp_path / "incomplete.zip"
    with ZipFile(ARCHIVE) as old, ZipFile(incomplete, "w") as new:
        for name in old.namelist():
            raw = old.read(name)
            if name.endswith("opls_system.in.settings"):
                raw = raw.replace(b"angle_coeff", b"removed_row", 1)
            new.writestr(name, raw)
    for path in (malformed, incomplete):
        options = replace(request.options, archive_path=path, expected_sha256=pin(path))
        p = prepare_forcefield(system, replace(request, options=options))
        assert p.payload["status"] == (
            "asset_unavailable" if path == malformed else "parameter_incomplete"
        )
        with pytest.raises(PDMSOPLSError):
            p.validate_integrity(system)


def test_actual_static_execution_and_no_md(workflow):
    _, _, _, root, result = workflow
    p = result["observed"]
    assert len(p["forces_kj_mol_angstrom"]) == 37
    assert len(p["virial_kj_mol_xx_yy_zz_xy_xz_yz"]) == 6
    assert np.isfinite(p["forces_kj_mol_angstrom"]).all()
    assert result["MD_executed"] is result["NPT_executed"] is False
    script = (root / "lammps/static.in").read_text()
    assert "kspace_style pppm 0.0001" in script
    assert "pair_modify table 0" not in script
    assert "minimize" not in script and "fix " not in script
    assert "run 0\n" in script
    assert not list(root.rglob("*.restart")) and not list(root.rglob("*.chk"))


def test_deterministic_export(assets, tmp_path):
    system, request = fixture(4)
    backend = prepare_forcefield(system, request).validate_integrity(system)
    assert export_lammps(system, backend, tmp_path / "a") == export_lammps(
        system, backend, tmp_path / "b"
    )


def test_md_guard_precedes_engine_and_files(assets, tmp_path, monkeypatch):
    import island.pdms_lammps as module

    system, request = fixture(2)
    p = prepare_forcefield(system, request)

    def forbidden(*a, **kw):
        raise AssertionError("must not instantiate/execute engine")

    monkeypatch.setattr(module, "_execute", forbidden)
    monkeypatch.setattr(module, "engine_record", forbidden)
    target = tmp_path / "MD-forbidden"
    with pytest.raises(PDMSOPLSError, match="lammps_md_acceptance_incomplete") as error:
        run_md(system, p, target, engine=assets, steps=10)
    assert error.value.result is p
    assert not target.exists()


@pytest.mark.parametrize(
    "field",
    [
        "system.json",
        "assignment.json",
        "source-preparation.json",
        "lammps/system.data",
        "lammps/coefficients.in",
        "lammps/static.in",
        "engine.json",
        "result.json",
    ],
)
def test_workflow_mutation(workflow, assets, tmp_path, field):
    _, request, _, root, result = workflow
    moved = tmp_path / "relocated"
    shutil.copytree(root, moved)
    path = moved / field
    path.write_bytes(path.read_bytes() + b"tamper")
    with pytest.raises(PDMSOPLSError, match="workflow_file_mutation"):
        reconstruct(
            moved,
            manifest_sha256=result["manifest_sha256"],
            request=request,
            engine=assets,
        )


@pytest.mark.parametrize(
    "field", ["lammps/static.in", "result.json", "assignment.json"]
)
def test_resigned_evidence_cannot_grant_authority(workflow, assets, tmp_path, field):
    _, request, _, root, _ = workflow
    moved = tmp_path / "re-signed"
    shutil.copytree(root, moved)
    path = moved / field
    if field.endswith(".in"):
        path.write_text(path.read_text() + "run 10\n")
    else:
        p = unpack(path.read_text())
        if field == "result.json":
            p["MD_executed"] = True
        else:
            p["status"] = "prepared"
        path.write_text(pack(p))
    new_pin = seal(moved)
    with pytest.raises(PDMSOPLSError):
        reconstruct(moved, manifest_sha256=new_pin, request=request, engine=assets)


def test_wrong_engine_hash_rejected_before_output(assets, tmp_path):
    system, request = fixture(0)
    prepared = prepare_forcefield(system, request)
    with pytest.raises(PDMSOPLSError, match="lammps_executable_hash_mismatch"):
        run_static(
            system,
            prepared,
            tmp_path / "forbidden",
            engine=replace(assets, expected_sha256="0" * 64),
        )
    assert not (tmp_path / "forbidden").exists()


def test_fresh_process_relocated_reconstruction(workflow, assets, tmp_path):
    _, request, prepared, root, result = workflow
    moved = tmp_path / "relocated"
    shutil.copytree(root, moved)
    shutil.copyfile(ARCHIVE, tmp_path / "source.zip")
    code = """import sys
from island.forcefields.pdms_lammps import *
from island.pdms_lammps import *
from island.forcefields.pdms_opls import PROFILE
r=ForceFieldRequest("opls",OPLSOptions(PROFILE,sys.argv[3],sys.argv[4],(int(sys.argv[7]),)))
s,p,v=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],request=r,engine=LAMMPSOptions(sys.argv[5],sys.argv[6]))
print(p.payload["identity"])
"""
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(moved),
            result["manifest_sha256"],
            str(tmp_path / "source.zip"),
            ARCHIVE_PIN,
            str(LMP),
            LMP_PIN,
            str(request.options.e_oxygen_site_ids[0]),
        ],
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src"), "OMP_NUM_THREADS": "1"},
    )
    assert proc.stdout.strip() == prepared.payload["identity"]


def test_immutable_thresholds():
    from island.pdms_lammps import ENERGY_BOUND, FORCE_BOUND
    from scripts.qualify_pdms_public import ENERGY_BOUND as old_e
    from scripts.qualify_pdms_public import FORCE_BOUND as old_f

    assert ENERGY_BOUND == old_e and FORCE_BOUND == old_f


def test_component_comparison_is_json_and_fail_closed(workflow, tmp_path):
    pytest.importorskip("openmm")
    import json

    from scripts.phase_4j47 import compare, write

    system, _, prepared, _, result = workflow
    report = compare(system, prepared, result["observed"])
    write(tmp_path / "comparison.json", report)
    saved = json.loads((tmp_path / "comparison.json").read_text())
    from island.forcefields.external_pdms import digest

    assert digest(report) == digest(saved)
    assert saved["acceptance_passed"] is False
    assert saved["identical_physical_inputs"] is True
    assert saved["identical_numerical_long_range_settings"] is False
    assert saved["MD_executed"] is False
    assert saved["component_comparison"]["virial"]["passed"] is False
    assert set(saved["component_comparison"]) == {
        "bond",
        "angle",
        "proper_torsion",
        "lj",
        "total",
        "coulomb",
        "coulomb_real_corrected",
        "coulomb_reciprocal_self_background",
        "improper",
        "lj14",
        "coulomb14",
        "forces",
        "virial",
        "pressure_proxy",
    }


def test_executable_relocation_keeps_build_identity(assets, tmp_path):
    from island.pdms_lammps import engine_record

    copied = tmp_path / "relocated-lammps"
    shutil.copyfile(LMP, copied)
    copied.chmod(0o755)
    assert engine_record(assets) == engine_record(LAMMPSOptions(copied, LMP_PIN))


def test_malformed_workflow_record_is_classified(workflow, assets, tmp_path):
    _, request, _, root, _ = workflow
    moved = tmp_path / "malformed"
    shutil.copytree(root, moved)
    (moved / "system.json").write_text('{"payload":"broken","sha256":"forged"}')
    checksum = seal(moved)
    with pytest.raises(PDMSOPLSError) as error:
        reconstruct(moved, manifest_sha256=checksum, request=request, engine=assets)
    assert error.value.status == "parameter_incomplete"


@pytest.mark.parametrize("mutation", ["charge", "row", "units"])
def test_independent_raw_auditor_rejects_resigned_mapping(workflow, tmp_path, mutation):
    from scripts.audit_phase_4j47_inputs import audit, digest

    _, _, _, root, _ = workflow
    moved = tmp_path / "dp2"
    shutil.copytree(root, moved)
    path = moved / "source-preparation.json"
    p = unpack(path.read_text())
    if mutation == "charge":
        p["mapping"][0]["charge_e"] = "-1.2500000000000000e+00"
    else:
        row = p["source"]["rows"][0]
        if mutation == "row":
            row["values"][0] = "999"
        else:
            row["units"] = "SI"
        row["evidence_checksum"] = digest(
            {k: v for k, v in row.items() if k != "evidence_checksum"}
        )
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    path.write_text(pack(p))
    with pytest.raises(ValueError, match="independent_"):
        audit(tmp_path, ARCHIVE, ARCHIVE_PIN, dps=(2,))


def test_missing_assets_cli_is_classified(tmp_path, capsys):
    import json

    from scripts.phase_4j47 import main

    rc = main(
        [
            "--archive",
            str(tmp_path / "absent.zip"),
            "--archive-sha256",
            ARCHIVE_PIN,
            "--lammps",
            str(LMP),
            "--lammps-sha256",
            LMP_PIN,
            "--lammps-source",
            str(tmp_path / "absent-src"),
            "--paper",
            str(tmp_path / "absent.pdf"),
            "--paper-sha256",
            "0" * 64,
            "--output",
            str(tmp_path / "forbidden"),
        ]
    )
    assert rc == 1
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "asset_unavailable" and report["MD_executed"] is False
    assert not (tmp_path / "forbidden").exists()


@pytest.mark.parametrize(
    "field,value",
    [
        ("profile_identity", "native_pcff"),
        ("assignment_identity", "forged"),
        ("executable_sha256", "0" * 64),
        ("status", "scientifically_validated"),
    ],
)
def test_resigned_result_identity_cannot_relabel_model(
    workflow, assets, tmp_path, field, value
):
    from island.forcefields.external_pdms import digest

    _, request, _, root, _ = workflow
    moved = tmp_path / "relabeled"
    shutil.copytree(root, moved)
    path = moved / "result.json"
    p = unpack(path.read_text())
    p[field] = value
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    path.write_text(pack(p))
    checksum = seal(moved)
    with pytest.raises(PDMSOPLSError):
        reconstruct(moved, manifest_sha256=checksum, request=request, engine=assets)
