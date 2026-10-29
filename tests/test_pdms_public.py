"""Real source-bound external PDMS preparation; no artificial physical rows."""

import os
import shutil
import subprocess
import sys
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.core import SimulationBox
from island.forcefields.pdms_opls import (
    PDMSOPLSPreparation,
    export_lammps,
    save_preparation,
)
from island.forcefields.research import (
    ForceFieldRequest,
    PDMSOPLSError,
    PDMSOPLSOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.pdms import build_openmm_model, finite_difference

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "island_ffref/pcff/PDMS_benchmark-main.zip"
PIN = "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4"


def molecule(dp):
    if dp == 0:
        system = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
        anchors = ()
    else:
        system = build_linear_polymer(
            "[*:1]O[Si](C)(C)[*:2]",
            dp=dp,
            head_end_group="[*][Si](C)(C)C",
            tail_end_group="[*]C",
            add_hydrogens=True,
            random_seed=2026,
        )
        anchors = (system.metadata["polymer"]["head_site_id"],)
    system.box = SimulationBox(40, 40, 40)
    return system, anchors


@pytest.fixture(scope="module")
def options():
    if not ARCHIVE.exists():
        pytest.skip("Explicit user-supplied local archive required")
    return PDMSOPLSOptions(ARCHIVE, PIN)


def prepare(system, options, anchors=()):
    return prepare_forcefield(
        system,
        ForceFieldRequest(
            "external_pdms_opls", replace(options, e_oxygen_site_ids=anchors)
        ),
    )


@pytest.mark.parametrize("dp,atoms", [(0, 27), (1, 27), (2, 37), (4, 57)])
def test_complete_graph_bound_assignment(options, dp, atoms):
    system, anchors = molecule(dp)
    result = prepare(system, options, anchors)
    result.require_executable(system)
    p = result.payload
    assert p["status"] == "prepared"
    assert len(p["mapping"]) == atoms
    assert len(p["interactions"]["bond"]) == atoms - 1
    assert len(p["source"]["rows"]) == 17
    assert p["model"]["paper_Qeq_model"] is False
    assert p["scientific_status"] == "scientifically_unqualified"
    assert not p["limitations"]["production_validated"]
    assert sum(Fraction(a["charge_e"]) for a in p["mapping"]) == Fraction(
        p["charge_policy"]["sum_e"]
    )


@pytest.mark.parametrize(
    "defect",
    [
        "missing_bond",
        "charge",
        "radical",
        "element",
        "mass",
        "sih",
        "sioh",
        "sioc",
        "wrong_anchor",
        "ambiguous",
    ],
)
def test_unsupported_graph_is_classified(options, defect, monkeypatch):
    system, anchors = molecule(4)
    if defect == "missing_bond":
        system.topology.remove_bond(*next(iter(system.topology.bonds)))
    elif defect in ("sih", "sioh", "sioc"):
        system = from_smiles(
            {"sih": "[SiH3]O[SiH3]", "sioh": "C[Si](C)(C)O", "sioc": "CO[Si](C)(C)C"}[
                defect
            ]
        )
    elif defect == "charge":
        system.topology.sites[1].formal_charge = 1
    elif defect == "radical":
        system.topology.sites[1].metadata["radical_electrons"] = 1
    elif defect == "element":
        si = next(a for a in system.topology.sites.values() if a.element == "Si")
        si.element, si.atomic_number = "P", 15
    elif defect == "mass":
        system.topology.sites[1].mass = 99
    elif defect == "wrong_anchor":
        anchors = (True,)
    else:
        anchors = ()
    result = prepare(system, options, anchors)
    assert result.payload["status"] in ("typing_blocked", "typing_ambiguous")
    assert result.payload["profile_identity"] is None
    monkeypatch.setattr("island.pdms._mm", lambda: pytest.fail("engine reached"))
    with pytest.raises(PDMSOPLSError) as e:
        create_evaluator(system, result)
    assert e.value.result is result


@pytest.mark.parametrize("defect", ["absent", "hash_missing", "hash_wrong", "scope"])
def test_asset_errors_preserve_classification(options, defect, monkeypatch):
    system, anchors = molecule(4)
    changes = {
        "absent": {"archive_path": "/absent/pdms.zip"},
        "hash_missing": {"expected_sha256": None},
        "hash_wrong": {"expected_sha256": "0" * 64},
        "scope": {"use_scope": "redistribute"},
    }[defect]
    result = prepare(system, replace(options, **changes), anchors)
    assert result.payload["status"] == "asset_unavailable"
    monkeypatch.setattr("island.pdms._mm", lambda: pytest.fail("engine reached"))
    with pytest.raises(PDMSOPLSError) as e:
        create_evaluator(system, result)
    assert e.value.status == "asset_unavailable"


@pytest.mark.parametrize("invalid_graph", [False, True])
def test_blocked_preparation_cannot_start_md(tmp_path, monkeypatch, invalid_graph):
    from island.pdms import run_pdms_md

    system, anchors = molecule(4)
    if invalid_graph:
        system.topology.remove_bond(*next(iter(system.topology.bonds)))
    result = prepare(system, PDMSOPLSOptions(), anchors)
    monkeypatch.setattr("island.pdms._mm", lambda: pytest.fail("engine reached"))
    with pytest.raises(PDMSOPLSError):
        run_pdms_md(system, result, tmp_path / "forbidden")
    assert not (tmp_path / "forbidden").exists()


@pytest.mark.parametrize(
    "field", ["charge", "row", "units", "profile", "mixing", "status"]
)
def test_resigned_forgery_rejected(options, field):
    from island.charge_references.records import pack
    from island.forcefields.external_pdms import digest

    system, anchors = molecule(4)
    result = prepare(system, options, anchors)
    p = result.payload
    if field == "charge":
        p["mapping"][0]["charge_e"] = "0"
    elif field == "row":
        p["source"]["rows"][0]["values"][0] = "0"
    elif field == "units":
        p["source"]["rows"][0]["units"] = "kJ/mol"
    elif field == "profile":
        p["profile_identity"] = "native_oplsaa"
    elif field == "mixing":
        p["model"]["mixing"] = "lorentz_berthelot"
    else:
        p["status"] = "scientifically_validated"
    p.pop("identity")
    p["identity"] = digest(p)
    forged = PDMSOPLSPreparation(pack(p), result.options)
    with pytest.raises(PDMSOPLSError):
        forged.require_executable(system)


@pytest.mark.parametrize("mutation", ["order", "molecule", "topology", "source"])
def test_stale_identity_rejected(options, mutation, tmp_path):
    system, anchors = molecule(4)
    if mutation == "source":
        target = tmp_path / "source.zip"
        shutil.copy2(ARCHIVE, target)
        options = replace(options, archive_path=target)
    result = prepare(system, options, anchors)
    if mutation == "source":
        with target.open("ab") as f:
            f.write(b"mutation")
    elif mutation == "order":
        system.topology.sites = dict(reversed(list(system.topology.sites.items())))
    elif mutation == "molecule":
        system.topology.sites[1].metadata["chain_id"] = "stale"
    else:
        system.topology.remove_bond(*next(iter(system.topology.bonds)))
    with pytest.raises(PDMSOPLSError):
        result.require_executable(system)


def test_deterministic_export_and_source_zeros(options, tmp_path):
    system, anchors = molecule(4)
    result = prepare(system, options, anchors)
    first = export_lammps(system, result, tmp_path / "a")
    assert export_lammps(system, result, tmp_path / "b") == first
    source_zeros = [
        r
        for r in result.payload["source"]["rows"]
        if r["family"] == "proper_torsion"
        and all(Fraction(x) == 0 for x in r["values"])
    ]
    assert len(source_zeros) == 3
    assert all(r["settings_evidence"]["raw_text"] for r in source_zeros)


def test_fresh_process_relocation(options, tmp_path):
    from island.charge_references.records import pack
    from island.workflows.bundle import system_data

    system, anchors = molecule(4)
    result = prepare(system, options, anchors)
    pin = save_preparation(result, tmp_path / "prepared.json")
    (tmp_path / "system.json").write_text(pack(system_data(system)))
    shutil.copy2(ARCHIVE, tmp_path / "source.zip")
    code = """
import sys
from island.charge_references.records import unpack
from island.forcefields.pdms_opls import load_preparation, PDMSOPLSOptions
from island.pdms import restore_pdms_system
from pathlib import Path
root=Path(sys.argv[1])
system=restore_pdms_system(unpack((root/'system.json').read_text()))
options=PDMSOPLSOptions(root/'source.zip',sys.argv[3],(system.metadata['polymer']['head_site_id'],))
p=load_preparation(root/'prepared.json',sys.argv[2],system,options)
print(p.payload['identity'])
"""
    completed = subprocess.run(
        [sys.executable, "-c", code, str(tmp_path), pin, PIN],
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        capture_output=True,
        text=True,
        check=True,
    )
    assert completed.stdout.strip() == result.payload["identity"]


def test_openmm_model_and_checkpoint_exact_continuation(options):
    mm = pytest.importorskip("openmm")
    from openmm import unit

    system, anchors = molecule(4)
    prepared = prepare(system, options, anchors)
    model = build_openmm_model(prepared.payload, (40, 40, 40))

    def make():
        integrator = mm.LangevinMiddleIntegrator(300, 1, 0.0005)
        integrator.setRandomNumberSeed(2026)
        return mm.Context(
            model, integrator, mm.Platform.getPlatformByName("Reference")
        ), integrator

    a, ia = make()
    xyz = (
        np.array(
            [system.coordinates.get(r["site_id"]) for r in prepared.payload["mapping"]]
        )
        * 0.1
    )
    a.setPositions(xyz)
    a.setVelocitiesToTemperature(300, 2026)
    ia.step(2)
    checkpoint = a.createCheckpoint()
    ia.step(2)
    final = a.createCheckpoint()
    b, ib = make()
    b.loadCheckpoint(checkpoint)
    ib.step(2)
    assert b.createCheckpoint() == final
    assert np.array_equal(
        a.getState(getPositions=True)
        .getPositions(asNumpy=True)
        .value_in_unit(unit.nanometer),
        b.getState(getPositions=True)
        .getPositions(asNumpy=True)
        .value_in_unit(unit.nanometer),
    )
    del a, b


def test_workflow_resume_and_thermodynamics(options, tmp_path):
    import json
    from hashlib import sha256

    mm = pytest.importorskip("openmm")
    from island.charge_references.records import pack
    from island.pdms import resume_pdms_md
    from island.workflows.bundle import system_data

    system, anchors = molecule(4)
    prepared = prepare(system, options, anchors)
    model = build_openmm_model(prepared.payload, (40, 40, 40))
    integrator = mm.LangevinMiddleIntegrator(300, 1, 0.0005)
    integrator.setRandomNumberSeed(2026)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(
        np.array(
            [system.coordinates.get(r["site_id"]) for r in prepared.payload["mapping"]]
        )
        * 0.1
    )
    context.setVelocitiesToTemperature(300, 2026)
    integrator.step(2)
    root = tmp_path / "start"
    root.mkdir()
    save_preparation(prepared, root / "preparation.json")
    (root / "system.json").write_text(pack(system_data(system)))
    (root / "openmm-system.xml").write_text(mm.XmlSerializer.serialize(model))
    (root / "integrator.xml").write_text(mm.XmlSerializer.serialize(integrator))
    (root / "checkpoint.chk").write_bytes(context.createCheckpoint())
    (root / "result.json").write_text(
        json.dumps({"openmm_version": mm.__version__, "timestep_fs": 0.5, "steps": 2})
    )
    files = {p.name: sha256(p.read_bytes()).hexdigest() for p in root.iterdir()}
    (root / "manifest.json").write_text(json.dumps({"files": files}))
    pin = sha256((root / "manifest.json").read_bytes()).hexdigest()
    # Finish the uninterrupted reference before creating another stochastic
    # Reference context: backend random-state ownership is process dependent.
    integrator.step(2)
    context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    uninterrupted = context.createCheckpoint()
    result = resume_pdms_md(
        root,
        tmp_path / "resumed",
        options=prepared.options,
        manifest_sha256=pin,
        steps=2,
    )
    assert (tmp_path / "resumed/checkpoint.chk").read_bytes() == uninterrupted
    assert result["thermodynamics"][0]["step"] == 4
    assert result["thermodynamics"][0]["time_ps"] == 0.002
    assert result["thermodynamics"][0]["temperature_kelvin"] > 0
    (root / "checkpoint.chk").write_bytes(b"tampered")
    with pytest.raises(PDMSOPLSError, match="checkpoint_asset_hash_mismatch"):
        resume_pdms_md(
            root,
            tmp_path / "forbidden",
            options=prepared.options,
            manifest_sha256=pin,
            steps=2,
        )
    assert not (tmp_path / "forbidden").exists()
    del context


@pytest.mark.skipif(
    os.environ.get("ISLAND_PDMS_NUMERICAL") != "1",
    reason="Explicit slow physical finite-difference qualification",
)
@pytest.mark.parametrize("dp", [0, 2, 4])
def test_actual_pdms_strict_finite_differences(options, dp):
    pytest.importorskip("openmm")
    system, anchors = molecule(dp)
    evaluator = create_evaluator(system, prepare(system, options, anchors))
    try:
        result = finite_difference(evaluator, all_components=True)
        assert result["passed"], result
    finally:
        evaluator.close()
