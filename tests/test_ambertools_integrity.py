"""Offline parsed Amber fixtures for Phase 4D2.1 integrity and restart checks."""

import copy
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.exceptions import AmberToolsStageError, InvalidAmberImportResultError
from island.forcefields.amber import import_amber_prmtop
from island.forcefields.ambertools.engine import _prmtop_lineage
from island.forcefields.ambertools.lineage import generated_atom_name
from island.forcefields.ambertools.models import (
    LEGACY_PREPARATION_SCHEMA as PREPARATION_SCHEMA,
)
from island.forcefields.ambertools.models import (
    AmberToolsPreparationResult,
    digest,
)

pmd = pytest.importorskip("parmed")
PHENOL = Path(__file__).parent / "fixtures/phenol_parmed_13239c2.prmtop"


def _synthetic_source(tmp_path):
    """Write a real parsed prmtop/restart, not a mocked ParmEd return value."""
    structure = pmd.Structure()
    atom_type = pmd.AtomType("c3", 1, 12.01, 6)
    atom_type.set_lj_params(0.1, 1.9)
    atoms = []
    for index in range(4):
        atom = pmd.Atom(
            name=f"C{index:03}", type="c3", atomic_number=6,
            mass=12.01, charge=0.0,
        )
        atom.atom_type = atom_type
        structure.add_atom(atom, "SYN", 1)
        atoms.append(atom)
    bond_type = pmd.BondType(100, 1.5)
    angle_type = pmd.AngleType(50, 109.5)
    dihedral_type = pmd.DihedralType(1, 3, 0, scee=1.2, scnb=2)
    structure.bond_types.append(bond_type)
    structure.angle_types.append(angle_type)
    structure.dihedral_types.append(dihedral_type)
    for left, right in ((0, 1), (1, 2), (2, 3)):
        structure.bonds.append(pmd.Bond(atoms[left], atoms[right], type=bond_type))
    for left, middle, right in ((0, 1, 2), (1, 2, 3)):
        structure.angles.append(pmd.Angle(
            atoms[left], atoms[middle], atoms[right], type=angle_type,
        ))
    structure.dihedrals.append(pmd.Dihedral(*atoms, type=dihedral_type))
    structure.bond_types.claim()
    structure.angle_types.claim()
    structure.dihedral_types.claim()
    parm = pmd.amber.AmberParm.from_structure(structure)
    xyz = np.array(((0, 0, 0), (1.5, 0, 0), (2.4, 1.1, 0), (3.4, 1.3, 1.0)), float)
    parm.coordinates = xyz
    prmtop = tmp_path / "source.prmtop"
    restart = tmp_path / "source.rst7"
    parm.save(str(prmtop))
    parm.save(str(restart))
    graph = Topology()
    site_ids = (101, 109, 117, 125)
    for index, site_id in enumerate(site_ids):
        graph.add_site(AtomSite(site_id, f"C{index}", 12.01,
                               element="C", atomic_number=6))
    for left, right in ((0, 1), (1, 2), (2, 3)):
        graph.add_bond(site_ids[left], site_ids[right], order=1)
    system = MolecularSystem(graph, Coordinates({
        site_ids[index]: xyz[index] for index in range(4)
    }))
    mapping = {index: site_ids[index] for index in range(4)}
    names = {generated_atom_name("C", index): site_ids[index]
             for index in range(4)}
    return system, prmtop, restart, mapping, names


def _valid_result(tmp_path, method="provided"):
    system, prmtop, _restart, mapping, names = _synthetic_source(tmp_path)
    imported = import_amber_prmtop(system, prmtop, mapping, source="synthetic parsed Amber")
    xyz = {str(site_id): system.coordinates.get(site_id).tolist()
           for site_id in sorted(system.topology.sites)}
    fake_hash = sha256(b"independent synthetic fixture payload").hexdigest()
    record = {
        "schema": PREPARATION_SCHEMA,
        "engine_version": "2",
        "requested_force_field": "gaff2",
        "charge_method": "provided",
        "charge_outcome": {"mode": "provided", "qm_run": False,
                           "serialization_tolerance_e": 1e-5},
        "charge_validation_tolerance_e": 1e-4,
        "provided_charge_tolerance_e": 1e-5,
        "input_coordinate_signature": digest(xyz),
        "input_coordinates_angstrom": xyz,
        "antechamber_coordinates_angstrom": copy.deepcopy(xyz),
        "expected_cip_by_site": {},
        "input_mol2_sha256": fake_hash,
        "input_lineage_sha256": fake_hash,
        "lineage": {
            "input_name_to_site_id": names,
            "typed_mol2_index_to_site_id": {i + 1: site_id
                                                 for i, site_id in mapping.items()},
            "prmtop_index_to_site_id": mapping,
        },
        "artifact_sha256": {
            name: imported.source_sha256 if name == "result.prmtop" else fake_hash
            for name in ("typed.mol2", "typed.frcmod", "result.prmtop",
                         "result.rst7", "leap.in", "leap.log")
        },
        "force_field_data": {"path": "/synthetic/dat/leap/parm/gaff2.dat", "sha256": fake_hash,
                             "header": "synthetic test source"},
        "leaprc": {"path": "/synthetic/dat/leap/cmd/leaprc.gaff2", "sha256": fake_hash},
        "amberhome": "/synthetic",
        "executables": {name: f"/synthetic/bin/{name}" for name in (
            "antechamber", "parmchk2", "tleap",
        )},
        "executable_sha256": {name: fake_hash for name in (
            "antechamber", "parmchk2", "tleap",
        )},
        "tool_versions": {name: "synthetic fixture: not run" for name in (
            "antechamber", "parmchk2", "tleap",
        )},
        "stages": [
            {"stage": "antechamber", "command": [
                "/synthetic/bin/antechamber", "-i", "input.mol2", "-fi", "mol2",
                "-o", "typed.mol2", "-fo", "mol2", "-at", "gaff2", "-c", "rc",
                "-nc", "0", "-m", "1", "-s", "2", "-j", "4", "-du", "yes",
                "-pf", "no", "-cf", "charges.txt",
            ], "returncode": 0},
            {"stage": "parmchk2", "command": [
                "/synthetic/bin/parmchk2", "-i", "typed.mol2", "-f", "mol2",
                "-o", "typed.frcmod", "-s", "2",
            ], "returncode": 0},
            {"stage": "tleap", "command": [
                "/synthetic/bin/tleap", "-f", "leap.in",
            ], "returncode": 0},
        ],
    }
    if method == "am1bcc":
        record["charge_method"] = method
        record["charge_outcome"] = {
            "mode": method, "qm_run": True, "sqm_out_sha256": fake_hash,
            "convergence_marker": "Calculation Completed", "sqm_warning_lines": [],
        }
        command = record["stages"][0]["command"]
        command[command.index("-c") + 1] = "bcc"
        del command[-2:]  # provided-charge file option
    imported = replace(imported, provenance={
        **dict(imported.provenance), "force_field": "gaff2",
        "charge_method": "provided" if method == "provided" else "AM1-BCC",
        "ambertools_preparation": record,
    }, result_signature="")
    imported = replace(imported, result_signature=imported.content_signature())
    outer = {**record, "imported_result_signature": imported.result_signature}
    return system, AmberToolsPreparationResult(imported, outer, digest(outer))


def _resign_both(result, change):
    """Alter both copies and recompute both signatures to test semantics, not hashes."""
    inner = copy.deepcopy(dict(result.imported_result.provenance["ambertools_preparation"]))
    change(inner)
    imported = replace(result.imported_result, provenance={
        **dict(result.imported_result.provenance),
        "ambertools_preparation": inner,
    }, result_signature="")
    imported = replace(imported, result_signature=imported.content_signature())
    outer = {**inner, "imported_result_signature": imported.result_signature}
    return replace(result, imported_result=imported, record=outer,
                   record_signature=digest(outer))


def test_actual_imported_result_reconstruction_and_snapshot(tmp_path):
    system, result = _valid_result(tmp_path)
    result.validate_integrity(system)
    copy.deepcopy(result).validate_integrity(system)
    replace(result).validate_integrity(system)
    snapshot = result.to_parameterized_system(system)
    assert snapshot.metadata["aggregate"]["production_validated"] is False
    assert snapshot.metadata["aggregate"]["simulation_readiness"] == "not_established"
    assert snapshot.aggregate_signature == result.imported_result.result_signature
    system.coordinates.set(101, (9, 9, 9))
    assert snapshot.system.coordinates.get(101)[0] != 9


@pytest.mark.parametrize("change,match", [
    (lambda r: r.update(requested_force_field="gaff"), "force-field|contradicts"),
    (lambda r: r.update(charge_method="am1bcc"), "charge-method|contradicts"),
    (lambda r: r["charge_outcome"].update(qm_run=True), "charge outcome|contradicts"),
    (lambda r: r["input_coordinates_angstrom"].pop("101"), "coordinate|contradicts"),
    (lambda r: r["lineage"]["prmtop_index_to_site_id"].update({0: 999}),
     "mapping|contradicts"),
    (lambda r: r["artifact_sha256"].update({"result.prmtop": "0" * 64}),
     "checksum|contradicts"),
    (lambda r: r["stages"].pop(), "stage|contradicts"),
])
def test_rehashed_outer_record_cannot_contradict_imported_provenance(tmp_path, change, match):
    system, result = _valid_result(tmp_path)
    original = system.to_dict()
    record = copy.deepcopy(dict(result.record))
    change(record)
    malformed = replace(result, record=record, record_signature=digest(record))
    with pytest.raises(InvalidAmberImportResultError, match=match):
        malformed.validate_integrity(system)
    with pytest.raises(InvalidAmberImportResultError):
        malformed.to_parameterized_system(system)
    assert system.to_dict() == original
    result.validate_integrity(system)


def test_missing_imported_result_raises_domain_error(tmp_path):
    system, result = _valid_result(tmp_path)
    malformed = replace(result, imported_result=None)
    with pytest.raises(InvalidAmberImportResultError, match="ImportedAmberResult"):
        malformed.validate_integrity(system)
    with pytest.raises(InvalidAmberImportResultError):
        malformed.to_parameterized_system(system)


@pytest.mark.parametrize("change,match", [
    (lambda r: r.update(requested_force_field="gaff"), "force-field"),
    (lambda r: r.update(charge_method="am1bcc"), "charge-method"),
    (lambda r: r["stages"].pop(), "stage"),
    (lambda r: r["stages"][0].update(command=[]), "command is missing or empty"),
    (lambda r: r["stages"][2].update(command=[]), "command is missing"),
    (lambda r: r.update(engine_version="1"), "engine version"),
    (lambda r: r["antechamber_coordinates_angstrom"].update({
        "101": [float("nan"), 0.0, 0.0],
    }), "antechamber coordinate"),
    (lambda r: r["charge_outcome"].update(qm_run=True), "charge outcome"),
    (lambda r: r["lineage"]["prmtop_index_to_site_id"].update({0: 999}),
     "prmtop_index_to_site_id"),
    (lambda r: r["artifact_sha256"].update({"result.prmtop": "0" * 64}),
     "prmtop artifact"),
    (lambda r: r.update(requested_force_field=[]), "force-field"),
    (lambda r: r.update(charge_method={}), "charge-method"),
    (lambda r: r["expected_cip_by_site"].update({101: []}), "CIP"),
    (lambda r: r["stages"][2]["command"].append({}), "nontext"),
    (lambda r: r["stages"][0]["command"].extend(["-at", "gaff"]), "duplicate"),
    (lambda r: r["stages"][1]["command"].extend(["-s", "1"]), "duplicate"),
    (lambda r: r["charge_outcome"].update(serialization_tolerance_e=100),
     "serialization tolerance"),
    (lambda r: r["charge_outcome"].pop("serialization_tolerance_e"),
     "serialization tolerance"),
    (lambda r: r.update(provided_charge_tolerance_e=100), "serialization tolerance"),
    (lambda r: r["charge_outcome"].update(convergence_marker="invented"), "SQM"),
])
def test_resigned_both_copies_still_reject_semantic_contradictions(
    tmp_path, change, match,
):
    system, result = _valid_result(tmp_path)
    original = system.to_dict()
    malformed = _resign_both(result, change)
    with pytest.raises(InvalidAmberImportResultError, match=match):
        malformed.validate_integrity(system)
    with pytest.raises(InvalidAmberImportResultError):
        malformed.to_parameterized_system(system)
    assert system.to_dict() == original


def test_preparation_record_owns_caller_data(tmp_path):
    system, result = _valid_result(tmp_path)
    caller_record = copy.deepcopy(dict(result.record))
    reconstructed = replace(result, record=caller_record)
    caller_record["lineage"]["prmtop_index_to_site_id"][0] = 999
    reconstructed.validate_integrity(system)


@pytest.mark.parametrize("bad_value", [float("nan"), float("inf")])
def test_real_parsed_restart_rejects_late_nonfinite_coordinate(tmp_path, bad_value):
    source = pmd.load_file(str(PHENOL))
    coordinates = np.arange(3 * len(source.atoms), dtype=float).reshape(-1, 3) / 10
    finite_coordinates = coordinates.copy()
    coordinates[-1, 1] = bad_value  # well after the finite format-identification records
    source.coordinates = coordinates
    bad_restart = tmp_path / "nonfinite.rst7"
    source.save(str(bad_restart))
    parsed = pmd.load_file(str(PHENOL), str(bad_restart))
    assert not np.isfinite(parsed.coordinates[-1, 1])
    graph = Topology()
    names = {}
    for atom in source.atoms:
        site_id = 201 + 7 * atom.idx
        symbol = pmd.periodic_table.Element[atom.atomic_number]
        graph.add_site(AtomSite(site_id, atom.name, atom.mass,
                               element=symbol, atomic_number=atom.atomic_number))
        names[atom.name] = site_id
    for bond in source.bonds:
        graph.add_bond(names[bond.atom1.name], names[bond.atom2.name], order=1)
    system = MolecularSystem(graph, Coordinates({
        201 + 7 * atom.idx: finite_coordinates[atom.idx]
        for atom in source.atoms
    }))
    original = system.to_dict()
    with pytest.raises(AmberToolsStageError, match="nonfinite_source_indices") as error:
        _prmtop_lineage(PHENOL, bad_restart, system, names, {}, tmp_path)
    assert error.value.stage == "tleap"
    assert system.to_dict() == original


def test_real_parsed_restart_passes_when_finite(tmp_path):
    system, prmtop, restart, mapping, names = _synthetic_source(tmp_path)
    observed_mapping, charges = _prmtop_lineage(
        prmtop, restart, system, names, {}, tmp_path,
    )
    assert observed_mapping == mapping
    assert charges == {site_id: 0 for site_id in mapping.values()}


def test_malformed_restart_shape_is_stage_error(tmp_path, monkeypatch):
    system, prmtop, restart, _mapping, names = _synthetic_source(tmp_path)
    parsed = pmd.load_file(str(prmtop), str(restart))
    monkeypatch.setattr(pmd, "load_file", lambda *_args: SimpleNamespace(
        atoms=parsed.atoms, coordinates=np.zeros((len(parsed.atoms), 2)),
    ))
    with pytest.raises(AmberToolsStageError, match="shape=") as error:
        _prmtop_lineage(prmtop, restart, system, names, {}, tmp_path)
    assert error.value.stage == "tleap"


@pytest.mark.parametrize("marker", [None, "", "invented", "not Calculation Completed"])
def test_am1bcc_record_requires_actual_completion_marker(tmp_path, marker):
    system, result = _valid_result(tmp_path, "am1bcc")
    result.to_parameterized_system(system)
    malformed = _resign_both(result, lambda r: r["charge_outcome"].update(
        convergence_marker=marker,
    ))
    for validate in (malformed.validate_integrity, malformed.to_parameterized_system):
        with pytest.raises(InvalidAmberImportResultError, match="completion marker"):
            validate(system)


@pytest.mark.parametrize("field", [
    "requested_force_field", "charge_method", "charge_outcome", "stages",
    "input_coordinates_angstrom", "input_coordinate_signature", "lineage",
    "input_mol2_sha256", "input_lineage_sha256", "artifact_sha256",
    "force_field_data", "leaprc", "amberhome", "executables", "executable_sha256",
    "tool_versions", "antechamber_coordinates_angstrom", "expected_cip_by_site",
])
def test_resigned_incomplete_preparation_is_rejected(tmp_path, field):
    system, result = _valid_result(tmp_path)
    malformed = _resign_both(result, lambda record: record.pop(field))
    for validate in (malformed.validate_integrity, malformed.to_parameterized_system):
        with pytest.raises(InvalidAmberImportResultError):
            validate(system)


@pytest.mark.parametrize("stage,flag,value", [
    (0, "-i", None), (0, "-nc", "2"), (0, "-m", "3"),
    (1, "-o", "unrelated.frcmod"),
])
def test_stage_commands_bind_inventory_and_artifacts(tmp_path, stage, flag, value):
    system, result = _valid_result(tmp_path)

    def change(record):
        command = record["stages"][stage]["command"]
        index = command.index(flag)
        if value is None:
            del command[index:index + 2]
        else:
            command[index + 1] = value

    malformed = _resign_both(result, change)
    for validate in (malformed.validate_integrity, malformed.to_parameterized_system):
        with pytest.raises(InvalidAmberImportResultError, match="command contradicts"):
            validate(system)
