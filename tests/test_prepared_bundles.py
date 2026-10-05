"""Data-only persistence regressions; synthetic native fixtures are software tests."""

import json
import shutil
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace

import pytest
from test_oplsaa_parameters import synthetic as opls_fixture  # noqa: F401
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.forcefields import (
    AmberBundleArtifacts,
    PreparedBundleError,
    PreparedForceFieldSources,
    adopt_forcefield,
    create_evaluator,
    load_prepared_forcefield,
    save_prepared_forcefield,
)
from island.workflows import storage
from island.workflows.bundle import system_data

# ruff: noqa: F811 -- pytest fixture injection


@pytest.fixture
def opls(opls_fixture, tmp_path):
    system, source, native = opls_fixture
    path = tmp_path / "source.xml"
    path.write_bytes(source.xml)
    return (
        system,
        adopt_forcefield(system, "oplsaa", native, source=source),
        PreparedForceFieldSources(opls_xml=path),
    )


def tree(path):
    return {
        str(p.relative_to(path)): p.read_bytes() for p in path.rglob("*") if p.is_file()
    }


def edit_manifest(root, change):
    envelope = storage.read_json(root / "manifest.json")
    change(envelope["payload"])
    envelope["sha256"] = storage.checksum(storage.json_bytes(envelope["payload"]))
    (root / "manifest.json").write_bytes(storage.json_bytes(envelope))


def edit_file(root, logical, transform):
    p = storage.read_json(root / "manifest.json")["payload"]
    path = root / p["files"][logical]["path"]
    path.write_bytes(transform(path.read_bytes()))
    edit_manifest(
        root,
        lambda p: p["files"][logical].update(
            sha256=storage.checksum(path.read_bytes())
        ),
    )


def test_roundtrip_relocation_ownership_and_coordinates(opls, tmp_path):
    system, prepared, sources = opls
    system.metadata["not_in_original"] = True
    with pytest.raises(PreparedBundleError, match="original"):
        save_prepared_forcefield(system, prepared, tmp_path / "bad", sources=sources)
    del system.metadata["not_in_original"]
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )
    before = tree(root)
    shutil.move(root, tmp_path / "relocated")
    root = tmp_path / "relocated"
    loaded = load_prepared_forcefield(root, sources=sources)
    assert tree(root) == before
    assert system_data(loaded.system) == system_data(system)
    assert loaded.prepared.identity == prepared.identity
    assert loaded.prepared.native_result.identity == prepared.native_result.identity
    changed = loaded.system
    changed.coordinates.set(39, (1.4, 0, 0))
    assert create_evaluator(
        changed, loaded.prepared
    ).evaluate().potential_energy == pytest.approx(8)
    changed.topology.sites[13].mass = 99
    with pytest.raises(Exception, match="binding"):
        create_evaluator(changed, loaded.prepared)
    assert loaded.system.topology.sites[13].mass == system.topology.sites[13].mass
    manifest = loaded.manifest
    manifest["payload"]["family"] = "pcff"
    assert loaded.manifest["payload"]["family"] == "oplsaa"
    assert tree(root) == before


@pytest.mark.parametrize(
    "kind",
    [
        "missing",
        "truncated",
        "checksum",
        "duplicate",
        "nan",
        "path",
        "symlink",
        "family",
        "schema",
        "metadata",
        "identity",
    ],
)
def test_bad_bundle(opls, tmp_path, kind):
    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )
    if kind == "missing":
        (root / "parameters.json").unlink()
    elif kind == "truncated":
        edit_file(root, "parameters", lambda _: b"{")
    elif kind == "checksum":
        (root / "parameters.json").write_text("{}")
    elif kind == "duplicate":
        edit_file(root, "parameters", lambda _: b'{"x":1,"x":2}')
    elif kind == "nan":
        edit_file(root, "system", lambda _: b'{"x":NaN}')
    elif kind == "path":
        edit_manifest(
            root, lambda p: p["files"]["system"].update(path="../system.json")
        )
    elif kind == "symlink":
        raw = (root / "system.json").read_bytes()
        (tmp_path / "outside").write_bytes(raw)
        (root / "system.json").unlink()
        (root / "system.json").symlink_to(tmp_path / "outside")
    elif kind == "metadata":
        edit_manifest(root, lambda p: p["prepared"].update(charge_method="am1bcc"))
    else:
        key = {"family": "family", "schema": "schema", "identity": "prepared_identity"}[
            kind
        ]
        edit_manifest(root, lambda p: p.update({key: "bad"}))
    with pytest.raises(PreparedBundleError):
        load_prepared_forcefield(root, sources=sources)


@pytest.mark.parametrize("change", ["mass", "stereo", "ids", "charge", "native"])
def test_semantic_tampering(opls, tmp_path, change):
    from island.charge_references.records import pack, unpack

    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )
    if change in ("mass", "stereo", "ids"):

        def modify(raw):
            p = storage.decode(json.loads(raw))
            if change == "mass":
                p["sites"][0]["mass"] = 9.0
            if change == "stereo":
                p["sites"][0]["metadata"]["cip_label"] = "R"
            if change == "ids":
                p["sites"][0]["id"] = 999
            return storage.json_bytes(storage.encode(p))

        edit_file(root, "system", modify)
    else:

        def modify(raw):
            p = unpack(raw.decode())
            if change == "charge":
                p["resolved"]["sites"][13]["charge_e"] = 0.3
            else:
                p["resolved"]["bonds"][0]["converted"]["k_kj_mol_angstrom2"] = 99
            return pack(p).encode()

        edit_file(root, "parameters", modify)
    with pytest.raises(PreparedBundleError):
        load_prepared_forcefield(root, sources=sources)


def test_sources_explicit_and_pinned(opls, tmp_path):
    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )
    for source in (
        None,
        PreparedForceFieldSources(pcff_frc=sources.opls_xml),
        PreparedForceFieldSources(opls_xml=tmp_path / "missing"),
    ):
        with pytest.raises(PreparedBundleError):
            load_prepared_forcefield(root, sources=source)
    sources.opls_xml.write_bytes(b"<wrong/>")
    with pytest.raises(PreparedBundleError, match="hash"):
        load_prepared_forcefield(root, sources=sources)


def test_overwrite_and_failed_publication(opls, tmp_path, monkeypatch):
    import island.forcefields.bundles as backend

    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )
    before = tree(root)
    with pytest.raises(PreparedBundleError, match="exists"):
        save_prepared_forcefield(system, prepared, root, sources=sources)
    assert tree(root) == before

    def fail(*args):
        raise OSError("injected rename failure")

    monkeypatch.setattr(backend.os, "replace", fail)
    with pytest.raises(PreparedBundleError, match="injected rename"):
        save_prepared_forcefield(system, prepared, tmp_path / "failed", sources=sources)
    assert not (tmp_path / "failed").exists()
    assert not list(tmp_path.glob(".failed.candidate-*"))
    assert tree(root) == before


def test_loading_runs_no_scientific_execution(opls, tmp_path, monkeypatch):
    import island.evaluation.oplsaa as evaluation
    import island.forcefields.oplsaa as native

    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )

    def forbidden(*args, **kwargs):
        pytest.fail("scientific execution during load")

    for name in ("parameterize_oplsaa", "type_atoms"):
        monkeypatch.setattr(native, name, forbidden)
    monkeypatch.setattr(evaluation, "OPLSSinglePointEvaluator", forbidden)
    monkeypatch.setattr(subprocess, "run", forbidden)
    for name in ("openmm", "foyer", "rdkit", "scipy", "parmed"):
        monkeypatch.setitem(sys.modules, name, None)
    loaded = load_prepared_forcefield(root, sources=sources)
    assert loaded.prepared.identity == prepared.identity


def test_import_is_optional():
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from island.forcefields import load_prepared_forcefield, save_prepared_forcefield; assert not {'openmm','foyer','parmed','rdkit','scipy'} & sys.modules.keys()",
        ],
        check=True,
    )


def amber_fixture(tmp_path, family):
    from test_ambertools_integrity import _valid_result

    from island.forcefields.ambertools.models import digest

    system, native = _valid_result(tmp_path)
    if family == "gaff":
        inner = deepcopy(
            dict(native.imported_result.provenance["ambertools_preparation"])
        )
        inner["requested_force_field"] = "gaff"
        for field in ("force_field_data", "leaprc"):
            inner[field]["path"] = inner[field]["path"].replace("gaff2", "gaff")
        cmd = inner["stages"][0]["command"]
        cmd[cmd.index("-at") + 1] = "gaff"
        cmd = inner["stages"][1]["command"]
        cmd[cmd.index("-s") + 1] = "1"
        imported = replace(
            native.imported_result,
            provenance={
                **dict(native.imported_result.provenance),
                "force_field": "gaff",
                "ambertools_preparation": inner,
            },
            result_signature="",
        )
        imported = replace(imported, result_signature=imported.content_signature())
        outer = {**inner, "imported_result_signature": imported.result_signature}
        native = replace(
            native,
            imported_result=imported,
            record=outer,
            record_signature=digest(outer),
        )
    prmtop = next(
        p
        for p in tmp_path.rglob("*.prmtop")
        if storage.checksum(p.read_bytes()) == native.imported_result.source_sha256
    )
    return (
        system,
        adopt_forcefield(system, family, native),
        AmberBundleArtifacts(prmtop),
    )


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
def test_amber_roundtrip_and_historical_relationship(tmp_path, monkeypatch, family):
    system, prepared, artifacts = amber_fixture(tmp_path, family)
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", artifacts=artifacts
    )
    artifacts.prmtop.unlink()  # only explicit bundled copy is operational
    from island.forcefields.ambertools import AmberToolsParameterizationEngine

    def forbidden(*args, **kwargs):
        pytest.fail("Amber execution on load")

    monkeypatch.setattr(AmberToolsParameterizationEngine, "parameterize", forbidden)
    monkeypatch.setattr(subprocess, "run", forbidden)
    loaded = load_prepared_forcefield(root)
    assert loaded.prepared.identity == prepared.identity
    assert system_data(loaded.system) == system_data(system)
    assert (
        create_evaluator(system, prepared).evaluate()
        == create_evaluator(loaded.system, loaded.prepared).evaluate()
    )
    moved = loaded.system
    moved.coordinates.set(min(moved.topology.sites), (0.01, 0, 0))
    create_evaluator(moved, loaded.prepared)
    with pytest.raises(PreparedBundleError, match="original"):
        save_prepared_forcefield(
            moved,
            loaded.prepared,
            tmp_path / "changed",
            artifacts=AmberBundleArtifacts(root / "result.prmtop"),
        )

    def modify(raw):
        p = storage.decode(json.loads(raw))
        p["record"]["charge_method"] = "am1bcc"
        return storage.json_bytes(storage.encode(p))

    edit_file(root, "preparation", modify)
    with pytest.raises(PreparedBundleError):
        load_prepared_forcefield(root)


def test_pcff_roundtrip_policy_and_no_repreparation(parameters, tmp_path, monkeypatch):
    import island.forcefields.pcff as backend
    from island import Coordinates
    from island.charge_references.records import pack, unpack
    from island.forcefields.pcff.automatic import graph_system

    system = graph_system(
        parameters.payload["charge_record"]["automatic_typing"]["graph"]
    )
    system.coordinates = Coordinates(
        {
            i: (float(n % 3), float(n // 3), 0.2 * (n % 2))
            for n, i in enumerate(sorted(system.topology.sites))
        }
    )
    native = backend.define_pcff_model(
        parameters,
        special_pairs=backend.special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    prepared = adopt_forcefield(system, "pcff", native)
    source = tmp_path / "source.frc"
    source.write_bytes(parameters.source.raw)
    sources = PreparedForceFieldSources(pcff_frc=source)
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )

    def forbidden(*args, **kwargs):
        pytest.fail("preparation during load")

    for name in (
        "type_pcff_atoms",
        "assign_automatic_pcff_charges",
        "assign_pcff_parameters",
        "define_pcff_model",
    ):
        monkeypatch.setattr(backend, name, forbidden)
    for name in ("openmm", "rdkit", "foyer", "scipy", "parmed"):
        monkeypatch.setitem(sys.modules, name, None)
    loaded = load_prepared_forcefield(root, sources=sources)
    assert loaded.prepared.identity == prepared.identity
    assert system_data(loaded.system) == system_data(system)

    def modify(raw):
        p = unpack(raw.decode())
        p["special_pairs"]["lj"][2] = 0.5
        return pack(p).encode()

    edit_file(root, "model", modify)
    with pytest.raises(PreparedBundleError):
        load_prepared_forcefield(root, sources=sources)


def test_typed_metadata_and_integer_inventory_roundtrip(opls, tmp_path):
    system, old, sources = opls
    system.metadata["coordinate_provenance_test"] = {
        17: ("stable", {"nested": [True, 3, 0.125]})
    }
    prepared = adopt_forcefield(system, "oplsaa", old.native_result, source=old.source)
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )
    loaded = load_prepared_forcefield(root, sources=sources)
    assert system_data(loaded.system) == system_data(system)
    assert type(loaded.system.metadata["coordinate_provenance_test"][17]) is tuple
    raw = loaded.prepared.native_result
    assert raw.identity == old.native_result.identity


@pytest.mark.parametrize("kind", ["missing", "wrong"])
def test_original_prmtop_required_before_publication(tmp_path, kind):
    system, prepared, artifacts = amber_fixture(tmp_path, "gaff2")
    if kind == "missing":
        artifacts.prmtop.unlink()
    else:
        artifacts.prmtop.write_bytes(b"wrong")
    with pytest.raises(PreparedBundleError):
        save_prepared_forcefield(
            system, prepared, tmp_path / "bundle", artifacts=artifacts
        )
    assert not (tmp_path / "bundle").exists()


def test_canonical_system_rejects_boolean_coordinates(opls, tmp_path):
    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )

    def modify(raw):
        p = storage.decode(json.loads(raw))
        p["coordinates"][13][0] = True
        return storage.json_bytes(storage.encode(p))

    edit_file(root, "system", modify)
    with pytest.raises(PreparedBundleError, match="canonical"):
        load_prepared_forcefield(root, sources=sources)


def test_concurrent_destination_wins_unchanged(opls, tmp_path, monkeypatch):
    from pathlib import Path

    system, prepared, sources = opls
    target = tmp_path / "bundle"
    original = Path.mkdir

    def race(path, *a, **kw):
        if path == target:
            original(path)
            (path / "existing-owner").write_text("do not change")
            raise FileExistsError("concurrent writer won")
        return original(path, *a, **kw)

    monkeypatch.setattr(Path, "mkdir", race)
    with pytest.raises(PreparedBundleError, match="concurrent writer"):
        save_prepared_forcefield(system, prepared, target, sources=sources)
    assert tree(target) == {"existing-owner": b"do not change"}
    assert not list(tmp_path.glob(".bundle.candidate-*"))


def test_cleanup_does_not_mask_original_failure(opls, tmp_path, monkeypatch):
    import island.forcefields.bundles as backend

    system, prepared, sources = opls

    def fail(*a, **kw):
        raise OSError("original publication failure")

    def cleanup(*a, **kw):
        raise OSError("secondary cleanup failure")

    monkeypatch.setattr(backend.os, "replace", fail)
    monkeypatch.setattr(backend.shutil, "rmtree", cleanup)
    with pytest.raises(PreparedBundleError, match="original publication") as error:
        save_prepared_forcefield(system, prepared, tmp_path / "bundle", sources=sources)
    assert error.value.__cause__ is not None
    assert not (tmp_path / "bundle").exists()


def test_amber_existing_nonperiodic_box_is_preserved(tmp_path):
    from island.core import SimulationBox

    system, old, artifacts = amber_fixture(tmp_path, "gaff2")
    system.box = SimulationBox(12.0, 13.0, 14.0, periodic=(False, False, False))
    prepared = adopt_forcefield(system, "gaff2", old.native_result)
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", artifacts=artifacts
    )
    loaded = load_prepared_forcefield(root)
    assert system_data(loaded.system) == system_data(system)
    assert loaded.prepared.identity == prepared.identity

    def modify(raw):
        p = storage.decode(json.loads(raw))
        p["box"]["periodic"] = (True, False, False)
        return storage.json_bytes(storage.encode(p))

    edit_file(root, "system", modify)
    with pytest.raises(PreparedBundleError, match="Periodic"):
        load_prepared_forcefield(root)


def test_source_can_be_resolved_at_a_new_explicit_location(opls, tmp_path):
    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "bundle", sources=sources
    )
    relocated = tmp_path / "moved-library.xml"
    sources.opls_xml.rename(relocated)
    loaded = load_prepared_forcefield(
        root, sources=PreparedForceFieldSources(opls_xml=relocated)
    )
    assert loaded.prepared.identity == prepared.identity
