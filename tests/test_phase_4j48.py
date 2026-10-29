"""Canonical OPLS selection uses exact chemistry and never grants blocked MD."""

import json
import os
import shutil
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from island.chemistry import from_smiles
from island.forcefields.opls_registry import (
    ForceFieldRequest,
    OPLSAssets,
    OPLSResolutionError,
    Preparation,
    canonical_registry,
    create_evaluator,
    prepare_forcefield,
    prepare_polymer,
    prepare_polymer_forcefield,
    use_opls_assets,
)
from island.forcefields.opls_registry.api import graph_binding, rank_profiles
from island.forcefields.opls_registry.registry import ARCHIVE_SHA256, NATIVE
from island.forcefields.pdms_opls import PROFILE, export_lammps
from island.forcefields.selection.inventory import canonical, digest
from island.opls_lammps import LAMMPSOptions, reconstruct, run_md, run_static
from scripts.qualify_pdms_public import molecule

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "island_ffref/pcff/PDMS_benchmark-main.zip"
LMP = ROOT.parent / "island-validation/phase4j22/lmp-kspace"
LMP_PIN = "6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191"
NATIVE_XML = (
    ROOT.parent
    / "island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml"
)


def fixture(dp=2):
    system, _ = molecule(dp)
    xyz = np.array([system.coordinates.get(i) for i in system.topology.sites])
    system.coordinates.translate(20 - xyz.mean(axis=0))
    return system


@pytest.fixture(scope="module")
def assets():
    if not ARCHIVE.is_file():
        pytest.skip("User-supplied local research ZIP unavailable")
    assert sha256(ARCHIVE.read_bytes()).hexdigest() == ARCHIVE_SHA256
    return OPLSAssets(ARCHIVE, ARCHIVE_SHA256)


def test_registry_owned_and_versioned():
    registry = canonical_registry()
    assert registry == json.loads(
        (ROOT / "docs/evidence/phase_4j48/canonical-registry.json").read_text()
    )
    assert registry["schema"] == "island_opls_registry_v1"
    assert {r["profile_identity"] for r in registry["entries"]} == {PROFILE, NATIVE}
    for row in registry["entries"]:
        assert row["family"] == "opls"
        assert row["source_hashes"]
        assert row["charge_origin"] and row["charge_method"]
        assert row["scientific_validation_status"] == "not_established"
    registry["entries"][0]["priority"] = 0
    assert canonical_registry() != registry


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_only_public_family_needed(assets, dp):
    system = fixture(dp)
    with use_opls_assets(assets):
        result = prepare_polymer_forcefield(system, forcefield="opls")
    p = result.payload
    assert p["status"] == "prepared_static_only", p
    assert p["selected_profile"] == PROFILE
    assert p["source_kind"] == "external_literature"
    assert p["charge_method"] == "author_dft_npa"
    assert p["charge_origin"] == "Figure 7 DFT/NPA; exact author serialized literals"
    assert not p["graph_charges"]["normalization"]
    assert not p["graph_charges"]["redistribution"]
    assert p["stages"]["typing_complete"] and p["stages"]["charge_complete"]
    assert p["stages"]["parameter_complete"] and not p["missing_terms"]
    assert p["stages"]["executable"] and not p["stages"]["lammps_MD_verified"]
    assert not p["stages"]["numerically_verified"]
    assert (
        not p["production_validated"] and p["simulation_readiness"] == "not_established"
    )
    backend = result.validate_integrity(system).validate_integrity(system).payload
    assert p["graph_charges"]["vector"] == [
        [r["site_id"], r["charge_e"]] for r in backend["mapping"]
    ]
    # Automatic selection preserves the original source-bound assignment verbatim.
    from island.forcefields.pdms_lammps import prepare_forcefield as old_prepare
    from tests.test_phase_4j47 import fixture as old_fixture

    old_system, old_request = old_fixture(dp)
    assert canonical(backend) == canonical(
        old_prepare(old_system, old_request).validate_integrity(old_system).payload
    )


def test_public_psmiles_constructor(assets):
    from island.core import SimulationBox

    with use_opls_assets(assets):
        result = prepare_polymer(
            psmiles="[*:1]O[Si](C)(C)[*:2]",
            forcefield="opls",
            dp=2,
            head_end_group="[*][Si](C)(C)C",
            tail_end_group="[*]C",
            box=SimulationBox(40, 40, 40),
        )
    assert result.payload["selected_profile"] == PROFILE
    assert result.validate_integrity().payload["status"] == "prepared_static_only"


def test_final_graph_and_optional_records(assets):
    system = fixture()
    p = prepare_forcefield(system, assets=assets).payload
    request = ForceFieldRequest(
        family="opls",
        final_graph=p["final_graph"],
        typed_graph=p["typed_graph"],
        graph_charges=p["graph_charges"],
    )
    result = prepare_forcefield(system, request, assets=assets)
    assert result.payload == p
    result.validate_integrity(system)


@pytest.mark.parametrize("record", ["typed_graph", "graph_charges", "final_graph"])
def test_exact_graph_record_mismatch(assets, record):
    system = fixture()
    p = prepare_forcefield(system, assets=assets).payload
    value = deepcopy(p[record])
    value["graph_identity"] = "wrong"
    result = prepare_forcefield(
        system, ForceFieldRequest(**{record: value}), assets=assets
    )
    assert result.payload["status"] in ("source_incompatible", "charge_incomplete")
    with pytest.raises(OPLSResolutionError):
        create_evaluator(system, result)


@pytest.mark.parametrize("mutation", ["profile", "charge_method", "family"])
def test_expert_overrides_are_constraints(assets, mutation):
    values = {"profile": "pcff", "charge_method": "qeq", "family": "pcff"}
    p = prepare_forcefield(
        fixture(), ForceFieldRequest(**{mutation: values[mutation]}), assets=assets
    ).payload
    assert p["status"] in ("source_incompatible", "charge_incomplete")
    assert not p["stages"]["executable"]


@pytest.mark.parametrize(
    "field",
    ["charge_method", "source_hashes", "lj_mixing", "coulomb_14", "profile_identity"],
)
def test_equal_priority_conflict_and_no_execution(assets, field):
    registry = canonical_registry()
    duplicate = deepcopy(
        next(e for e in registry["entries"] if e["profile_identity"] == PROFILE)
    )
    duplicate[field] = "conflict"
    registry["entries"].append(duplicate)
    system = fixture()
    p = prepare_forcefield(system, assets=assets, registry=registry)
    assert p.payload["status"] == "profile_conflict", p.payload
    assert p._backend is None
    assert (
        p.payload["diagnostics"][0]["stage"] == "equal_priority_opls_profiles_disagree"
    )
    with pytest.raises(OPLSResolutionError):
        p.validate_integrity(system)


def test_rank_order_independent_and_specificity_first():
    rows = canonical_registry()["entries"]
    assert rank_profiles(rows) == rank_profiles(rows[::-1])
    assert rank_profiles(rows)["profile_identity"] == PROFILE


@pytest.mark.parametrize("mutation", ["source", "priority", "profile"])
def test_registry_tampering(assets, mutation):
    registry = canonical_registry()
    row = next(e for e in registry["entries"] if e["profile_identity"] == PROFILE)
    row[
        {
            "source": "source_hashes",
            "priority": "priority",
            "profile": "profile_identity",
        }[mutation]
    ] = "tampered"
    assert (
        prepare_forcefield(fixture(), assets=assets, registry=registry).payload[
            "status"
        ]
        == "source_incompatible"
    )


@pytest.mark.parametrize(
    "kind", ["absent", "no_hash", "wrong_hash", "malformed", "mutated_source"]
)
def test_asset_failures(assets, tmp_path, kind):
    source = tmp_path / "source.zip"
    options = assets
    if kind == "absent":
        options = OPLSAssets()
    elif kind == "no_hash":
        options = replace(assets, archive_sha256=None)
    elif kind == "wrong_hash":
        options = replace(assets, archive_sha256="0" * 64)
    else:
        source.write_bytes(
            b"not a ZIP" if kind == "malformed" else ARCHIVE.read_bytes() + b"mutation"
        )
        options = replace(assets, archive_path=source)
    p = prepare_forcefield(fixture(), assets=options)
    assert p.payload["status"] in ("asset_unavailable", "source_incompatible")
    assert p.payload["attempted_profiles"] and p.payload["missing_charges"]
    assert p._backend is None


@pytest.mark.parametrize(
    "smiles", ["[SiH3]O[SiH3]", "CO[Si](C)(C)C", "CP(=O)(O)O", "C[SiH3]"]
)
def test_unsupported_domain_no_family_fallback(assets, smiles):
    system = from_smiles(smiles, random_seed=2026)
    p = prepare_forcefield(system, assets=assets)
    assert p.payload["status"] in ("typing_blocked", "typing_ambiguous")
    assert not p.payload["stages"]["executable"]
    assert p.payload["missing_terms"]
    assert p.payload["selected_profile"] is None
    with pytest.raises(OPLSResolutionError):
        create_evaluator(system, p)


def test_unanchored_symmetric_chain_is_ambiguous(assets):
    system = fixture()
    system.metadata.pop("polymer")
    p = prepare_forcefield(system, assets=assets)
    assert p.payload["status"] == "typing_ambiguous"
    assert "explicit_E_oxygen_required" in canonical(p.payload["diagnostics"])


def test_missing_interaction_exact_diagnostics(assets, monkeypatch):
    from island.forcefields import pdms_opls

    original = pdms_opls.load_source

    def missing(options):
        source = original(options)
        source["rows"] = [r for r in source["rows"] if r["family"] != "angle"]
        return source

    monkeypatch.setattr(pdms_opls, "load_source", missing)
    result = prepare_forcefield(fixture(), assets=assets)
    assert result.payload["status"] == "parameter_incomplete"
    assert "missing_interaction_rows" in canonical(result.payload["diagnostics"])
    assert any(r["family"] == "angle" for r in result.payload["missing_terms"])
    assert result._backend is None


def test_missing_source_charge_is_not_created(assets, monkeypatch):
    from island.forcefields import pdms_opls

    original = pdms_opls.load_source

    def missing(options):
        source = original(options)
        source["charge_literals"]["Si"] = []
        return source

    monkeypatch.setattr(pdms_opls, "load_source", missing)
    result = prepare_forcefield(fixture(), assets=assets)
    assert result.payload["status"] in ("parameter_incomplete", "charge_incomplete")
    assert "source_charge_literal_missing" in canonical(result.payload["diagnostics"])
    assert result.payload["missing_charges"]
    assert result._backend is None


@pytest.mark.parametrize("mutation", ["charge", "profile", "source", "status", "units"])
def test_resigned_assignment_mutations(assets, mutation):
    system = fixture()
    result = prepare_forcefield(system, assets=assets)
    p = result.payload
    if mutation == "charge":
        p["graph_charges"]["vector"][0][1] = "0"
    elif mutation == "units":
        p["graph_charges"]["unit"] = "coulomb"
    elif mutation == "status":
        p["stages"]["lammps_MD_verified"] = True
    else:
        p[{"profile": "selected_profile", "source": "source_identity"}[mutation]] = (
            "forged"
        )
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    forged = Preparation(canonical(p), result.request, assets)
    with pytest.raises(OPLSResolutionError):
        create_evaluator(system, forged)


@pytest.mark.parametrize("mutation", ["coordinate", "box", "order", "bond"])
def test_stale_graph_execution_binding(assets, mutation):
    from island.core import SimulationBox

    system = fixture()
    result = prepare_forcefield(system, assets=assets)
    graph_id = result.payload["graph_identity"]
    if mutation == "coordinate":
        system.coordinates.translate((0.01, 0, 0))
    elif mutation == "box":
        system.box = SimulationBox(41, 40, 40)
    elif mutation == "order":
        sites = system.topology.sites
        atom = sites.pop(next(iter(sites)))
        sites[atom.id] = atom
    else:
        system.topology.bonds.pop(next(iter(system.topology.bonds)))
    if mutation in ("coordinate", "box"):
        assert digest(graph_binding(system)) == graph_id
        assert (
            prepare_forcefield(system, assets=assets).payload["selected_profile"]
            == PROFILE
        )
    with pytest.raises(OPLSResolutionError):
        result.validate_integrity(system)


def test_name_and_coordinate_do_not_select_profile(assets):
    system = fixture()
    first = prepare_forcefield(system, assets=assets).payload
    system.metadata["name"] = "this is not PDMS"
    system.coordinates.translate((0.001, 0.003, 0.007))
    second = prepare_forcefield(system, assets=assets).payload
    assert first["selected_profile"] == second["selected_profile"]
    assert first["graph_identity"] == second["graph_identity"]
    assert first["graph_charges"] == second["graph_charges"]


@pytest.mark.parametrize("smiles", ["CC", "c1ccccc1", "COC"])
def test_native_opls_charge_resolution(assets, smiles):
    if not NATIVE_XML.is_file():
        pytest.skip("Pinned native XML unavailable")
    options = replace(
        assets,
        native_source_path=NATIVE_XML,
        native_source_sha256=sha256(NATIVE_XML.read_bytes()).hexdigest(),
    )
    p = prepare_forcefield(
        from_smiles(smiles, random_seed=2026), assets=options
    ).payload
    assert p["selected_profile"] == NATIVE
    assert p["charge_method"] == "native"
    assert p["stages"]["charge_complete"]
    assert p["stages"]["parameter_complete"], p["diagnostics"]


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_deterministic_lammps_export_and_md_block(assets, tmp_path, dp):
    system = fixture(dp)
    result = prepare_forcefield(system, assets=assets)
    backend = result.validate_integrity(system).validate_integrity(system)
    first = export_lammps(system, backend, tmp_path / "first")
    assert first == export_lammps(system, backend, tmp_path / "second")
    with pytest.raises(OPLSResolutionError) as error:
        run_md(system, result, tmp_path / "must-not-exist")
    assert error.value.status == "scientifically_unqualified"
    assert not (tmp_path / "must-not-exist").exists()
    assert not list(tmp_path.rglob("*.chk"))


@pytest.fixture(scope="module")
def workflow(assets, tmp_path_factory):
    if not LMP.is_file():
        pytest.skip("Pinned local LAMMPS unavailable")
    system = fixture()
    result = prepare_forcefield(system, assets=assets)
    root = tmp_path_factory.mktemp("j48") / "workflow"
    engine = LAMMPSOptions(LMP, LMP_PIN)
    observed = run_static(system, result, root, engine=engine)
    assert observed["status"] == "static_executed_unqualified"
    return system, result, root, engine, observed


def test_static_relocation_and_fresh_process(assets, workflow, tmp_path):
    system, result, root, engine, observed = workflow
    moved = tmp_path / "moved"
    shutil.copytree(root, moved)
    source = tmp_path / "renamed.zip"
    shutil.copyfile(ARCHIVE, source)
    moved_assets = replace(assets, archive_path=source)
    restored, prepared, _ = reconstruct(
        moved,
        manifest_sha256=observed["manifest_sha256"],
        assets=moved_assets,
        engine=engine,
    )
    assert graph_binding(system) == graph_binding(restored)
    assert prepared.json_text == result.json_text
    code = """import sys
from island.opls_lammps import *
from island.forcefields.opls_registry import OPLSAssets
s,p,r=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],assets=OPLSAssets(sys.argv[3],sys.argv[4]),engine=LAMMPSOptions(sys.argv[5],sys.argv[6]))
print(p.identity)
"""
    for seed in ("1", "98317"):
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                code,
                str(moved),
                observed["manifest_sha256"],
                str(source),
                ARCHIVE_SHA256,
                str(LMP),
                LMP_PIN,
            ],
            cwd=tmp_path,
            env={**os.environ, "PYTHONPATH": str(ROOT / "src"), "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
        )
        assert proc.stdout.strip() == result.identity


@pytest.mark.parametrize(
    "kind", ["asset_unavailable", "typing_ambiguous", "source_incompatible"]
)
def test_blocked_fresh_process(assets, tmp_path, kind):
    source = tmp_path / "renamed.zip"
    shutil.copyfile(ARCHIVE, source)
    code = """import sys
from island.forcefields.opls_registry import *
from island.chemistry import from_smiles
from island.core import SimulationBox
s=from_smiles("C[Si](C)(C)O[Si](C)(C)C" if sys.argv[3]!="typing_ambiguous" else "CP(=O)(O)O", random_seed=2026)
s.box=SimulationBox(40,40,40)
a=OPLSAssets() if sys.argv[3]=="asset_unavailable" else OPLSAssets(sys.argv[1],sys.argv[2] if sys.argv[3]!="source_incompatible" else "0"*64)
p=prepare_forcefield(s,assets=a)
try: create_evaluator(s,p)
except OPLSResolutionError as e:
 assert e.status==p.payload["status"]
 print(e.status)
"""
    proc = subprocess.run(
        [sys.executable, "-c", code, str(source), ARCHIVE_SHA256, kind],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.stdout.strip() == kind


@pytest.mark.parametrize("flag", [("--forcefield", "opls"), ("--forcefield=opls",)])
def test_cli_no_profile_or_charge_method(assets, tmp_path, flag):
    config = tmp_path / "assets.json"
    config.write_text(
        json.dumps({"archive_path": str(ARCHIVE), "archive_sha256": ARCHIVE_SHA256})
    )
    output = tmp_path / "report.json"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "island.inspect_forcefield",
            "--psmiles",
            "[*:1]O[Si](C)(C)[*:2]",
            "--dp",
            "2",
            "--head-end-group",
            "[*][Si](C)(C)C",
            "--tail-end-group",
            "[*]C",
            *flag,
            "--json",
            str(output),
        ],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHONPATH": str(ROOT / "src"),
            "ISLAND_OPLS_ASSETS": str(config),
            "ISLAND_OPLS_ASSETS_SHA256": sha256(config.read_bytes()).hexdigest(),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    p = json.loads(output.read_text())
    assert p["selected_profile"] == PROFILE and p["charge_method"] == "author_dft_npa"
    assert p["source_kind"] == "external_literature"


def test_registry_file_order_does_not_resolve_differently(assets):
    system = fixture()
    registry = canonical_registry()
    registry["entries"].reverse()
    assert (
        prepare_forcefield(system, assets=assets, registry=registry).json_text
        == prepare_forcefield(system, assets=assets).json_text
    )


def test_incomplete_supplied_charge_vector_rejected(assets):
    system = fixture()
    p = prepare_forcefield(system, assets=assets).payload
    charges = deepcopy(p["graph_charges"])
    charges["vector"].pop()
    result = prepare_forcefield(
        system, ForceFieldRequest(graph_charges=charges), assets=assets
    )
    assert result.payload["status"] == "charge_incomplete"
    assert not result.payload["stages"]["charge_complete"]
    assert result.payload["missing_charges"]


def test_source_mutation_after_preparation(assets, tmp_path):
    path = tmp_path / "source.zip"
    shutil.copyfile(ARCHIVE, path)
    local = replace(assets, archive_path=path)
    system = fixture()
    result = prepare_forcefield(system, assets=local)
    path.write_bytes(path.read_bytes() + b"mutation")
    with pytest.raises(OPLSResolutionError):
        create_evaluator(system, result)


def test_builder_partition_mutation_is_graph_mutation(assets):
    system = fixture()
    result = prepare_forcefield(system, assets=assets)
    system.metadata["polymer"]["head_site_id"] = system.metadata["polymer"][
        "tail_site_id"
    ]
    assert digest(graph_binding(system)) != result.payload["graph_identity"]
    with pytest.raises(OPLSResolutionError):
        result.validate_integrity(system)


def test_no_engine_before_blocked_exception(assets, monkeypatch, tmp_path):
    from island import pdms, pdms_lammps

    def forbidden(*a, **k):
        raise AssertionError("Blocked preparation constructed a numerical engine")

    monkeypatch.setattr(pdms, "_mm", forbidden)
    monkeypatch.setattr(pdms_lammps, "_execute", forbidden)
    system = fixture()
    blocked = prepare_forcefield(system, assets=OPLSAssets())
    with pytest.raises(OPLSResolutionError) as error:
        create_evaluator(system, blocked)
    assert error.value.status == "asset_unavailable"
    assert error.value.result is blocked
    prepared = prepare_forcefield(system, assets=assets)
    with pytest.raises(OPLSResolutionError):
        run_md(system, prepared, tmp_path / "no-MD")
    assert not (tmp_path / "no-MD").exists()


def test_automatic_openmm_evaluator_keeps_source_assignment(assets):
    system = fixture(0)
    result = prepare_forcefield(system, assets=assets)
    evaluator = create_evaluator(system, result)
    value = evaluator.evaluate()
    assert np.isfinite(value.potential_energy)
    assert not result.payload["stages"]["scientifically_validated"]


@pytest.mark.parametrize(
    "mutation", ["canonical_charge", "cutoff", "mixing", "scaling", "reciprocal", "box"]
)
def test_resigned_workflow_tampering(assets, workflow, tmp_path, mutation):
    from island.charge_references.records import pack, unpack
    from island.forcefields.pdms_opls import PDMSOPLSError
    from island.pdms_lammps import seal

    _, _, root, engine, _ = workflow
    moved = tmp_path / mutation
    shutil.copytree(root, moved)
    if mutation == "canonical_charge":
        value = json.loads((moved / "canonical-assignment.json").read_text())
        value["graph_charges"]["vector"][0][1] = "0"
        value["identity"] = digest({k: v for k, v in value.items() if k != "identity"})
        (moved / "canonical-assignment.json").write_text(canonical(value))
    elif mutation == "box":
        value = unpack((moved / "system.json").read_text())
        value["box"]["lengths"] = (41, 40, 40)
        (moved / "system.json").write_text(pack(value))
    else:
        path = moved / "lammps/system.in"
        old, new = {
            "cutoff": ("9 11", "9 12"),
            "mixing": ("mix geometric", "mix arithmetic"),
            "scaling": ("lj/coul 0 0 0.5", "lj/coul 0 0 1"),
            "reciprocal": ("pppm 0.0001", "pppm 0.00001"),
        }[mutation]
        assert old in path.read_text()
        path.write_text(path.read_text().replace(old, new))
    with pytest.raises((OPLSResolutionError, PDMSOPLSError)):
        reconstruct(moved, manifest_sha256=seal(moved), assets=assets, engine=engine)


def test_receipt_validator_rejects_independent_pin_mutation(tmp_path):
    from scripts.validate_phase_4j48_receipt import validate

    (tmp_path / "receipt.json").write_text("{}")
    with pytest.raises(ValueError, match="Independent receipt hash mismatch"):
        validate(tmp_path, repo=ROOT, expected_sha256="0" * 64)


def test_receipt_validator_rejects_resigned_schema(tmp_path):
    from scripts.validate_phase_4j48_receipt import validate

    for name in ("receipt", "manifest"):
        p = {"schema": "forged"}
        (tmp_path / (name + ".json")).write_text(
            canonical({**p, "identity": digest(p)})
        )
    with pytest.raises(ValueError, match="Wrong schema"):
        validate(tmp_path, repo=ROOT)
