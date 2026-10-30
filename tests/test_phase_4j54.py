"""Separate Class-II source domains, exact missing couplings and J55 claims."""

import json
import os
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from island.chemistry import from_smiles
from island.forcefield_policy import ClaimRequirements
from island.forcefield_policy.claims import assess
from island.forcefields.selection.inventory import canonical, digest
from island.pcff_pdms import (
    CANDIDATE,
    REGISTERED_PROFILES,
    Class2AuditError,
    PCFFPDMSAssets,
    build_openmm,
    export_lammps,
    prepare_polymer,
    prepare_polymer_forcefield,
    reconstruct,
    save,
    use_pcff_pdms_assets,
)
from scripts.audit_pcff_pdms_class2_sources import archive_audit, independent_frc
from scripts.phase_4j54 import fixtures, repeat_audit

ROOT = Path(__file__).resolve().parents[1]
DISTRIBUTION = (
    ROOT.parent / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
)
ARCHIVE = ROOT / "island_ffref/pcff/PDMS_benchmark-main.zip"


@pytest.fixture(scope="module")
def assets():
    if not (DISTRIBUTION / "pcff.frc").is_file():
        pytest.skip("pinned local PCFF distribution unavailable")
    return PCFFPDMSAssets(DISTRIBUTION)


@pytest.fixture(scope="module")
def hmds(assets):
    s = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
    with use_pcff_pdms_assets(assets):
        p = prepare_polymer_forcefield(s, forcefield="pcff")
    return s, p


@pytest.fixture(scope="module", params=list(fixtures()), ids=lambda v: v[0])
def bounded(request, assets):
    name, system = request.param
    return (
        name,
        system,
        prepare_polymer_forcefield(
            system, forcefield="pcff", assets=assets, audit_si_h=name == "disiloxane"
        ),
    )


def test_family_only_never_registers_incomplete_profile(hmds):
    s, prepared = hmds
    p = prepared.validate_integrity(s)
    assert p["family"] == "pcff" and p["candidate_profile"] == CANDIDATE
    assert p["selected_profile"] is None and not REGISTERED_PROFILES
    assert p["typing_status"] == p["charge_status"] == "complete"
    assert p["parameter_completeness"] == "incomplete" and not p["parameter_complete"]
    assert len(p["missing_terms"]) == 378
    assert p["lammps_export_status"] == p["openmm_build_status"] == "blocked"
    assert p["cross_engine_equivalence"] == "not_claimed"
    assert p["scientific_validation"] == p["production_readiness"] == "not_established"
    assert not p["MD_authorized"] and not p["production_validated"]
    assert not p["normalization"] and not p["redistribution"]
    with pytest.raises(Exception, match="parameter_complete"):
        prepared.require(ClaimRequirements(parameter_complete=True))
    with pytest.raises(Exception, match="cross_engine_equivalence"):
        prepared.require(ClaimRequirements(cross_engine_equivalence=True))


def test_bounded_independent_source_and_graph_reconstruction(bounded, assets):
    name, s, p = bounded
    r = p.validate_integrity(s)
    independent = independent_frc(s, r, DISTRIBUTION / "pcff.frc")
    assert independent["status"] == "passed"
    assert (
        len(r["missing_terms"])
        == {"disiloxane": 54, "hmds": 378, "dp2": 528, "dp4": 828, "dp6": 1128}[name]
    )
    assert all(
        c["residual"] in ("0.0000", "0.000", "0.00", "0")
        for c in independent["component_totals"]
    )
    assert r["domain"]["SiH_control"] == (name == "disiloxane")


@pytest.mark.parametrize(
    "family",
    [
        "bond-bond",
        "bond-angle",
        "angle-angle",
        "bond-bond_1_3",
        "middle_bond-torsion_3",
        "end_bond-torsion_3",
        "angle-torsion_3",
        "angle-angle-torsion_1",
    ],
)
def test_missing_class2_family_rejects_completeness(hmds, family):
    _, p = hmds
    r = p.payload
    required = [
        {k: v for k, v in row.items() if k not in ("assignment", "status")}
        for row in r["interaction_inventory"]
    ]
    assigned = [
        row["assignment"]
        for row in r["interaction_inventory"]
        if row["assignment"] is not None and row["family"] != family
    ]
    result = assess(
        required,
        assigned,
        r["conventions"],
        typing_complete=True,
        charge_complete=True,
        source_valid=True,
    )
    assert result["parameter_completeness"] == "incomplete"
    assert any(row["family"] == family for row in result["missing_terms"])
    assert all(
        row["assignment"] is None
        for row in result["missing_terms"]
        if row["family"] == family
    )


def test_dependency_loss_propagates_without_zero_filling(hmds):
    _, p = hmds
    r = p.payload
    required = r["interaction_inventory"]
    assignments = [
        x["assignment"]
        for x in required
        if x["assignment"] and x["family"] != "quartic_bond"
    ]
    q = assess(
        required,
        assignments,
        r["conventions"],
        typing_complete=True,
        charge_complete=True,
        source_valid=True,
    )
    for row in q["interaction_inventory"]:
        if row["family"] in ("bond-bond", "bond-angle"):
            assert (
                row["status"] == "incomplete" and row["reason"] == "missing_dependency"
            )


@pytest.mark.parametrize(
    "mutation",
    [
        "functional_form",
        "unit",
        "source_hash",
        "charge",
        "row_value",
        "zero_filling",
        "bond_bond",
        "bond_angle",
        "angle_angle",
        "torsion_coupling",
        "typing",
        "profile",
        "readiness",
    ],
)
def test_resigned_report_cannot_change_source_assignment(hmds, mutation):
    s, p = hmds
    r = p.payload
    row = next(x for x in r["source_rows"].values() if x["family"] == "quartic_bond")
    if mutation == "functional_form":
        row["functional_form"] = "harmonic"
    elif mutation == "unit":
        row["original_units"][0] = "nanometer"
    elif mutation == "source_hash":
        row["source_sha256"] = "0" * 64
    elif mutation == "charge":
        r["charge_vector"][next(iter(r["charge_vector"]))] += 0.01
    elif mutation == "row_value":
        row["converted_values"][0] = "1"
    elif mutation == "zero_filling":
        r["missing_terms"][0]["assignment"] = {"status": "assigned", "values": [0]}
    elif mutation in ("bond_bond", "bond_angle", "angle_angle", "torsion_coupling"):
        family = {
            "bond_bond": "bond-bond",
            "bond_angle": "bond-angle",
            "angle_angle": "angle-angle",
            "torsion_coupling": "torsion_3",
        }[mutation]
        r["interaction_inventory"] = [
            x for x in r["interaction_inventory"] if x["family"] != family
        ]
    elif mutation == "typing":
        r["typing"][next(iter(r["typing"]))] = "si"
    elif mutation == "profile":
        r["selected_profile"] = CANDIDATE
    else:
        r["production_validated"] = True
    r.pop("identity")
    r["identity"] = digest(r)
    with pytest.raises(Class2AuditError, match="mutation"):
        replace(p, json_text=canonical(r)).validate_integrity(s)


@pytest.mark.parametrize("mutation", ["coordinates", "order", "bond", "box"])
def test_topology_and_execution_input_mutation(hmds, mutation):
    from island.charge_references.records import pack, unpack
    from island.workflows.bundle import system_data, system_from

    s, p = hmds
    data = unpack(pack(system_data(s)))
    if mutation == "coordinates":
        data["coordinates"][next(iter(data["coordinates"]))] = (1.0, 2.0, 3.0)
    elif mutation == "order":
        data["sites"].reverse()
    elif mutation == "bond":
        data["bonds"].pop()
    changed = system_from(data)
    if mutation == "box":
        from island.core import SimulationBox

        changed.box = SimulationBox(40, 40, 40)
    with pytest.raises((ValueError, RuntimeError)):
        p.validate_integrity(changed)


@pytest.mark.parametrize("missing", ["pcff.frc", "pcff.rlb", "pcff_templates.dat"])
def test_incomplete_distribution_rejected(assets, tmp_path, missing):
    shutil.copytree(assets.directory, tmp_path / "sources")
    (tmp_path / "sources" / missing).unlink()
    with pytest.raises(Class2AuditError, match="incomplete_source_archive"):
        PCFFPDMSAssets(tmp_path / "sources").load()


def test_replaced_source_rejected(assets, tmp_path):
    shutil.copytree(assets.directory, tmp_path / "sources")
    with (tmp_path / "sources/pcff.frc").open("ab") as f:
        f.write(b"\n")
    with pytest.raises(Class2AuditError, match="source_hash_mutation"):
        PCFFPDMSAssets(tmp_path / "sources").load()


@pytest.mark.parametrize(
    "smiles",
    [
        "CO[Si](C)(C)C",
        "O[Si](C)(C)O",
        "C[Si]1(C)O[Si](C)(C)O1",
        "CP(C)C",
        "[SiH3]O[SiH3]",
    ],
)
def test_unsupported_terminals_and_chemistry_rejected(smiles, assets):
    s = from_smiles(smiles, random_seed=2026)
    with pytest.raises(Class2AuditError, match="outside_declared_domain"):
        prepare_polymer_forcefield(s, forcefield="pcff", assets=assets)


@pytest.mark.parametrize(
    "family,profile",
    [
        ("opls", CANDIDATE),
        ("gaff", None),
        ("pcff", "island_external_opls_pdms_xiang2025_dft_local_v1"),
    ],
)
def test_wrong_profile_family_rejected(hmds, family, profile):
    s, _ = hmds
    with pytest.raises(Class2AuditError, match="wrong_profile_family"):
        prepare_polymer_forcefield(s, forcefield=family, profile=profile)


def test_pcff_profile_rejected_by_canonical_opls(hmds):
    from island.forcefields.opls_registry import ForceFieldRequest, prepare_forcefield

    s, _ = hmds
    result = prepare_forcefield(s, ForceFieldRequest(family="opls", profile=CANDIDATE))
    assert result.payload["status"] not in ("prepared", "prepared_static_only")
    assert result.payload["selected_profile"] is None


def test_no_backend_for_incomplete_class2(hmds, tmp_path):
    s, p = hmds
    for run in (
        lambda: export_lammps(s, p, tmp_path / "export"),
        lambda: build_openmm(s, p),
    ):
        with pytest.raises(Exception, match="parameter_complete"):
            run()
    assert not (tmp_path / "export").exists()


def test_repeat_graph_and_explicit_caps_are_distinct(assets):
    p = repeat_audit()
    assert p["parameter_completeness"] == p["typing_status"] == "incomplete"
    assert p["missing_terms"][0]["attachment_maps"] == [1, 2]
    s, q = prepare_polymer(
        psmiles="[*:1]O[Si](C)(C)[*:2]",
        forcefield="pcff",
        dp=2,
        head_end_group="[*][Si](C)(C)C",
        tail_end_group="[*]C",
        assets=assets,
    )
    assert q.payload["typing_status"] == "complete" and s.number_of_sites == 37


def test_author_compass_input_is_separate_and_incomplete():
    if not ARCHIVE.is_file():
        pytest.skip("local author archive unavailable")
    a = archive_audit(ARCHIVE)
    assert (
        a["parameter_completeness"] == "incomplete" and a["registered_profiles"] == []
    )
    assert len(a["active_coefficient_rows"]) == 53
    assert a["improper_center_types"] == {6: 26112}
    assert a["silicon_neighbor_triples_without_source_improper_or_AA"] == 12800
    assert any(r["explicit_zero_cells"] for r in a["active_coefficient_rows"])
    assert all("/COMPASS/" in r["source_file"] for r in a["active_coefficient_rows"])


@pytest.mark.parametrize("mode", ["source_hash", "incomplete"])
def test_author_archive_replacement_rejected(tmp_path, mode):
    if not ARCHIVE.is_file():
        pytest.skip("local author archive unavailable")
    raw = ARCHIVE.read_bytes()
    target = tmp_path / "source.zip"
    target.write_bytes(raw + b"mutation" if mode == "source_hash" else raw[:100])
    with pytest.raises(ValueError, match="archive_hash_mutation"):
        archive_audit(target)


def test_relocated_fresh_process_reconstruction(hmds, assets, tmp_path):
    s, p = hmds
    pin = save(s, p, tmp_path / "bundle")
    moved = tmp_path / "relocated"
    shutil.copytree(tmp_path / "bundle", moved / "bundle")
    shutil.copytree(assets.directory, moved / "sources")
    _, q = reconstruct(moved / "bundle", pin, assets=PCFFPDMSAssets(moved / "sources"))
    assert q.json_text == p.json_text
    for seed in (1, 777):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "scripts.phase_4j54",
                "reconstruct",
                "--directory",
                str(moved / "bundle"),
                "--distribution",
                str(moved / "sources"),
                "--manifest-sha256",
                pin,
            ],
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            check=True,
        )
        assert json.loads(result.stdout)["identity"] == p.payload["identity"]


@pytest.mark.parametrize(
    "claim", ["--require-parameter-complete", "--require-cross-engine-equivalence"]
)
def test_cli_claims_are_explicit_and_independent(assets, tmp_path, claim):
    command = [
        sys.executable,
        "-m",
        "island.pcff_pdms",
        "--smiles",
        "C[Si](C)(C)O[Si](C)(C)C",
        "--forcefield",
        "pcff",
        "--distribution",
        str(assets.directory),
        "--json",
        str(tmp_path / "report.json"),
    ]
    assert subprocess.run(command, capture_output=True, check=False).returncode == 0
    assert (
        subprocess.run(command + [claim], capture_output=True, check=False).returncode
        == 1
    )
    p = json.loads((tmp_path / "report.json").read_text())
    assert p["parameter_completeness"] == "incomplete"
    assert p["cross_engine_equivalence"] == "not_claimed"


def test_source_rows_and_conversions_are_fully_pinned(hmds):
    _, p = hmds
    for row in p.payload["source_rows"].values():
        assert row["source_file"] == "pcff.frc"
        assert row["evidence_checksum"] == digest(
            {k: v for k, v in row.items() if k != "evidence_checksum"}
        )
        assert len(row["converted_values"]) == len(row["converted_units"])
        assert row["citation"]["revision"] and row["location"]["line"] > 0
        assert row["doi"] is None  # never relabel FRC rows as a different paper version
    assert p.payload["conventions"]["lj_mixing"]["value"]["lines"]
    assert p.payload["conventions"]["coulomb_14"]["value"] is None


def test_author_rows_retain_exact_roles_units_and_explicit_zero_evidence():
    if not ARCHIVE.is_file():
        pytest.skip("local author archive unavailable")
    from scripts.audit_pcff_pdms_class2_sources import signed

    a = archive_audit(ARCHIVE)
    for row in a["active_coefficient_rows"]:
        assert len(row["original_values"]) == len(row["original_units"])
        assert row["converted_units"] == row["original_units"]
        assert row["converted_values"] == row["original_values"]
        assert row["labels"] and all(isinstance(p, list) for p in row["labels"])
        assert signed({k: v for k, v in row.items() if k != "evidence_checksum"}) == row
    assert any("dihedral" in s["raw_text"] for s in a["settings"])
