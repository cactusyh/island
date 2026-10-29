"""New graph mapping, source permission separation and longrange controls."""

import json
import os
import shutil
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.forcefields.external_pdms import digest
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.opls_external.pdms_mapping import (
    load_mapping,
    map_author_pdms,
    source_rights,
)
from island.forcefields.oplsaa.pdms_author import load_author_opls_audit
from island.forcefields.oplsaa.pdms_literature import (
    OPLSPDMSLiteratureError,
    build_opls_pdms_literature,
)
from island.workflows.bundle import system_data
from scripts.assess_opls_pdms_conventions import fixtures
from scripts.qualify_pdms_longrange_backend import BOUNDS, qualify, switch

ROOT = Path(__file__).resolve().parents[1]
PIN = "33fb0a3f5be078b29532123c9563b47068fd54ca827fc0df0f0160454e1d0e49"


@pytest.fixture(scope="module")
def author():
    sources = os.environ.get("ISLAND_PDMS_AUTHOR_SOURCES")
    docs = os.environ.get("ISLAND_LAMMPS_DOCS")
    if not sources or not docs:
        pytest.skip(
            "Explicit local author sources and engine docs required; no discovery"
        )
    supplement = build_opls_pdms_literature(
        DomainSelection(
            (ROOT / "docs/evidence/phase_4j45/selection.json").read_text(),
            (ROOT / "docs/evidence/phase_4j44/source-audit.json").read_text(),
        )
    )
    return load_author_opls_audit(
        ROOT / "docs/evidence/opls_pdms_author_v3/record.json",
        PIN,
        sources=Path(sources),
        supplement=supplement,
        references=ROOT / "island_ffref",
        engine_docs=Path(docs),
    )


@pytest.fixture(scope="module")
def cases():
    return {name: (system, oxygen) for name, system, oxygen in fixtures()}


def test_publication_grant_does_not_license_other_assets():
    rights = source_rights(ROOT / "island_ffref")
    assert rights["paper"]["license"] == "CC-BY-4.0"
    assert rights["author_repository"]["use_permission"] == "not_established"
    assert (
        rights["supporting_information"]["paper_grant_automatically_transferred"]
        is False
    )


@pytest.mark.parametrize("name", ["hmds", "tms_capped_dp2", "tms_capped_dp4"])
def test_explicit_graph_partition_preserves_every_source_charge(author, cases, name):
    system, oxygen = cases[name]
    result = map_author_pdms(
        system, author, references=ROOT / "island_ffref", e_oxygen_site_id=oxygen
    )
    result.validate_integrity(system)
    p = result.payload
    assert p["status"] == "asset_unavailable"
    assert p["atom_order"] == list(system.topology.sites)
    assert len([r for r in p["mapping"] if r["fragment"] == "E"]) == 27
    assert len(p["terminal_silicon"]) == 2
    assert sum(Fraction(r["printed_Figure7_DFT_charge"]) for r in p["mapping"]) == 0
    actual = sum(Fraction(r["exact_author_serialized_charge"]) for r in p["mapping"])
    assert actual == Fraction(p["charge_policy"]["exact_serialized_sum_fraction"])
    expected_residual = {
        "hmds": Fraction(3, 10**16),
        "tms_capped_dp2": Fraction(1, 5 * 10**15),
        "tms_capped_dp4": Fraction(0),
    }
    assert (
        actual == expected_residual[name]
    )  # Preserve source sums, including exact zero.
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        result.require_executable(system)
    assert error.value.status == "asset_unavailable" and error.value.result is result
    assert p["execution_profile"] is None and not p["numerical_MD_run"]


@pytest.mark.parametrize(
    "name", ["tms_capped_dp2_unoriented", "tms_capped_dp4_unoriented"]
)
def test_ambiguous_orientation_never_uses_site_order(author, cases, name):
    system, _ = cases[name]
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        map_author_pdms(system, author, references=ROOT / "island_ffref")
    assert error.value.status == "typing_ambiguous"
    assert error.value.diagnostics[0]["stage"] == "explicit_E_oxygen_required"


@pytest.mark.parametrize(
    "mutation",
    [
        "valence",
        "formal_charge",
        "radical",
        "phosphorus",
        "wrong_anchor",
        "sih3",
        "mass",
    ],
)
def test_unsupported_graphs_are_classified(author, cases, mutation):
    system, oxygen = deepcopy(cases["tms_capped_dp4"])
    if mutation == "valence":
        system.topology.remove_bond(*next(iter(system.topology.bonds)))
    elif mutation == "formal_charge":
        system.topology.sites[1].formal_charge = 1
    elif mutation == "radical":
        system.topology.sites[1].metadata["radical_electrons"] = 1
    elif mutation == "phosphorus":
        atom = next(a for a in system.topology.sites.values() if a.element == "Si")
        atom.element, atom.atomic_number = "P", 15
    elif mutation == "mass":
        system.topology.sites[1].mass = 99.0
    elif mutation == "wrong_anchor":
        oxygen = True  # bool is not an exact integer stable atom ID.
    else:
        system, oxygen = cases["unsupported_sih3"]
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        map_author_pdms(
            system, author, references=ROOT / "island_ffref", e_oxygen_site_id=oxygen
        )
    assert error.value.status == "typing_blocked"


@pytest.mark.parametrize(
    "field",
    [
        "status",
        "execution_authorized",
        "execution_profile",
        "graph_identity",
        "charge",
        "rights",
        "units",
        "source_hashes",
    ],
)
def test_resigned_forged_records_cannot_promote_execution(author, cases, field):
    system, oxygen = cases["tms_capped_dp2"]
    result = map_author_pdms(
        system, author, references=ROOT / "island_ffref", e_oxygen_site_id=oxygen
    )
    p = result.payload
    if field == "charge":
        p["mapping"][0]["exact_author_serialized_charge"] = "0"
    elif field == "rights":
        p["source_rights"]["author_repository"]["use_permission"] = "authorized"
    elif field == "units":
        p["charge_policy"]["unit"] = "coulomb"
    elif field == "source_hashes":
        p["source_hashes"]["README.md"] = "0" * 64
    else:
        p[field] = {
            "status": "prepared",
            "execution_authorized": True,
            "execution_profile": "fake",
            "graph_identity": "0" * 64,
        }[field]
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    forged = replace(result, json_text=pack(p))
    with pytest.raises(OPLSPDMSLiteratureError):
        forged.require_executable(system)


@pytest.mark.parametrize("mutation", ["order", "component", "graph"])
def test_stale_final_graph_rejected(author, cases, mutation):
    system, oxygen = deepcopy(cases["tms_capped_dp4"])
    result = map_author_pdms(
        system, author, references=ROOT / "island_ffref", e_oxygen_site_id=oxygen
    )
    if mutation == "order":
        system.topology.sites = dict(reversed(list(system.topology.sites.items())))
    elif mutation == "component":
        system.topology.sites[1].metadata["chain_id"] = "stale"
    else:
        system.topology.remove_bond(*next(iter(system.topology.bonds)))
    with pytest.raises(OPLSPDMSLiteratureError):
        result.validate_integrity(system)


@pytest.mark.parametrize("name", ["hmds", "tms_capped_dp2", "tms_capped_dp4"])
def test_relocated_fresh_process_cannot_execute(author, cases, tmp_path, name):
    system, oxygen = cases[name]
    result = map_author_pdms(
        system, author, references=ROOT / "island_ffref", e_oxygen_site_id=oxygen
    )
    record = tmp_path / "record.json"
    record.write_text(result.json_text)
    system_file = tmp_path / "system.json"
    system_file.write_text(pack(system_data(system)))
    refs = tmp_path / "relocated-references"
    for relative in ("pcff/jp4c08471.pdf", "pcff/jp4c08471_si_001.pdf"):
        p = refs / relative
        p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "island_ffref" / relative, p)
    pin = sha256(record.read_bytes()).hexdigest()
    assert (
        load_mapping(
            record, pin, system, author, references=refs, e_oxygen_site_id=oxygen
        ).json_text
        == result.json_text
    )
    command = [
        sys.executable,
        str(ROOT / "scripts/audit_pdms_execution_boundary.py"),
        "--sources",
        str(author.sources),
        "--references",
        str(refs),
        "--engine-docs",
        str(author.engine_docs),
        "--selection",
        str(ROOT / "docs/evidence/phase_4j45/selection.json"),
        "--catalog",
        str(ROOT / "docs/evidence/phase_4j44/source-audit.json"),
        "--author-record",
        str(ROOT / "docs/evidence/opls_pdms_author_v3/record.json"),
        "--author-record-sha256",
        PIN,
        "--reconstruct",
        str(record),
        "--sha256",
        pin,
        "--system",
        str(system_file),
        "--require-executable",
    ]
    if oxygen is not None:
        command += ["--e-oxygen-site-id", str(oxygen)]
    process = subprocess.run(
        command, cwd=ROOT, capture_output=True, text=True, check=False
    )
    assert process.returncode == 1, process.stderr
    assert json.loads(process.stdout)["status"] == "asset_unavailable"
    assert list(tmp_path.glob("*.chk")) + list(tmp_path.glob("*.dcd")) == []


def test_charmm_switch_is_not_openmm_builtin_switch():
    assert switch(9) == 1 and switch(11) == 0
    assert switch(10) == pytest.approx(0.53746875)
    assert switch(10) != pytest.approx(0.5)
    for r in (9, 11):
        derivative = (switch(r + 1e-5) - switch(r - 1e-5)) / 2e-5
        assert abs(derivative) < 1e-5


def test_wrong_executable_pin_prevents_launch(tmp_path, monkeypatch):
    exe = tmp_path / "exe"
    exe.write_text("never run")
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: pytest.fail("Executable launched before hash validation"),
    )
    with pytest.raises(ValueError):
        qualify(exe, "0" * 64, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_actual_longrange_equations_with_independent_forces(tmp_path):
    exe = os.environ.get("ISLAND_LAMMPS_EXECUTABLE")
    pin = os.environ.get("ISLAND_LAMMPS_EXECUTABLE_SHA256")
    if not exe or not pin:
        pytest.skip("Explicit pinned LAMMPS executable required")
    p = qualify(exe, pin, tmp_path / "control")
    assert p["software_equations_verified"]
    assert not p["author_default_accuracy_control_verified"]
    assert not p["physical_PDMS_executable"] and not p["PDMS_MD_run"]
    assert p["controls"][-1]["status"] == "verified"
    assert (
        max(p["controls"][-1]["component_absolute_errors"].values())
        < BOUNDS["tight_energy_kj_mol"]
    )
