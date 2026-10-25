"""Opt-in source completeness, independent Si mapping and persistence negatives."""

import copy
import os
import tarfile
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    PCFFSiPAssessment,
    assess_pcff_sip,
    load_pcff_sip_assessment,
    load_pcff_source,
    prepare_pcff_sip,
)
from island.forcefields.pcff.source import PCFFSource
from scripts.pcff_j41_reference import audit_queries, inventory


@pytest.fixture(scope="module")
def source():
    path = Path(
        os.environ.get(
            "ISLAND_J41_SOURCE", "../lammps/lammps/tools/msi2lmp/frc_files/pcff.frc"
        )
    )
    if path.is_file():
        return load_pcff_source(path)
    archive = Path("docs/evidence/phase_4j41/execution.tar.xz")
    if archive.is_file():
        with tarfile.open(archive) as tar:
            raw = tar.extractfile("study/pcff.frc").read()
        return PCFFSource(raw, sha256(raw).hexdigest())
    pytest.skip("Reviewed PCFF source asset unavailable")


@pytest.fixture(scope="module")
def silane(source):
    system = from_smiles("[SiH4]", random_seed=2026)
    assessment = assess_pcff_sip(system, source, profile="island_pcff_silane_source_v1")
    return system, assessment


def test_silane_complete_source_only_model(silane):
    system, result = silane
    model = prepare_pcff_sip(system, result)
    assert result.payload["complete"]
    assert model.assignment.payload["parameter_coverage_complete"]
    terms = model.payload["terms"]
    assert len([t for t in terms if t["family"] == "angle-angle"]) == 12
    assert all(t["origin"] == "source_row" and t["source_rows"] for t in terms)
    assert result.payload["charge_origin"] == "PCFF/source_native_bond_increments"
    assert result.payload["production_validated"] is False
    assert result.payload["simulation_readiness"] == "not_established"


def test_independent_inventory_and_source_reconstruction(source, silane):
    system, assessment = silane
    raw = inventory(source.raw)
    assert raw["primary_labels"] == ["p", "p=", "si", "sio", "sz"]
    assert {"hsi", "osi", "oas", "osh", "oss"} <= set(raw["descriptive_labels"])
    assert raw["selected_counts"]["bond_increments"] > 0
    checked = audit_queries(system, assessment, source.raw)
    assert checked["typing_equal"] and checked["inventory_equal"]
    assert checked["charges"]["components"][0]["residual"] == "0.0000"


@pytest.mark.parametrize(
    "smiles,profile,missing",
    [
        ("C[SiH3]", "organosilicon_source", "bond-bond_1_3"),
        ("C[Si](C)(C)O[Si](C)(C)C", "siloxane_source", "angle-angle"),
    ],
)
def test_missing_couplings_reject_publication(source, smiles, profile, missing):
    system = from_smiles(smiles, random_seed=2026)
    result = assess_pcff_sip(system, source, profile="island_pcff_" + profile + "_v1")
    assert not result.payload["complete"]
    assert result.payload["inspection"]["charges_complete"]
    queries = result.payload["inspection"]["interaction_queries"]
    assert any(
        q["family"] == missing and q["resolution"]["status"] == "missing"
        for q in queries
    )
    assert audit_queries(system, result, source.raw)["inventory_equal"]
    with pytest.raises(PCFFError, match="Incomplete Si/P"):
        prepare_pcff_sip(system, result)


@pytest.mark.parametrize(
    "smiles,profile",
    [
        ("COP(=O)(OC)OC", "phosphate_ester_audit"),
        ("CP(=O)(OC)OC", "phosphonate_audit"),
        ("[PH3]", "phosphonate_audit"),
    ],
)
def test_p_oxidation_template_conflict_never_silently_typed(source, smiles, profile):
    system = from_smiles(smiles, random_seed=2026)
    result = assess_pcff_sip(system, source, profile="island_pcff_" + profile + "_v1")
    assert not result.payload["complete"]
    assert result.payload["typing"] is None
    assert any(
        "phosphorus_oxidation" in d["reason"] for d in result.payload["diagnostics"]
    )
    with pytest.raises(PCFFError, match="Incomplete Si/P"):
        prepare_pcff_sip(system, result)


@pytest.mark.parametrize(
    "smiles,reason",
    [
        ("[SiH2]", "unsupported_silicon_valence"),
        ("[SiH3-]", "charged_ionic"),
        ("C", "outside_selected_source_profile"),
    ],
)
def test_unsupported_domains(source, smiles, reason):
    system = from_smiles(smiles, random_seed=2026)
    result = assess_pcff_sip(system, source, profile="island_pcff_silane_source_v1")
    assert not result.payload["complete"]
    assert any(reason in d["reason"] for d in result.payload["diagnostics"])


@pytest.mark.parametrize(
    "field",
    [
        "source",
        "graph_identity",
        "atom_order",
        "coordinates",
        "charge_origin",
        "typing_identity",
        "complete",
        "diagnostics",
        "special_pairs",
        "inspection",
        "charges",
        "component_totals",
        "source_rows",
    ],
)
def test_persistence_mutations_fail_closed(silane, field, tmp_path):
    system, assessment = silane
    payload = copy.deepcopy(assessment.payload)
    if field == "source":
        payload[field]["sha256"] = "0" * 64
    elif field == "atom_order":
        payload[field].reverse()
    elif field == "coordinates":
        payload[field][0][1][0] += 0.1
    elif field == "special_pairs":
        payload[field]["coulomb"][2] = 0.5
    elif field == "charges":
        payload["inspection"]["native_charge_record"]["native_charge_record"][
            "partial_charges"
        ][1] += 0.01
    elif field == "component_totals":
        payload["inspection"]["native_charge_record"]["native_charge_record"][
            "components"
        ][0]["charge"] = 1.0
    elif field == "source_rows":
        payload["inspection"]["source_rows"] = {}
    elif field == "inspection":
        payload[field]["interaction_queries"].pop()
    elif field == "complete":
        payload[field] = False
    elif field == "diagnostics":
        payload[field] = [{"reason": "invented"}]
    else:
        payload[field] = "mutated"
    path = tmp_path / "record.json"
    path.write_text(pack(payload))
    with pytest.raises(PCFFError, match="Contradictory"):
        load_pcff_sip_assessment(path, assessment.source, system=system)


def test_final_graph_mutation_rejected(silane):
    system, assessment = silane
    changed = system.copy()
    changed.topology.remove_site(5)
    with pytest.raises(PCFFError):
        prepare_pcff_sip(changed, assessment)


def test_wrong_source_and_incomplete_interactions_rejected(silane):
    system, assessment = silane
    raw = assessment.source.raw.replace(b"silicon atom", b"SILICON atom", 1)
    wrong = PCFFSource(raw, sha256(raw).hexdigest())
    with pytest.raises(PCFFError, match="not audited|hash mismatch"):
        assess_pcff_sip(system, wrong, profile="island_pcff_silane_source_v1")
    model = prepare_pcff_sip(system, assessment)
    p = model.assignment.payload
    p["assignments"].pop()
    incomplete = replace(model.assignment, json_text=pack(p))
    with pytest.raises(PCFFError):
        incomplete.validate_integrity(system)


def test_unregistered_profile_and_roundtrip(silane, tmp_path):
    system, assessment = silane
    path = tmp_path / "assessment.json"
    path.write_text(assessment.json_text)
    loaded = load_pcff_sip_assessment(path, assessment.source, system=system)
    assert loaded.identity == assessment.identity
    assert loaded.json_text == assessment.json_text
    with pytest.raises(PCFFError, match="Unregistered"):
        assess_pcff_sip(system, assessment.source, profile="universal_Si_P")
    with pytest.raises(PCFFError):
        replace(
            PCFFSiPAssessment(assessment.json_text, assessment.source), json_text="{}"
        ).validate_integrity(system)
