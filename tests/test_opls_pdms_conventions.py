"""Figure 7 partitioning must not manufacture missing execution conventions."""

import json
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import FrozenInstanceError
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.forcefields.external_pdms import digest
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_conventions import (
    OPLSPDMSConventionAudit,
    OPLSPDMSConventionsBlockedError,
    assess_opls_pdms_conventions,
    load_opls_pdms_conventions,
    source_evidence,
)
from island.forcefields.oplsaa.pdms_literature import (
    SOURCE_FILES,
    OPLSPDMSLiteratureError,
    build_opls_pdms_literature,
)
from island.forcefields.preparation import ForceFieldRequest, ForceFieldRequestError
from island.workflows.bundle import system_from
from scripts.assess_opls_pdms_conventions import fixtures

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "docs/evidence/phase_4j44/source-audit.json"
SELECTION = ROOT / "docs/evidence/phase_4j45/selection.json"


@pytest.fixture
def supplement():
    return build_opls_pdms_literature(
        DomainSelection(SELECTION.read_text(), CATALOG.read_text())
    )


@pytest.fixture
def references(tmp_path):
    target = tmp_path / "sources"
    for name in SOURCE_FILES:
        original = ROOT / "island_ffref" / name
        if not original.is_file():
            pytest.skip("Requires the two user-supplied local Xiang PDFs")
        dest = target / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, dest)
    return target


@pytest.fixture
def cases():
    return {name: (system, oxygen) for name, system, oxygen in fixtures()}


def assess(case, supplement, references):
    system, oxygen = case
    return assess_opls_pdms_conventions(
        system, supplement, references=references, e_oxygen_site_id=oxygen
    )


@pytest.mark.parametrize("name", ["hmds", "tms_capped_dp2", "tms_capped_dp4"])
def test_source_partition_is_exact_without_charge_repair(
    name, cases, supplement, references
):
    system, oxygen = cases[name]
    result = assess((system, oxygen), supplement, references)
    p = result.payload
    result.validate_integrity(system, references=references)
    assert p["status"] == "parameter_incomplete"
    q = p["charge_partition"]
    assert Decimal(q["exact_charge_sum"]) == 0
    assert sum(map(Decimal, p["exact_charge_vector"])) == 0
    assert not q["normalization_or_redistribution"]
    assert not q["author_atom_mapping_verified"]
    assert Counter(r["element"] for r in q["rows"] if r["fragment"] == "E") == {
        "Si": 2,
        "O": 1,
        "C": 6,
        "H": 18,
    }
    assert [r["site_id"] for r in q["rows"]] == list(system.topology.sites)
    source_charges = {
        e: supplement.lookup("partial_charge", [e])[0]["original_values"]
        for e in ("O", "Si", "C", "H")
    }
    for r in q["rows"]:
        assert r["charge"] == source_charges[r["element"]][int(r["fragment"] == "E")]
    assert p["stages"]["source_consistent_charge_partition"]
    assert not p["stages"]["coefficient_units_complete"]
    assert (
        supplement.identity
        == "511948798d014d4f9bd2521e50761644ccf49f135a713b49de5cb4d17be22538"
    )
    assert all(r["source_units"] is None for r in supplement.payload["rows"])


def test_assigning_both_terminal_oxygens_is_not_the_source_partition(
    cases, supplement, references
):
    p = assess(cases["tms_capped_dp4"], supplement, references).payload
    rows = p["charge_partition"]["rows"]
    wrong = [
        Decimal("-0.578")
        if r["element"] == "O" and r["fragment"] == "M" and r["site_id"] == 13
        else Decimal(r["charge"])
        for r in rows
    ]
    assert sum(wrong) == Decimal("0.022")
    assert sum(map(Decimal, p["exact_charge_vector"])) == 0


@pytest.mark.parametrize(
    "name", ["tms_capped_dp2_unoriented", "tms_capped_dp4_unoriented"]
)
def test_longer_chain_requires_explicit_orientation(
    name, cases, supplement, references
):
    p = assess(cases[name], supplement, references).payload
    assert p["status"] == "typing_ambiguous"
    assert p["exact_charge_vector"] is None
    assert p["diagnostics"][0]["stage"] == "explicit_E_oxygen_required"
    assert len(p["diagnostics"][0]["candidate_site_ids"]) == 2


def test_explicit_opposite_orientation_is_a_distinct_graph_bound_record(
    cases, supplement, references
):
    system, oxygen = cases["tms_capped_dp4"]
    first = assess((system, oxygen), supplement, references).payload
    second = assess((system, 13), supplement, references).payload
    assert first["graph_identity"] == second["graph_identity"]
    assert first["identity"] != second["identity"]
    assert first["exact_charge_vector"] != second["exact_charge_vector"]
    assert sum(map(Decimal, second["exact_charge_vector"])) == 0
    assert (
        second["charge_partition"]["selection_origin"]
        == "caller_explicit_graph_bound_E_oxygen"
    )


@pytest.mark.parametrize("oxygen", [True, "1", [1, 13], 2, 5, -1])
def test_no_arbitrary_or_internal_E_oxygen(oxygen, cases, supplement, references):
    system, _ = cases["tms_capped_dp4"]
    p = assess((system, oxygen), supplement, references).payload
    assert p["status"] == "typing_blocked"
    assert p["diagnostics"][0]["stage"] == "invalid_E_oxygen_site"
    assert p["exact_charge_vector"] is None


@pytest.mark.parametrize("change", ["valence", "charge", "radical", "element", "si_h"])
def test_unsupported_graphs_remain_classified(change, cases, supplement, references):
    system, oxygen = cases["hmds"]
    if change == "si_h":
        system, oxygen = cases["unsupported_sih3"]
    elif change == "valence":
        bond = next(iter(system.topology.bonds))
        system.topology.remove_bond(*bond)
        system.topology.rebuild_derived_interactions()
    else:
        atom = next(iter(system.topology.sites.values()))
        if change == "charge":
            atom.formal_charge = 1
        elif change == "radical":
            atom.metadata["radical_electrons"] = 1
        else:
            atom.element, atom.atomic_number = "P", 15
    result = assess((system, oxygen), supplement, references)
    assert result.payload["status"] == "typing_blocked"
    assert result.payload["exact_charge_vector"] is None
    result.validate_integrity(system, references=references)
    with pytest.raises(OPLSPDMSConventionsBlockedError) as error:
        result.require_executable()
    assert error.value.result is result
    assert error.value.status == "typing_blocked"


@pytest.mark.parametrize(
    "change",
    [
        "charge",
        "fragment",
        "unit",
        "mixing",
        "scaling",
        "evidence",
        "source",
        "authority",
        "profile",
        "status",
        "bool_type",
        "missing_gap",
    ],
)
def test_resigning_does_not_forge_conventions(change, cases, supplement, references):
    result = assess(cases["hmds"], supplement, references)
    p = result.payload
    if change == "charge":
        p["exact_charge_vector"][0] = "0.0"
    elif change == "fragment":
        p["charge_partition"]["rows"][0]["fragment"] = "M"
    elif change == "unit":
        p["source_evidence"]["unit_evidence"]["coefficient_units"] = "real"
    elif change == "mixing":
        p["missing_dependencies"]["mixing"] = "geometric assumed"
    elif change == "scaling":
        p["missing_dependencies"]["special_pairs"] = [0, 0, 0.5]
    elif change == "evidence":
        p["source_evidence"]["figure7"]["image_sha256"] = "0" * 64
    elif change == "source":
        p["source_evidence"]["source_hashes"][SOURCE_FILES[0]] = "0" * 64
    elif change == "authority":
        p["execution_authorized"] = True
    elif change == "profile":
        p["profile_identity"] = "native_oplsaa"
    elif change == "status":
        p["status"] = "prepared"
    elif change == "bool_type":
        p["execution_authorized"] = 0
    else:
        p["missing_dependencies"].pop("improper_applicability")
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        OPLSPDMSConventionAudit(pack(p), supplement)
    assert error.value.status == "source_incompatible"


@pytest.mark.parametrize(
    "change", ["order", "component", "graph", "polymer", "orientation"]
)
def test_stale_final_graph_or_order_rejected(change, cases, supplement, references):
    system, oxygen = cases["tms_capped_dp4"]
    result = assess((system, oxygen), supplement, references)
    if change == "order":
        system.topology.sites = dict(reversed(list(system.topology.sites.items())))
    elif change == "component":
        next(iter(system.topology.sites.values())).metadata["chain_id"] = "stale"
    elif change == "graph":
        next(iter(system.topology.sites.values())).formal_charge = 1
    elif change == "polymer":
        system.metadata["polymer"]["end_groups"]["head"] = "[*]O"
    else:
        system.metadata["polymer"]["head_site_id"] = 13
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        result.validate_integrity(system, references=references)
    assert error.value.status == "source_incompatible"


def test_only_coordinates_can_change(cases, supplement, references):
    system, _ = cases["hmds"]
    result = assess(cases["hmds"], supplement, references)
    for i in system.topology.sites:
        system.coordinates.set(i, system.coordinates.get(i) + 1.5)
    result.validate_integrity(system, references=references)


def test_immutability_and_no_backend_authority(
    cases, supplement, references, monkeypatch
):
    def forbidden(*args, **kwargs):
        pytest.fail(
            "Data audit must not reach a numerical or native force-field adapter"
        )

    monkeypatch.setattr("island.forcefields.preparation.prepare_forcefield", forbidden)
    monkeypatch.setattr(
        "island.forcefields.oplsaa.adapter.check_foyer_installation", forbidden
    )
    result = assess(cases["hmds"], supplement, references)
    with pytest.raises(FrozenInstanceError):
        result.json_text = "forged"
    result.payload["charge_partition"]["rows"].clear()
    source_evidence(supplement)["figure7"]["Qeq_Charge"]["Si"][0] = "99"
    assert assess(cases["hmds"], supplement, references) == result
    assert result.payload["charge_partition"]["rows"]
    with pytest.raises(ForceFieldRequestError):
        ForceFieldRequest("oplsaa", result)
    with pytest.raises(OPLSPDMSConventionsBlockedError) as error:
        result.require_executable()
    assert error.value.result is result
    assert any(d.get("dependency") == "special_pairs" for d in error.value.diagnostics)


@pytest.mark.parametrize("change", ["pin", "source", "missing"])
def test_source_hashes_and_independent_file_pin_required(
    change, cases, supplement, references, tmp_path
):
    result = assess(cases["hmds"], supplement, references)
    path = tmp_path / "record.json"
    path.write_text(result.json_text)
    pin = sha256(path.read_bytes()).hexdigest()
    if change == "pin":
        pin = None
    elif change == "source":
        (references / SOURCE_FILES[0]).write_bytes(b"changed source")
    else:
        (references / SOURCE_FILES[0]).unlink()
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        load_opls_pdms_conventions(
            path, pin, supplement=supplement, references=references
        )
    assert error.value.status == (
        "asset_unavailable" if change == "missing" else "source_incompatible"
    )


@pytest.mark.parametrize(
    "name", ["hmds", "tms_capped_dp2_unoriented", "unsupported_sih3"]
)
def test_relocation_and_fresh_process_keep_all_blocked_statuses(
    name, cases, supplement, references, tmp_path
):
    result = assess(cases[name], supplement, references)
    moved = tmp_path / "relocated"
    moved.mkdir()
    path = moved / "record.json"
    path.write_text(result.json_text)
    selection, catalog = moved / "selection.json", moved / "catalog.json"
    selection.write_text(supplement.selection.json_text)
    catalog.write_text(supplement.selection.catalog_text)
    shutil.move(str(references), moved / "references")
    args = [
        sys.executable,
        "scripts/assess_opls_pdms_conventions.py",
        "--selection",
        str(selection),
        "--catalog",
        str(catalog),
        "--references",
        str(moved / "references"),
        "--reconstruct",
        str(path),
        "--sha256",
        sha256(path.read_bytes()).hexdigest(),
        "--require-executable",
    ]
    run = subprocess.run(args, text=True, capture_output=True, check=False)
    assert run.returncode == 1 and not run.stderr
    output = json.loads(run.stdout)
    assert output["status"] == result.payload["status"]
    assert output["identity"] == result.payload["identity"]
    assert output["backend_executions"] == output["execution_profiles"] == []
    loaded = load_opls_pdms_conventions(
        path,
        sha256(path.read_bytes()).hexdigest(),
        supplement=supplement,
        references=moved / "references",
    )
    loaded.validate_integrity(
        system_from(result.payload["system"]), references=moved / "references"
    )
    assert not any(p.suffix in (".data", ".chk", ".dcd") for p in moved.rglob("*"))
