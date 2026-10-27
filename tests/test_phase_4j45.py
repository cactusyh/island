"""Source-bound closure controls; no invented numerical parameter fixtures."""

import json
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError
from hashlib import sha256
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields.external_pdms import digest
from island.forcefields.literature_closure import (
    STAGES,
    DomainClosure,
    DomainClosureBlockedError,
    DomainSelection,
    close_domain,
)
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j45 import execute, fixtures, reconstruct

EVIDENCE = Path(__file__).resolve().parents[1] / "docs/evidence"
SELECTION = EVIDENCE / "phase_4j45/selection.json"
CATALOG = EVIDENCE / "phase_4j44/source-audit.json"


@pytest.fixture
def selection():
    return DomainSelection(SELECTION.read_text(), CATALOG.read_text())


def resign(p):
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    return DomainClosure(pack(p))


@pytest.mark.parametrize(
    "name,system", list(fixtures()), ids=lambda x: x if isinstance(x, str) else "graph"
)
def test_every_fixture_is_graph_bound_and_classified(name, system, selection):
    result = close_domain(system, selection)
    result.validate_integrity(system, selection)
    p = result.payload
    assert p["missing_rows"] and p["exact_charge_vector"] is None
    assert p["profile_identity"] is p["parameter_assignment"] is None
    assert set(p["stages"]) == set(STAGES)
    assert all(v is False for v in p["stages"].values())
    assert set(p["missing_rows"]) == set(p["requirements"]["requests"])
    assert p["limitations"]["production_validated"] is False
    assert p["limitations"]["simulation_readiness"] == "not_established"
    with pytest.raises(DomainClosureBlockedError) as exc:
        result.require_executable()
    assert exc.value.result is result and exc.value.status == p["status"]
    assert exc.value.diagnostics == p["missing_rows"]
    if name == "hexamethyldisiloxane":
        assert p["status"] == "parameter_incomplete"
        assert not p["unsupported_environments"]
        assert all(
            v
            for k, v in p["literal_row_candidates"].items()
            if k.startswith(("bonds:", "angles:", "rb:"))
        )
        assert sum(bool(v) for v in p["literal_row_candidates"].values()) > 50
    else:
        assert p["status"] == "typing_ambiguous"
        assert p["unsupported_environments"]


def test_selection_reuses_j44_without_promoting_candidates(selection):
    p = selection.payload
    assert len(selection.catalog["files"]) == 52
    assert len(selection.catalog["extracted_row_candidates"]) == 136
    assert len(p["candidate_rows"]) == 17
    assert len(p["alternatives"]) == 8
    assert p["accepted_rows"] == p["execution_profiles"] == []
    assert p["profile_identity"] is None
    assert all(r["converted_values"] is None for r in p["candidate_rows"])
    assert any("Si-H" in v["reason"] for v in p["alternatives"])
    # Explicit source zeros remain evidence; they confer no assignment authority.
    assert any("0" in r["original_values"] for r in p["candidate_rows"])


@pytest.mark.parametrize(
    "field",
    [
        "candidate_rows",
        "required_unresolved_policies",
        "profile_identity",
        "alternatives",
        "citations",
    ],
)
def test_selection_mutations_cannot_supply_missing_authority(selection, field):
    p = selection.payload
    p[field] = [] if isinstance(p[field], list) else "forged"
    with pytest.raises(ValueError, match="selection_evidence"):
        DomainSelection(json.dumps(p), selection.catalog_text)


def test_catalog_hash_is_independent(selection):
    with pytest.raises(ValueError, match="catalog_mutated"):
        DomainSelection(selection.json_text, selection.catalog_text + " ")


@pytest.mark.parametrize(
    "mutation",
    [
        "graph",
        "atom_order",
        "component",
        "units",
        "rows",
        "charges",
        "profile",
        "status",
        "stage",
        "source",
        "zero",
        "checksum",
    ],
)
def test_mutations_cannot_promote_or_rebind_blocked_record(selection, mutation):
    system = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
    original = close_domain(system, selection)
    p = original.payload
    if mutation == "graph":
        p["requirements"]["graph"]["bonds"].pop()
    elif mutation == "atom_order":
        p["requirements"]["atom_order"].reverse()
    elif mutation == "component":
        first = next(iter(p["requirements"]["graph"]["molecule_membership"]))
        p["requirements"]["graph"]["molecule_membership"][first] = "stale"
    elif mutation == "units":
        p["unresolved_policies"]["units"] = "units real inferred"
    elif mutation == "rows":
        p["missing_rows"].pop(next(iter(p["missing_rows"])))
    elif mutation == "charges":
        p["exact_charge_vector"] = [0] * system.number_of_sites
    elif mutation == "profile":
        p["profile_identity"] = "native_oplsaa"
    elif mutation == "status":
        p["status"] = "executable"
    elif mutation == "stage":
        p["stages"]["source_complete"] = True
    elif mutation == "source":
        p["selection_sha256"] = "f" * 64
    elif mutation == "zero":
        p["accepted_rows"] = [{"missing_term": 0}]
    else:
        p["identity"] = "a" * 64
    with pytest.raises(ValueError):
        altered = DomainClosure(pack(p)) if mutation == "checksum" else resign(p)
        altered.validate_integrity(system, selection)


@pytest.mark.parametrize(
    "smiles",
    [
        "[SiH3]O[SiH3]",
        "C[SiH3]",
        "[SiH2]O[SiH3]",
        "[O-][Si](C)(C)C",
        "COP(=O)(O)OC",
        "CO[Si](C)(C)C",
    ],
)
def test_unsupported_valence_charge_radical_and_foreign_types(selection, smiles):
    p = close_domain(from_smiles(smiles), selection).payload
    assert p["status"] == "typing_ambiguous"
    assert p["unsupported_environments"]
    assert not p["execution_profiles"]


def test_wrong_graph_and_molecule_membership_rejected(selection):
    system = from_smiles("C[Si](C)(C)O[Si](C)(C)C")
    result = close_domain(system, selection)
    other = from_smiles("CO[Si](C)(C)C")
    with pytest.raises(ValueError, match="stale_graph"):
        result.validate_integrity(other, selection)
    system.topology.sites[1].metadata["chain_id"] = "mutated"
    with pytest.raises(ValueError, match="stale_graph"):
        result.validate_integrity(system, selection)


def test_sources_missing_and_mutated_are_classified(selection, tmp_path):
    errors = selection.verify_sources(tmp_path)
    assert len(errors) == 52 and all(e["status"] == "asset_unavailable" for e in errors)
    path = tmp_path / selection.catalog["files"][0]["filename"]
    path.parent.mkdir(parents=True)
    path.write_text("not the cited source")
    errors = selection.verify_sources(tmp_path)
    assert any(e["status"] == "source_incompatible" for e in errors)
    system = from_smiles("C[Si](C)(C)O[Si](C)(C)C")
    result = close_domain(system, selection, source_errors=errors)
    assert result.payload["status"] == "source_incompatible"
    with pytest.raises(DomainClosureBlockedError):
        result.require_executable()


def test_immutable_record_and_payload_isolation(selection):
    result = close_domain(from_smiles("[SiH3]O[SiH3]"), selection)
    with pytest.raises(FrozenInstanceError):
        result.json_text = "forged"
    result.payload["missing_rows"].clear()
    assert result.payload["missing_rows"]


@pytest.mark.parametrize(
    "status",
    [
        "parameter_incomplete",
        "typing_ambiguous",
        "asset_unavailable",
        "source_incompatible",
    ],
)
def test_fresh_process_reconstructs_each_blocked_status(selection, status):
    system = from_smiles(
        "[SiH3]O[SiH3]" if status == "typing_ambiguous" else "C[Si](C)(C)O[Si](C)(C)C"
    )
    errors = (
        []
        if status in ("parameter_incomplete", "typing_ambiguous")
        else [{"status": status, "reason": "explicit source diagnostic"}]
    )
    result = close_domain(system, selection, source_errors=errors)
    code = "from island.forcefields.literature_closure import DomainClosure,DomainClosureBlockedError; import sys; r=DomainClosure(sys.stdin.read());\ntry: r.require_executable()\nexcept DomainClosureBlockedError as e: print(e.status)"
    output = subprocess.run(
        [sys.executable, "-c", code],
        input=result.json_text,
        text=True,
        capture_output=True,
        check=True,
    )
    assert output.stdout.strip() == status


def test_relocation_fresh_process_and_coherent_graph_mutation(tmp_path):
    study = tmp_path / "study"
    missing = tmp_path / "missing-references"
    report = execute(study, missing, SELECTION, CATALOG)
    assert len(report["cases"]) == 6
    moved = tmp_path / "relocated"
    shutil.copytree(study, moved)
    expected = reconstruct(moved)
    process = subprocess.run(
        [sys.executable, "scripts/phase_4j45.py", "--reconstruct", str(moved)],
        text=True,
        capture_output=True,
        check=True,
    )
    assert json.loads(process.stdout) == expected
    assert {
        str(p.relative_to(study)): sha256(p.read_bytes()).hexdigest()
        for p in study.rglob("*")
        if p.is_file()
    } == {
        str(p.relative_to(moved)): sha256(p.read_bytes()).hexdigest()
        for p in moved.rglob("*")
        if p.is_file()
    }
    name = report["cases"][0]["name"]
    data = read(moved / name / "system.json")
    data["sites"][0]["metadata"]["chain_id"] = "forged"
    write(moved / name / "system.json", data)
    with pytest.raises(ValueError):
        reconstruct(moved)
    assert not list(study.rglob("*.data"))
    assert not list(study.rglob("*.chk"))
    assert not list(study.rglob("*.dcd"))
