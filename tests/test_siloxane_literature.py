"""Real paper rows and fail-closed graph audits; no synthetic force-field values."""

import json
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import FrozenInstanceError
from hashlib import sha256
from pathlib import Path

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields.external_pdms import digest
from island.forcefields.siloxane_literature import (
    CATALOG_IDENTITY,
    SiloxaneParameterAssessment,
    SiloxaneParameterError,
    SiloxaneParametersBlockedError,
    SiloxaneParameterSupplement,
    assess_siloxane_parameters,
    load_siloxane_parameters,
)
from scripts.assess_sun_siloxane_parameters import fixtures, run

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "docs/evidence/sun_siloxane_parameters/parameters.json"
PIN = "d23244f11d6be8fb683ab01b1f586ee466f60c3b65d10b8da9379eda111961ca"
REFERENCES = ROOT / "island_ffref"


@pytest.fixture
def supplement():
    return SiloxaneParameterSupplement(pack(json.loads(CATALOG.read_text())))


@pytest.fixture
def local_sources(supplement):
    if not all(
        (REFERENCES / s["file"]).is_file()
        for s in supplement.payload["sources"].values()
    ):
        pytest.skip("explicit user-supplied local source PDFs are unavailable")
    return REFERENCES


def test_real_paper_rows_and_printed_blanks(supplement):
    p = supplement.payload
    assert sha256(CATALOG.read_bytes()).hexdigest() == PIN
    assert digest(p) == CATALOG_IDENTITY
    assert Counter(r["source"] for r in p["rows"]) == {"sun1995": 77, "sun1997": 28}
    partial = [r for r in p["rows"] if r["missing_fields"]]
    assert len(partial) == 8
    assert all(r["missing_fields"] == ["k4"] for r in partial)
    assert all(r["original_values"][-1] is None for r in partial)
    assert any(r["explicit_zero_fields"] for r in p["rows"])
    for row in p["rows"]:
        assert row["evidence_checksum"] == digest(
            {k: v for k, v in row.items() if k != "evidence_checksum"}
        )
        assert row["conversion"]["converted_values"] == row["original_values"]
        assert row["execution_authorized"] is False


@pytest.mark.parametrize(
    "family,labels,values",
    [
        ("quartic_bond", "o2-si4", ["1.640", "350.1232", "-517.3424", "673.7067"]),
        ("quartic_angle", "si4-o2-si4", ["159.0", "8.5000", "-13.4188", "-4.1785"]),
        ("bond_increment", "o2-si4", ["-0.2225"]),
        ("nonbond(9-6)", "si4c", ["4.290", "0.131"]),
        ("bond-bond", "si4-o2-si4", ["41.1143"]),
        ("bond-angle", "si4-o2-si4", ["28.6686", "28.6686"]),
    ],
)
def test_exact_siloxane_rows(supplement, family, labels, values):
    (row,) = supplement.lookup("sun1997", family, labels)
    assert row["original_values"] == values
    assert row["source_evidence"]["pdf_page"] in (21, 22)


def test_source_bound_literal_lookup_and_payload_ownership(supplement):
    assert not supplement.lookup("sun1997", "quartic_bond", "si4-o2")
    assert not supplement.lookup("sun1997", "quartic_bond", "Si-C")
    assert supplement.lookup("sun1995", "torsion_3", "*-Si-C-*")
    assert not supplement.lookup("sun1995", "torsion_3", "H-Si-C-H")
    assert supplement.lookup("sun1995", "nonbond(9-6)", "Si")[0]["original_values"] == [
        "4.450",
        "0.190",
    ]
    assert supplement.lookup("sun1997", "nonbond(9-6)", "si4")[0][
        "original_values"
    ] == ["4.405", "0.198"]
    supplement.payload["rows"].clear()
    assert len(supplement.payload["rows"]) == 105
    with pytest.raises(FrozenInstanceError):
        supplement.json_text = "changed"


@pytest.mark.parametrize(
    "mutation", ["zero", "units", "source", "charges", "authority", "profile", "rows"]
)
def test_resigned_catalog_changes_rejected(supplement, mutation):
    p = supplement.payload
    if mutation == "zero":
        next(r for r in p["rows"] if r["missing_fields"])["original_values"][-1] = "0.0"
    elif mutation == "units":
        p["rows"][0]["source_units"] = ["unknown"]
    elif mutation == "source":
        p["sources"]["sun1997"]["sha256"] = "a" * 64
    elif mutation == "charges":
        p["exact_charge_vector"] = [0.0]
    elif mutation == "authority":
        p["execution_authorized"] = True
    elif mutation == "profile":
        p["execution_profiles"] = ["native_pcff"]
    else:
        p["rows"].pop()
    p["checksum"] = digest({k: v for k, v in p.items() if k != "checksum"})
    with pytest.raises(
        SiloxaneParameterError, match="parameter_rows_or_authority_changed"
    ):
        SiloxaneParameterSupplement(pack(p))


def test_missing_and_changed_assets_classified(supplement, local_sources, tmp_path):
    with pytest.raises(SiloxaneParameterError) as exc:
        load_siloxane_parameters(
            tmp_path / "absent.json", PIN, references=local_sources
        )
    assert exc.value.status == "asset_unavailable"
    with pytest.raises(SiloxaneParameterError) as exc:
        load_siloxane_parameters(CATALOG, "a" * 64, references=local_sources)
    assert exc.value.status == "source_incompatible"
    with pytest.raises(SiloxaneParameterError) as exc:
        supplement.verify_sources(tmp_path)
    assert exc.value.status == "asset_unavailable"
    source = next(iter(supplement.payload["sources"].values()))
    changed = tmp_path / source["file"]
    changed.parent.mkdir(parents=True)
    changed.write_text("not the source paper")
    with pytest.raises(SiloxaneParameterError) as exc:
        supplement.verify_sources(tmp_path)
    assert exc.value.status == "source_incompatible"


def test_loader_requires_numeric_evidence(local_sources, tmp_path):
    relocated = tmp_path / "parameters.json"
    shutil.copyfile(CATALOG, relocated)
    with pytest.raises(SiloxaneParameterError) as exc:
        load_siloxane_parameters(relocated, PIN, references=local_sources)
    assert exc.value.status == "asset_unavailable"
    (tmp_path / "source-rows.tsv").write_text("invented values")
    with pytest.raises(SiloxaneParameterError) as exc:
        load_siloxane_parameters(relocated, PIN, references=local_sources)
    assert exc.value.status == "source_incompatible"


def test_graph_has_real_rows_and_exact_unresolved_dependencies(
    supplement, local_sources
):
    system = from_smiles("C[Si](C)(C)O[Si](C)(C)C")
    result = assess_siloxane_parameters(system, supplement, references=local_sources)
    result.validate_integrity(system, supplement, references=local_sources)
    p = result.payload
    assert p["status"] == "parameter_incomplete"
    assert not p["unsupported_sites"]
    associations = p["row_associations"]
    assert all(
        not r["prior_silane_candidates"]
        for key, r in associations.items()
        if key.startswith("nonbonded:")
    )
    assert any(r["comparison_only_source_rows"] for r in associations.values())
    assert any(
        "sun1997:quartic_bond:o2-si4" in r["primary_source_rows"]
        for r in associations.values()
    )
    assert any(
        "sun1997:nonbond(9-6):si4c" in r["primary_source_rows"]
        for r in associations.values()
    )
    assert any(
        "sun1995:quartic_bond:Si-C" in r["prior_silane_candidates"]
        for r in associations.values()
    )
    assert any(
        "sun1995:quartic_angle:Si-C-H" in r["partial_source_rows"]
        for r in associations.values()
    )
    assert any(
        k.startswith("quartic_bond:")
        and r["labels"] == ["c4", "h1"]
        and r["state"] == "source_row_absent"
        for k, r in associations.items()
    )
    assert set(p["missing_dependencies"]) == set(associations)
    assert all(r["selected_row"] is None for r in associations.values())
    assert (
        p["parameter_assignment"]
        is p["profile_identity"]
        is p["exact_charge_vector"]
        is None
    )
    with pytest.raises(SiloxaneParametersBlockedError) as exc:
        result.require_executable()
    assert exc.value.result is result
    assert exc.value.diagnostics == p["missing_dependencies"]


@pytest.mark.parametrize(
    "smiles",
    ["[SiH2]O[SiH3]", "[O-][Si](C)(C)C", "CO[Si](C)(C)C", "COP(=O)(O)OC", "C[SiH3]"],
)
def test_unsupported_chemical_environments_remain_classified(
    supplement, local_sources, smiles
):
    result = assess_siloxane_parameters(
        from_smiles(smiles), supplement, references=local_sources
    )
    assert result.payload["status"] == "typing_ambiguous"
    assert result.payload["typing_diagnostics"]
    with pytest.raises(SiloxaneParametersBlockedError):
        result.require_executable()


@pytest.mark.parametrize(
    "mutation",
    ["status", "inventory", "order", "component", "graph", "charge", "profile"],
)
def test_assessment_mutations_cannot_authorize_or_rebind(
    supplement, local_sources, mutation
):
    system = from_smiles("[SiH3]O[SiH3]")
    result = assess_siloxane_parameters(system, supplement, references=local_sources)
    p = result.payload
    if mutation == "status":
        p["status"] = "prepared"
    elif mutation == "inventory":
        p["missing_dependencies"].clear()
    elif mutation == "order":
        p["requirements"]["atom_order"].reverse()
    elif mutation == "component":
        p["requirements"]["graph"]["component_membership"]["0"] = 7
    elif mutation == "graph":
        p["requirements"]["graph"]["bonds"].pop()
    elif mutation == "charge":
        p["exact_charge_vector"] = [0] * system.number_of_sites
    else:
        p["profile_identity"] = "forged"
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(SiloxaneParameterError):
        mutated = SiloxaneParameterAssessment(pack(p))
        mutated.validate_integrity(system, supplement, references=local_sources)


def test_stale_graph_and_molecule_identity(supplement, local_sources):
    system = from_smiles("[SiH3]O[SiH3]")
    result = assess_siloxane_parameters(system, supplement, references=local_sources)
    next(iter(system.topology.sites.values())).metadata["chain_id"] = "stale"
    with pytest.raises(SiloxaneParameterError, match="stale_or_mutated"):
        result.validate_integrity(system, supplement, references=local_sources)


def test_relocated_sources_and_fresh_process(supplement, local_sources, tmp_path):
    refs = tmp_path / "references"
    for source in supplement.payload["sources"].values():
        target = refs / source["file"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(local_sources / source["file"], target)
    catalog = tmp_path / "catalog/parameters.json"
    catalog.parent.mkdir()
    shutil.copyfile(CATALOG, catalog)
    shutil.copyfile(
        CATALOG.parent / "source-rows.tsv", catalog.parent / "source-rows.tsv"
    )
    output = tmp_path / "study"
    report = run(catalog, PIN, refs, output)
    assert len(report["fixtures"]) == len(list(fixtures())) == 6
    assert not report["execution_profiles"] and not report["backend_executions"]
    assert report["acceptance_gate"] == report["native_pcff_acceptance_gate"] == 1
    command = [
        sys.executable,
        str(ROOT / "scripts/assess_sun_siloxane_parameters.py"),
        "--catalog",
        str(catalog),
        "--sha256",
        PIN,
        "--references",
        str(refs),
        "--output",
        str(output),
        "--reconstruct",
        "--require-executable",
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    assert completed.returncode == 1, completed.stderr
    assert json.loads(completed.stdout) == report
    saved = output / "report.json"
    altered = json.loads(saved.read_text())
    altered["acceptance_gate"] = 0
    saved.write_text(json.dumps(altered))
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    assert completed.returncode == 1
    assert json.loads(completed.stdout)["status"] == "source_incompatible"
