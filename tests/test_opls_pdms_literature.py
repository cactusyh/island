"""Data staging must never confer OPLS typing, charge or execution authority."""

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
from island.forcefields.literature_closure import (
    DomainClosureBlockedError,
    DomainSelection,
)
from island.forcefields.oplsaa.pdms_literature import (
    SOURCE_FILES,
    OPLSPDMSLiteratureError,
    OPLSPDMSLiteratureSupplement,
    build_opls_pdms_literature,
    load_opls_pdms_literature,
)
from island.forcefields.preparation import ForceFieldRequest, ForceFieldRequestError

ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / "docs/evidence/phase_4j45/selection.json"
CATALOG = ROOT / "docs/evidence/phase_4j44/source-audit.json"


@pytest.fixture
def supplement():
    return build_opls_pdms_literature(
        DomainSelection(SELECTION.read_text(), CATALOG.read_text())
    )


@pytest.fixture
def references(tmp_path):
    root = tmp_path / "references"
    for name in SOURCE_FILES:
        original = ROOT / "island_ffref" / name
        if not original.is_file():
            pytest.skip("Requires the two explicitly supplied local Xiang PDFs")
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, target)
    return root


def test_source_rows_preserved_without_physical_assignment(supplement):
    p = supplement.payload
    assert Counter(row["kind"] for row in p["rows"]) == {
        "nonbonded": 4,
        "bond": 3,
        "angle": 6,
        "proper_torsion": 4,
        "partial_charge": 4,
    }
    original = supplement.selection.payload["candidate_rows"]
    assert [row["source_evidence"] for row in p["rows"][:17]] == original
    for row in p["rows"]:
        tokens = row["source_evidence"]["raw_text"].split()
        assert row["original_values"] == [
            tokens[i] for i in row["numeric_columns_zero_based"]
        ]
        assert row["source_evidence"]["page"] == 6
        assert row["source_units"] is row["converted_values"] is None
        assert row["accepted_assignment"] is None
        assert row["evidence_checksum"] == digest(
            {k: v for k, v in row.items() if k != "evidence_checksum"}
        )
    charges = supplement.lookup("partial_charge", ["Si"])
    assert charges[0]["original_values"] == ["0.9", "0.7"]
    assert charges[0]["column_names"] == ["inside monomer", "end monomer"]
    assert charges[0]["charge_method_assignment"] is None
    assert p["exact_charge_vector"] is p["profile_identity"] is None
    assert not p["execution_authorized"] and p["acceptance_gate"] == 1


def test_literal_lookup_does_not_type_reverse_or_borrow(supplement):
    assert supplement.lookup("bond", ["o", "si"])[0]["original_values"] == [
        "374",
        "1.64",
    ]
    for query in (["si", "o"], ["O", "Si"], ["*", "si"], ["si", "h"]):
        assert supplement.lookup("bond", query) == ()
    assert supplement.lookup("class2_cross_term", ["c", "si", "o"]) == ()


def test_explicit_source_zeros_are_preserved_only_for_their_rows(supplement):
    row = supplement.lookup("proper_torsion", ["si", "o", "si", "o"])[0]
    assert row["original_values"] == ["0", "0", "0", "0"]
    assert row["explicit_zero_value_indices"] == [0, 1, 2, 3]
    assert supplement.lookup("proper_torsion", ["h", "si", "o", "si"]) == ()
    assert row["converted_values"] is None


def test_immutable_and_payload_lookup_isolation(supplement):
    with pytest.raises(FrozenInstanceError):
        supplement.json_text = "forged"
    supplement.payload["rows"].clear()
    supplement.lookup("bond", ["o", "si"])[0]["original_values"].clear()
    assert len(supplement.payload["rows"]) == 21
    assert supplement.lookup("bond", ["o", "si"])[0]["original_values"]


@pytest.mark.parametrize(
    "mutation",
    [
        "value",
        "unit",
        "charge",
        "source",
        "family",
        "zero",
        "missing_row",
        "duplicate",
        "profile",
        "authority",
        "limits",
        "type",
        "label",
        "column",
        "evidence",
        "mixing",
        "charge_assignment",
        "converted",
        "status",
        "identity",
    ],
)
def test_resigning_cannot_change_data_or_infer_missing_authority(supplement, mutation):
    p = supplement.payload
    first = p["rows"][0]
    if mutation == "value":
        first["original_values"][0] = "0.15"
    elif mutation == "unit":
        first["source_units"] = "kcal/mol"
    elif mutation == "charge":
        p["rows"][-1]["original_values"][0] = "0.06"
    elif mutation == "source":
        first["source_evidence"]["source_sha256"] = "f" * 64
    elif mutation == "family":
        p["source_family"] = "native:oplsaa"
    elif mutation == "zero":
        nonzero_torsion = next(
            row
            for row in p["rows"]
            if row["kind"] == "proper_torsion" and "0.18" in row["original_values"]
        )
        nonzero_torsion["original_values"] = ["0"] * 4
    elif mutation == "missing_row":
        p["rows"].pop()
    elif mutation == "duplicate":
        p["rows"].append(dict(first))
    elif mutation == "profile":
        p["profile_identity"] = "executable"
    elif mutation == "authority":
        p["execution_authorized"] = True
    elif mutation == "limits":
        p["limitations"]["production_validated"] = True
    elif mutation == "type":
        p["execution_authorized"] = 0
    elif mutation == "label":
        first["literal_labels"] = ["Si", "Si"]
    elif mutation == "column":
        first["numeric_columns_zero_based"].reverse()
    elif mutation == "evidence":
        first["source_evidence"]["raw_text"] = "pair_coeff O O 0 0"
    elif mutation == "mixing":
        p["unresolved_policies"]["mixing"] = "geometric assumed"
    elif mutation == "charge_assignment":
        p["exact_charge_vector"] = [0.0]
    elif mutation == "converted":
        first["converted_values"] = ["0.58576", "0.29"]
    elif mutation == "status":
        p["status"] = "prepared"
    else:
        p["identity"] = "f" * 64
    if mutation != "identity":
        p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        OPLSPDMSLiteratureSupplement(pack(p), supplement.selection)
    assert error.value.status == "source_incompatible"


@pytest.mark.parametrize("pin", [None, "0" * 64, "invalid"])
def test_missing_wrong_independent_pin_is_classified(supplement, tmp_path, pin):
    path = tmp_path / "parameters.json"
    path.write_text(supplement.json_text)
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        load_opls_pdms_literature(
            path, pin, selection=supplement.selection, references=tmp_path
        )
    assert error.value.status == "source_incompatible"


def test_malformed_export_and_missing_sources_classified(supplement, tmp_path):
    path = tmp_path / "parameters.json"
    path.write_text('{"invalid":')
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        load_opls_pdms_literature(
            path,
            sha256(path.read_bytes()).hexdigest(),
            selection=supplement.selection,
            references=tmp_path,
        )
    assert error.value.status == "source_incompatible"
    path.write_text(supplement.json_text)
    with pytest.raises(OPLSPDMSLiteratureError) as error:
        load_opls_pdms_literature(
            path,
            sha256(path.read_bytes()).hexdigest(),
            selection=supplement.selection,
            references=tmp_path,
        )
    assert error.value.status == "asset_unavailable"
    assert len(error.value.diagnostics) == 2


def test_native_request_rejects_a_literature_record(supplement):
    with pytest.raises(ForceFieldRequestError):
        ForceFieldRequest("oplsaa", supplement)


@pytest.mark.parametrize("kind", ["missing", "mutated"])
def test_cli_missing_or_mutated_historical_input_is_classified(tmp_path, kind):
    path = tmp_path / "selection.json"
    if kind == "mutated":
        path.write_text("{}")
    output = tmp_path / "output"
    run = subprocess.run(
        [
            sys.executable,
            "scripts/supplement_opls_pdms_literature.py",
            "--selection",
            str(path),
            "--catalog",
            str(CATALOG),
            "--references",
            str(tmp_path),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert run.returncode == 1 and not run.stderr
    assert json.loads(run.stdout)["status"] == (
        "asset_unavailable" if kind == "missing" else "source_incompatible"
    )
    assert not output.exists()


def test_graph_assessment_preserves_boundary_without_numerical_engine(
    supplement, references, monkeypatch
):
    def forbidden(*args, **kwargs):
        pytest.fail("Literature lookup must not prepare or instantiate an engine")

    monkeypatch.setattr(
        "island.forcefields.oplsaa.adapter.check_foyer_installation", forbidden
    )
    monkeypatch.setattr("island.forcefields.preparation.prepare_forcefield", forbidden)
    system = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
    result = supplement.assess(system, references=references)
    assert result.payload["status"] == "parameter_incomplete"
    assert result.payload["missing_rows"]
    with pytest.raises(DomainClosureBlockedError) as error:
        result.require_executable()
    assert error.value.result is result
    other = from_smiles("CO[Si](C)(C)C")
    assert supplement.assess(other, references=references).payload["status"] == (
        "typing_ambiguous"
    )
    with pytest.raises(ValueError, match="stale_graph"):
        result.validate_integrity(other, supplement.selection)
    (references / SOURCE_FILES[0]).write_bytes(b"mutated paper")
    result = supplement.assess(system, references=references)
    assert result.payload["status"] == "source_incompatible"
    with pytest.raises(DomainClosureBlockedError):
        result.require_executable()


def test_relocated_export_reconstructs_in_a_fresh_process(
    supplement, references, tmp_path
):
    relocated = tmp_path / "moved"
    relocated.mkdir()
    path = relocated / "parameters.json"
    path.write_text(supplement.json_text)
    selection_path, catalog_path = (
        relocated / "selection.json",
        relocated / "catalog.json",
    )
    selection_path.write_text(supplement.selection.json_text)
    catalog_path.write_text(supplement.selection.catalog_text)
    relocated_references = relocated / "references"
    shutil.move(str(references), relocated_references)
    pin = sha256(path.read_bytes()).hexdigest()
    loaded = load_opls_pdms_literature(
        path, pin, selection=supplement.selection, references=relocated_references
    )
    assert loaded == supplement
    args = [
        sys.executable,
        "scripts/supplement_opls_pdms_literature.py",
        "--selection",
        str(selection_path),
        "--catalog",
        str(catalog_path),
        "--references",
        str(relocated_references),
        "--reconstruct",
        str(path),
        "--sha256",
        pin,
    ]
    run = subprocess.run(args, text=True, capture_output=True, check=True)
    report = json.loads(run.stdout)
    assert report["identity"] == supplement.identity
    assert report["row_counts"]["partial_charge"] == 4
    assert len(report["cases"]) == 6
    assert report["acceptance_gate"] == 1
    assert report["execution_profiles"] == report["backend_executions"] == []
    blocked = subprocess.run(
        args + ["--require-executable"], text=True, capture_output=True, check=False
    )
    assert blocked.returncode == 1 and json.loads(blocked.stdout) == report
    assert not any(p.suffix in (".data", ".chk", ".dcd") for p in relocated.rglob("*"))
