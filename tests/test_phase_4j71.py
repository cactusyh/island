"""Literature controls are separate from retained QM evidence; no new QM/MD."""

import copy
import json
import os
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from island.pcff_completion import CompletionError, canonical, signed
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.literature_intake import (
    LiteratureRecord,
    candidate_sites,
    coherent_bundle_gate,
    compatibility,
    converter_gate,
    make_record,
    normalize,
    record_object,
    semantic_declarations,
)
from scripts.extract_phase_4j71 import (
    BASE,
    BRANCH,
    DECLARATIONS,
    REQUIRED,
    ROOT,
    inputs,
    reports,
)
from scripts.replay_phase_4j71 import replay
from scripts.validate_phase_4j71_receipt import semantics

REPO = Path(__file__).resolve().parents[1]
E = REPO / ROOT


def resign(value):
    return signed({k: v for k, v in value.items() if k != "identity"})


@pytest.fixture(scope="module")
def retained():
    _, cells, reading, sources = inputs(REPO, E)
    return cells, reading, sources, reports(REPO, E)


def test_actual_inventory_and_complete_reading_coverage(retained):
    reading = retained[1]
    assert len(reading["files"]) == 25
    original = read(E / "source-inventory.json")
    assert {r["sha256"] for r in original["files"]} == {
        r["sha256"] for r in reading["files"]
    }
    bad = [
        r
        for r in reading["files"]
        if r.get("reading_status") == "unreadable_truncated_checkpoint"
    ]
    assert len(bad) == 1 and not bad[0]["actually_read_successfully"]
    opls = next(
        r for r in reading["files"] if r["relative_path"].endswith("/ja11225.pdf")
    )
    assert opls["page_coverage"] == list(range(1, 14))
    assert opls["full_text_kind"] == "supporting_information_only_main_absent"
    assert opls["read_mode"] == "manual_rendered_no_OCR"
    assert reading["OCR_used"] is False
    for r in reading["files"]:
        if r["actually_read_successfully"]:
            assert r["page_coverage"] == list(range(1, r["page_count"] + 1))
            assert len(r["page_content_sha256"]) == r["page_count"]


@pytest.mark.parametrize(
    "field",
    [
        "source_file_sha256",
        "citation",
        "location",
        "raw_value",
        "raw_unit",
        "scalar_energy_equation",
        "analytic_force_equation",
        "equilibrium_values",
        "phase_convention",
        "multiplicity",
        "orientation_reversal_rule",
        "cross_term_dependencies",
        "provenance_category",
        "source_force_field_family",
        "source_model_resolution",
        "extraction_confidence",
        "execution_eligible",
        "normalized",
        "ordered_sites",
    ],
)
def test_resigned_parameter_mutation_rejected(retained, field):
    cell = retained[0][0]
    source = retained[2][cell["source"]]
    r = make_record(cell, source)
    r[field] = "forged"
    with pytest.raises(CompletionError, match="mutation"):
        LiteratureRecord(canonical(resign(r))).validate(cell, source)


@pytest.mark.parametrize(
    "field",
    [
        "pdf_page",
        "table",
        "equation",
        "row_label",
        "cell_coordinates",
        "bbox_pdf_points",
    ],
)
def test_resigned_citation_table_cell_location_mutation(retained, field):
    cell = retained[0][0]
    source = retained[2][cell["source"]]
    r = copy.deepcopy(make_record(cell, source))
    r["location"][field] = "forged"
    with pytest.raises(CompletionError):
        LiteratureRecord(canonical(resign(r))).validate(cell, source)


@pytest.mark.parametrize(
    "name",
    [
        "source-inventory.json",
        "reading-status.json",
        "extraction-ledger.json",
        "archive-members.json",
        "contract.json",
    ],
)
def test_review_ledger_hash_resigning_cannot_change_authority(tmp_path, name):
    for path in E.iterdir():
        if path.is_file():
            shutil.copyfile(path, tmp_path / path.name)
    p = read(tmp_path / name)
    p["forged_replacement"] = True
    (tmp_path / name).write_text(json.dumps(resign(p)))
    with pytest.raises(CompletionError, match="mutation"):
        inputs(REPO, tmp_path)


def test_immutable_record_and_payload_copy(retained):
    c = retained[0][0]
    obj = record_object(c, retained[2][c["source"]])
    p = obj.payload
    p["raw_value"][0] = "1000"
    assert obj.payload["raw_value"][0] != "1000"
    with pytest.raises(FrozenInstanceError):
        obj.json_text = "forged"


def test_original_file_hash_checker_rejects_changed_source(monkeypatch):
    import scripts.read_phase_4j71_sources as checker

    original_pin = checker.pin
    original_hashes = {
        str(REPO / r["relative_path"]): r["sha256"]
        for r in read(E / "reading-status.json")["files"]
    }
    monkeypatch.setattr(
        checker,
        "pin",
        lambda p: (
            "f" * 64
            if str(p).endswith("/ja11225.pdf")
            else (
                original_hashes[str(p)]
                if str(p) in original_hashes
                else original_pin(p)
            )
        ),
    )
    with pytest.raises(CompletionError, match="source_hash_mutation"):
        checker.verify_originals(REPO, E)


def test_printed_zeros_blanks_and_no_inferred_decimal(retained):
    records = retained[3]["parameter-records.json"]["records"]
    blank = next(r for r in records if None in r["raw_value"])
    for v, n in zip(blank["raw_value"], blank["normalized"], strict=True):
        if v is None:
            assert n["value"] is n["factor"] is None
    zero = next(r for r in records if r["explicit_printed_zero_fields"])
    assert not zero["zero_filled"]
    unresolved = next(
        r for r in records if r["extraction_issue"] and "0.945" in r["extraction_issue"]
    )
    assert unresolved["raw_value"] == [".945"]
    assert unresolved["extraction_confidence"] == "extraction_unresolved"
    assert not unresolved["execution_eligible"]


def test_normalization_preserves_raw_physical_factors():
    assert normalize(".1640", "nm")["value"] == "1.6400"
    assert float(normalize("4.184", "kJ/mol")["value"]) == 1.0
    assert float(normalize("180", "degree")["value"]) == pytest.approx(
        3.141592653589793
    )
    assert normalize(None, "kcal/mol")["value"] is None


def test_actual_primary_numeric_values_and_source_typing(retained):
    records = retained[3]["parameter-records.json"]["records"]
    sun95 = next(
        x
        for x in records
        if x["source_file_sha256"].startswith("9c6d6315b91e")
        and x["interaction_family"] == "quartic_bond"
        and x["ordered_sites"] == ["Si", "C"]
    )
    assert sun95["raw_value"] == ["1.899", "189.65", "-279.42", "307.51"]
    sun97 = next(
        x
        for x in records
        if x["source_file_sha256"].startswith("65bda3a5b14c")
        and x["interaction_family"] == "nonbond(9-6)"
        and x["ordered_sites"] == ["si4"]
    )
    assert sun97["raw_value"] == ["4.405", "0.198"]
    assert not sun97["native_pcff_authority"]
    assert sun97["provenance_category"] == "literature_external"


@pytest.mark.parametrize(
    "family,expected",
    [
        ("PCFF", "exact_literature_external_match"),
        ("COMPASS", "compatible_but_non_native"),
        ("OPLS", "different_force_function"),
        ("GAFF", "different_force_function"),
    ],
)
def test_manufactured_family_separation_same_element_pattern(family, expected):
    # Synthetic source-contract control, not a discovered PCFF row.
    r = {
        "interaction_family": "angle-angle",
        "source_model_resolution": "all_atom",
        "source_force_field_family": family,
        "missing_convention_fields": [],
        "extraction_confidence": "direct_printed_numeric_cells",
    }
    t = {"family": "angle-angle", "ordered_types": ["hc", "c3", "sio", "hc"]}
    assert (
        compatibility(
            r, t, ["H", "C", "Si", "H"], same_environment=True, complete_contract=True
        )
        == expected
    )
    r["source_model_resolution"] = "united_atom"
    assert (
        compatibility(
            r, t, ["H", "C", "Si", "H"], same_environment=True, complete_contract=True
        )
        == "different_resolution"
    )


def test_missing_convention_and_ambiguous_ocr_rejected():
    r = {
        "interaction_family": "angle-angle",
        "source_model_resolution": "all_atom",
        "source_force_field_family": "PCFF",
        "missing_convention_fields": ["equilibrium"],
        "extraction_confidence": "direct_printed_numeric_cells",
    }
    t = {"family": "angle-angle", "ordered_types": ["hc", "c3", "sio", "hc"]}
    assert (
        compatibility(
            r, t, ["H", "C", "Si", "H"], same_environment=True, complete_contract=True
        )
        == "ambiguous"
    )
    r.update(
        missing_convention_fields=[], extraction_confidence="extraction_unresolved"
    )
    assert (
        compatibility(
            r, t, ["H", "C", "Si", "H"], same_environment=True, complete_contract=True
        )
        == "ambiguous"
    )


def test_actual_conflicts_and_equal_duplicates_preserved(retained):
    c = retained[3]["source-conflicts.json"]
    assert len(c["unresolved_conflicts"]) == 2
    for conflict in c["unresolved_conflicts"]:
        a, b = conflict["records"]
        assert a["raw_value"] != b["raw_value"]
        assert a["source_file_sha256"] != b["source_file_sha256"]
    assert len(c["duplicate_sources"]) == 16
    assert not c["values_averaged"] and not c["preferred_values_selected"]


def test_all_378_stable_instances_and_no_exact_assignments(retained):
    m = retained[3]["mapping-report.json"]
    parent = read(REPO / "docs/evidence/phase_4j70/completeness-matrix.json")
    assert [r["id"] for r in m["instances"]] == [r["id"] for r in parent["instances"]]
    assert (
        m["graph_identity"] == parent["graph_identity"]
        and m["atom_order"] == parent["atom_order"]
    )
    assert len(m["instances"]) == 378
    assert all(not r["resolved"] and r["exact_matches"] == [] for r in m["instances"])
    assert m["family_counts"]["angle-angle"] == 78
    assert set(m["candidate_status_counts"]) == {
        "ambiguous",
        "compatible_but_non_native",
    }


def test_source_topology_order_and_equilibrium_reversal_not_silently_fixed(retained):
    cells = retained[0]
    c = next(
        c
        for c in cells
        if c.get("model") == "COMPASS"
        and c["source"] == "8b0e8b74303b"
        and c.get("command") == "angle_coeff"
        and c["location"]["row_label"] == "102"
        and c.get("substyle") is None
    )
    assert candidate_sites(c) == ["O", "Si", "C"]
    aat = next(
        c
        for c in cells
        if c.get("model") == "COMPASS"
        and c["source"] == "8b0e8b74303b"
        and c.get("substyle") == "aat"
        and c["location"]["row_label"] == "24"
    )
    assert aat["raw_values"][1:] == ["112.0355", "114.9060"]


def test_partial_bundle_converter_and_static_closure_not_executed(retained):
    r = retained[3]
    bundle = r["partial-literature-bundle.json"]
    assert (
        len(bundle["remaining_row_ids"]) == 378
        and bundle["external_executable_rows"] == []
    )
    with pytest.raises(CompletionError, match="partial_bundle"):
        converter_gate(bundle)
    assert not r["strict-converter-report.json"]["executed"]
    assert not r["static-closure-report.json"]["executed"]
    assert r["strict-converter-report.json"]["attempts"] == 0


@pytest.mark.parametrize(
    "mutation", ["zero", "wildcard", "duplicate_record", "duplicate_target"]
)
def test_manufactured_zero_wildcard_duplicate_rejection(mutation):
    matrix = [{"id": "one"}]
    records = [{"record_id": "a", "zero_filled": False, "wildcard_resolved": False}]
    if mutation == "zero":
        records[0]["zero_filled"] = True
    elif mutation == "wildcard":
        records[0]["wildcard_resolved"] = True
    elif mutation == "duplicate_record":
        records += copy.deepcopy(records)
    else:
        matrix += copy.deepcopy(matrix)
    with pytest.raises(CompletionError):
        coherent_bundle_gate(matrix, [], records)


def test_manufactured_complete_label_without_source_evidence_cannot_activate_converter():
    with pytest.raises(CompletionError, match="source_mapping_evidence"):
        converter_gate(
            {
                "parameter_completeness": "complete",
                "coherent_contract": True,
                "resolved_instances": 378,
                "remaining_row_ids": [],
            }
        )
    with pytest.raises(CompletionError, match="partial_bundle"):
        converter_gate(
            {
                "parameter_completeness": "complete",
                "coherent_contract": True,
                "resolved_instances": 377,
                "remaining_row_ids": ["missing"],
            }
        )


@pytest.mark.parametrize("field", list(DECLARATIONS))
def test_resigned_scientific_metadata_rejected(field):
    r = signed(DECLARATIONS)
    v = r[field]
    r[field] = not v if type(v) is bool else v + 1 if type(v) is int else "forged"
    with pytest.raises(CompletionError, match="declaration"):
        semantic_declarations(resign(r))


def make_receipt():
    r = signed(
        {
            **DECLARATIONS,
            "schema": "island.j71-receipt.v1",
            "base_commit": BASE,
            "required_J70": REQUIRED,
            "branch": BRANCH,
            "manifest_sha256": None,
            "software_counts": {},
            "terminal_unique_tests": 0,
            "historical_files_preserved": 0,
            "offline_reconstruction": replay(REPO, E),
        }
    )
    m = signed(
        {
            "schema": "island.j71-evidence-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {},
            "historical_files": {},
            "artifacts": {},
            "scientific_scope": DECLARATIONS,
        }
    )
    return r, m


@pytest.mark.parametrize("field", list(DECLARATIONS))
def test_resigned_manifest_scientific_declaration_rejected(field):
    r, m = make_receipt()
    m = copy.deepcopy(m)
    m["scientific_scope"][field] = "forged"
    with pytest.raises(CompletionError, match="scientific_semantics"):
        semantics(r, resign(m))


@pytest.mark.parametrize(
    "kind",
    [
        "valid",
        "receipt_schema",
        "receipt_extra",
        "manifest_schema",
        "manifest_extra",
        "ancestry",
        "nested_boolean",
        "nested_count",
    ],
)
def test_receipt_schema_and_nested_semantics(kind):
    r, m = make_receipt()
    if kind == "receipt_schema":
        r["schema"] = "wrong"
    elif kind == "receipt_extra":
        r["unexpected"] = True
    elif kind == "manifest_schema":
        m["schema"] = "wrong"
    elif kind == "manifest_extra":
        m["unexpected"] = True
    elif kind == "ancestry":
        r["required_J70"] = BASE
    elif kind == "nested_boolean":
        r["offline_reconstruction"]["validation_untouched"] = 1
    elif kind == "nested_count":
        r["offline_reconstruction"]["new_QM"] = 1
    if kind == "valid":
        semantics(r, m)
    else:
        with pytest.raises(CompletionError):
            semantics(resign(r), resign(m))


def test_canonical_selector_actual_fail_closed():
    from island.pcff_pdms import PCFFPDMSAssets, prepare_polymer

    assets = (
        REPO.parent
        / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
    )
    if not assets.is_dir():
        pytest.skip(
            "declared source unavailable; offline capability report remains testable"
        )
    system, prepared = prepare_polymer(
        psmiles="[*:1]O[Si](C)(C)[*:2]",
        forcefield="pcff",
        dp=1,
        head_end_group="[*][Si](C)(C)C",
        tail_end_group="[*]C",
        assets=PCFFPDMSAssets(assets),
    )
    p = prepared.validate_integrity(system)
    assert len(p["missing_terms"]) == 378 and p["selected_profile"] is None
    assert not p["parameter_complete"]
    with pytest.raises(Exception, match="parameter_complete"):
        prepared.require()


def test_no_scientific_execution_and_validation_nonactivation(retained, monkeypatch):
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: pytest.fail(
            "literature replay cannot launch scientific subprocesses"
        ),
    )
    assert reports(REPO, E) == retained[3]
    d = retained[3]["decision.json"]
    assert d["new_QM"] == d["MD_steps"] == d["QM_fits"] == 0
    assert (
        not d["reserved_validation_activated"] and not d["reserved_validation_consumed"]
    )


def test_fresh_process_relocation_offline_without_pdfs(tmp_path):
    for name in ("src", "scripts"):
        shutil.copytree(
            REPO / name, tmp_path / name, ignore=shutil.ignore_patterns("__pycache__")
        )
    shutil.copytree(E, tmp_path / ROOT)
    parent = tmp_path / "docs/evidence/phase_4j70"
    parent.mkdir(parents=True)
    for name in (
        "receipt.json",
        "manifest.json",
        "completeness-matrix.json",
        "canonical-capability.json",
        "qm-fit-targets.json",
    ):
        shutil.copyfile(REPO / "docs/evidence/phase_4j70" / name, parent / name)
    assert not (tmp_path / "island_ffref").exists()
    env = {**os.environ, "PYTHONPATH": str(tmp_path) + ":" + str(tmp_path / "src")}
    p = subprocess.run(
        [sys.executable, "scripts/replay_phase_4j71.py", "--repo", str(tmp_path)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert p.returncode == 0, p.stdout + p.stderr
    assert "source_conflict" in p.stdout and "'new_QM': 0" in p.stdout


def test_offline_reconstruction_rejects_resigned_semantic_mutation(tmp_path):
    shutil.copytree(E, tmp_path / "j71")
    p = read(tmp_path / "j71/decision.json")
    p["accepted_coefficients"] = 378
    (tmp_path / "j71/decision.json").write_text(json.dumps(resign(p)))
    with pytest.raises(CompletionError, match="semantic_reconstruction"):
        replay(REPO, tmp_path / "j71")


def test_manufactured_AA_mapping_preserves_center_and_shared_neighbor():
    r = {
        "interaction_family": "angle-angle",
        "source_model_resolution": "all_atom",
        "source_force_field_family": "PCFF",
        "missing_convention_fields": [],
        "extraction_confidence": "direct_printed_numeric_cells",
    }
    t = {"family": "angle-angle", "ordered_types": ["hc", "c3", "sio", "c3"]}
    for sites in (["H", "C", "Si", "C"], ["C", "C", "Si", "H"]):
        assert (
            compatibility(r, t, sites, same_environment=True, complete_contract=True)
            == "exact_literature_external_match"
        )
    assert (
        compatibility(
            r, t, ["Si", "C", "H", "C"], same_environment=True, complete_contract=True
        )
        == "not_applicable"
    )
