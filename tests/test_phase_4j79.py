"""J79 software contracts. No synthetic control is scientific qualification."""

import json
import os
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import pytest

from scripts import phase_4j79_literature as j79
from scripts import validate_phase_4j79_receipt as validator

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / j79.ROOT


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j79.scope()


@pytest.fixture(scope="module")
def docs():
    return j79.reconstruct(ROOT)


def test_provenance_is_metadata_only(docs):
    evidence = docs["literature-evidence.json"]
    assert evidence["bibliography"]["DOI"] == "10.1016/S0032-3861(96)00600-3"
    assert evidence["primary_artifact_sha256"] is None
    assert evidence["parameter_rows_source_extracted"] == 0
    assert "not_read" in evidence["reading_status"]
    access = j79.read(ROOT / "source-access.json")
    assert access["full_text_read"] is False
    assert access["publisher_XML_contains_full_text"] is False
    assert len(access["crossref_response_sha256"]) == 64
    assert access["crossref_identity"]["page"] == "1035-1045"
    assert any(
        r.get("error") == "HTTP Error 403: Forbidden" for r in access["access_attempts"]
    )


TRANSCRIPTION = (
    "partial_charge:Si +0.640",
    "partial_charge:O -0.440",
    "partial_charge:C -0.259",
    "partial_charge:H +0.053",
    "quartic_bond:Si-O 1.6562 306.1232 -517.3424 673.7067",
    "quartic_bond:Si-C 1.9073 157.0049 -237.7023 356.0328",
    "quartic_bond:C-H 1.1010 345.0000 -691.8900 844.6000",
    "quartic_angle:Si-O-Si 157.0260 9.0740 -19.5576 8.5000",
    "quartic_angle:O-Si-O 110.6930 70.3069 -6.9375 0",
    "quartic_angle:C-Si-O 114.9060 23.0218 -31.3993 -924.9814",
    "quartic_angle:H-C-Si 111.5360 30.2481 -15.5255 0",
    "quartic_angle:C-Si-C 114.9060 23.0218 -31.3993 24.9814",
    "torsion_3:Si-O-Si-O 0.3000 0.3658 0",
    "torsion_3:Si-O-Si-C 0 0 -0.1300",
    "torsion_3:H-C-Si-O -1.3513 0 -0.0580",
    "torsion_3:H-C-Si-C -1.3513 0 -0.0580",
    "lj_9_6:Si 4.2840 0.07000",
    "lj_9_6:O 3.3500 0.24000",
    "lj_9_6:C 4.0100 0.05400",
    "lj_9_6:H 2.9950 0.02000",
)


@pytest.mark.parametrize("line", TRANSCRIPTION)
def test_exact_reported_decimal_transcription(line):
    key, *expected = line.split()
    row = next(r for r in j79.reported_rows() if r["id"] == key)
    assert list(row["raw_reported_values"].values()) == expected


@pytest.mark.parametrize(
    "field",
    (
        "source_hash",
        "source_page",
        "source_table",
        "source_equation",
        "source_units",
        "phase_convention",
        "orientation_rule",
    ),
)
def test_missing_primary_fields_are_explicit(docs, field):
    assert all(
        r[field] is None for r in docs["literature-evidence.json"]["reported_rows"]
    )


def test_repeat_neutrality_does_not_authorize_caps(docs):
    q = docs["convention-compatibility.json"]["charge_arithmetic_only"]
    assert Decimal(q["repeat_SiOC2H6_charge_e"]) == 0
    for name in ("hmds", "dp2", "dp4", "dp6"):
        assert Decimal(q[name]["charge_if_repeat_values_transferred_e"]) == Decimal(
            "0.240"
        )
        assert q[name]["transfer_allowed"] is False


def test_9_6_and_12_6_are_distinct_models(docs):
    c = docs["convention-compatibility.json"]
    assert c["generic_parameter_model"]["compatible_with_reported_paper_9_6"] is False
    assert c["J78_proposal"]["symbolic_form_match_to_reported_paper"] is True
    assert c["isolated_PCFF_model"]["called"] is False
    assert c["isolated_PCFF_model"]["paper_mixing_compatibility"].startswith(
        "unresolved"
    )
    assert c["paper_conventions"]["mixing_rule"]["value"] is None
    assert c["paper_conventions"]["dielectric"]["value"] is None
    assert c["paper_conventions"]["cutoff"]["value"] is None


def test_unit_conversion_is_arithmetic_not_reinterpretation(docs):
    assert j79.convert("4.2840", "angstrom_to_nm") == "0.42840"
    assert Decimal(j79.convert("4.184", "kJ_to_kcal")) == 1
    row = next(
        r
        for r in docs["literature-evidence.json"]["reported_rows"]
        if r["id"] == "lj_9_6:Si"
    )
    assert row["raw_reported_values"]["epsilon"] == "0.07000"
    assert row["reported_equation"] == j79.VDW
    assert (
        "IF source"
        in row["derived_conditional_conversions"]["epsilon_energy_kcal"]["condition"]
    )
    assert (
        docs["convention-compatibility.json"]["angle"]["applied_to_coefficients"]
        is False
    )
    with pytest.raises(ValueError, match="explicit_unit"):
        j79.convert("1", "LJ9_to_LJ12")


def test_torsion_sign_and_unknown_phase_are_retained(docs):
    row = next(
        r
        for r in docs["literature-evidence.json"]["reported_rows"]
        if r["id"] == "torsion_3:H-C-Si-O"
    )
    assert row["raw_reported_values"]["V1"] == "-1.3513"
    assert row["raw_reported_values"]["V3"] == "-0.0580"
    with pytest.raises(ValueError, match="torsion_sign_phase"):
        j79.torsion_contract(row)


@pytest.mark.parametrize(
    "key",
    (
        "source_verified",
        "source_equation",
        "prefactor",
        "phase_convention",
        "angle_origin",
        "orientation_rule",
        "multiplicities",
    ),
)
def test_incomplete_synthetic_torsion_contract_rejected(key):
    # Pure schema control; synthetic metadata must never enter source evidence.
    r = {
        "source_verified": True,
        "source_equation": "synthetic",
        "prefactor": 1,
        "phase_convention": "explicit",
        "angle_origin": "explicit",
        "orientation_rule": "explicit",
        "multiplicities": [1, 2, 3],
    }
    assert j79.torsion_contract(r) == r
    r[key] = None
    with pytest.raises(ValueError):
        j79.torsion_contract(r)


def test_all_terms_and_site_identities_mapped(docs):
    walk = docs["j78-literature-crosswalk.json"]
    inv = j79.read(ROOT / "j78-input.json")
    assert walk["unique_terms"] == 129
    assert walk["site_instances"] == 5135
    assert {r["J78_parameter_id"] for r in walk["rows"]} == {
        r["id"] for r in inv["terms"]
    }
    assert sum(walk["status_counts"].values()) == 129
    assert walk["status_counts"]["literature_candidate_exact_form"] == 0
    assert {i for r in walk["rows"] for i in r["J78_instance_ids"]} == {
        i["id"] for i in inv["instances"]
    }
    assert all(
        not r["extension_term_resolved"] and not r["native_row_resolved"]
        for r in walk["rows"]
    )


def test_element_reversal_is_only_candidate_matching():
    term = {
        "family": "torsion_3",
        "ordered_types": ["ext_C_D", "ext_Si_D", "ext_O_bridge", "ext_Si_D"],
    }
    assert j79.candidates(term, j79.reported_rows()) == ["torsion_3:Si-O-Si-C"]
    term["ordered_types"] = ["ext_C_D", "ext_O_bridge", "ext_Si_D", "ext_Si_D"]
    assert j79.candidates(term, j79.reported_rows()) == []


def test_cross_terms_and_cap_blocking(docs):
    rows = docs["j78-literature-crosswalk.json"]["rows"]
    assert all(r["status"] == "cross_term_unresolved" for r in rows if r["cross_term"])
    assert all(
        r["status"] == "end_group_unresolved"
        for r in rows
        if r["end_group_involved"] and not r["cross_term"]
    )
    classes = docs["unresolved-terms.json"]["required_classes"]
    assert len(classes) == 7 and all(r["status"] == "unresolved" for r in classes)
    assert all("not verifiable" in r["reason"] for r in classes)


@pytest.mark.parametrize("mutation", ("missing", "duplicate"))
def test_incomplete_and_duplicate_term_coverage_rejected(mutation):
    inv = j79.read(ROOT / "j78-input.json")
    if mutation == "missing":
        inv["terms"].pop()
    else:
        inv["terms"][-1] = inv["terms"][0]
    with pytest.raises(ValueError):
        j79.crosswalk(inv, j79.literature())


@pytest.mark.parametrize(
    "mutation",
    ("raw_value", "page", "table", "unit", "source_hash", "phase", "native", "profile"),
)
def test_resigned_scientific_metadata_mutation_rejected(tmp_path, mutation):
    shutil.copytree(ROOT, tmp_path / "evidence")
    root = tmp_path / "evidence"
    p = root / "literature-evidence.json"
    d = j79.read(p)
    row = d["reported_rows"][0]
    if mutation == "raw_value":
        row["raw_reported_values"]["q"] = "+0.641"
    else:
        field, value = {
            "page": ("source_page", 1036),
            "table": ("source_table", "I"),
            "unit": ("source_units", "kcal/mol"),
            "source_hash": ("source_hash", "0" * 64),
            "phase": ("phase_convention", "guessed"),
            "native": ("native_row", True),
            "profile": ("usable_for_runtime_assignment", True),
        }[mutation]
        row[field] = value
    row.pop("identity")
    d["reported_rows"][0] = j79.signed(row)
    d.pop("identity")
    j79.write(p, j79.signed(d))
    with pytest.raises(ValueError, match="semantic_mutation"):
        j79.replay(root)


@pytest.mark.parametrize(
    "name",
    (
        "j78-input.json",
        "source-access.json",
        "model-inputs.json",
        "local-source-search.json",
    ),
)
def test_captured_source_hash_mutation_rejected(tmp_path, name):
    shutil.copytree(ROOT, tmp_path / "evidence")
    p = tmp_path / "evidence" / name
    p.write_bytes(p.read_bytes() + b" ")
    with pytest.raises(ValueError, match="source_hash_mutation"):
        j79.reconstruct(p.parent)


@pytest.mark.parametrize(
    "field",
    (
        "accepted_coefficients",
        "native_exact_rows_resolved",
        "runtime_registered",
        "automatic_assignment_enabled",
        "production_validated",
        "J79_MD_steps",
    ),
)
def test_resigned_receipt_readiness_rejected(field):
    receipt = {
        "schema": "island.j79-receipt.v1",
        "base": j79.BASE,
        "branch": j79.BRANCH,
        "scope": dict(j79.SCOPE),
        "decision": "literature_reconciliation_partial",
    }
    manifest = {
        "schema": "island.j79-manifest.v1",
        "base": j79.BASE,
        "branch": j79.BRANCH,
        "scope": dict(j79.SCOPE),
        "freeze_sha256": j79.FREEZE_SHA,
    }
    receipt["scope"][field] = True if isinstance(receipt["scope"][field], bool) else 1
    with pytest.raises(ValueError, match="readiness_or_scope"):
        validator.semantic(j79.signed(receipt), j79.signed(manifest))


def test_no_runtime_registration_or_automatic_assignment(docs):
    assert j79.SCOPE["runtime_registered"] is False
    assert j79.SCOPE["automatic_assignment_enabled"] is False
    assert j79.SCOPE["accepted_coefficients"] == 0
    assert docs["reconciliation-report.json"]["scope"] == j79.SCOPE
    native = j79.read(ROOT / "model-inputs.json")
    assert "island_pcff_pdms_derived_v1" not in native[j79.MODEL_FILES[3]]


@pytest.mark.parametrize(
    "name", ("sqm", "antechamber", "lammps", "msi2lmp", "phase_4j79_fake_fit")
)
def test_fake_child_scientific_executable_never_runs(tmp_path, name):
    marker = tmp_path / "launched"
    fake = tmp_path / name
    fake.write_text(
        f"#!{sys.executable}\nfrom pathlib import Path\nPath({str(marker)!r}).touch()\n"
    )
    fake.chmod(0o755)
    with pytest.raises((PermissionError, RuntimeError)):
        subprocess.run([str(fake)], check=True)
    assert not marker.exists()


def test_python_child_inherits_fence(tmp_path):
    code = "from scripts.phase_4j76_source_resolution import scope; scope(); import importlib.util; assert importlib.util.find_spec('json'); import openmm"
    p = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )
    assert p.returncode != 0
    assert "J73_scope" in p.stderr or "scope" in p.stderr.lower()


def test_fresh_process_relocated_offline_data_replay(tmp_path):
    root = tmp_path / "relocated"
    shutil.copytree(ROOT, root)
    p = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j79_literature",
            "replay",
            "--directory",
            str(root),
        ],
        cwd=REPO,
        check=False,
        capture_output=True,
        text=True,
    )
    assert p.returncode == 0, p.stderr + p.stdout
    assert json.loads(p.stdout)["decision"] == "literature_reconciliation_partial"
    assert j79.replay(root) == j79.replay(ROOT)


def test_original_J78_and_native_evidence_preserved():
    assert validator.history(REPO, ROOT) > 10154


@pytest.mark.parametrize(
    "control",
    ("real_QM", "real_MD", "converter", "fitting", "evaluator", "reserved_validation"),
)
def test_prohibited_scientific_controls_explicitly_skipped(control):
    pytest.skip(
        "J79 evidence-only scope: " + control + " is not authorized or launched"
    )
