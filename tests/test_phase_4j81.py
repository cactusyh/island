"""Synthetic parser/receipt controls never establish scientific parameter data."""

import copy
import io
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from island.pcff_pdms.pilot_intake import gates
from island.pcff_pdms.source_ingestion import parsing, pipeline
from island.pcff_pdms.source_ingestion.schema import PARSER_VERSION, validate
from scripts import phase_4j81_ingest as j81
from scripts import validate_phase_4j81_receipt as validator
from tests.test_phase_4j80 import fixture_contract

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j81.scope()


def bundle(tmp_path, kind="original_pcff_frc"):
    raw = b"! synthetic parser control only\r\n#quartic_bond synthetic\r\n1.0 1 t_si t_o 1.5 2.0 -3.0 4.0 ! original comment\r\n"
    (tmp_path / "fixture.frc").write_bytes(raw)
    _, _, _, contract = fixture_contract()
    fields = contract["fields"]
    fields["internal_and_end_group_typing"]["value"]["source_type_ids"] = [
        "t_si",
        "t_o",
    ]
    row = {
        "id": "synthetic-candidate",
        "source_id": "synthetic-source",
        "source_location": {
            "line_start": 3,
            "line_end": 3,
            "section": "quartic_bond",
            "parameter_record_id": "synthetic-record-1",
        },
        "raw_text": raw.decode().splitlines(keepends=True)[2],
        "atom_type_sequence": ["t_si", "t_o"],
        "element_sequence": ["Si", "O"],
        "parameter_family": "bond",
        "numerical_values": ["1.5", "2.0", "-3.0", "4.0"],
        "units": copy.deepcopy(fields["energy_and_length_units"]["value"]),
        "functional_form": fields["bond_and_angle_polynomials"]["value"][
            "bond_equation"
        ],
        "mixing_rule": copy.deepcopy(fields["combination_rules"]["value"]),
        "exclusion_and_1_4_policy": copy.deepcopy(
            fields["pair_exclusions_and_1_4"]["value"]
        ),
        "charge_convention": copy.deepcopy(fields["charges_and_neutrality"]["value"]),
        "torsion_convention": copy.deepcopy(fields["torsion"]["value"]),
        "chain_scope": "internal",
        "reporting_mode": "directly_reported",
        "chemical_environment": {"synthetic_exact_graph": "not scientific"},
        "native_identity": {},
    }
    row["units"]["per_value"] = [
        "angstrom",
        "kJ/mol/angstrom^2",
        "kJ/mol/angstrom^3",
        "kJ/mol/angstrom^4",
    ]
    artifact = {
        "id": "synthetic-source",
        "source_path": "fixture.frc",
        "source_type": kind,
        "sha256": parsing.sha(raw),
        "file_size": len(raw),
        "bibliography": {
            "title": "Synthetic software control",
            "authors": ["Fixture"],
            "year": 2000,
            "identifier": "synthetic:fixture",
        },
        "source_version": "synthetic-v1",
        "acquisition_method": "explicit_local_file",
        "authorization_status": "authorized_local_use",
        "license_or_redistribution_status": "synthetic_only",
        "extraction_tool": "island.j81-text",
        "extraction_version": "1",
        "extraction_timestamp": "2026-11-05T00:00:00Z",
        "parser_version": PARSER_VERSION,
        "authority": "original_author",
        "source_family": "published_PDMS",
    }
    result = j81.empty_bundle()
    result.update(artifacts=[artifact], candidates=[row], convention_contract=contract)
    return result


def audited(snapshot):
    mapped = pipeline.bridge(snapshot)
    artifact = mapped["artifacts"][0]
    return {
        artifact["id"]: {
            "artifact_metadata": copy.deepcopy(artifact),
            "capture_identity": pipeline.digest(snapshot),
            "verified_rows": {
                r["candidate_id"]: copy.deepcopy(r) for r in mapped["row_evidence"]
            },
            "convention_assertions": copy.deepcopy(
                snapshot["manifest"]["convention_contract"]["fields"]
            ),
            "evidence_kind": "independently_reviewed_primary_artifact",
        }
    }


def run(tmp_path, manifest=None, trust=True):
    m = bundle(tmp_path) if manifest is None else manifest
    snap = pipeline.capture(m, tmp_path)
    audits = audited(snap) if trust else {}
    return snap, audits, pipeline.normalize(snap, audits, ["synthetic-candidate"])


def test_missing_source_stops_conventions():
    snap = pipeline.capture(j81.empty_bundle(), REPO)
    r = pipeline.normalize(snap, {}, ["missing-candidate"])
    assert r["decision"] == "source_ingestion_blocked"
    assert r["convention_comparison"]["status"] == "not_evaluated"
    assert r["verified_literature_rows"] == 0
    assert r["candidate_rows"] == []
    assert r["J80_preserved"]["unresolved_cross_terms"] == 87


def test_raw_frc_format_locations_and_decimal_values_preserved(tmp_path):
    snap, _, r = run(tmp_path)
    row = r["candidate_rows"][0]
    assert row["raw_text"].endswith("! original comment\r\n")
    assert row["raw_values"] == ["1.5", "2.0", "-3.0", "4.0"]
    assert row["source_location"]["line_start"] == 3
    assert row["verification_status"] == "literature_candidate_verified"
    assert r["decision"] == "source_parsed_ready_for_qm_authorization"
    assert r["native_exact_rows"] == 0
    assert not r["qm_authorization_granted"]
    assert snap["artifacts"][0]["status"] == "source_parsed"


def test_unattested_source_is_only_parsed(tmp_path):
    _, _, r = run(tmp_path, trust=False)
    assert r["decision"] == "source_parsed_convention_blocked"
    assert r["verified_literature_rows"] == 0
    assert r["convention_comparison"]["status"] == "not_evaluated"


@pytest.mark.parametrize("field", ("sha256", "file_size"))
def test_source_hash_size_mutation(tmp_path, field):
    m = bundle(tmp_path)
    m["artifacts"][0][field] = "0" * 64 if field == "sha256" else 1
    s = pipeline.capture(m, tmp_path)
    assert s["artifacts"][0]["status"] == "source_present_unparsed"
    assert "mutation" in s["artifacts"][0]["errors"][0]
    assert (
        pipeline.normalize(s, {}, ["synthetic-candidate"])["decision"]
        == "source_ingestion_blocked"
    )


@pytest.mark.parametrize(
    "field",
    (
        "source_version",
        "authorization_status",
        "extraction_timestamp",
        "license_or_redistribution_status",
        "bibliography",
        "parser_version",
    ),
)
def test_required_provenance_metadata(tmp_path, field):
    m = bundle(tmp_path)
    del m["artifacts"][0][field]
    with pytest.raises(ValueError, match="schema_missing"):
        pipeline.capture(m, tmp_path)


@pytest.mark.parametrize(
    "change", ("line", "section", "family", "raw", "value", "types")
)
def test_exact_frc_location_and_transcription(tmp_path, change):
    m = bundle(tmp_path)
    c = m["candidates"][0]
    if change == "line":
        c["source_location"]["line_start"] = 2
    elif change == "section":
        c["source_location"]["section"] = "quartic_angle"
    elif change == "family":
        c["parameter_family"] = "angle"
    elif change == "raw":
        c["raw_text"] = c["raw_text"].replace("\r\n", "\n")
    elif change == "value":
        c["numerical_values"][0] = "99.9"
    else:
        c["atom_type_sequence"].reverse()
    _, _, r = run(tmp_path, m, False)
    assert r["decision"] == "source_ingestion_blocked"
    assert not r["candidate_rows"][0]["extracted_fact"]


@pytest.mark.parametrize("change", ("source", "candidate", "location"))
def test_duplicate_identities_cannot_hide_conflicts(tmp_path, change):
    m = bundle(tmp_path)
    if change == "source":
        m["artifacts"].append(copy.deepcopy(m["artifacts"][0]))
    else:
        c = copy.deepcopy(m["candidates"][0])
        if change == "location":
            c["id"] = "different-id-same-cell"
            c["numerical_values"] = ["99"]
        m["candidates"].append(c)
    with pytest.raises(ValueError, match="duplicate"):
        pipeline.capture(m, tmp_path)


@pytest.mark.parametrize("field", list(gates.CONVENTIONS))
def test_every_missing_convention_blocks_even_audited(tmp_path, field):
    m = bundle(tmp_path)
    m["convention_contract"]["fields"][field]["value"] = None
    _, _, r = run(tmp_path, m)
    assert not r["convention_comparison"]["passed"]
    assert (
        r["candidate_rows"][0]["verification_status"]
        == "literature_candidate_convention_blocked"
    )


@pytest.mark.parametrize(
    "field",
    (
        "units",
        "functional_form",
        "mixing_rule",
        "exclusion_and_1_4_policy",
        "torsion_convention",
    ),
)
def test_row_convention_mismatch_blocks(tmp_path, field):
    m = bundle(tmp_path)
    m["candidates"][0][field] = "wrong"
    _, _, r = run(tmp_path, m)
    assert r["verified_literature_rows"] == 0
    assert r["decision"] == "source_parsed_convention_blocked"


def test_no_9_6_to_12_6_conversion(tmp_path):
    m = bundle(tmp_path)
    v = m["convention_contract"]["fields"]["vdw_function"]["value"]
    v["target_family"] = "lj_12_6"
    _, _, r = run(tmp_path, m)
    assert (
        "no_9_6_to_12_6_coercion_or_analog_form"
        in r["convention_comparison"]["blocking_reasons"]
    )
    assert r["candidate_rows"][0]["converted_values"] is None


@pytest.mark.parametrize(
    "family", ("COMPASS", "OPLS", "GAFF", "CFF91", "CVFF", "analog")
)
def test_cross_family_analog_not_native(tmp_path, family):
    m = bundle(tmp_path)
    m["artifacts"][0]["source_family"] = family
    _, _, r = run(tmp_path, m)
    assert (
        r["candidate_rows"][0]["verification_status"]
        == "inferred_or_analog_value_blocked"
    )
    assert r["native_exact_rows"] == 0


@pytest.mark.parametrize("mode", ("inferred", "reconstructed"))
def test_non_reported_value_blocked(tmp_path, mode):
    m = bundle(tmp_path)
    m["candidates"][0]["reporting_mode"] = mode
    _, _, r = run(tmp_path, m)
    assert (
        r["candidate_rows"][0]["verification_status"]
        == "inferred_or_analog_value_blocked"
    )


def test_caps_and_missing_cross_terms_remain_unresolved(tmp_path):
    m = bundle(tmp_path)
    m["candidates"][0]["chain_scope"] = "capped_end_group"
    _, _, r = run(tmp_path, m, False)
    assert r["candidate_rows"][0]["verification_status"] == "end_group_scope_blocked"
    assert len(r["missing_cross_term_classes"]) == 9
    assert r["unresolved_cross_terms"] == 87 and r["unresolved_end_group_terms"] == 69


def test_source_and_audit_mutation_cannot_be_self_signed(tmp_path):
    snap, audit, _ = run(tmp_path)
    snap["manifest"]["candidates"][0]["source_location"]["parameter_record_id"] = (
        "invented"
    )
    r = pipeline.normalize(snap, audit, ["synthetic-candidate"])
    assert not r["source_gate"]["passed"]
    assert r["verified_literature_rows"] == 0


def make_pdf(tmp_path):
    from pypdf import PdfReader
    from reportlab.pdfgen.canvas import Canvas

    buffer = io.BytesIO()
    canvas = Canvas(buffer, invariant=1)
    snippet = "t_si t_o 1.5 2.0 -3.0 4.0"
    canvas.drawString(30, 700, snippet)
    canvas.save()
    raw = buffer.getvalue()
    m = bundle(tmp_path, "primary_pdf")
    (tmp_path / "synthetic.pdf").write_bytes(raw)
    m["artifacts"][0].update(
        source_path="synthetic.pdf",
        sha256=parsing.sha(raw),
        file_size=len(raw),
        extraction_tool="pypdf",
        extraction_version="6.1.1",
    )
    text = PdfReader(io.BytesIO(raw)).pages[0].extract_text()
    m["candidates"][0].update(
        raw_text=snippet,
        source_location={
            "page": 1,
            "table": "synthetic-table",
            "row": "1",
            "column": "types,values",
            "figure": None,
            "appendix": None,
            "extracted_text_sha256": parsing.sha(text.encode()),
            "ocr_status": "not_used",
            "ocr_confidence": None,
        },
    )
    return m


def test_pdf_exact_page_text_checksum(tmp_path):
    _, _, r = run(tmp_path, make_pdf(tmp_path))
    row = r["candidate_rows"][0]
    assert row["verification_status"] == "literature_candidate_verified"
    assert row["source_location"]["page"] == 1
    assert row["extracted_fact"]["parser"]["version"] == "6.1.1"


@pytest.mark.parametrize("mutation", ("page", "checksum", "row", "ocr"))
def test_pdf_location_checksum_and_ocr_fail_closed(tmp_path, mutation):
    m = make_pdf(tmp_path)
    loc = m["candidates"][0]["source_location"]
    if mutation == "page":
        loc["page"] = 10
    elif mutation == "checksum":
        loc["extracted_text_sha256"] = "0" * 64
    elif mutation == "row":
        loc["row"] = None
    else:
        loc.update(ocr_status="required", ocr_confidence=0.99)
    _, _, r = run(tmp_path, m, False)
    assert r["decision"] == "source_ingestion_blocked"


def test_archive_member_hash_no_extraction_or_execution(tmp_path):
    m = bundle(tmp_path, "authorized_archive")
    content = (tmp_path / "fixture.frc").read_bytes()
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        z.writestr("original.frc", content)
        z.writestr("do_not_run.py", "raise RuntimeError('MUST NEVER RUN')")
    raw = buffer.getvalue()
    (tmp_path / "fixture.zip").write_bytes(raw)
    m["artifacts"][0].update(
        source_path="fixture.zip", sha256=parsing.sha(raw), file_size=len(raw)
    )
    m["candidates"][0]["source_location"].update(
        archive_member="original.frc", member_sha256=parsing.sha(content)
    )
    snap, _, r = run(tmp_path, m)
    assert r["decision"] == "source_parsed_ready_for_qm_authorization"
    assert len(snap["artifacts"][0]["extraction"]["archive_members"]) == 2
    assert not (tmp_path / "original.frc").exists()


@pytest.mark.parametrize("name", ("../escape", "/absolute", "a\\b"))
def test_unsafe_archives_rejected(name):
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w") as z:
        z.writestr(name, "synthetic")
    with pytest.raises(ValueError, match="unsafe_archive"):
        parsing.archive(b.getvalue())


def test_explicit_local_only_symlink_and_missing_source(tmp_path):
    m = bundle(tmp_path)
    (tmp_path / "fixture.frc").unlink()
    s = pipeline.capture(m, tmp_path)
    assert s["artifacts"][0]["status"] == "source_missing"
    (tmp_path / "fixture.frc").symlink_to("elsewhere")
    assert "symlink" in pipeline.capture(m, tmp_path)["artifacts"][0]["errors"][0]
    with pytest.raises(ValueError, match="explicit_local"):
        parsing.local_bytes(tmp_path, "https://publisher.example/forbidden.pdf")


def test_native_requires_exact_types_and_complete_identity(tmp_path):
    m = bundle(tmp_path)
    m["artifacts"][0]["source_family"] = "native_PCFF"
    _, _, r = run(tmp_path, m)
    assert r["native_exact_rows"] == 0
    assert (
        r["candidate_rows"][0]["verification_status"]
        == "inferred_or_analog_value_blocked"
    )


def test_schema_matches_generated_bundle(tmp_path):
    validate(bundle(tmp_path))
    assert (
        j81.schema()["properties"]["artifacts"]["items"]["additionalProperties"]
        is False
    )


def test_pinned_manifest_and_audit_required(tmp_path):
    p = tmp_path / "bundle.json"
    j81.write(p, bundle(tmp_path))
    with pytest.raises(ValueError, match="external_hash"):
        j81.prepare(REPO, p)
    with pytest.raises(ValueError, match="checksum"):
        j81.prepare(REPO, p, "0" * 64)


def test_snapshot_and_resigned_reports_cannot_grant_verification(tmp_path):
    inputs = j81.prepare(REPO)
    out = tmp_path / "out"
    j81.emit(out, inputs)
    expected = j81.pin(out / "captured-input.json")
    assert j81.replay(out, expected)["decision"] == "source_ingestion_blocked"
    data = j81.read(out / "result.json")
    data["verified_literature_rows"] = 20
    data = j81.signed({k: v for k, v in data.items() if k != "identity"})
    j81.write(out / "result.json", data)
    with pytest.raises(ValueError, match="semantic_mutation"):
        j81.replay(out, expected)
    j81.write(out / "captured-input.json", {"fabricated": True})
    with pytest.raises(ValueError, match="external_pin"):
        j81.replay(out, expected)


def test_relocated_fresh_process_without_sources(tmp_path):
    out = tmp_path / "bundle with spaces"
    j81.emit(out, j81.prepare(REPO))
    old = tmp_path / "relocated"
    shutil.copytree(out, old)
    shutil.rmtree(out)
    code = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j81_ingest",
            "replay",
            "--directory",
            str(old),
            "--snapshot-sha256",
            j81.pin(old / "captured-input.json"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert code.returncode == 0, code.stderr
    assert json.loads(code.stdout)["decision"] == "source_ingestion_blocked"


@pytest.mark.parametrize(
    "name", ("sqm", "msi2lmp", "lammps", "antechamber", "sander", "pmemd")
)
def test_forbidden_synthetic_child_never_executes(tmp_path, name):
    marker = tmp_path / "ran"
    executable = tmp_path / name
    executable.write_text("#!/bin/sh\ntouch " + str(marker) + "\n")
    executable.chmod(0o755)
    child = subprocess.run(
        [
            sys.executable,
            "-c",
            "import subprocess,sys; from scripts.phase_4j76_source_resolution import scope; scope(); subprocess.run([sys.argv[1]])",
            str(executable),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert child.returncode != 0
    assert not marker.exists()
    assert "scope" in child.stderr.lower()


def test_no_runtime_registration_or_assignment_from_ready_control(tmp_path):
    _, _, r = run(tmp_path)
    assert r["decision"] == "source_parsed_ready_for_qm_authorization"
    assert r["scope"] == pipeline.SCOPE
    assert not r["scope"]["runtime_registration"]
    assert not r["scope"]["automatic_assignment"]
    assert r["scope"]["accepted_parameters"] == 0
    assert not r["qm_authorization_granted"]


@pytest.mark.parametrize(
    "control",
    ("QM", "physical fitting", "MD", "converter", "evaluator", "reserved validation"),
)
def test_scientific_controls_explicitly_not_run(control):
    pytest.skip("J81 evidence-only scope: " + control + " is not launched")


@pytest.mark.parametrize("kind", ("materials_studio_export", "parameter_table"))
def test_other_authorized_text_sources(tmp_path, kind):
    _, _, result = run(tmp_path, bundle(tmp_path, kind))
    assert (
        result["candidate_rows"][0]["verification_status"]
        == "literature_candidate_verified"
    )


def test_corrupt_pdf_and_archive_are_blocked(tmp_path):
    for kind in ("primary_pdf", "authorized_archive"):
        m = bundle(tmp_path, kind)
        raw = b"%PDF-1.7\ncorrupt" if kind == "primary_pdf" else b"not an archive"
        (tmp_path / "fixture.frc").write_bytes(raw)
        m["artifacts"][0].update(sha256=parsing.sha(raw), file_size=len(raw))
        s = pipeline.capture(m, tmp_path)
        assert s["artifacts"][0]["status"] == "source_present_unparsed"
        assert (
            pipeline.normalize(s, {}, ["synthetic-candidate"])["decision"]
            == "source_ingestion_blocked"
        )


def test_synthetic_exact_native_does_not_register(tmp_path):
    m = bundle(tmp_path)
    m["artifacts"][0]["source_family"] = "native_PCFF"
    c = m["candidates"][0]
    c["native_identity"] = {
        "ordered_source_types": list(c["atom_type_sequence"]),
        "upstream_release": "synthetic-not-a-native-release",
        "source_row_id": "synthetic",
        "equation": c["functional_form"],
        "ordering": "exact_source_order",
        "exact_environment": copy.deepcopy(c["chemical_environment"]),
        "units_and_phase": {"synthetic": "documented"},
        "exclusion_and_mixing_policy": {"synthetic": "documented"},
        "provenance": {"synthetic": "not scientific"},
    }
    snap, audit, r = run(tmp_path, m)
    assert r["candidate_rows"][0]["verification_status"] == "exact_native_row"
    assert not r["scope"]["runtime_registration"]
    snap["manifest"]["candidates"][0]["atom_type_sequence"] = ["Si", "O"]
    r = pipeline.normalize(snap, audit, ["synthetic-candidate"])
    assert r["native_exact_rows"] == 0


@pytest.mark.parametrize("ready", (False, True))
def test_blocked_and_verified_status_receipts(tmp_path, monkeypatch, ready):
    inputs = j81.prepare(REPO)
    if ready:
        snap, audit, _ = run(tmp_path)
        inputs.update(
            capture=snap,
            trusted_audits=audit,
            expected_J79_candidate_ids=["synthetic-candidate"],
        )
    root = tmp_path / "sealed"
    j81.emit(root, inputs)
    pin = j81.pin(root / "captured-input.json")
    # Small manufactured receipt control; actual historical audit runs separately.
    monkeypatch.setattr(validator, "historical", lambda repo: {})
    monkeypatch.setattr(validator, "verification", lambda repo, root: {})
    validator.seal(REPO, root, pin)
    r = validator.validate(REPO, root, pin)
    expected = (
        "source_parsed_ready_for_qm_authorization"
        if ready
        else "source_ingestion_blocked"
    )
    assert r["decision"] == expected
    data = j81.read(root / "receipt.json")
    data["scope"]["runtime_registration"] = True
    j81.write(
        root / "receipt.json",
        j81.signed({k: v for k, v in data.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="resigned_scope"):
        validator.validate(REPO, root, pin)


def test_historical_target_identities_and_counts():
    inputs = j81.prepare(REPO)
    docs = j81.reconstruct(inputs)
    targets = docs["unresolved-terms.json"]
    assert len(targets["extension_targets"]) == 129
    assert len(targets["cross_term_ids"]) == 87
    assert len(targets["end_group_term_ids"]) == 69
    assert docs["result.json"]["blocked_rows"] == 20
    assert (
        j81.pin(REPO / "docs/evidence/phase_4j80/receipt.json") == j81.J80_RECEIPT_SHA
    )


@pytest.mark.parametrize(
    "driver",
    (
        "fit_hmds.py",
        "run_qm.py",
        "run_md.py",
        "run_converter.py",
        "evaluate_energy.py",
        "reserved_validation.py",
    ),
)
def test_forbidden_python_scientific_drivers_rejected(tmp_path, driver):
    marker = tmp_path / "ran"
    script = tmp_path / driver
    script.write_text(
        "from pathlib import Path; Path(" + repr(str(marker)) + ").touch()"
    )
    code = "import subprocess,sys; from scripts.phase_4j76_source_resolution import scope; scope(); subprocess.run([sys.executable,sys.argv[1]],check=True)"
    p = subprocess.run(
        [sys.executable, "-c", code, str(script)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode != 0 and not marker.exists()
    assert "scope" in p.stderr.lower()


def test_fake_backend_discovery_allowed_import_blocked(tmp_path):
    marker = tmp_path / "imported"
    (tmp_path / "psi4.py").write_text(
        "from pathlib import Path; Path(" + repr(str(marker)) + ").touch()"
    )
    code = "import sys,importlib.util; from scripts.phase_4j76_source_resolution import scope; scope(); sys.path.insert(0,sys.argv[1]); assert importlib.util.find_spec('psi4') is not None; import psi4"
    p = subprocess.run(
        [sys.executable, "-c", code, str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode != 0 and not marker.exists()
    assert "scope" in p.stderr.lower()


def test_unaudited_claim_diagnostics_never_verify(tmp_path):
    m = bundle(tmp_path)
    m["candidates"][0]["functional_form"] = "wrong"
    _, _, r = run(tmp_path, m, False)
    assert (
        "row_equation_mismatch:bond"
        in r["candidate_rows"][0]["declared_field_diagnostics"]
    )
    assert r["convention_comparison"]["status"] == "not_evaluated"
    assert r["verified_literature_rows"] == 0


@pytest.mark.parametrize(
    "change", ("native_identity", "contract_fields", "contract_value")
)
def test_malformed_metadata_fails_schema(tmp_path, change):
    m = bundle(tmp_path)
    if change == "native_identity":
        m["candidates"][0][change] = "not an object"
    elif change == "contract_fields":
        m["convention_contract"]["fields"] = []
    else:
        m["convention_contract"]["fields"]["torsion"]["value"] = "not an object"
    with pytest.raises(ValueError, match="schema_type"):
        pipeline.capture(m, tmp_path)


def test_missing_per_coefficient_units_block_verification(tmp_path):
    m = bundle(tmp_path)
    m["candidates"][0]["units"].pop("per_value")
    _, _, r = run(tmp_path, m)
    assert r["verified_literature_rows"] == 0
    assert (
        "explicit_per_value_units_required"
        in r["candidate_rows"][0]["blocking_reasons"]
    )


def test_intake_always_emits_a_status_receipt(tmp_path):
    root = tmp_path / "out"
    j81.emit(root, j81.prepare(REPO))
    r = j81.read(root / "intake-receipt.json")
    assert r["decision"] == "source_ingestion_blocked"
    assert r["verified_literature_rows"] == 0
    assert not r["scope"]["automatic_assignment"]


def test_untyped_literature_is_parsed_but_never_exact(tmp_path):
    m = bundle(tmp_path)
    m["candidates"][0]["atom_type_sequence"] = []
    snap, _, r = run(tmp_path, m)
    assert snap["artifacts"][0]["status"] == "source_parsed"
    assert r["verified_literature_rows"] == 0
    assert (
        "exact_source_atom_types_unresolved"
        in r["candidate_rows"][0]["blocking_reasons"]
    )


def test_unknown_forcefield_family_cannot_qualify(tmp_path):
    m = bundle(tmp_path)
    m["artifacts"][0]["source_family"] = "unknown"
    _, _, r = run(tmp_path, m)
    assert not r["source_gate"]["passed"]
    assert r["verified_literature_rows"] == 0
