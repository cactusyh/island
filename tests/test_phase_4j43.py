"""Literature contract controls only: no PDMS coefficients or positive asset."""

import copy
import json
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace
from hashlib import sha256

import pytest

from island.charge_references.records import pack
from island.exceptions import IslandError
from island.forcefields import ExternalPDMSBlockedError, create_evaluator
from island.forcefields.external_pdms import ROWS_SCHEMA, UNITS, bundle_template
from island.forcefields.external_pdms_adapter import ExternalPDMSPreparation
from island.forcefields.pdms_literature import (
    CLASSIFICATIONS,
    SOURCE_SCHEMA,
    LiteratureIngestion,
    ingest_literature,
    inspect_source,
    reconstruct_rows,
)
from island.workflows.bundle import system_data
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j41_c5 import CASES, fixture
from scripts.phase_4j43 import execute, reconstruct


@pytest.fixture
def control(tmp_path):
    """A parser sentinel, explicitly no literature authority or force-field row."""
    raw = b"parser_sentinel 0\n"
    row = {
        "family": "parser_sentinel",
        "types": ["not_an_atom_type"],
        "values": ["0"],
        "units": ["dimensionless"],
        "provenance": "Contract control only; not a parameter or publication",
        "citations": [
            {"file": "raw.txt", "line": 1, "raw": raw.decode().strip(), "columns": [1]}
        ],
    }
    document = {
        "schema": ROWS_SCHEMA,
        "source_family": "external:contract_control",
        "rows": {"sentinel": row},
        "provenance": "No physical coefficients",
    }
    files = {
        "raw.txt": raw,
        "rows.json": json.dumps(document).encode(),
        "license.txt": b"Locally generated test control only; no literature asset.\n",
    }
    manifest = {
        "schema": SOURCE_SCHEMA,
        "bibliography": {
            "doi": "not-a-publication-contract-control",
            "title": "Parser control",
            "authors": ["Test control"],
            "year": 2026,
        },
        "classification": "incomplete_literature_source",
        "source_family": "external:contract_control",
        "model_name": "No physical model",
        "license": {
            "identifier": "test-control-only",
            "authorized_use": True,
            "redistribution": "allowed",
            "evidence_file": "license.txt",
        },
        "files": {
            n: {
                "sha256": sha256(b).hexdigest(),
                "family": "external:contract_control",
                "role": {
                    "raw.txt": "raw_table",
                    "rows.json": "row_table",
                    "license.txt": "license",
                }[n],
                "provenance": "Generated contract control",
            }
            for n, b in files.items()
        },
        "primary_source": "rows.json",
        "parameterization": {
            "method": "No parameterization; parser test",
            "quantum_chemistry": None,
            "condensed_phase_targets": [],
        },
        "functional_form": "class2",
        "units": copy.deepcopy(UNITS),
        "charge_method": "No charges",
        "atom_type_definitions": {"sentinel": "Not a physical atom"},
        "reported_validation_properties": [],
    }
    for name, data in files.items():
        (tmp_path / name).write_bytes(data)
    return tmp_path, manifest, document


def save_source(control):
    root, manifest, _ = control
    path = root / "manifest.json"
    path.write_text(json.dumps(manifest))
    return path, sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("case", CASES, ids=lambda c: c[0])
def test_no_assets_preserves_every_graph_bound_requirement(case, monkeypatch):
    def forbidden(*a, **kw):
        pytest.fail("Absent assets must not initiate file or numerical discovery")

    import island.forcefields.pdms_literature as module

    monkeypatch.setattr(module, "pinned_bytes", forbidden)
    system = fixture(case)
    result = ingest_literature(system)
    p = result.payload
    assert p["status"] == "asset_unavailable" and p["source_record"] is None
    assert p["missing_interactions"] == p["preparation"]["requirements"]["interactions"]
    assert p["unresolved_dependencies"]
    assert p["execution_profile"] is None and p["backend_executions"] == []
    preparation = ExternalPDMSPreparation(pack(p["preparation"]))
    with pytest.raises(ExternalPDMSBlockedError):
        create_evaluator(system, preparation)


def test_no_independent_pin_does_not_open_source(control, monkeypatch):
    from pathlib import Path

    path, _ = save_source(control)

    def forbidden(*a):
        pytest.fail("Unpinned file was opened")

    monkeypatch.setattr(Path, "read_bytes", forbidden)
    result = inspect_source(path, None).payload
    assert result["manifest"] is None
    assert result["diagnostics"][0]["context"] == "independent_sha256_required"


def test_exact_raw_decimal_reconstruction_is_not_a_pdms_parameter(control):
    path, pin = save_source(control)
    record = inspect_source(path, pin)
    assert record.payload["diagnostics"] == []
    citation = record.payload["extracted_rows"]["sentinel"]["citations"][0]
    assert citation == {
        "source_file": "raw.txt",
        "line": 1,
        "numeric_columns": [1],
        "raw_text": "parser_sentinel 0",
        "decimal_values": ["0"],
    }
    assert record.payload["classification"] == "incomplete_literature_source"
    record.payload["extracted_rows"].clear()
    assert record.payload["extracted_rows"]
    with pytest.raises(FrozenInstanceError):
        record.json_text = "changed"


@pytest.mark.parametrize(
    "mutation", ["raw", "line", "columns", "value", "inserted_zero", "pdf"]
)
def test_row_mutations_rejected_by_independent_reader(control, mutation):
    _, _, document = control
    row = document["rows"]["sentinel"]
    assets = {"raw.txt": b"parser_sentinel 0\n"}
    if mutation == "raw":
        row["citations"][0]["raw"] = "changed"
    elif mutation == "line":
        row["citations"][0]["line"] = 0
    elif mutation == "columns":
        row["citations"][0]["columns"] = [99]
    elif mutation == "value":
        row["values"] = ["1"]
    elif mutation == "inserted_zero":
        row["values"].append("0")
    else:
        assets["raw.txt"] = b"%PDF fake parser control"
    extracted, errors = reconstruct_rows(document, assets)
    assert extracted == {} and errors[0]["stage"] == "source_row_incomplete"


@pytest.mark.parametrize(
    "mutation,context",
    [
        ("hash", "asset_sha256_mismatch"),
        ("license", "license_or_authorization_missing"),
        ("redistribution", "license_or_authorization_missing"),
        ("doi", "bibliography_incomplete"),
        ("units", "functional_form_or_units_unsupported"),
        ("functional_form", "functional_form_or_units_unsupported"),
        ("family", "mixed_source_or_file_fields_invalid"),
        ("path", "Bundle paths must be checked relative paths"),
        ("pdf_only", "reconstructable_parameter_table_missing"),
        ("atom_types", "chemical_provenance_missing"),
        ("charge", "chemical_provenance_missing"),
        ("method", "parameterization_metadata_missing"),
        ("license_file", "license_evidence_missing"),
    ],
)
def test_source_mutations_keep_classified_missing_evidence(control, mutation, context):
    root, m, _ = control
    if mutation == "hash":
        m["files"]["raw.txt"]["sha256"] = "0" * 64
    elif mutation == "license":
        m["license"]["authorized_use"] = False
    elif mutation == "redistribution":
        m["license"]["redistribution"] = "prohibited"
    elif mutation == "doi":
        m["bibliography"]["doi"] = ""
    elif mutation == "units":
        m["units"]["energy"] = "eV"
    elif mutation == "functional_form":
        m["functional_form"] = "opls"
    elif mutation == "family":
        m["files"]["raw.txt"]["family"] = "COMPASS"
    elif mutation == "path":
        m["files"]["../escape"] = m["files"].pop("raw.txt")
    elif mutation == "pdf_only":
        m["primary_source"] = None
    elif mutation == "atom_types":
        m["atom_type_definitions"] = {}
    elif mutation == "charge":
        m["charge_method"] = ""
    elif mutation == "method":
        m["parameterization"]["method"] = ""
    else:
        m["license"]["evidence_file"] = "missing"
    path, pin = save_source(control)
    result = ingest_literature(fixture(CASES[0]), source=path, source_sha256=pin)
    p = result.payload
    assert p["status"] == "parameter_incomplete"
    assert any(d["context"] == context for d in p["diagnostics"])
    assert p["source_record"]["classification"] in CLASSIFICATIONS
    assert p["execution_profile"] is None and p["acceptance_gate"] == 1
    assert not (root / "trajectory").exists()


@pytest.mark.parametrize("raw", [b"{", b"[]", b"null", b"\xff"])
def test_malformed_manifest_is_classified(tmp_path, raw):
    path = tmp_path / "manifest.json"
    path.write_bytes(raw)
    result = ingest_literature(
        fixture(CASES[0]), source=path, source_sha256=sha256(raw).hexdigest()
    )
    assert result.payload["status"] == "parameter_incomplete"
    assert result.payload["diagnostics"]


def test_no_complete_bundle_from_parser_rows(control):
    root, _, _ = control
    path, pin = save_source(control)
    system = fixture(CASES[0])
    bundle = root / "incomplete-bundle.json"
    bundle.write_text(pack(bundle_template(system)))
    p = ingest_literature(
        system,
        source=path,
        source_sha256=pin,
        bundle=bundle,
        bundle_sha256=sha256(bundle.read_bytes()).hexdigest(),
    ).payload
    assert p["status"] == "parameter_incomplete" and p["data_status"] == "incomplete"
    assert p["bundle_validation"]["status"] == "parameter_incomplete"
    assert p["missing_interactions"] and p["unresolved_dependencies"]


@pytest.mark.parametrize(
    "mutation", ["readiness", "authority", "checksum", "graph", "rows"]
)
def test_mutated_relocated_study_cannot_promote_evidence(tmp_path, mutation):
    execute(tmp_path / "study")
    if mutation in ("readiness", "authority"):
        path = tmp_path / "study/report.json"
        data = read(path)
        data["acceptance_gate" if mutation == "readiness" else "execution_profiles"] = (
            0 if mutation == "readiness" else ["forged"]
        )
        write(path, data)
    elif mutation == "rows":
        write(
            tmp_path / "study/extracted-parameter-rows.json", {"rows": {"forged": "0"}}
        )
    else:
        path = tmp_path / "study" / CASES[0][0] / "ingestion.json"
        record = LiteratureIngestion(path.read_text())
        p = record.payload
        if mutation == "graph":
            p["preparation"]["requirements"]["graph_identity"] = "0" * 64
            path.write_text(pack(p))
        else:
            raw = json.loads(record.json_text)
            raw["payload"]["status"] = "prepared"
            path.write_text(json.dumps(raw))
    with pytest.raises((ValueError, IslandError)):
        reconstruct(tmp_path / "study")


def test_pinned_source_and_incomplete_bundle_relocate_without_discovery(control):
    root, _, _ = control
    path, pin = save_source(control)
    system = fixture(CASES[0])
    bundle = root / "incomplete-bundle.json"
    bundle.write_text(pack(bundle_template(system)))
    bundle_pin = sha256(bundle.read_bytes()).hexdigest()
    before = ingest_literature(
        system, source=path, source_sha256=pin, bundle=bundle, bundle_sha256=bundle_pin
    )
    target = root / "relocated-assets"
    target.mkdir()
    for name in (
        "manifest.json",
        "rows.json",
        "raw.txt",
        "license.txt",
        "incomplete-bundle.json",
    ):
        shutil.copyfile(root / name, target / name)
    after = ingest_literature(
        system,
        source=target / "manifest.json",
        source_sha256=pin,
        bundle=target / "incomplete-bundle.json",
        bundle_sha256=bundle_pin,
    )
    assert before == after
    write(target / "system.json", system_data(system))
    child = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys; from pathlib import Path; "
                "from scripts.phase_4j41_c2 import read; "
                "from island.workflows.bundle import system_from; "
                "from island.forcefields.pdms_literature import ingest_literature; "
                "p=Path(sys.argv[1]); "
                "sys.stdout.write(ingest_literature(system_from(read(p/'system.json')), "
                "source=p/'manifest.json',source_sha256=sys.argv[2], "
                "bundle=p/'incomplete-bundle.json',bundle_sha256=sys.argv[3]).json_text)"
            ),
            str(target),
            pin,
            bundle_pin,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr
    assert LiteratureIngestion(child.stdout) == before
    (target / "raw.txt").write_text("tampered citation")
    changed = ingest_literature(
        system,
        source=target / "manifest.json",
        source_sha256=pin,
        bundle=target / "incomplete-bundle.json",
        bundle_sha256=bundle_pin,
    )
    assert changed.payload["status"] == "parameter_incomplete"
    assert changed.payload["bundle_validation"] is None
    assert changed.payload["source_record"] != before.payload["source_record"]


@pytest.mark.parametrize("pin", [None, "0" * 64])
def test_bundle_pin_required_even_with_intact_source(control, pin):
    root, _, _ = control
    source, source_pin = save_source(control)
    system = fixture(CASES[0])
    bundle = root / "incomplete-bundle.json"
    bundle.write_text(pack(bundle_template(system)))
    p = ingest_literature(
        system,
        source=source,
        source_sha256=source_pin,
        bundle=bundle,
        bundle_sha256=pin,
    ).payload
    assert p["status"] == "parameter_incomplete" and p["bundle_validation"] is None
    assert p["diagnostics"][0]["stage"] == "literature_bundle_incomplete"


def test_result_cannot_grant_execution():
    result = ingest_literature(fixture(CASES[0]))
    p = result.payload
    p["execution_profile"] = "forged"
    with pytest.raises(ValueError, match="literature_cannot_promote_readiness"):
        replace(result, json_text=pack(p))


def test_relocation_and_fresh_process_have_no_numerical_authority(tmp_path):
    execute(tmp_path / "study")
    shutil.copytree(tmp_path / "study", tmp_path / "relocated")
    child = subprocess.run(
        [
            sys.executable,
            "scripts/phase_4j43.py",
            "--reconstruct",
            str(tmp_path / "relocated"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr
    assert json.loads(child.stdout) == reconstruct(tmp_path / "study")
    for p in (tmp_path / "study").rglob("*"):
        if p.is_file():
            assert (
                p.read_bytes()
                == (
                    tmp_path / "relocated" / p.relative_to(tmp_path / "study")
                ).read_bytes()
            )


def test_cli_acceptance_remains_nonzero(tmp_path):
    child = subprocess.run(
        [sys.executable, "scripts/phase_4j43.py", "--output", str(tmp_path / "study")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 1, child.stderr
    assert json.loads(child.stdout)["acceptance_gate"] == 1
