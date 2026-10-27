"""Bounded audit controls; parser sentinels are not physical parameter assets."""

import copy
import json
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace
from hashlib import sha256

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields.external_pdms import digest
from island.forcefields.literature_audit import (
    BUNDLE_SCHEMA,
    ROW_SCHEMA,
    STAGES,
    LiteratureBlockedError,
    audit_parameter_bundle,
    fixture_requirements,
    inventory_literature,
    validate_evidence_row,
)
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j44 import execute, reconstruct


@pytest.fixture
def control(tmp_path):
    root = tmp_path / "references"
    root.mkdir()
    raw = b"audit_parser_only 1.000\naudit_parser_only 2.000\n"
    (root / "sentinel.txt").write_bytes(raw)
    row = {
        "schema": ROW_SCHEMA,
        "kind": "atom_typing",
        "model_family": "external:opls_aa",
        "functional_form": "opls_aa",
        "source": {
            "file": "sentinel.txt",
            "sha256": sha256(raw).hexdigest(),
            "doi": "parser-control-not-a-publication",
            "citation": "Nonphysical test control",
            "location": "line 1; no literature table",
        },
        "row": {
            "family": "atom_type",
            "types": ["audit_parser_only"],
            "values": ["1.000"],
            "units": ["dalton"],
            "provenance": "No physical parameter or chemical authority",
            "citations": [
                {
                    "file": "sentinel.txt",
                    "line": 1,
                    "raw": "audit_parser_only 1.000",
                    "columns": [1],
                }
            ],
        },
        "conventions": {
            "mixing_rule": "geometric_sigma_and_epsilon",
            "one_four_scaling": "lj=0.5;coulomb=0.5",
            "charge_method": "no physical charges supplied",
            "torsion_convention": "opls_fourier_v1",
            "state_point": "not_applicable_parser_test",
        },
        "conversion": {
            "policy": "identity_v1",
            "original_values": ["1.000"],
            "converted_values": ["1.000"],
            "original_units": ["dalton"],
            "converted_units": ["dalton"],
            "code": "island.identity_decimal_v1",
            "version": 1,
        },
    }
    row["evidence_checksum"] = digest(row)
    return root, row


def resign(row):
    row["evidence_checksum"] = digest(
        {k: v for k, v in row.items() if k != "evidence_checksum"}
    )
    return row


def test_recursive_inventory_includes_hidden_binary_and_symlink_files(tmp_path):
    root = tmp_path / "refs"
    (root / ".hidden").mkdir(parents=True)
    (root / ".hidden/table.tsv").write_bytes(b"not a parameter table\n")
    (root / "paper.pdf").write_bytes(b"%PDF parser signature only")
    (root / "source.zip").write_bytes(b"PK\x03\x04 parser signature only")
    (root / "external-link").symlink_to(tmp_path / "unavailable-external-target")
    audit = inventory_literature(root)
    records = {f["filename"]: f for f in audit.payload["files"]}
    assert set(records) == {
        ".hidden/table.tsv",
        "paper.pdf",
        "source.zip",
        "external-link",
    }
    for name in (".hidden/table.tsv", "paper.pdf", "source.zip"):
        assert records[name]["sha256"] == sha256((root / name).read_bytes()).hexdigest()
        assert records[name]["size"] == (root / name).stat().st_size
    assert records["paper.pdf"]["file_type"] == "pdf"
    assert records["source.zip"]["file_type"] == "zip_archive"
    assert records["external-link"]["hash_basis"] == "symlink_text"
    assert all(
        f["bibliography"]["doi"] is None and f["license"] == "unknown"
        for f in records.values()
    )
    assert all(
        f["numerical_parameter_rows"] == "not_independently_audited"
        for f in records.values()
    )
    assert audit.payload["stages"]["source_inventory_completed"]
    assert not audit.payload["stages"]["citation_audit_completed"]


@pytest.mark.parametrize("present", [True, False])
def test_empty_or_absent_directory_is_asset_unavailable(tmp_path, present):
    root = tmp_path / "refs"
    if present:
        root.mkdir()
    p = inventory_literature(root).payload
    assert p["files"] == [] and p["status"] == "asset_unavailable"
    assert p["stages"]["source_inventory_completed"] is present
    assert p["execution_profiles"] == [] and p["acceptance_gate"] == 1
    assert all(p["stages"][k] is False for k in STAGES[1:])


def test_inventory_is_immutable_and_source_mutation_rejected(control):
    root, _ = control
    audit = inventory_literature(root)
    audit.payload["files"].clear()
    assert audit.payload["files"]
    with pytest.raises(FrozenInstanceError):
        audit.json_text = "changed"
    (root / "sentinel.txt").write_text("changed source")
    with pytest.raises(ValueError, match="stale_or_mutated"):
        audit.validate_integrity(root)


def test_raw_token_reconstruction_has_no_citation_or_execution_authority(control):
    root, row = control
    evidence = validate_evidence_row(row, root=root, audit=inventory_literature(root))
    assert evidence["extracted"]["citations"][0]["decimal_values"] == ["1.000"]
    assert evidence["execution_authorized"] is False
    assert (
        evidence["citation_metadata_status"] == "supplied_not_independently_qualified"
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "source_hash",
        "source_file",
        "doi",
        "location",
        "wildcard",
        "kind",
        "functional_form",
        "family",
        "units",
        "wrong_mass_units",
        "mixing",
        "scaling",
        "torsion",
        "charge_method",
        "state_point",
        "original",
        "converted",
        "conversion_code",
        "conversion_version",
        "raw",
        "line",
        "columns",
        "inserted_zero",
        "checksum",
    ],
)
def test_mutated_evidence_rows_are_rejected(control, mutation):
    root, row = control
    if mutation in ("source_hash", "source_file", "doi", "location"):
        key = {"source_hash": "sha256", "source_file": "file"}.get(mutation, mutation)
        row["source"][key] = "0" * 64 if mutation == "source_hash" else ""
    elif mutation == "wildcard":
        row["row"]["types"] = ["*"]
    elif mutation == "kind":
        row["kind"] = "class2_cross_term"
    elif mutation == "functional_form":
        row["functional_form"] = "harmonic_substitute"
    elif mutation == "family":
        row["model_family"] = "native:OPLS-AA"
    elif mutation == "units":
        row["row"]["units"] = ["unresolved"]
        row["conversion"]["original_units"] = ["unresolved"]
        row["conversion"]["converted_units"] = ["unresolved"]
    elif mutation == "wrong_mass_units":
        row["row"]["units"] = ["kJ/mol"]
        row["conversion"]["original_units"] = ["kJ/mol"]
        row["conversion"]["converted_units"] = ["kJ/mol"]
    elif mutation in ("mixing", "scaling", "torsion", "charge_method", "state_point"):
        key = {
            "mixing": "mixing_rule",
            "scaling": "one_four_scaling",
            "torsion": "torsion_convention",
        }.get(mutation, mutation)
        row["conventions"][key] = ""
    elif mutation in ("original", "converted"):
        row["conversion"][mutation + "_values"] = ["0"]
    elif mutation in ("conversion_code", "conversion_version"):
        row["conversion"][mutation.removeprefix("conversion_")] = "forged"
    elif mutation in ("raw", "line", "columns"):
        row["row"]["citations"][0][mutation] = {
            "raw": "invented",
            "line": 0,
            "columns": [99],
        }[mutation]
    elif mutation == "inserted_zero":
        row["row"]["values"].append("0")
        row["conversion"]["original_values"] = row["row"]["values"][:]
        row["conversion"]["converted_values"] = row["row"]["values"][:]
        row["row"]["units"].append("dalton")
        row["conversion"]["original_units"] = row["row"]["units"][:]
        row["conversion"]["converted_units"] = row["row"]["units"][:]
    if mutation != "checksum":
        resign(row)
    else:
        row["evidence_checksum"] = "0" * 64
    with pytest.raises(ValueError):
        validate_evidence_row(row, root=root, audit=inventory_literature(root))


@pytest.mark.parametrize("form", ["class2", "opls_aa"])
def test_complete_graph_inventory_does_not_imply_complete_parameters(tmp_path, form):
    root = tmp_path / "empty"
    root.mkdir()
    system = from_smiles("CCO", random_seed=2026)
    audit = inventory_literature(root)
    result = audit_parameter_bundle(system, root=root, audit=audit, form=form)
    p = result.payload
    assert p["missing_rows"] == p["requirements"]["requests"]
    assert p["unresolved_dependencies"]
    if form == "class2":
        assert any(k.startswith("bond-angle:") for k in p["missing_rows"])
    else:
        assert any(k.startswith("special_pair:") for k in p["missing_rows"])
        assert "policy:torsion_convention" in p["missing_rows"]
    with pytest.raises(LiteratureBlockedError) as caught:
        result.require_executable()
    assert caught.value.result is result and caught.value.status == "asset_unavailable"


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_row",
        "stale_graph",
        "wrong_order",
        "conflict",
        "mixed_family",
        "mixed_paper",
        "bad_pin",
    ],
)
def test_incomplete_candidate_bundles_cannot_hide_missing_dependencies(
    control, mutation
):
    root, row = control
    form = "class2" if mutation == "mixed_family" else "opls_aa"
    system = from_smiles("CC", random_seed=2026)
    required = fixture_requirements(system, form)
    data = {
        "schema": BUNDLE_SCHEMA,
        "graph_identity": required["graph_identity"],
        "atom_order": required["atom_order"],
        "rows": {},
        "assignments": {},
    }
    if mutation == "missing_row":
        data["assignments"] = {next(iter(required["requests"])): "absent"}
    elif mutation == "stale_graph":
        data["graph_identity"] = "0" * 64
    elif mutation == "wrong_order":
        data["atom_order"] = required["atom_order"][::-1]
    elif mutation in ("conflict", "mixed_family", "mixed_paper"):
        second = copy.deepcopy(row)
        if mutation == "conflict":
            second["row"]["values"] = ["2.000"]
            second["row"]["citations"][0].update(line=2, raw="audit_parser_only 2.000")
            second["conversion"]["original_values"] = ["2.000"]
            second["conversion"]["converted_values"] = ["2.000"]
        elif mutation == "mixed_paper":
            second["source"]["doi"] = "other-paper-not-authorized"
        else:
            row.update(functional_form="class2", model_family="external:pcff")
            second.update(functional_form="class2", model_family="external:compass")
        data["rows"] = {"a": resign(row), "b": resign(second)}
    path = root.parent / "candidate.json"
    path.write_text(json.dumps(data))
    result = audit_parameter_bundle(
        system,
        root=root,
        audit=inventory_literature(root),
        form=form,
        bundle=path,
        bundle_sha256="0" * 64
        if mutation == "bad_pin"
        else sha256(path.read_bytes()).hexdigest(),
    )
    assert result.payload["diagnostics"]
    assert result.payload["missing_rows"] == required["requests"]
    assert result.payload["execution_profile"] is None


@pytest.mark.parametrize(
    "mutation", ["stages", "profile", "source", "topology", "evidence"]
)
def test_reconstruction_rejects_tampered_study(tmp_path, mutation):
    refs = tmp_path / "refs"
    refs.mkdir()
    study = tmp_path / "study"
    execute(study, refs)
    if mutation in ("stages", "profile"):
        report = read(study / "report.json")
        if mutation == "stages":
            report["stages"]["scientific_validation_complete"] = True
        else:
            report["execution_profiles"] = ["forged"]
        write(study / "report.json", report)
    elif mutation == "source":
        (refs / "new.txt").write_text("new source")
    elif mutation == "topology":
        path = study / "siosi_oligomer/system.json"
        data = read(path)
        data["sites"][0]["metadata"]["chain_id"] = "mutated"
        write(path, data)
    else:
        write(study / "extracted-parameter-rows.json", {"rows": {"invented": "0"}})
    with pytest.raises((ValueError, KeyError)):
        reconstruct(study, refs)


def test_relocation_fresh_process_and_nonzero_acceptance(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    study = tmp_path / "study"
    child = subprocess.run(
        [
            sys.executable,
            "scripts/phase_4j44.py",
            "--references",
            str(refs),
            "--output",
            str(study),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 1, child.stderr
    assert json.loads(child.stdout)["source_files"] == 0
    relocated, new_refs = tmp_path / "relocated", tmp_path / "relocated-refs"
    shutil.copytree(study, relocated)
    shutil.copytree(refs, new_refs)
    child = subprocess.run(
        [
            sys.executable,
            "scripts/phase_4j44.py",
            "--references",
            str(new_refs),
            "--reconstruct",
            str(relocated),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr
    assert json.loads(child.stdout) == reconstruct(study, refs)
    for path in study.rglob("*"):
        if path.is_file():
            assert (
                path.read_bytes() == (relocated / path.relative_to(study)).read_bytes()
            )


def test_promoted_inventory_cannot_be_constructed(tmp_path):
    audit = inventory_literature(tmp_path)
    p = audit.payload
    p["stages"]["typing_complete"] = True
    with pytest.raises(ValueError, match="inventory_is_not_qualification"):
        replace(audit, json_text=pack(p))


def test_reviewed_supplied_inventory_and_partial_parameter_evidence():
    from pathlib import Path

    from island.forcefields.literature_audit import LiteratureSourceReview

    path = (
        Path(__file__).resolve().parents[1]
        / "docs/evidence/phase_4j44/source-audit.json"
    )
    review = LiteratureSourceReview(path.read_text())
    data = review.payload
    assert len(data["files"]) == 52
    assert len(data["extracted_row_candidates"]) == 136
    assert data["execution_profiles"] == data["accepted_parameter_rows"] == []
    assert {r["model_scope"] for r in data["extracted_row_candidates"]} >= {
        "COMPASS_Table_S3",
        "OPLS_Table_S4",
        "improved_UA_Table_S1",
    }
    assert any(
        "equilibrium dependency" in f["interaction"]
        for f in data["dependency_findings"]
    )
    assert any(
        "Buckingham" in f.get("diagnostic", "") for f in data["dependency_findings"]
    )
    assert any("Si-H" in f["interaction"] for f in data["dependency_findings"])
    assert any(
        f["bibliography"]["year"] == 1997 and "1999" in f["filename"]
        for f in data["files"]
    )
    # Records remain immutable copies, not mutable scientific authority.
    data["execution_profiles"].append("forged")
    assert review.payload["execution_profiles"] == []


@pytest.mark.parametrize(
    "field",
    [
        "source_sha256",
        "original_values",
        "raw_text",
        "model_scope",
        "evidence_checksum",
    ],
)
def test_review_rejects_extraction_mutation(field):
    from pathlib import Path

    from island.forcefields.literature_audit import LiteratureSourceReview

    path = (
        Path(__file__).resolve().parents[1]
        / "docs/evidence/phase_4j44/source-audit.json"
    )
    data = json.loads(path.read_text())
    data["extracted_row_candidates"][0][field] = "mutated"
    with pytest.raises(ValueError):
        LiteratureSourceReview(json.dumps(data))


def test_review_inventory_is_hash_bound(control):
    from pathlib import Path

    from island.forcefields.literature_audit import LiteratureSourceReview

    root, _ = control
    review = LiteratureSourceReview(
        (
            Path(__file__).resolve().parents[1]
            / "docs/evidence/phase_4j44/source-audit.json"
        ).read_text()
    )
    with pytest.raises(ValueError, match="review_inventory_changed"):
        review.validate_inventory(inventory_literature(root))


def test_atom_labels_must_be_present_in_cited_raw_row(control):
    root, row = control
    row["row"]["types"] = ["uncited_chemical_type"]
    with pytest.raises(ValueError, match="atom_labels_not_reconstructed"):
        validate_evidence_row(resign(row), root=root, audit=inventory_literature(root))
