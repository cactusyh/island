"""Network-free archive decoding/tampering tests using real synthetic Amber data."""

import copy
import importlib.util
import json
import shutil
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import pytest

pytest.importorskip("parmed")
from test_ambertools_integrity import _valid_result

from island.exceptions import EvaluationInputError, InvalidAmberImportResultError
from island.forcefields import import_amber_prmtop
from island.forcefields.ambertools.models import AmberToolsPreparationResult, digest

SCRIPT = Path(__file__).parents[1] / "scripts/validate_singlepoint_references.py"
spec = importlib.util.spec_from_file_location("singlepoint_archive", SCRIPT)
archive = importlib.util.module_from_spec(spec)
spec.loader.exec_module(archive)


@pytest.fixture
def synthetic_archive(tmp_path):
    directory = tmp_path / "synthetic"
    directory.mkdir()
    system, template = _valid_result(directory)
    record = copy.deepcopy(dict(template.record))
    record.pop("imported_result_signature")
    shutil.copyfile(directory / "source.prmtop", directory / "result.prmtop")
    shutil.copyfile(directory / "source.rst7", directory / "result.rst7")
    for name in ("typed.frcmod", "leap.in", "leap.log", "lineage.json"):
        (directory / name).write_text("synthetic offline test artifact\n")
    mol2 = (
        "@<TRIPOS>ATOM\n"
        + "\n".join(
            f"{index} {name} 0 0 0 C.3 1 SYN 0.0"
            for index, name in enumerate(record["lineage"]["input_name_to_site_id"], 1)
        )
        + "\n@<TRIPOS>BOND\n"
    )
    for name in ("input.mol2", "typed.mol2"):
        (directory / name).write_text(mol2)
    (directory / "input_system.json").write_text(json.dumps(system.to_dict()))
    record["artifact_sha256"] = {
        name: sha256((directory / name).read_bytes()).hexdigest()
        for name in record["artifact_sha256"]
    }
    record["input_mol2_sha256"] = sha256(
        (directory / "input.mol2").read_bytes()
    ).hexdigest()
    record["input_lineage_sha256"] = sha256(
        (directory / "lineage.json").read_bytes()
    ).hexdigest()
    imported = import_amber_prmtop(
        system,
        directory / "result.prmtop",
        dict(template.imported_result.mapping),
        source="AmberTools gaff2 generated prmtop",
        force_field="gaff2",
        charge_method="provided",
    )
    imported = replace(
        imported,
        provenance={**dict(imported.provenance), "ambertools_preparation": record},
        result_signature="",
    )
    imported = replace(imported, result_signature=imported.content_signature())
    outer = {**record, "imported_result_signature": imported.result_signature}
    original = AmberToolsPreparationResult(imported, outer, digest(outer))
    original.validate_integrity(system)
    row = {
        "record": outer,
        "record_signature": original.record_signature,
        "imported_signature": imported.result_signature,
        "source_sha256": imported.source_sha256,
        "artifact_subdirectory": directory.name,
        "retained_file_sha256": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in directory.iterdir()
            if p.is_file()
        },
    }
    return tmp_path, json.loads(json.dumps(row)), original


def test_roundtrip_restores_integer_keys_and_historical_signatures(synthetic_archive):
    root, row, original = synthetic_archive
    before = copy.deepcopy(row)
    system, loaded, directory = archive.load_archived_preparation(root, row)
    assert row == before
    assert loaded.record_signature == original.record_signature
    assert (
        loaded.imported_result.result_signature
        == original.imported_result.result_signature
    )
    assert dict(loaded.record) == dict(original.record)
    assert all(
        type(key) is int for key in loaded.record["lineage"]["prmtop_index_to_site_id"]
    )
    assert all(type(key) is str for key in loaded.record["input_coordinates_angstrom"])
    assert (
        archive.charge_report(system, loaded, directory)["final_charge_residual_e"] == 0
    )
    # Archive can move; the signed historical paths must not be rewritten.
    relocated = root / "relocated"
    shutil.copytree(directory, relocated / directory.name)
    assert (
        archive.load_archived_preparation(relocated, row)[1].record_signature
        == original.record_signature
    )


@pytest.mark.parametrize(
    "field,value",
    [
        ("charge_validation_tolerance_e", 0.002),
        ("requested_force_field", "gaff"),
    ],
)
def test_changed_signed_preparation_is_rejected(synthetic_archive, field, value):
    root, row, _ = synthetic_archive
    row["record"][field] = value
    with pytest.raises(EvaluationInputError, match="record signature"):
        archive.load_archived_preparation(root, row)


def test_rehashed_preparation_cannot_replace_historical_import_signature(
    synthetic_archive,
):
    root, row, _ = synthetic_archive
    row["record"]["warnings"] = ["new contradictory history"]
    row["record_signature"] = digest(archive.restore_preparation_record(row["record"]))
    with pytest.raises(InvalidAmberImportResultError, match="signature"):
        archive.load_archived_preparation(root, row)


@pytest.mark.parametrize("name", ["result.prmtop", "typed.mol2", "input_system.json"])
def test_changed_archive_artifact_is_rejected(synthetic_archive, name):
    root, row, _ = synthetic_archive
    (root / row["artifact_subdirectory"] / name).write_text("tampered")
    with pytest.raises(EvaluationInputError, match="checksum"):
        archive.load_archived_preparation(root, row)


def test_noncanonical_mapping_key_rejected(synthetic_archive):
    root, row, _ = synthetic_archive
    mapping = row["record"]["lineage"]["prmtop_index_to_site_id"]
    mapping["00"] = mapping.pop("0")
    with pytest.raises(EvaluationInputError, match="Noncanonical"):
        archive.load_archived_preparation(root, row)
