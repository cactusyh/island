"""J93 metadata-gate controls; synthetic files are not scientific evidence."""

import shutil
from pathlib import Path

import pytest

from island.forcefields.opls_providers.intake93 import PRIOR_HASHES, intake
from scripts.phase_4j73_source_intake import pin, read, signed, write
from scripts.phase_4j93_metadata_gate import emit, replay


def source(root):
    root.mkdir()
    (root / "LIG.lammps.lmp").write_text("pinned synthetic data control\n")
    (root / "LIG.lammps.in").write_text("pinned synthetic settings control\n")
    return {
        "source_dir": root,
        "data_sha256": pin(root / "LIG.lammps.lmp"),
        "settings_sha256": pin(root / "LIG.lammps.in"),
    }


def test_missing_structure_and_metadata_stop_before_graph_inference(tmp_path):
    result = intake(**source(tmp_path / "external")).report
    assert result["decision"] == "parameter_bundle_incomplete"
    assert result["exact_graph_identity"] is None
    assert result["formal_charge"] is None
    assert result["terminal_caps"] is None
    assert result["pair14_count"] == result["exclusion_count"] == 0
    assert result["source_files_parsed_this_step"] == 0
    assert "original_input_structure_missing" in result["blockers"]
    assert "terminal_cap_information_missing" in result["blockers"]
    assert result["J92_receipt_preserved"] is True
    assert result["J92_artifact_hashes"] == PRIOR_HASHES


def test_structure_presence_does_not_fill_missing_request_or_caps(tmp_path):
    args = source(tmp_path / "external")
    (args["source_dir"] / "original.sdf").write_text("unparsed control\n")
    result = intake(**args).report
    assert result["original_structure_files"] == ["original.sdf"]
    assert result["exact_graph_identity"] is None
    assert result["formal_charge"] is None
    assert result["terminal_caps"] is None
    assert "request.json_missing" in result["blockers"]


def test_source_hash_mutation_is_blocked(tmp_path):
    args = source(tmp_path / "external")
    (args["source_dir"] / "LIG.lammps.lmp").write_text("tampered\n")
    result = intake(**args).report
    assert result["decision"] == "parameter_bundle_incomplete"
    assert result["blockers"] == ["source_SHA256_mutation:LIG.lammps.lmp"]


def test_prior_receipt_mutation_is_blocked(tmp_path):
    args = source(tmp_path / "external")
    prior = tmp_path / "prior"
    shutil.copytree(Path("docs/evidence/phase_4j92_real_pe5"), prior)
    receipt = read(prior / "qualification-receipt.json")
    receipt["decision"] = "audited_exact_graph"
    write(
        prior / "qualification-receipt.json",
        signed({key: value for key, value in receipt.items() if key != "identity"}),
    )
    result = intake(**args, prior_root=prior).report
    assert result["decision"] == "parameter_bundle_incomplete"
    assert "J92_historical_artifact_mutation" in result["blockers"][0]


def test_deterministic_and_relocated_replay(tmp_path):
    args = source(tmp_path / "external")
    output = tmp_path / "output"
    result = emit(Path.cwd(), output, **args)
    snapshot = pin(output / "captured-input.json")
    assert replay(Path.cwd(), output, snapshot, **args) == result
    relocated = tmp_path / "relocated"
    shutil.copytree(args["source_dir"], relocated)
    relocated_args = {**args, "source_dir": relocated}
    assert replay(Path.cwd(), output, snapshot, **relocated_args) == result


def test_receipt_and_source_tamper_rejection(tmp_path):
    args = source(tmp_path / "external")
    output = tmp_path / "output"
    emit(Path.cwd(), output, **args)
    snapshot = pin(output / "captured-input.json")
    (args["source_dir"] / "LIG.lammps.in").write_text("tampered\n")
    with pytest.raises(ValueError, match="semantic_mutation"):
        replay(Path.cwd(), output, snapshot, **args)
    receipt = read(output / "receipt.json")
    receipt["scope"]["runtime_registration"] = True
    write(
        output / "receipt.json",
        signed({key: value for key, value in receipt.items() if key != "identity"}),
    )
    with pytest.raises(ValueError):
        replay(Path.cwd(), output, snapshot, **args)


def test_source_inside_worktree_is_rejected():
    result = intake(
        Path.cwd() / "bundle",
        data_sha256="0" * 64,
        settings_sha256="0" * 64,
    ).report
    assert "licensed_files_must_remain_outside_Git_worktree" in result["blockers"]


def test_real_external_metadata_gate_when_available():
    root = Path(
        "/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-opls-bundles/pe5"
    )
    if not root.is_dir():
        pytest.skip("External PE5 source directory is unavailable.")
    result = intake(
        root,
        data_sha256="19aa11e34253a6838319cdbf183aacae35ace320563e3e9bf57c8f9b8e4a79f6",
        settings_sha256="89203559d4717fa50d1b1681c62db2fe0d2ddc264a7d18257584cf24e4cb3a63",
    ).report
    assert result["decision"] == "parameter_bundle_incomplete"
    assert result["provider"] == "LigParGen"
    assert result["provider_version"] == "3.0.1"
    assert result["exact_graph_identity"] is None
    assert result["formal_charge"] is None
    assert result["terminal_caps"] is None
    assert result["J92_receipt_preserved"] is True
