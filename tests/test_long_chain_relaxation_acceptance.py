"""Software guards for retained-input acceptance tooling; no charge generation."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from island.workflows import storage
from scripts import validate_long_chain_relaxation as acceptance


def test_missing_retained_inputs_write_failure_and_exit_nonzero(tmp_path):
    source = tmp_path / "missing-inputs"
    source.mkdir()
    output = tmp_path / "output"
    assert acceptance.main(["--source", str(source), "--output", str(output)]) == 1
    report = storage.read_json(output / "report.json")
    assert report["acceptance"] == "unmet"
    assert "WorkflowError" in report["failure"]
    assert "input-system.json" in report["failure"]
    assert report["retained_source_unchanged"]
    assert not (output / "workflow").exists()


def test_existing_output_and_source_are_preserved(tmp_path):
    path = tmp_path / "source"
    path.mkdir()
    saved = path / "historical.json"
    saved.write_bytes(b"unchanged")
    with pytest.raises(ValueError, match="separate"):
        acceptance.run(path, path / "experiment")
    existing = tmp_path / "existing-output"
    existing.mkdir()
    with pytest.raises(FileExistsError):
        acceptance.run(path, existing)
    assert saved.read_bytes() == b"unchanged"


def state(force=1.0):
    return SimpleNamespace(
        coordinates={9: (1.0, 2.0, 3.0)},
        velocities={9: (4.0, 5.0, 6.0)},
        evaluation=SimpleNamespace(
            forces={9: (force, 0.0, 0.0)},
            potential_energy=2.0,
            energy_components={"bond": 2.0},
        ),
        kinetic_energy=3.0,
        total_energy=5.0,
    )


def test_comparison_rejects_outside_predeclared_tolerance():
    assert not any(acceptance.compare_states(state(), state()).values())
    with pytest.raises(AssertionError):
        acceptance.compare_states(state(), state(force=1.001))


def test_reject_changed_artifact_before_reconstruction(tmp_path, monkeypatch):
    case = tmp_path / "peo20"
    artifacts = case / "artifacts"
    artifacts.mkdir(parents=True)
    (artifacts / "result.prmtop").write_bytes(b"corrupted")
    payload = {
        "record": {
            "artifact_dir": "/historical/artifacts",
            "artifact_sha256": {"result.prmtop": "bad"},
            "input_mol2_sha256": "unused",
            "input_lineage_sha256": "unused",
        }
    }
    storage.publish(case / "input-system.json", storage.json_bytes(storage.encode({})))
    storage.publish(
        case / "preparation.json", storage.json_bytes(storage.encode(payload))
    )
    monkeypatch.setattr(acceptance, "system_from", lambda _: object())
    monkeypatch.setattr(
        acceptance, "preparation_from", lambda *_: pytest.fail("must reject first")
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        acceptance.load_retained(Path(tmp_path))
