import pytest

from island.exceptions import AnalysisError
from island.trajectory_store import (
    TrajectoryStoreConfig,
    TrajectoryStoreWriter,
    load_restart,
    read_trajectory_store,
    relocate_trajectory_store,
)


def _ids():
    return {name: f"{index:064x}" for index, name in enumerate(("workflow_identity", "backend_identity", "graph_identity", "packing_identity"), 1)}


def test_chunked_store_separates_frames_thermodynamics_and_restart(tmp_path):
    writer = TrajectoryStoreWriter(tmp_path / "store", identities=_ids(), config=TrajectoryStoreConfig(trajectory_interval=2, thermodynamic_interval=2, checkpoint_interval=2))
    writer.append_trajectory_frame({"step": 1, "positions": [[0.0, 0.0, 0.0]]})
    writer.append_trajectory_frame({"step": 2, "positions": [[1.0, 0.0, 0.0]]})
    writer.append_thermodynamic_record({"step": 1, "target_temperature_kelvin": 300.0, "measured_temperature_kelvin": 299.5})
    writer.append_thermodynamic_record({"step": 2, "target_temperature_kelvin": 300.0, "measured_temperature_kelvin": 300.1})
    record = writer.publish_restart(step=2, time_ps=0.4, checkpoint=b"checkpoint", openmm={"version": "8.6.1", "platform": "Reference"}, configuration_identity="c" * 64)
    assert record["trajectory_chunks"]
    assert record["thermodynamic_chunks"]
    loaded = load_restart(tmp_path / "store", expected_identity=record["identity"])
    frames, thermo = read_trajectory_store(tmp_path / "store")
    assert next(iter(frames))["step"] == 1
    assert list(thermo)[1]["measured_temperature_kelvin"] == 300.1
    assert loaded["checkpoint_sha256"] == record["checkpoint_sha256"]


def test_failed_or_tampered_restart_does_not_validate(tmp_path):
    writer = TrajectoryStoreWriter(tmp_path / "store", identities=_ids())
    with pytest.raises(AnalysisError):
        writer.publish_restart(step=0, time_ps=0.0, checkpoint=b"", openmm={}, configuration_identity="d" * 64)
    writer.publish_restart(step=0, time_ps=0.0, checkpoint=b"checkpoint", openmm={}, configuration_identity="d" * 64)
    checkpoint = next((tmp_path / "store").glob("checkpoint-*.bin"))
    checkpoint.write_bytes(b"tampered")
    with pytest.raises(AnalysisError, match="checkpoint hash"):
        load_restart(tmp_path / "store")


def test_relocation_preserves_restart_identity(tmp_path):
    writer = TrajectoryStoreWriter(tmp_path / "store", identities=_ids())
    writer.publish_restart(step=0, time_ps=0.0, checkpoint=b"checkpoint", openmm={}, configuration_identity="e" * 64)
    destination = relocate_trajectory_store(tmp_path / "store", tmp_path / "moved")
    assert load_restart(destination)["step"] == 0


def test_store_rejects_nonmonotonic_frames(tmp_path):
    writer = TrajectoryStoreWriter(tmp_path / "store", identities=_ids())
    writer.append_trajectory_frame({"step": 2})
    with pytest.raises(AnalysisError, match="increase strictly"):
        writer.append_trajectory_frame({"step": 1})
