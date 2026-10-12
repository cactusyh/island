import base64
import json
from hashlib import sha256

import numpy as np
import pytest

from island import periodic_md, trajectory_store
from island.analysis import (
    TrajectoryFrame,
    TrajectoryMetadata,
    VerifiedTrajectory,
    calculate_observables,
    load_analysis_record,
    make_analysis_record,
    write_analysis_record,
)
from island.exceptions import AnalysisError
from island.trajectory_store import (
    TrajectoryStoreConfig,
    TrajectoryStoreWriter,
    load_restart,
)


def _metadata(count=4):
    digest = "a" * 64
    return TrajectoryMetadata(
        final_graph_identity=digest,
        packing_identity="b" * 64,
        force_field_family="GAFF2",
        force_field_source_identity="c" * 64,
        typed_graph_identity="d" * 64,
        charge_identity="e" * 64,
        workflow_identity="f" * 64,
        checkpoint_identity="1" * 64,
        atom_ids=tuple(range(1, count + 1)),
        masses_dalton=tuple([1.0] * count),
        atom_mapping={index: index - 1 for index in range(1, count + 1)},
    )


def _trajectory(coordinates, *, groups, bonds):
    frames = tuple(
        TrajectoryFrame(np.asarray(coordinates, dtype=float), np.diag([10.0, 10.0, 10.0]), step, float(step))
        for step in range(2)
    )
    return VerifiedTrajectory(_metadata(len(coordinates)), frames, "trajectory", {"trajectory": "9" * 64}, tuple(bonds), tuple(groups))


def test_boundary_crossing_uses_unwrapped_chain_geometry():
    trajectory = _trajectory([[9.0, 1.0, 1.0], [1.0, 1.0, 1.0]], groups=((1, 2),), bonds=((1, 2),))
    values = calculate_observables(trajectory, ("end_to_end_distance", "re_over_rg"), block_size=1)
    assert values["end_to_end_distance"]["value"] == pytest.approx(2.0)
    assert values["re_over_rg"]["value"] == pytest.approx(2.0)


def test_multiple_chains_are_aggregated_per_chain():
    trajectory = _trajectory(
        [[9.0, 1.0, 1.0], [1.0, 1.0, 1.0], [5.0, 1.0, 1.0], [6.0, 1.0, 1.0]],
        groups=((1, 2), (3, 4)),
        bonds=((1, 2), (3, 4)),
    )
    values = calculate_observables(trajectory, ("re_over_rg",), block_size=1)
    assert values["re_over_rg"]["value"] == pytest.approx(2.0, rel=1e-6)
    assert len(values["re_over_rg"]["per_chain"])


def test_branched_periodic_unwrapping_is_blocked():
    trajectory = _trajectory(
        [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [3.0, 1.0, 1.0], [2.0, 2.0, 1.0]],
        groups=((1, 2, 3, 4),),
        bonds=((1, 2), (2, 3), (2, 4)),
    )
    values = calculate_observables(trajectory, ("re_over_rg",), block_size=1)
    assert values["re_over_rg"]["status"] == "blocked"
    assert "branched" in values["re_over_rg"]["failure"]


def test_nested_binding_tamper_requires_trusted_identity(tmp_path):
    trajectory = _trajectory([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]], groups=((1, 2),), bonds=((1, 2),))
    record = make_analysis_record(trajectory, calculate_observables(trajectory, ("density",), block_size=1))
    path = write_analysis_record(record, tmp_path / "record.json")
    payload = json.loads(path.read_text())
    payload["bindings"]["final_graph_identity"] = "0" * 64
    body = dict(payload)
    body.pop("identity")
    from hashlib import sha256
    payload["identity"] = sha256((json.dumps(body, sort_keys=True, separators=(",", ":")) + "\n").encode()).hexdigest()
    path.write_text(json.dumps(payload))
    with pytest.raises(AnalysisError, match="final_graph_identity"):
        load_analysis_record(path, expected_bindings=record.payload["bindings"])


def _ids():
    return {name: f"{index:064x}" for index, name in enumerate(("workflow_identity", "backend_identity", "graph_identity", "packing_identity"), 1)}


def test_checkpoint_interval_is_enforced(tmp_path):
    writer = TrajectoryStoreWriter(tmp_path / "store", identities=_ids(), config=TrajectoryStoreConfig(checkpoint_interval=4))
    with pytest.raises(AnalysisError, match="checkpoint interval"):
        writer.publish_restart(step=2, time_ps=0.0, checkpoint=b"state", openmm={}, configuration_identity="a" * 64)


def test_manifest_accepted_restart_wins_and_extra_restart_is_rejected(tmp_path):
    root = tmp_path / "store"
    writer = TrajectoryStoreWriter(root, identities=_ids(), config=TrajectoryStoreConfig(checkpoint_interval=4))
    first = writer.publish_restart(step=0, time_ps=0.0, checkpoint=b"zero", openmm={}, configuration_identity="a" * 64)
    writer.publish_restart(step=4, time_ps=1.0, checkpoint=b"four", openmm={}, configuration_identity="a" * 64)
    (root / "restart-999999999999.json").write_text(json.dumps(first.payload if hasattr(first, "payload") else first))
    with pytest.raises(AnalysisError, match="extra restart"):
        load_restart(root)


def test_failed_publication_preserves_previous_accepted_boundary(tmp_path, monkeypatch):
    root = tmp_path / "store"
    writer = TrajectoryStoreWriter(root, identities=_ids(), config=TrajectoryStoreConfig(checkpoint_interval=4))
    first = writer.publish_restart(step=0, time_ps=0.0, checkpoint=b"zero", openmm={}, configuration_identity="a" * 64)
    writer.append_trajectory_frame({"step": 1, "atom_ids": [1], "positions": [[0.0, 0.0, 0.0]]})
    original = trajectory_store._atomic_write

    def fail_manifest(path, raw):
        if path.name == "manifest.json" and path.parent == root:
            raise OSError("simulated manifest publication failure")
        return original(path, raw)

    monkeypatch.setattr(trajectory_store, "_atomic_write", fail_manifest)
    with pytest.raises(OSError, match="simulated"):
        writer.publish_restart(step=4, time_ps=1.0, checkpoint=b"four", openmm={}, configuration_identity="a" * 64)
    assert writer.accepted_restart["identity"] == first["identity"]
    assert load_restart(root)["identity"] == first["identity"]


def test_tampered_committed_trajectory_chunk_is_rejected(tmp_path):
    root = tmp_path / "store"
    writer = TrajectoryStoreWriter(root, identities=_ids(), config=TrajectoryStoreConfig(trajectory_interval=1, checkpoint_interval=1))
    writer.append_trajectory_frame({"step": 0, "atom_ids": [1], "positions": [[0.0, 0.0, 0.0]]})
    writer.publish_restart(step=0, time_ps=0.0, checkpoint=b"zero", openmm={}, configuration_identity="a" * 64)
    chunk = next(root.glob("trajectory-*.jsonl"))
    chunk.write_text(chunk.read_text().replace("0.0", "1.0", 1))
    with pytest.raises(AnalysisError, match="chunk hash"):
        load_restart(root)


def test_opt_in_j23_workflow_export_to_store(tmp_path, monkeypatch):
    def checkpoint(value):
        raw = value.encode()
        return base64.b64encode(raw).decode(), sha256(raw).hexdigest()

    first, first_hash = checkpoint("zero")
    second, second_hash = checkpoint("one")
    manifest = {
        "identity": "1" * 64,
        "backend": {"identity": "2" * 64, "final_graph_identity": "3" * 64},
        "packing_identity": "4" * 64,
        "atom_ids": [1],
        "engine": {"openmm_version": "8.6.1", "platform": "Reference"},
        "config": {"identity": "5" * 64, "temperature_kelvin": 300.0, "pressure_bar": None},
    }
    frames = []
    for step, encoded, checksum in ((0, first, first_hash), (1, second, second_hash)):
        frames.append({"step": step, "time_ps": float(step), "positions_nm": [[0.0, 0.0, 0.0]], "box_vectors_nm": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "temperature_kelvin": 300.0, "pressure_bar": None, "potential_energy_kj_mol": 0.0, "checkpoint_base64": encoded, "checkpoint_sha256": checksum})
    fake = periodic_md.PeriodicMDWorkflow(json.dumps({"manifest": manifest, "frames": frames}))
    monkeypatch.setattr(periodic_md, "load_periodic_md_workflow", lambda path: fake)
    record = periodic_md.export_periodic_md_workflow_to_trajectory_store(tmp_path / "workflow", tmp_path / "store", trajectory_interval=2, thermodynamic_interval=2, checkpoint_interval=1)
    assert record["step"] == 1
    assert load_restart(tmp_path / "store")["step"] == 1
