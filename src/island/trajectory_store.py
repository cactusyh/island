"""Bounded-memory, opt-in trajectory and restart persistence.

This store is additive to the short J23 in-memory workflow. It keeps trajectory
frames, thermodynamic records and accepted binary restart checkpoints in separate
files and never places a full stochastic checkpoint in a trajectory frame.
"""

from __future__ import annotations

import json
import os
import shutil
from collections.abc import Iterator
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

from island.exceptions import AnalysisError

SCHEMA = "island.trajectory-store.v1"


def _digest(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def _file_hash(path: Path) -> str:
    return _digest(path.read_bytes())


def _atomic_write(path: Path, raw: bytes) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_bytes(raw)
    with temporary.open("rb") as stream:
        os.fsync(stream.fileno())
    os.replace(temporary, path)


@dataclass(frozen=True)
class TrajectoryStoreConfig:
    trajectory_interval: int = 100
    thermodynamic_interval: int = 10
    checkpoint_interval: int = 1000

    def __post_init__(self) -> None:
        for name in ("trajectory_interval", "thermodynamic_interval", "checkpoint_interval"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise AnalysisError(f"{name} must be a positive integer")

    def to_dict(self) -> dict[str, int]:
        return {
            "trajectory_interval": self.trajectory_interval,
            "thermodynamic_interval": self.thermodynamic_interval,
            "checkpoint_interval": self.checkpoint_interval,
        }


class TrajectoryStoreWriter:
    """Incrementally write sampled records and accepted restart boundaries."""

    def __init__(self, root: str | Path, *, identities: dict[str, str], config: TrajectoryStoreConfig | None = None) -> None:
        self.root = Path(root)
        if self.root.exists() and any(self.root.iterdir()):
            raise AnalysisError(f"Trajectory store must be a new directory: {self.root}")
        self.root.mkdir(parents=True, exist_ok=True)
        required = ("workflow_identity", "backend_identity", "graph_identity", "packing_identity")
        if any(type(identities.get(name)) is not str for name in required):
            raise AnalysisError("Trajectory store identities are incomplete")
        self.identities = dict(identities)
        self.config = config or TrajectoryStoreConfig()
        self._frame_buffer: list[dict[str, Any]] = []
        self._thermo_buffer: list[dict[str, Any]] = []
        self._last_step = -1
        self._accepted_restart: dict[str, Any] | None = None
        _atomic_write(self.root / "manifest.json", _json_bytes({"schema": SCHEMA, "identities": self.identities, "config": self.config.to_dict(), "accepted_restart": None}))

    @property
    def accepted_restart(self) -> dict[str, Any] | None:
        return None if self._accepted_restart is None else dict(self._accepted_restart)

    def append_trajectory_frame(self, frame: dict[str, Any]) -> None:
        step = frame.get("step")
        if type(step) is not int or step <= self._last_step:
            raise AnalysisError("Trajectory frame steps must increase strictly")
        self._last_step = step
        self._frame_buffer.append(dict(frame))
        if len(self._frame_buffer) >= self.config.trajectory_interval:
            self.flush_trajectory()

    def append_thermodynamic_record(self, record: dict[str, Any]) -> None:
        if type(record) is not dict or type(record.get("step")) is not int:
            raise AnalysisError("Thermodynamic records require an integer step")
        self._thermo_buffer.append(dict(record))
        if len(self._thermo_buffer) >= self.config.thermodynamic_interval:
            self.flush_thermodynamics()

    def flush_trajectory(self) -> str | None:
        if not self._frame_buffer:
            return None
        index = len(list(self.root.glob("trajectory-*.jsonl")))
        path = self.root / f"trajectory-{index:08d}.jsonl"
        raw = b"".join(_json_bytes(item) for item in self._frame_buffer)
        _atomic_write(path, raw)
        self._frame_buffer.clear()
        return _file_hash(path)

    def flush_thermodynamics(self) -> str | None:
        if not self._thermo_buffer:
            return None
        index = len(list(self.root.glob("thermodynamics-*.jsonl")))
        path = self.root / f"thermodynamics-{index:08d}.jsonl"
        raw = b"".join(_json_bytes(item) for item in self._thermo_buffer)
        _atomic_write(path, raw)
        self._thermo_buffer.clear()
        return _file_hash(path)

    def publish_restart(self, *, step: int, time_ps: float, checkpoint: bytes, openmm: dict[str, Any], configuration_identity: str) -> dict[str, Any]:
        if type(step) is not int or step < 0 or step < (self._accepted_restart or {}).get("step", -1):
            raise AnalysisError("Restart step must not move backwards")
        if type(checkpoint) is not bytes or not checkpoint:
            raise AnalysisError("Restart checkpoint must be nonempty bytes")
        self.flush_trajectory()
        self.flush_thermodynamics()
        checkpoint_path = self.root / f"checkpoint-{step:012d}.bin"
        if checkpoint_path.exists():
            raise AnalysisError("Restart checkpoint boundary already exists")
        temporary = checkpoint_path.with_name(f".{checkpoint_path.name}.tmp")
        try:
            temporary.write_bytes(checkpoint)
            with temporary.open("rb") as stream:
                os.fsync(stream.fileno())
            os.replace(temporary, checkpoint_path)
            chunks = {str(path.name): _file_hash(path) for path in sorted(self.root.glob("trajectory-*.jsonl"))}
            thermo_chunks = {str(path.name): _file_hash(path) for path in sorted(self.root.glob("thermodynamics-*.jsonl"))}
            record = {"schema": "island.trajectory-restart.v1", "identities": self.identities, "configuration_identity": configuration_identity, "step": step, "time_ps": time_ps, "checkpoint": checkpoint_path.name, "checkpoint_sha256": _file_hash(checkpoint_path), "trajectory_chunks": chunks, "thermodynamic_chunks": thermo_chunks, "previous_accepted_checkpoint": None if self._accepted_restart is None else self._accepted_restart["checkpoint_sha256"], "openmm": dict(openmm)}
            record["identity"] = _digest(_json_bytes(record))
            _atomic_write(self.root / f"restart-{step:012d}.json", _json_bytes(record))
            self._accepted_restart = record
            manifest = {"schema": SCHEMA, "identities": self.identities, "config": self.config.to_dict(), "accepted_restart": record["identity"]}
            _atomic_write(self.root / "manifest.json", _json_bytes(manifest))
            return dict(record)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise


def _records(root: Path, prefix: str) -> Iterator[dict[str, Any]]:
    for path in sorted(root.glob(f"{prefix}-*.jsonl")):
        for line in path.read_text().splitlines():
            yield json.loads(line)


def read_trajectory_store(root: str | Path) -> tuple[Iterator[dict[str, Any]], Iterator[dict[str, Any]]]:
    """Return lazy frame and thermodynamic iterators after manifest validation."""
    path = Path(root)
    manifest = json.loads((path / "manifest.json").read_text())
    if manifest.get("schema") != SCHEMA:
        raise AnalysisError("Unsupported trajectory store schema")
    return _records(path, "trajectory"), _records(path, "thermodynamics")


def load_restart(root: str | Path, *, expected_identity: str | None = None) -> dict[str, Any]:
    path = Path(root)
    manifest = json.loads((path / "manifest.json").read_text())
    identity = manifest.get("accepted_restart")
    if not identity:
        raise AnalysisError("Trajectory store has no accepted restart")
    records = sorted(path.glob("restart-*.json"))
    record = json.loads(records[-1].read_text())
    body = dict(record)
    stored = body.pop("identity", None)
    if stored != _digest(_json_bytes(body)) or stored != identity or (expected_identity is not None and stored != expected_identity):
        raise AnalysisError("Restart record identity mismatch")
    checkpoint = path / record["checkpoint"]
    if not checkpoint.exists() or _file_hash(checkpoint) != record["checkpoint_sha256"]:
        raise AnalysisError("Restart checkpoint hash mismatch")
    for name, checksum in record["trajectory_chunks"].items():
        chunk = path / name
        if not chunk.exists() or _file_hash(chunk) != checksum:
            raise AnalysisError("Trajectory chunk hash mismatch")
    for name, checksum in record["thermodynamic_chunks"].items():
        chunk = path / name
        if not chunk.exists() or _file_hash(chunk) != checksum:
            raise AnalysisError("Thermodynamic chunk hash mismatch")
    return record


def relocate_trajectory_store(source: str | Path, destination: str | Path) -> Path:
    source_path, destination_path = Path(source), Path(destination)
    if destination_path.exists():
        raise AnalysisError(f"Relocation target already exists: {destination_path}")
    shutil.copytree(source_path, destination_path)
    load_restart(destination_path) if json.loads((destination_path / "manifest.json").read_text()).get("accepted_restart") else None
    return destination_path
