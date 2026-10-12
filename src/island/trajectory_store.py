"""Transactional bounded-memory trajectory and restart persistence."""

from __future__ import annotations

import json
import os
import shutil
import uuid
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

from island.exceptions import AnalysisError

SCHEMA_V1 = "island.trajectory-store.v1"
SCHEMA = "island.trajectory-store.v2"
RESTART_SCHEMA = "island.trajectory-restart.v2"


def _digest(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _identity(value: Any, label: str) -> str:
    if type(value) is not str or len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
        raise AnalysisError(f"{label} must be a lowercase SHA-256 identity")
    return value


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def _file_hash(path: Path) -> str:
    return _digest(path.read_bytes())


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _atomic_write(path: Path, raw: bytes) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("xb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        _fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _writer_lock(root: Path):
    lock = root / ".writer.lock"
    try:
        descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as error:
        raise AnalysisError("Trajectory store writer is already active") from error
    try:
        os.close(descriptor)
        yield
    finally:
        lock.unlink(missing_ok=True)
        _fsync_directory(root)


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
        return {name: getattr(self, name) for name in ("trajectory_interval", "thermodynamic_interval", "checkpoint_interval")}


def _read_manifest(root: Path) -> dict[str, Any]:
    try:
        manifest = json.loads((root / "manifest.json").read_text())
    except Exception as error:
        raise AnalysisError(f"Cannot read trajectory-store manifest: {error}") from error
    if manifest.get("schema") not in {SCHEMA, SCHEMA_V1}:
        raise AnalysisError("Unsupported trajectory store schema")
    return manifest


def _validate_frame_records(root: Path, names: dict[str, str], identities: dict[str, str]) -> None:
    previous_step = -1
    expected_atom_ids = identities.get("atom_ids")
    if expected_atom_ids is not None:
        expected_atom_ids = tuple(expected_atom_ids)
    for name, checksum in names.items():
        chunk = root / name
        if not chunk.exists() or _file_hash(chunk) != checksum:
            raise AnalysisError(f"Trajectory chunk hash mismatch: {name}")
        for line in chunk.read_text().splitlines():
            try:
                frame = json.loads(line)
            except Exception as error:
                raise AnalysisError(f"Invalid trajectory frame in {name}: {error}") from error
            step = frame.get("step")
            if type(step) is not int or step <= previous_step:
                raise AnalysisError("Trajectory frame steps must increase strictly")
            previous_step = step
            atom_ids = frame.get("atom_ids")
            if atom_ids is not None:
                if tuple(atom_ids) != tuple(dict.fromkeys(atom_ids)):
                    raise AnalysisError("Trajectory frame atom IDs are not unique")
                if expected_atom_ids is not None and tuple(atom_ids) != expected_atom_ids:
                    raise AnalysisError("Trajectory frame stable atom IDs do not match the store identity")
                for coordinate_key in ("positions", "coordinates"):
                    if coordinate_key in frame and len(frame[coordinate_key]) != len(atom_ids):
                        raise AnalysisError("Trajectory atom count does not match stable atom IDs")
            lineage = frame.get("trajectory_lineage")
            if lineage is not None and identities.get("trajectory_identity") is not None and lineage != identities["trajectory_identity"]:
                raise AnalysisError("Trajectory lineage identity mismatch")


def _validate_thermodynamic_records(root: Path, names: dict[str, str]) -> None:
    previous_step = -1
    for name, checksum in names.items():
        chunk = root / name
        if not chunk.exists() or _file_hash(chunk) != checksum:
            raise AnalysisError(f"Thermodynamic chunk hash mismatch: {name}")
        for line in chunk.read_text().splitlines():
            record = json.loads(line)
            step = record.get("step")
            if type(step) is not int or step <= previous_step:
                raise AnalysisError("Thermodynamic record steps must increase strictly")
            previous_step = step


def _record_body(record: dict[str, Any]) -> dict[str, Any]:
    body = dict(record)
    stored = body.pop("identity", None)
    if stored != _digest(_json_bytes(body)):
        raise AnalysisError("Restart record identity mismatch")
    return body


def _validate_v2(root: Path, manifest: dict[str, Any]) -> dict[str, Any] | None:
    identities = manifest.get("identities", {})
    for name in ("workflow_identity", "backend_identity", "graph_identity", "packing_identity"):
        _identity(identities.get(name), name)
    accepted = manifest.get("accepted_restart")
    history = manifest.get("restart_history", [])
    if type(history) is not list or len(set(history)) != len(history):
        raise AnalysisError("Restart history is invalid")
    restart_files = sorted(path.name for path in root.glob("restart-*.json"))
    if set(restart_files) != set(history):
        raise AnalysisError("extra restart record is present outside the manifest history")
    if not accepted:
        if history or list(root.glob("checkpoint-*.bin")):
            raise AnalysisError("Uncommitted restart artifacts exist without an accepted boundary")
        return None
    if not history:
        raise AnalysisError("Accepted restart is absent from restart history")
    previous_checkpoint = None
    previous_step = None
    interval = manifest.get("config", {}).get("checkpoint_interval")
    if type(interval) is not int or interval < 1:
        raise AnalysisError("Checkpoint interval policy is invalid")
    accepted_record = None
    checkpoint_names = set()
    for filename in history:
        record = json.loads((root / filename).read_text())
        body = _record_body(record)
        if body.get("schema") != RESTART_SCHEMA:
            raise AnalysisError("Unsupported restart record schema")
        if body.get("identities") != identities:
            raise AnalysisError("Restart identity lineage mismatch")
        if body.get("previous_accepted_checkpoint") != previous_checkpoint:
            raise AnalysisError("Restart previous-boundary link mismatch")
        if type(body.get("step")) is not int or (body["step"] != 0 and body["step"] % interval != 0) or (previous_step is not None and body["step"] - previous_step != interval):
            raise AnalysisError("Committed restart violates checkpoint interval")
        _identity(body.get("configuration_identity"), "configuration_identity")
        _identity(body.get("checkpoint_sha256"), "checkpoint_sha256")
        checkpoint = root / body["checkpoint"]
        checkpoint_names.add(body["checkpoint"])
        if not checkpoint.exists() or _file_hash(checkpoint) != body["checkpoint_sha256"]:
            raise AnalysisError("Restart checkpoint hash mismatch")
        _validate_frame_records(root, body.get("trajectory_chunks", {}), identities)
        _validate_thermodynamic_records(root, body.get("thermodynamic_chunks", {}))
        if record.get("identity") == accepted:
            accepted_record = record
        previous_checkpoint = body["checkpoint_sha256"]
        previous_step = body["step"]
    if accepted_record is None:
        raise AnalysisError("Manifest accepted restart does not match restart history")
    body = _record_body(accepted_record)
    if {path.name for path in root.glob("checkpoint-*.bin")} != checkpoint_names:
        raise AnalysisError("Checkpoint file set does not match restart history")
    if {path.name for path in root.glob("trajectory-*.jsonl")} != set(body.get("trajectory_chunks", {})):
        raise AnalysisError("Trajectory chunk set does not match the accepted restart")
    if {path.name for path in root.glob("thermodynamics-*.jsonl")} != set(body.get("thermodynamic_chunks", {})):
        raise AnalysisError("Thermodynamic chunk set does not match the accepted restart")
    return accepted_record


class TrajectoryStoreWriter:
    """Single-writer, staged trajectory and restart publication."""

    def __init__(self, root: str | Path, *, identities: dict[str, str], config: TrajectoryStoreConfig | None = None) -> None:
        self.root = Path(root)
        if self.root.exists() and any(self.root.iterdir()):
            raise AnalysisError(f"Trajectory store must be a new directory: {self.root}")
        self.root.mkdir(parents=True, exist_ok=True)
        required = ("workflow_identity", "backend_identity", "graph_identity", "packing_identity")
        for name in required:
            _identity(identities.get(name), name)
        self.identities = dict(identities)
        self.config = config or TrajectoryStoreConfig()
        self._frame_buffer: list[dict[str, Any]] = []
        self._thermo_buffer: list[dict[str, Any]] = []
        self._last_step = -1
        self._last_thermo_step = -1
        self._accepted_restart: dict[str, Any] | None = None
        self._history: list[str] = []
        _atomic_write(self.root / "manifest.json", _json_bytes({"schema": SCHEMA, "identities": self.identities, "config": self.config.to_dict(), "accepted_restart": None, "restart_history": []}))

    @property
    def accepted_restart(self) -> dict[str, Any] | None:
        return None if self._accepted_restart is None else dict(self._accepted_restart)

    def append_trajectory_frame(self, frame: dict[str, Any]) -> None:
        if type(frame) is not dict or type(frame.get("step")) is not int or frame["step"] <= self._last_step:
            raise AnalysisError("Trajectory frame steps must increase strictly")
        atom_ids = frame.get("atom_ids")
        if atom_ids is not None and tuple(atom_ids) != tuple(dict.fromkeys(atom_ids)):
            raise AnalysisError("Trajectory frame atom IDs are not unique")
        self._last_step = frame["step"]
        self._frame_buffer.append(dict(frame))
        if len(self._frame_buffer) >= self.config.trajectory_interval:
            self.flush_trajectory()

    def append_thermodynamic_record(self, record: dict[str, Any]) -> None:
        if type(record) is not dict or type(record.get("step")) is not int or record["step"] <= self._last_thermo_step:
            raise AnalysisError("Thermodynamic record steps must increase strictly")
        self._last_thermo_step = record["step"]
        self._thermo_buffer.append(dict(record))
        if len(self._thermo_buffer) >= self.config.thermodynamic_interval:
            self.flush_thermodynamics()

    def flush_trajectory(self) -> str | None:
        if not self._frame_buffer:
            return None
        pending = self.root / ".pending"
        pending.mkdir(exist_ok=True)
        index = len(list(self.root.glob("trajectory-*.jsonl"))) + len(list(pending.glob("trajectory-*.jsonl")))
        path = pending / f"trajectory-{index:08d}.jsonl"
        _atomic_write(path, b"".join(_json_bytes(item) for item in self._frame_buffer))
        self._frame_buffer.clear()
        return _file_hash(path)

    def flush_thermodynamics(self) -> str | None:
        if not self._thermo_buffer:
            return None
        pending = self.root / ".pending"
        pending.mkdir(exist_ok=True)
        index = len(list(self.root.glob("thermodynamics-*.jsonl"))) + len(list(pending.glob("thermodynamics-*.jsonl")))
        path = pending / f"thermodynamics-{index:08d}.jsonl"
        _atomic_write(path, b"".join(_json_bytes(item) for item in self._thermo_buffer))
        self._thermo_buffer.clear()
        return _file_hash(path)

    def publish_restart(self, *, step: int, time_ps: float, checkpoint: bytes, openmm: dict[str, Any], configuration_identity: str) -> dict[str, Any]:
        if type(step) is not int or step < 0:
            raise AnalysisError("Restart step must be a nonnegative integer")
        previous_step = -1 if self._accepted_restart is None else self._accepted_restart["step"]
        if step <= previous_step:
            raise AnalysisError("Restart step must increase strictly")
        if step != 0 and step % self.config.checkpoint_interval != 0:
            raise AnalysisError("Restart step violates checkpoint interval")
        if previous_step >= 0 and step - previous_step != self.config.checkpoint_interval:
            raise AnalysisError("Restart publication skipped a required checkpoint interval")
        if type(checkpoint) is not bytes or not checkpoint:
            raise AnalysisError("Restart checkpoint must be nonempty bytes")
        _identity(configuration_identity, "configuration_identity")
        with _writer_lock(self.root):
            stage = self.root / f".stage-{uuid.uuid4().hex}"
            stage.mkdir()
            moved: list[Path] = []
            old_trajectory = set(self.root.glob("trajectory-*.jsonl"))
            old_thermodynamic = set(self.root.glob("thermodynamics-*.jsonl"))
            try:
                self.flush_trajectory()
                self.flush_thermodynamics()
                pending = self.root / ".pending"
                if pending.exists():
                    for source in sorted(pending.glob("*.jsonl")):
                        shutil.copyfile(source, stage / source.name)
                checkpoint_name = f"checkpoint-{step:012d}.bin"
                restart_name = f"restart-{step:012d}.json"
                checkpoint_path = stage / checkpoint_name
                _atomic_write(checkpoint_path, checkpoint)
                trajectory_chunks = {path.name: _file_hash(path) for path in sorted([*self.root.glob("trajectory-*.jsonl"), *stage.glob("trajectory-*.jsonl")])}
                thermo_chunks = {path.name: _file_hash(path) for path in sorted([*self.root.glob("thermodynamics-*.jsonl"), *stage.glob("thermodynamics-*.jsonl")])}
                record = {"schema": RESTART_SCHEMA, "identities": self.identities, "configuration_identity": configuration_identity, "step": step, "time_ps": float(time_ps), "checkpoint": checkpoint_name, "checkpoint_sha256": _file_hash(checkpoint_path), "trajectory_chunks": trajectory_chunks, "thermodynamic_chunks": thermo_chunks, "previous_accepted_checkpoint": None if self._accepted_restart is None else self._accepted_restart["checkpoint_sha256"], "openmm": dict(openmm)}
                record["identity"] = _digest(_json_bytes(record))
                _atomic_write(stage / restart_name, _json_bytes(record))
                manifest = {"schema": SCHEMA, "identities": self.identities, "config": self.config.to_dict(), "accepted_restart": record["identity"], "restart_history": self._history + [restart_name]}
                for source in sorted(stage.iterdir()):
                    destination = self.root / source.name
                    os.replace(source, destination)
                    moved.append(destination)
                _validate_v2(self.root, manifest)
                _atomic_write(self.root / "manifest.json", _json_bytes(manifest))
                self._accepted_restart = record
                self._history.append(restart_name)
                if pending.exists():
                    shutil.rmtree(pending)
                return dict(record)
            except Exception:
                for path in moved:
                    path.unlink(missing_ok=True)
                for path in set(self.root.glob("trajectory-*.jsonl")) - old_trajectory:
                    path.unlink(missing_ok=True)
                for path in set(self.root.glob("thermodynamics-*.jsonl")) - old_thermodynamic:
                    path.unlink(missing_ok=True)
                pending = self.root / ".pending"
                if pending.exists():
                    shutil.rmtree(pending)
                raise
            finally:
                shutil.rmtree(stage, ignore_errors=True)


def _records(root: Path, prefix: str) -> Iterator[dict[str, Any]]:
    for path in sorted(root.glob(f"{prefix}-*.jsonl")):
        for line in path.read_text().splitlines():
            yield json.loads(line)


def read_trajectory_store(root: str | Path) -> tuple[Iterator[dict[str, Any]], Iterator[dict[str, Any]]]:
    path = Path(root)
    manifest = _read_manifest(path)
    if manifest.get("schema") == SCHEMA:
        _validate_v2(path, manifest)
    return _records(path, "trajectory"), _records(path, "thermodynamics")


def _load_v1(root: Path, manifest: dict[str, Any], expected_identity: str | None) -> dict[str, Any]:
    accepted = manifest.get("accepted_restart")
    records = [json.loads(path.read_text()) for path in sorted(root.glob("restart-*.json"))]
    record = next((item for item in records if item.get("identity") == accepted), None)
    if record is None:
        raise AnalysisError("Manifest accepted restart does not match v1 records")
    if expected_identity is not None and record["identity"] != expected_identity:
        raise AnalysisError("Restart record identity mismatch")
    body = _record_body(record)
    checkpoint = root / body["checkpoint"]
    if not checkpoint.exists() or _file_hash(checkpoint) != body["checkpoint_sha256"]:
        raise AnalysisError("Restart checkpoint hash mismatch")
    _validate_frame_records(root, body.get("trajectory_chunks", {}), manifest.get("identities", {}))
    _validate_thermodynamic_records(root, body.get("thermodynamic_chunks", {}))
    return record


def load_restart(root: str | Path, *, expected_identity: str | None = None) -> dict[str, Any]:
    path = Path(root)
    manifest = _read_manifest(path)
    if manifest.get("schema") == SCHEMA_V1:
        return _load_v1(path, manifest, expected_identity)
    record = _validate_v2(path, manifest)
    if record is None:
        raise AnalysisError("Trajectory store has no accepted restart")
    if expected_identity is not None and record["identity"] != expected_identity:
        raise AnalysisError("Restart record identity mismatch")
    return record


def relocate_trajectory_store(source: str | Path, destination: str | Path) -> Path:
    source_path, destination_path = Path(source), Path(destination)
    if destination_path.exists():
        raise AnalysisError(f"Relocation target already exists: {destination_path}")
    stage = destination_path.parent / f".{destination_path.name}.{uuid.uuid4().hex}.copy"
    try:
        shutil.copytree(source_path, stage)
        manifest = _read_manifest(stage)
        if manifest.get("accepted_restart"):
            load_restart(stage)
        os.replace(stage, destination_path)
        _fsync_directory(destination_path.parent)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return destination_path
