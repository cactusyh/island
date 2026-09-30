"""Strict data-only records, checked relative paths and durable publication."""

import hashlib
import json
import os
import socket
import tempfile
from collections.abc import Mapping
from contextlib import contextmanager
from pathlib import Path

from island.dynamics._checkpoint_data import strict_load
from island.exceptions import WorkflowBusyError, WorkflowError


def checksum(raw):
    return hashlib.sha256(raw).hexdigest()


def json_bytes(value):
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    ).encode()


def read_json(path):
    try:
        return strict_load(Path(path).read_text())
    except Exception as error:
        raise WorkflowError(f"Invalid JSON record {path}: {error}") from error


def child(root, name):
    if (
        not isinstance(name, str)
        or not name
        or Path(name).is_absolute()
        or any(p in (".", "..") for p in name.split("/"))
    ):
        raise WorkflowError("Bundle paths must be checked relative paths")
    path = Path(root) / name
    if not path.resolve().is_relative_to(Path(root).resolve()):
        raise WorkflowError("Bundle path escapes run directory")
    return path


def sync_directory(directory):
    fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def publish(path, raw, *, replace=False):
    """Fsync bytes, atomically publish, then fsync directory (same-host POSIX)."""
    path = Path(path)
    temporary = None
    try:
        fd, temporary = tempfile.mkstemp(prefix=".publish-", dir=path.parent)
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        if replace:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
        sync_directory(path.parent)
    except OSError as error:
        raise WorkflowError(f"Publication failed for {path}: {error}") from error
    finally:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)


@contextmanager
def writer(root):
    path = Path(root) / ".writer.lock"
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as error:
        raise WorkflowBusyError(
            f"Writer lock exists: {path}; inspect owner before recovery"
        ) from error
    except OSError as error:
        raise WorkflowError(f"Cannot acquire workflow writer lock: {error}") from error
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump({"pid": os.getpid(), "host": socket.gethostname()}, stream)
            stream.flush()
            os.fsync(stream.fileno())
        yield
    finally:
        path.unlink()
        sync_directory(root)


# Explicit container tags preserve integer metadata keys and tuples losslessly.
# No dynamic class lookup or executable deserialization is supported.
def encode(value):
    if value is None or type(value) in (str, int, float, bool):
        return value
    if isinstance(value, Mapping):
        if any(type(k) not in (str, int) for k in value):
            raise WorkflowError("Only string/integer metadata keys supported")
        return {"map": [[k, encode(v)] for k, v in value.items()]}
    if isinstance(value, (tuple, list)):
        return {
            "tuple" if isinstance(value, tuple) else "list": [encode(v) for v in value]
        }
    raise WorkflowError(f"Unsupported data-only value: {type(value).__name__}")


def decode(value):
    if value is None or type(value) in (str, int, float, bool):
        return value
    if type(value) is not dict or len(value) != 1:
        raise WorkflowError("Malformed tagged data")
    kind, rows = next(iter(value.items()))
    if type(rows) is not list:
        raise WorkflowError("Malformed container")
    if kind == "map":
        result = {}
        for row in rows:
            if (
                type(row) is not list
                or len(row) != 2
                or type(row[0]) not in (str, int)
                or row[0] in result
            ):
                raise WorkflowError("Invalid/duplicate mapping key")
            result[row[0]] = decode(row[1])
        return result
    if kind in ("tuple", "list"):
        values = [decode(v) for v in rows]
        return tuple(values) if kind == "tuple" else values
    raise WorkflowError("Unknown container tag")
