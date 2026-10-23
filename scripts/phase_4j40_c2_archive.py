"""Lossless, hash-verified archive parts for retained long trajectories."""

import json
import tarfile
import tempfile
from contextlib import contextmanager
from hashlib import sha256
from pathlib import Path

from island.analysis.equilibration import require


@contextmanager
def open_artifacts(directory):
    directory = Path(directory)
    index = directory / "artifact-parts.json"
    if not index.exists():
        with tarfile.open(directory / "artifacts.tar.xz") as archive:
            yield archive
        return
    require(
        not (directory / "artifacts.tar.xz").exists(),
        "Ambiguous archive representation",
    )
    manifest = json.loads(index.read_text())
    require(
        manifest["schema"] == "island.lossless-archive-parts.v1" and manifest["parts"],
        "Invalid archive parts",
    )
    with tempfile.TemporaryDirectory(prefix="j40-c2-parts-") as temporary:
        path = Path(temporary) / "artifacts.tar.xz"
        digest = sha256()
        size = 0
        with path.open("wb") as stream:
            for i, part in enumerate(manifest["parts"]):
                require(
                    part["name"] == f"artifacts.tar.xz.part-{i:03d}",
                    "Changed archive part order/path",
                )
                raw = (directory / part["name"]).read_bytes()
                require(
                    len(raw) == part["size_bytes"]
                    and sha256(raw).hexdigest() == part["sha256"],
                    "Changed archive part",
                )
                stream.write(raw)
                digest.update(raw)
                size += len(raw)
        require(
            digest.hexdigest() == manifest["sha256"] and size == manifest["size_bytes"],
            "Changed combined archive",
        )
        with tarfile.open(path) as archive:
            yield archive
