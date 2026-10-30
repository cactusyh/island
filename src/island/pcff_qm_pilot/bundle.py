"""Deterministic small numerical artifacts; scratch and wavefunctions stay external."""

import io
import lzma
import tarfile
from pathlib import Path

from island.pcff_completion import require


def pack(directory, output):
    directory = Path(directory)
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w", format=tarfile.USTAR_FORMAT) as archive:
        for path in sorted(directory.rglob("*")):
            require(not path.is_symlink(), "bundle_symlink_forbidden")
            if path.is_file():
                data = path.read_bytes()
                info = tarfile.TarInfo(path.relative_to(directory).as_posix())
                info.size, info.mode, info.mtime = len(data), 0o644, 0
                archive.addfile(info, io.BytesIO(data))
    Path(output).write_bytes(lzma.compress(raw.getvalue(), preset=6))


def unpack(path, destination):
    destination = Path(destination)
    require(not destination.exists(), "bundle_destination_exists")
    with tarfile.open(path, mode="r:xz") as archive:
        names = set()
        for info in archive.getmembers():
            require(info.isfile() and info.name not in names, "unsafe_bundle_member")
            require(
                (destination / info.name)
                .resolve()
                .is_relative_to(destination.resolve()),
                "bundle_path_escape",
            )
            names.add(info.name)
        archive.extractall(destination, filter="data")
