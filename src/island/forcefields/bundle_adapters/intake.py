"""Read explicitly supplied local files only. Never discover, fetch or execute sources."""

from pathlib import Path

from island.forcefields.assignment_core.schema import is_hash, require

from .adapter import adapt
from .schema import CAPTURE, loads, relative_path, sha, timestamp, validate_manifest


def _read(path, expected):
    require(is_hash(expected), "external_source_pin_required")
    path = Path(path)
    require(
        not path.is_symlink() and path.is_file(), "source_file_unavailable_or_symlink"
    )
    raw = path.read_bytes()
    require(sha(raw) == expected, "source_SHA256_mutation:" + path.name)
    return raw


def capture(
    manifest=None,
    manifest_sha256=None,
    audit=None,
    audit_sha256=None,
    *,
    capture_timestamp,
):
    """CLI installs the inherited fence before calling this read-only function."""
    timestamp(capture_timestamp)
    require(
        (manifest is None) == (manifest_sha256 is None),
        "manifest_and_external_hash_required",
    )
    require(
        (audit is None) == (audit_sha256 is None), "audit_and_external_hash_required"
    )
    require(audit is None or manifest is not None, "audit_without_bundle")
    output = {
        "schema": CAPTURE,
        "manifest": None,
        "manifest_raw": None,
        "manifest_sha256": manifest_sha256,
        "files": [],
        "audit": None,
        "audit_raw": None,
        "audit_sha256": audit_sha256,
        "capture_timestamp": capture_timestamp,
    }
    if manifest is None:
        return output
    path = Path(manifest)
    output["manifest_raw"] = _read(path, manifest_sha256).decode("utf-8")
    m = loads(output["manifest_raw"])
    validate_manifest(m)
    root = path.parent.resolve()
    for f in sorted(m["files"], key=lambda f: f["path"]):
        relative = relative_path(f["path"])
        candidate = root
        for part in relative.parts:
            candidate = candidate / part
            require(not candidate.is_symlink(), "source_path_symlink_forbidden")
        require(candidate.resolve().is_relative_to(root), "source_root_escape")
        raw = _read(candidate, f["sha256"])
        require(len(raw) == f["size"], "source_size_mutation")
        output["files"].append(
            {
                "path": f["path"],
                "sha256": sha(raw),
                "size": len(raw),
                "raw_text": raw.decode("utf-8")
                if f["role"] == "record_source"
                else None,
            }
        )
    output["manifest"] = m
    if audit is not None:
        output["audit_raw"] = _read(audit, audit_sha256).decode("utf-8")
        output["audit"] = loads(output["audit_raw"])
    return output


def ingest(*args, **kwargs):
    return adapt(capture(*args, **kwargs))
