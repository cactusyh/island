"""Transactional, data-only native preparation bundles. Source libraries stay external."""

import os
import shutil
import tempfile
from copy import deepcopy
from dataclasses import dataclass
from functools import wraps
from math import isfinite
from pathlib import Path

from island.exceptions import ForceFieldError

SCHEMA = "island_prepared_forcefield_bundle_v1"
FILES = {
    "gaff": {
        "system": "system.json",
        "preparation": "preparation.json",
        "prmtop": "result.prmtop",
    },
    "gaff2": {
        "system": "system.json",
        "preparation": "preparation.json",
        "prmtop": "result.prmtop",
    },
    "oplsaa": {"system": "system.json", "parameters": "parameters.json"},
    "pcff": {
        "system": "system.json",
        "assignment": "assignment.json",
        "model": "model.json",
    },
}


class PreparedBundleError(ForceFieldError):
    """Invalid, unavailable or incompatible bundle, artifact, source or publication."""


def _require(ok, message):
    if not ok:
        raise PreparedBundleError(message)


def _boundary(fn):
    @wraps(fn)
    def call(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except PreparedBundleError:
            raise
        except Exception as exc:
            raise PreparedBundleError(f"{fn.__name__}: {exc}") from exc

    return call


@dataclass(frozen=True)
class AmberBundleArtifacts:
    """Explicit permission to copy this original prmtop into the bundle."""

    prmtop: Path

    def __post_init__(self):
        _require(isinstance(self.prmtop, (str, Path)), "Original prmtop path required")
        object.__setattr__(self, "prmtop", Path(self.prmtop))


@dataclass(frozen=True)
class PreparedForceFieldSources:
    """Explicit local library resolution, never embedded in a bundle or searched."""

    opls_xml: Path | None = None
    pcff_frc: Path | None = None

    def __post_init__(self):
        for field in ("opls_xml", "pcff_frc"):
            value = getattr(self, field)
            _require(
                value is None or isinstance(value, (str, Path)), f"Invalid {field} path"
            )
            if value is not None:
                object.__setattr__(self, field, Path(value))


@dataclass(frozen=True)
class LoadedPreparedForceField:
    """Accessors return owned copies of the original system, facade and manifest."""

    _system: object
    _prepared: object
    _manifest: dict

    def __post_init__(self):
        for field in ("_system", "_prepared", "_manifest"):
            object.__setattr__(self, field, deepcopy(getattr(self, field)))

    @property
    def system(self):
        return deepcopy(self._system)

    @property
    def prepared(self):
        return deepcopy(self._prepared)

    @property
    def manifest(self):
        return deepcopy(self._manifest)


def _source(family, sources):
    if sources is None:
        sources = PreparedForceFieldSources()
    _require(
        type(sources) is PreparedForceFieldSources, "PreparedForceFieldSources required"
    )
    if family in ("gaff", "gaff2"):
        _require(
            sources.opls_xml is None and sources.pcff_frc is None,
            "Amber reconstruction uses the original prmtop, not OPLS/PCFF sources",
        )
        return None
    if family == "oplsaa":
        _require(
            sources.opls_xml is not None and sources.pcff_frc is None,
            "OPLS reconstruction requires only explicit opls_xml",
        )
        from .oplsaa import load_oplsaa_source

        return load_oplsaa_source(sources.opls_xml)
    _require(
        sources.pcff_frc is not None and sources.opls_xml is None,
        "PCFF reconstruction requires only explicit pcff_frc",
    )
    from .pcff import load_pcff_source

    return load_pcff_source(sources.pcff_frc)


def _read(root, name):
    from island.workflows import storage

    path = storage.child(root, name)
    _require(
        not path.is_symlink() and path.is_file(),
        f"Missing or unsafe bundle artifact: {name}",
    )
    return path.read_bytes()


def _json(raw):
    from island.dynamics._checkpoint_data import strict_load

    return strict_load(raw.decode("utf-8"))


def _remove(path, original=None):
    """Best-effort cleanup never obscures the original failure."""
    try:
        shutil.rmtree(path)
    except OSError as exc:
        if original is not None and hasattr(original, "add_note"):
            original.add_note(f"Temporary bundle cleanup failed at {path}: {exc}")


def _system_record(encoded, family):
    from island.workflows import storage
    from island.workflows.bundle import system_from

    payload = storage.decode(encoded)
    box = payload.get("box")
    if box is None:
        return system_from(payload)
    # Amber already permits a nonperiodic bounding box. Preserve it as metadata;
    # this does not introduce cell interactions or broaden OPLS/PCFF support.
    _require(
        family in ("gaff", "gaff2")
        and type(box) is dict
        and set(box) == {"lengths", "periodic"},
        "Unsupported system box",
    )
    lengths, periodic = box["lengths"], box["periodic"]
    _require(
        isinstance(lengths, (tuple, list))
        and len(lengths) == 3
        and all(type(v) in (int, float) and isfinite(v) and v > 0 for v in lengths),
        "Invalid nonperiodic box lengths",
    )
    _require(
        isinstance(periodic, (tuple, list))
        and len(periodic) == 3
        and all(type(v) is bool and not v for v in periodic),
        "Periodic boxes are outside bundle scope",
    )
    from island.core import SimulationBox

    system = system_from({**payload, "box": None})
    system.box = SimulationBox(*lengths, periodic=deepcopy(periodic))
    return system


def _reconstruct(root, sources):
    from island.workflows import storage
    from island.workflows.bundle import preparation_from, system_data

    from .preparation import adopt_forcefield

    envelope = _json(_read(root, "manifest.json"))
    _require(
        type(envelope) is dict and set(envelope) == {"payload", "sha256"},
        "Malformed manifest envelope",
    )
    p = envelope["payload"]
    _require(
        envelope["sha256"] == storage.checksum(storage.json_bytes(p)),
        "Manifest checksum mismatch",
    )
    _require(
        type(p) is dict
        and set(p) == {"schema", "family", "files", "prepared", "prepared_identity"},
        "Malformed manifest fields",
    )
    _require(
        p["schema"] == SCHEMA and type(p["family"]) is str and p["family"] in FILES,
        "Unsupported bundle schema/family",
    )
    family = p["family"]
    _require(
        type(p["files"]) is dict and set(p["files"]) == set(FILES[family]),
        "Missing/unexpected logical artifacts",
    )
    raw = {}
    for logical, filename in FILES[family].items():
        entry = p["files"][logical]
        _require(
            type(entry) is dict
            and set(entry) == {"path", "sha256"}
            and entry["path"] == filename,
            f"Unsafe/noncanonical artifact reference: {logical}",
        )
        raw[logical] = _read(root, filename)
        _require(
            storage.checksum(raw[logical]) == entry["sha256"],
            f"Artifact checksum mismatch: {logical}",
        )
    encoded = _json(raw["system"])
    system = _system_record(encoded, family)
    _require(
        storage.json_bytes(storage.encode(system_data(system)))
        == storage.json_bytes(encoded),
        "System record is not lossless/canonical",
    )
    source = _source(family, sources)
    if family in ("gaff", "gaff2"):
        data = storage.decode(_json(raw["preparation"]))
        _require(
            type(data) is dict
            and set(data)
            == {"record", "record_signature", "import_source", "import_provenance"},
            "Malformed Amber reconstruction record",
        )
        # Parse the verified bytes, not a path that could change between hash and use.
        temporary = Path(tempfile.mkdtemp(prefix="island-prepared-read-"))
        original = None
        try:
            prmtop = temporary / "result.prmtop"
            prmtop.write_bytes(raw["prmtop"])
            native = preparation_from(system, data, prmtop)
        except BaseException as exc:
            original = exc
            raise
        finally:
            _remove(temporary, original)
    elif family == "oplsaa":
        from .oplsaa import OPLSParameterizationResult

        native = OPLSParameterizationResult(raw["parameters"].decode("utf-8"))
        native.validate_integrity(system, source)
    else:
        from .pcff import PCFFClass2Result, PCFFModelSpecification

        assignment = PCFFClass2Result(raw["assignment"].decode("utf-8"), source)
        assignment.validate_integrity(system)
        native = PCFFModelSpecification(raw["model"].decode("utf-8"), assignment)
        native.validate_integrity(system)
    prepared = adopt_forcefield(
        system, family, native, source=source if family == "oplsaa" else None
    )
    _require(
        prepared.metadata == p["prepared"]
        and prepared.identity == p["prepared_identity"],
        "Reconstructed facade/native identities contradict manifest",
    )
    return LoadedPreparedForceField(system, prepared, envelope)


@_boundary
def load_prepared_forcefield(path, *, sources=None):
    """Load without preparation/evaluation; source libraries must be explicitly supplied.

    Returns LoadedPreparedForceField. Its .system and .prepared are independent
    owned copies suitable for create_evaluator(). Amber requires ParmEd to reparse
    the bundled original prmtop; OPLS/PCFF need no typing/execution dependencies.
    """
    return _reconstruct(Path(path), sources)


@_boundary
def save_prepared_forcefield(system, prepared, path, *, artifacts=None, sources=None):
    """Validate, stage, reconstruct, then exclusively publish a new bundle directory.

    No overwrite mode. AmberBundleArtifacts explicitly authorizes copying the
    original prmtop. Pinned OPLS/PCFF libraries remain external; supply their local
    paths via PreparedForceFieldSources for candidate validation and later loading.
    """
    from island.workflows import storage
    from island.workflows.bundle import system_data

    from .preparation import PreparedForceField

    target = Path(path)
    _require(
        not target.exists() and not target.is_symlink(),
        f"Destination already exists: {target}",
    )
    _require(target.parent.is_dir(), "Bundle parent directory must already exist")
    _require(
        type(prepared) is PreparedForceField, "Validated PreparedForceField required"
    )
    prepared.validate_integrity(system)
    system_raw = storage.json_bytes(storage.encode(system_data(system)))
    _require(
        system_raw == storage.json_bytes(storage.encode(system_data(prepared._system))),
        "Save requires the original preparation system/coordinates/provenance",
    )
    metadata = prepared.metadata
    family = metadata["family"]
    _source(
        family, sources
    )  # Explicit source pin/dependency validation before staging.
    native = prepared.native_result
    raw = {"system": system_raw}
    if family in ("gaff", "gaff2"):
        _require(
            type(artifacts) is AmberBundleArtifacts,
            "AmberBundleArtifacts(original prmtop) required",
        )
        raw["prmtop"] = artifacts.prmtop.read_bytes()
        _require(
            storage.checksum(raw["prmtop"]) == native.imported_result.source_sha256,
            "Original prmtop checksum differs from signed import",
        )
        raw["preparation"] = storage.json_bytes(
            storage.encode(
                {
                    "record": dict(native.record),
                    "record_signature": native.record_signature,
                    "import_source": native.imported_result.source,
                    "import_provenance": dict(native.imported_result.provenance),
                }
            )
        )
    else:
        _require(artifacts is None, "Amber artifacts are inappropriate for OPLS/PCFF")
        if family == "oplsaa":
            raw["parameters"] = native.json_text.encode("utf-8")
        else:
            raw["assignment"] = native.assignment.json_text.encode("utf-8")
            raw["model"] = native.json_text.encode("utf-8")
    p = {
        "schema": SCHEMA,
        "family": family,
        "prepared": metadata,
        "prepared_identity": prepared.identity,
        "files": {
            key: {"path": FILES[family][key], "sha256": storage.checksum(value)}
            for key, value in raw.items()
        },
    }
    staging = Path(
        tempfile.mkdtemp(prefix=f".{target.name}.candidate-", dir=target.parent)
    )
    reserved = False
    published = False
    original = None
    try:
        for key, value in raw.items():
            storage.publish(staging / FILES[family][key], value)
        storage.publish(
            staging / "manifest.json",
            storage.json_bytes(
                {"payload": p, "sha256": storage.checksum(storage.json_bytes(p))}
            ),
        )
        _reconstruct(
            staging, sources
        )  # Full candidate semantic reconstruction, no backend execution.
        target.mkdir()  # Atomic exclusive reservation; concurrent/existing targets win unchanged.
        reserved = True
        os.replace(
            staging, target
        )  # Replace only our empty reservation, never a prior bundle.
        published = True
        storage.sync_directory(target.parent)
    except BaseException as exc:
        original = exc
        raise
    finally:
        if staging.exists():
            _remove(staging, original)
        if reserved and not published:
            try:
                target.rmdir()  # Do not recursively remove an unexpectedly modified destination.
            except OSError as exc:
                if original is not None and hasattr(original, "add_note"):
                    original.add_note(f"Empty reservation cleanup failed: {exc}")
    return target
