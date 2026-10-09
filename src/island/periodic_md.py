"""Append-only Reference-platform periodic MD over validated J22 backends.

Offline inspection is data-only. Continuation restores owned OpenMM force XML
and binary checkpoints, never a chemistry/typing/parameter preparation path.
"""

import base64
import ctypes
import json
import math
import os
import platform
import shutil
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from functools import wraps
from hashlib import sha256
from pathlib import Path

from island.exceptions import ValidationError
from island.graph import (
    FinalChemicalGraph,
    UnifiedGraphCharges,
    UnifiedParameterAssignment,
    UnifiedTypedGraph,
)
from island.graph.final import pack
from island.periodic import (
    BACKEND_SCHEMA,
    FLAGS,
    PeriodicForceFieldConfig,
    _payload,
    _require,
    _signed,
)

CONFIG = "island_periodic_md_config_v1"
MANIFEST = "island_periodic_md_manifest_v1"
GENERATION = "island_periodic_md_generation_v1"
FRAME = "island_periodic_md_frame_v1"


def _hash(raw):
    return sha256(raw).hexdigest()


def _number(x, label, positive=False):
    _require(
        type(x) in (float, int) and math.isfinite(x),
        f"{label} must be finite and nonboolean",
    )
    if positive:
        _require(x > 0, f"{label} must be positive")


@dataclass(frozen=True)
class PeriodicMDWorkflowConfig:
    family: str
    mode: str = "nvt"
    steps: int = 10
    temperature_kelvin: float = 300.0
    timestep_fs: float = 0.2
    friction_per_ps: float = 1.0
    seed: int = 2026
    velocity_seed: int = 2027
    pressure_bar: float | None = None
    barostat_interval: int = 5
    barostat_seed: int = 2028
    schedule: tuple = ()
    minimize: bool = True
    minimization_tolerance_kj_mol_nm: float = 10.0
    minimization_max_iterations: int = 1000
    platform: str = "Reference"

    def __post_init__(self):
        _require(
            self.family in {"PCFF", "OPLS-AA", "GAFF", "GAFF2"},
            "Unknown workflow family",
        )
        _require(
            self.mode in {"nvt", "npt", "annealing"}, "Unsupported periodic MD mode"
        )
        _require(
            self.platform == "Reference",
            "Unsupported checkpoint platform: only Reference exact restoration is supported",
        )
        for key in (
            "steps",
            "seed",
            "velocity_seed",
            "barostat_seed",
            "barostat_interval",
            "minimization_max_iterations",
        ):
            value = getattr(self, key)
            _require(
                type(value) is int and 0 < value < 2**31,
                f"{key} must be a positive integer below 2^31",
            )
        for key in (
            "temperature_kelvin",
            "timestep_fs",
            "friction_per_ps",
            "minimization_tolerance_kj_mol_nm",
        ):
            _number(getattr(self, key), key, True)
        _require(type(self.minimize) is bool, "minimize must be boolean")
        if self.mode == "npt":
            _number(self.pressure_bar, "pressure_bar", True)
        else:
            _require(
                self.pressure_bar is None, "Pressure is supported only in NPT mode"
            )
        schedule = tuple(tuple(row) for row in self.schedule)
        object.__setattr__(self, "schedule", schedule)
        if self.mode == "annealing":
            _require(
                bool(schedule),
                "Annealing requires explicit (end_step, temperature_kelvin) schedule",
            )
            prior = 0
            for row in schedule:
                _require(
                    len(row) == 2
                    and type(row[0]) is int
                    and prior < row[0] <= self.steps,
                    "Schedule boundaries must strictly increase",
                )
                _number(row[1], "schedule temperature", True)
                prior = row[0]
            _require(
                prior == self.steps and schedule[0][1] == self.temperature_kelvin,
                "Schedule must cover steps and match initial temperature",
            )
        else:
            _require(not schedule, "Schedule requires annealing mode")

    @property
    def payload(self):
        return json.loads(_signed({"schema": CONFIG, **self.__dict__}))["payload"]

    @property
    def identity(self):
        return self.payload["identity"]

    def temperature(self, step):
        if self.mode != "annealing":
            return self.temperature_kelvin
        for end, value in self.schedule:
            if step < end:
                return value
        return self.schedule[-1][1]


def _config(p):
    raw = _payload(pack(p), CONFIG)
    cfg = PeriodicMDWorkflowConfig(
        **{k: v for k, v in raw.items() if k not in {"schema", "identity"}}
    )
    _require(cfg.payload == raw, "Workflow configuration identity mismatch")
    return cfg


def _sync(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _write(path, raw):
    with open(path, "xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())


def _publish_directory(stage, target):
    """Linux atomic, exclusive directory publication (RENAME_NOREPLACE)."""
    _sync(stage)
    libc = ctypes.CDLL(None, use_errno=True)
    _require(
        hasattr(libc, "renameat2"),
        "Unsupported publication platform: renameat2 required",
    )
    rename = libc.renameat2
    rename.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(stage), -100, os.fsencode(target), 1) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), str(target))
    _sync(target.parent)


@contextmanager
def _writer(root):
    path = root / ".writer.lock"
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise ValidationError("Periodic workflow writer already active") from exc
    try:
        os.close(fd)
        yield
    finally:
        path.unlink()
        _sync(root)


def _backend(p):
    _payload(pack(p), BACKEND_SCHEMA)
    graph = FinalChemicalGraph(pack(p["final_graph"]))
    from island.graph.final import _system_from_payload

    graph.validate_integrity(_system_from_payload(p["system"]))
    typed = UnifiedTypedGraph(pack(p["typed_graph"]))
    charges = UnifiedGraphCharges(pack(p["graph_charges"]))
    assignment = UnifiedParameterAssignment(pack(p["parameter_assignment"]))
    typed.validate_integrity(graph)
    charges.validate_integrity(graph, typed)
    assignment.validate_integrity(graph, typed, charges)
    cfg = PeriodicForceFieldConfig.from_json(pack(p["config"]))
    for key, value in {
        "final_graph_identity": graph.identity,
        "typed_graph_identity": typed.identity,
        "charge_identity": charges.identity,
        "assignment_identity": assignment.identity,
        "config_identity": cfg.identity,
    }.items():
        _require(p[key] == value, f"Stale backend {key}")
    _require(
        cfg.payload["family"] == typed.payload["force_field"]
        and cfg.payload["source"] == typed.payload["source"],
        "Backend family/source mismatch",
    )
    _require(
        p["system"]["box"]["lengths"] == cfg.payload["box_lengths"]
        and p["system"]["box"]["periodic"] == [True] * 3,
        "Backend initial box mismatch",
    )
    return graph


def _vectors(value, n, label):
    _require(type(value) is list and len(value) == n, f"Invalid {label} coverage")
    for row in value:
        _require(type(row) is list and len(row) == 3, f"Invalid {label} shape")
        for x in row:
            _number(x, label)


def _validate_frame(p, manifest, prior_step):
    _payload(pack(p), FRAME)
    cfg = _config(manifest["config"])
    _require(
        p["manifest_identity"] == manifest["identity"], "Stale checkpoint manifest"
    )
    _require(
        type(p["step"]) is int and p["step"] == prior_step + 1,
        "Checkpoint step lineage mismatch",
    )
    _require(p["step"] <= cfg.steps, "Checkpoint exceeds configured steps")
    _number(p["time_ps"], "time")
    _require(
        abs(p["time_ps"] - p["step"] * cfg.timestep_fs / 1000) < 1e-10,
        "Checkpoint time mismatch",
    )
    for key in ("positions_nm", "velocities_nm_ps"):
        _vectors(p[key], len(manifest["atom_ids"]), key)
    _vectors(p["box_vectors_nm"], 3, "box vectors")
    for i, row in enumerate(p["box_vectors_nm"]):
        _require(
            row[i] > 0 and all(x == 0 for j, x in enumerate(row) if j != i),
            "Invalid orthorhombic checkpoint box",
        )
        _require(
            row[i] * 10 > 2 * manifest["backend"]["config"]["nonbonded_cutoff"],
            "Dynamic box violates periodic cutoff",
        )
    if cfg.mode != "npt":
        _require(
            p["box_vectors_nm"] == manifest["initial_box_vectors_nm"],
            "Unexpected box change outside NPT",
        )
    _require(
        p["temperature_kelvin"] == cfg.temperature(p["step"])
        and p["pressure_bar"] == cfg.pressure_bar,
        "Checkpoint thermodynamic settings changed",
    )
    _require(
        p["integrator"] == manifest["integrator"]
        and p["barostat"] == manifest["barostat"],
        "Checkpoint integrator/barostat mismatch",
    )
    _require(
        p["backend_identity"] == manifest["backend"]["identity"]
        and p["graph_identity"] == manifest["backend"]["final_graph_identity"],
        "Checkpoint graph/backend mismatch",
    )
    _require(
        p["packing_identity"] == manifest["packing_identity"],
        "Checkpoint packing identity mismatch",
    )
    for key in ("potential_energy_kj_mol", "kinetic_energy_kj_mol"):
        _number(p[key], key)
    raw = base64.b64decode(p["checkpoint_base64"], validate=True)
    _require(_hash(raw) == p["checkpoint_sha256"], "Tampered binary checkpoint")


@dataclass(frozen=True)
class PeriodicMDWorkflow:
    json_text: str

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["manifest"]["identity"]


def _record_boundary(function):
    @wraps(function)
    def checked(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except ValidationError:
            raise
        except (ValueError, TypeError, KeyError, OSError) as error:
            raise ValidationError(
                f"Malformed periodic workflow record: {error}"
            ) from error

    return checked


@_record_boundary
def load_periodic_md_workflow(
    path,
    *,
    expected_identity=None,
    expected_config=None,
    expected_backend_identity=None,
    expected_checkpoint_identity=None,
):
    """Inspect data and all checkpoints without importing any scientific backend."""
    root = Path(path)
    manifest = _payload((root / "manifest.json").read_text(), MANIFEST)
    cfg = _config(manifest["config"])
    graph = _backend(manifest["backend"])
    _require(
        manifest["config"]["family"] == manifest["backend"]["config"]["family"],
        "Wrong-family workflow backend",
    )
    _require(
        manifest["atom_ids"] == sorted(s["id"] for s in graph.payload["sites"]),
        "Workflow atom ordering mismatch",
    )
    _require(
        manifest["initial_box_vectors_nm"]
        == [
            [length / 10 if i == j else 0 for j in range(3)]
            for i, length in enumerate(manifest["backend"]["config"]["box_lengths"])
        ],
        "Workflow initial box mismatch",
    )
    _require(
        manifest["packing_identity"]
        == manifest["backend"]["system"]["metadata"]
        .get("packing", {})
        .get("plan_identity"),
        "Workflow packing identity mismatch",
    )
    _require(
        manifest["integrator"] == _integrator_config(cfg)
        and manifest["barostat"] == _barostat_config(cfg),
        "Workflow integrator/barostat mismatch",
    )
    _require(
        _hash((root / "system.xml").read_bytes()) == manifest["system_xml_sha256"],
        "Tampered OpenMM system XML",
    )
    _require(all(manifest[k] == v for k, v in FLAGS.items()), "Changed readiness flags")
    if expected_identity is not None:
        _require(
            manifest["identity"] == expected_identity,
            "Unexpected workflow manifest identity",
        )
    if expected_config is not None:
        _require(
            cfg.identity == expected_config.identity,
            "Changed workflow configuration: temperature/pressure/timestep/seed/schedule",
        )
    if expected_backend_identity is not None:
        _require(
            manifest["backend"]["identity"] == expected_backend_identity,
            "Changed workflow backend identity",
        )
    generations = sorted(root.glob("generation-*"))
    _require(bool(generations), "No committed checkpoint")
    frames = []
    previous = None
    step = -1
    for n, directory in enumerate(generations):
        _require(
            directory.name == f"generation-{n:08d}" and directory.is_dir(),
            "Broken generation lineage",
        )
        p = _payload((directory / "record.json").read_text(), GENERATION)
        _require(
            p["manifest_identity"] == manifest["identity"]
            and p["previous"] == previous
            and p["generation"] == n,
            "Stale generation lineage",
        )
        if n == 0:
            diagnostic_raw = (root / "minimization.json").read_bytes()
            diagnostic = _payload(
                diagnostic_raw.decode(), "island_periodic_minimization_diagnostic_v1"
            )
            _require(
                _hash(diagnostic_raw) == p["minimization_sha256"]
                and diagnostic["manifest_identity"] == manifest["identity"],
                "Minimization diagnostic identity mismatch",
            )
        else:
            _require(
                p["minimization_sha256"] is None,
                "Unexpected minimization during resume",
            )
        _require(bool(p["frames"]), "Empty checkpoint generation")
        for frame in p["frames"]:
            _validate_frame(frame, manifest, step)
            step = frame["step"]
            frames.append(frame)
        previous = p["identity"]
    if expected_checkpoint_identity is not None:
        _require(
            previous == expected_checkpoint_identity,
            "Stale workflow checkpoint identity",
        )
    return PeriodicMDWorkflow(
        json.dumps(
            {
                "manifest": manifest,
                "frames": frames,
                "generation": len(generations) - 1,
                "checkpoint_identity": previous,
            },
            sort_keys=True,
            allow_nan=False,
        )
    )


def periodic_md_workflow_status(path, **expectations):
    p = load_periodic_md_workflow(path, **expectations).payload
    return {
        "manifest_identity": p["manifest"]["identity"],
        "checkpoint_identity": p["checkpoint_identity"],
        "step": p["frames"][-1]["step"],
        "status": "complete"
        if p["frames"][-1]["step"] == p["manifest"]["config"]["steps"]
        else "ready_to_resume",
        **FLAGS,
    }


def read_periodic_md_frames(path, **expectations):
    return tuple(load_periodic_md_workflow(path, **expectations).payload["frames"])


def _mm():
    try:
        import openmm
        from openmm import unit
    except ImportError as exc:
        raise ValidationError(
            "Periodic MD requires optional dependency OpenMM"
        ) from exc
    return openmm, unit


def _engine(mm):
    return {
        "openmm_version": mm.__version__,
        "platform": "Reference",
        "machine": platform.machine(),
        "system": platform.system(),
        "byteorder": __import__("sys").byteorder,
    }


def _integrator_config(cfg):
    return {
        "type": "LangevinMiddleIntegrator",
        "timestep_fs": cfg.timestep_fs,
        "friction_per_ps": cfg.friction_per_ps,
        "seed": cfg.seed,
        "velocity_seed": cfg.velocity_seed,
    }


def _barostat_config(cfg):
    return (
        None
        if cfg.mode != "npt"
        else {
            "type": "MonteCarloBarostat",
            "pressure_bar": cfg.pressure_bar,
            "interval": cfg.barostat_interval,
            "seed": cfg.barostat_seed,
        }
    )


def _model(backend, cfg, mm):
    from island._periodic_openmm import build

    model = build(backend["numerical_data"], backend["config"], mm)
    if cfg.mode == "npt":
        bar = mm.MonteCarloBarostat(
            cfg.pressure_bar, cfg.temperature_kelvin, cfg.barostat_interval
        )
        bar.setRandomNumberSeed(cfg.barostat_seed)
        model.addForce(bar)
    return model


def _context(model, cfg, mm):
    integrator = mm.LangevinMiddleIntegrator(
        cfg.temperature_kelvin, cfg.friction_per_ps, cfg.timestep_fs / 1000
    )
    integrator.setRandomNumberSeed(cfg.seed)
    return mm.Context(
        model, integrator, mm.Platform.getPlatformByName("Reference")
    ), integrator


def _temperature(ctx, integrator, cfg, step, mm):
    value = cfg.temperature(step)
    integrator.setTemperature(value)
    if cfg.mode == "npt":
        ctx.setParameter(mm.MonteCarloBarostat.Temperature(), value)


def _frame(ctx, manifest, step, mm, unit):
    state = ctx.getState(getPositions=True, getVelocities=True, getEnergy=True)
    raw = ctx.createCheckpoint()
    cfg = _config(manifest["config"])
    p = {
        "schema": FRAME,
        "manifest_identity": manifest["identity"],
        "backend_identity": manifest["backend"]["identity"],
        "graph_identity": manifest["backend"]["final_graph_identity"],
        "packing_identity": manifest["packing_identity"],
        "step": step,
        "time_ps": state.getTime().value_in_unit(unit.picosecond),
        "positions_nm": state.getPositions(asNumpy=True)
        .value_in_unit(unit.nanometer)
        .tolist(),
        "velocities_nm_ps": state.getVelocities(asNumpy=True)
        .value_in_unit(unit.nanometer / unit.picosecond)
        .tolist(),
        "box_vectors_nm": state.getPeriodicBoxVectors(asNumpy=True)
        .value_in_unit(unit.nanometer)
        .tolist(),
        "potential_energy_kj_mol": state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        ),
        "kinetic_energy_kj_mol": state.getKineticEnergy().value_in_unit(
            unit.kilojoule_per_mole
        ),
        "temperature_kelvin": cfg.temperature(step),
        "pressure_bar": cfg.pressure_bar,
        "integrator": manifest["integrator"],
        "barostat": manifest["barostat"],
        "checkpoint_base64": base64.b64encode(raw).decode(),
        "checkpoint_sha256": _hash(raw),
    }
    p = json.loads(_signed(p))["payload"]
    _validate_frame(p, manifest, step - 1)
    _require(ctx.getStepCount() == step, "OpenMM checkpoint step mismatch")
    return p


def _generation(stage, manifest, frames, n, previous, minimization_sha256=None):
    record = {
        "schema": GENERATION,
        "manifest_identity": manifest["identity"],
        "generation": n,
        "minimization_sha256": minimization_sha256,
        "previous": previous,
        "frames": frames,
    }
    _write(stage / "record.json", _signed(record).encode())


def _advance(ctx, integrator, manifest, start, count, mm, unit):
    cfg = _config(manifest["config"])
    frames = []
    for step in range(start + 1, start + count + 1):
        _temperature(ctx, integrator, cfg, step - 1, mm)
        integrator.step(1)
        _temperature(ctx, integrator, cfg, step, mm)
        frames.append(_frame(ctx, manifest, step, mm, unit))
    return frames


def _budget(steps, remaining):
    count = remaining if steps is None else steps
    _require(
        type(count) is int and 0 <= count <= remaining, "Invalid workflow step budget"
    )
    return count


def _minimize(ctx, model, cfg, mm, unit):
    import numpy as np

    initial = (
        ctx.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )
    box = ctx.getState().getPeriodicBoxVectors()
    mm.LocalEnergyMinimizer.minimize(
        ctx,
        cfg.minimization_tolerance_kj_mol_nm / math.sqrt(model.getNumParticles()),
        cfg.minimization_max_iterations,
    )
    state = ctx.getState(getPositions=True, getEnergy=True, getForces=True)
    verify, vi = _context(
        mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(model)), cfg, mm
    )
    try:
        verify.setPeriodicBoxVectors(*box)
        verify.setPositions(state.getPositions())
        fresh = verify.getState(getEnergy=True, getForces=True)
        e = fresh.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces = fresh.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        largest = float(np.max(np.linalg.norm(forces, axis=1)))
        _require(
            np.isfinite(forces).all() and math.isfinite(e),
            "Nonfinite minimized candidate",
        )
        _require(
            e <= initial + 1e-8 and largest <= cfg.minimization_tolerance_kj_mol_nm,
            "Minimization did not meet independent force/energy criterion",
        )
        _require(
            np.allclose(
                forces,
                state.getForces(asNumpy=True).value_in_unit(
                    unit.kilojoule_per_mole / unit.nanometer
                ),
                atol=1e-9,
                rtol=1e-12,
            ),
            "Independent minimization force mismatch",
        )
        _require(
            state.getPeriodicBoxVectors() == box, "Minimization changed periodic box"
        )
        return {
            "initial_energy_kj_mol": initial,
            "final_energy_kj_mol": e,
            "maximum_force_kj_mol_nm": largest,
            "independent_context_verified": True,
        }
    finally:
        del verify, vi


def start_periodic_md_workflow(backend, config, path, *, steps=None):
    """Validate an existing backend, minimize optionally, initialize and publish."""
    from island.periodic import PeriodicParameterizedSystem

    _require(
        type(backend) is PeriodicParameterizedSystem,
        "PeriodicParameterizedSystem required",
    )
    _require(
        type(config) is PeriodicMDWorkflowConfig, "PeriodicMDWorkflowConfig required"
    )
    backend.validate_integrity()
    p = backend.payload
    _require(p["config"]["family"] == config.family, "Wrong-family periodic backend")
    target = Path(path)
    _require(not target.exists(), "Workflow destination exists")
    count = _budget(steps, config.steps)
    mm, unit = _mm()
    model = _model(p, config, mm)
    xml = mm.XmlSerializer.serialize(model)
    manifest = {
        "schema": MANIFEST,
        "backend": p,
        "config": config.payload,
        "engine": _engine(mm),
        "atom_ids": sorted(s["id"] for s in p["final_graph"]["sites"]),
        "initial_box_vectors_nm": [
            [length / 10 if i == j else 0 for j in range(3)]
            for i, length in enumerate(p["config"]["box_lengths"])
        ],
        "packing_identity": p["system"]["metadata"]
        .get("packing", {})
        .get("plan_identity"),
        "system_xml_sha256": _hash(xml.encode()),
        "integrator": _integrator_config(config),
        "barostat": _barostat_config(config),
        **FLAGS,
    }
    manifest = json.loads(_signed(manifest))["payload"]
    stage = Path(tempfile.mkdtemp(prefix=".periodic-md-", dir=target.parent))
    ctx, integrator = _context(model, config, mm)
    try:
        ctx.setPositions(
            [p["system"]["coordinates"][str(i)] for i in manifest["atom_ids"]]
            * unit.angstrom
        )
        minimization = (
            _minimize(ctx, model, config, mm, unit)
            if config.minimize
            else {"status": "not_requested"}
        )
        ctx.setVelocitiesToTemperature(config.temperature_kelvin, config.velocity_seed)
        frames = [_frame(ctx, manifest, 0, mm, unit)]
        frames.extend(_advance(ctx, integrator, manifest, 0, count, mm, unit))
        _write(stage / "manifest.json", pack(manifest).encode())
        _write(stage / "system.xml", xml.encode())
        _write(
            stage / "minimization.json",
            _signed(
                {
                    "schema": "island_periodic_minimization_diagnostic_v1",
                    "manifest_identity": manifest["identity"],
                    **minimization,
                }
            ).encode(),
        )
        gen = stage / "generation-00000000"
        gen.mkdir()
        _generation(
            gen,
            manifest,
            frames,
            0,
            None,
            _hash((stage / "minimization.json").read_bytes()),
        )
        load_periodic_md_workflow(stage, expected_identity=manifest["identity"])
        _publish_directory(stage, target)
    except BaseException as error:
        # Failed candidates are retained independently; no valid run is replaced.
        try:
            state = ctx.getState(getPositions=True)
            diagnostic = {
                "schema": "island_periodic_failed_candidate_v1",
                "error": str(error),
                "manifest_identity": manifest["identity"],
                "positions_nm": state.getPositions(asNumpy=True)
                .value_in_unit(unit.nanometer)
                .tolist(),
            }
            _write(
                stage / "failure.json", json.dumps(diagnostic, allow_nan=False).encode()
            )
            error.add_note(f"Failed candidate retained at {stage}")
        except Exception as diagnostic_error:  # noqa: BLE001
            error.add_note(
                f"Could not retain candidate coordinates: {diagnostic_error}"
            )
        raise
    finally:
        del ctx, integrator
    return periodic_md_workflow_status(target)


def resume_periodic_md_workflow(
    path,
    *,
    steps=None,
    expected_identity=None,
    expected_config=None,
    expected_backend_identity=None,
    expected_checkpoint_identity=None,
):
    root = Path(path)
    with _writer(root):
        loaded = load_periodic_md_workflow(
            root,
            expected_identity=expected_identity,
            expected_config=expected_config,
            expected_backend_identity=expected_backend_identity,
            expected_checkpoint_identity=expected_checkpoint_identity,
        ).payload
        manifest = loaded["manifest"]
        cfg = _config(manifest["config"])
        last = loaded["frames"][-1]
        count = _budget(steps, cfg.steps - last["step"])
        if count == 0:
            return periodic_md_workflow_status(root)
        mm, unit = _mm()
        _require(
            _engine(mm) == manifest["engine"],
            "Unsupported OpenMM version/platform/architecture for exact checkpoint restoration",
        )
        model = _model(manifest["backend"], cfg, mm)
        _require(
            _hash(mm.XmlSerializer.serialize(model).encode())
            == manifest["system_xml_sha256"],
            "Reconstructed model contradicts checkpoint force field",
        )
        ctx, integrator = _context(model, cfg, mm)
        stage = Path(tempfile.mkdtemp(prefix=".candidate-", dir=root))
        try:
            ctx.loadCheckpoint(
                base64.b64decode(last["checkpoint_base64"], validate=True)
            )
            _temperature(ctx, integrator, cfg, last["step"], mm)
            restored = _frame(ctx, manifest, last["step"], mm, unit)
            _require(
                restored == last,
                "Unsupported exact restoration: checkpoint state differs",
            )
            frames = _advance(ctx, integrator, manifest, last["step"], count, mm, unit)
            n = loaded["generation"] + 1
            _generation(stage, manifest, frames, n, loaded["checkpoint_identity"])
            candidate = _payload((stage / "record.json").read_text(), GENERATION)
            _require(
                candidate["previous"] == loaded["checkpoint_identity"]
                and candidate["frames"] == frames
                and candidate["manifest_identity"] == manifest["identity"]
                and candidate["generation"] == n,
                "Candidate checkpoint publication validation failed",
            )
            for f in frames:
                _validate_frame(f, manifest, f["step"] - 1)
            _publish_directory(stage, root / f"generation-{n:08d}")
        except BaseException as error:
            _write(
                stage / "failure.json",
                json.dumps(
                    {"error": str(error), "manifest_identity": manifest["identity"]}
                ).encode(),
            )
            error.add_note(
                f"Previous checkpoint preserved; candidate diagnostic: {stage}"
            )
            raise
        finally:
            del ctx, integrator
    return periodic_md_workflow_status(root)


def save_periodic_md_workflow(path, destination, *, expected_identity=None):
    """Validate and transactionally relocate/copy every published checkpoint."""
    source = Path(path)
    target = Path(destination)
    _require(not target.exists(), "Workflow destination exists")
    with _writer(source):
        current = load_periodic_md_workflow(source, expected_identity=expected_identity)
        stage = Path(tempfile.mkdtemp(prefix=".periodic-copy-", dir=target.parent))
        try:
            for name in ("manifest.json", "system.xml", "minimization.json"):
                _write(stage / name, (source / name).read_bytes())
            for directory in sorted(source.glob("generation-*")):
                out = stage / directory.name
                out.mkdir()
                _write(out / "record.json", (directory / "record.json").read_bytes())
                _sync(out)
            load_periodic_md_workflow(stage, expected_identity=current.identity)
            _publish_directory(stage, target)
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    return periodic_md_workflow_status(target)


def minimize_periodic_energy(backend, config, path):
    """Publish an independently checked minimized step-zero periodic workflow."""
    _require(
        config.minimize is True, "Explicit minimization-enabled configuration required"
    )
    return start_periodic_md_workflow(backend, config, path, steps=0)


def compare_periodic_md_workflows(first, second):
    """Compare complete trajectories exactly, including opaque engine RNG bytes.

    Generation layout may differ; configuration, initial backend, each frame and
    every checkpoint byte must agree. No tolerance-based equivalence is implied.
    """
    a = load_periodic_md_workflow(first).payload
    b = load_periodic_md_workflow(second).payload
    _require(
        a["manifest"] == b["manifest"], "Cannot compare different workflow contracts"
    )
    complete = all(
        p["frames"][-1]["step"] == p["manifest"]["config"]["steps"] for p in (a, b)
    )
    states = a["frames"] == b["frames"]
    checkpoints = [f["checkpoint_base64"] for f in a["frames"]] == [
        f["checkpoint_base64"] for f in b["frames"]
    ]
    return {
        "manifest_identity": a["manifest"]["identity"],
        "complete": complete,
        "exact_frames": states,
        "exact_checkpoint_and_rng_bytes": checkpoints,
        "passed": complete and states and checkpoints,
        **FLAGS,
    }
