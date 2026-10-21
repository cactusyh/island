"""Reproducible, bounded periodic NPT workflow records and diagnostics.

This module is deliberately independent of a force-field preparation path.  It
accepts a validated periodic backend and records the states produced by the
Reference OpenMM engine.  Density is always computed from the retained box and
the authoritative total mass; target thermodynamic settings are never used as
measurements.  The records are JSON plus opaque checkpoint bytes, so they can be
validated and relocated without importing OpenMM.
"""

from __future__ import annotations

import base64
import json
import math
import os
import shutil
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

import numpy as np

from island.exceptions import AnalysisError, ValidationError

CONFIG_SCHEMA = "island.periodic-npt-config.v1"
MANIFEST_SCHEMA = "island.periodic-npt-manifest.v1"
FRAME_SCHEMA = "island.periodic-npt-frame.v1"
ANALYSIS_SCHEMA = "island.periodic-npt-analysis.v1"
FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
AVOGADRO_DALTON_ANGSTROM3 = 1.66053906660


def _json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def _digest(value: Any) -> str:
    return sha256(_json(value)).hexdigest()


def _signed(value: dict[str, Any]) -> dict[str, Any]:
    body = json.loads(json.dumps(value, sort_keys=True, allow_nan=False))
    body["identity"] = _digest(body)
    return body


def _body(value: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in value.items() if k != "identity"}


def _identity(value: Any, label: str) -> None:
    if type(value) is not str or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValidationError(f"{label} must be a lowercase SHA-256 identity")


def _positive(value: Any, label: str) -> None:
    if type(value) not in (int, float) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
        raise ValidationError(f"{label} must be a positive finite number")


def _integer(value: Any, label: str) -> None:
    if type(value) is not int or value < 1:
        raise ValidationError(f"{label} must be a positive integer")


@dataclass(frozen=True)
class NPTConfiguration:
    """Immutable identity-bound settings for the J33 staged workflow."""

    family: str
    source_identity: str
    final_graph_identity: str
    packing_identity: str
    box_lengths_nm: tuple[float, float, float]
    temperature_kelvin: float = 300.0
    pressure_bar: float = 1.0
    timestep_fs: float = 0.5
    thermostat_type: str = "LangevinMiddleIntegrator"
    barostat_type: str = "MonteCarloBarostat"
    temperature_interval: int = 1
    pressure_interval: int = 1
    seed: int = 2033
    velocity_seed: int = 2034
    barostat_seed: int = 2035
    nvt_steps: int = 8
    npt_steps: int = 32
    production_steps: int = 0
    cutoff_nm: float = 0.9
    electrostatics: str = "ewald"
    platform: str = "Reference"
    box_scaling_policy: str = "orthorhombic_monte_carlo_isotropic"
    minimization_policy: str = "independent_force_checked"
    checkpoint_interval: int = 1
    schema_version: int = 1

    def __post_init__(self) -> None:
        if self.family not in {"PCFF", "OPLS-AA", "GAFF", "GAFF2"}:
            raise ValidationError("Unsupported NPT force-field family")
        for name in ("source_identity", "final_graph_identity", "packing_identity"):
            _identity(getattr(self, name), name)
        if len(self.box_lengths_nm) != 3 or not all(type(x) in (int, float) and not isinstance(x, bool) for x in self.box_lengths_nm):
            raise ValidationError("NPT requires three orthorhombic box lengths")
        for x in self.box_lengths_nm:
            _positive(x, "box length")
        for name in ("temperature_kelvin", "pressure_bar", "timestep_fs", "cutoff_nm"):
            _positive(getattr(self, name), name)
        for name in ("temperature_interval", "pressure_interval", "seed", "velocity_seed", "barostat_seed", "nvt_steps", "npt_steps", "checkpoint_interval"):
            _integer(getattr(self, name), name)
        if type(self.production_steps) is not int or self.production_steps < 0:
            raise ValidationError("production_steps must be a nonnegative integer")
        if self.thermostat_type != "LangevinMiddleIntegrator":
            raise ValidationError("Unsupported thermostat; Reference LangevinMiddleIntegrator required")
        if self.barostat_type != "MonteCarloBarostat":
            raise ValidationError("barostat must be the actual MonteCarloBarostat")
        if self.platform != "Reference":
            raise ValidationError("Unsupported NPT platform; only Reference is reproducible")
        if self.box_scaling_policy != "orthorhombic_monte_carlo_isotropic":
            raise ValidationError("Unsupported box scaling policy")
        if self.electrostatics not in {"ewald", "pme"}:
            raise ValidationError("NPT electrostatics must be explicit Ewald or PME")
        if self.schema_version != 1:
            raise ValidationError("Unsupported NPT configuration schema version")

    @property
    def payload(self) -> dict[str, Any]:
        values = {"schema": CONFIG_SCHEMA, **self.__dict__}
        values["box_lengths_nm"] = list(self.box_lengths_nm)
        return _signed(values)

    @property
    def identity(self) -> str:
        return self.payload["identity"]

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> NPTConfiguration:
        if payload.get("schema") != CONFIG_SCHEMA or payload.get("identity") != _digest(_body(payload)):
            raise ValidationError("NPT configuration identity mismatch")
        values = _body(payload)
        values.pop("schema")
        values["box_lengths_nm"] = tuple(values["box_lengths_nm"])
        return cls(**values)


def _box_matrix(box: Any) -> np.ndarray:
    array = np.asarray(box, dtype=float)
    if array.shape != (3, 3) or not np.isfinite(array).all():
        raise AnalysisError("Every NPT frame requires finite 3x3 box vectors")
    diagonal = np.diag(array)
    if np.any(diagonal <= 0) or not np.allclose(array, np.diag(diagonal), atol=1e-10, rtol=0):
        raise AnalysisError("NPT supports positive orthorhombic boxes only")
    return array


def frame_measurements(frames: Iterable[dict[str, Any]], total_mass_dalton: float) -> list[dict[str, Any]]:
    """Derive volume and density from actual frame boxes.

    The returned records retain target and measured quantities separately.  A
    pressure value is accepted only when the backend explicitly supplies one;
    this avoids turning a target pressure into a claimed measurement.
    """

    _positive(total_mass_dalton, "total mass")
    result = []
    for number, frame in enumerate(frames):
        box = _box_matrix(frame.get("box_vectors_nm"))
        volume_nm3 = float(np.prod(np.diag(box)))
        density = total_mass_dalton * AVOGADRO_DALTON_ANGSTROM3 / (volume_nm3 * 1000.0)
        out = dict(frame)
        out["frame_index"] = number
        out["volume_nm3"] = volume_nm3
        out["density_g_cm3"] = density
        measured = frame.get("measured_pressure_bar")
        if measured is not None:
            if type(measured) not in (int, float) or isinstance(measured, bool) or not math.isfinite(measured):
                raise AnalysisError("Measured pressure must be finite")
            if frame.get("pressure_observation_source") != "backend":
                raise AnalysisError("Measured pressure must identify a backend observation")
        else:
            out["pressure_observation_source"] = "unavailable"
        result.append(out)
    if not result:
        raise AnalysisError("NPT analysis requires at least one retained frame")
    return result


def _blocks(values: np.ndarray, block_count: int) -> dict[str, Any]:
    if values.ndim != 1 or len(values) < block_count * 2:
        return {"status": "insufficient_sampling", "blocks": 0, "reason": "at least two frames per declared block are required"}
    size = len(values) // block_count
    if size < 2:
        return {"status": "insufficient_sampling", "blocks": 0, "reason": "at least two frames per declared block are required"}
    chunks = values[: size * block_count].reshape(block_count, size)
    means = chunks.mean(axis=1)
    return {
        "status": "ok",
        "blocks": block_count,
        "block_size_frames": size,
        "mean": float(values.mean()),
        "standard_deviation": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
        "standard_error": float(means.std(ddof=1) / math.sqrt(block_count)) if block_count > 1 else 0.0,
        "block_means": means.tolist(),
        "first_half_second_half_difference": float(values[len(values) // 2 :].mean() - values[: len(values) // 2].mean()),
        "rolling_mean": np.convolve(values, np.ones(size) / size, mode="valid").tolist(),
    }


def analyze_npt_frames(frames: Iterable[dict[str, Any]], *, total_mass_dalton: float, target_temperature_kelvin: float, target_pressure_bar: float, block_count: int = 4) -> dict[str, Any]:
    """Produce stage-aware finite-run density diagnostics.

    ``equilibrium_status`` is always ``not_established``.  A successful block
    calculation means only that the declared finite-run software criteria were
    computable.
    """

    _positive(target_temperature_kelvin, "target temperature")
    _positive(target_pressure_bar, "target pressure")
    _integer(block_count, "block_count")
    measured = frame_measurements(frames, total_mass_dalton)
    npt = [f for f in measured if f.get("stage") == "npt"]
    if not npt:
        return {"schema": ANALYSIS_SCHEMA, "status": "not_executed", "reason": "No retained NPT frames", "equilibrium_status": "not_established", **FLAGS}
    density = np.asarray([f["density_g_cm3"] for f in npt], dtype=float)
    volume = np.asarray([f["volume_nm3"] for f in npt], dtype=float)
    temperature = np.asarray([f["measured_temperature_kelvin"] for f in npt if f.get("measured_temperature_kelvin") is not None], dtype=float)
    fixed_box = bool(np.allclose(volume, volume[0], atol=1e-12, rtol=0))
    density_blocks = _blocks(density, block_count)
    status = "ok" if density_blocks["status"] == "ok" and not fixed_box else ("blocked" if fixed_box else "insufficient_sampling")
    return {
        "schema": ANALYSIS_SCHEMA,
        "status": status,
        "equilibrium_status": "not_established",
        "npt_executed": True,
        "target_temperature_kelvin": target_temperature_kelvin,
        "target_pressure_bar": target_pressure_bar,
        "frame_count": len(npt),
        "density_units": "g/cm^3",
        "density": density_blocks,
        "volume": _blocks(volume, block_count),
        "temperature": _blocks(temperature, block_count) if len(temperature) else {"status": "blocked", "reason": "backend did not retain kinetic temperature"},
        "initial_density_g_cm3": float(density[0]),
        "final_density_g_cm3": float(density[-1]),
        "fixed_box": fixed_box,
        "density_relaxation_status": "not_demonstrated" if fixed_box else status,
        "pressure_measurement_status": "backend" if any(f.get("measured_pressure_bar") is not None for f in npt) else "unavailable",
        "minimum_distance_violation_count": int(sum((f.get("minimum_distance_angstrom") is not None and f["minimum_distance_angstrom"] < 0.5) for f in npt)),
        "observables": {
            "chain_radius_of_gyration": {"status": "blocked", "reason": "J33 record requires an explicit chain selection; no selection was supplied"},
            "end_to_end_distance": {"status": "blocked", "reason": "Chain ancestry selection is not implicit in a packed backend"},
            "orientation_and_nematic_order": {"status": "blocked", "reason": "Orientation selection is not implicit in a packed backend"},
            "rdf": {"status": "blocked", "reason": "Atom-class selection is required"},
            "contact_fraction": {"status": "blocked", "reason": "Contact cutoff and atom selection are required"},
            "network_connectivity": {"status": "blocked", "reason": "Final graph/network selection is required"},
            "crosslink_density": {"status": "blocked", "reason": "Crosslink provenance and denominator are not supplied"},
        },
        "frame_content_identity": _digest(measured),
        "frames": measured,
        **FLAGS,
    }


def validate_npt_analysis(record: dict[str, Any]) -> dict[str, Any]:
    if record.get("schema") != ANALYSIS_SCHEMA:
        raise AnalysisError("Unsupported NPT analysis schema")
    if record.get("equilibrium_status") != "not_established":
        raise AnalysisError("NPT short run cannot claim thermodynamic equilibrium")
    if record.get("production_validated") is not False or record.get("simulation_readiness") != "not_established":
        raise AnalysisError("NPT readiness flags were promoted")
    frames = frame_measurements(record.get("frames", []), record.get("total_mass_dalton", 1.0)) if record.get("frames") else []
    if frames and record.get("frame_content_identity") != _digest(frames):
        raise AnalysisError("NPT frame content identity mismatch")
    return record


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_bytes(_json(value))


def write_npt_record(root: str | Path, record: dict[str, Any]) -> dict[str, Any]:
    """Publish an analysis record atomically into a new directory."""
    root = Path(root)
    if root.exists():
        raise AnalysisError("NPT record destination exists")
    validate_npt_analysis(record)
    stage = Path(tempfile.mkdtemp(prefix=".npt-record-", dir=root.parent))
    try:
        _write_json(stage / "analysis.json", _signed(record))
        os.replace(stage, root)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return record


def load_npt_record(root: str | Path) -> dict[str, Any]:
    root = Path(root)
    try:
        record = json.loads((root / "analysis.json").read_text())
    except (OSError, ValueError) as error:
        raise AnalysisError(f"Cannot read NPT record: {error}") from error
    if record.get("identity") != _digest(_body(record)):
        raise AnalysisError("NPT record identity mismatch")
    return validate_npt_analysis(_body(record))


def relocate_npt_record(source: str | Path, destination: str | Path) -> dict[str, Any]:
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise AnalysisError("NPT relocation destination exists")
    record = load_npt_record(source)
    shutil.copytree(source, destination)
    try:
        loaded = load_npt_record(destination)
    except BaseException:
        shutil.rmtree(destination, ignore_errors=True)
        raise
    if loaded != record:
        shutil.rmtree(destination, ignore_errors=True)
        raise AnalysisError("Relocated NPT record differs")
    return loaded


def run_npt_workflow(backend: Any, config: NPTConfiguration, destination: str | Path, *, total_mass_dalton: float | None = None) -> dict[str, Any]:
    """Execute minimize → NVT → NPT → optional NVT on one packed backend.

    The stage transition creates a fresh OpenMM context from the exact state of
    the preceding stage and carries positions, velocities and the actual box;
    no coordinates are copied from another trajectory.  The persisted bundle is
    intentionally compact but retains every frame's checkpoint, box and
    thermodynamic observation.  OpenMM is imported only when this function is
    called; offline loading remains data-only.
    """
    from island.periodic_md import _engine, _minimize, _model

    if destination is None:
        raise ValidationError("NPT destination is required")
    root = Path(destination)
    if root.exists():
        raise ValidationError("NPT workflow destination exists")
    backend.validate_integrity()
    payload = backend.payload
    _identity(payload["final_graph_identity"], "final graph identity")
    if config.family != payload["config"]["family"]:
        raise ValidationError("NPT force-field family differs from packed backend")
    if config.final_graph_identity != payload["final_graph_identity"]:
        raise ValidationError("NPT final graph identity differs from packed backend")
    if config.source_identity != payload["config"]["source"].get("sha256"):
        raise ValidationError("NPT source identity differs from packed backend")
    packed_identity = payload["system"]["metadata"].get("packing", {}).get("plan_identity")
    if config.packing_identity != packed_identity:
        raise ValidationError("NPT packing identity differs from packed backend")
    if payload["system"]["box"]["periodic"] != [True, True, True]:
        raise ValidationError("NPT requires fully periodic boundaries")
    lengths_angstrom = payload["config"]["box_lengths"]
    expected_nm = tuple(float(x) / 10.0 for x in lengths_angstrom)
    if not np.allclose(expected_nm, config.box_lengths_nm, atol=1e-12, rtol=0):
        raise ValidationError("NPT configuration box differs from packed backend")
    if total_mass_dalton is None:
        total_mass_dalton = float(sum(site["mass"] for site in payload["final_graph"]["sites"]))
    _positive(total_mass_dalton, "total mass")
    try:
        import openmm
        from openmm import unit
    except ImportError as error:
        raise ValidationError("NPT requires optional dependency OpenMM") from error

    class _StageConfig:
        mode = "nvt"
        pressure_bar = None
        temperature_kelvin = config.temperature_kelvin
        timestep_fs = config.timestep_fs
        friction_per_ps = 1.0
        seed = config.seed
        velocity_seed = config.velocity_seed
        barostat_seed = config.barostat_seed
        barostat_interval = config.pressure_interval
        minimize = True
        minimization_tolerance_kj_mol_nm = 100.0
        minimization_max_iterations = 10000
        schedule = ()

        def temperature(self, _step):
            return self.temperature_kelvin

    nvt_cfg = _StageConfig()
    npt_cfg = _StageConfig()
    npt_cfg.mode = "npt"
    npt_cfg.pressure_bar = config.pressure_bar

    def make_model(stage_cfg):
        return _model(payload, stage_cfg, openmm)

    def make_context(model, stage_cfg):
        integrator = openmm.LangevinMiddleIntegrator(
            stage_cfg.temperature_kelvin, stage_cfg.friction_per_ps, stage_cfg.timestep_fs / 1000
        )
        integrator.setRandomNumberSeed(stage_cfg.seed)
        return openmm.Context(model, integrator, openmm.Platform.getPlatformByName("Reference")), integrator

    atom_ids = sorted(site["id"] for site in payload["final_graph"]["sites"])
    coordinates = [payload["system"]["coordinates"][str(i)] for i in atom_ids]
    frames: list[dict[str, Any]] = []
    model = make_model(nvt_cfg)
    ctx, integrator = make_context(model, nvt_cfg)
    try:
        ctx.setPositions(np.asarray(coordinates, dtype=float) * unit.angstrom)
        box = [length / 10 for length in lengths_angstrom]
        ctx.setPeriodicBoxVectors(
            openmm.Vec3(box[0], 0, 0) * unit.nanometer,
            openmm.Vec3(0, box[1], 0) * unit.nanometer,
            openmm.Vec3(0, 0, box[2]) * unit.nanometer,
        )
        minimization = _minimize(ctx, model, nvt_cfg, openmm, unit)
        ctx.setVelocitiesToTemperature(config.temperature_kelvin, config.velocity_seed)

        def capture(stage: str, step: int, context, checkpoint: bytes | None = None) -> None:
            state = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            if checkpoint is None:
                checkpoint = context.createCheckpoint()
            positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer).tolist()
            vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer).tolist()
            kinetic = float(state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
            measured_temperature = 2 * kinetic / (3 * len(atom_ids) * 0.00831446261815324)
            frames.append({
                "schema": FRAME_SCHEMA,
                "stage": stage,
                "step": step,
                "time_ps": float(state.getTime().value_in_unit(unit.picosecond)),
                "positions_nm": positions,
                "velocities_nm_ps": state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond).tolist(),
                "box_vectors_nm": vectors,
                "volume_nm3": float(np.linalg.det(np.asarray(vectors))),
                "potential_energy_kj_mol": float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)),
                "kinetic_energy_kj_mol": kinetic,
                "target_temperature_kelvin": config.temperature_kelvin,
                "target_pressure_bar": config.pressure_bar if stage == "npt" else None,
                "measured_temperature_kelvin": measured_temperature,
                "measured_temperature_source": "kinetic_energy",
                "measured_pressure_bar": None,
                "pressure_observation_source": "unavailable",
                "minimum_distance_angstrom": None,
                "checkpoint_sha256": sha256(checkpoint).hexdigest(),
                "checkpoint_base64": base64.b64encode(checkpoint).decode(),
            })

        capture("minimize", 0, ctx)
        for i in range(1, config.nvt_steps + 1):
            integrator.step(1)
            if i % config.checkpoint_interval == 0 or i == config.nvt_steps:
                capture("nvt", i, ctx)
        state = ctx.getState(getPositions=True, getVelocities=True)
        npt_model = make_model(npt_cfg)
        npt_ctx, npt_integrator = make_context(npt_model, npt_cfg)
        try:
            npt_ctx.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
            npt_ctx.setPositions(state.getPositions())
            npt_ctx.setVelocities(state.getVelocities())
            for i in range(1, config.npt_steps + 1):
                npt_integrator.step(1)
                if i % config.checkpoint_interval == 0 or i == config.npt_steps:
                    capture("npt", config.nvt_steps + i, npt_ctx)
            if config.production_steps:
                final_state = npt_ctx.getState(getPositions=True, getVelocities=True)
                production_model = make_model(nvt_cfg)
                production_ctx, production_integrator = make_context(production_model, nvt_cfg)
                production_ctx.setPeriodicBoxVectors(*final_state.getPeriodicBoxVectors())
                production_ctx.setPositions(final_state.getPositions())
                production_ctx.setVelocities(final_state.getVelocities())
                for i in range(1, config.production_steps + 1):
                    production_integrator.step(1)
                    if i % config.checkpoint_interval == 0 or i == config.production_steps:
                        capture("production_nvt", config.nvt_steps + config.npt_steps + i, production_ctx)
                del production_ctx, production_integrator
        finally:
            del npt_ctx, npt_integrator
    finally:
        del ctx, integrator

    manifest = _signed({
        "schema": MANIFEST_SCHEMA,
        "config": config.payload,
        "backend_identity": payload["identity"],
        "final_graph_identity": payload["final_graph_identity"],
        "packing_identity": payload["system"]["metadata"].get("packing", {}).get("plan_identity"),
        "source_identity": payload["config"]["source"].get("sha256"),
        "atom_ids": atom_ids,
        "total_mass_dalton": total_mass_dalton,
        "engine": _engine(openmm),
        "minimization": minimization,
        "stage_counts": {stage: sum(f["stage"] == stage for f in frames) for stage in {f["stage"] for f in frames}},
        **FLAGS,
    })
    record = analyze_npt_frames(frames, total_mass_dalton=total_mass_dalton, target_temperature_kelvin=config.temperature_kelvin, target_pressure_bar=config.pressure_bar, block_count=4)
    record.update({"total_mass_dalton": total_mass_dalton, "bindings": {"backend_identity": payload["identity"], "final_graph_identity": payload["final_graph_identity"], "packing_identity": manifest["packing_identity"], "source_identity": manifest["source_identity"], "npt_config_identity": config.identity, "engine": manifest["engine"]}, "manifest_identity": manifest["identity"]})
    validate_npt_analysis(record)
    stage_dir = Path(tempfile.mkdtemp(prefix=".npt-workflow-", dir=root.parent))
    try:
        _write_json(stage_dir / "manifest.json", manifest)
        _write_json(stage_dir / "frames.json", _signed({"schema": "island.periodic-npt-frames.v1", "frames": frames}))
        _write_json(stage_dir / "analysis.json", _signed(record))
        os.replace(stage_dir, root)
    except BaseException:
        shutil.rmtree(stage_dir, ignore_errors=True)
        raise
    return {"manifest": manifest, "frames": frames, "analysis": record}


def load_npt_workflow(root: str | Path) -> dict[str, Any]:
    """Load a relocated workflow while checking manifest, frame and analysis identities."""
    root = Path(root)
    try:
        manifest = json.loads((root / "manifest.json").read_text())
        frames_record = json.loads((root / "frames.json").read_text())
        analysis_record = json.loads((root / "analysis.json").read_text())
    except (OSError, ValueError) as error:
        raise AnalysisError(f"Cannot read NPT workflow: {error}") from error
    for value, label in ((manifest, "manifest"), (frames_record, "frames"), (analysis_record, "analysis")):
        if value.get("identity") != _digest(_body(value)):
            raise AnalysisError(f"NPT {label} identity mismatch")
    if manifest.get("schema") != MANIFEST_SCHEMA or frames_record.get("schema") != "island.periodic-npt-frames.v1":
        raise AnalysisError("Unsupported NPT workflow schema")
    config = NPTConfiguration.from_payload(manifest["config"])
    if manifest.get("source_identity") != config.source_identity or manifest.get("final_graph_identity") != config.final_graph_identity or manifest.get("packing_identity") != config.packing_identity:
        raise AnalysisError("NPT manifest/config identity binding mismatch")
    frames = frames_record["frames"]
    if not frames:
        raise AnalysisError("NPT workflow has no retained frames")
    allowed = {"minimize", "nvt", "npt", "production_nvt"}
    previous_step = -1
    initial_box = [[config.box_lengths_nm[0], 0, 0], [0, config.box_lengths_nm[1], 0], [0, 0, config.box_lengths_nm[2]]]
    for frame in frames:
        if frame.get("schema") != FRAME_SCHEMA or frame.get("stage") not in allowed:
            raise AnalysisError("Unsupported or missing NPT stage label")
        if type(frame.get("step")) is not int or frame["step"] <= previous_step:
            raise AnalysisError("NPT frame steps must increase strictly")
        previous_step = frame["step"]
        raw = base64.b64decode(frame.get("checkpoint_base64", ""), validate=True)
        if sha256(raw).hexdigest() != frame.get("checkpoint_sha256"):
            raise AnalysisError("NPT checkpoint bytes changed")
        box = _box_matrix(frame.get("box_vectors_nm"))
        if frame["stage"] in {"minimize", "nvt"} and not np.allclose(box, initial_box, atol=1e-10, rtol=0):
            raise AnalysisError("NVT stage box changed before NPT")
        if frame.get("measured_temperature_source") not in {None, "kinetic_energy"}:
            raise AnalysisError("Measured temperature source is not a retained kinetic observation")
        if frame.get("stage") != "npt" and frame.get("target_pressure_bar") is not None:
            raise AnalysisError("Pressure target leaked into a non-NPT stage")
    observed_counts = {stage: sum(frame["stage"] == stage for frame in frames) for stage in {frame["stage"] for frame in frames}}
    if manifest.get("stage_counts") != observed_counts:
        raise AnalysisError("NPT stage counts do not match retained frames")
    if analysis_record.get("manifest_identity") != manifest["identity"]:
        raise AnalysisError("NPT analysis/manifest binding mismatch")
    bindings = analysis_record.get("bindings", {})
    expected_bindings = {
        "backend_identity": manifest.get("backend_identity"),
        "final_graph_identity": manifest.get("final_graph_identity"),
        "packing_identity": manifest.get("packing_identity"),
        "source_identity": manifest.get("source_identity"),
        "npt_config_identity": config.identity,
    }
    if any(bindings.get(key) != value for key, value in expected_bindings.items()):
        raise AnalysisError("NPT analysis binding was changed")
    if analysis_record.get("target_temperature_kelvin") != config.temperature_kelvin or analysis_record.get("target_pressure_bar") != config.pressure_bar:
        raise AnalysisError("NPT target thermodynamic metadata does not match configuration")
    if analysis_record.get("frame_content_identity") != _digest(frame_measurements(frames, manifest["total_mass_dalton"])):
        raise AnalysisError("NPT frame content changed")
    validate_npt_analysis(analysis_record)
    return {"manifest": manifest, "frames": frames, "analysis": _body(analysis_record)}


def compare_npt_workflows(first: str | Path, second: str | Path) -> dict[str, Any]:
    a, b = load_npt_workflow(first), load_npt_workflow(second)
    if a["manifest"] != b["manifest"]:
        raise AnalysisError("Cannot compare different NPT contracts")
    return {"same_frames": a["frames"] == b["frames"], "same_analysis": a["analysis"] == b["analysis"], "passed": a["frames"] == b["frames"] and a["analysis"] == b["analysis"]}


__all__ = [
    "ANALYSIS_SCHEMA",
    "CONFIG_SCHEMA",
    "FLAGS",
    "NPTConfiguration",
    "NPTWorkflowConfig",
    "analyze_npt_frames",
    "compare_npt_workflows",
    "frame_measurements",
    "load_npt_record",
    "load_npt_workflow",
    "relocate_npt_record",
    "run_npt_workflow",
    "validate_npt_analysis",
    "write_npt_record",
]

# Descriptive alias used by workflow callers; the frozen dataclass remains the
# single implementation and therefore has one identity representation.
NPTWorkflowConfig = NPTConfiguration
