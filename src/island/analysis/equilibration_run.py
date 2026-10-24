"""Additive J40 staged heating executor using the J33 checkpoint architecture.

The integration body intentionally preserves the bounded J33 Reference model.
Historical executor, record schemas and numerical force-field settings are unchanged.
Heating steps are NVT in the underlying bundle; the signed schedule distinguishes
heating from the subsequent constant-temperature NVT pre-equilibration.
"""

import base64
import os
import shutil
import tempfile
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

import numpy as np

from island.analysis.npt import (
    FLAGS,
    FRAME_SCHEMA,
    MANIFEST_SCHEMA,
    NPTConfiguration,
    _identity,
    _positive,
    _signed,
    _write_json,
    analyze_npt_frames,
    validate_npt_analysis,
)
from island.exceptions import ValidationError


@dataclass(frozen=True)
class HeatingSchedule:
    stages: tuple[tuple[int, float], ...]
    nvt_steps: int
    final_temperature: float

    def __post_init__(self):
        if (
            type(self.stages) is not tuple
            or not self.stages
            or any(type(x) is not tuple or len(x) != 2 for x in self.stages)
            or type(self.nvt_steps) is not int
        ):
            raise ValidationError(
                "Heating schedule requires stages and integer NVT steps"
            )
        previous = 0
        temperature = 0
        for end, target in self.stages:
            _positive(target, "heating temperature")
            if (
                type(end) is not int
                or not previous < end < self.nvt_steps
                or target < temperature
            ):
                raise ValidationError(
                    "Heating boundaries must increase before NVT pre-equilibration"
                )
            previous, temperature = end, target
        if temperature != self.final_temperature:
            raise ValidationError("Heating must end at the NPT target temperature")

    def temperature(self, step):
        for end, target in self.stages:
            if step <= end:
                return target
        return self.final_temperature

    @property
    def payload(self):
        return _signed(
            {
                "schema": "island.heating-schedule.v1",
                "stages": self.stages,
                "nvt_steps": self.nvt_steps,
                "final_temperature": self.final_temperature,
            }
        )


def run_heated_workflow(
    backend: Any,
    config: NPTConfiguration,
    destination: str | Path,
    *,
    total_mass_dalton: float | None = None,
    stop_npt_step: int | None = None,
    continuation: bool = False,
    heating=((64, 100.0), (128, 200.0), (192, 300.0)),
) -> dict[str, Any]:
    """Execute minimize → NVT → NPT → optional NVT on one packed backend.

    The stage transition creates a fresh OpenMM context from the exact state of
    the preceding stage and carries positions, velocities and the actual box;
    no coordinates are copied from another trajectory.  The persisted bundle is
    intentionally compact but retains every frame's checkpoint, box and
    thermodynamic observation.  OpenMM is imported only when this function is
    called; offline loading remains data-only.
    """
    from island.periodic_md import _engine, _minimize, _model

    schedule = HeatingSchedule(
        tuple(tuple(x) for x in heating), config.nvt_steps, config.temperature_kelvin
    )
    capture_temperature = schedule.stages[0][1]

    if stop_npt_step is not None:
        if (
            type(stop_npt_step) is not int
            or not 0 < stop_npt_step < config.npt_steps
            or stop_npt_step % config.checkpoint_interval
        ):
            raise ValidationError(
                "Split boundary must be an interior retained NPT checkpoint"
            )
        if config.production_steps:
            raise ValidationError("Bounded continuation requires production_steps=0")
        continuation = True

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
    packed_identity = (
        payload["system"]["metadata"].get("packing", {}).get("plan_identity")
    )
    if config.packing_identity != packed_identity:
        raise ValidationError("NPT packing identity differs from packed backend")
    if payload["system"]["box"]["periodic"] != [True, True, True]:
        raise ValidationError("NPT requires fully periodic boundaries")
    lengths_angstrom = payload["config"]["box_lengths"]
    expected_nm = tuple(float(x) / 10.0 for x in lengths_angstrom)
    if not np.allclose(expected_nm, config.box_lengths_nm, atol=1e-12, rtol=0):
        raise ValidationError("NPT configuration box differs from packed backend")
    if total_mass_dalton is None:
        total_mass_dalton = float(
            sum(site["mass"] for site in payload["final_graph"]["sites"])
        )
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
    nvt_cfg.temperature_kelvin = schedule.stages[0][1]
    npt_cfg = _StageConfig()
    npt_cfg.temperature_kelvin = config.temperature_kelvin
    npt_cfg.mode = "npt"
    npt_cfg.pressure_bar = config.pressure_bar

    def make_model(stage_cfg):
        return _model(payload, stage_cfg, openmm)

    def make_context(model, stage_cfg):
        integrator = openmm.LangevinMiddleIntegrator(
            stage_cfg.temperature_kelvin,
            stage_cfg.friction_per_ps,
            stage_cfg.timestep_fs / 1000,
        )
        integrator.setRandomNumberSeed(stage_cfg.seed)
        return openmm.Context(
            model, integrator, openmm.Platform.getPlatformByName("Reference")
        ), integrator

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
        ctx.setVelocitiesToTemperature(capture_temperature, config.velocity_seed)

        def capture(
            stage: str, step: int, context, checkpoint: bytes | None = None
        ) -> None:
            state = context.getState(
                getPositions=True, getVelocities=True, getEnergy=True
            )
            if checkpoint is None:
                checkpoint = context.createCheckpoint()
            positions = (
                state.getPositions(asNumpy=True).value_in_unit(unit.nanometer).tolist()
            )
            vectors = (
                state.getPeriodicBoxVectors(asNumpy=True)
                .value_in_unit(unit.nanometer)
                .tolist()
            )
            kinetic = float(
                state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            )
            measured_temperature = (
                2 * kinetic / (3 * len(atom_ids) * 0.00831446261815324)
            )
            frames.append(
                {
                    "schema": FRAME_SCHEMA,
                    "stage": stage,
                    "step": step,
                    "time_ps": float(state.getTime().value_in_unit(unit.picosecond)),
                    "positions_nm": positions,
                    "velocities_nm_ps": state.getVelocities(asNumpy=True)
                    .value_in_unit(unit.nanometer / unit.picosecond)
                    .tolist(),
                    "box_vectors_nm": vectors,
                    "volume_nm3": float(np.linalg.det(np.asarray(vectors))),
                    "potential_energy_kj_mol": float(
                        state.getPotentialEnergy().value_in_unit(
                            unit.kilojoule_per_mole
                        )
                    ),
                    "kinetic_energy_kj_mol": kinetic,
                    "target_temperature_kelvin": capture_temperature
                    if stage != "npt"
                    else config.temperature_kelvin,
                    "target_pressure_bar": config.pressure_bar
                    if stage == "npt"
                    else None,
                    "measured_temperature_kelvin": measured_temperature,
                    "measured_temperature_source": "kinetic_energy",
                    "measured_pressure_bar": None,
                    "pressure_observation_source": "unavailable",
                    "minimum_distance_angstrom": None,
                    "checkpoint_sha256": sha256(checkpoint).hexdigest(),
                    "checkpoint_base64": base64.b64encode(checkpoint).decode(),
                }
            )

        capture("minimize", 0, ctx)
        for i in range(1, config.nvt_steps + 1):
            capture_temperature = schedule.temperature(i)
            integrator.setTemperature(capture_temperature)
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
            initial_npt_checkpoint = (
                npt_ctx.createCheckpoint() if continuation else None
            )
            for i in range(1, (stop_npt_step or config.npt_steps) + 1):
                npt_integrator.step(1)
                if i % config.checkpoint_interval == 0 or i == config.npt_steps:
                    capture("npt", config.nvt_steps + i, npt_ctx)
            if continuation:
                from island.analysis.npt_continuation import restart_assets

                restart = restart_assets(
                    payload,
                    config,
                    npt_model,
                    npt_integrator,
                    openmm,
                    initial_npt_checkpoint,
                )
            if config.production_steps:
                final_state = npt_ctx.getState(getPositions=True, getVelocities=True)
                production_model = make_model(nvt_cfg)
                production_ctx, production_integrator = make_context(
                    production_model, nvt_cfg
                )
                production_ctx.setPeriodicBoxVectors(
                    *final_state.getPeriodicBoxVectors()
                )
                production_ctx.setPositions(final_state.getPositions())
                production_ctx.setVelocities(final_state.getVelocities())
                for i in range(1, config.production_steps + 1):
                    production_integrator.step(1)
                    if (
                        i % config.checkpoint_interval == 0
                        or i == config.production_steps
                    ):
                        capture(
                            "production_nvt",
                            config.nvt_steps + config.npt_steps + i,
                            production_ctx,
                        )
                del production_ctx, production_integrator
        finally:
            del npt_ctx, npt_integrator
    except Exception as error:
        # Preserve actual partial trajectory/checkpoints; never publish it as a
        # complete or resumable success bundle.
        root.mkdir(exist_ok=False)
        _write_json(
            root / "failure.json",
            _signed(
                {
                    "schema": "island.heated-failure.v1",
                    "config": config.payload,
                    "heating": schedule.payload,
                    "backend_identity": payload["identity"],
                    "frames": frames,
                    "exception": type(error).__name__,
                    "diagnostic": str(error),
                    **FLAGS,
                }
            ),
        )
        raise
    finally:
        del ctx, integrator

    manifest = _signed(
        {
            "schema": MANIFEST_SCHEMA,
            "heating": schedule.payload,
            "config": config.payload,
            "backend_identity": payload["identity"],
            "final_graph_identity": payload["final_graph_identity"],
            "packing_identity": payload["system"]["metadata"]
            .get("packing", {})
            .get("plan_identity"),
            "source_identity": payload["config"]["source"].get("sha256"),
            "atom_ids": atom_ids,
            "total_mass_dalton": total_mass_dalton,
            "engine": _engine(openmm),
            "minimization": minimization,
            "stage_counts": {
                stage: sum(f["stage"] == stage for f in frames)
                for stage in {f["stage"] for f in frames}
            },
            **FLAGS,
        }
    )
    record = analyze_npt_frames(
        frames,
        total_mass_dalton=total_mass_dalton,
        target_temperature_kelvin=config.temperature_kelvin,
        target_pressure_bar=config.pressure_bar,
        block_count=4,
    )
    record.update(
        {
            "total_mass_dalton": total_mass_dalton,
            "bindings": {
                "backend_identity": payload["identity"],
                "final_graph_identity": payload["final_graph_identity"],
                "packing_identity": manifest["packing_identity"],
                "source_identity": manifest["source_identity"],
                "npt_config_identity": config.identity,
                "engine": manifest["engine"],
            },
            "manifest_identity": manifest["identity"],
        }
    )
    validate_npt_analysis(record)
    stage_dir = Path(tempfile.mkdtemp(prefix=".npt-workflow-", dir=root.parent))
    try:
        _write_json(stage_dir / "manifest.json", manifest)
        _write_json(
            stage_dir / "frames.json",
            _signed({"schema": "island.periodic-npt-frames.v1", "frames": frames}),
        )
        _write_json(stage_dir / "analysis.json", _signed(record))
        if continuation:
            from island.analysis.npt_continuation import write_restart

            write_restart(stage_dir, restart, manifest, frames, stop_npt_step)
        os.replace(stage_dir, root)
    except BaseException:
        shutil.rmtree(stage_dir, ignore_errors=True)
        raise
    return {"manifest": manifest, "frames": frames, "analysis": record}
