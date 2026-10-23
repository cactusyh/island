"""Offline observations and fail-closed analysis of complete or failed real runs."""

import base64
import json
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.analysis.equilibration import EquilibrationPolicy, diagnose, require, series
from island.analysis.equilibration_protocol import CorrectiveProtocol
from island.analysis.npt import (
    FLAGS,
    NPTConfiguration,
    _body,
    _digest,
    _signed,
    frame_measurements,
    load_npt_workflow,
)
from island.periodic_md import _backend


def read(path):
    return json.loads(Path(path).read_text())


def completed_history(root):
    """Strict frame/restart import for the additive C2 statistics schema.

    The immutable v1 summary stays hash-bound, not reinterpreted as a C2
    statistic. V1 NumPy convolution is not bit portable across installations;
    C2 statistics use the unchanged J40 fsum implementation on actual frames.
    No trajectory values or scientific bounds are adjusted.
    """
    workflow = load_npt_workflow(root)
    restart = read(root / "restart.json")
    require(restart["identity"] == _digest(_body(restart)), "Changed restart identity")
    expected_files = {
        "npt-initial-checkpoint.txt",
        "backend.json",
        "system.xml",
        "integrator.xml",
        "manifest.json",
        "frames.json",
        "analysis.json",
    }
    require(set(restart["files"]) == expected_files, "Changed restart file inventory")
    require(
        all(
            sha256((root / name).read_bytes()).hexdigest() == digest
            for name, digest in restart["files"].items()
        ),
        "Changed restart source/frame assets",
    )
    backend = read(root / "backend.json")
    _backend(backend)
    m, frames = workflow["manifest"], workflow["frames"]
    require(
        restart["backend_identity"] == backend["identity"] == m["backend_identity"],
        "Changed complete backend",
    )
    for key in (
        "assignment_identity",
        "typed_graph_identity",
        "charge_identity",
        "final_graph_identity",
    ):
        require(restart[key] == backend[key], "Changed restart " + key)
    require(
        restart["frame_content_identity"] == _digest(frames),
        "Changed complete frame identity",
    )
    require(
        restart["checkpoint_identity"]
        == restart["rng_state_sha256"]
        == frames[-1]["checkpoint_sha256"],
        "Changed final checkpoint/RNG",
    )
    from island.analysis.npt_continuation import parent_identity

    require(
        restart["parent_workflow_identity"] == parent_identity(m)
        and restart["manifest_identity"] == m["identity"]
        and restart["npt_config_identity"] == m["config"]["identity"]
        and restart["source_identity"] == m["source_identity"],
        "Changed complete restart bindings",
    )
    require(
        restart["split_boundary_npt_step"] is None,
        "Incomplete history relabeled complete",
    )
    require(
        frames[-1]["step"] == m["config"]["nvt_steps"] + m["config"]["npt_steps"],
        "NPT budget incomplete",
    )
    return workflow


def trajectory(root, label="uninterrupted"):
    """A failure's retained frames are observations, never a completed workflow."""
    root = Path(root)
    backend = read(root / "backend.json")
    _backend(backend)
    config = read(root / "config.json")
    NPTConfiguration.from_payload(config)
    failure = root / label / "failure.json"
    if failure.exists():
        result = read(failure)
        require(result["identity"] == _digest(_body(result)), "Changed failure record")
        require(
            result["config"] == config
            and result["backend_identity"] == backend["identity"],
            "Changed failure configuration/backend",
        )
        for k, v in FLAGS.items():
            require(
                type(result.get(k)) is type(v) and result[k] == v, "Changed readiness"
            )
        row = read(root / "result.json")
        if "retained_failure_identity" in row:
            require(
                row["retained_failure_identity"] == result["identity"],
                "Changed retained failure/frame content",
            )
        frames = result["frames"]
        complete = False
    else:
        workflow = completed_history(root / label)
        require(workflow["manifest"]["config"] == config, "Changed workflow config")
        require(
            workflow["manifest"]["backend_identity"] == backend["identity"],
            "Changed backend",
        )
        frames = workflow["frames"]
        complete = True
    require(
        frames and frames[0]["stage"] == "minimize", "Missing minimization observation"
    )
    previous = -1
    for frame in frames:
        require(
            type(frame["step"]) is int and frame["step"] > previous, "Frame ordering"
        )
        previous = frame["step"]
        raw = base64.b64decode(frame["checkpoint_base64"], validate=True)
        require(
            sha256(raw).hexdigest() == frame["checkpoint_sha256"],
            "Changed checkpoint/RNG",
        )
        for key in ("positions_nm", "velocities_nm_ps", "box_vectors_nm"):
            a = np.asarray(frame[key], dtype=float)
            expected = (
                (3, 3)
                if key == "box_vectors_nm"
                else (len(backend["final_graph"]["sites"]), 3)
            )
            require(
                a.shape == expected and np.isfinite(a).all(),
                "Invalid actual frame " + key,
            )
        require(
            frame["measured_pressure_bar"] is None
            and frame["pressure_observation_source"] == "unavailable",
            "Pressure target substituted for observation",
        )
        temperature = (
            2
            * frame["kinetic_energy_kj_mol"]
            / (3 * len(backend["final_graph"]["sites"]) * 0.00831446261815324)
        )
        require(
            frame["measured_temperature_kelvin"] == temperature,
            "Changed measured temperature",
        )
        require(
            frame["volume_nm3"] == float(np.linalg.det(frame["box_vectors_nm"])),
            "Changed actual volume",
        )
        expected_stage = (
            "minimize"
            if frame["step"] == 0
            else "nvt"
            if frame["step"] <= config["nvt_steps"]
            else "npt"
        )
        require(frame["stage"] == expected_stage, "Changed stage")
    return backend, config, frames, complete


def observations(root):
    backend, config, frames, complete = trajectory(root)
    mass = sum(s["mass"] for s in backend["final_graph"]["sites"])
    measured = frame_measurements(frames, mass)
    for frame in measured:
        frame["elapsed_protocol_time_ps"] = frame["step"] * config["timestep_fs"] / 1000
    npt = [f for f in measured if f["stage"] == "npt"]
    policy = EquilibrationPolicy()
    metrics = {
        name: series([f[key] for f in npt], policy)
        for name, key in (
            ("density", "density_g_cm3"),
            ("volume", "volume_nm3"),
            ("temperature", "measured_temperature_kelvin"),
        )
    }
    return _signed(
        dict(
            schema="island.j40-c2-observations.v1",
            trajectory_identity=_digest(frames),
            config_identity=config["identity"],
            backend_identity=backend["identity"],
            complete=complete,
            last_retained_step=frames[-1]["step"],
            npt_frames=len(npt),
            total_retained_frames=len(frames),
            npt_duration_ps=max(0, frames[-1]["step"] - config["nvt_steps"])
            * config["timestep_fs"]
            / 1000,
            observations=measured,
            metrics=metrics,
            pressure_status="pressure_unavailable",
            policy=policy.payload,
            protocol=CorrectiveProtocol().payload,
            scientific_acceptance=False,
            **FLAGS,
        )
    )


def components(root):
    """Single-point component decomposition on each actual retained geometry/box.

    Independent contexts do not advance or change the trajectory's RNG/state.
    Pressure is deliberately not approximated from an incomplete force virial.
    """
    import openmm as mm
    from openmm import unit

    from island._periodic_openmm import build

    backend, _, frames, _ = trajectory(root)
    model = build(backend["numerical_data"], backend["config"], mm)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    groups = {}
    for force in model.getForces():
        groups.setdefault(force.getForceGroup(), []).append(force.getName())
    rows = []
    try:
        for frame in frames:
            context.setPeriodicBoxVectors(
                *[mm.Vec3(*v) for v in frame["box_vectors_nm"]]
            )
            context.setPositions(frame["positions_nm"])
            energies = {
                "+".join(names): float(
                    context.getState(getEnergy=True, groups=1 << group)
                    .getPotentialEnergy()
                    .value_in_unit(unit.kilojoule_per_mole)
                )
                for group, names in groups.items()
            }
            error = abs(sum(energies.values()) - frame["potential_energy_kj_mol"])
            rows.append(
                {
                    "step": frame["step"],
                    "components_kj_mol": energies,
                    "sum_error_kj_mol": error,
                    "passed": error <= 1e-5,
                }
            )
    finally:
        del context, integrator
    return _signed(
        dict(
            schema="island.j40-c2-component-energies.v1",
            trajectory_identity=_digest(frames),
            backend_identity=backend["identity"],
            engine=mm.version.full_version,
            platform="Reference",
            bound_kj_mol=1e-5,
            rows=rows,
            passed=all(x["passed"] for x in rows),
            **FLAGS,
        )
    )


def replica_analysis(roots):
    workflows = []
    for root in roots:
        backend, config, frames, _ = trajectory(root)
        if not any(f["stage"] == "npt" for f in frames):
            return {
                "status": "blocked",
                "reason": "NPT did not execute in every replica",
            }
        workflows.append(
            {
                "manifest": {
                    "config": config,
                    "backend_identity": backend["identity"],
                    "total_mass_dalton": sum(
                        s["mass"] for s in backend["final_graph"]["sites"]
                    ),
                },
                "frames": frames,
            }
        )
    # Independently regenerated AmberTools sources have distinct source/assignment
    # receipts. Never replace those identities to trick the historical comparator.
    if len({w["manifest"]["backend_identity"] for w in workflows}) != 1:
        from math import hypot

        records = [observations(root) for root in roots]
        agreement = {}
        for name in ("density", "volume", "temperature"):
            ms = [r["metrics"][name] for r in records]
            pairs = []
            okay = all(m["status"] == "ok" for m in ms)
            for i, a in enumerate(ms):
                for b in ms[i + 1 :]:
                    delta = abs(a.get("mean", 0) - b.get("mean", 0))
                    se = hypot(a.get("standard_error", 0), b.get("standard_error", 0))
                    pairs.append(
                        {"mean_difference": delta, "combined_standard_error": se}
                    )
                    okay &= delta <= EquilibrationPolicy().replica_z * se
            agreement[name] = {
                "status": "ok" if okay else "insufficient_sampling",
                "comparisons": pairs,
            }
        return _signed(
            dict(
                schema="island.j40-c2-distinct-source-replicas.v1",
                status="blocked",
                reason="independently_regenerated_source_identities_differ",
                backend_identities=[
                    w["manifest"]["backend_identity"] for w in workflows
                ],
                replicas=[r["metrics"] for r in records],
                replica_agreement=agreement,
                scientific_acceptance=False,
                **FLAGS,
            )
        )
    # Partial histories are explicitly marked in the containing matrix.
    return diagnose(workflows)
