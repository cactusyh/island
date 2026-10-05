"""Relationships between independently valid saved records; no force evaluation."""

from dataclasses import fields
from math import isclose

from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import _checkpoint_data as data
from island.dynamics.integrity import verify_final
from island.dynamics.thermal import mass_inventory
from island.exceptions import WorkflowError


def compatible_boundary(saved, recomputed):
    """Exact synchronized state; existing fresh-verification tolerance for evaluations."""
    try:
        saved.validate_integrity()
        recomputed.validate_integrity()
        if type(saved) is not type(recomputed):
            raise ValueError("Different frame types")
        for field in fields(saved):
            if field.name not in {
                "evaluation",
                "total_energy",
                "energy_deviation",
                "kinetic_energy",
                "instantaneous_temperature_kelvin",
            } and getattr(saved, field.name) != getattr(recomputed, field.name):
                raise ValueError(f"Changed boundary {field.name}")
        # Same arithmetic tolerances used for kinetic/thermal frame validation.
        for name in ("kinetic_energy", "instantaneous_temperature_kelvin"):
            if hasattr(saved, name) and not isclose(
                getattr(saved, name),
                getattr(recomputed, name),
                rel_tol=1e-12,
                abs_tol=1e-10,
            ):
                raise ValueError(f"Changed boundary {name}")
        for name in ("coordinate_fingerprint", "evaluation_fingerprint"):
            if getattr(saved.evaluation, name) != getattr(recomputed.evaluation, name):
                raise ValueError(f"Changed evaluation {name}")
        verify_final(saved.evaluation, recomputed.evaluation)
    except Exception as error:
        raise WorkflowError(
            f"Conflicting adjacent segment boundary: {error}"
        ) from error


def validate_trajectory(config, segments, setup):
    """Bind original state to durable initialization, minimum and parameter identity."""
    system, preparation, initialization, minimum = setup
    validate_prepared_trajectory(
        config,
        segments,
        system,
        initialization,
        minimum,
        preparation.imported_result.content_signature(),
    )


def validate_prepared_trajectory(
    config, segments, system, initialization, minimum, parameter
):
    """Backend-independent setup/trajectory relationships; native parameter identity supplied."""
    masses = mass_inventory(system)
    if dict(initialization.masses) != masses:
        raise WorkflowError("Initialization masses differ from starting system")
    if minimum.final_evaluation.parameter_fingerprint != parameter:
        raise WorkflowError("Minimum parameter identity differs from preparation")
    coords = {s: tuple(system.coordinates.get(s)) for s in sorted(masses)}
    expected_physical = data.physical(config.langevin(config.segment_steps))
    origin = None
    previous_payload = None
    for segment in segments:
        p = segment.payload
        if (
            p["physical"] != expected_physical
            or p["masses"] != data.rows(masses)
            or p["system_fingerprint"] != data.compatibility(system)
        ):
            raise WorkflowError(
                "Trajectory physical/system/mass identity differs from setup"
            )
        initial = data.frame_from(p["origin"]["initial_state"], True)
        if (
            dict(initial.coordinates) != coords
            or initial.coordinate_fingerprint != coordinate_hash(coords)
            or dict(initial.velocities) != dict(initialization.velocities)
            or initial.velocity_fingerprint != initialization.velocity_fingerprint
        ):
            raise WorkflowError(
                "Trajectory origin differs from saved starting coordinates/initialization"
            )
        verify_final(minimum.final_evaluation, initial.evaluation)
        if origin is not None and p["origin"] != origin:
            raise WorkflowError("Original trajectory identity changed")
        if previous_payload is not None:
            if (
                p["lineage"][:-1] != previous_payload["lineage"]
                or p["environment"] != previous_payload["environment"]
            ):
                raise WorkflowError("Trajectory lineage/environment changed")
        else:
            if len(p["lineage"]) != 1:
                raise WorkflowError("First workflow segment has prior lineage")
            compatible_boundary(initial, segment.frames[0])
        origin = p["origin"]
        previous_payload = p


def validate_prepared_attempt(config, segment, system, initialization, minimum, parent):
    """Bind an attempt's authoritative inputs, without accepting its failed observations.

    Native segment integrity validates local accounting and the resulting RNG state.
    A rejected attempt does not become the parent of any later attempt.
    """
    from island.dynamics.integrity import same_model
    from island.dynamics.thermal import _rng

    p = segment.payload
    prior = None if parent is None else parent.payload
    start = 0 if parent is None else parent.absolute_step
    row = p["lineage"][-1]
    expected_options = {
        "steps": min(config.segment_steps, config.total_steps - start),
        "max_evaluations": config.max_evaluations_per_segment,
        "max_frames": config.max_frames_per_segment,
        "recording_interval": config.recording_interval,
    }
    if (
        start >= config.total_steps
        or row["start_step"] != start
        or row["options"] != expected_options
    ):
        raise WorkflowError(
            "Dynamics attempt start/budgets differ from authoritative boundary"
        )
    masses = mass_inventory(system)
    if (
        p["physical"] != data.physical(config.langevin(config.segment_steps))
        or p["masses"] != data.rows(masses)
        or p["system_fingerprint"] != data.compatibility(system)
    ):
        raise WorkflowError(
            "Dynamics attempt physical/system/mass identity differs from setup"
        )
    first = segment.frames[0]
    if first.step != start or first.time_ps != start * (config.timestep_fs * 0.001):
        raise WorkflowError("Dynamics attempt boundary step/time differs")
    if prior is None:
        if len(p["lineage"]) != 1 or row["parent_checksum"] is not None:
            raise WorkflowError("First attempt has unrelated parent lineage")
        coords = {s: tuple(system.coordinates.get(s)) for s in sorted(masses)}
        origin = data.frame_from(p["origin"]["initial_state"], True)
        for frame in (origin, first):
            if dict(frame.coordinates) != coords or dict(frame.velocities) != dict(
                initialization.velocities
            ):
                raise WorkflowError(
                    "First attempt origin differs from saved initialization/setup"
                )
            # No numerical comparison with the minimum here: a failed backend
            # observation may disagree. Accepted segments retain verify_final checks.
            if frame.evaluation is not None:
                same_model(minimum.final_evaluation, frame.evaluation)
        counters = {"evaluations": 0, "random_steps": 0, "normal_draws": 0}
        initial_rng = _rng(config.thermostat_seed).bit_generator.state
    else:
        if (
            p["origin"] != prior["origin"]
            or p["lineage"][:-1] != prior["lineage"]
            or row["parent_checksum"] != parent.content_checksum
            or p["environment"] != prior["environment"]
        ):
            raise WorkflowError(
                "Dynamics attempt origin/parent lineage differs from accepted boundary"
            )
        # Startup failure restores the saved frame; successful startup may have
        # a numerically tolerated fresh evaluation. Failed final/trial observations
        # are deliberately not compared with the input boundary.
        compatible_boundary(data.frame_from(prior["state"], True), first)
        counters = prior["counters"]
        initial_rng = prior["rng"]["state"]
        if any(p["rng"][k] != prior["rng"][k] for k in ("algorithm", "numpy_version")):
            raise WorkflowError("Dynamics attempt RNG identity differs from parent")
    steps = segment.final_state.step - start
    random_steps = steps + int(p["diagnostic"]["failure_stage"] == "trial")
    expected_counters = {
        "evaluations": counters["evaluations"] + segment.evaluations,
        "random_steps": counters["random_steps"] + random_steps,
        "normal_draws": counters["normal_draws"] + 3 * len(masses) * random_steps,
    }
    if p["counters"] != expected_counters:
        raise WorkflowError(
            "Dynamics attempt cumulative counters differ from accepted boundary"
        )
    # No draws is the case where the stored input/output RNG states must be equal.
    # For consumed draws the existing native segment contract validates the complete
    # output PCG64 state and draw counts; no trajectory or RNG replay is introduced.
    if random_steps == 0 and p["rng"]["state"] != initial_rng:
        raise WorkflowError("Dynamics attempt changed RNG without consuming draws")
