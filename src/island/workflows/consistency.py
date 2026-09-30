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
    masses = mass_inventory(system)
    if dict(initialization.masses) != masses:
        raise WorkflowError("Initialization masses differ from starting system")
    parameter = preparation.imported_result.content_signature()
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
