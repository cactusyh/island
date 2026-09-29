"""Reconstructed frames must not bypass stable-ID, arithmetic or state checks."""

from copy import deepcopy
from dataclasses import replace

import pytest
from test_dynamics import VELOCITIES, Harmonic, molecule, options

from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import DynamicsOptions, run_nve
from island.dynamics.models import velocity_hash
from island.exceptions import DynamicsInputError, InvalidDynamicsResultError


@pytest.mark.parametrize(
    "changes",
    [
        {"timestep_fs": 0},
        {"timestep_fs": float("nan")},
        {"steps": True},
        {"steps": 0},
        {"max_evaluations": 0},
        {"max_frames": 0},
        {"recording_interval": -1},
        {"max_energy_deviation": float("inf")},
    ],
)
def test_invalid_options(changes):
    fields = {"timestep_fs": 0.1, "steps": 2, "max_evaluations": 4, "max_frames": 3}
    fields.update(changes)
    with pytest.raises(DynamicsInputError):
        DynamicsOptions(**fields)


@pytest.mark.parametrize(
    "field,value",
    [
        ("time_ps", 999.0),
        ("kinetic_energy", 999.0),
        ("total_energy", float("nan")),
        ("energy_deviation", 999.0),
        ("coordinate_fingerprint", "wrong"),
        ("velocity_fingerprint", "wrong"),
        ("evaluation", None),
        ("step", True),
        ("velocities", {17: (0, 0, 0)}),
    ],
)
def test_bad_frame_fields(field, value):
    system = molecule()
    result = run_nve(system, Harmonic(), VELOCITIES, options(2))
    with pytest.raises(InvalidDynamicsResultError):
        altered = replace(result.frames[-1], **{field: value})
        replace(result, frames=(*result.frames[:-1], altered)).to_system(
            system, allow_incomplete=True
        )


def test_missing_coordinates_and_forces_with_matching_fingerprints():
    system = molecule()
    before = deepcopy(system.to_dict())
    result = run_nve(system, Harmonic(), VELOCITIES, options(1))
    last = result.final_state
    reduced_x, reduced_v = {17: last.coordinates[17]}, {17: last.velocities[17]}
    record = replace(
        last.evaluation,
        forces={17: last.evaluation.forces[17]},
        coordinate_fingerprint=coordinate_hash(reduced_x),
    )
    altered = replace(
        last,
        coordinates=reduced_x,
        velocities=reduced_v,
        evaluation=record,
        coordinate_fingerprint=coordinate_hash(reduced_x),
        velocity_fingerprint=velocity_hash(reduced_v),
    )
    with pytest.raises(InvalidDynamicsResultError):
        replace(
            result, frames=(result.initial_state, altered), final_evaluation=record
        ).to_system(system, allow_incomplete=True)
    assert system.to_dict() == before


def test_bad_order_duplicate_frames_and_velocity_identity():
    result = run_nve(molecule(), Harmonic(), VELOCITIES, options(2))
    for frames in (
        result.frames[::-1],
        (result.frames[0], *result.frames),
        (result.frames[0], result.frames[-1]),
    ):
        with pytest.raises(InvalidDynamicsResultError):
            replace(result, frames=frames).validate_integrity()
    for changes in (
        {"masses": {17: 3.0, 91: 8.0}},
        {"options": replace(result.options, timestep_fs=0.25)},
        {
            "final_evaluation": replace(
                result.final_evaluation, model_fingerprint="different"
            )
        },
    ):
        with pytest.raises(InvalidDynamicsResultError):
            replace(result, **changes).validate_integrity()


def test_changed_final_verification_does_not_replace_accepted_state():
    class FinalChanged(Harmonic):
        def evaluate(self, *args, **kwargs):
            record = super().evaluate(*args, **kwargs)
            if self.calls == 4:
                return replace(
                    record,
                    potential_energy=record.potential_energy + 1.0,
                    energy_components={"harmonic": record.potential_energy + 1.0},
                )
            return record

    regular = run_nve(molecule(), Harmonic(), VELOCITIES, options(2))
    failed = run_nve(molecule(), FinalChanged(), VELOCITIES, options(2))
    assert not failed.completed and failed.completed_steps == 2
    assert failed.termination_reason == "final_evaluation_failed"
    assert failed.final_state.coordinates == regular.final_state.coordinates
    assert failed.final_state.velocities == regular.final_state.velocities
    assert failed.final_state.total_energy == regular.final_state.total_energy


def test_unavailable_backend_has_domain_exception():
    from island.exceptions import DynamicsUnavailableError, EvaluationUnavailableError

    class Unavailable(Harmonic):
        def evaluate(self, *args, **kwargs):
            if self.calls == 1:
                raise EvaluationUnavailableError("backend dependency disappeared")
            return super().evaluate(*args, **kwargs)

    with pytest.raises(DynamicsUnavailableError, match="backend dependency"):
        run_nve(molecule(), Unavailable(), VELOCITIES, options(2))
