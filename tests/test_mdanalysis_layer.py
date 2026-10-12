import json

import numpy as np
import pytest

from island.analysis import (
    TrajectoryFrame,
    TrajectoryMetadata,
    VerifiedTrajectory,
    calculate_observables,
    load_analysis_record,
    make_analysis_record,
    relocate_analysis_record,
    validate_physical_gates,
    write_analysis_record,
)
from island.exceptions import AnalysisError


def _trajectory():
    digest = "a" * 64
    metadata = TrajectoryMetadata(
        final_graph_identity=digest,
        packing_identity="b" * 64,
        force_field_family="GAFF2",
        force_field_source_identity="c" * 64,
        typed_graph_identity="d" * 64,
        charge_identity="e" * 64,
        workflow_identity="f" * 64,
        checkpoint_identity="1" * 64,
        atom_ids=(10, 20, 30),
        masses_dalton=(12.0, 12.0, 16.0),
        atom_mapping={10: 0, 20: 1, 30: 2},
    )
    frames = tuple(
        TrajectoryFrame(
            np.array([[1.0 + shift, 1.0, 1.0], [2.0 + shift, 1.0, 1.0], [1.0 + shift, 2.0, 1.0]]),
            np.diag([20.0, 20.0, 20.0]),
            shift,
            float(shift),
        )
        for shift in range(4)
    )
    return VerifiedTrajectory(metadata, frames, "workflow", {"frame": digest}, chain_groups=((10, 30),))


def test_analysis_record_binds_observable_and_relocates(tmp_path):
    trajectory = _trajectory()
    values = calculate_observables(trajectory, ("density", "radius_of_gyration", "glass_transition"), block_size=2)
    assert values["density"]["status"] == "ok"
    assert values["glass_transition"]["status"] == "blocked"
    record = make_analysis_record(trajectory, values)
    path = write_analysis_record(record, tmp_path / "record.json")
    moved = relocate_analysis_record(path, tmp_path / "moved" / "record.json")
    assert load_analysis_record(moved, expected_identity=record.identity).identity == record.identity


def test_re_over_rg_can_be_requested_without_endpoint_observable():
    values = calculate_observables(_trajectory(), ("re_over_rg",), block_size=2)
    assert set(values) == {"re_over_rg"}
    assert values["re_over_rg"]["status"] == "ok"
    assert values["re_over_rg"]["value"] > 0


def test_re_over_rg_missing_endpoints_is_explicitly_blocked():
    trajectory = _trajectory()
    trajectory = VerifiedTrajectory(trajectory.metadata, trajectory.frames, "workflow", {"frame": "a" * 64})
    values = calculate_observables(trajectory, ("re_over_rg",), block_size=2)
    assert values["re_over_rg"]["status"] == "blocked"
    assert "chain_groups" in values["re_over_rg"]["failure"]


def test_re_over_rg_zero_radius_is_explicitly_blocked():
    trajectory = _trajectory()
    frames = tuple(
        TrajectoryFrame(np.ones((3, 3)), frame.box_vectors_angstrom, frame.step, frame.time_ps)
        for frame in trajectory.frames
    )
    collapsed = VerifiedTrajectory(trajectory.metadata, frames, "workflow", {"frame": "a" * 64}, chain_groups=((10, 30),))
    values = calculate_observables(collapsed, ("re_over_rg",), block_size=2)
    assert values["re_over_rg"]["status"] == "blocked"
    assert "zero" in values["re_over_rg"]["failure"]


def test_end_to_end_aggregates_multiple_chains_without_mutation():
    trajectory = _trajectory()
    original = trajectory.frames[0].coordinates_angstrom.copy()
    multi = VerifiedTrajectory(
        trajectory.metadata,
        trajectory.frames,
        "workflow",
        {"frame": "a" * 64},
        chain_groups=((10, 20), (20, 30)),
    )
    values = calculate_observables(multi, ["end_to_end_distance", "re_over_rg"], block_size=2)
    assert values["end_to_end_distance"]["status"] == "ok"
    assert values["re_over_rg"]["status"] == "ok"
    np.testing.assert_array_equal(trajectory.frames[0].coordinates_angstrom, original)


def test_nonfinite_frame_and_empty_chain_are_rejected_or_blocked():
    trajectory = _trajectory()
    with pytest.raises(AnalysisError, match="finite"):
        TrajectoryFrame(np.array([[np.nan, 0.0, 0.0]]), np.diag([20.0, 20.0, 20.0]), 0, 0.0)
    empty = VerifiedTrajectory(
        trajectory.metadata,
        trajectory.frames,
        "workflow",
        {"frame": "a" * 64},
        chain_groups=((),),
    )
    values = calculate_observables(empty, ("end_to_end_distance",), block_size=2)
    assert values["end_to_end_distance"]["status"] == "blocked"


def test_analysis_serialization_is_deterministic_and_preserves_blocked_values():
    trajectory = _trajectory()
    values = calculate_observables(trajectory, ("re_over_rg", "glass_transition"), block_size=2)
    record = make_analysis_record(trajectory, values)
    assert record.to_json() == record.to_json()
    assert record.payload["observables"]["glass_transition"]["status"] == "blocked"


def test_physical_gates_require_explicit_convergence_and_readiness():
    trajectory = _trajectory()
    observables = calculate_observables(trajectory, ("density_convergence",), block_size=2)
    gates = validate_physical_gates(trajectory, observables)
    assert gates["density_plateau"]["status"] == "blocked"
    assert gates["energy_stability"]["status"] == "blocked"
    assert gates["topology_preservation"]["status"] == "blocked"
    assert gates["aggregate"]["status"] == "blocked"


def test_target_thermostat_values_are_not_reported_as_measured_temperature():
    trajectory = _trajectory()
    frame = TrajectoryFrame(
        trajectory.frames[0].coordinates_angstrom,
        trajectory.frames[0].box_vectors_angstrom,
        0,
        0.0,
        300.0,
        1.0,
    )
    second = TrajectoryFrame(
        trajectory.frames[1].coordinates_angstrom,
        trajectory.frames[1].box_vectors_angstrom,
        1,
        1.0,
        300.0,
        1.0,
    )
    measured = VerifiedTrajectory(trajectory.metadata, (frame, second), "workflow", {"frame": "a" * 64})
    values = calculate_observables(measured, ("temperature_density_volume",), block_size=1)
    assert values["temperature_density_volume"]["status"] == "blocked"
    assert "target thermostat" in values["temperature_density_volume"]["failure"]


def test_network_connectivity_and_nematic_order_use_explicit_metadata():
    trajectory = _trajectory()
    metadata = TrajectoryMetadata(
        **{**trajectory.metadata.to_dict(), "atom_mapping": {10: 0, 20: 1, 30: 2}, "source_metadata": {"crosslink_count": 1, "crosslink_density_convention": "per_mass_dalton"}}
    )
    enriched = VerifiedTrajectory(metadata, trajectory.frames, "workflow", {"frame": "a" * 64}, topology_bonds=((10, 20), (20, 30)), chain_groups=((10, 30),))
    values = calculate_observables(enriched, ("network_connectivity", "crosslink_density", "nematic_order", "chain_orientation"), block_size=2)
    assert values["network_connectivity"]["status"] == "ok"
    assert values["crosslink_density"]["status"] == "ok"
    assert values["nematic_order"]["status"] == "ok"
    assert values["chain_orientation"]["status"] == "ok"


def test_triclinic_rdf_is_explicitly_blocked():
    trajectory = _trajectory()
    frames = tuple(TrajectoryFrame(frame.coordinates_angstrom, np.array([[20.0, 1.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]]), frame.step, frame.time_ps) for frame in trajectory.frames)
    triclinic = VerifiedTrajectory(trajectory.metadata, frames, "workflow", {"frame": "a" * 64})
    values = calculate_observables(triclinic, ("radial_distribution_function",), block_size=2)
    assert values["radial_distribution_function"]["status"] == "blocked"
    assert "orthorhombic" in values["radial_distribution_function"]["failure"]


def test_analysis_record_tamper_rejected(tmp_path):
    trajectory = _trajectory()
    record = make_analysis_record(trajectory, calculate_observables(trajectory, ("density",), block_size=2))
    path = write_analysis_record(record, tmp_path / "record.json")
    payload = json.loads(path.read_text())
    payload["observables"]["density"]["value"] += 1
    path.write_text(json.dumps(payload))
    with pytest.raises(AnalysisError, match="checksum/identity"):
        load_analysis_record(path)


def test_mismatched_mapping_rejected():
    with pytest.raises(AnalysisError, match="bijection"):
        TrajectoryMetadata(
            final_graph_identity="a" * 64,
            packing_identity="b" * 64,
            force_field_family="GAFF2",
            force_field_source_identity="c" * 64,
            typed_graph_identity="d" * 64,
            charge_identity="e" * 64,
            workflow_identity="f" * 64,
            checkpoint_identity="1" * 64,
            atom_ids=(1, 2),
            masses_dalton=(1.0, 1.0),
            atom_mapping={1: 0, 2: 2},
        )
