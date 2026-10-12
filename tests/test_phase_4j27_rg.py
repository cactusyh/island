from dataclasses import replace

import numpy as np
import pytest

from island.analysis import (
    TrajectoryFrame,
    TrajectoryMetadata,
    VerifiedTrajectory,
    calculate_observables,
)
from island.exceptions import AnalysisError


def test_periodic_linear_chain_rg_uses_minimum_image_geometry():
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
        atom_ids=(1, 2),
        masses_dalton=(1.0, 1.0),
        atom_mapping={1: 0, 2: 1},
    )
    frames = tuple(
        TrajectoryFrame(
            np.array([[9.0, 1.0, 1.0], [1.0, 1.0, 1.0]]),
            np.diag([10.0, 10.0, 10.0]),
            step,
            float(step),
        )
        for step in range(2)
    )
    trajectory = VerifiedTrajectory(metadata, frames, "trajectory", {"trajectory": "9" * 64}, ((1, 2),), ((1, 2),))
    values = calculate_observables(trajectory, ("radius_of_gyration", "end_to_end_distance", "re_over_rg"), block_size=1)
    assert values["end_to_end_distance"]["value"] == 2.0
    assert values["re_over_rg"]["value"] == 2.0
    assert values["radius_of_gyration"]["value"] == 1.0


def _make_trajectory(coordinates, groups, bonds, boxes=None, second_coordinates=None):
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
        atom_ids=tuple(range(1, len(coordinates) + 1)),
        masses_dalton=tuple([1.0] * len(coordinates)),
        atom_mapping={index: index - 1 for index in range(1, len(coordinates) + 1)},
    )
    boxes = boxes or [np.diag([10.0, 10.0, 10.0])] * 2
    second_coordinates = coordinates if second_coordinates is None else second_coordinates
    frames = tuple(TrajectoryFrame(np.asarray(frame), np.asarray(boxes[index]), index, float(index)) for index, frame in enumerate((coordinates, second_coordinates)))
    return VerifiedTrajectory(metadata, frames, "trajectory", {"trajectory": "9" * 64}, tuple(bonds), tuple(groups))


def test_translated_chains_and_different_lengths_use_equal_chain_weight():
    trajectory = _make_trajectory(
        [[9.0, 1.0, 1.0], [1.0, 1.0, 1.0], [4.0, 4.0, 1.0], [5.0, 4.0, 1.0], [6.0, 4.0, 1.0]],
        ((1, 2), (3, 4, 5)),
        ((1, 2), (3, 4), (4, 5)),
    )
    values = calculate_observables(trajectory, ("radius_of_gyration", "re_over_rg"), block_size=1)
    assert values["radius_of_gyration"]["status"] == "ok"
    assert len(values["radius_of_gyration"]["per_chain"]) == 2
    assert values["re_over_rg"]["status"] == "ok"
    assert len(values["re_over_rg"]["per_chain"]) == 2


def test_multiple_frames_recompute_chain_geometry_before_blocks():
    trajectory = _make_trajectory(
        [[9.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
        ((1, 2),),
        ((1, 2),),
        second_coordinates=[[8.0, 1.0, 1.0], [2.0, 1.0, 1.0]],
    )
    values = calculate_observables(trajectory, ("radius_of_gyration", "end_to_end_distance", "re_over_rg"), block_size=1)
    assert values["radius_of_gyration"]["value"] == pytest.approx(1.5)
    assert values["end_to_end_distance"]["value"] == pytest.approx(3.0)
    assert values["re_over_rg"]["value"] == pytest.approx(2.0)


def test_periodic_system_without_chain_groups_is_blocked_for_chain_rg():
    trajectory = _make_trajectory([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]], (), ((1, 2),))
    values = calculate_observables(trajectory, ("radius_of_gyration",), block_size=1)
    assert values["radius_of_gyration"]["status"] == "blocked"
    assert "chain_groups" in values["radius_of_gyration"]["failure"]


def test_nonperiodic_whole_system_rg_remains_supported():
    trajectory = _make_trajectory([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], (), ((1, 2),))
    nonperiodic = VerifiedTrajectory(replace(trajectory.metadata, box_periodic=(False, False, False)), trajectory.frames, trajectory.source_path, trajectory.trajectory_file_hashes, trajectory.topology_bonds, ())
    values = calculate_observables(nonperiodic, ("radius_of_gyration",), block_size=1)
    assert values["radius_of_gyration"]["status"] == "ok"
    assert values["radius_of_gyration"]["value"] == pytest.approx(1.0)


def test_zero_rg_and_triclinic_rg_are_blocked():
    collapsed = _make_trajectory([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], ((1, 2),), ((1, 2),))
    collapsed_values = calculate_observables(collapsed, ("radius_of_gyration", "re_over_rg"), block_size=1)
    assert collapsed_values["radius_of_gyration"]["status"] == "blocked"
    assert collapsed_values["re_over_rg"]["status"] == "blocked"
    triclinic = _make_trajectory([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]], ((1, 2),), ((1, 2),), boxes=[np.array([[10.0, 1.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])] * 2)
    values = calculate_observables(triclinic, ("radius_of_gyration", "re_over_rg"), block_size=1)
    assert values["radius_of_gyration"]["status"] == "blocked"
    assert "orthorhombic" in values["radius_of_gyration"]["failure"]


def test_nonfinite_coordinates_are_rejected_before_rg():
    with pytest.raises(AnalysisError, match="finite"):
        TrajectoryFrame(np.array([[np.nan, 0.0, 0.0], [1.0, 0.0, 0.0]]), np.diag([10.0, 10.0, 10.0]), 0, 0.0)
