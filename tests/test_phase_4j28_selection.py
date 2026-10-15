from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from island.analysis import (
    AnalysisSelection,
    TrajectoryFrame,
    TrajectoryMetadata,
    VerifiedTrajectory,
    calculate_observables,
    load_analysis_record,
    load_verified_trajectory,
    make_analysis_record,
)
from island.exceptions import AnalysisError


def _trajectory(*, periodic=True, bonds=((1, 2),), chain_groups=((1, 2),)):
    metadata = TrajectoryMetadata(
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
        atom_mapping={1: 0, 2: 1},
        box_periodic=(True, True, True) if periodic else (False, False, False),
    )
    frames = tuple(TrajectoryFrame(np.array([[9.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), np.diag([10.0, 10.0, 10.0]), step, float(step)) for step in range(2))
    return VerifiedTrajectory(metadata, frames, "trajectory", {"trajectory": "9" * 64}, bonds, chain_groups)


def test_explicit_h_backbone_selection_is_not_branches_from_hydrogen_degree():
    graph = {
        "graph_identity": "a" * 64,
        "sites": [
            {"id": 1, "atomic_number": 6, "metadata": {"chain_id": "A"}},
            {"id": 2, "atomic_number": 6, "metadata": {"chain_id": "A"}},
            {"id": 3, "atomic_number": 1, "metadata": {"chain_id": "A"}},
            {"id": 4, "atomic_number": 1, "metadata": {"chain_id": "A"}},
        ],
        "bonds": [{"sites": [1, 2]}, {"sites": [1, 3]}, {"sites": [1, 4]}],
        "molecule_membership": {"1": "molecule-A", "2": "molecule-A", "3": "molecule-A", "4": "molecule-A"},
    }
    selection = AnalysisSelection.from_final_graph(graph, atom_mapping_identity="b" * 64)
    assert selection.chains[0]["backbone_path"] == [1, 2]
    assert selection.chains[0]["rg_atom_ids"] == [1, 2, 3, 4]


def test_nonperiodic_six_angstrom_re_does_not_apply_minimum_image():
    trajectory = _trajectory(periodic=False)
    frames = tuple(TrajectoryFrame(np.array([[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]]), frame.box_vectors_angstrom, frame.step, frame.time_ps) for frame in trajectory.frames)
    trajectory = VerifiedTrajectory(trajectory.metadata, frames, trajectory.source_path, trajectory.trajectory_file_hashes, ((1, 2),), ((1, 2),))
    values = calculate_observables(trajectory, ("end_to_end_distance",), block_size=1)
    assert values["end_to_end_distance"]["value"] == pytest.approx(6.0)


def test_loader_retains_final_graph_bonds_and_chain_selection():
    workflow = next(Path("../island-validation/phase4j23/acceptance").glob("test_four_family_minimization_0/min"))
    trajectory = load_verified_trajectory(workflow)
    assert trajectory.topology_bonds
    assert trajectory.molecule_membership
    assert trajectory.selection is not None


def test_selection_identity_is_persisted_and_changes_with_selection():
    trajectory = _trajectory()
    first = AnalysisSelection.from_trajectory(trajectory, rg_atom_ids=(1, 2), endpoint_ids=(1, 2), backbone_path=(1, 2))
    second = AnalysisSelection.from_trajectory(trajectory, rg_atom_ids=(1, 2), endpoint_ids=(1, 2), backbone_path=(1, 2))
    second = replace(second, provenance={"variant": "second"})
    assert first.identity != second.identity
    record = make_analysis_record(trajectory, calculate_observables(trajectory, ("radius_of_gyration",), block_size=1), selection=first)
    assert record.payload["selection"]["identity"] == first.identity


def test_periodic_geometry_requires_authoritative_bonds():
    trajectory = _trajectory(bonds=(), chain_groups=((1, 2),))
    selection = AnalysisSelection.from_trajectory(trajectory, rg_atom_ids=(1, 2), endpoint_ids=(1, 2), backbone_path=(1, 2))
    with pytest.raises(AnalysisError, match="authoritative backbone bonds"):
        VerifiedTrajectory(trajectory.metadata, trajectory.frames, trajectory.source_path, trajectory.trajectory_file_hashes, (), (), {}, selection)


def test_explicit_h_rg_uses_full_molecular_selection_and_backbone_endpoints():
    metadata = replace(_trajectory().metadata, atom_ids=(1, 2, 3, 4), masses_dalton=(12.0, 12.0, 1.0, 1.0), atom_mapping={1: 0, 2: 1, 3: 2, 4: 3})
    frames = tuple(TrajectoryFrame(np.array([[9.0, 0.0, 0.0], [1.0, 0.0, 0.0], [8.5, 0.0, 0.0], [1.5, 0.0, 0.0]]), np.diag([10.0, 10.0, 10.0]), step, float(step)) for step in range(2))
    trajectory = VerifiedTrajectory(metadata, frames, "trajectory", {"trajectory": "9" * 64}, ((1, 2), (1, 3), (2, 4)))
    selection = AnalysisSelection.from_trajectory(trajectory, rg_atom_ids=(1, 2, 3, 4), endpoint_ids=(1, 2), backbone_path=(1, 2))
    values = calculate_observables(trajectory, ("radius_of_gyration", "end_to_end_distance", "re_over_rg"), block_size=1, selection=selection)
    assert values["end_to_end_distance"]["value"] == pytest.approx(2.0)
    assert values["radius_of_gyration"]["value"] == pytest.approx(1.0469736606781692)
    assert values["re_over_rg"]["value"] == pytest.approx(2.0 / 1.0469736606781692)


def test_selection_binding_rejects_recomputed_nested_mutation(tmp_path):
    trajectory = _trajectory()
    selection = AnalysisSelection.from_trajectory(trajectory, rg_atom_ids=(1, 2), endpoint_ids=(1, 2), backbone_path=(1, 2))
    record = make_analysis_record(trajectory, calculate_observables(trajectory, ("radius_of_gyration",), block_size=1), selection=selection)
    path = tmp_path / "analysis.json"
    path.write_bytes(record.to_json())
    assert load_analysis_record(path, expected_selection=selection).identity == record.identity
    changed = dict(record.payload)
    changed["selection"] = dict(changed["selection"])
    changed["selection"]["chains"] = [dict(changed["selection"]["chains"][0], rg_atom_ids=[1])]
    changed["identity"] = None
    import json
    from hashlib import sha256
    body = dict(changed)
    body.pop("identity")
    changed["identity"] = sha256((json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()).hexdigest()
    path.write_text(json.dumps(changed))
    with pytest.raises(AnalysisError):
        load_analysis_record(path, expected_selection=selection)
