import json
import subprocess
import sys
from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest

from island.analysis import (
    calculate_observables,
    make_analysis_record,
    relocate_analysis_record,
    write_analysis_record,
)
from island.exceptions import AnalysisError
from scripts.phase_4j30_fixture import build_multichain_trajectory

FIXTURE = Path("../island-validation/phase4j23/acceptance/test_four_family_exact_relocat0/full")


def _trajectory():
    if not (FIXTURE / "manifest.json").exists():
        pytest.skip(f"retained J23 fixture is unavailable: {FIXTURE}")
    return build_multichain_trajectory(FIXTURE)


def _independent_chain(trajectory, chain, frame_index=0):
    frame = trajectory.frames[frame_index]
    mapping = trajectory.metadata.atom_mapping
    box = np.diag(frame.box_vectors_angstrom)
    coordinates = frame.coordinates_angstrom
    adjacency = {atom_id: set() for atom_id in trajectory.metadata.atom_ids}
    for left, right in trajectory.topology_bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
    path = [int(atom_id) for atom_id in chain["backbone_path"]]
    placed = {path[0]: coordinates[mapping[path[0]]].copy()}
    for left, right in pairwise(path):
        delta = coordinates[mapping[right]] - coordinates[mapping[left]]
        delta -= np.rint(delta / box) * box
        placed[right] = placed[left] + delta
    pending = {int(atom_id) for atom_id in chain["rg_atom_ids"]} - set(placed)
    while pending:
        progressed = False
        for atom_id in sorted(pending):
            parents = [parent for parent in adjacency[atom_id] if parent in placed]
            if len(parents) == 1:
                delta = coordinates[mapping[atom_id]] - coordinates[mapping[parents[0]]]
                delta -= np.rint(delta / box) * box
                placed[atom_id] = placed[parents[0]] + delta
                progressed = True
        if not progressed:
            raise AssertionError("independent reference could not place chain")
        pending = {int(atom_id) for atom_id in chain["rg_atom_ids"]} - set(placed)
    atom_ids = [int(atom_id) for atom_id in chain["rg_atom_ids"]]
    masses = np.asarray([trajectory.metadata.masses_dalton[mapping[atom_id]] for atom_id in atom_ids])
    points = np.asarray([placed[atom_id] for atom_id in atom_ids])
    center = np.average(points, axis=0, weights=masses)
    rg = float(np.sqrt(np.average(np.sum((points - center) ** 2, axis=1), weights=masses)))
    endpoints = chain["endpoint_ids"]
    re = float(np.linalg.norm(placed[int(endpoints[1])] - placed[int(endpoints[0])]))
    return rg, re, re / rg, float(masses.sum())


def test_real_derived_multichain_all_observables_and_independent_values():
    trajectory = _trajectory()
    assert len(trajectory.selection.chains) == 2
    assert len(trajectory.selection.chains[0]["rg_atom_ids"]) != len(trajectory.selection.chains[1]["rg_atom_ids"])
    names = ("radius_of_gyration", "end_to_end_distance", "re_over_rg", "gyration_tensor", "chain_anisotropy", "chain_orientation", "nematic_order", "density", "network_connectivity")
    values = calculate_observables(trajectory, names, block_size=2)
    assert all(values[name]["status"] == "ok" for name in names)
    references = [[_independent_chain(trajectory, chain, frame_index) for frame_index in range(trajectory.frame_count)] for chain in trajectory.selection.chains]
    assert values["radius_of_gyration"]["per_chain"] == pytest.approx([np.mean([reference[0] for reference in chain_references]) for chain_references in references])
    assert values["end_to_end_distance"]["value"] == pytest.approx(np.mean([reference[1] for chain_references in references for reference in chain_references[:6]]))
    assert values["re_over_rg"]["value"] == pytest.approx(np.mean([reference[2] for chain_references in references for reference in chain_references[:6]]))


def test_translation_and_mass_weighted_aggregation():
    trajectory = _trajectory()
    baseline = calculate_observables(trajectory, ("radius_of_gyration", "end_to_end_distance", "re_over_rg"), block_size=2)
    translated_frames = []
    second_ids = set(trajectory.selection.chains[1]["rg_atom_ids"])
    for frame in trajectory.frames:
        coordinates = frame.coordinates_angstrom.copy()
        for atom_id in second_ids:
            coordinates[trajectory.metadata.atom_mapping[atom_id]] += np.asarray((30.0, 0.0, 0.0))
        translated_frames.append(replace(frame, coordinates_angstrom=coordinates))
    translated = replace(trajectory, frames=tuple(translated_frames))
    translated_values = calculate_observables(translated, ("radius_of_gyration", "end_to_end_distance", "re_over_rg"), block_size=2)
    for name in baseline:
        assert translated_values[name]["value"] == pytest.approx(baseline[name]["value"])
    mass_selection = replace(trajectory.selection, chain_aggregation="mass_weighted")
    weighted = calculate_observables(trajectory, ("radius_of_gyration", "re_over_rg"), block_size=2, selection=mass_selection)
    references = [[_independent_chain(trajectory, chain, frame_index) for frame_index in range(trajectory.frame_count)] for chain in trajectory.selection.chains]
    weights = np.asarray([chain_references[0][3] for chain_references in references])
    weights /= weights.sum()
    expected = np.mean([float(np.dot(weights, [references[chain_index][frame_index][0] for chain_index in range(len(references))])) for frame_index in range(6)])
    assert weighted["radius_of_gyration"]["value"] == pytest.approx(expected)


def test_multichain_integrity_and_separate_process_reconstruction(tmp_path):
    trajectory = _trajectory()
    values = calculate_observables(trajectory, ("radius_of_gyration", "end_to_end_distance", "re_over_rg", "density"), block_size=2)
    record = make_analysis_record(trajectory, values, provenance={"phase": "4J30"})
    source = write_analysis_record(record, tmp_path / "analysis.json")
    relocated = relocate_analysis_record(source, tmp_path / "relocated.json")
    code = """
import json, sys
from island.analysis import load_analysis_record
r = load_analysis_record(sys.argv[1], expected_selection_identity=sys.argv[2])
print(json.dumps(r.payload, sort_keys=True, separators=(\",\", \":\")))
"""
    result = subprocess.run([sys.executable, "-c", code, str(relocated), trajectory.selection.identity], cwd=Path.cwd(), check=True, capture_output=True, text=True)
    assert result.stdout.strip() == json.dumps(record.payload, sort_keys=True, separators=(",", ":"))
    reordered = replace(trajectory.selection, chains=tuple(reversed(trajectory.selection.chains)))
    assert reordered.identity != trajectory.selection.identity
    missing_pair = frozenset((trajectory.selection.chains[1]["backbone_path"][0], trajectory.selection.chains[1]["backbone_path"][1]))
    with pytest.raises(AnalysisError, match="authoritative backbone bonds"):
        replace(trajectory, topology_bonds=tuple(bond for bond in trajectory.topology_bonds if frozenset(bond) != missing_pair))
