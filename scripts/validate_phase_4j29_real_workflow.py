"""Run the real J23-to-analysis acceptance path with an independent control."""

from __future__ import annotations

import argparse
import json
from itertools import pairwise
from pathlib import Path

import numpy as np

from island.analysis import (
    calculate_observables,
    load_verified_trajectory,
    make_analysis_record,
    relocate_analysis_record,
    write_analysis_record,
)


def _independent_first_frame(trajectory):
    selection = trajectory.selection
    chain = selection.chains[0]
    mapping = trajectory.metadata.atom_mapping
    frame = trajectory.frames[0]
    box = np.diag(frame.box_vectors_angstrom)
    coordinates = frame.coordinates_angstrom
    bonds = {int(atom): set() for atom in trajectory.metadata.atom_ids}
    for left, right in trajectory.topology_bonds:
        bonds[left].add(right)
        bonds[right].add(left)
    path = [int(atom) for atom in chain["backbone_path"]]
    placed = {path[0]: coordinates[mapping[path[0]]].copy()}
    for left, right in pairwise(path):
        delta = coordinates[mapping[right]] - coordinates[mapping[left]]
        delta -= np.rint(delta / box) * box
        placed[right] = placed[left] + delta
    pending = {int(atom) for atom in chain["rg_atom_ids"]} - set(placed)
    while pending:
        progressed = False
        for atom in sorted(pending):
            parents = [parent for parent in bonds[atom] if parent in placed]
            if len(parents) != 1:
                continue
            parent = parents[0]
            delta = coordinates[mapping[atom]] - coordinates[mapping[parent]]
            delta -= np.rint(delta / box) * box
            placed[atom] = placed[parent] + delta
            progressed = True
        if not progressed:
            raise RuntimeError("independent reference could not place all selected atoms")
        pending = {int(atom) for atom in chain["rg_atom_ids"]} - set(placed)
    atom_ids = [int(atom) for atom in chain["rg_atom_ids"]]
    masses = np.asarray([trajectory.metadata.masses_dalton[mapping[atom]] for atom in atom_ids], dtype=float)
    points = np.asarray([placed[atom] for atom in atom_ids], dtype=float)
    center = np.average(points, axis=0, weights=masses)
    rg = float(np.sqrt(np.average(np.sum((points - center) ** 2, axis=1), weights=masses)))
    endpoints = chain["endpoint_ids"]
    re = float(np.linalg.norm(placed[int(endpoints[1])] - placed[int(endpoints[0])]))
    return rg, re, re / rg


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", type=Path)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--relocated-record", type=Path, required=True)
    parser.add_argument("--block-size", type=int, default=2)
    args = parser.parse_args()
    trajectory = load_verified_trajectory(args.workflow)
    names = ("radius_of_gyration", "end_to_end_distance", "re_over_rg", "gyration_tensor", "chain_anisotropy", "chain_orientation", "nematic_order", "density", "network_connectivity")
    observables = calculate_observables(trajectory, names, block_size=args.block_size)
    required = {name for name in names if observables[name]["status"] != "ok"}
    if required:
        raise RuntimeError(f"real J23 workflow did not produce requested observables: {sorted(required)}")
    record = make_analysis_record(trajectory, observables, provenance={"phase": "4J29", "workflow": str(args.workflow)}, evidence_references=("docs/phase_4j29.md",))
    write_analysis_record(record, args.record)
    relocate_analysis_record(args.record, args.relocated_record)
    reference = _independent_first_frame(trajectory)
    result = {
        "record_identity": record.identity,
        "selection_identity": trajectory.selection.identity,
        "trajectory_identity": trajectory.identity,
        "workflow_identity": trajectory.metadata.workflow_identity,
        "observables": observables,
        "independent_first_frame": {"radius_of_gyration": reference[0], "end_to_end_distance": reference[1], "re_over_rg": reference[2]},
        "fixture": {"final_graph_identity": trajectory.metadata.final_graph_identity, "packing_identity": trajectory.metadata.packing_identity, "force_field_family": trajectory.metadata.force_field_family, "force_field_source_identity": trajectory.metadata.force_field_source_identity, "typed_graph_identity": trajectory.metadata.typed_graph_identity, "charge_identity": trajectory.metadata.charge_identity, "checkpoint_identity": trajectory.metadata.checkpoint_identity},
    }
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
