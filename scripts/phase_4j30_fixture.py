"""Build a bounded multi-chain analysis fixture from a verified J23 workflow."""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.analysis import (
    AnalysisSelection,
    TrajectoryFrame,
    TrajectoryMetadata,
    VerifiedTrajectory,
    load_verified_trajectory,
)


def _digest(value) -> str:
    return sha256((json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()).hexdigest()


def build_multichain_trajectory(base_path: str | Path) -> VerifiedTrajectory:
    base = load_verified_trajectory(base_path)
    selected_second = (1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15)
    second_ids = {old: 100 + index for index, old in enumerate(selected_second, 1)}
    chain1_ids = tuple(base.metadata.atom_ids)
    second_id_values = tuple(second_ids.values())
    atom_ids = tuple(value for pair in zip(chain1_ids[: len(second_id_values)], second_id_values) for value in pair) + chain1_ids[len(second_id_values):]
    base_mapping = base.metadata.atom_mapping
    chain1_masses = {old: base.metadata.masses_dalton[base_mapping[old]] for old in chain1_ids}
    chain2_masses = {new: chain1_masses[old] for old, new in second_ids.items()}
    masses = tuple(chain1_masses[atom_id] if atom_id in chain1_masses else chain2_masses[atom_id] for atom_id in atom_ids)
    mapping = {atom_id: index for index, atom_id in enumerate(atom_ids)}
    bonds = list(base.topology_bonds)
    bonds.extend((second_ids[left], second_ids[right]) for left, right in base.topology_bonds if left in second_ids and right in second_ids)
    membership = {atom_id: "chain-0" for atom_id in chain1_ids}
    membership.update({atom_id: "chain-1" for atom_id in second_ids.values()})
    frames = []
    for frame in base.frames:
        old_coordinates = {atom_id: frame.coordinates_angstrom[base_mapping[atom_id]] for atom_id in base.metadata.atom_ids}
        second_coordinates = {second_ids[old]: old_coordinates[old] + np.asarray((0.0, 15.0, 0.0)) for old in selected_second}
        coordinates = np.asarray([old_coordinates[atom_id] if atom_id in old_coordinates else second_coordinates[atom_id] for atom_id in atom_ids])
        frames.append(TrajectoryFrame(coordinates, frame.box_vectors_angstrom, frame.step, frame.time_ps, frame.temperature_kelvin, frame.pressure_bar, frame.potential_energy_kj_mol, frame.volume_angstrom3, frame.target_temperature_kelvin, frame.target_pressure_bar, frame.measured_temperature_kelvin, frame.measured_pressure_bar))
    metadata = TrajectoryMetadata(
        final_graph_identity=_digest({"base": base.metadata.final_graph_identity, "membership": membership, "bonds": bonds}),
        packing_identity=_digest({"base": base.metadata.packing_identity, "chains": ["chain-0", "chain-1"]}),
        force_field_family=base.metadata.force_field_family,
        force_field_source_identity=base.metadata.force_field_source_identity,
        typed_graph_identity=_digest({"base": base.metadata.typed_graph_identity, "atom_ids": list(atom_ids)}),
        charge_identity=_digest({"base": base.metadata.charge_identity, "atom_ids": list(atom_ids)}),
        workflow_identity=_digest({"base": base.metadata.workflow_identity, "fixture": "phase-4j30"}),
        checkpoint_identity=_digest({"base": base.metadata.checkpoint_identity, "fixture": "phase-4j30"}),
        assignment_identity=_digest({"base": base.metadata.assignment_identity, "fixture": "phase-4j30"}),
        atom_ids=atom_ids, masses_dalton=masses, atom_mapping=mapping,
        box_periodic=base.metadata.box_periodic, temperature_kelvin=base.metadata.temperature_kelvin,
        pressure_bar=base.metadata.pressure_bar, timestep_fs=base.metadata.timestep_fs,
        ensemble=base.metadata.ensemble, source_files=base.metadata.source_files,
        source_metadata={**base.metadata.source_metadata, "derived_multichain_fixture": True},
    )
    trajectory = VerifiedTrajectory(metadata, tuple(frames), f"derived:{base.source_path}", {"base_workflow": _digest(base.trajectory_file_hashes)}, tuple(bonds), (), membership)
    selection = AnalysisSelection(
        final_graph_identity=metadata.final_graph_identity,
        atom_mapping_identity=trajectory.mapping_identity,
        molecule_membership=membership,
        chains=(
            {"chain_id": "chain-0", "molecule_id": "chain-0", "rg_atom_ids": list(chain1_ids), "backbone_atom_ids": list(range(1, 7)), "backbone_path": list(range(1, 7)), "endpoint_ids": [1, 6]},
            {"chain_id": "chain-1", "molecule_id": "chain-1", "rg_atom_ids": list(second_ids.values()), "backbone_atom_ids": [second_ids[index] for index in (1, 2, 3, 4)], "backbone_path": [second_ids[index] for index in (1, 2, 3, 4)], "endpoint_ids": [second_ids[1], second_ids[4]]},
        ),
        provenance={"source": "retained_j23_workflow", "packing_contract": "J21-derived bounded fixture"},
    )
    return VerifiedTrajectory(metadata, tuple(frames), trajectory.source_path, trajectory.trajectory_file_hashes, tuple(bonds), tuple((chain["endpoint_ids"]) for chain in selection.chains), membership, selection)
