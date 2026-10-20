"""J31 topology, geometry, lineage, persistence and negative controls."""

import json
import os
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from island.analysis import (
    calculate_network_observables,
    load_network_analysis_bundle,
    load_network_analysis_record,
    make_network_analysis_record,
    relocate_network_analysis_record,
    save_network_analysis_bundle,
    unwrap_network,
    write_network_analysis_record,
)
from island.analysis.network import _bytes, _image_geometry, _signed
from island.core import AtomSite
from island.crosslinking import apply_crosslink_plan, identify_reactive_sites
from island.exceptions import AnalysisError, IslandError
from island.graph import final_graph
from scripts.phase_4j31_network_fixture import (
    RULE,
    build_network_fixture,
    graph_control_trajectory,
    independent_topology,
    select_network,
)


@pytest.fixture(scope="module", params=("uncrosslinked", "branched", "crosslinked"))
def fixture(request):
    return build_network_fixture(request.param)


@pytest.fixture(scope="module")
def cyclic():
    return build_network_fixture("crosslinked")


def test_real_j20_j21_lineage_and_independent_graph_reconstruction(fixture):
    trajectory, selection = graph_control_trajectory(fixture)
    values = calculate_network_observables(trajectory, selection)
    reference = independent_topology(fixture.packed.graph)
    assert len(set(trajectory.molecule_membership.values())) == 2
    assert (
        values["network_connectivity"]["value"]["component_count"]
        == reference["component_count"]
    )
    assert values["degree_distribution"]["value"] == reference["degree_counts"]
    assert values["cycle_rank"]["value"] == reference["cycle_rank"]
    assert (
        values["crosslink_density"]["crosslink_count"] == reference["crosslink_count"]
    )
    assert values["crosslink_density"]["value"] == pytest.approx(
        [reference["crosslink_density"]] * trajectory.frame_count
    )
    assert values["strand_length_distribution"]["value"] == reference["strand_lengths"]
    assert values["gel_fraction"]["value"] == max(reference["component_sizes"]) / len(
        trajectory.metadata.atom_ids
    )
    assert values["periodic_geometry"]["status"] == "ok"
    fixture.crosslinked.plan.validate_integrity(fixture.input_graph)
    fixture.packed.plan.validate_integrity(
        [("two-chain-network", fixture.crosslinked.system, fixture.crosslinked.graph)]
    )
    labels = {
        site["id"]: site["element"] for site in fixture.input_graph.payload["sites"]
    }
    if fixture.kind != "uncrosslinked":
        assert any(
            {labels[a], labels[b]} == {"C", "O"}
            for a, b in fixture.crosslinked.plan.payload["selected_crosslink_bonds"]
        )
        assert values["branch_point_statistics"]["value"]["branch_point_count"] > 0
    assert values["gel_fraction"]["policy"].startswith("largest_atom_count")


def test_determinism_and_read_only_inputs(fixture):
    before = deepcopy(fixture.packed.system.to_dict())
    other = build_network_fixture(fixture.kind)
    assert other.crosslinked.plan.identity == fixture.crosslinked.plan.identity
    assert other.packed.graph.identity == fixture.packed.graph.identity
    trajectory, selection = graph_control_trajectory(fixture)
    first = calculate_network_observables(trajectory, selection)
    assert first == calculate_network_observables(trajectory, selection)
    assert before == fixture.packed.system.to_dict()
    assert (
        selection.payload["molecule_membership"]
        != selection.payload["component_membership"]
    )


def _least_squares_rg(trajectory, graph):
    """Independent incidence-matrix solve, not a spanning-tree traversal."""
    p = graph.payload
    component = min(p["components"], key=lambda c: (-len(c), min(c)))
    indices = {sid: n for n, sid in enumerate(component)}
    frame = trajectory.frames[0]
    lengths = np.diag(frame.box_vectors_angstrom)
    rows, displacements = [], []
    for bond in p["bonds"]:
        a, b = bond["sites"]
        if a not in indices:
            continue
        row = np.zeros(len(component))
        row[indices[b]], row[indices[a]] = 1, -1
        rows.append(row)
        delta = (
            frame.coordinates_angstrom[trajectory.metadata.atom_mapping[b]]
            - frame.coordinates_angstrom[trajectory.metadata.atom_mapping[a]]
        )
        displacements.append(delta - np.floor(delta / lengths + 0.5) * lengths)
    points, *_ = np.linalg.lstsq(
        np.asarray(rows), np.asarray(displacements), rcond=None
    )
    assert np.max(np.abs(np.asarray(rows) @ points - displacements)) < 1e-10
    masses = [
        trajectory.metadata.masses_dalton[trajectory.metadata.atom_mapping[i]]
        for i in component
    ]
    center = np.average(points, weights=masses, axis=0)
    return float(
        np.sqrt(np.average(np.sum((points - center) ** 2, axis=1), weights=masses))
    )


def test_network_rg_independent_incidence_reference(fixture):
    trajectory, selection = graph_control_trajectory(fixture)
    values = calculate_network_observables(trajectory, selection)
    assert values["network_radius_of_gyration"]["value"] == pytest.approx(
        [_least_squares_rg(trajectory, fixture.packed.graph)] * trajectory.frame_count
    )
    assert all(
        row["value"] == [False, False, False]
        for row in values["periodic_spanning"]["per_frame"]
    )


def test_integer_box_translation_invariance(fixture):
    trajectory, selection = graph_control_trajectory(fixture)
    baseline = calculate_network_observables(trajectory, selection)
    frames = []
    for frame in trajectory.frames:
        coordinates = frame.coordinates_angstrom.copy()
        for sid, molecule in trajectory.molecule_membership.items():
            if molecule.endswith("chain-1"):
                coordinates[trajectory.metadata.atom_mapping[sid]] += [
                    80.0,
                    -40.0,
                    120.0,
                ]
        frames.append(replace(frame, coordinates_angstrom=coordinates))
    translated = replace(trajectory, frames=tuple(frames))
    with pytest.raises(AnalysisError, match="identity mismatch"):
        selection.validate(translated)
    moved = select_network(translated, fixture)
    values = calculate_network_observables(translated, moved)
    assert values["network_radius_of_gyration"]["value"] == pytest.approx(
        baseline["network_radius_of_gyration"]["value"]
    )
    assert values["contact_fraction"]["value"] == baseline["contact_fraction"]["value"]
    assert values["periodic_spanning"] == baseline["periodic_spanning"]
    for name in (
        "degree_distribution",
        "strand_length_distribution",
        "crosslink_density",
        "cycle_rank",
    ):
        assert values[name] == baseline[name]


@pytest.mark.parametrize(
    "name",
    (
        "final_graph_identity",
        "packing_identity",
        "force_field_source_identity",
        "typed_graph_identity",
        "charge_identity",
        "assignment_identity",
        "workflow_identity",
        "checkpoint_identity",
        "force_field_family",
    ),
)
def test_stale_metadata_rejected(cyclic, name):
    trajectory, selection = graph_control_trajectory(cyclic)
    changed = replace(
        trajectory, metadata=replace(trajectory.metadata, **{name: "f" * 64})
    )
    with pytest.raises(AnalysisError):
        selection.validate(changed)


def test_changed_plan_mapping_bonds_and_frame_content_rejected(cyclic):
    trajectory, selection = graph_control_trajectory(cyclic)
    meta = replace(
        trajectory.metadata,
        source_metadata={
            **trajectory.metadata.source_metadata,
            "crosslink_plan_identity": "f" * 64,
        },
    )
    with pytest.raises(AnalysisError, match="plan"):
        selection.validate(replace(trajectory, metadata=meta))
    mapping = dict(trajectory.metadata.atom_mapping)
    a, b = tuple(mapping)[:2]
    mapping[a], mapping[b] = mapping[b], mapping[a]
    with pytest.raises(AnalysisError):
        selection.validate(
            replace(
                trajectory, metadata=replace(trajectory.metadata, atom_mapping=mapping)
            )
        )
    with pytest.raises(AnalysisError, match="bond"):
        selection.validate(
            replace(trajectory, topology_bonds=trajectory.topology_bonds[:-1])
        )
    frame = trajectory.frames[-1]
    coords = frame.coordinates_angstrom.copy()
    coords[0] += [0.01, 0, 0]
    changed = replace(
        trajectory,
        frames=trajectory.frames[:-1] + (replace(frame, coordinates_angstrom=coords),),
    )
    with pytest.raises(AnalysisError, match="identity mismatch"):
        selection.validate(changed)


def test_invalid_valence_duplicate_bond_stale_plan(cyclic):
    system = cyclic.input_system.copy()
    target = next(
        i
        for i, a in system.topology.sites.items()
        if a.element == "C" and not a.metadata.get("j31_reactive")
    )
    system.topology.sites[target].metadata["j31_reactive"] = True
    with pytest.raises(IslandError, match="capacity"):
        identify_reactive_sites(final_graph(system), RULE)
    system.topology.add_site(AtomSite(9999, "H", 1.008, element="H", atomic_number=1))
    system.coordinates.set(9999, [0, 0, 0])
    system.topology.add_bond(target, 9999, order=1)
    with pytest.raises(IslandError, match="valence"):
        final_graph(system)
    a, b = next(iter(cyclic.input_system.topology.bonds))
    with pytest.raises(IslandError):
        cyclic.input_system.copy().topology.add_bond(a, b)
    changed = cyclic.input_system.copy()
    changed.topology.sites[target].metadata["j31_changed"] = True
    with pytest.raises(IslandError):
        apply_crosslink_plan(changed, final_graph(changed), cyclic.crosslinked.plan)
    with pytest.raises(IslandError):
        apply_crosslink_plan(
            cyclic.crosslinked.system, cyclic.crosslinked.graph, cyclic.crosslinked.plan
        )


def test_non_tree_cycle_is_not_hidden_by_modulo():
    from island.analysis import TrajectoryFrame

    frame = TrajectoryFrame(
        np.array([[0.0, 0, 0], [4.0, 0, 0], [8.0, 0, 0]]),
        np.diag([10.0, 10.0, 10.0]),
        0,
        0.0,
    )
    geometry = _image_geometry(
        frame, {1: 0, 2: 1, 3: 2}, {1: {2, 3}, 2: {1, 3}, 3: {1, 2}}, 1e-6
    )
    assert not geometry["consistent"]
    assert geometry["winding_cycles"] == [{"bond": [2, 3], "winding": [-1, 0, 0]}]


@pytest.mark.parametrize("distance", (5.0, -5.0, 15.0, 5.0 + 1e-8))
def test_half_box_ambiguity_rejected(distance):
    from island.analysis import TrajectoryFrame

    frame = TrajectoryFrame(
        np.array([[0.0, 0, 0], [distance, 0, 0]]), np.diag([10.0, 10.0, 10.0]), 0, 0.0
    )
    with pytest.raises(AnalysisError, match="Ambiguous"):
        _image_geometry(frame, {1: 0, 2: 1}, {1: {2}, 2: {1}}, 1e-6)


def test_second_frame_cycle_blocks_rg_but_reports_percolation(cyclic):
    trajectory, _ = graph_control_trajectory(cyclic)
    graph = cyclic.packed.graph.payload
    heavy = {s["id"] for s in graph["sites"] if s["element"] != "H"}
    adjacency = {i: set() for i in heavy}
    for bond in graph["bonds"]:
        a, b = bond["sites"]
        if a in heavy and b in heavy:
            adjacency[a].add(b)
            adjacency[b].add(a)
    # Independently peel leaves; the remaining induced graph is the sole ring.
    remaining = set(heavy)
    while True:
        leaves = {i for i in remaining if len(adjacency[i] & remaining) < 2}
        if not leaves:
            break
        remaining -= leaves
    start = min(remaining)
    cycle, previous, current = [start], None, start
    while True:
        nxt = next(i for i in sorted(adjacency[current] & remaining) if i != previous)
        if nxt == start:
            break
        cycle.append(nxt)
        previous, current = current, nxt
    xyz = trajectory.frames[-1].coordinates_angstrom.copy()
    for n, sid in enumerate(cycle):
        xyz[trajectory.metadata.atom_mapping[sid], 0] = 40.0 * n / len(cycle)
    changed = replace(
        trajectory,
        frames=(
            trajectory.frames[0],
            replace(trajectory.frames[-1], coordinates_angstrom=xyz),
        ),
    )
    selection = select_network(changed, cyclic)
    values = calculate_network_observables(changed, selection)
    assert values["network_radius_of_gyration"]["per_frame"][0]["status"] == "ok"
    assert values["network_radius_of_gyration"]["per_frame"][1]["status"] == "blocked"
    assert values["periodic_spanning"]["per_frame"][1]["value"][0]
    with pytest.raises(AnalysisError, match="cycle closure"):
        unwrap_network(changed, selection, frame_index=1)


@pytest.mark.parametrize(
    "defect", ("triclinic", "nonperiodic", "contact_cutoff", "missing_assignment")
)
def test_unsupported_geometry_and_assignment_rejected(cyclic, defect):
    trajectory, _ = graph_control_trajectory(cyclic)
    if defect == "triclinic":
        box = trajectory.frames[0].box_vectors_angstrom.copy()
        box[0, 1] = 0.1
        trajectory = replace(
            trajectory,
            frames=(replace(trajectory.frames[0], box_vectors_angstrom=box),),
        )
    elif defect == "nonperiodic":
        trajectory = replace(
            trajectory,
            metadata=replace(trajectory.metadata, box_periodic=(True, True, False)),
        )
    elif defect == "missing_assignment":
        trajectory = replace(
            trajectory, metadata=replace(trajectory.metadata, assignment_identity=None)
        )
    with pytest.raises(AnalysisError):
        select_network(
            trajectory,
            cyclic,
            contact_cutoff_angstrom=25.0 if defect == "contact_cutoff" else 4.0,
        )


def test_contacts_declared_denominator_and_unsupported_observables(cyclic):
    trajectory, selection = graph_control_trajectory(cyclic)
    values = calculate_network_observables(trajectory, selection)
    row = values["contact_fraction"]["per_frame"][0]
    assert 0 <= row["contacts"] <= row["denominator"]
    assert row["value"] == row["contacts"] / row["denominator"]
    assert (
        values["dangling_end_fraction"]["value"]
        == len(values["dangling_end_fraction"]["dangling_end_ids"]) / 12
    )
    no_cutoff = select_network(trajectory, cyclic, contact_cutoff_angstrom=None)
    assert (
        calculate_network_observables(trajectory, no_cutoff)["contact_fraction"][
            "status"
        ]
        == "blocked"
    )
    assert (
        calculate_network_observables(trajectory, selection, observables=("Tg",))["Tg"][
            "status"
        ]
        == "blocked"
    )


def test_record_relocation_fresh_process_and_owned_payload(cyclic, tmp_path):
    trajectory, selection = graph_control_trajectory(cyclic)
    record = make_network_analysis_record(trajectory, selection)
    path = write_network_analysis_record(record, tmp_path / "record.json")
    moved = relocate_network_analysis_record(path, tmp_path / "moved.json")
    loaded = load_network_analysis_record(
        moved, expected_identity=record.identity, trajectory=trajectory
    )
    assert loaded.payload == record.payload
    payload = loaded.payload
    payload["selection"]["assignment_identity"] = "f" * 64
    assert loaded.payload == record.payload
    code = "from island.analysis import load_network_analysis_record; import sys; print(load_network_analysis_record(sys.argv[1],expected_identity=sys.argv[2]).to_json(),end='')"
    result = subprocess.run(
        [sys.executable, "-c", code, str(moved), record.identity],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout == record.to_json()
    with pytest.raises(IslandError):
        write_network_analysis_record(record, path)
    with pytest.raises(AnalysisError):
        load_network_analysis_record(moved, expected_selection_identity="f" * 64)


def test_transactional_bundle_relocation_and_fresh_reconstruction(cyclic, tmp_path):
    trajectory, selection = graph_control_trajectory(cyclic)
    path = save_network_analysis_bundle(trajectory, selection, tmp_path / "original")
    moved = tmp_path / "moved"
    path.rename(moved)
    restored, chosen, record = load_network_analysis_bundle(moved)
    assert chosen.identity == selection.identity
    assert (
        calculate_network_observables(restored, chosen) == record.payload["observables"]
    )
    code = "from island.analysis import load_network_analysis_bundle; import sys; print(load_network_analysis_bundle(sys.argv[1],expected_identity=sys.argv[2])[2].to_json(),end='')"
    result = subprocess.run(
        [sys.executable, "-c", code, str(moved), record.identity],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout == record.to_json()
    with pytest.raises(AnalysisError, match="exists"):
        save_network_analysis_bundle(trajectory, selection, moved)
    assert not list(tmp_path.glob(".network-stage-*"))


def test_transactional_failure_cleans_stage(cyclic, tmp_path, monkeypatch):
    import island.periodic_md

    trajectory, selection = graph_control_trajectory(cyclic)

    def fail(*args):
        raise OSError("injected publication failure")

    monkeypatch.setattr(island.periodic_md, "_publish_directory", fail)
    with pytest.raises(OSError, match="publication"):
        save_network_analysis_bundle(trajectory, selection, tmp_path / "unpublished")
    assert not (tmp_path / "unpublished").exists()
    assert not list(tmp_path.glob(".network-stage-*"))


@pytest.mark.parametrize(
    "defect", ("frame", "selection", "observables", "readiness", "schema")
)
def test_resigned_bundle_mutations_rejected(cyclic, tmp_path, defect):
    trajectory, selection = graph_control_trajectory(cyclic)
    path = (
        save_network_analysis_bundle(trajectory, selection, tmp_path / "bundle")
        / "bundle.json"
    )
    p = json.loads(path.read_text())
    if defect == "frame":
        p["frames"][-1]["coordinates"][0][0] += 0.01
    elif defect == "selection":
        sel = p["record"]["selection"]
        sel["assignment_identity"] = "f" * 64
        p["record"]["selection"] = _signed(
            {k: v for k, v in sel.items() if k != "identity"}
        )
    elif defect == "observables":
        p["record"]["observables"]["cycle_rank"]["value"] = 999
    elif defect == "readiness":
        p["record"]["production_validated"] = True
    else:
        p["record"]["schema"] = "unknown"
    p["record"] = _signed({k: v for k, v in p["record"].items() if k != "identity"})
    path.write_bytes(_bytes(_signed({k: v for k, v in p.items() if k != "identity"})))
    with pytest.raises(AnalysisError):
        load_network_analysis_bundle(path.parent)


def test_historical_receipts_unchanged():
    import hashlib
    import subprocess

    base = "aa3014648ddda4b032d6b07cdea381792330f99d"
    for path in Path("docs/evidence").glob("phase_4j*.json"):
        if path.name == "phase_4j31.json":
            continue
        original = subprocess.check_output(["git", "show", f"{base}:{path}"])
        assert (
            hashlib.sha256(path.read_bytes()).digest()
            == hashlib.sha256(original).digest()
        )


def test_no_forcefield_or_chemistry_import_for_network_module():
    code = "import island.analysis.network,sys; assert not set(('rdkit','openmm','foyer','parmed')) & sys.modules.keys()"
    subprocess.run([sys.executable, "-c", code], check=True, env={**os.environ})
