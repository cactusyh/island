"""Bounded graph-bound periodic network geometry and transactional records.

No force-field terms or trajectory labels are inferred here. The caller supplies
an authoritative final graph and actual assignment identities (or explicitly
blocked graph-only evidence). Geometry is evaluated independently in every frame.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import tempfile
from collections import deque
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np

from island.exceptions import AnalysisError
from island.graph import FinalChemicalGraph
from island.graph.final import pack
from island.graph.network import network_topology

from .mdanalysis import TrajectoryFrame, TrajectoryMetadata, VerifiedTrajectory

NETWORK_SELECTION_SCHEMA = "island.network-analysis-selection.v1"
NETWORK_RECORD_SCHEMA = "island.network-analysis-record.v1"
BUNDLE_SCHEMA = "island.network-analysis-bundle.v1"
FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
POLICY = "heavy_skeleton_strands_all_atom_clusters_v1"
GEOMETRY = "orthorhombic_minimum_image_tree_cycle_v1"
IDENTITIES = (
    "final_graph_identity",
    "crosslink_plan_identity",
    "atom_mapping_identity",
    "packing_identity",
    "force_field_source_identity",
    "typed_graph_identity",
    "charge_identity",
    "assignment_identity",
    "workflow_identity",
    "checkpoint_identity",
    "trajectory_identity",
    "frame_content_identity",
)


def _bytes(value):
    return (
        json.dumps(
            json.loads(json.dumps(value, allow_nan=False)),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode()


def _digest(value):
    return sha256(_bytes(value)).hexdigest()


def _require(ok, message):
    if not ok:
        raise AnalysisError(message)


def _identity(value):
    _require(
        type(value) is str
        and len(value) == 64
        and all(c in "0123456789abcdef" for c in value),
        "Malformed network identity",
    )


def _frame_data(trajectory):
    return [
        {
            "coordinates": f.coordinates_angstrom.tolist(),
            "box": f.box_vectors_angstrom.tolist(),
            "step": f.step,
            "time_ps": f.time_ps,
        }
        for f in trajectory.frames
    ]


def _graph(raw):
    return FinalChemicalGraph(pack(raw))


def _body(raw):
    return {k: v for k, v in raw.items() if k != "identity"}


def _signed(value):
    return {**value, "identity": _digest(value)}


def _validate_lineage(graph, plan, packing_plan, lineage):
    from island.crosslinking import apply_crosslink_plan
    from island.graph.final import _system_from_payload
    from island.packing import apply_periodic_packing_plan

    source_graph = _graph(lineage["graph"])
    source_system = _system_from_payload(lineage["system"])
    source_graph.validate_integrity(source_system)
    plan.validate_integrity(source_graph)
    system, crosslinked, _ = apply_crosslink_plan(source_system, source_graph, plan)
    from island.core import Coordinates

    system.coordinates = Coordinates(
        {int(i): xyz for i, xyz in lineage["post_crosslink_coordinates"].items()}
    )
    crosslinked.validate_integrity(system)
    labels = packing_plan.payload["unit_labels"]
    _require(
        len(labels) == 1,
        "Bounded J31 lineage requires one prelinked network packing unit",
    )
    packed = apply_periodic_packing_plan(
        [(labels[0], system, crosslinked)], packing_plan
    )
    _require(
        packed.graph.payload == graph.payload,
        "Final graph contradicts replayed J20/J21 lineage",
    )


def _selection_payload(
    trajectory,
    graph,
    plan,
    packing_plan,
    *,
    lineage,
    contact_cutoff_angstrom,
    tolerance_angstrom,
    provenance,
):
    _require(type(trajectory) is VerifiedTrajectory, "Expected VerifiedTrajectory")
    _validate_lineage(graph, plan, packing_plan, lineage)
    facts = network_topology(graph)
    p = graph.payload
    _require(
        graph.identity == trajectory.metadata.final_graph_identity,
        "Final graph identity mismatch",
    )
    _require(
        set(trajectory.metadata.atom_ids) == {s["id"] for s in p["sites"]},
        "Final graph atom IDs mismatch",
    )
    bonds = sorted(tuple(sorted(b)) for b in trajectory.topology_bonds)
    _require(
        bonds == sorted(tuple(b["sites"]) for b in p["bonds"]),
        "Final graph bond mismatch",
    )
    membership = {int(i): v for i, v in p["molecule_membership"].items()}
    _require(
        membership == trajectory.molecule_membership, "Molecule membership mismatch"
    )
    _require(
        len(set(membership.values())) >= 2, "At least two molecular chains required"
    )
    for s in p["sites"]:
        _require(
            s["mass"]
            == trajectory.metadata.masses_dalton[
                trajectory.metadata.atom_mapping[s["id"]]
            ],
            "Final graph mass mismatch",
        )
    _require(
        trajectory.metadata.box_periodic == (True, True, True),
        "Fully periodic trajectory required",
    )
    _require(
        type(tolerance_angstrom) in (int, float)
        and math.isfinite(tolerance_angstrom)
        and tolerance_angstrom > 0,
        "Positive geometry tolerance required",
    )
    if contact_cutoff_angstrom is not None:
        _require(
            type(contact_cutoff_angstrom) in (int, float)
            and math.isfinite(contact_cutoff_angstrom)
            and contact_cutoff_angstrom > 0,
            "Positive contact cutoff required",
        )
    for frame in trajectory.frames:
        box = np.diag(frame.box_vectors_angstrom)
        _require(
            np.all(box > 0)
            and np.allclose(
                frame.box_vectors_angstrom, np.diag(box), rtol=0, atol=1e-10
            ),
            "Orthorhombic positive box required",
        )
        _require(
            tolerance_angstrom < min(box) / 4, "Geometry tolerance exceeds box budget"
        )
        if contact_cutoff_angstrom is not None:
            _require(
                2 * contact_cutoff_angstrom < min(box),
                "Contact cutoff must be below half the shortest box",
            )
    _identity(trajectory.metadata.assignment_identity)
    plan.validate_integrity()
    crosslinks = p["crosslink_provenance"]
    _require(
        all(row["plan_identity"] == plan.identity for row in crosslinks),
        "Crosslink plan provenance mismatch",
    )
    _require(
        len(crosslinks) == len(plan.payload["selected_crosslink_bonds"]),
        "Crosslink count disagrees with plan",
    )
    _require(
        trajectory.metadata.source_metadata.get("crosslink_plan_identity")
        == plan.identity,
        "Trajectory crosslink plan mismatch",
    )
    packing = [
        tr
        for tr in p["transformations"]
        if tr["operation"] == "periodic_multichain_packing"
    ]
    _require(
        len(packing) == 1
        and packing[0]["parameters"]["plan_identity"]
        == trajectory.metadata.packing_identity,
        "Packing identity mismatch",
    )
    # J21 explicitly records the remapping. Verify each final crosslink against
    # the J20 plan; no bonds may be created by analysis or by a string label.
    packing_plan.validate_integrity()
    _require(
        packing_plan.identity == trajectory.metadata.packing_identity,
        "Packing plan witness mismatch",
    )
    remaps = packing_plan.payload["site_id_remapping"]
    expected = []
    for mapping in remaps:
        if all(
            str(i) in mapping
            for b in plan.payload["selected_crosslink_bonds"]
            for i in b
        ):
            expected.extend(
                sorted((mapping[str(a)], mapping[str(b)]))
                for a, b in plan.payload["selected_crosslink_bonds"]
            )
    _require(
        sorted(expected) == sorted(row["sites"] for row in crosslinks),
        "Packed crosslinks disagree with J20 plan",
    )
    metadata = trajectory.metadata
    return {
        "schema": NETWORK_SELECTION_SCHEMA,
        "final_graph": p,
        "lineage": lineage,
        "crosslink_plan": plan.payload,
        "packing_plan": packing_plan.payload,
        "final_graph_identity": graph.identity,
        "crosslink_plan_identity": plan.identity,
        "molecule_membership": [[i, membership[i]] for i in sorted(membership)],
        "component_membership": facts["component_membership"],
        "atom_ids": list(metadata.atom_ids),
        "atom_mapping_identity": trajectory.mapping_identity,
        "packing_identity": metadata.packing_identity,
        "force_field_family": metadata.force_field_family,
        "force_field_source_identity": metadata.force_field_source_identity,
        "typed_graph_identity": metadata.typed_graph_identity,
        "charge_identity": metadata.charge_identity,
        "assignment_identity": metadata.assignment_identity,
        "workflow_identity": metadata.workflow_identity,
        "checkpoint_identity": metadata.checkpoint_identity,
        "trajectory_identity": trajectory.identity,
        "frame_content_identity": _digest(_frame_data(trajectory)),
        "periodic_geometry_policy": GEOMETRY,
        "network_policy": POLICY,
        "tolerance_angstrom": tolerance_angstrom,
        "contact_cutoff_angstrom": contact_cutoff_angstrom,
        "gel_policy": "largest_atom_count_component_minimum_stable_id_tiebreak",
        "provenance": dict(provenance or {}),
    }


@dataclass(frozen=True)
class NetworkAnalysisSelection:
    """Owned canonical JSON, including the authoritative graph and J20 plan."""

    json_text: str

    def __post_init__(self):
        p = json.loads(self.json_text)
        _require(
            p.get("schema") == NETWORK_SELECTION_SCHEMA
            and p.get("identity") == _digest(_body(p)),
            "Network selection checksum/schema mismatch",
        )
        for name in IDENTITIES:
            _identity(p.get(name))
        _require(
            p.get("network_policy") == POLICY
            and p.get("periodic_geometry_policy") == GEOMETRY,
            "Unsupported network policy",
        )
        facts = network_topology(_graph(p["final_graph"]))
        _require(
            _graph(p["final_graph"]).identity == p["final_graph_identity"]
            and facts["component_membership"] == p["component_membership"],
            "Network graph/component binding mismatch",
        )
        from island.crosslinking import CrosslinkNetworkPlan

        plan = CrosslinkNetworkPlan(pack(p["crosslink_plan"]))
        from island.packing import PeriodicPackingPlan

        plan.validate_integrity()
        _validate_lineage(
            _graph(p["final_graph"]),
            plan,
            PeriodicPackingPlan(pack(p["packing_plan"])),
            p["lineage"],
        )
        _require(
            plan.identity == p["crosslink_plan_identity"],
            "Network plan binding mismatch",
        )

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    @classmethod
    def from_trajectory(
        cls,
        trajectory,
        *,
        graph,
        crosslink_plan,
        packing_plan,
        input_system,
        input_graph,
        crosslinked_coordinates,
        contact_cutoff_angstrom=None,
        tolerance_angstrom=1e-6,
        provenance=None,
    ):
        p = _selection_payload(
            trajectory,
            graph,
            crosslink_plan,
            packing_plan,
            lineage={
                "system": input_system.to_dict(),
                "graph": input_graph.payload,
                "post_crosslink_coordinates": crosslinked_coordinates,
            },
            contact_cutoff_angstrom=contact_cutoff_angstrom,
            tolerance_angstrom=tolerance_angstrom,
            provenance=provenance,
        )
        return cls(_bytes(_signed(p)).decode())

    def to_dict(self):
        return self.payload

    def validate(self, trajectory):
        from island.crosslinking import CrosslinkNetworkPlan
        from island.packing import PeriodicPackingPlan

        p = self.payload
        expected = _selection_payload(
            trajectory,
            _graph(p["final_graph"]),
            CrosslinkNetworkPlan(pack(p["crosslink_plan"])),
            PeriodicPackingPlan(pack(p["packing_plan"])),
            lineage=p["lineage"],
            contact_cutoff_angstrom=p["contact_cutoff_angstrom"],
            tolerance_angstrom=p["tolerance_angstrom"],
            provenance=p["provenance"],
        )
        _require(
            _signed(expected) == p,
            "Network selection identity mismatch (stale graph, plan, mapping, assignment, workflow, checkpoint or trajectory)",
        )


def _image_geometry(frame, mapping, adjacency, tolerance):
    """Sorted BFS roots/neighbours; verify every non-tree edge without modulo.

    Nonzero integer closure vectors are winding cycles, not a valid finite
    embedding. They establish graph percolation axes but block finite-network Rg.
    Half-box ties have two equally near images and are rejected explicitly.
    """
    lengths = np.diag(frame.box_vectors_angstrom)
    xyz = frame.coordinates_angstrom
    images, points, windings, trees = {}, {}, [], []
    for root in sorted(adjacency):
        if root in images:
            continue
        images[root] = np.zeros(3, dtype=int)
        points[root] = xyz[mapping[root]].copy()
        queue = deque([root])
        while queue:
            a = queue.popleft()
            for b in sorted(adjacency[a]):
                delta = xyz[mapping[b]] - xyz[mapping[a]]
                frac = delta / lengths
                _require(
                    not np.any(
                        np.abs(np.abs(frac - np.rint(frac)) - 0.5) * lengths
                        <= tolerance
                    ),
                    "Ambiguous half-box periodic image",
                )
                edge_image = -np.rint(frac).astype(int)
                candidate = images[a] + edge_image
                if b not in images:
                    images[b] = candidate
                    points[b] = xyz[mapping[b]] + candidate * lengths
                    trees.append([a, b])
                    queue.append(b)
                elif a < b:
                    winding = images[b] - candidate
                    residual = points[b] - points[a] - (delta + edge_image * lengths)
                    _require(
                        np.allclose(
                            residual, winding * lengths, atol=tolerance, rtol=0
                        ),
                        "Inconsistent periodic image arithmetic",
                    )
                    if np.any(winding):
                        windings.append({"bond": [a, b], "winding": winding.tolist()})
    return {
        "positions": points,
        "images": {str(i): images[i].tolist() for i in sorted(images)},
        "tree_edges": trees,
        "winding_cycles": windings,
        "consistent": not windings,
    }


def unwrap_network(trajectory, selection, *, frame_index=0):
    """Strict finite embedding; reject winding cycles and ambiguous images."""
    selection.validate(trajectory)
    p = selection.payload
    adjacency = {i: set() for i in p["atom_ids"]}
    for bond in p["final_graph"]["bonds"]:
        a, b = bond["sites"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    result = _image_geometry(
        trajectory.frames[frame_index],
        trajectory.metadata.atom_mapping,
        adjacency,
        p["tolerance_angstrom"],
    )
    _require(
        result["consistent"], "Inconsistent periodic cycle closure: nonzero winding"
    )
    return result


def _ok(value, units, **extra):
    return {"status": "ok", "value": value, "units": units, **extra}


def _blocked(reason):
    return {"status": "blocked", "value": None, "failure": reason}


def calculate_network_observables(trajectory, selection, *, observables=None):
    selection.validate(trajectory)
    p = selection.payload
    facts = network_topology(_graph(p["final_graph"]))
    clusters = facts["components"]
    largest = min(clusters, key=lambda c: (-len(c), min(c)))
    adjacency = {i: set() for i in p["atom_ids"]}
    for bond in p["final_graph"]["bonds"]:
        a, b = bond["sites"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    n = len(adjacency)
    heavy = [s["id"] for s in p["final_graph"]["sites"] if s["element"] != "H"]
    result = {
        "network_connectivity": _ok(
            {
                "component_count": len(clusters),
                "components": clusters,
                "gel_component_ids": largest,
            },
            "count/stable_atom_id",
        ),
        "degree_distribution": _ok(
            {str(k): v for k, v in facts["degree_distribution"].items()},
            "atom_count",
            policy="all_atoms",
        ),
        "branch_point_statistics": _ok(
            {
                "branch_point_ids": facts["branch_points"],
                "branch_point_count": len(facts["branch_points"]),
            },
            "count",
            policy="heavy_atom_degree_gt_two",
        ),
        "cluster_size_distribution": _ok(
            {str(k): v for k, v in facts["cluster_size_distribution"].items()},
            "cluster_count_by_atom_count",
        ),
        "gel_fraction": _ok(len(largest) / n, "atom_fraction", policy=p["gel_policy"]),
        "cycle_rank": _ok(facts["cycle_rank"], "count"),
        "strand_length_distribution": _ok(
            facts["strand_lengths_bonds"],
            "heavy_bond_count",
            policy="maximal_edge_disjoint_paths_crosslink_endpoints_are_junctions",
        ),
        "dangling_end_fraction": _ok(
            len(facts["dangling_end_ids"]) / len(heavy) if heavy else 0.0,
            "heavy_atom_fraction",
            dangling_end_ids=facts["dangling_end_ids"],
        ),
        "crosslink_density": _ok(
            [len(facts["crosslink_bonds"]) / f.volume for f in trajectory.frames],
            "crosslink_count/angstrom^3",
            crosslink_count=len(facts["crosslink_bonds"]),
            denominator="per_frame_periodic_cell_volume_angstrom3",
            volumes=[f.volume for f in trajectory.frames],
        ),
    }
    rg, percolation, audits, contacts = [], [], [], []
    eligible = [
        (a, b) for a, b in combinations(sorted(heavy), 2) if b not in adjacency[a]
    ]
    for frame in trajectory.frames:
        try:
            geometry = _image_geometry(
                frame,
                trajectory.metadata.atom_mapping,
                adjacency,
                p["tolerance_angstrom"],
            )
            axes = [
                any(row["winding"][axis] != 0 for row in geometry["winding_cycles"])
                for axis in range(3)
            ]
            percolation.append(
                _ok(
                    axes,
                    "boolean_by_xyz",
                    policy="nonzero_cycle_lattice_winding",
                    winding_cycles=geometry["winding_cycles"],
                )
            )
            audits.append({k: v for k, v in geometry.items() if k != "positions"})
            if geometry["consistent"]:
                points = np.asarray([geometry["positions"][i] for i in largest])
                masses = np.asarray(
                    [
                        trajectory.metadata.masses_dalton[
                            trajectory.metadata.atom_mapping[i]
                        ]
                        for i in largest
                    ]
                )
                center = np.average(points, axis=0, weights=masses)
                rg.append(
                    _ok(
                        float(
                            np.sqrt(
                                np.average(
                                    np.sum((points - center) ** 2, axis=1),
                                    weights=masses,
                                )
                            )
                        ),
                        "angstrom",
                    )
                )
            else:
                rg.append(
                    _blocked(
                        "Inconsistent periodic cycle closure: winding network has no finite unique Rg"
                    )
                )
        except AnalysisError as error:
            audits.append({"consistent": False, "failure": str(error)})
            percolation.append(_blocked(str(error)))
            rg.append(_blocked(str(error)))
        cutoff = p["contact_cutoff_angstrom"]
        if cutoff is None:
            contacts.append(_blocked("Contact cutoff is not declared"))
        elif not eligible:
            contacts.append(_blocked("No eligible nonbonded heavy-atom pairs"))
        else:
            lengths = np.diag(frame.box_vectors_angstrom)
            hits = 0
            for a, b in eligible:
                delta = (
                    frame.coordinates_angstrom[trajectory.metadata.atom_mapping[b]]
                    - frame.coordinates_angstrom[trajectory.metadata.atom_mapping[a]]
                )
                delta -= np.rint(delta / lengths) * lengths
                hits += np.linalg.norm(delta) <= cutoff
            contacts.append(
                _ok(
                    hits / len(eligible),
                    "eligible_pair_fraction",
                    contacts=int(hits),
                    denominator=len(eligible),
                    cutoff_angstrom=cutoff,
                    policy="all_nonbonded_heavy_pairs_minimum_image",
                )
            )
    for name, values in (
        ("network_radius_of_gyration", rg),
        ("periodic_spanning", percolation),
        ("contact_fraction", contacts),
    ):
        result[name] = {
            "status": "ok" if all(v["status"] == "ok" for v in values) else "blocked",
            "per_frame": values,
            "value": [v["value"] for v in values],
            "frame_steps": [f.step for f in trajectory.frames],
        }
    result["periodic_geometry"] = {
        "status": "ok" if all(g["consistent"] for g in audits) else "blocked",
        "per_frame": audits,
    }
    if observables is not None:
        return {
            name: result.get(name, _blocked(f"Unsupported network observable: {name}"))
            for name in observables
        }
    return result


@dataclass(frozen=True)
class NetworkAnalysisRecord:
    json_text: str

    def __post_init__(self):
        p = json.loads(self.json_text)
        _require(
            p.get("schema") == NETWORK_RECORD_SCHEMA
            and p.get("identity") == _digest(_body(p)),
            "Network record checksum/schema mismatch",
        )
        selection = NetworkAnalysisSelection(_bytes(p["selection"]).decode())
        _require(
            p["bindings"] == {key: selection.payload[key] for key in IDENTITIES},
            "Network record binding mismatch",
        )
        _require(
            all(p.get(k) == v and type(p.get(k)) is type(v) for k, v in FLAGS.items()),
            "Network readiness flags cannot be promoted",
        )

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    def to_json(self):
        return self.json_text


def make_network_analysis_record(
    trajectory, selection, observables=None, *, provenance=None
):
    values = calculate_network_observables(trajectory, selection)
    _require(
        observables is None or values == observables,
        "Network observable payload contradicts independent recomputation",
    )
    return NetworkAnalysisRecord(
        _bytes(
            _signed(
                {
                    "schema": NETWORK_RECORD_SCHEMA,
                    "selection": selection.payload,
                    "bindings": {k: selection.payload[k] for k in IDENTITIES},
                    "observables": values,
                    "provenance": dict(provenance or {}),
                    **FLAGS,
                }
            )
        ).decode()
    )


def write_network_analysis_record(record, path):
    from island.periodic import _publish_file

    record = NetworkAnalysisRecord(record.to_json())
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    _publish_file(target, record.to_json().encode())
    return target


def load_network_analysis_record(
    path, *, expected_selection_identity=None, expected_identity=None, trajectory=None
):
    record = NetworkAnalysisRecord(Path(path).read_text())
    p = record.payload
    if expected_identity is not None:
        _require(
            record.identity == expected_identity, "Unexpected network record identity"
        )
    if expected_selection_identity is not None:
        _require(
            p["selection"]["identity"] == expected_selection_identity,
            "Unexpected network selection identity",
        )
    if trajectory is not None:
        selection = NetworkAnalysisSelection(_bytes(p["selection"]).decode())
        _require(
            calculate_network_observables(trajectory, selection) == p["observables"],
            "Network record observables mismatch",
        )
    return record


def relocate_network_analysis_record(source, destination):
    record = load_network_analysis_record(source)
    return write_network_analysis_record(record, destination)


def save_network_analysis_bundle(trajectory, selection, path):
    """Stage, fsync and atomically publish a complete, relocatable analysis bundle."""
    from island.periodic_md import _publish_directory

    selection.validate(trajectory)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    _require(
        not target.exists() and not target.is_symlink(),
        "Network bundle destination exists",
    )
    stage = Path(tempfile.mkdtemp(prefix=".network-stage-", dir=target.parent))
    try:
        record = make_network_analysis_record(trajectory, selection)
        data = {
            "schema": BUNDLE_SCHEMA,
            "metadata": trajectory.metadata.to_dict(),
            "frames": _frame_data(trajectory),
            "source_path": trajectory.source_path,
            "trajectory_file_hashes": trajectory.trajectory_file_hashes,
            "topology_bonds": trajectory.topology_bonds,
            "molecule_membership": [
                [k, v] for k, v in sorted(trajectory.molecule_membership.items())
            ],
            "record": record.payload,
            **FLAGS,
        }
        raw = _bytes(_signed(data))
        with (stage / "bundle.json").open("xb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        load_network_analysis_bundle(stage, expected_identity=record.identity)
        _publish_directory(stage, target)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return target


def load_network_analysis_bundle(path, *, expected_identity=None):
    p = json.loads((Path(path) / "bundle.json").read_text())
    _require(
        p.get("schema") == BUNDLE_SCHEMA and p.get("identity") == _digest(_body(p)),
        "Network bundle checksum/schema mismatch",
    )
    _require(
        all(p.get(k) == v and type(p.get(k)) is type(v) for k, v in FLAGS.items()),
        "Network bundle readiness mismatch",
    )
    metadata = dict(p["metadata"])
    metadata["atom_mapping"] = dict(metadata["atom_mapping"])
    trajectory = VerifiedTrajectory(
        TrajectoryMetadata(**metadata),
        tuple(
            TrajectoryFrame(
                np.asarray(f["coordinates"]),
                np.asarray(f["box"]),
                f["step"],
                f["time_ps"],
            )
            for f in p["frames"]
        ),
        p["source_path"],
        p["trajectory_file_hashes"],
        tuple(tuple(b) for b in p["topology_bonds"]),
        molecule_membership=dict(p["molecule_membership"]),
    )
    record = NetworkAnalysisRecord(_bytes(p["record"]).decode())
    if expected_identity is not None:
        _require(
            record.identity == expected_identity,
            "Unexpected network bundle record identity",
        )
    selection = NetworkAnalysisSelection(_bytes(record.payload["selection"]).decode())
    _require(
        calculate_network_observables(trajectory, selection)
        == record.payload["observables"],
        "Network bundle observable mismatch",
    )
    return trajectory, selection, record
