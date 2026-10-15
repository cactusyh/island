"""Provenance-bound post-processing for J23 periodic trajectories.

The numerical kernels in this module deliberately use NumPy arrays so that a
verified J23 workflow can be analysed without changing its graph, typing, or
charge identities.  MDAnalysis is an optional adapter for external trajectory
formats; it is never used to infer a topology or atom mapping.
"""

from __future__ import annotations

import json
import math
import shutil
from dataclasses import dataclass, field
from hashlib import sha256
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np

from island.exceptions import AnalysisError

SCHEMA = "island.periodic-analysis.v2"
LEGACY_SCHEMAS = {"island.periodic-analysis.v1"}
_MISSING = object()


def _finite(value: Any, label: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise AnalysisError(f"{label} must be finite")
    return float(value)


def _digest_bytes(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _digest_file(path: Path) -> str:
    h = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def _identity(value: Any, label: str) -> str:
    if type(value) is not str or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise AnalysisError(f"{label} must be a lowercase SHA-256 identity")
    return value


@dataclass(frozen=True)
class TrajectoryFrame:
    """One frame in Å, ps, and bar, with an orthorhombic or triclinic box."""

    coordinates_angstrom: np.ndarray
    box_vectors_angstrom: np.ndarray
    step: int
    time_ps: float
    temperature_kelvin: float | None = None
    pressure_bar: float | None = None
    potential_energy_kj_mol: float | None = None
    volume_angstrom3: float | None = None
    target_temperature_kelvin: float | None = None
    target_pressure_bar: float | None = None
    measured_temperature_kelvin: float | None = None
    measured_pressure_bar: float | None = None

    def __post_init__(self):
        xyz = np.asarray(self.coordinates_angstrom, dtype=float)
        box = np.asarray(self.box_vectors_angstrom, dtype=float)
        if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
            raise AnalysisError("Frame coordinates must be finite N x 3")
        if box.shape != (3, 3) or not np.isfinite(box).all():
            raise AnalysisError("Frame box must be finite 3 x 3")
        if abs(float(np.linalg.det(box))) <= 0:
            raise AnalysisError("Frame box must have positive volume")
        if type(self.step) is not int or self.step < 0:
            raise AnalysisError("Frame step must be a nonnegative integer")
        _finite(self.time_ps, "Frame time")
        for name in (
            "temperature_kelvin", "pressure_bar", "potential_energy_kj_mol",
            "volume_angstrom3", "target_temperature_kelvin", "target_pressure_bar",
            "measured_temperature_kelvin", "measured_pressure_bar",
        ):
            value = getattr(self, name)
            if value is not None:
                _finite(value, name)
        object.__setattr__(self, "coordinates_angstrom", xyz.copy())
        object.__setattr__(self, "box_vectors_angstrom", box.copy())
        if self.target_temperature_kelvin is None:
            object.__setattr__(self, "target_temperature_kelvin", self.temperature_kelvin)
        if self.target_pressure_bar is None:
            object.__setattr__(self, "target_pressure_bar", self.pressure_bar)

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.box_vectors_angstrom)))


@dataclass(frozen=True)
class TrajectoryMetadata:
    """Required identity and protocol metadata for scientific post-processing."""

    final_graph_identity: str
    packing_identity: str
    force_field_family: str
    force_field_source_identity: str
    typed_graph_identity: str
    charge_identity: str
    workflow_identity: str
    checkpoint_identity: str
    atom_ids: tuple[int, ...]
    masses_dalton: tuple[float, ...]
    atom_mapping: dict[int, int]
    box_periodic: tuple[bool, bool, bool] = (True, True, True)
    temperature_kelvin: float | None = None
    pressure_bar: float | None = None
    timestep_fs: float | None = None
    ensemble: str | None = None
    source_files: dict[str, str] = field(default_factory=dict)
    source_metadata: dict[str, Any] = field(default_factory=dict)
    target_temperature_kelvin: float | None = None
    target_pressure_bar: float | None = None
    measured_temperature_kelvin: float | None = None
    measured_pressure_bar: float | None = None
    assignment_identity: str | None = None

    def __post_init__(self):
        for name in (
            "final_graph_identity",
            "packing_identity",
            "force_field_source_identity",
            "typed_graph_identity",
            "charge_identity",
            "workflow_identity",
            "checkpoint_identity",
        ):
            _identity(getattr(self, name), name)
        if self.assignment_identity is not None:
            _identity(self.assignment_identity, "assignment_identity")
        if type(self.force_field_family) is not str or not self.force_field_family:
            raise AnalysisError("force_field_family must be nonempty")
        ids = tuple(self.atom_ids)
        if not ids or any(type(i) is not int for i in ids) or len(set(ids)) != len(ids):
            raise AnalysisError("atom_ids must be unique integer stable IDs")
        masses = tuple(_finite(x, "mass") for x in self.masses_dalton)
        if len(masses) != len(ids) or any(x <= 0 for x in masses):
            raise AnalysisError("masses_dalton must cover all atom IDs positively")
        mapping = dict(self.atom_mapping)
        if set(mapping) != set(ids) or sorted(mapping.values()) != list(range(len(ids))):
            raise AnalysisError("atom_mapping must be a bijection from stable IDs to trajectory indices")
        periodic = tuple(self.box_periodic)
        if len(periodic) != 3 or any(type(value) is not bool for value in periodic):
            raise AnalysisError("box_periodic must contain three boolean values")
        object.__setattr__(self, "atom_ids", ids)
        object.__setattr__(self, "masses_dalton", masses)
        object.__setattr__(self, "atom_mapping", mapping)
        object.__setattr__(self, "box_periodic", periodic)
        object.__setattr__(self, "source_files", dict(self.source_files))
        object.__setattr__(self, "source_metadata", dict(self.source_metadata))
        if self.target_temperature_kelvin is None:
            object.__setattr__(self, "target_temperature_kelvin", self.temperature_kelvin)
        if self.target_pressure_bar is None:
            object.__setattr__(self, "target_pressure_bar", self.pressure_bar)
        for name in (
            "temperature_kelvin", "pressure_bar", "timestep_fs",
            "target_temperature_kelvin", "target_pressure_bar",
            "measured_temperature_kelvin", "measured_pressure_bar",
        ):
            value = getattr(self, name)
            if value is not None:
                _finite(value, name)

    @property
    def mass_total_dalton(self) -> float:
        return float(sum(self.masses_dalton))

    def to_dict(self) -> dict[str, Any]:
        return {
            "final_graph_identity": self.final_graph_identity,
            "packing_identity": self.packing_identity,
            "force_field_family": self.force_field_family,
            "force_field_source_identity": self.force_field_source_identity,
            "typed_graph_identity": self.typed_graph_identity,
            "charge_identity": self.charge_identity,
            "assignment_identity": self.assignment_identity,
            "workflow_identity": self.workflow_identity,
            "checkpoint_identity": self.checkpoint_identity,
            "atom_ids": list(self.atom_ids),
            "masses_dalton": list(self.masses_dalton),
            "atom_mapping": [[k, self.atom_mapping[k]] for k in self.atom_ids],
            "box_periodic": list(self.box_periodic),
            "temperature_kelvin": self.temperature_kelvin,
            "pressure_bar": self.pressure_bar,
            "target_temperature_kelvin": self.target_temperature_kelvin,
            "target_pressure_bar": self.target_pressure_bar,
            "measured_temperature_kelvin": self.measured_temperature_kelvin,
            "measured_pressure_bar": self.measured_pressure_bar,
            "timestep_fs": self.timestep_fs,
            "ensemble": self.ensemble,
            "source_files": dict(self.source_files),
            "source_metadata": dict(self.source_metadata),
        }


SELECTION_SCHEMA = "island.analysis-selection.v1"


def _graph_longest_path(adjacency: dict[int, set[int]], nodes: set[int]) -> tuple[int, ...]:
    if not nodes:
        raise AnalysisError("Selection molecule has no heavy atoms")
    best: tuple[int, ...] = ()
    for start in sorted(nodes):
        queue = [(start, (start,))]
        while queue:
            current, path = queue.pop(0)
            if len(path) > len(best):
                best = path
            for neighbour in sorted(adjacency[current] & nodes):
                if neighbour not in path:
                    queue.append((neighbour, path + (neighbour,)))
    return best


@dataclass(frozen=True)
class AnalysisSelection:
    final_graph_identity: str
    atom_mapping_identity: str
    molecule_membership: dict[int, str]
    chains: tuple[dict[str, Any], ...]
    mass_weighting: str = "atomic_mass"
    chain_aggregation: str = "equal_weight"
    pbc_geometry_policy: str = "orthorhombic_bonded_path_v1"
    uncertainty: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    evidence: tuple[str, ...] = ()

    def __post_init__(self):
        _identity(self.final_graph_identity, "final_graph_identity")
        _identity(self.atom_mapping_identity, "atom_mapping_identity")
        if self.mass_weighting != "atomic_mass" or self.chain_aggregation != "equal_weight":
            raise AnalysisError("Only atomic-mass and equal-weight selection conventions are supported")
        if self.pbc_geometry_policy != "orthorhombic_bonded_path_v1":
            raise AnalysisError("Unsupported PBC geometry policy")
        membership = {int(atom_id): str(molecule) for atom_id, molecule in self.molecule_membership.items()}
        if not membership or not self.chains:
            raise AnalysisError("Analysis selection requires molecule membership and chains")
        for chain in self.chains:
            required = ("chain_id", "molecule_id", "rg_atom_ids", "backbone_path", "endpoint_ids")
            if any(key not in chain for key in required):
                raise AnalysisError("Analysis selection chain record is incomplete")
            rg_ids = tuple(int(value) for value in chain["rg_atom_ids"])
            path = tuple(int(value) for value in chain["backbone_path"])
            endpoints = tuple(int(value) for value in chain["endpoint_ids"])
            if len(path) < 2 or endpoints != (path[0], path[-1]):
                raise AnalysisError("Analysis selection has invalid backbone or endpoints")
            if any(atom_id not in membership for atom_id in rg_ids):
                raise AnalysisError("Analysis selection references an unknown atom")
        object.__setattr__(self, "molecule_membership", membership)
        object.__setattr__(self, "chains", tuple(dict(chain) for chain in self.chains))
        object.__setattr__(self, "uncertainty", dict(self.uncertainty))
        object.__setattr__(self, "provenance", dict(self.provenance))
        object.__setattr__(self, "evidence", tuple(self.evidence))

    def _body(self) -> dict[str, Any]:
        return {"schema": SELECTION_SCHEMA, "final_graph_identity": self.final_graph_identity, "atom_mapping_identity": self.atom_mapping_identity, "molecule_membership": [[atom_id, self.molecule_membership[atom_id]] for atom_id in sorted(self.molecule_membership)], "chains": list(self.chains), "mass_weighting": self.mass_weighting, "chain_aggregation": self.chain_aggregation, "pbc_geometry_policy": self.pbc_geometry_policy, "uncertainty": self.uncertainty, "provenance": self.provenance, "evidence": list(self.evidence)}

    @property
    def identity(self) -> str:
        return _digest_bytes(_json_bytes(self._body()))

    def to_dict(self) -> dict[str, Any]:
        return {**self._body(), "identity": self.identity}

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> AnalysisSelection:
        if payload.get("schema") != SELECTION_SCHEMA:
            raise AnalysisError("Unsupported analysis selection schema")
        body = dict(payload)
        stored = body.pop("identity", None)
        body.pop("schema", None)
        membership = dict(body.pop("molecule_membership", []))
        selection = cls(molecule_membership=membership, **body)
        if stored != selection.identity:
            raise AnalysisError("Analysis selection checksum/identity mismatch")
        return selection

    @classmethod
    def from_final_graph(cls, graph: dict[str, Any], *, atom_mapping_identity: str, provenance: dict[str, Any] | None = None) -> AnalysisSelection:
        graph_identity = graph.get("graph_identity") or graph.get("final_graph_identity")
        _identity(graph_identity, "final_graph_identity")
        sites = {int(site["id"]): site for site in graph.get("sites", [])}
        membership = {int(atom_id): str(molecule) for atom_id, molecule in graph.get("molecule_membership", {}).items()}
        if not sites or set(membership) != set(sites):
            raise AnalysisError("Final graph molecule membership is incomplete")
        adjacency = {atom_id: set() for atom_id in sites}
        for bond in graph.get("bonds", []):
            left, right = (int(value) for value in bond["sites"])
            adjacency[left].add(right)
            adjacency[right].add(left)
        chains = []
        for molecule in sorted(set(membership.values())):
            atom_ids = sorted(atom_id for atom_id, value in membership.items() if value == molecule)
            heavy = {atom_id for atom_id in atom_ids if sites[atom_id].get("atomic_number") != 1}
            path = _graph_longest_path(adjacency, heavy)
            if set(path) != heavy:
                raise AnalysisError("Final graph contains unsupported branched heavy-atom topology")
            chains.append({"chain_id": molecule, "molecule_id": molecule, "rg_atom_ids": atom_ids, "backbone_atom_ids": list(path), "backbone_path": list(path), "endpoint_ids": [path[0], path[-1]]})
        return cls(final_graph_identity=graph_identity, atom_mapping_identity=atom_mapping_identity, molecule_membership=membership, chains=tuple(chains), provenance=provenance or {})

    @classmethod
    def from_trajectory(cls, trajectory: VerifiedTrajectory, *, rg_atom_ids: tuple[int, ...], endpoint_ids: tuple[int, int], backbone_path: tuple[int, ...]) -> AnalysisSelection:
        membership = trajectory.molecule_membership or {atom_id: "molecule-0" for atom_id in rg_atom_ids}
        return cls(final_graph_identity=trajectory.metadata.final_graph_identity, atom_mapping_identity=trajectory.mapping_identity, molecule_membership=membership, chains=({"chain_id": "molecule-0", "molecule_id": "molecule-0", "rg_atom_ids": list(rg_atom_ids), "backbone_atom_ids": list(backbone_path), "backbone_path": list(backbone_path), "endpoint_ids": list(endpoint_ids)},))

    def validate(self, trajectory: VerifiedTrajectory) -> None:
        if self.final_graph_identity != trajectory.metadata.final_graph_identity:
            raise AnalysisError("Analysis selection final graph identity mismatch")
        if self.atom_mapping_identity != trajectory.mapping_identity:
            raise AnalysisError("Analysis selection atom mapping identity mismatch")
        actual_membership = trajectory.molecule_membership or {atom_id: molecule for atom_id, molecule in self.molecule_membership.items()}
        if any(actual_membership.get(atom_id) != molecule for atom_id, molecule in self.molecule_membership.items()):
            raise AnalysisError("Analysis selection molecule membership mismatch")
        bonds = {frozenset(bond) for bond in trajectory.topology_bonds}
        for chain in self.chains:
            if any(actual_membership.get(int(atom_id)) != chain["molecule_id"] for atom_id in chain["rg_atom_ids"]):
                raise AnalysisError("Analysis selection references atoms from the wrong molecule")
            if any(int(atom_id) not in actual_membership for atom_id in chain["backbone_path"]):
                raise AnalysisError("Analysis selection backbone contains an unknown atom")
            path = chain["backbone_path"]
            if any(frozenset((left, right)) not in bonds for left, right in pairwise(path)):
                raise AnalysisError("Analysis selection requires authoritative backbone bonds")


@dataclass(frozen=True)
class VerifiedTrajectory:
    metadata: TrajectoryMetadata
    frames: tuple[TrajectoryFrame, ...]
    source_path: str
    trajectory_file_hashes: dict[str, str]
    topology_bonds: tuple[tuple[int, int], ...] = ()
    chain_groups: tuple[tuple[int, ...], ...] = ()
    molecule_membership: dict[int, str] = field(default_factory=dict)
    selection: AnalysisSelection | None = None

    def __post_init__(self):
        if not self.frames:
            raise AnalysisError("Verified trajectory contains no frames")
        if tuple(sorted(self.frames, key=lambda f: f.step)) != self.frames:
            raise AnalysisError("Trajectory frames must be ordered by step")
        if any(left.step >= right.step for left, right in zip(self.frames, self.frames[1:])):
            raise AnalysisError("Trajectory frame steps must increase strictly")
        n = len(self.metadata.atom_ids)
        for frame in self.frames:
            if frame.coordinates_angstrom.shape[0] != n:
                raise AnalysisError("Frame atom count does not match atom mapping")
        for bond in self.topology_bonds:
            if len(bond) != 2 or bond[0] not in self.metadata.atom_mapping or bond[1] not in self.metadata.atom_mapping:
                raise AnalysisError("Topology bond references unknown stable ID")
        membership = {int(atom_id): str(value) for atom_id, value in self.molecule_membership.items()}
        if membership and set(membership) != set(self.metadata.atom_ids):
            raise AnalysisError("Molecule membership does not cover stable atom IDs")
        object.__setattr__(self, "molecule_membership", membership)
        if self.selection is not None:
            self.selection.validate(self)
        object.__setattr__(self, "trajectory_file_hashes", dict(self.trajectory_file_hashes))

    @property
    def frame_count(self) -> int:
        return len(self.frames)

    @property
    def mapping_identity(self) -> str:
        payload = {"atom_ids": list(self.metadata.atom_ids), "masses_dalton": list(self.metadata.masses_dalton), "atom_mapping": [[atom_id, self.metadata.atom_mapping[atom_id]] for atom_id in self.metadata.atom_ids]}
        return _digest_bytes(_json_bytes(payload))

    @property
    def identity(self) -> str:
        payload = {"source_path": self.source_path, "file_hashes": dict(self.trajectory_file_hashes), "frame_steps": [frame.step for frame in self.frames], "mapping_identity": self.mapping_identity, "final_graph_identity": self.metadata.final_graph_identity}
        return _digest_bytes(_json_bytes(payload))


def _workflow_hashes(root: Path) -> dict[str, str]:
    result = {}
    for path in sorted(p for p in root.rglob("*") if p.is_file() and ".writer.lock" not in p.name):
        result[str(path.relative_to(root))] = _digest_file(path)
    return result


def load_verified_trajectory(path: str | Path, *, expected_identity: str | None = None) -> VerifiedTrajectory:
    """Load and validate a J23 workflow, preserving all source identities."""
    from island.periodic_md import load_periodic_md_workflow

    root = Path(path).resolve()
    try:
        workflow = load_periodic_md_workflow(root, expected_identity=expected_identity)
        payload = workflow.payload
        manifest = payload["manifest"]
        backend = manifest["backend"]
        graph = backend["final_graph"]
        sites = graph["sites"]
        ids = tuple(manifest["atom_ids"])
        masses = tuple(float(next(site["mass"] for site in sites if site["id"] == sid)) for sid in ids)
        mapping = {sid: i for i, sid in enumerate(ids)}
        frames = []
        for item in payload["frames"]:
            frames.append(
                TrajectoryFrame(
                    np.asarray(item["positions_nm"], dtype=float) * 10.0,
                    np.asarray(item["box_vectors_nm"], dtype=float) * 10.0,
                    item["step"], item["time_ps"], item.get("temperature_kelvin"),
                    item.get("pressure_bar"), item.get("potential_energy_kj_mol"),
                )
            )
        source = (
            backend.get("source")
            or backend.get("config", {}).get("source")
            or backend.get("typed_graph", {}).get("source")
            or {}
        )
        if isinstance(source, dict):
            source_identity = source.get("sha256") or source.get("identity")
        else:
            source_identity = source
        if not source_identity:
            source_identity = backend.get("config", {}).get("source_sha256")
        if not source_identity:
            raise AnalysisError("J23 backend lacks an explicit force-field source identity")
        topology_bonds = tuple(tuple(int(value) for value in bond["sites"]) for bond in graph.get("bonds", []))
        molecule_membership = {int(atom_id): str(value) for atom_id, value in graph.get("molecule_membership", {}).items()}
        metadata = TrajectoryMetadata(
            final_graph_identity=backend["final_graph_identity"],
            packing_identity=manifest["packing_identity"],
            force_field_family=backend["config"]["family"],
            force_field_source_identity=source_identity,
            typed_graph_identity=backend["typed_graph_identity"],
            charge_identity=backend["charge_identity"],
            assignment_identity=backend.get("assignment_identity"),
            workflow_identity=manifest["identity"],
            checkpoint_identity=payload["checkpoint_identity"],
            atom_ids=ids,
            masses_dalton=masses,
            atom_mapping=mapping,
            temperature_kelvin=manifest["config"].get("temperature_kelvin"),
            pressure_bar=manifest["config"].get("pressure_bar"),
            timestep_fs=manifest["config"].get("timestep_fs"),
            ensemble=manifest["config"].get("mode"),
            source_files={"workflow_root": _digest_bytes(_json_bytes(_workflow_hashes(root)))},
            source_metadata={"production_validated": manifest["production_validated"], "simulation_readiness": manifest["simulation_readiness"]},
        )
        trajectory_hashes = _workflow_hashes(root)
        provisional = VerifiedTrajectory(metadata, tuple(frames), str(root), trajectory_hashes, topology_bonds, ())
        selection = AnalysisSelection.from_final_graph(graph, atom_mapping_identity=provisional.mapping_identity, provenance={"source": "verified_final_graph", "workflow_identity": manifest["identity"]})
        return VerifiedTrajectory(metadata, tuple(frames), str(root), trajectory_hashes, topology_bonds, tuple(tuple(chain["endpoint_ids"]) for chain in selection.chains), molecule_membership, selection)
    except AnalysisError:
        raise
    except Exception as error:
        raise AnalysisError(f"Cannot load verified J23 trajectory: {error}") from error


def load_mdanalysis_trajectory(
    topology: str | Path,
    trajectory: str | Path,
    metadata: TrajectoryMetadata,
    *,
    bonds: tuple[tuple[int, int], ...] = (),
    chain_groups: tuple[tuple[int, ...], ...] = (),
    atom_id_order: tuple[int, ...] | None = None,
) -> VerifiedTrajectory:
    """Adapt an MDAnalysis-readable trajectory with an explicit stable mapping.

    The optional dependency is imported only here.  The caller must supply all
    identities and masses; no PSMILES or topology is regenerated.
    """
    expected_order = tuple(atom_id for atom_id, _ in sorted(metadata.atom_mapping.items(), key=lambda pair: pair[1]))
    if atom_id_order is None or tuple(atom_id_order) != expected_order:
        raise AnalysisError("External trajectory atom order does not match the explicit stable atom mapping")
    for name, path in (("topology", topology), ("trajectory", trajectory)):
        expected_hash = metadata.source_files.get(name)
        if expected_hash is None:
            raise AnalysisError(f"External {name} requires a trusted source hash")
        if _digest_file(Path(path)) != expected_hash:
            raise AnalysisError(f"External {name} source hash mismatch")
    try:
        import MDAnalysis as mda
    except ImportError as error:
        raise AnalysisError("MDAnalysis is required for external trajectory formats") from error
    try:
        universe = mda.Universe(str(topology), str(trajectory))
        if len(universe.atoms) != len(metadata.atom_ids):
            raise AnalysisError("MDAnalysis atom count disagrees with atom mapping")
        if atom_id_order is None:
            raise AnalysisError("External trajectories require an explicit stable atom-ID order")
        expected_order = tuple(atom_id for atom_id, _ in sorted(metadata.atom_mapping.items(), key=lambda pair: pair[1]))
        if tuple(atom_id_order) != expected_order:
            raise AnalysisError("External trajectory atom order does not match the stable atom mapping")
        frames = []
        for ts in universe.trajectory:
            dimensions = np.asarray(ts.dimensions, dtype=float)
            if dimensions.shape != (6,) or np.any(dimensions[3:] != 90):
                raise AnalysisError("Only orthorhombic MDAnalysis unit cells are supported")
            box = np.diag(dimensions[:3])
            frames.append(TrajectoryFrame(universe.atoms.positions.copy(), box, int(ts.frame), float(ts.time)))
        hashes = {str(Path(trajectory)): _digest_file(Path(trajectory)), str(Path(topology)): _digest_file(Path(topology))}
        return VerifiedTrajectory(metadata, tuple(frames), str(Path(trajectory).resolve()), hashes, bonds, chain_groups)
    except AnalysisError:
        raise
    except Exception as error:
        raise AnalysisError(f"Cannot read MDAnalysis trajectory: {error}") from error


def _block(values: np.ndarray, block_size: int) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    if not values.size or not np.isfinite(values).all():
        raise AnalysisError("Observable values must be finite and nonempty")
    if values.size < 2:
        raise AnalysisError("At least two frames are required for uncertainty")
    if type(block_size) is not int or block_size < 1:
        raise AnalysisError("block_size must be a positive integer")
    nblocks = values.size // block_size
    if nblocks < 2:
        raise AnalysisError("At least two complete blocks are required for uncertainty")
    means = values[: nblocks * block_size].reshape(nblocks, block_size).mean(axis=1)
    return float(means.mean()), float(means.std(ddof=1) / math.sqrt(nblocks)), nblocks


def _coords(traj: VerifiedTrajectory) -> np.ndarray:
    return np.stack([f.coordinates_angstrom for f in traj.frames])


def _gyration(xyz: np.ndarray, weights: np.ndarray) -> tuple[float, np.ndarray, float]:
    w = weights / weights.sum()
    center = np.einsum("i,ij->j", w, xyz)
    centered = xyz - center
    tensor = np.einsum("i,ij,ik->jk", w, centered, centered)
    eigen = np.linalg.eigvalsh(tensor)
    rg = float(np.sqrt(max(0.0, np.trace(tensor))))
    fractions = eigen / eigen.sum() if eigen.sum() else np.zeros(3)
    anisotropy = float(1.5 * np.dot(fractions, fractions) - 0.5) if eigen.sum() else None
    return rg, tensor, anisotropy


def _orthorhombic_lengths(frame: TrajectoryFrame) -> np.ndarray:
    box = np.asarray(frame.box_vectors_angstrom, dtype=float)
    diagonal = np.diag(np.diag(box))
    if not np.allclose(box, diagonal, rtol=0.0, atol=1e-10):
        raise AnalysisError("Periodic chain unwrapping requires an orthorhombic box")
    lengths = np.diag(box)
    if np.any(lengths <= 0):
        raise AnalysisError("Periodic box lengths must be positive")
    return lengths


def _unwrap_chain(trajectory: VerifiedTrajectory, frame: TrajectoryFrame, chain: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    if len(chain) < 2:
        raise AnalysisError("A chain must contain at least two stable atom IDs")
    lengths = _orthorhombic_lengths(frame) if trajectory.metadata.box_periodic != (False, False, False) else None
    mapping = trajectory.metadata.atom_mapping
    if any(atom_id not in mapping for atom_id in chain):
        raise AnalysisError("Chain membership references an unknown stable atom ID")
    if trajectory.topology_bonds:
        adjacency = {atom_id: set() for atom_id in trajectory.metadata.atom_ids}
        for left, right in trajectory.topology_bonds:
            adjacency[left].add(right)
            adjacency[right].add(left)
        if len(chain) == 2 and not frozenset(chain) in {frozenset((left, right)) for left, right in trajectory.topology_bonds}:
            paths = [[chain[0]]]
            found = []
            while paths:
                path = paths.pop(0)
                if path[-1] == chain[1]:
                    found.append(path)
                    continue
                for neighbour in sorted(adjacency[path[-1]]):
                    if neighbour not in path:
                        paths.append(path + [neighbour])
            if len(found) != 1:
                raise AnalysisError("Chain endpoints do not identify one linear bonded path")
            chain = tuple(found[0])
        selected = set(chain)
        # Branching is assessed on the declared backbone subgraph.  Hydrogen
        # and pendant bonds are valid all-atom neighbours and must not turn a
        # linear polymer backbone into a false branched topology.
        backbone_adjacency = {
            atom_id: {neighbour for neighbour in adjacency[atom_id] if neighbour in selected}
            for atom_id in selected
        }
        if any(len(backbone_adjacency[atom_id]) > 2 for atom_id in selected):
            raise AnalysisError("branched or crosslinked periodic chains are unsupported")
        required = {frozenset((left, right)) for left, right in pairwise(chain)}
        available = {frozenset((left, right)) for left, right in trajectory.topology_bonds if left in selected and right in selected}
        if required != available or len(available) != len(chain) - 1:
            raise AnalysisError("Chain membership is not a simple linear bonded path")
    coordinates = frame.coordinates_angstrom[[mapping[atom_id] for atom_id in chain]].copy()
    unwrapped = np.empty_like(coordinates)
    unwrapped[0] = coordinates[0]
    for index in range(1, len(coordinates)):
        displacement = coordinates[index] - coordinates[index - 1]
        if lengths is not None:
            displacement -= np.rint(displacement / lengths) * lengths
        unwrapped[index] = unwrapped[index - 1] + displacement
    masses = np.asarray([trajectory.metadata.masses_dalton[mapping[atom_id]] for atom_id in chain], dtype=float)
    return unwrapped, masses


def _chain_metrics(trajectory: VerifiedTrajectory, frame: TrajectoryFrame, chain: tuple[int, ...]) -> tuple[float, float, np.ndarray, np.ndarray]:
    coordinates, masses = _unwrap_chain(trajectory, frame, chain)
    rg, tensor, _ = _gyration(coordinates, masses)
    return float(np.linalg.norm(coordinates[-1] - coordinates[0])), rg, tensor, coordinates[-1] - coordinates[0]


def _selection_chains(trajectory: VerifiedTrajectory, selection: AnalysisSelection | None) -> tuple[dict[str, Any], ...]:
    if selection is not None:
        selection.validate(trajectory)
        return selection.chains
    return tuple({"rg_atom_ids": list(chain), "backbone_path": list(chain), "endpoint_ids": [chain[0], chain[-1]]} for chain in trajectory.chain_groups if chain)


def _unwrap_selected_atoms(
    trajectory: VerifiedTrajectory,
    frame: TrajectoryFrame,
    chain_record: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return all-atom coordinates, endpoint vector, and authoritative masses."""
    path = tuple(int(value) for value in chain_record["backbone_path"])
    rg_ids = tuple(int(value) for value in chain_record["rg_atom_ids"])
    endpoints = tuple(int(value) for value in chain_record["endpoint_ids"])
    if len(path) < 2 or endpoints != (path[0], path[-1]):
        raise AnalysisError("selection endpoint IDs do not match the ordered backbone path")
    mapping = trajectory.metadata.atom_mapping
    if any(atom_id not in mapping for atom_id in (*path, *rg_ids)):
        raise AnalysisError("selection references an unknown stable atom ID")
    bonds = {frozenset((int(left), int(right))) for left, right in trajectory.topology_bonds}
    if trajectory.metadata.box_periodic != (False, False, False) and not bonds:
        raise AnalysisError("periodic graph-bound analysis requires authoritative bonds")
    # Start with the bonded backbone, then place hydrogens and pendants from
    # their already placed neighbours.  This deliberately rejects ambiguous
    # periodic graph geometry rather than inventing an image.
    backbone, _ = _unwrap_chain(trajectory, frame, path)
    positions = {atom_id: coordinate.copy() for atom_id, coordinate in zip(path, backbone)}
    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in rg_ids}
    for bond in trajectory.topology_bonds:
        left, right = bond
        if left in adjacency and right in adjacency:
            adjacency[left].add(right)
            adjacency[right].add(left)
    pending = set(rg_ids) - set(positions)
    lengths = _orthorhombic_lengths(frame) if trajectory.metadata.box_periodic != (False, False, False) else None
    coordinates = frame.coordinates_angstrom
    while pending:
        candidates = []
        for atom_id in sorted(pending):
            neighbours = [neighbour for neighbour in adjacency.get(atom_id, ()) if neighbour in positions]
            if len(neighbours) == 1:
                candidates.append((atom_id, neighbours[0]))
            elif len(neighbours) > 1:
                raise AnalysisError("ambiguous periodic images in all-atom selection")
        if not candidates:
            raise AnalysisError("selection atoms lack a connected authoritative bond path")
        for atom_id, neighbour in candidates:
            displacement = coordinates[mapping[atom_id]] - coordinates[mapping[neighbour]]
            if lengths is not None:
                displacement -= np.rint(displacement / lengths) * lengths
            positions[atom_id] = positions[neighbour] + displacement
            pending.remove(atom_id)
    ordered = np.asarray([positions[atom_id] for atom_id in rg_ids], dtype=float)
    masses = np.asarray([trajectory.metadata.masses_dalton[mapping[atom_id]] for atom_id in rg_ids], dtype=float)
    endpoint_vector = positions[endpoints[1]] - positions[endpoints[0]]
    return ordered, endpoint_vector, masses


def _selected_chain_metrics(trajectory: VerifiedTrajectory, frame: TrajectoryFrame, chain_record: dict[str, Any]) -> tuple[float, float, np.ndarray, np.ndarray]:
    coordinates, endpoint_vector, masses = _unwrap_selected_atoms(trajectory, frame, chain_record)
    rg, tensor, _ = _gyration(coordinates, masses)
    return float(np.linalg.norm(endpoint_vector)), rg, tensor, endpoint_vector


def _rdf(traj: VerifiedTrajectory, bins: int, r_max: float) -> dict[str, Any]:
    if bins < 2 or r_max <= 0:
        raise AnalysisError("RDF bins and r_max must be positive")
    edges = np.linspace(0.0, r_max, bins + 1)
    counts = np.zeros(bins, dtype=float)
    n = len(traj.metadata.atom_ids)
    for frame in traj.frames:
        if not np.allclose(frame.box_vectors_angstrom, np.diag(np.diag(frame.box_vectors_angstrom))):
            raise AnalysisError("RDF requires orthorhombic periodic boxes")
        xyz = frame.coordinates_angstrom
        box = np.diag(frame.box_vectors_angstrom)
        for i in range(n - 1):
            delta = xyz[i + 1 :] - xyz[i]
            delta -= np.rint(delta / box) * box
            distances = np.linalg.norm(delta, axis=1)
            counts += np.histogram(distances[distances < r_max], edges)[0]
    radii = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 * math.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    density = n / np.mean([f.volume for f in traj.frames])
    normalizer = len(traj.frames) * n * density * shell / 2.0
    values = np.divide(counts, normalizer, out=np.zeros_like(counts), where=normalizer > 0)
    return {"r_angstrom": radii.tolist(), "g_r": values.tolist(), "counts": counts.tolist(), "bins": bins, "r_max_angstrom": r_max}


def calculate_observables(
    trajectory: VerifiedTrajectory,
    names: tuple[str, ...] | list[str],
    *,
    block_size: int = 5,
    rdf_bins: int = 100,
    rdf_r_max_angstrom: float = 10.0,
    selection: AnalysisSelection | None = None,
) -> dict[str, dict[str, Any]]:
    """Calculate requested descriptors with block uncertainty and failure states."""
    if type(trajectory) is not VerifiedTrajectory:
        raise AnalysisError("Expected VerifiedTrajectory")
    selection = selection or trajectory.selection
    if selection is not None:
        selection.validate(trajectory)
    chain_records = _selection_chains(trajectory, selection)
    names = tuple(names)
    known = set(observable_definitions())
    unknown = set(names) - known
    if unknown:
        raise AnalysisError(f"Unknown observables: {sorted(unknown)}")
    xyz = _coords(trajectory)
    results: dict[str, dict[str, Any]] = {}

    def scalar(name: str, values: np.ndarray, units: str, definition: str, requirements: str):
        mean, uncertainty, blocks = _block(values, block_size)
        results[name] = {"status": "ok", "observable": name, "definition": definition, "units": units, "value": mean, "uncertainty": uncertainty, "uncertainty_method": "standard error of complete block means", "blocks": blocks, "requirements": requirements}

    if "density" in names:
        values = np.asarray([trajectory.metadata.mass_total_dalton / f.volume * 1.66053906660 for f in trajectory.frames])
        scalar("density", values, "g/cm^3", "total system mass divided by instantaneous periodic volume", "periodic box and authoritative atomic masses")
    if "density_convergence" in names:
        values = np.asarray([trajectory.metadata.mass_total_dalton / f.volume * 1.66053906660 for f in trajectory.frames])
        mean, uncertainty, blocks = _block(values, block_size)
        results["density_convergence"] = {"status": "ok", "observable": "density_convergence", "definition": "block means, drift, and plateau diagnostic for density", "units": "g/cm^3", "value": mean, "uncertainty": uncertainty, "blocks": blocks, "block_means": values[: blocks * block_size].reshape(blocks, block_size).mean(axis=1).tolist(), "requirements": "at least two complete blocks; plateau decision remains an explicit gate"}
    if "radius_of_gyration" in names or "gyration_tensor" in names or "chain_anisotropy" in names:
        if not chain_records and trajectory.metadata.box_periodic == (False, False, False):
            values = []
            tensors = []
            anisotropy = []
            for frame_xyz in xyz:
                rg, tensor, kappa = _gyration(frame_xyz, np.asarray(trajectory.metadata.masses_dalton, dtype=float))
                values.append(rg)
                tensors.append(tensor)
                anisotropy.append(kappa)
            if "radius_of_gyration" in names:
                scalar("radius_of_gyration", np.asarray(values), "angstrom", "mass-weighted sqrt(trace(gyration tensor)) for the nonperiodic whole system", "nonperiodic coordinates and authoritative masses")
            if "gyration_tensor" in names:
                tensor_array = np.asarray(tensors)
                results["gyration_tensor"] = {"status": "ok", "observable": "gyration_tensor", "definition": "mass-weighted covariance of nonperiodic whole-system coordinates", "units": "angstrom^2", "value": tensor_array.mean(axis=0).tolist(), "uncertainty": tensor_array.std(axis=0, ddof=1).tolist() if len(tensor_array) > 1 else None, "uncertainty_method": "frame standard deviation", "requirements": "nonperiodic coordinates and authoritative masses"}
            if "chain_anisotropy" in names:
                scalar("chain_anisotropy", np.asarray([value for value in anisotropy if value is not None]), "dimensionless", "relative anisotropy of the nonperiodic whole-system gyration tensor", "nonzero gyration tensor")
        elif not chain_records:
            for name in ("radius_of_gyration", "gyration_tensor", "chain_anisotropy"):
                if name in names:
                    results[name] = {"status": "blocked", "observable": name, "definition": observable_definitions()[name]["definition"], "units": observable_definitions()[name]["units"], "value": None, "uncertainty": None, "failure": "periodic Rg analysis requires explicit chain_groups", "requirements": observable_definitions()[name]["requirements"]}
        else:
            radii = []
            tensors = []
            anisotropy = []
            try:
                for frame in trajectory.frames:
                    frame_radii = []
                    frame_tensors = []
                    frame_anisotropy = []
                    for chain in chain_records:
                        if selection is not None:
                            _, rg, tensor, _ = _selected_chain_metrics(trajectory, frame, chain)
                        else:
                            _, rg, tensor, _ = _chain_metrics(trajectory, frame, tuple(chain["backbone_path"]))
                        frame_radii.append(rg)
                        frame_tensors.append(tensor)
                        eigenvalues = np.linalg.eigvalsh(tensor)
                        total = float(eigenvalues.sum())
                        frame_anisotropy.append(float(1.5 * np.dot(eigenvalues / total, eigenvalues / total) - 0.5) if total > 0 else None)
                    radii.append(frame_radii)
                    tensors.append(frame_tensors)
                    anisotropy.append(frame_anisotropy)
                radii_array = np.asarray(radii, dtype=float)
                if not np.isfinite(radii_array).all() or np.any(radii_array <= 0):
                    raise AnalysisError("radius of gyration contains a zero or nonfinite value")
                if "radius_of_gyration" in names:
                    mean, uncertainty, blocks = _block(radii_array.mean(axis=1), block_size)
                    results["radius_of_gyration"] = {"status": "ok", "observable": "radius_of_gyration", "definition": "mass-weighted square root of the gyration-tensor trace for each unwrapped chain", "units": "angstrom", "value": mean, "uncertainty": uncertainty, "uncertainty_method": "standard error of complete block means over equal-weight per-chain frame means", "blocks": blocks, "per_chain": radii_array.mean(axis=0).tolist(), "requirements": "explicit linear chain groups and orthorhombic bonded-path unwrapping"}
                if "gyration_tensor" in names:
                    tensor_array = np.asarray(tensors, dtype=float)
                    frame_tensor_means = tensor_array.mean(axis=1)
                    results["gyration_tensor"] = {"status": "ok", "observable": "gyration_tensor", "definition": "mass-weighted covariance of unwrapped coordinates for each chain, equal-weight averaged per frame", "units": "angstrom^2", "value": frame_tensor_means.mean(axis=0).tolist(), "uncertainty": frame_tensor_means.std(axis=0, ddof=1).tolist() if len(frame_tensor_means) > 1 else None, "uncertainty_method": "frame standard deviation", "requirements": "explicit linear chain groups and orthorhombic bonded-path unwrapping"}
                if "chain_anisotropy" in names:
                    anisotropy_array = np.asarray([[value for value in frame if value is not None] for frame in anisotropy], dtype=float)
                    if anisotropy_array.size:
                        scalar("chain_anisotropy", anisotropy_array.mean(axis=1), "dimensionless", "relative anisotropy of each unwrapped chain, equal-weight averaged per frame", "nonzero unwrapped gyration tensor")
                    else:
                        raise AnalysisError("gyration tensor has zero spread")
            except AnalysisError as error:
                for name in ("radius_of_gyration", "gyration_tensor", "chain_anisotropy"):
                    if name in names and name not in results:
                        results[name] = {"status": "blocked", "observable": name, "definition": observable_definitions()[name]["definition"], "units": observable_definitions()[name]["units"], "value": None, "uncertainty": None, "failure": str(error), "requirements": observable_definitions()[name]["requirements"]}
    if "end_to_end_distance" in names:
        if not chain_records:
            results["end_to_end_distance"] = {"status": "blocked", "observable": "end_to_end_distance", "definition": "mean explicit chain endpoint distance", "units": "angstrom", "value": None, "uncertainty": None, "failure": "explicit chain_groups are absent", "requirements": "one or more nonempty stable-ID chain endpoint mappings"}
        else:
            values = []
            try:
                per_frame = []
                for frame in trajectory.frames:
                    chain_values = [(_selected_chain_metrics(trajectory, frame, chain)[0] if selection is not None else _chain_metrics(trajectory, frame, tuple(chain["backbone_path"]))[0]) for chain in chain_records]
                    per_frame.append(chain_values)
                values = np.asarray(per_frame, dtype=float).mean(axis=1)
                scalar("end_to_end_distance", np.asarray(values), "angstrom", "mean explicit chain endpoint distance", "explicit stable endpoint mapping")
            except (AnalysisError, KeyError, IndexError) as error:
                results["end_to_end_distance"] = {"status": "blocked", "observable": "end_to_end_distance", "definition": "mean explicit chain endpoint distance", "units": "angstrom", "value": None, "uncertainty": None, "failure": str(error), "requirements": "one or more nonempty stable-ID chain endpoint mappings"}
    if "re_over_rg" in names:
        if not chain_records:
            results["re_over_rg"] = {"status": "blocked", "observable": "re_over_rg", "definition": "mean of per-chain Re/Rg values over frames", "units": "dimensionless", "value": None, "uncertainty": None, "failure": "explicit chain_groups are absent", "requirements": "linear chain endpoint and bonded-path mappings"}
        else:
            try:
                ratios = []
                for frame in trajectory.frames:
                    frame_ratios = []
                    for chain in chain_records:
                        if selection is not None:
                            end_to_end, rg, _, _ = _selected_chain_metrics(trajectory, frame, chain)
                        else:
                            end_to_end, rg, _, _ = _chain_metrics(trajectory, frame, tuple(chain["backbone_path"]))
                        if not np.isfinite(rg) or rg <= 0:
                            raise AnalysisError("radius of gyration contains a zero or nonfinite denominator")
                        frame_ratios.append(end_to_end / rg)
                    ratios.append(frame_ratios)
                ratios = np.asarray(ratios, dtype=float)
                mean, uncertainty, blocks = _block(ratios.mean(axis=1), block_size)
                results["re_over_rg"] = {"status": "ok", "observable": "re_over_rg", "definition": "mean of per-chain Re/Rg values over frames; chains receive equal weight", "units": "dimensionless", "value": mean, "uncertainty": uncertainty, "uncertainty_method": "standard error of complete block means over per-frame chain means", "blocks": blocks, "per_chain": ratios.mean(axis=0).tolist(), "requirements": "linear chain endpoint and bonded-path mappings"}
            except AnalysisError as error:
                results["re_over_rg"] = {"status": "blocked", "observable": "re_over_rg", "definition": "mean of per-chain Re/Rg values over frames", "units": "dimensionless", "value": None, "uncertainty": None, "failure": str(error), "requirements": "linear chain endpoint and bonded-path mappings"}
    if "radial_distribution_function" in names:
        try:
            value = _rdf(trajectory, rdf_bins, rdf_r_max_angstrom)
            results["radial_distribution_function"] = {"status": "ok", "observable": "radial_distribution_function", "definition": "minimum-image pair distribution normalized by ideal-gas shell counts", "units": "dimensionless", "value": value, "uncertainty": None, "uncertainty_method": "pair-count accumulation; block bootstrap required for publication", "requirements": "orthorhombic periodic cells and sufficient independent frames"}
        except AnalysisError as error:
            results["radial_distribution_function"] = {"status": "blocked", "observable": "radial_distribution_function", "definition": "minimum-image pair distribution normalized by ideal-gas shell counts", "units": "dimensionless", "value": None, "uncertainty": None, "failure": str(error), "requirements": "orthorhombic periodic cells and sufficient independent frames"}
    if "mean_squared_displacement" in names:
        if not trajectory.metadata.source_metadata.get("unwrapped_coordinates", False):
            results["mean_squared_displacement"] = {"status": "blocked", "observable": "mean_squared_displacement", "definition": "mean squared displacement from validated unwrapped coordinates", "units": "angstrom^2", "value": None, "uncertainty": None, "failure": "unwrapped coordinates or an explicit minimum-image unwrapping protocol are absent", "requirements": "validated unwrapped coordinates and at least three frames"}
        else:
            if len(trajectory.frames) < 3:
                raise AnalysisError("mean_squared_displacement requires at least three frames")
            origin = xyz[0]
            values = np.mean(np.sum((xyz - origin[None, :, :]) ** 2, axis=2), axis=1)
            results["mean_squared_displacement"] = {"status": "ok", "observable": "mean_squared_displacement", "definition": "mean squared displacement relative to first retained frame", "units": "angstrom^2", "value": values.tolist(), "uncertainty": None, "requirements": "validated unwrapped coordinates and at least three frames"}
    if "temperature_density_volume" in names:
        temperatures = [f.measured_temperature_kelvin for f in trajectory.frames]
        if any(value is None for value in temperatures):
            results["temperature_density_volume"] = {"status": "blocked", "observable": "temperature_density_volume", "definition": "temperature-dependent density and volume profile using measured frame temperatures", "units": {"temperature": "K", "density": "g/cm^3", "volume": "angstrom^3"}, "value": None, "uncertainty": None, "failure": "measured temperature is absent from one or more frames; target thermostat values are not measurements", "requirements": "trajectory with measured temperature per frame"}
        else:
            results["temperature_density_volume"] = {"status": "ok", "observable": "temperature_density_volume", "definition": "temperature-dependent density and volume profile using measured frame temperatures", "units": {"temperature": "K", "density": "g/cm^3", "volume": "angstrom^3"}, "value": [{"measured_temperature_kelvin": f.measured_temperature_kelvin, "density_g_cm3": trajectory.metadata.mass_total_dalton / f.volume * 1.66053906660, "volume_angstrom3": f.volume} for f in trajectory.frames], "uncertainty": None, "requirements": "trajectory with measured temperature per frame"}
    if "nematic_order" in names:
        if not chain_records:
            results["nematic_order"] = {"status": "blocked", "observable": "nematic_order", "definition": observable_definitions()["nematic_order"]["definition"], "units": "dimensionless", "value": None, "uncertainty": None, "failure": "chain orientation mapping absent", "requirements": observable_definitions()["nematic_order"]["requirements"]}
        else:
            order_values = []
            for frame in trajectory.frames:
                axes = []
                for chain in chain_records:
                    if len(chain) < 2:
                        continue
                    try:
                        axis = (_selected_chain_metrics(trajectory, frame, chain)[3] if selection is not None else _chain_metrics(trajectory, frame, tuple(chain["backbone_path"]))[3])
                    except AnalysisError:
                        axes = []
                        break
                    norm = np.linalg.norm(axis)
                    if norm > 0:
                        axes.append(axis / norm)
                if not axes:
                    break
                q = np.mean([1.5 * np.outer(axis, axis) - 0.5 * np.eye(3) for axis in axes], axis=0)
                order_values.append(float(np.linalg.eigvalsh(q)[-1]))
            if len(order_values) < 2:
                results["nematic_order"] = {"status": "blocked", "observable": "nematic_order", "definition": observable_definitions()["nematic_order"]["definition"], "units": "dimensionless", "value": None, "uncertainty": None, "failure": "chain orientation vectors are unavailable", "requirements": observable_definitions()["nematic_order"]["requirements"]}
            else:
                scalar("nematic_order", np.asarray(order_values), "dimensionless", observable_definitions()["nematic_order"]["definition"], observable_definitions()["nematic_order"]["requirements"])
    if "chain_orientation" in names:
        if not chain_records:
            results["chain_orientation"] = {"status": "blocked", "observable": "chain_orientation", "definition": observable_definitions()["chain_orientation"]["definition"], "units": "unit vector", "value": None, "uncertainty": None, "failure": "chain orientation mapping absent", "requirements": observable_definitions()["chain_orientation"]["requirements"]}
        else:
            orientations = []
            for frame in trajectory.frames:
                vectors = []
                for chain in chain_records:
                    if len(chain) < 2:
                        continue
                    try:
                        axis = (_selected_chain_metrics(trajectory, frame, chain)[3] if selection is not None else _chain_metrics(trajectory, frame, tuple(chain["backbone_path"]))[3])
                    except AnalysisError:
                        vectors = []
                        break
                    norm = np.linalg.norm(axis)
                    if norm:
                        vectors.append((axis / norm).tolist())
                orientations.append(vectors)
            results["chain_orientation"] = {"status": "ok" if any(orientations) else "blocked", "observable": "chain_orientation", "definition": observable_definitions()["chain_orientation"]["definition"], "units": "unit vector", "value": orientations if any(orientations) else None, "uncertainty": None, "failure": None if any(orientations) else "chain orientation vectors are unavailable", "requirements": observable_definitions()["chain_orientation"]["requirements"]}
    if "contact_fraction" in names:
        cutoff = trajectory.metadata.source_metadata.get("contact_cutoff_angstrom")
        if type(cutoff) not in (int, float) or cutoff <= 0:
            results["contact_fraction"] = {"status": "blocked", "observable": "contact_fraction", "definition": observable_definitions()["contact_fraction"]["definition"], "units": "dimensionless", "value": None, "uncertainty": None, "failure": "contact_cutoff_angstrom is absent", "requirements": observable_definitions()["contact_fraction"]["requirements"]}
        else:
            bonded = {tuple(sorted((trajectory.metadata.atom_mapping[a], trajectory.metadata.atom_mapping[b]))) for a, b in trajectory.topology_bonds}
            pairs = [(i, j) for i in range(xyz.shape[1]) for j in range(i + 1, xyz.shape[1]) if (i, j) not in bonded]
            if not pairs:
                results["contact_fraction"] = {"status": "blocked", "observable": "contact_fraction", "definition": observable_definitions()["contact_fraction"]["definition"], "units": "dimensionless", "value": None, "uncertainty": None, "failure": "no nonbonded pairs remain after exclusions", "requirements": observable_definitions()["contact_fraction"]["requirements"]}
            else:
                values = []
                for frame, frame_xyz in zip(trajectory.frames, xyz):
                    box = np.diag(frame.box_vectors_angstrom)
                    count = 0
                    for i, j in pairs:
                        delta = frame_xyz[j] - frame_xyz[i]
                        delta -= np.rint(delta / box) * box
                        count += float(np.linalg.norm(delta) < cutoff)
                    values.append(count / len(pairs))
                scalar("contact_fraction", np.asarray(values), "dimensionless", observable_definitions()["contact_fraction"]["definition"], observable_definitions()["contact_fraction"]["requirements"])
    if "crosslink_density" in names:
        count = trajectory.metadata.source_metadata.get("crosslink_count")
        convention = trajectory.metadata.source_metadata.get("crosslink_density_convention")
        if type(count) is int and count >= 0 and convention in {"per_volume_angstrom3", "per_mass_dalton"}:
            denominator = trajectory.frames[-1].volume if convention == "per_volume_angstrom3" else trajectory.metadata.mass_total_dalton
            results["crosslink_density"] = {"status": "ok", "observable": "crosslink_density", "definition": observable_definitions()["crosslink_density"]["definition"], "units": convention, "value": count / denominator, "uncertainty": None, "requirements": observable_definitions()["crosslink_density"]["requirements"]}
    if "network_connectivity" in names and trajectory.topology_bonds:
            adjacency = {atom_id: set() for atom_id in trajectory.metadata.atom_ids}
            for left, right in trajectory.topology_bonds:
                adjacency[left].add(right); adjacency[right].add(left)
            components = []
            unseen = set(adjacency)
            while unseen:
                stack = [unseen.pop()]; component = []
                while stack:
                    current = stack.pop(); component.append(current)
                    for neighbour in adjacency[current] & unseen:
                        unseen.remove(neighbour); stack.append(neighbour)
                components.append(component)
            degrees = [len(neighbours) for neighbours in adjacency.values()]
            results["network_connectivity"] = {"status": "ok", "observable": "network_connectivity", "definition": observable_definitions()["network_connectivity"]["definition"], "units": "count/fraction", "value": {"component_count": len(components), "component_sizes": sorted((len(component) for component in components), reverse=True), "mean_degree": float(np.mean(degrees))}, "uncertainty": None, "requirements": observable_definitions()["network_connectivity"]["requirements"]}
    blocked = {"torsional_state_statistics", "cohesive_energy_density", "fractional_free_volume", "crosslink_density", "network_connectivity", "relaxation_time", "glass_transition", "mechanical_stress"}
    for name in names:
        if name in blocked and name not in results:
            results[name] = {"status": "blocked", "observable": name, "definition": observable_definitions()[name]["definition"], "units": observable_definitions()[name]["units"], "value": None, "uncertainty": None, "failure": observable_definitions()[name]["failure"], "requirements": observable_definitions()[name]["requirements"]}
    return results


def validate_physical_gates(
    trajectory: VerifiedTrajectory,
    observables: dict[str, dict[str, Any]] | None = None,
    *,
    replica_trajectories: tuple[VerifiedTrajectory, ...] = (),
) -> dict[str, dict[str, Any]]:
    """Evaluate explicit validation gates without promoting readiness.

    A gate can pass as an internal consistency check while the aggregate remains
    blocked when the source workflow retains its conservative readiness flags.
    """
    if type(trajectory) is not VerifiedTrajectory:
        raise AnalysisError("Expected VerifiedTrajectory")
    gates: dict[str, dict[str, Any]] = {}
    gates["identity"] = {"status": "passed", "reason": "graph, mapping, source and checkpoint identities validated on load"}
    finite = all(np.isfinite(frame.coordinates_angstrom).all() and np.isfinite(frame.box_vectors_angstrom).all() for frame in trajectory.frames)
    gates["finite_state"] = {"status": "passed" if finite else "blocked", "reason": "all retained coordinates and boxes are finite" if finite else "nonfinite retained state"}
    ordered = all(a.step < b.step for a, b in zip(trajectory.frames, trajectory.frames[1:]))
    gates["trajectory_order"] = {"status": "passed" if ordered else "blocked", "reason": "strictly increasing accepted steps" if ordered else "duplicate or reversed accepted steps"}
    density = (observables or {}).get("density_convergence")
    plateau = density is not None and density.get("status") == "ok" and density.get("blocks", 0) >= 2 and density.get("plateau") is True
    gates["density_plateau"] = {"status": "passed" if plateau else "blocked", "reason": "declared plateau criterion passed" if plateau else "a block estimate exists only; an explicit plateau criterion is required"}
    source_flags = trajectory.metadata.source_metadata
    for gate_name, field_name, label in (
        ("energy_stability", "energy_stability_validated", "energy"),
        ("temperature_stability", "temperature_stability_validated", "temperature"),
        ("pressure_stability", "pressure_stability_validated", "pressure"),
        ("volume_stability", "volume_stability_validated", "volume"),
        ("chain_dimension_convergence", "chain_dimension_convergence_validated", "chain dimensions"),
        ("rdf_consistency", "rdf_consistency_validated", "RDF"),
        ("topology_preservation", "topology_validated", "bond and valence topology"),
        ("crosslink_preservation", "crosslink_topology_validated", "crosslink topology"),
    ):
        passed = source_flags.get(field_name) is True
        gates[gate_name] = {"status": "passed" if passed else "blocked", "reason": f"{label} validation is explicitly recorded" if passed else f"independent {label} validation is absent"}
    replica_ok = True
    if replica_trajectories:
        replica_ok = all(t.metadata.final_graph_identity == trajectory.metadata.final_graph_identity and t.metadata.atom_ids == trajectory.metadata.atom_ids for t in replica_trajectories)
    else:
        replica_ok = False
    gates["replica_agreement"] = {"status": "passed" if replica_ok else "blocked", "reason": "independent replica identities supplied" if replica_ok else "independent replicas are required"}
    readiness_ok = source_flags.get("production_validated") is True and source_flags.get("simulation_readiness") == "established"
    gates["aggregate"] = {"status": "passed" if all(g["status"] == "passed" for name, g in gates.items() if name != "aggregate") and readiness_ok else "blocked", "reason": "all gates and source readiness flags support physical interpretation" if readiness_ok else "physical production readiness is not established"}
    return gates


def observable_definitions() -> dict[str, dict[str, str]]:
    """Definitions, units, requirements and explicit failure states."""
    return {
        "density": {"definition": "mass / periodic volume", "units": "g/cm^3", "requirements": "periodic box and masses", "failure": "missing box or masses"},
        "density_convergence": {"definition": "block density means and plateau diagnostic", "units": "g/cm^3", "requirements": "two or more complete blocks", "failure": "insufficient frames or no plateau evidence"},
        "radius_of_gyration": {"definition": "mass-weighted sqrt(trace(gyration tensor)) per unwrapped chain and frame, then equal-weight chain average; nonperiodic systems may use the whole system", "units": "angstrom", "requirements": "periodic linear chain groups and bonded-path unwrapping, or nonperiodic coordinates", "failure": "explicit chain groups absent for periodic data or valid coordinates absent"},
        "end_to_end_distance": {"definition": "minimum-image distance between explicitly mapped chain endpoints per chain and frame, then equal-weight chain average", "units": "angstrom", "requirements": "linear chain endpoint mapping and orthorhombic bonded path", "failure": "endpoint provenance absent or periodic graph is branched"},
        "re_over_rg": {"definition": "per-chain Re/Rg from unwrapped coordinates, then equal-weight average across chains and frames", "units": "dimensionless", "requirements": "linear chain endpoints, bonded paths and strictly positive per-chain Rg", "failure": "endpoint provenance absent, periodic graph is branched, or Rg is zero/nonfinite"},
        "gyration_tensor": {"definition": "mass-weighted coordinate covariance", "units": "angstrom^2", "requirements": "stable atom IDs and coordinates", "failure": "missing coordinates or masses"},
        "chain_anisotropy": {"definition": "relative eigenvalue anisotropy of each unwrapped chain gyration tensor, equal-weight averaged per frame", "units": "dimensionless", "requirements": "periodic linear chain groups and nonzero unwrapped gyration tensors", "failure": "zero-spread selection"},
        "chain_orientation": {"definition": "normalized first-to-last vector for each chain from unwrapped bonded-path coordinates, retained per chain and frame", "units": "unit vector", "requirements": "stable-ID chain groups with at least two atoms", "failure": "chain orientation mapping absent or graph is branched"},
        "radial_distribution_function": {"definition": "minimum-image pair distribution g(r)", "units": "dimensionless", "requirements": "orthorhombic PBC and independent frames", "failure": "triclinic or insufficient trajectory"},
        "mean_squared_displacement": {"definition": "mean squared displacement from first retained frame", "units": "angstrom^2", "requirements": "unwrapped coordinates or explicit unwrapping protocol", "failure": "insufficient frames or unvalidated wrapping"},
        "nematic_order": {"definition": "largest eigenvalue of the second-rank tensor formed from per-chain unwrapped orientation vectors, averaged per frame", "units": "dimensionless", "requirements": "chain axes and orientation mapping", "failure": "chain orientation mapping absent or graph is branched"},
        "contact_fraction": {"definition": "fraction of nonbonded stable-ID pairs below declared cutoff", "units": "dimensionless", "requirements": "cutoff, exclusions, and independent frames", "failure": "contact policy absent"},
        "torsional_state_statistics": {"definition": "declared torsion bins and transition counts", "units": "radian/count", "requirements": "authoritative torsion mapping and trajectory", "failure": "torsion mapping absent"},
        "cohesive_energy_density": {"definition": "physically defined cohesive energy per volume", "units": "energy/volume", "requirements": "validated energy decomposition and convention", "failure": "potential energy alone is insufficient"},
        "fractional_free_volume": {"definition": "1 - occupied probe volume / cell volume under explicit algorithm", "units": "dimensionless", "requirements": "probe radius, grid/sampling rule, convergence", "failure": "algorithm and uncertainty not declared"},
        "crosslink_density": {"definition": "crosslink count per declared mass or volume convention", "units": "declared convention", "requirements": "immutable crosslink record", "failure": "crosslink provenance absent"},
        "network_connectivity": {"definition": "component and degree statistics of final network", "units": "count/fraction", "requirements": "final graph and crosslink topology", "failure": "network topology absent"},
        "relaxation_time": {"definition": "fit to a declared autocorrelation or stress relaxation model", "units": "ps", "requirements": "long stationary trajectory and fit diagnostics", "failure": "trajectory length/convergence insufficient"},
        "temperature_density_volume": {"definition": "temperature-dependent density and volume profile from measured frame temperatures", "units": "K, g/cm^3, angstrom^3", "requirements": "measured temperature per frame and equilibrated windows", "failure": "measured temperature scan absent; target values are insufficient"},
        "glass_transition": {"definition": "piecewise fit to converged density or volume versus temperature", "units": "K", "requirements": "multiple equilibrated temperatures and fit uncertainty", "failure": "J23 does not establish a converged temperature scan"},
        "mechanical_stress": {"definition": "stress or modulus from a declared deformation protocol", "units": "pressure", "requirements": "nonequilibrium deformation and virial convention", "failure": "equilibrium trajectory alone is insufficient"},
    }


@dataclass(frozen=True)
class AnalysisRecord:
    payload: dict[str, Any]

    def __post_init__(self):
        if self.payload.get("schema") not in {SCHEMA, *LEGACY_SCHEMAS}:
            raise AnalysisError("Unsupported analysis record schema")
        if self.payload.get("production_validated") is not False or self.payload.get("simulation_readiness") != "not_established":
            raise AnalysisError("Analysis records cannot promote scientific readiness")
        if self.payload.get("schema") == SCHEMA:
            bindings = self.payload.get("bindings")
            if not isinstance(bindings, dict):
                raise AnalysisError("Analysis record bindings are missing")
            for name in ("final_graph_identity", "packing_identity", "force_field_source_identity", "typed_graph_identity", "charge_identity", "workflow_identity", "checkpoint_identity", "atom_mapping_identity", "trajectory_identity"):
                _identity(bindings.get(name), name)
            if bindings.get("assignment_identity") is not None:
                _identity(bindings["assignment_identity"], "assignment_identity")
            for name, value in bindings.get("source_hashes", {}).items():
                _identity(value, f"source hash {name}")
            metadata = self.payload["trajectory"]["metadata"]
            for name in ("final_graph_identity", "packing_identity", "force_field_source_identity", "typed_graph_identity", "charge_identity", "workflow_identity", "checkpoint_identity"):
                if bindings[name] != metadata.get(name):
                    raise AnalysisError(f"Analysis binding mismatch: {name}")
            if bindings.get("assignment_identity") != metadata.get("assignment_identity"):
                raise AnalysisError("Analysis binding mismatch: assignment_identity")
            mapping_payload = {name: metadata[name] for name in ("atom_ids", "masses_dalton", "atom_mapping")}
            if _digest_bytes(_json_bytes(mapping_payload)) != bindings["atom_mapping_identity"]:
                raise AnalysisError("Analysis binding mismatch: atom_mapping_identity")
            for name, observable in self.payload.get("observables", {}).items():
                if not observable.get("units"):
                    raise AnalysisError(f"Observable units are absent: {name}")
            _validate_finite_payload(self.payload.get("observables", {}))
            selection_payload = self.payload.get("selection")
            if selection_payload is not None:
                selection = AnalysisSelection.from_dict(selection_payload)
                if self.payload.get("bindings", {}).get("selection_identity") != selection.identity:
                    raise AnalysisError("Analysis binding mismatch: selection_identity")

    @property
    def identity(self) -> str:
        body = dict(self.payload)
        body.pop("identity", None)
        return _digest_bytes(_json_bytes(body))

    def to_json(self) -> bytes:
        payload = dict(self.payload)
        payload["identity"] = self.identity
        return _json_bytes(payload)


def _validate_finite_payload(value: Any, path: str = "payload") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise AnalysisError(f"Nonfinite numeric value at {path}")
    if isinstance(value, dict):
        for key, child in value.items():
            _validate_finite_payload(child, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            _validate_finite_payload(child, f"{path}[{index}]")


def _record_bindings(trajectory: VerifiedTrajectory, selection: AnalysisSelection | None = None) -> dict[str, Any]:
    metadata = trajectory.metadata
    result = {"final_graph_identity": metadata.final_graph_identity, "packing_identity": metadata.packing_identity, "force_field_source_identity": metadata.force_field_source_identity, "typed_graph_identity": metadata.typed_graph_identity, "charge_identity": metadata.charge_identity, "assignment_identity": metadata.assignment_identity, "workflow_identity": metadata.workflow_identity, "checkpoint_identity": metadata.checkpoint_identity, "atom_mapping_identity": trajectory.mapping_identity, "trajectory_identity": trajectory.identity, "source_hashes": dict(metadata.source_files)}
    if selection is not None:
        result["selection_identity"] = selection.identity
    return result


def make_analysis_record(trajectory: VerifiedTrajectory, observables: dict[str, dict[str, Any]], *, provenance: dict[str, Any] | None = None, evidence_references: tuple[str, ...] = (), selection: AnalysisSelection | None = None) -> AnalysisRecord:
    selection = selection or trajectory.selection
    if selection is not None:
        selection.validate(trajectory)
    metadata = trajectory.metadata.to_dict()
    payload = {"schema": SCHEMA, "implementation_version": "2", "identity": None, "bindings": _record_bindings(trajectory, selection), "trajectory": {"source_path": trajectory.source_path, "file_hashes": dict(trajectory.trajectory_file_hashes), "frame_count": trajectory.frame_count, "metadata": metadata}, "selection": selection.to_dict() if selection is not None else None, "observables": observables, "provenance": dict(provenance or {}), "evidence_references": list(evidence_references), "production_validated": False, "simulation_readiness": "not_established"}
    record = AnalysisRecord(payload)
    return AnalysisRecord({**payload, "identity": record.identity})


def write_analysis_record(record: AnalysisRecord, path: str | Path) -> Path:
    if type(record) is not AnalysisRecord:
        raise AnalysisError("Expected AnalysisRecord")
    target = Path(path)
    if target.exists():
        raise AnalysisError(f"Analysis record already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    raw = record.to_json()
    target.write_bytes(raw)
    return target


def load_analysis_record(path: str | Path, *, expected_identity: str | None = None, expected_bindings: dict[str, Any] | None = None, expected_trajectory: VerifiedTrajectory | None = None, expected_selection_identity: str | None = None, expected_selection: AnalysisSelection | None = None) -> AnalysisRecord:
    target = Path(path)
    try:
        payload = json.loads(target.read_text())
        stored = payload.get("identity")
        body = dict(payload); body.pop("identity", None)
        calculated = _digest_bytes(_json_bytes(body))
        if stored != calculated or (expected_identity is not None and stored != expected_identity):
            raise AnalysisError("Analysis record checksum/identity mismatch")
        record = AnalysisRecord(payload)
        trusted = expected_bindings
        if expected_trajectory is not None:
            trusted = _record_bindings(expected_trajectory, expected_trajectory.selection)
        if expected_selection is not None:
            expected_selection_identity = expected_selection.identity
        if expected_selection_identity is not None:
            trusted = dict(trusted or {})
            trusted["selection_identity"] = expected_selection_identity
        if trusted is not None:
            actual = payload.get("bindings", {})
            for name, expected in trusted.items():
                if actual.get(name) != expected:
                    raise AnalysisError(f"Analysis binding mismatch: {name}")
        return record
    except AnalysisError:
        raise
    except Exception as error:
        raise AnalysisError(f"Invalid analysis record: {error}") from error


def relocate_analysis_record(source: str | Path, destination: str | Path) -> Path:
    record = load_analysis_record(source)
    target = Path(destination)
    if target.exists():
        raise AnalysisError(f"Relocation target already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    loaded = load_analysis_record(target, expected_identity=record.identity)
    if loaded.identity != record.identity:
        raise AnalysisError("Relocated analysis identity changed")
    return target
