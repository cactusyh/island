"""Read one pinned workflow publication; never execute scientific engines."""

import json
import math
import tempfile
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

from island.exceptions import AnalysisError, IslandError

from .geometry import GeometryOptions, geometry_metrics


@dataclass(frozen=True)
class AnalysisOptions:
    geometry: GeometryOptions = field(default_factory=GeometryOptions)
    start_step: int | None = None
    end_step: int | None = None
    start_time_ps: float | None = None
    end_time_ps: float | None = None
    stride: int = 1

    def __post_init__(self):
        if type(self.geometry) is not GeometryOptions:
            raise AnalysisError("Expected GeometryOptions")
        if type(self.stride) is not int or self.stride < 1:
            raise AnalysisError("stride must be a positive retained-sample integer")
        for name in ("start_step", "end_step"):
            x = getattr(self, name)
            if x is not None and (type(x) is not int or x < 0):
                raise AnalysisError(f"Invalid {name}")
        for name in ("start_time_ps", "end_time_ps"):
            x = getattr(self, name)
            if x is None:
                continue
            try:
                if type(x) not in (int, float) or x < 0 or not math.isfinite(x):
                    raise AnalysisError(f"Invalid {name}")
            except (OverflowError, TypeError, ValueError) as error:
                raise AnalysisError(f"Unrepresentable {name}") from error
        for lo, hi in (
            (self.start_step, self.end_step),
            (self.start_time_ps, self.end_time_ps),
        ):
            if lo is not None and hi is not None and lo > hi:
                raise AnalysisError("Reversed analysis window")


@dataclass(frozen=True)
class AnalysisReport:
    """Owned JSON data; payload returns a fresh copy, not mutable shared state."""

    _json: str

    def validate_integrity(self):
        """Check internal v1 consistency without requiring source files."""
        from .integrity import validate_report

        validate_report(self._json)

    @property
    def payload(self):
        from .integrity import validate_report

        return validate_report(self._json)


def resolve_endpoints(system):
    """Validate existing linear-repeat provenance; never invent endpoints."""
    from island.conformations.random_walk import _validate_linear_polymer

    try:
        layout = _validate_linear_polymer(system)
        ids = (layout.head_site_id, layout.tail_site_id)
        for s in ids:
            site = system.topology.sites[s]
            if (
                type(s) is not int
                or site.atomic_number <= 1
                or type(site.metadata.get("source_repeat_atom_index")) is not int
                or site.metadata.get("generated_hydrogen", False)
            ):
                raise ValueError("Endpoints lack original heavy-atom repeat provenance")
        return ids, "builder_head_tail", None
    except (IslandError, ValueError, TypeError, KeyError, AttributeError) as error:
        return None, "unavailable", f"endpoints_unavailable: {error}"


def _pin(source, target):
    from island.dynamics._checkpoint_data import checksum, strict_load
    from island.workflows import storage

    raw = (source / "manifest.json").read_bytes()
    envelope = strict_load(raw.decode())
    if (
        set(envelope) != {"payload", "sha256"}
        or checksum(envelope["payload"]) != envelope["sha256"]
    ):
        raise AnalysisError("Manifest envelope/checksum mismatch")
    m = envelope["payload"]
    # Bound the private copy by the declared storage budget before reading files.
    limit = m["config"]["max_artifact_bytes"]
    total = 0
    for name, descriptor in m["files"].items():
        path = storage.child(source, name)
        size = descriptor["bytes"]
        if type(size) is not int or size < 0 or path.stat().st_size != size:
            raise AnalysisError("Invalid snapshot artifact size")
        total += size
        if total > limit:
            raise AnalysisError("Snapshot exceeds declared artifact budget")
        data = path.read_bytes()
        if len(data) != size or storage.checksum(data) != descriptor["sha256"]:
            raise AnalysisError(f"Snapshot artifact checksum mismatch: {name}")
        destination = storage.child(target, name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
    (target / "manifest.json").write_bytes(raw)
    return envelope["sha256"], storage.checksum(raw)


def analyze_workflow(directory, options=None):
    """Analyze inclusive windows, then every stride-th retained frame.

    Deep saved-bundle validation requires ParmEd (and RDKit for assigned stereo).
    No OpenMM Context or other scientific calculation is created.
    """
    options = AnalysisOptions() if options is None else options
    if type(options) is not AnalysisOptions:
        raise AnalysisError("Expected AnalysisOptions")
    from island.workflows import read_workflow_frames, workflow_status
    from island.workflows.chain import _load_bundle, _load_segment

    try:
        source = Path(directory).resolve()
        with tempfile.TemporaryDirectory(prefix="island-analysis-") as temporary:
            root = Path(temporary)
            manifest_checksum, file_checksum = _pin(source, root)
            manifest = workflow_status(root)
            frames = read_workflow_frames(root)
            if not frames:
                raise AnalysisError("No accepted retained frames available")
            system, _, _ = _load_bundle(root, manifest)
            sites = system.topology.sites
            masses = {s: a.mass for s, a in sites.items()}
            numbers = {s: a.atomic_number for s, a in sites.items()}
            geometry = options.geometry
            convention, reason = "explicit_stable_ids", None
            if geometry.endpoints is None:
                endpoints, convention, reason = resolve_endpoints(system)
                geometry = replace(geometry, endpoints=endpoints)
            segments = [_load_segment(root, name) for name in manifest["segments"]]
            frame_sources = {}
            for name, segment in zip(manifest["segments"], segments, strict=True):
                for frame in segment.frames:
                    frame_sources.setdefault(frame.step, name)
            selected = [
                f
                for f in frames
                if (
                    (options.start_step is None or f.step >= options.start_step)
                    and (options.end_step is None or f.step <= options.end_step)
                    and (
                        options.start_time_ps is None
                        or f.time_ps >= options.start_time_ps
                    )
                    and (
                        options.end_time_ps is None or f.time_ps <= options.end_time_ps
                    )
                )
            ][:: options.stride]
            if not selected:
                raise AnalysisError("No retained frames in requested selection")
            rows = []
            for frame in selected:
                g = geometry_metrics(
                    frame.coordinates, masses, atomic_numbers=numbers, options=geometry
                )
                g = replace(
                    g,
                    endpoint_convention=convention,
                    diagnostics=g.diagnostics + ((reason,) if reason else ()),
                )
                rows.append(
                    {
                        "step": frame.step,
                        "time_ps": frame.time_ps,
                        "coordinate_fingerprint": frame.coordinate_fingerprint,
                        "source_segment": frame_sources[frame.step],
                        "evaluation_fingerprint": frame.evaluation.evaluation_fingerprint,
                        "potential_energy_kj_mol": frame.potential_energy,
                        "kinetic_energy_kj_mol": frame.kinetic_energy,
                        "total_energy_kj_mol": frame.total_energy,
                        "temperature_kelvin": frame.instantaneous_temperature_kelvin,
                        "geometry": g.to_dict(),
                    }
                )
            payload = {
                "schema": "island.trajectory-analysis.v1",
                "implementation_version": "1",
                "source": {
                    "directory": str(source),
                    "manifest_checksum": manifest_checksum,
                    "manifest_file_sha256": file_checksum,
                    "status": manifest["status"],
                    "accepted_step": manifest["accepted_step"],
                    "run_id": manifest["run_id"],
                    "evidence": manifest["evidence"],
                    "production_validated": manifest["production_validated"],
                    "simulation_readiness": manifest["simulation_readiness"],
                    "model_fingerprint": frames[0].evaluation.model_fingerprint,
                    "parameter_fingerprint": frames[0].evaluation.parameter_fingerprint,
                    "backend": {
                        "name": frames[0].evaluation.backend_name,
                        "version": frames[0].evaluation.backend_version,
                        "platform": frames[0].evaluation.platform,
                        "settings": dict(frames[0].evaluation.settings),
                    },
                    "artifacts": manifest["files"],
                    "checkpoint": manifest["checkpoint"],
                    "segments": [
                        {"file": n, "checksum": s.content_checksum}
                        for n, s in zip(manifest["segments"], segments, strict=True)
                    ],
                    "origin": segments[0].payload["origin"],
                },
                "options": asdict(options),
                "frames": rows,
                "sample_count": len(rows),
                "sampling": "inclusive windows then stride over retained samples; no interpolation",
                "boundary_policy": "retain earlier accepted frame; original records unchanged",
                "conventions": {
                    "weighting": geometry.weighting,
                    "center": "sum(w*r)/sum(w)",
                    "tensor": "sum(w*(r-center)*(r-center)^T)/sum(w)",
                    "Rg": "sqrt(trace(G))",
                    "eigenvalues": "ascending angstrom^2",
                    "kappa_squared": "1.5*sum(lambda^2)/sum(lambda)^2-0.5",
                    "temperature_dof": "3N; translation and rotation included",
                    "missing_values": "JSON null; CSV empty field; see diagnostics",
                },
                "production_validated": False,
                "simulation_readiness": "not_established",
                "interpretation": "retained samples only; no equilibrium or independence claim",
            }
            result = AnalysisReport(
                json.dumps(payload, sort_keys=True, allow_nan=False)
            )
            result.validate_integrity()
            return result
    except AnalysisError:
        raise
    except Exception as error:
        raise AnalysisError(f"Workflow analysis failed: {error}") from error
