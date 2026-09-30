"""Bounded orchestration; manifests alone select published recovery boundaries."""

from copy import deepcopy
from importlib.util import find_spec
from pathlib import Path
from uuid import uuid4

from island.dynamics import (
    DynamicsSegment,
    DynamicsSegmentOptions,
    create_dynamics_checkpoint,
    initialize_velocities,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    save_dynamics_checkpoint,
)
from island.dynamics._checkpoint_data import canonical
from island.dynamics._checkpoint_data import checksum as payload_checksum
from island.evaluation import OpenMMSinglePointEvaluator
from island.exceptions import WorkflowError
from island.forcefields import AmberToolsParameterizationEngine
from island.minimization import minimize_geometry

from . import bundle, storage
from .config import WorkflowConfig

SCHEMA = "island_single_chain_manifest_v1"
BUNDLE_SCHEMA = "island_single_chain_bundle_v1"


def preflight(config, *, operation="start"):
    """Start checks construction/tools; resume only needs saved-input consumers."""
    if type(config) is not WorkflowConfig or operation not in (
        "start",
        "resume",
        "inspect",
        "prepared_start",
    ):
        raise WorkflowError(
            "Expected workflow config and start/prepared_start/resume/inspect operation"
        )
    names = ["rdkit"] if operation == "inspect" else ["openmm", "parmed", "rdkit"]
    if operation in ("start", "prepared_start"):
        names += ["scipy"]
    missing = [name for name in names if find_spec(name) is None]
    if missing:
        raise WorkflowError(
            f"Unavailable {operation} dependencies: {', '.join(missing)}"
        )
    if operation == "start":
        from island.forcefields.ambertools.engine import _discover

        _discover(config.amber_options())
    return {"operation": operation, "dependencies": names, "available": True}


def inspect_chain(config):
    """Build an owned local-template chain for inspecting the exact charge IDs."""
    preflight(config, operation="inspect")
    from island.builders import build_linear_polymer

    system = build_linear_polymer(
        config.psmiles,
        dp=config.dp,
        coordinate_method="local_templates",
        template_seed=config.template_seed,
        assembly_seed=config.assembly_seed,
        stereo_seed=config.stereo_seed,
        tacticity=config.tacticity,
        atactic_fraction=config.atactic_fraction,
    )
    if system.number_of_sites > 100:
        raise WorkflowError("AmberTools preparation is limited to 100 sites")
    return system


def _charges(config, system):
    if config.charge_method == "provided":
        from island.forcefields import ProvidedChargeEngine

        ProvidedChargeEngine().assign(
            system,
            dict(config.provided_charges),
            source="workflow explicit charges",
            tolerance=config.charge_tolerance,
        )


def _new_manifest(config, evidence):
    return {
        "schema": SCHEMA,
        "run_id": uuid4().hex,
        "config": config.to_dict(),
        "evidence": evidence,
        "status": "starting",
        "stages": [],
        "files": {},
        "bundle": None,
        "segments": [],
        "checkpoint": None,
        "accepted_step": 0,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def _manifest(root, manifest):
    storage.publish(
        Path(root) / "manifest.json",
        storage.json_bytes(
            {
                "payload": manifest,
                "sha256": payload_checksum(manifest),
            }
        ),
        replace=True,
    )


def _put(root, manifest, label, raw):
    # UUID names allow explicit resume despite unpublished orphan files.
    suffix = "dat" if label == "source" else "json"
    name = f"{label}-{uuid4().hex}.{suffix}"
    if (
        sum(v["bytes"] for v in manifest["files"].values()) + len(raw)
        > manifest["config"]["max_artifact_bytes"]
    ):
        raise WorkflowError("Published artifact byte budget exhausted")
    storage.publish(storage.child(root, name), raw)
    manifest["files"][name] = {"sha256": storage.checksum(raw), "bytes": len(raw)}
    return name


def _outcome(root, manifest, stage, status, **details):
    manifest["stages"].append({"stage": stage, "status": status, **details})
    _manifest(root, manifest)


def _failure(root, manifest, stage, error):
    manifest["status"] = "stage_failed"
    status = (
        "unavailable"
        if "Unavailable" in type(error).__name__ or "Unavailable" in str(error)
        else "failed"
    )
    _outcome(
        root,
        manifest,
        stage,
        status,
        error_type=type(error).__name__,
        message=str(error),
        details={
            name: getattr(error, name)
            for name in (
                "stage",
                "artifact_dir",
                "command",
                "stdout",
                "stderr",
                "returncode",
                "evaluations",
            )
            if hasattr(error, name)
        },
    )
    return deepcopy(manifest)


def _prepare_bundle(root, manifest, config, system, preparation, artifact_directory):
    preparation.validate_integrity(system)
    if system.number_of_sites > 100:
        raise WorkflowError("AmberTools preparation is limited to 100 sites")
    _charges(config, system)
    r = preparation.record
    if (
        r["requested_force_field"] != config.force_field
        or r["charge_method"] != config.charge_method
        or r["charge_validation_tolerance_e"] != config.charge_tolerance
    ):
        raise WorkflowError("Configuration contradicts historical preparation")
    if config.provided_charges is not None:
        for site, charge in config.provided_charges.items():
            actual = preparation.imported_result.charge_result.assignments[site].charge
            if abs(actual - charge) > r["provided_charge_tolerance_e"]:
                raise WorkflowError("Provided charge differs from prepared value")
    required = {
        **r["artifact_sha256"],
        "input.mol2": r["input_mol2_sha256"],
        "lineage.json": r["input_lineage_sha256"],
    }
    if r["charge_outcome"]["qm_run"]:
        required["sqm.out"] = r["charge_outcome"]["sqm_out_sha256"]
    artifacts = {}
    # Copy all retained regular artifacts, including actual SQM output/commands.
    directory = Path(artifact_directory)
    for name, expected in required.items():
        if storage.checksum(storage.child(directory, name).read_bytes()) != expected:
            raise WorkflowError(f"Historical artifact checksum mismatch: {name}")
    for path in sorted(directory.iterdir()):
        if path.is_file() and not path.is_symlink():
            raw = path.read_bytes()
            artifacts[path.name] = _put(root, manifest, "source", raw)
    record = {
        "record": dict(r),
        "record_signature": preparation.record_signature,
        "import_source": preparation.imported_result.source,
        "import_provenance": dict(preparation.imported_result.provenance),
    }
    # Verify operational reconstruction now, before accepting a durable bundle.
    bundle.preparation_from(
        system, record, storage.child(root, artifacts["result.prmtop"])
    )
    return artifacts, record


def _finish_setup(
    root, manifest, config, system, preparation, artifacts, preparation_record
):
    stage = "minimization"
    try:
        evaluator = OpenMMSinglePointEvaluator(
            system, preparation.imported_result, platform="Reference"
        )
        with evaluator.open_session() as session:
            minimum = minimize_geometry(system, session, config.minimization)
        minimum.validate_integrity()
        report = _put(
            root,
            manifest,
            "minimization",
            storage.json_bytes(storage.encode(bundle.record(minimum))),
        )
        _outcome(
            root,
            manifest,
            stage,
            "passed"
            if minimum.converged and minimum.final_evaluation_verified
            else "failed",
            record=report,
            evaluations=minimum.evaluations,
            reason=minimum.termination_reason,
        )
        if not minimum.converged or not minimum.final_evaluation_verified:
            manifest["status"] = "stage_failed"
            _manifest(root, manifest)
            return False
        starting = minimum.to_system(system)
        stage = "initialization"
        initialized = initialize_velocities(
            starting,
            temperature_kelvin=config.temperature_kelvin,
            seed=config.velocity_seed,
        )
        payload = {
            "schema": BUNDLE_SCHEMA,
            "config": config.to_dict(),
            "input_system": bundle.system_data(system),
            "starting_system": bundle.system_data(starting),
            "preparation": preparation_record,
            "artifacts": artifacts,
            "minimum": bundle.record(minimum),
            "initialization": bundle.record(initialized),
        }
        manifest["bundle"] = _put(
            root, manifest, "bundle", storage.json_bytes(storage.encode(payload))
        )
        manifest["status"] = "ready"
        _outcome(
            root,
            manifest,
            stage,
            "passed",
            initialization_fingerprint=initialized.initialization_fingerprint,
        )
        return True
    except Exception as error:  # noqa: BLE001 -- preserve stage failure diagnostics
        _failure(root, manifest, stage, error)
        return False


def start_workflow(config, *, segments=1):
    """Create a new run; execute at most the explicitly requested segment count."""
    return _start(config, segments, None)


def start_prepared_workflow(
    config, system, preparation, artifact_directory, *, segments=1, evidence
):
    """Explicit archived/synthetic acceptance path; never substitutes for live start."""
    from island.core import MolecularSystem
    from island.forcefields.ambertools.models import AmberToolsPreparationResult

    if not isinstance(system, MolecularSystem) or not isinstance(
        preparation, AmberToolsPreparationResult
    ):
        raise WorkflowError(
            "Prepared workflow requires MolecularSystem and AmberToolsPreparationResult"
        )
    if evidence not in ("archived_parameters", "synthetic_software_test"):
        raise WorkflowError("Prepared input requires explicit evidence classification")
    return _start(
        config,
        segments,
        (system.copy(), preparation, Path(artifact_directory), evidence),
    )


def _segment_count(segments):
    if type(segments) is not int or segments <= 0:
        raise WorkflowError("segments must be a positive integer")


def _start(config, segments, prepared):
    if type(config) is not WorkflowConfig:
        raise WorkflowError("Expected WorkflowConfig")
    _segment_count(segments)
    root = Path(config.output_directory)
    try:
        root.mkdir(parents=True, exist_ok=False)
    except OSError as error:
        raise WorkflowError(f"New workflow requires a new directory: {root}") from error
    with storage.writer(root):
        manifest = _new_manifest(
            config, "live_parameterization" if prepared is None else prepared[3]
        )
        _manifest(root, manifest)
        stage = "preflight"
        try:
            preflight(
                config, operation="start" if prepared is None else "prepared_start"
            )
            _outcome(root, manifest, stage, "passed")
            stage = "construction"
            if prepared is None:
                system = inspect_chain(config)
                _charges(config, system)
                input_record = _put(
                    root,
                    manifest,
                    "input-system",
                    storage.json_bytes(storage.encode(bundle.system_data(system))),
                )
                _outcome(
                    root,
                    manifest,
                    stage,
                    "passed",
                    sites=system.number_of_sites,
                    record=input_record,
                )
                stage = "parameterization"
                preparation = AmberToolsParameterizationEngine().parameterize(
                    system,
                    config.amber_options(root / "preparation-work"),
                )
                directory = Path(preparation.record["artifact_dir"])
            else:
                system, preparation, directory, _ = prepared
                _outcome(
                    root, manifest, stage, "supplied", sites=system.number_of_sites
                )
                stage = "parameterization"
            artifacts, record = _prepare_bundle(
                root, manifest, config, system, preparation, directory
            )
            _outcome(
                root,
                manifest,
                stage,
                "passed" if prepared is None else "supplied",
                preparation_signature=preparation.record_signature,
                imported_signature=preparation.imported_result.result_signature,
            )
        except Exception as error:  # noqa: BLE001 -- preserve stage failure diagnostics
            return _failure(root, manifest, stage, error)
        if not _finish_setup(
            root, manifest, config, system, preparation, artifacts, record
        ):
            return deepcopy(manifest)
        return _advance(root, manifest, segments)


def workflow_status(directory):
    """Read only the published manifest; verify every referenced file checksum."""
    root = Path(directory)
    try:
        envelope = storage.read_json(root / "manifest.json")
        if set(envelope) != {"payload", "sha256"}:
            raise ValueError("Invalid manifest envelope")
        m = envelope["payload"]
        expected = {
            "schema",
            "run_id",
            "config",
            "evidence",
            "status",
            "stages",
            "files",
            "bundle",
            "segments",
            "checkpoint",
            "accepted_step",
            "production_validated",
            "simulation_readiness",
        }
        if (
            set(m) != expected
            or m["schema"] != SCHEMA
            or payload_checksum(m) != envelope["sha256"]
        ):
            raise ValueError("Manifest schema/checksum mismatch")
        config = WorkflowConfig.from_dict(m["config"])
        if (
            m["production_validated"] is not False
            or m["simulation_readiness"] != "not_established"
            or type(m["accepted_step"]) is not int
            or not 0 <= m["accepted_step"] <= config.total_steps
            or m["status"]
            not in (
                "starting",
                "ready",
                "paused",
                "completed",
                "budget_exhausted",
                "stage_failed",
            )
        ):
            raise ValueError("Invalid workflow status")
        if m["evidence"] not in (
            "live_parameterization",
            "archived_parameters",
            "synthetic_software_test",
        ):
            raise ValueError("Unknown evidence classification")
        if (
            type(m["stages"]) is not list
            or type(m["segments"]) is not list
            or len(m["segments"]) > config.max_segments
        ):
            raise ValueError("Invalid stage/segment records")
        if len(set(m["segments"])) != len(m["segments"]):
            raise ValueError("Duplicate segment reference")
        for stage in m["stages"]:
            if (
                type(stage) is not dict
                or stage.get("stage")
                not in (
                    "preflight",
                    "construction",
                    "parameterization",
                    "minimization",
                    "initialization",
                    "dynamics",
                )
                or stage.get("status")
                not in (
                    "passed",
                    "supplied",
                    "failed",
                    "unavailable",
                    "paused",
                    "completed",
                    "budget_exhausted",
                    "stage_failed",
                )
            ):
                raise ValueError("Malformed stage outcome")
            if "record" in stage and stage["record"] not in m["files"]:
                raise ValueError("Unlisted stage record")
            if "evaluations" in stage and (
                type(stage["evaluations"]) is not int or stage["evaluations"] < 0
            ):
                raise ValueError("Malformed stage call count")
        if sum(s["stage"] == "dynamics" for s in m["stages"]) > config.max_segments:
            raise ValueError("Segment attempt budget exceeded")
        size = 0
        for name, row in m["files"].items():
            if (
                set(row) != {"sha256", "bytes"}
                or type(row["bytes"]) is not int
                or row["bytes"] < 0
            ):
                raise ValueError("Malformed artifact descriptor")
            path = storage.child(root, name)
            if path.stat().st_size != row["bytes"]:
                raise ValueError(f"Artifact size mismatch: {name}")
            size += row["bytes"]
            if (
                size > config.max_artifact_bytes
                or storage.checksum(path.read_bytes()) != row["sha256"]
            ):
                raise ValueError(f"Artifact checksum/storage limit mismatch: {name}")
        for name in (m["bundle"], m["checkpoint"], *m["segments"]):
            if name is not None and name not in m["files"]:
                raise ValueError("Unlisted operational artifact")
        previous = None
        accepted = None
        for name in m["segments"]:
            segment = _load_segment(root, name)
            if previous is not None and segment.frames[0] != previous:
                raise ValueError("Conflicting adjacent segment boundary")
            parent = segment.payload["lineage"][-1]["parent_checksum"]
            if parent != (None if accepted is None else accepted.content_checksum):
                raise ValueError("Segment lineage does not extend published boundary")
            if accepted is None and segment.frames[0].step != 0:
                raise ValueError("Workflow must start at step zero")
            previous = segment.final_state
            if segment.payload["diagnostic"][
                "final_verified"
            ] and segment.termination_reason in (
                "completed",
                "maximum_frames",
                "maximum_evaluations",
            ):
                accepted = create_dynamics_checkpoint(segment)
        if m["checkpoint"] is not None:
            saved = load_dynamics_checkpoint(storage.child(root, m["checkpoint"]))
            if saved != accepted or saved.absolute_step != m["accepted_step"]:
                raise ValueError("Manifest checkpoint/progress contradicts segments")
        elif m["accepted_step"] != 0 or accepted is not None:
            raise ValueError("Progress without published checkpoint")
        if (m["status"] == "completed") != (m["accepted_step"] == config.total_steps):
            raise ValueError("False workflow completion")
        return m
    except WorkflowError:
        raise
    except Exception as error:
        raise WorkflowError(f"Invalid workflow manifest: {error}") from error


def _load_segment(root, name):
    envelope = storage.read_json(storage.child(root, name))
    if set(envelope) != {"payload", "sha256"}:
        raise WorkflowError("Invalid segment envelope")
    return DynamicsSegment(canonical(envelope["payload"]), envelope["sha256"])


def _load_bundle(root, manifest):
    if manifest["bundle"] is None:
        raise WorkflowError(
            "No prepared bundle was published; start a new run after resolving the failed stage"
        )
    try:
        p = storage.decode(storage.read_json(storage.child(root, manifest["bundle"])))
        if (
            set(p)
            != {
                "schema",
                "config",
                "input_system",
                "starting_system",
                "preparation",
                "artifacts",
                "minimum",
                "initialization",
            }
            or p["schema"] != BUNDLE_SCHEMA
        ):
            raise ValueError("Invalid bundle schema")
        if p["config"] != manifest["config"]:
            raise ValueError("Configuration differs from immutable prepared bundle")
        original = bundle.system_from(p["input_system"])
        for logical, path in p["artifacts"].items():
            if Path(logical).name != logical or path not in manifest["files"]:
                raise ValueError("Unlisted source artifact")
        r = p["preparation"]["record"]
        for name, sha in {
            **r["artifact_sha256"],
            "input.mol2": r["input_mol2_sha256"],
            "lineage.json": r["input_lineage_sha256"],
        }.items():
            if manifest["files"][p["artifacts"][name]]["sha256"] != sha:
                raise ValueError("Artifact contradicts historical signature")
        if (
            r["charge_outcome"]["qm_run"]
            and manifest["files"][p["artifacts"]["sqm.out"]]["sha256"]
            != r["charge_outcome"]["sqm_out_sha256"]
        ):
            raise ValueError("SQM artifact contradicts historical signature")
        prepared = bundle.preparation_from(
            original,
            p["preparation"],
            storage.child(root, p["artifacts"]["result.prmtop"]),
        )
        minimum = bundle.minimum_from(p["minimum"])
        starting = minimum.to_system(original)
        if bundle.system_data(starting) != p["starting_system"]:
            raise ValueError("Starting system differs from verified minimum")
        initialization = bundle.initialization_from(p["initialization"])
        from island.minimization.models import system_identity

        if initialization.system_fingerprint != system_identity(starting):
            raise ValueError("Initialization system mismatch")
        config = WorkflowConfig.from_dict(manifest["config"])
        if (
            initialization.temperature_kelvin != config.temperature_kelvin
            or initialization.velocity_seed != config.velocity_seed
            or minimum.options != config.minimization
            or r["requested_force_field"] != config.force_field
            or r["charge_method"] != config.charge_method
            or r["charge_validation_tolerance_e"] != config.charge_tolerance
        ):
            raise ValueError("Bundle/config identity mismatch")
        return starting, prepared, initialization
    except Exception as error:
        raise WorkflowError(f"Invalid prepared bundle: {error}") from error


def read_workflow_frames(directory):
    """Owned retained endpoint frames, deduplicating identical shared boundaries."""
    m = workflow_status(directory)
    frames = []
    for name in m["segments"]:
        for frame in _load_segment(directory, name).frames:
            if frames and frame.step == frames[-1].step:
                if frame != frames[-1]:
                    raise WorkflowError("Conflicting duplicate frame")
                continue
            if frames and frame.step < frames[-1].step:
                raise WorkflowError("Out-of-order frame")
            frames.append(frame)
    return tuple(frames)


def resume_workflow(directory, *, segments=1):
    """Explicit continuation from the manifest; no build, QM, minimization or reseed."""
    _segment_count(segments)
    root = Path(directory)
    with storage.writer(root):
        manifest = workflow_status(root)
        if manifest["status"] == "completed":
            return manifest
        config = WorkflowConfig.from_dict(manifest["config"])
        preflight(config, operation="resume")
        return _advance(root, manifest, segments)


def _advance(root, manifest, count):
    config = WorkflowConfig.from_dict(manifest["config"])
    system, prepared, initialization = _load_bundle(root, manifest)
    evaluator = OpenMMSinglePointEvaluator(
        system, prepared.imported_result, platform="Reference"
    )
    checkpoint = (
        None
        if manifest["checkpoint"] is None
        else load_dynamics_checkpoint(storage.child(root, manifest["checkpoint"]))
    )
    if checkpoint is not None:
        physical = checkpoint.payload["physical"]
        for name in (
            "temperature_kelvin",
            "friction_per_ps",
            "timestep_fs",
            "thermostat_seed",
        ):
            if physical[name] != getattr(config, name):
                raise WorkflowError(
                    f"Checkpoint/config physical setting differs: {name}"
                )
    for _ in range(count):
        if manifest["accepted_step"] == config.total_steps:
            break
        if (
            sum(s["stage"] == "dynamics" for s in manifest["stages"])
            >= config.max_segments
        ):
            manifest["status"] = "budget_exhausted"
            _manifest(root, manifest)
            break
        steps = min(
            config.segment_steps, config.total_steps - manifest["accepted_step"]
        )
        try:
            with evaluator.open_session() as session:
                if checkpoint is None:
                    result = run_dynamics_segment(
                        system,
                        session,
                        initialization.velocities,
                        config.langevin(steps),
                    )
                else:
                    result = resume_dynamics(
                        checkpoint,
                        system,
                        session,
                        DynamicsSegmentOptions(
                            steps,
                            config.max_evaluations_per_segment,
                            config.max_frames_per_segment,
                            config.recording_interval,
                        ),
                    )
        except Exception as error:  # noqa: BLE001 -- preserve stage outcome, no retry
            return _failure(root, manifest, "dynamics", error)
        # Publish all immutable data before exposing progress. Any exception
        # leaves the previous manifest/boundary authoritative; orphan files ignored.
        candidate = deepcopy(manifest)
        name = _put(
            root,
            candidate,
            "segment",
            storage.json_bytes(
                {"payload": result.payload, "sha256": result.content_checksum}
            ),
        )
        segment_name = name
        eligible = result.payload["diagnostic"][
            "final_verified"
        ] and result.termination_reason in (
            "completed",
            "maximum_frames",
            "maximum_evaluations",
        )
        if eligible:
            candidate["segments"].append(segment_name)
            next_checkpoint = create_dynamics_checkpoint(result)
            name = f"checkpoint-{uuid4().hex}.json"
            raw = storage.json_bytes(
                {
                    "payload": next_checkpoint.payload,
                    "sha256": next_checkpoint.content_checksum,
                }
            )
            if (
                sum(row["bytes"] for row in candidate["files"].values()) + len(raw)
                > config.max_artifact_bytes
            ):
                raise WorkflowError(
                    "Published artifact byte budget exhausted; prior boundary retained"
                )
            save_dynamics_checkpoint(next_checkpoint, storage.child(root, name))
            storage.sync_directory(root)
            candidate["files"][name] = {
                "sha256": storage.checksum(storage.child(root, name).read_bytes()),
                "bytes": storage.child(root, name).stat().st_size,
            }
            candidate["checkpoint"] = name
            candidate["accepted_step"] = next_checkpoint.absolute_step
        candidate["status"] = (
            "completed"
            if candidate["accepted_step"] == config.total_steps
            else "stage_failed"
            if not eligible
            else "paused"
            if result.completed
            else "budget_exhausted"
        )
        candidate["stages"].append(
            {
                "stage": "dynamics",
                "status": candidate["status"],
                "reason": result.termination_reason,
                "evaluations": result.evaluations,
                "record": segment_name,
            }
        )
        _manifest(root, candidate)
        manifest = candidate
        if not eligible or not result.completed:
            break
        checkpoint = next_checkpoint
    return deepcopy(manifest)
