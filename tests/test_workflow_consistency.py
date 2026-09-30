"""Cross-record workflow consistency, including controlled restart offsets."""

from dataclasses import replace

import pytest
from test_workflow import config, start
from test_workflow import prepared as prepared_fixture

prepared = prepared_fixture

from island.dynamics import (
    create_dynamics_checkpoint,
    initialize_velocities,
    run_dynamics_segment,
    save_dynamics_checkpoint,
)
from island.evaluation import OpenMMEvaluationSession, OpenMMSinglePointEvaluator
from island.exceptions import WorkflowError
from island.workflows import (
    read_workflow_frames,
    resume_workflow,
    storage,
    workflow_status,
)
from island.workflows.chain import _load_bundle, _manifest, _put


@pytest.mark.parametrize("offset", [5e-9, -5e-9])
def test_tolerated_startup_offset_remains_readable(
    tmp_path, prepared, monkeypatch, offset
):
    root = tmp_path / "run"
    first = start(root, prepared)
    old_bytes = (root / first["segments"][0]).read_bytes()
    original = OpenMMEvaluationSession.evaluate_fresh
    calls = 0

    def shifted(self, *a, **kw):
        nonlocal calls
        calls += 1
        r = original(self, *a, **kw)
        if calls == 1:
            components = dict(r.energy_components)
            key = next(iter(components))
            components[key] += offset
            return replace(
                r,
                potential_energy=r.potential_energy + offset,
                energy_components=components,
            )
        return r

    with monkeypatch.context() as patch:
        patch.setattr(OpenMMEvaluationSession, "evaluate_fresh", shifted)
        result = resume_workflow(root)
    assert result["accepted_step"] == 4
    assert workflow_status(root)["accepted_step"] == 4
    frames = read_workflow_frames(root)
    assert [f.step for f in frames] == list(range(5))
    from island.workflows.chain import _load_segment

    assert frames[2] == _load_segment(root, first["segments"][0]).final_state
    assert (
        frames[2].evaluation
        != _load_segment(root, result["segments"][1]).frames[0].evaluation
    )
    assert (root / first["segments"][0]).read_bytes() == old_bytes
    assert resume_workflow(root)["status"] == "completed"


def alternative(root, manifest, mode):
    system, preparation, initialization = _load_bundle(root, manifest)
    velocities = {
        s: tuple(2 * x for x in v) for s, v in initialization.velocities.items()
    }
    if mode == "seed":
        velocities = initialize_velocities(
            system, temperature_kelvin=300, seed=99
        ).velocities
    with OpenMMSinglePointEvaluator(
        system, preparation.imported_result
    ).open_session() as session:
        return run_dynamics_segment(
            system, session, velocities, config(root).langevin(2)
        )


def publish_alternative(root, manifest, result):
    name = _put(
        root,
        manifest,
        "inconsistent-segment",
        storage.json_bytes(
            {
                "payload": result.payload,
                "sha256": result.content_checksum,
            }
        ),
    )
    cp = create_dynamics_checkpoint(result)
    cpname = "inconsistent-checkpoint.json"
    save_dynamics_checkpoint(cp, root / cpname)
    raw = (root / cpname).read_bytes()
    manifest["files"][cpname] = {"sha256": storage.checksum(raw), "bytes": len(raw)}
    manifest["segments"] = [name]
    manifest["checkpoint"] = cpname
    manifest["stages"][-1]["record"] = name
    _manifest(root, manifest)


@pytest.mark.parametrize("mode", ["twice", "seed"])
def test_valid_but_different_initialization_rejected(tmp_path, prepared, mode):
    root = tmp_path / "run"
    m = start(root, prepared)
    other = alternative(root, m, mode)
    publish_alternative(root, m, other)
    for read in (workflow_status, read_workflow_frames, resume_workflow):
        with pytest.raises(WorkflowError):
            read(root)


@pytest.mark.parametrize("completed", [False, True])
def test_config_disagreement_rejected_by_every_reader(tmp_path, prepared, completed):
    root = tmp_path / "run"
    m = start(root, prepared, segments=3 if completed else 1)
    m["config"]["temperature_kelvin"] = 777
    _manifest(root, m)
    for read in (workflow_status, read_workflow_frames, resume_workflow):
        with pytest.raises(WorkflowError):
            read(root)


@pytest.mark.parametrize("kind", ["energy", "force"])
def test_outside_tolerance_retains_saved_boundary(
    tmp_path, prepared, monkeypatch, kind
):
    root = tmp_path / "run"
    first = start(root, prepared)
    original = OpenMMEvaluationSession.evaluate_fresh
    calls = 0

    def shifted(self, *a, **kw):
        nonlocal calls
        calls += 1
        r = original(self, *a, **kw)
        if kind == "energy":
            components = dict(r.energy_components)
            components[next(iter(components))] += 0.01
            return replace(
                r,
                potential_energy=r.potential_energy + 0.01,
                energy_components=components,
            )
        return replace(
            r, forces={s: tuple(v + 0.01 for v in f) for s, f in r.forces.items()}
        )

    with monkeypatch.context() as patch:
        patch.setattr(OpenMMEvaluationSession, "evaluate_fresh", shifted)
        result = resume_workflow(root)
    assert calls == 1
    assert result["status"] == "stage_failed"
    assert result["checkpoint"] == first["checkpoint"]
    assert workflow_status(root)["accepted_step"] == 2
    assert len(read_workflow_frames(root)) == 3


def test_invalid_publication_is_transactional(tmp_path, prepared, monkeypatch):
    root = tmp_path / "run"
    first = start(root, prepared)
    other = alternative(root, first, "twice")
    before = (root / "manifest.json").read_bytes()
    checkpoint = (root / first["checkpoint"]).read_bytes()
    with monkeypatch.context() as patch:
        patch.setattr("island.workflows.chain.resume_dynamics", lambda *a, **kw: other)
        with pytest.raises(WorkflowError):
            resume_workflow(root)
    assert (root / "manifest.json").read_bytes() == before
    assert (root / first["checkpoint"]).read_bytes() == checkpoint
    assert workflow_status(root)["accepted_step"] == 2
    assert resume_workflow(root)["accepted_step"] == 4


@pytest.mark.parametrize("field", ["coordinates", "velocities", "model"])
def test_changed_boundary_state_rejected(tmp_path, prepared, field):
    from island.core.coordinate_provenance import coordinate_hash
    from island.dynamics.models import kinetic_energy, velocity_hash
    from island.dynamics.thermal import GAS_CONSTANT
    from island.workflows.chain import _load_segment
    from island.workflows.consistency import compatible_boundary

    root = tmp_path / "run"
    m = start(root, prepared)
    frame = _load_segment(root, m["segments"][0]).final_state
    system, prep, init = _load_bundle(root, m)
    if field == "model":
        changed = replace(
            frame, evaluation=replace(frame.evaluation, model_fingerprint="0" * 64)
        )
        changed.validate_integrity()
        with pytest.raises(WorkflowError, match="model"):
            compatible_boundary(frame, changed)
        return
    vectors = dict(getattr(frame, field))
    site = next(iter(vectors))
    vectors[site] = tuple(v + 1e-5 for v in vectors[site])
    if field == "coordinates":
        evaluation = OpenMMSinglePointEvaluator(system, prep.imported_result).evaluate(
            vectors
        )
        total = evaluation.potential_energy + frame.kinetic_energy
        changed = replace(
            frame,
            coordinates=vectors,
            coordinate_fingerprint=coordinate_hash(vectors),
            evaluation=evaluation,
            total_energy=total,
            energy_deviation=frame.energy_deviation + total - frame.total_energy,
        )
    else:
        kinetic = kinetic_energy(init.masses, vectors)
        total = frame.potential_energy + kinetic
        changed = replace(
            frame,
            velocities=vectors,
            velocity_fingerprint=velocity_hash(vectors),
            kinetic_energy=kinetic,
            total_energy=total,
            instantaneous_temperature_kelvin=2
            * kinetic
            / (frame.degrees_of_freedom * GAS_CONSTANT),
            energy_deviation=frame.energy_deviation + total - frame.total_energy,
        )
    changed.validate_integrity()
    with pytest.raises(WorkflowError, match=field):
        compatible_boundary(frame, changed)


def test_completed_readers_do_not_execute_engines(tmp_path, prepared, monkeypatch):
    root = tmp_path / "run"
    start(root, prepared, segments=3)
    from island.workflows import chain

    def forbidden(*a, **kw):
        raise AssertionError("inspection must not run a scientific engine")

    for name in (
        "initialize_velocities",
        "minimize_geometry",
        "run_dynamics_segment",
        "resume_dynamics",
    ):
        monkeypatch.setattr(chain, name, forbidden)
    monkeypatch.setattr(OpenMMSinglePointEvaluator, "open_session", forbidden)
    monkeypatch.setattr(OpenMMSinglePointEvaluator, "evaluate", forbidden)
    monkeypatch.setattr(
        chain.AmberToolsParameterizationEngine, "parameterize", forbidden
    )
    assert workflow_status(root)["status"] == "completed"
    assert len(read_workflow_frames(root)) == 7
    assert resume_workflow(root)["status"] == "completed"


def test_derived_kinetic_values_use_existing_arithmetic_tolerance(tmp_path, prepared):
    from island.dynamics.thermal import GAS_CONSTANT
    from island.workflows.chain import _load_segment
    from island.workflows.consistency import compatible_boundary

    root = tmp_path / "run"
    m = start(root, prepared)
    frame = _load_segment(root, m["segments"][0]).final_state
    kinetic = frame.kinetic_energy + 1e-12
    changed = replace(
        frame,
        kinetic_energy=kinetic,
        total_energy=frame.total_energy + (kinetic - frame.kinetic_energy),
        energy_deviation=frame.energy_deviation + (kinetic - frame.kinetic_energy),
        instantaneous_temperature_kelvin=2
        * kinetic
        / (frame.degrees_of_freedom * GAS_CONSTANT),
    )
    changed.validate_integrity()
    compatible_boundary(frame, changed)
