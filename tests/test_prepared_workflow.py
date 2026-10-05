"""Synthetic software interoperability, not force-field scientific acceptance."""

from dataclasses import replace

import pytest
from test_pcff_evaluation import bound, parameters, synthetic  # noqa: F401

# ruff: noqa: F811 -- fixture injection
from test_prepared_bundles import opls, opls_fixture, tree  # noqa: F401

from island.exceptions import DynamicsInputError, WorkflowError
from island.forcefields import save_prepared_forcefield
from island.minimization import MinimizationOptions
from island.workflows import (
    PreparedWorkflowConfig,
    prepared_workflow_status,
    read_prepared_workflow_frames,
    resume_prepared_workflow,
    start_prepared_bundle_workflow,
    storage,
)
from island.workflows import prepared as workflow


@pytest.fixture
def case(opls, tmp_path):
    system, prepared, sources = opls
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "input", sources=sources
    )
    config = PreparedWorkflowConfig(
        total_steps=4,
        segment_steps=2,
        recording_interval=1,
        max_evaluations_per_segment=4,
        max_frames_per_segment=3,
        minimization=MinimizationOptions(max_iterations=100, max_evaluations=200),
    )
    return root, tmp_path / "run", sources, config


def test_split_complete_and_inspection_without_context(case, monkeypatch):
    root, run, sources, config = case
    before = tree(root)
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert m["status"] == "paused", m
    assert m["accepted_step"] == 2
    m = resume_prepared_workflow(run, sources=sources)
    assert m["status"] == "completed"
    import openmm

    monkeypatch.setattr(
        openmm, "Context", lambda *a, **k: pytest.fail("inspection Context")
    )
    assert prepared_workflow_status(run, sources=sources) == m
    assert [f.step for f in read_prepared_workflow_frames(run, sources=sources)] == [
        0,
        1,
        2,
        3,
        4,
    ]
    assert resume_prepared_workflow(run, sources=sources) == m
    assert tree(root) == before


@pytest.mark.parametrize("completed", [False, True])
@pytest.mark.parametrize(
    "key,value",
    [
        ("temperature_kelvin", 777),
        ("velocity_seed", 2),
        ("timestep_fs", 0.2),
        ("friction_per_ps", 3),
        ("thermostat_seed", 5),
    ],
)
def test_rechecksummed_config_rejected(case, completed, key, value):
    root, run, sources, config = case
    m = start_prepared_bundle_workflow(
        root, run, config, sources=sources, segments=2 if completed else 1
    )
    m["config"][key] = value
    workflow._manifest(run, m)
    for operation in (
        prepared_workflow_status,
        read_prepared_workflow_frames,
        resume_prepared_workflow,
    ):
        with pytest.raises(WorkflowError, match="Configuration"):
            operation(run, sources=sources)


def test_candidate_failure_preserves_boundary(case, monkeypatch):
    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    previous = (run / "manifest.json").read_bytes()
    cp = (run / m["checkpoint"]).read_bytes()
    original = workflow._publish_manifest

    def fail(root, candidate, sources):
        if candidate["accepted_step"] > 2:
            raise WorkflowError("injected candidate failure")
        return original(root, candidate, sources)

    monkeypatch.setattr(workflow, "_publish_manifest", fail)
    with pytest.raises(WorkflowError, match="candidate"):
        resume_prepared_workflow(run, sources=sources)
    assert (run / "manifest.json").read_bytes() == previous
    assert (run / m["checkpoint"]).read_bytes() == cp
    assert prepared_workflow_status(run, sources=sources)["accepted_step"] == 2


def test_failed_minimum_inspectable(case):
    root, run, sources, config = case
    config = replace(config, minimization=MinimizationOptions(max_evaluations=1))
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert m["status"] == "stage_failed"
    assert prepared_workflow_status(run, sources=sources) == m
    assert read_prepared_workflow_frames(run, sources=sources) == ()
    with pytest.raises(WorkflowError, match="Failed stage"):
        resume_prepared_workflow(run, sources=sources)


def test_writer_exclusion_and_missing_input(case):
    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    with storage.writer(run), pytest.raises(WorkflowError):
        resume_prepared_workflow(run, sources=sources)
    (run / m["input_files"][1]).unlink()
    with pytest.raises(WorkflowError):
        prepared_workflow_status(run, sources=sources)


@pytest.mark.parametrize("offset", [5e-9, -5e-9, 1.0])
def test_fresh_boundary_tolerance_and_failed_startup(case, monkeypatch, offset):
    from island.evaluation.oplsaa_session import OPLSEvaluationSession

    root, run, sources, config = case
    start_prepared_bundle_workflow(root, run, config, sources=sources)
    original = OPLSEvaluationSession.evaluate_fresh
    calls = []

    def modified(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        calls.append(1)
        if len(calls) == 1:
            components = dict(result.energy_components)
            components[next(iter(components))] += offset
            return replace(
                result,
                potential_energy=result.potential_energy + offset,
                energy_components=components,
            )
        return result

    monkeypatch.setattr(OPLSEvaluationSession, "evaluate_fresh", modified)
    result = resume_prepared_workflow(run, sources=sources)
    if abs(offset) < 1e-8:
        assert result["status"] == "completed"
        assert len(read_prepared_workflow_frames(run, sources=sources)) == 5
        assert resume_prepared_workflow(run, sources=sources) == result
    else:
        assert result["status"] == "stage_failed"
        assert result["accepted_step"] == 2
        assert len(read_prepared_workflow_frames(run, sources=sources)) == 3


def rewrite_setup(run, m, transform):
    path = run / m["bundle"]
    p = storage.decode(storage.read_json(path))
    transform(p)
    raw = storage.json_bytes(storage.encode(p))
    path.write_bytes(raw)
    m["files"][m["bundle"]] = {"bytes": len(raw), "sha256": storage.checksum(raw)}
    workflow._manifest(run, m)


@pytest.mark.parametrize(
    "field",
    ["mass", "element", "stereo", "coordinate", "velocity", "minimum", "family"],
)
def test_semantically_changed_setup(case, field):
    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)

    def alter(p):
        site = p["starting_system"]["sites"][0]
        if field == "mass":
            site["mass"] += 1
        elif field == "element":
            site["element"] = "F"
        elif field == "stereo":
            site["metadata"]["cip_label"] = "R"
        elif field == "coordinate":
            p["starting_system"]["coordinates"][site["id"]] = (9.0, 0.0, 0.0)
        elif field == "velocity":
            from island.dynamics import initialize_velocities

            system, _, _ = workflow._load_setup(run, m, sources)
            p["initialization"] = workflow.bundle.record(
                initialize_velocities(system, temperature_kelvin=300, seed=12)
            )
        elif field == "minimum":
            p["minimum"]["initial_coordinates"][site["id"]] = (8.0, 0.0, 0.0)
        else:
            p["prepared_identity"] = "0" * 64

    rewrite_setup(run, m, alter)
    for op in (
        prepared_workflow_status,
        read_prepared_workflow_frames,
        resume_prepared_workflow,
    ):
        with pytest.raises(WorkflowError):
            op(run, sources=sources)


@pytest.mark.parametrize("budget", ["evaluations", "frames", "segments"])
def test_budget_exhaustion_is_not_completion(case, budget):
    root, run, sources, config = case
    if budget == "evaluations":
        config = replace(config, max_evaluations_per_segment=3)
    elif budget == "frames":
        config = replace(config, max_frames_per_segment=2)
    else:
        config = replace(config, max_segments=1)
    m = start_prepared_bundle_workflow(root, run, config, sources=sources, segments=2)
    assert m["status"] == "budget_exhausted"
    assert m["accepted_step"] < 4
    assert m["checkpoint"] is not None
    assert prepared_workflow_status(run, sources=sources) == m


def test_offline_inspection_without_scientific_dependencies(case, monkeypatch):
    import sys

    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources, segments=2)
    for name in ("openmm", "scipy", "rdkit", "parmed", "foyer"):
        monkeypatch.setitem(sys.modules, name, None)
    assert prepared_workflow_status(run, sources=sources) == m
    assert len(read_prepared_workflow_frames(run, sources=sources)) == 5
    assert resume_prepared_workflow(run, sources=sources) == m


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
def test_amber_dispatch_without_parameterization(tmp_path, monkeypatch, family):
    from test_prepared_bundles import amber_fixture

    from island.forcefields import AmberToolsParameterizationEngine

    system, prepared, artifacts = amber_fixture(tmp_path, family)
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "input", artifacts=artifacts
    )

    def forbidden(*a, **k):
        pytest.fail("parameterization during workflow")

    monkeypatch.setattr(AmberToolsParameterizationEngine, "parameterize", forbidden)
    config = PreparedWorkflowConfig(
        total_steps=2,
        segment_steps=1,
        recording_interval=1,
        max_evaluations_per_segment=3,
        max_frames_per_segment=2,
    )
    m = start_prepared_bundle_workflow(root, tmp_path / "run", config)
    assert m["status"] == "paused", m
    assert resume_prepared_workflow(tmp_path / "run")["status"] == "completed"


def test_independently_valid_wrong_origin_rejected_before_publication(case):
    from island.dynamics import create_dynamics_checkpoint, run_dynamics_segment
    from island.forcefields import create_evaluator

    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    prior = (run / "manifest.json").read_bytes()
    system, prepared, init = workflow._load_setup(run, m, sources)
    with create_evaluator(system, prepared).open_session() as session:
        segment = run_dynamics_segment(
            system,
            session,
            {s: tuple(2 * v for v in xyz) for s, xyz in init.velocities.items()},
            config.langevin(2),
        )
    cp = create_dynamics_checkpoint(segment)
    name = workflow._put(
        run,
        m,
        "rogue",
        storage.json_bytes(
            {"payload": segment.payload, "sha256": segment.content_checksum}
        ),
    )
    checkpoint = workflow._put(
        run,
        m,
        "rogue-cp",
        storage.json_bytes({"payload": cp.payload, "sha256": cp.content_checksum}),
    )
    old = m["segments"][0]
    m["segments"] = [name]
    m["checkpoint"] = checkpoint
    for stage in m["stages"]:
        if stage.get("record") == old:
            stage["record"] = name
    with pytest.raises(WorkflowError, match="origin"):
        workflow._publish_manifest(run, m, sources)
    assert (run / "manifest.json").read_bytes() == prior
    assert prepared_workflow_status(run, sources=sources)["accepted_step"] == 2


@pytest.mark.parametrize("failure", ["backend", "final"])
def test_failed_segment_retains_diagnostic_and_prior_boundary(
    case, monkeypatch, failure
):
    from island.evaluation.oplsaa_session import OPLSEvaluationSession
    from island.exceptions import EvaluationError

    root, run, sources, config = case
    before = start_prepared_bundle_workflow(root, run, config, sources=sources)
    calls = []
    original = OPLSEvaluationSession.evaluate_fresh

    def fresh(self, *args, **kwargs):
        calls.append(1)
        result = original(self, *args, **kwargs)
        if len(calls) == 2:
            components = dict(result.energy_components)
            components[next(iter(components))] += 1
            return replace(
                result,
                potential_energy=result.potential_energy + 1,
                energy_components=components,
            )
        return result

    def broken(self, *args, **kwargs):
        raise EvaluationError("injected backend error")

    if failure == "backend":
        monkeypatch.setattr(OPLSEvaluationSession, "evaluate", broken)
    else:
        monkeypatch.setattr(OPLSEvaluationSession, "evaluate_fresh", fresh)
    result = resume_prepared_workflow(run, sources=sources)
    assert result["status"] == "stage_failed"
    assert result["checkpoint"] == before["checkpoint"]
    assert result["accepted_step"] == 2
    diagnostic = workflow._load_segment(run, result["stages"][-1]["record"])
    from island.dynamics import create_dynamics_checkpoint

    with pytest.raises(DynamicsInputError):
        create_dynamics_checkpoint(diagnostic)
    assert prepared_workflow_status(run, sources=sources) == result


@pytest.mark.parametrize("field", ["evaluations", "reason"])
def test_rechecksummed_accounting_contradiction(case, field):
    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    m["stages"][-1][field] = 999 if field == "evaluations" else "invented"
    workflow._manifest(run, m)
    with pytest.raises(WorkflowError, match="accounting"):
        prepared_workflow_status(run, sources=sources)


def test_pcff_dispatch_and_no_preparation(bound, tmp_path, monkeypatch):
    import island.forcefields.pcff as backend
    from island.forcefields import PreparedForceFieldSources, adopt_forcefield

    system, spec, _ = bound
    path = tmp_path / "source.frc"
    path.write_bytes(spec.assignment.source.raw)
    sources = PreparedForceFieldSources(pcff_frc=path)
    prepared = adopt_forcefield(system, "pcff", spec)
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "input", sources=sources
    )

    def forbidden(*a, **kw):
        pytest.fail("scientific preparation during workflow")

    for name in (
        "type_pcff_atoms",
        "assign_automatic_pcff_charges",
        "assign_pcff_parameters",
        "define_pcff_model",
    ):
        monkeypatch.setattr(backend, name, forbidden)
    config = PreparedWorkflowConfig(
        total_steps=2,
        segment_steps=1,
        recording_interval=1,
        max_evaluations_per_segment=3,
        max_frames_per_segment=2,
        minimization=MinimizationOptions(max_iterations=500, max_evaluations=1000),
    )
    m = start_prepared_bundle_workflow(root, tmp_path / "run", config, sources=sources)
    assert m["status"] == "paused", m
    assert (
        resume_prepared_workflow(tmp_path / "run", sources=sources)["status"]
        == "completed"
    )


def test_no_final_check_budget_is_not_resumable(case):
    root, run, sources, config = case
    config = replace(config, max_evaluations_per_segment=1)
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert m["status"] == "budget_exhausted"
    assert m["checkpoint"] is None and m["accepted_step"] == 0
    assert read_prepared_workflow_frames(run, sources=sources) == ()
    diagnostic = workflow._load_segment(run, m["stages"][-1]["record"])
    assert diagnostic.evaluations == 1
    from island.dynamics import create_dynamics_checkpoint

    with pytest.raises(DynamicsInputError):
        create_dynamics_checkpoint(diagnostic)


@pytest.mark.parametrize("change", ["stage", "manifest"])
def test_rechecksummed_status_contradiction(case, change):
    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    if change == "stage":
        m["stages"][-1]["status"] = "completed"
    else:
        m["status"] = "ready"
    workflow._manifest(run, m)
    with pytest.raises(WorkflowError, match="status"):
        prepared_workflow_status(run, sources=sources)


def append_diagnostic(run, manifest, segment):
    from copy import deepcopy

    candidate = deepcopy(manifest)
    name = workflow._put(
        run,
        candidate,
        "failed-attempt",
        storage.json_bytes(
            {"payload": segment.payload, "sha256": segment.content_checksum}
        ),
    )
    candidate["stages"].append(
        {
            "stage": "dynamics",
            "status": "stage_failed",
            "record": name,
            "evaluations": segment.evaluations,
            "reason": segment.termination_reason,
        }
    )
    candidate["status"] = "stage_failed"
    return candidate


def failed_attempt(system, prepared, velocities, options, monkeypatch, checkpoint=None):
    from island.dynamics import resume_dynamics, run_dynamics_segment
    from island.exceptions import EvaluationError
    from island.forcefields import create_evaluator

    with create_evaluator(system, prepared).open_session() as session:

        def fail(*a, **k):
            raise EvaluationError("injected trial failure")

        monkeypatch.setattr(session, "evaluate", fail)
        segment = (
            run_dynamics_segment(system, session, velocities, options)
            if checkpoint is None
            else resume_dynamics(checkpoint, system, session, options)
        )
    segment.validate_integrity()
    assert segment.termination_reason == "evaluation_failed"
    return segment


def test_unrelated_failed_attempt_rejected_before_publication(case, monkeypatch):
    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    system, prepared, init = workflow._load_setup(run, m, sources)
    segment = failed_attempt(
        system,
        prepared,
        {s: tuple(2 * v for v in xyz) for s, xyz in init.velocities.items()},
        config.langevin(2),
        monkeypatch,
    )
    assert segment.frames[0].step == 0 and m["accepted_step"] == 2
    candidate = append_diagnostic(run, m, segment)
    prior = (run / "manifest.json").read_bytes()
    checkpoint = (run / m["checkpoint"]).read_bytes()
    with pytest.raises(WorkflowError):
        workflow._publish_manifest(run, candidate, sources)
    assert (run / "manifest.json").read_bytes() == prior
    assert (run / m["checkpoint"]).read_bytes() == checkpoint
    # Even a correctly rechecksummed externally installed contradiction is rejected.
    workflow._manifest(run, candidate)
    for operation in (
        prepared_workflow_status,
        read_prepared_workflow_frames,
        resume_prepared_workflow,
    ):
        with pytest.raises(WorkflowError, match="attempt|boundary|origin|lineage"):
            operation(run, sources=sources)


@pytest.mark.parametrize("change", ["origin", "parent", "physical", "rng"])
def test_semantically_valid_failed_resume_contradictions(case, monkeypatch, change):
    from island.dynamics import (
        DynamicsSegment,
        DynamicsSegmentOptions,
        create_dynamics_checkpoint,
        load_dynamics_checkpoint,
        run_dynamics_segment,
    )
    from island.dynamics import _checkpoint_data as data
    from island.forcefields import create_evaluator

    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    system, prepared, init = workflow._load_setup(run, m, sources)
    cp = load_dynamics_checkpoint(run / m["checkpoint"])
    if change == "origin":
        with create_evaluator(system, prepared).open_session() as session:
            other = run_dynamics_segment(
                system,
                session,
                {s: tuple(2 * v for v in xyz) for s, xyz in init.velocities.items()},
                config.langevin(2),
            )
        cp = create_dynamics_checkpoint(other)
    options = DynamicsSegmentOptions(2, 4, 3, 1)
    segment = failed_attempt(
        system, prepared, init.velocities, options, monkeypatch, cp
    )
    p = segment.payload
    if change == "parent":
        p["lineage"][-1]["parent_checksum"] = "0" * 64
    elif change == "physical":
        p["physical"]["friction_per_ps"] = 9.0
        p["origin"]["trajectory_fingerprint"] = data.trajectory_identity(p)
    elif change == "rng":
        # A genuine startup failure has consumed no normals; changing a valid
        # PCG64 state keeps native local integrity but breaks the parent relation.
        from island.dynamics import resume_dynamics
        from island.exceptions import EvaluationError

        with create_evaluator(system, prepared).open_session() as session:

            def fail(*a, **k):
                raise EvaluationError("startup observation failed")

            monkeypatch.setattr(session, "evaluate_fresh", fail)
            p = resume_dynamics(cp, system, session, options).payload
        p["rng"]["state"]["state"]["state"] += 1
    segment = DynamicsSegment(data.canonical(p), data.checksum(p))
    segment.validate_integrity()
    candidate = append_diagnostic(run, m, segment)
    prior = (run / "manifest.json").read_bytes()
    with pytest.raises(WorkflowError):
        workflow._publish_manifest(run, candidate, sources)
    assert (run / "manifest.json").read_bytes() == prior
    workflow._manifest(run, candidate)
    for operation in (
        prepared_workflow_status,
        read_prepared_workflow_frames,
        resume_prepared_workflow,
    ):
        with pytest.raises(WorkflowError, match="attempt|boundary|origin|lineage|RNG"):
            operation(run, sources=sources)


@pytest.mark.parametrize("failure", ["backend", "startup", "final"])
def test_legitimate_first_attempt_failure_remains_inspectable(
    case, monkeypatch, failure
):
    from island.evaluation.oplsaa_session import OPLSEvaluationSession
    from island.exceptions import EvaluationError

    root, run, sources, config = case
    execute = workflow.run_dynamics_segment

    def injected(*args, **kwargs):
        session = args[1]
        fresh = session.evaluate_fresh
        calls = []

        def fail(*a, **k):
            raise EvaluationError("legitimate first attempt backend failure")

        def checked(*a, **k):
            calls.append(1)
            result = fresh(*a, **k)
            if len(calls) == 2:
                components = dict(result.energy_components)
                components[next(iter(components))] += 1.0
                return replace(
                    result,
                    potential_energy=result.potential_energy + 1.0,
                    energy_components=components,
                )
            return result

        assert isinstance(session, OPLSEvaluationSession)
        monkeypatch.setattr(
            session,
            "evaluate" if failure == "backend" else "evaluate_fresh",
            checked if failure == "final" else fail,
        )
        return execute(*args, **kwargs)

    monkeypatch.setattr(workflow, "run_dynamics_segment", injected)
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert (
        m["status"] == "stage_failed"
        and m["accepted_step"] == 0
        and m["checkpoint"] is None
    )
    assert prepared_workflow_status(run, sources=sources) == m
    assert read_prepared_workflow_frames(run, sources=sources) == ()


def test_budget_exhausted_history_can_continue(case):
    root, run, sources, config = case
    config = replace(config, max_evaluations_per_segment=3)
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert m["status"] == "budget_exhausted" and m["accepted_step"] == 1
    m = resume_prepared_workflow(run, sources=sources, segments=5)
    # Each call stops at the next budget boundary; explicit continuation is safe.
    while m["status"] == "budget_exhausted":
        m = resume_prepared_workflow(run, sources=sources)
    assert m["status"] == "completed" and m["accepted_step"] == 4
    assert [
        f.step for f in read_prepared_workflow_frames(run, sources=sources)
    ] == list(range(5))


def test_first_failed_attempt_must_use_saved_velocities(case, monkeypatch):
    root, run, sources, config = case
    monkeypatch.setattr(workflow, "_advance", lambda root, m, count, sources: m)
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert m["status"] == "ready"
    system, prepared, init = workflow._load_setup(run, m, sources)
    segment = failed_attempt(
        system,
        prepared,
        {s: tuple(2 * v for v in xyz) for s, xyz in init.velocities.items()},
        config.langevin(2),
        monkeypatch,
    )
    candidate = append_diagnostic(run, m, segment)
    with pytest.raises(WorkflowError, match="origin"):
        workflow._publish_manifest(run, candidate, sources)
    assert prepared_workflow_status(run, sources=sources) == m


def test_rejected_attempts_do_not_advance_parent_or_counters(case, monkeypatch):
    from island.dynamics import DynamicsSegmentOptions, load_dynamics_checkpoint

    root, run, sources, config = case
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    system, prepared, init = workflow._load_setup(run, m, sources)
    cp = load_dynamics_checkpoint(run / m["checkpoint"])
    for _ in range(2):
        segment = failed_attempt(
            system,
            prepared,
            init.velocities,
            DynamicsSegmentOptions(2, 4, 3, 1),
            monkeypatch,
            cp,
        )
        m = append_diagnostic(run, m, segment)
        workflow._publish_manifest(run, m, sources)
    assert prepared_workflow_status(run, sources=sources) == m
    assert m["accepted_step"] == 2
    assert load_dynamics_checkpoint(run / m["checkpoint"]) == cp
    assert len(read_prepared_workflow_frames(run, sources=sources)) == 3
