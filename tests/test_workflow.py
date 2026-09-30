"""Synthetic software integration uses actual ParmEd/OpenMM, never live AmberTools."""

import shutil
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace

import pytest

from island.dynamics import load_dynamics_checkpoint
from island.exceptions import WorkflowBusyError, WorkflowError
from island.workflows import (
    WorkflowConfig,
    inspect_chain,
    read_workflow_frames,
    resume_workflow,
    start_prepared_workflow,
    start_workflow,
    storage,
    workflow_status,
)


@pytest.fixture
def prepared(tmp_path):
    pytest.importorskip("openmm")
    pytest.importorskip("scipy")
    from test_ambertools_integrity import _valid_result

    source = tmp_path / "source"
    source.mkdir()
    system, preparation = _valid_result(source)
    for name in (*preparation.record["artifact_sha256"], "input.mol2", "lineage.json"):
        (source / name).write_bytes(
            (source / "source.prmtop").read_bytes()
            if name == "result.prmtop"
            else b"independent synthetic fixture payload"
        )
    return system, preparation, source


def config(path, **kwargs):
    return WorkflowConfig(
        "[*]CC[*]",
        3,
        str(path),
        "gaff2",
        "provided",
        provided_charges={s: 0.0 for s in (101, 109, 117, 125)},
        total_steps=6,
        segment_steps=2,
        max_evaluations_per_segment=4,
        max_frames_per_segment=3,
        recording_interval=1,
        **kwargs,
    )


def start(path, prepared, **kw):
    return start_prepared_workflow(
        config(path), *prepared, evidence="synthetic_software_test", **kw
    )


def test_real_synthetic_relocation_cross_process_and_uninterrupted(tmp_path, prepared):
    system_before = deepcopy(prepared[0].to_dict())
    signature = prepared[1].record_signature
    path = tmp_path / "split"
    first = start(path, prepared)
    assert first["status"] == "paused", first
    assert first["accepted_step"] == 2
    assert (
        prepared[0].to_dict() == system_before
        and prepared[1].record_signature == signature
    )
    # New process cannot access the old bundle or the source artifact directory.
    moved = tmp_path / "relocated"
    shutil.move(path, moved)
    shutil.rmtree(prepared[2])
    code = """
import sys
from island.workflows import resume_workflow
from island.workflows import chain
# Resuming must not call any setup stage.
def forbidden(*a, **kw):
    raise AssertionError("setup repeated")
chain.inspect_chain = forbidden
chain.AmberToolsParameterizationEngine.parameterize = forbidden
chain.minimize_geometry = forbidden
chain.initialize_velocities = forbidden
print(resume_workflow(sys.argv[1])["status"])
"""
    for status in ("paused", "completed"):
        out = subprocess.run(
            [sys.executable, "-c", code, str(moved)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert out.returncode == 0, out.stderr
        assert out.stdout.strip() == status
    final = workflow_status(moved)
    assert final["accepted_step"] == 6
    frames = read_workflow_frames(moved)
    assert [f.step for f in frames] == list(range(7))
    # Independent uninterrupted integration from durable inputs, no setup rerun.
    from island.dynamics import run_dynamics_segment
    from island.evaluation import OpenMMSinglePointEvaluator
    from island.workflows.chain import _load_bundle

    starting, prep, init = _load_bundle(moved, final)
    with OpenMMSinglePointEvaluator(
        starting, prep.imported_result
    ).open_session() as session:
        whole = run_dynamics_segment(
            starting,
            session,
            init.velocities,
            replace(config(moved).langevin(6), max_evaluations=8, max_frames=7),
        )
    cp = load_dynamics_checkpoint(moved / final["checkpoint"])
    assert cp.state == whole.final_state
    assert cp.rng_state == whole.payload["rng"]["state"]
    assert cp.payload["counters"]["evaluations"] == 12
    assert whole.evaluations == 8
    assert cp.payload["origin"] == whole.payload["origin"]
    assert cp.state.time_ps == 6 * 0.1 * 0.001


def test_publication_failure_preserves_boundary(tmp_path, prepared, monkeypatch):
    path = tmp_path / "run"
    initial = start(path, prepared)
    assert initial["status"] == "paused", initial
    saved = (path / "manifest.json").read_bytes()
    original = storage.os.replace

    def fail(source, target):
        if str(target).endswith("manifest.json"):
            raise OSError("injected atomic publication interruption")
        return original(source, target)

    with monkeypatch.context() as patch:
        patch.setattr(storage.os, "replace", fail)
        with pytest.raises(WorkflowError, match="Publication failed"):
            resume_workflow(path)
    assert (path / "manifest.json").read_bytes() == saved
    assert workflow_status(path)["accepted_step"] == 2
    assert resume_workflow(path)["accepted_step"] == 4


def test_lock_corruption_and_existing_directory(tmp_path, prepared):
    path = tmp_path / "run"
    m = start(path, prepared)
    assert m["status"] == "paused", m
    with storage.writer(path), pytest.raises(WorkflowBusyError):
        resume_workflow(path)
    with pytest.raises(WorkflowError, match="new directory"):
        start(path, prepared)
    file = path / m["bundle"]
    file.write_bytes(file.read_bytes() + b" ")
    with pytest.raises(WorkflowError, match="size mismatch"):
        resume_workflow(path)


@pytest.mark.parametrize(
    "changes",
    [
        {"dp": 0},
        {"temperature_kelvin": float("nan")},
        {"force_field": "opls"},
        {"charge_method": "am1bcc"},
        {"template_seed": -1},
        {"max_segments": True},
        {"minimization": {}},
        {"schema": "v9"},
        {"provided_charges": {True: 0.0}},
    ],
)
def test_config_rejects(changes, tmp_path):
    with pytest.raises(WorkflowError):
        replace(config(tmp_path), **changes)


def test_data_codec_and_corruption(tmp_path):
    original = {17: {"nested": (3, 4)}, "17": [2**127, 1.2345678901234567]}
    assert (
        storage.decode(
            storage.read_json(
                _write(tmp_path, storage.json_bytes(storage.encode(original)))
            )
        )
        == original
    )
    for raw in (b'{"a":1,"a":2}', b"NaN", b"{"):
        with pytest.raises(WorkflowError):
            storage.read_json(_write(tmp_path, raw))
    with pytest.raises(WorkflowError):
        storage.decode({"map": [[17, 1], [17, 2]]})
    with pytest.raises(WorkflowError):
        storage.child(tmp_path, "../outside")


def _write(tmp_path, raw):
    p = tmp_path / "data.json"
    p.write_bytes(raw)
    return p


def test_missing_dependency_recorded(tmp_path, monkeypatch):
    monkeypatch.setattr("island.workflows.chain.find_spec", lambda name: None)
    m = start_workflow(config(tmp_path / "run"))
    assert m["status"] == "stage_failed"
    assert m["stages"][-1]["status"] == "unavailable"
    assert m["checkpoint"] is None


def test_inspect_inventory_and_limit(tmp_path):
    pytest.importorskip("rdkit")
    cfg = replace(config(tmp_path), provided_charges={}, dp=3)
    built = inspect_chain(cfg)
    assert built.number_of_sites == 20
    with pytest.raises(WorkflowError, match="100 sites"):
        inspect_chain(replace(cfg, dp=20))


def test_optional_import_isolation():
    code = """
import sys
class Block:
    def find_spec(self, fullname, *args):
        if fullname.split(".")[0] in {"openmm", "parmed", "rdkit", "scipy"}:
            raise AssertionError(fullname)
sys.meta_path.insert(0, Block())
import island.workflows
"""
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr


def test_minimization_nonconvergence_stops_before_initialization(
    tmp_path, prepared, monkeypatch
):
    from island.minimization import MinimizationOptions

    cfg = replace(
        config(tmp_path / "run"), minimization=MinimizationOptions(max_evaluations=1)
    )

    def forbidden(*a, **kw):
        raise AssertionError("initialization must not run")

    monkeypatch.setattr("island.workflows.chain.initialize_velocities", forbidden)
    m = start_prepared_workflow(cfg, *prepared, evidence="synthetic_software_test")
    assert m["status"] == "stage_failed"
    assert m["stages"][-1]["stage"] == "minimization"
    assert m["bundle"] is None and m["checkpoint"] is None
    assert workflow_status(cfg.output_directory)["status"] == "stage_failed"


def test_failed_dynamics_retains_published_boundary_and_diagnostics(
    tmp_path, prepared, monkeypatch
):
    path = tmp_path / "run"
    initial = start(path, prepared)
    from island.evaluation import OpenMMEvaluationSession

    def fail(*a, **kw):
        raise RuntimeError("injected propagation evaluation failure")

    with monkeypatch.context() as patch:
        patch.setattr(OpenMMEvaluationSession, "evaluate", fail)
        m = resume_workflow(path)
    assert m["status"] == "stage_failed"
    assert m["accepted_step"] == initial["accepted_step"]
    assert m["checkpoint"] == initial["checkpoint"]
    assert m["segments"] == initial["segments"]
    assert m["stages"][-1]["reason"] == "evaluation_failed"
    assert m["stages"][-1]["evaluations"] == 3  # startup, failed trial, fresh final
    assert workflow_status(path)["status"] == "stage_failed"
    assert resume_workflow(path)["accepted_step"] == 4
    assert [f.step for f in read_workflow_frames(path)] == list(range(5))


@pytest.mark.parametrize(
    "budget", ["max_evaluations_per_segment", "max_frames_per_segment"]
)
def test_budget_stop_is_saved_not_completed(tmp_path, prepared, budget):
    cfg = replace(config(tmp_path / "run"), **{budget: 2})
    m = start_prepared_workflow(cfg, *prepared, evidence="synthetic_software_test")
    assert m["status"] == "budget_exhausted"
    assert m["checkpoint"] is not None
    assert m["accepted_step"] < 2
    assert workflow_status(cfg.output_directory)["status"] == "budget_exhausted"


def test_provided_coverage_and_source_checksum_errors(tmp_path, prepared):
    cfg = replace(config(tmp_path / "coverage"), provided_charges={101: 0.0})
    m = start_prepared_workflow(cfg, *prepared, evidence="synthetic_software_test")
    assert m["status"] == "stage_failed" and m["bundle"] is None
    (prepared[2] / "result.prmtop").write_text("corrupt")
    m = start(tmp_path / "checksum", prepared)
    assert m["status"] == "stage_failed"
    assert "checksum" in m["stages"][-1]["message"]


def test_parameter_failure_records_stage_without_fallback(tmp_path, monkeypatch):
    from island.exceptions import AmberToolsStageError
    from island.workflows import chain

    monkeypatch.setattr(chain, "preflight", lambda *a, **kw: {})
    from test_dynamics import molecule

    monkeypatch.setattr(chain, "inspect_chain", lambda *a: molecule())
    calls = []

    def fail(self, system, options):
        calls.append(options.charge_method)
        raise AmberToolsStageError("sqm", "not converged", artifact_dir="retained")

    monkeypatch.setattr(chain.AmberToolsParameterizationEngine, "parameterize", fail)
    cfg = replace(
        config(tmp_path / "run"), charge_method="am1bcc", provided_charges=None
    )
    m = start_workflow(cfg)
    assert calls == ["am1bcc"]
    assert m["status"] == "stage_failed" and m["checkpoint"] is None
    assert m["stages"][-1]["stage"] == "parameterization"
    assert "not converged" in m["stages"][-1]["message"]


def test_incompatible_config_rejected_without_propagation(tmp_path, prepared):
    path = tmp_path / "run"
    m = start(path, prepared)
    m["config"]["temperature_kelvin"] = 301
    from island.workflows.chain import _manifest

    _manifest(path, m)
    with pytest.raises(WorkflowError, match="immutable prepared bundle"):
        resume_workflow(path)


def test_byte_budget_and_missing_artifact(tmp_path, prepared):
    cfg = replace(config(tmp_path / "small"), max_artifact_bytes=1)
    m = start_prepared_workflow(cfg, *prepared, evidence="synthetic_software_test")
    assert m["status"] == "stage_failed"
    path = tmp_path / "run"
    m = start(path, prepared)
    (path / m["checkpoint"]).unlink()
    with pytest.raises(WorkflowError):
        resume_workflow(path)


def test_signature_and_duplicate_boundary_rejection(tmp_path, prepared):
    from island.workflows.chain import _manifest

    path = tmp_path / "run"
    m = start(path, prepared)
    # Valid segment repeated in the manifest is not another continuation.
    m["segments"].append(m["segments"][0])
    _manifest(path, m)
    with pytest.raises(WorkflowError, match="Duplicate segment"):
        read_workflow_frames(path)
    m["segments"].pop()
    _manifest(path, m)
    # Replaying a separately named valid segment creates a conflicting boundary.
    duplicate = "replayed-segment.json"
    (path / duplicate).write_bytes((path / m["segments"][0]).read_bytes())
    m["files"][duplicate] = dict(m["files"][m["segments"][0]])
    m["segments"].append(duplicate)
    _manifest(path, m)
    with pytest.raises(WorkflowError, match="Conflicting adjacent segment boundary"):
        read_workflow_frames(path)
    m["segments"].pop()
    _manifest(path, m)
    # Reconstructing changed preparation data is checked against historical signatures.
    bundle_path = path / m["bundle"]
    payload = storage.decode(storage.read_json(bundle_path))
    payload["preparation"]["record"]["input_coordinates_angstrom"]["101"][0] += 1
    raw = storage.json_bytes(storage.encode(payload))
    bundle_path.write_bytes(raw)
    m["files"][m["bundle"]] = {"sha256": storage.checksum(raw), "bytes": len(raw)}
    _manifest(path, m)
    with pytest.raises(WorkflowError, match="signature|provenance"):
        resume_workflow(path)


def test_unsupported_build_and_attempt_budget(tmp_path, prepared):
    from island.exceptions import IslandError

    with pytest.raises(IslandError):
        inspect_chain(replace(config(tmp_path), psmiles="[*]=CC[*]"))
    cfg = replace(config(tmp_path / "run"), max_segments=1)
    m = start_prepared_workflow(
        cfg, *prepared, evidence="synthetic_software_test", segments=2
    )
    assert m["status"] == "budget_exhausted" and m["accepted_step"] == 2
    assert resume_workflow(cfg.output_directory)["accepted_step"] == 2
