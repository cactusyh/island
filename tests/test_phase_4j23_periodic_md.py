"""J23 real retained-source software contracts; no production acceptance implied."""

import builtins
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from island.exceptions import ValidationError
from island.periodic import _signed
from island.periodic_md import (
    PeriodicMDWorkflowConfig,
    compare_periodic_md_workflows,
    load_periodic_md_workflow,
    periodic_md_workflow_status,
    read_periodic_md_frames,
    resume_periodic_md_workflow,
    save_periodic_md_workflow,
    start_periodic_md_workflow,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from periodic_j22_fixtures import retained_fixture


@pytest.fixture(scope="module", params=("pcff", "oplsaa", "gaff", "gaff2"))
def backend(request):
    return retained_fixture(request.param, method="ewald", tolerance=1e-5)[0]


@pytest.fixture(scope="module")
def gaff():
    return retained_fixture("gaff", method="ewald", tolerance=1e-5)[0]


def child(code, *arguments):
    env = {**os.environ, "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}
    result = subprocess.run(
        [sys.executable, "-c", code, *map(str, arguments)],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


@pytest.mark.parametrize("mode", ("nvt", "npt", "annealing"))
def test_four_family_exact_relocated_child_resume(backend, mode, tmp_path):
    cfg = PeriodicMDWorkflowConfig(
        backend.config.payload["family"],
        mode=mode,
        steps=6,
        minimize=False,
        barostat_interval=1,
        pressure_bar=1.0 if mode == "npt" else None,
        schedule=((2, 300.0), (4, 350.0), (6, 310.0)) if mode == "annealing" else (),
    )
    full = tmp_path / "full"
    split = tmp_path / "split"
    moved = tmp_path / "moved"
    start_periodic_md_workflow(backend, cfg, full)
    start_periodic_md_workflow(backend, cfg, split, steps=3)
    boundary = read_periodic_md_frames(split)[-1]
    assert boundary == read_periodic_md_frames(full)[3]
    save_periodic_md_workflow(split, moved)
    manifest_id = load_periodic_md_workflow(moved).identity
    child(
        """
import builtins,sys
old=builtins.__import__
def guard(name, globals=None, locals=None, fromlist=(), level=0):
 if level==0 and name.split('.')[0] in {'rdkit','foyer','parmed','scipy'}:raise RuntimeError('forbidden scientific reconstruction '+name)
 return old(name,globals,locals,fromlist,level)
builtins.__import__=guard
from island.periodic_md import resume_periodic_md_workflow
resume_periodic_md_workflow(sys.argv[1],expected_identity=sys.argv[2])
""",
        moved,
        manifest_id,
    )
    actual = read_periodic_md_frames(moved)
    expected = read_periodic_md_frames(full)
    assert (
        actual == expected
    )  # Includes ALL binary checkpoint bytes, random state, box, velocities, time.
    assert compare_periodic_md_workflows(full, moved)["passed"]
    assert periodic_md_workflow_status(moved)["status"] == "complete"
    assert all(f["graph_identity"] == backend.final_graph.identity for f in actual)
    if mode == "npt":
        assert actual[-1]["box_vectors_nm"] != actual[0]["box_vectors_nm"]


def test_four_family_minimization(backend, tmp_path):
    cfg = PeriodicMDWorkflowConfig(
        backend.config.payload["family"], steps=1, minimize=True
    )
    start_periodic_md_workflow(backend, cfg, tmp_path / "min", steps=0)
    diagnostic = json.loads((tmp_path / "min/minimization.json").read_text())["payload"]
    assert diagnostic["independent_context_verified"]
    assert diagnostic["maximum_force_kj_mol_nm"] <= cfg.minimization_tolerance_kj_mol_nm
    assert (
        diagnostic["final_energy_kj_mol"] <= diagnostic["initial_energy_kj_mol"] + 1e-8
    )


@pytest.fixture
def run(gaff, tmp_path):
    path = tmp_path / "run"
    cfg = PeriodicMDWorkflowConfig(
        "GAFF",
        mode="npt",
        steps=4,
        minimize=False,
        pressure_bar=1.0,
        barostat_interval=1,
    )
    start_periodic_md_workflow(gaff, cfg, path, steps=1)
    return path, cfg


@pytest.mark.parametrize(
    "field,value",
    [
        ("temperature_kelvin", 310.0),
        ("pressure_bar", 2.0),
        ("timestep_fs", 0.3),
        ("seed", 99),
        ("barostat_seed", 99),
        ("barostat_interval", 2),
        ("family", "PCFF"),
    ],
)
def test_changed_configuration_rejects(run, field, value):
    path, cfg = run
    before = read_periodic_md_frames(path)
    with pytest.raises(ValidationError, match="Changed workflow configuration"):
        resume_periodic_md_workflow(
            path, expected_config=replace(cfg, **{field: value})
        )
    assert read_periodic_md_frames(path) == before


def test_changed_schedule_rejects(gaff, tmp_path):
    cfg = PeriodicMDWorkflowConfig(
        "GAFF",
        mode="annealing",
        steps=4,
        minimize=False,
        schedule=((2, 300.0), (4, 350.0)),
    )
    path = tmp_path / "anneal"
    start_periodic_md_workflow(gaff, cfg, path, steps=1)
    with pytest.raises(ValidationError, match="configuration"):
        resume_periodic_md_workflow(
            path, expected_config=replace(cfg, schedule=((2, 300.0), (4, 360.0)))
        )


def test_offline_inspection_no_scientific_import(run):
    path, _ = run
    child(
        """
import builtins,sys
old=builtins.__import__
def guard(name, globals=None, locals=None, fromlist=(), level=0):
 if level==0 and name.split('.')[0] in {'openmm','rdkit','foyer','parmed','scipy'}:raise RuntimeError('forbidden offline import '+name)
 return old(name,globals,locals,fromlist,level)
builtins.__import__=guard
from island.periodic_md import periodic_md_workflow_status,read_periodic_md_frames
assert periodic_md_workflow_status(sys.argv[1])['step']==1
assert len(read_periodic_md_frames(sys.argv[1]))==2
""",
        path,
    )


def test_missing_openmm(run, monkeypatch):
    path, _ = run
    original = builtins.__import__

    def block(name, *args, **kw):
        if name == "openmm":
            raise ImportError("blocked")
        return original(name, *args, **kw)

    monkeypatch.setattr(builtins, "__import__", block)
    with pytest.raises(ValidationError, match="optional dependency OpenMM"):
        resume_periodic_md_workflow(path)
    assert periodic_md_workflow_status(path)["step"] == 1


@pytest.mark.parametrize(
    "key",
    (
        "checkpoint_sha256",
        "step",
        "positions_nm",
        "box_vectors_nm",
        "temperature_kelvin",
        "backend_identity",
    ),
)
def test_tampered_checkpoint(run, key):
    path, _ = run
    p = json.loads((path / "generation-00000000/record.json").read_text())["payload"]
    f = p["frames"][-1]
    if key == "positions_nm":
        f[key][0][0] = float("nan")
    elif key == "box_vectors_nm":
        f[key][0][0] = -1
    elif key in ("step", "temperature_kelvin"):
        f[key] += 1
    else:
        f[key] = "bad"
    (path / "generation-00000000/record.json").write_text(
        json.dumps({"payload": p, "sha256": "bad"})
    )
    with pytest.raises(ValidationError):
        load_periodic_md_workflow(path)


def test_semantically_resigned_state_must_match_binary(run):
    path, _ = run
    file = path / "generation-00000000/record.json"
    p = json.loads(file.read_text())["payload"]
    f = p["frames"][-1]
    f["positions_nm"][0][0] += 0.01
    p["frames"][-1] = json.loads(
        _signed({k: v for k, v in f.items() if k != "identity"})
    )["payload"]
    file.write_text(_signed({k: v for k, v in p.items() if k != "identity"}))
    with pytest.raises(ValidationError, match="restoration"):
        resume_periodic_md_workflow(path)


def test_failed_publication_preserves_checkpoint(run, monkeypatch):
    import island.periodic_md as workflow

    path, _ = run
    before = load_periodic_md_workflow(path).json_text

    def fail(*args):
        raise OSError("injected publication failure")

    monkeypatch.setattr(workflow, "_publish_directory", fail)
    with pytest.raises(OSError, match="injected"):
        resume_periodic_md_workflow(path)
    assert load_periodic_md_workflow(path).json_text == before


@pytest.mark.parametrize(
    "kwargs",
    [
        {"platform": "CPU"},
        {"mode": "npt"},
        {"seed": True},
        {"temperature_kelvin": float("inf")},
        {"mode": "annealing", "schedule": ((1, 300.0),)},
    ],
)
def test_unsupported_configuration(kwargs):
    with pytest.raises(ValidationError):
        PeriodicMDWorkflowConfig("PCFF", **kwargs)


def test_wrong_family(gaff, tmp_path):
    with pytest.raises(ValidationError, match="Wrong-family"):
        start_periodic_md_workflow(
            gaff, PeriodicMDWorkflowConfig("PCFF"), tmp_path / "wrong"
        )


@pytest.mark.parametrize("key", ("graph", "source", "backend", "boundary"))
def test_stale_manifest_bindings(run, key):
    path, _ = run
    file = path / "manifest.json"
    p = json.loads(file.read_text())["payload"]
    if key == "graph":
        p["backend"]["final_graph"]["sites"][0]["formal_charge"] = 1
    elif key == "source":
        p["backend"]["typed_graph"]["source"]["sha256"] = "b" * 64
    elif key == "backend":
        p["backend"]["identity"] = "b" * 64
    else:
        p["backend"]["system"]["box"]["periodic"] = [False, True, True]
    file.write_text(_signed({k: v for k, v in p.items() if k != "identity"}))
    with pytest.raises(ValidationError):
        load_periodic_md_workflow(path)


def test_stale_trusted_checkpoint(run):
    path, _ = run
    with pytest.raises(ValidationError, match="identity"):
        load_periodic_md_workflow(path, expected_identity="b" * 64)
    with pytest.raises(ValidationError, match="checkpoint"):
        load_periodic_md_workflow(path, expected_checkpoint_identity="b" * 64)
    with pytest.raises(ValidationError, match="backend"):
        load_periodic_md_workflow(path, expected_backend_identity="b" * 64)


def test_pinned_lammps_remains_blocked():
    p = json.loads(
        (
            Path(__file__).resolve().parents[1] / "docs/evidence/phase_4j22.json"
        ).read_text()
    )
    assert (
        p["lammps_periodic_runtime_comparison"]["status"]
        == "blocked: KSPACE package unavailable"
    )
    assert not p["lammps_periodic_runtime_comparison"]["kspace_package"]
