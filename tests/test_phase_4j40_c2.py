"""Real long-run failures, provenance mutations and fresh-process integrity."""

import copy
import json
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

from island.analysis.equilibration import require
from island.analysis.equilibration_protocol import CorrectiveProtocol
from island.analysis.npt import _signed
from island.exceptions import AnalysisError, IslandError
from scripts.phase_4j40_c2_analysis import observations, trajectory
from scripts.validate_phase_4j40 import CASES
from scripts.validate_phase_4j40_c2_receipt import validate_run


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    configured = os.environ.get("ISLAND_J40_C2_RUN")
    if configured:
        return Path(configured)
    root = tmp_path_factory.mktemp("j40-c2")
    with tarfile.open("docs/evidence/phase_4j40_c2/artifacts.tar.xz") as archive:
        archive.extractall(root, filter="data")
    return root / "run"


@pytest.mark.parametrize("case", CASES)
def test_actual_long_protocol_attempted_for_three_replicas(run, case):
    for i in range(3):
        root = run / case / f"replica-{i}"
        backend, cfg, frames, complete = trajectory(root)
        assert cfg["nvt_steps"] == 40000
        assert cfg["npt_steps"] == 2000000
        assert len(frames) > 400
        assert frames[-1]["stage"] == "npt"
        assert frames[-1]["step"] > 40000
        assert cfg["checkpoint_interval"] == 100
        assert backend["final_graph_identity"]
        if complete:
            assert frames[-1]["step"] == cfg["nvt_steps"] + cfg["npt_steps"]
        else:
            failure = json.loads((root / "uninterrupted/failure.json").read_text())
            assert "less than twice the nonbonded cutoff" in failure["diagnostic"]
        record = observations(root)
        assert not record["scientific_acceptance"]
        assert record["pressure_status"] == "pressure_unavailable"


@pytest.mark.parametrize("case", CASES)
def test_fresh_process_partial_continuation_equals_actual_history(run, case):
    for i in range(3):
        root = run / case / f"replica-{i}"
        _, _, frames, complete = trajectory(root)
        attempt = json.loads((root / "resumed-attempt/attempt.json").read_text())
        assert attempt["frames"] == frames
        assert attempt["status"] == ("completed" if complete else "failed")
        assert json.loads((root / "resume-attempt.command.json").read_text())[
            "exit_code"
        ] == int(not complete)
        assert json.loads((root / "partial-comparison.json").read_text())[
            "byte_exact_retained_history"
        ]


@pytest.mark.parametrize(
    "mutation",
    (
        "coordinate",
        "velocity",
        "box",
        "checkpoint",
        "temperature",
        "pressure",
        "stage",
        "config",
        "source",
        "readiness",
    ),
)
def test_tampered_failed_history_rejected(run, tmp_path, mutation):
    source = run / "opls-network/replica-0"
    shutil.copyfile(source / "result.json", tmp_path / "result.json")
    shutil.copyfile(source / "backend.json", tmp_path / "backend.json")
    shutil.copyfile(source / "config.json", tmp_path / "config.json")
    target = tmp_path / "uninterrupted"
    target.mkdir()
    record = json.loads((source / "uninterrupted/failure.json").read_text())
    # Content changes must be detected even after the outer failure hash is renewed.
    f = record["frames"][-1]
    if mutation == "coordinate":
        f["positions_nm"][0][0] += 0.01
    elif mutation == "velocity":
        f["velocities_nm_ps"][0][0] += 0.01
    elif mutation == "box":
        f["box_vectors_nm"][0][0] += 0.01
    elif mutation == "checkpoint":
        f["checkpoint_sha256"] = "f" * 64
    elif mutation == "temperature":
        f["measured_temperature_kelvin"] = 300.0
    elif mutation == "pressure":
        f["measured_pressure_bar"] = 1.0
    elif mutation == "stage":
        f["stage"] = "nvt"
    elif mutation == "config":
        record["config"]["pressure_bar"] = 2.0
    elif mutation == "source":
        b = json.loads((tmp_path / "backend.json").read_text())
        b["source_identity"] = "f" * 64
        (tmp_path / "backend.json").write_text(json.dumps(b))
    else:
        record["production_validated"] = True
    record.pop("identity")
    (target / "failure.json").write_text(json.dumps(_signed(record)))
    with pytest.raises((IslandError, ValueError)):
        trajectory(tmp_path)


def test_offline_matrix_reconstruction(run):
    assert not validate_run(run)["scientific_acceptance"]


@pytest.mark.parametrize(
    "mutation", ("origin", "replica", "scientific_acceptance", "policy")
)
def test_matrix_mutations_rejected(run, tmp_path, mutation):
    # Symlink retained directories read-only; mutate only the new matrix.
    for p in run.iterdir():
        if p.name != "matrix.json":
            (tmp_path / p.name).symlink_to(p)
    matrix = copy.deepcopy(json.loads((run / "matrix.json").read_text()))
    if mutation == "origin":
        matrix["cases"]["opls-peo"][0]["charge_origin"] = "RESP/Psi4"
    elif mutation == "replica":
        matrix["cases"]["opls-peo"].pop()
    elif mutation == "scientific_acceptance":
        matrix["scientific_acceptance"] = True
    else:
        matrix["protocol"]["maximum_windows"] = 1
    matrix.pop("identity")
    (tmp_path / "matrix.json").write_text(json.dumps(_signed(matrix)))
    with pytest.raises(AnalysisError):
        validate_run(tmp_path)


def test_child_reconstruction_and_nonzero_scientific_gate(run):
    base = [
        sys.executable,
        "-m",
        "scripts.validate_phase_4j40_c2_receipt",
        str(run),
        "--bundle",
    ]
    for extra, code in ((["--diagnostic-only"], 0), ([], 1)):
        result = subprocess.run(
            [*base, *extra], capture_output=True, text=True, check=False
        )
        assert result.returncode == code, result.stdout + result.stderr


def test_unchanged_policy_and_scientific_boundaries():
    p = CorrectiveProtocol()
    assert p.replicas == 3
    assert p.maximum_windows * p.npt_window_steps == 2000000
    with pytest.raises(AnalysisError):
        require(False, "Pressure target is not an observed pressure")
