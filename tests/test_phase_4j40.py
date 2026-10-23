"""J40 independent diagnostics, false-success and persistence controls."""

import copy
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import numpy as np
import pytest

from island.analysis.equilibration import (
    EquilibrationPolicy,
    diagnose,
    series,
    validate_analysis,
)
from island.analysis.equilibration_run import HeatingSchedule
from island.analysis.npt import _body, _signed
from island.analysis.npt_continuation import load_continuation
from island.exceptions import AnalysisError, ValidationError


@pytest.mark.parametrize("values", [np.zeros(1024), np.arange(1024), np.arange(4)])
def test_insufficient_series(values):
    assert series(values, EquilibrationPolicy())["status"] == "insufficient_sampling"


def test_independent_white_noise():
    result = series(
        100 + np.random.default_rng(42).normal(0, 0.1, 8192), EquilibrationPolicy()
    )
    assert result["status"] == "ok"
    assert result["effective_sample_size"] > 1000
    assert result["standard_error"] > 0
    assert len(result["block_means"]) == 4


def test_correlated_noise_reduces_ess():
    rng = np.random.default_rng(20)
    x = np.zeros(4096)
    for i in range(1, len(x)):
        x[i] = 0.99 * x[i - 1] + rng.normal()
    assert series(x + 100, EquilibrationPolicy())["effective_sample_size"] < 100


@pytest.mark.parametrize(
    "kwargs",
    [
        {"burn_fraction": 0.8},
        {"blocks": 3},
        {"minimum_frames": True},
        {"minimum_ess": float("nan")},
        {"maximum_relative_drift": -1},
    ],
)
def test_invalid_policy(kwargs):
    with pytest.raises(AnalysisError):
        EquilibrationPolicy(**kwargs)


@pytest.mark.parametrize(
    "stages",
    [
        (),
        ((64, float("nan")),),
        ((64, 100), (32, 300)),
        ((64, 300), (128, 100)),
        ((500, 300),),
        ((64, 200),),
    ],
)
def test_invalid_heating(stages):
    with pytest.raises(ValidationError):
        HeatingSchedule(stages, 448, 300)


def test_heating_boundaries():
    schedule = HeatingSchedule(((64, 100), (128, 200), (192, 300)), 448, 300)
    assert [schedule.temperature(i) for i in (0, 64, 65, 129, 192, 448)] == [
        100,
        100,
        200,
        300,
        300,
        300,
    ]


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    declared = os.environ.get("ISLAND_J40_RUN")
    if declared:
        return Path(declared)
    root = tmp_path_factory.mktemp("j40")
    with tarfile.open("docs/evidence/phase_4j40/artifacts.tar.gz") as archive:
        archive.extractall(root, filter="data")
    return root / "run"


@pytest.fixture(scope="module")
def workflows(bundle):
    return [
        load_continuation(
            bundle / "opls-resp-aromatic" / f"replica-{i}" / "uninterrupted"
        )[0]
        for i in range(2)
    ]


def test_actual_diagnostics(workflows):
    record = diagnose(workflows)
    validate_analysis(record, workflows)
    assert not record["acceptance_met"]
    assert record["equilibration_status"] == "insufficient_sampling"
    assert record["replicas"][0]["metrics"]["pressure"]["status"] == "blocked"


@pytest.mark.parametrize("key", ["seed", "velocity_seed", "barostat_seed"])
def test_replica_independence(workflows, key):
    changed = copy.deepcopy(workflows)
    changed[1]["manifest"]["config"][key] = changed[0]["manifest"]["config"][key]
    with pytest.raises(AnalysisError, match="distinct"):
        diagnose(changed)


@pytest.mark.parametrize("mutation", ["equilibrium", "policy", "density", "trajectory"])
def test_fabricated_equilibration_rejected(workflows, mutation):
    record = _body(diagnose(workflows))
    if mutation == "equilibrium":
        record["acceptance_met"] = True
        record["equilibration_status"] = "converged_under_declared_finite_run_criteria"
    elif mutation == "policy":
        record["policy"]["minimum_ess"] = 1
    elif mutation == "density":
        record["replicas"][0]["metrics"]["density"]["mean"] += 1
    else:
        record["replicas"][0]["trajectory_identity"] = "0" * 64
    with pytest.raises(AnalysisError):
        validate_analysis(_signed(record), workflows)


def test_missing_npt(workflows):
    changed = copy.deepcopy(workflows)
    changed[0]["frames"] = [f for f in changed[0]["frames"] if f["stage"] != "npt"]
    with pytest.raises(AnalysisError, match="did not execute"):
        diagnose(changed)


def test_missing_replica(workflows):
    with pytest.raises(AnalysisError, match="replicas"):
        diagnose(workflows[:1])


def test_child_offline_validation(bundle):
    command = [
        sys.executable,
        "-m",
        "scripts.validate_phase_4j40_receipt",
        str(bundle),
        "--bundle",
        "--diagnostic-only",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["acceptance_met"] is False


def test_native_origins_and_cross_engine(bundle):
    matrix = json.loads((bundle / "matrix.json").read_text())
    assert (
        matrix["unsupported_origins"]["native_aromatic_oplsaa"] == "source_incompatible"
    )
    assert matrix["unsupported_origins"]["native_am1bcc"] == "blocked_precision"
    assert matrix["unsupported_origins"]["unsupported_pcff"] == "typing_rejection"
    for key, row in matrix["cases"].items():
        assert row["charge_origin"] == (
            "native_oplsaa" if key in ("opls-peo", "opls-network") else "RESP/Psi4"
        )
        reference = json.loads((bundle / key / "cross-engine.json").read_text())
        assert reference["lammps_exit_code"] == 0
        assert reference["executable_identity"] == matrix["engine_identity"]
