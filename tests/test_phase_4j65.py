"""Protocol intake, bounded coupled fitting and failure/uncertainty controls."""

import copy
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, digest, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hmds_baseline import extended_parts
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_fit import (
    accepts,
    features,
    infer,
    performance,
    prediction_uncertainty,
    predictions,
    uncertainty,
    weighted_problem,
)
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_paths import equivalent, leakage
from island.pcff_qm_pilot.hmds_precision import protocol_contrast, read
from scripts.analyze_phase_4j63 import sample
from scripts.analyze_phase_4j65 import observation, validation_preflight
from scripts.phase_4j64 import run_job

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "docs/evidence/phase_4j65"


@pytest.fixture(scope="module")
def plan():
    return read(ROOT / "plan/scientific-plan.json")


@pytest.fixture(scope="module")
def manufactured(plan):
    model = load_model(REPO)
    rows = plan["primary"] + plan["joints"]
    records = []
    values = np.r_[np.random.default_rng(65).uniform(-2, 2, 56), 23.4, 18.7, 137.2]
    for row in rows:
        p = extended_parts(model, row["coordinates"])
        records.append(
            {
                "id": row["id"],
                "block": row["block"],
                "coordinates": row["coordinates"],
                "forces": p["fixed_forces"] + p["force_jacobian"] @ values[:-1],
                "relative_energy_kcal_mol": p["fixed_energy"]
                + p["energy_jacobian"] @ values[:-1]
                + values[-1],
            }
        )
    return features(model, records), records, values


def test_exact_reuse_inventory_and_original_producers(plan, tmp_path):
    assert len(plan["reuse"]) == 7
    assert len(plan["primary"]) == 17
    row = next(r for r in plan["primary"] if r["id"] == "center")
    item = observation(REPO, tmp_path, plan, row)
    old = read(REPO / "docs/evidence/phase_4j64/normalized-reference.json")["samples"][
        "exact/center"
    ]
    assert item["original_normalized_sample"] == old
    assert item["producer"]["kind"] == "reused_original_execution"
    assert item["producer"]["new_native_attempt"] is False
    assert item["relative_energy_kcal_mol"] == 0


@pytest.mark.parametrize(
    "change", ["basis", "grid", "thresholds", "auxiliary", "coordinates"]
)
def test_protocol_reuse_rejects_mixed_observations(plan, tmp_path, change):
    other = copy.deepcopy(plan)
    row = copy.deepcopy(plan["primary"][0])
    if change == "basis":
        other["protocol"]["basis"] = "def2-qzvp"
    elif change == "grid":
        other["protocol"]["grid"] = [150, 974]
    elif change == "thresholds":
        other["protocol"]["thresholds"] = [1e-10, 1e-8]
    elif change == "auxiliary":
        other["protocol"]["settings"]["df_basis_scf"] = "cc-pvtz-jkfit"
    else:
        row["coordinates"][0][0] += 0.001
    with pytest.raises(CompletionError, match="exact_geometry_and_protocol"):
        observation(REPO, tmp_path, other, row)


def test_cross_site_predictions_frozen_without_parameter_refitting(plan):
    for row in plan["cross_site_predictions"]:
        assert row["predicted_QZVP_total_curvature"] == pytest.approx(
            row["fixed_total_curvature"] + 2 * row["fixed_correction"]
        )
        assert row["correction_refitted"] is False
        assert row["same_site_path_equivalent_at1e_3A"] is True
        assert row["same_site_path_equivalent_at1e_6A"] is False
        assert len(row["unknown_null"]["unknown_columns"]) == 56


def test_qzvp_grid_contrast_keeps_scf_basis_and_auxiliary(plan):
    original = request(
        plan,
        plan["center_coordinates"],
        basis="def2-qzvp",
        grid=[150, 974],
        thresholds=[1e-10, 1e-8],
    )
    dense = read(ROOT / "plan/requests.json")["controls/QZVP-dense-center"]
    assert protocol_contrast(original, dense, ["grid"])["changed_fields"] == ["grid"]
    bad = copy.deepcopy(dense)
    bad["settings"]["df_basis_scf"] = "other"
    with pytest.raises(CompletionError, match="confounded"):
        protocol_contrast(original, bad, ["grid"])


def test_manufactured_all58_recovery_with_nonzero_offset(manufactured, plan):
    data, rows, truth = manufactured
    result = infer(data, plan["parameter_specifications"], plan["fit"])
    assert result["solver_success"] and result["identifiable"]
    assert not result["active_bound_parameters"]
    assert np.max(abs(np.array(result["values"]) - truth)) < 1e-6
    report = performance(data, result["values"], [r["id"] for r in rows])
    assert report["energy_RMSE"] < 1e-8 and report["force_RMS"] < 1e-8
    assert result["actual_weighted_design"]["columns"] == 59
    assert result["actual_weighted_design"]["scaled_condition"] > 1000


def test_objective_uses_actual_weights_and_offset(manufactured, plan):
    data, _, truth = manufactured
    a, y = weighted_problem(data, plan["fit"])
    assert a.shape == (22 * 82, 59)
    assert np.linalg.norm(a @ truth - y) < 1e-9
    assert np.all(a[22:, -1] == 0)
    assert np.all(a[:22, -1] > 0)
    shift = copy.deepcopy(data)
    shift["QM_energy"] += 2.6
    result = infer(shift, plan["parameter_specifications"], plan["fit"])
    assert np.max(abs(np.array(result["values"][:-1]) - truth[:-1])) < 1e-6
    assert result["values"][-1] == pytest.approx(truth[-1] + 2.6, abs=1e-7)


def test_validation_offset_never_refits(manufactured, plan):
    data, rows, truth = manufactured
    validation = copy.deepcopy(data)
    validation["QM_energy"] += 0.4
    report = performance(validation, truth, [r["id"] for r in rows])
    assert report["energy_RMSE"] == pytest.approx(0.4)
    assert report["force_RMS"] < 1e-12
    assert not report["validation_offsets_refitted"]
    assert not accepts(report, plan["fit"])


def test_bounds_are_enforced_and_activity_visible(manufactured, plan):
    data, _, truth = manufactured
    forced = copy.deepcopy(data)
    outside = truth.copy()
    outside[-3] = 2 * plan["parameter_specifications"][-2]["upper"]
    e, f = predictions(data, outside)
    forced["QM_energy"], forced["QM_force"] = e, f
    result = infer(forced, plan["parameter_specifications"], plan["fit"])
    assert result["active_bound_parameters"]
    for s, v in zip(
        plan["parameter_specifications"], result["values"][:-1], strict=True
    ):
        assert s["lower"] - 1e-8 <= v <= s["upper"] + 1e-8


def test_rank_loss_is_not_optimizer_or_scientific_success(manufactured, plan):
    data, _, _ = manufactured
    result = infer(data, plan["parameter_specifications"], plan["fit"], [0, 0])
    assert result["actual_weighted_design"]["rank"] < 59
    assert not result["identifiable"]


def test_bootstrap_failures_retained_without_replacement(manufactured, plan):
    data, rows, truth = manufactured
    one = {k: (v[:1] if isinstance(v, np.ndarray) else v) for k, v in data.items()}
    settings = copy.deepcopy(plan["fit"])
    settings["uncertainty"]["replicates"] = 4
    result = uncertainty(
        one, rows[:1], truth, plan["parameter_specifications"], settings
    )
    assert len(result["attempts"]) == 4
    assert result["eligible_replicates"] == 0 and not result["qualified"]
    assert (
        result["failed_replicates_retained"] == 4
        and result["replacement_resamples"] == 0
    )
    assert result["parameter_intervals"] is None


@pytest.mark.parametrize("bad", [None, np.zeros(58), np.r_[np.zeros(58), np.nan]])
def test_incomplete_fitted_vector_rejected(manufactured, bad):
    with pytest.raises((CompletionError, TypeError)):
        predictions(manufactured[0], bad)


def test_missing_force_vector_rejected(plan):
    row = {
        "coordinates": plan["center_coordinates"],
        "forces": np.zeros((26, 3)),
        "relative_energy_kcal_mol": 0,
    }
    with pytest.raises(CompletionError, match="finite_complete"):
        features(load_model(REPO), [row])


def test_reserved_validation_leakage_and_nonactivation(plan, tmp_path):
    model = load_model(REPO)
    geo = np.array(plan["validation_coordinates"][0])
    with pytest.raises(CompletionError):
        leakage(model, [geo], [geo + np.array([1, 2, 3])])
    assert equivalent(model, geo, geo + np.array([1, 2, 3]))
    write(tmp_path / "plan/scientific-plan.json", plan)
    write(tmp_path / "trial-fit.json", {"validation_authorized": False})
    with pytest.raises(CompletionError, match="untouched_validation_not_activated"):
        validation_preflight(REPO, tmp_path, tmp_path / "work")
    assert not (tmp_path / "validation-activation.json").exists()


def test_startup_failure_retains_process_diagnostics_and_counts(tmp_path, monkeypatch):
    root = tmp_path / "evidence"
    work = tmp_path / "work"
    work.mkdir()
    plan = {
        "execution_sources": {},
        "limits": {"maximum_new_gradients": 2, "campaign_wall_seconds": 36000},
    }
    for d in (root / "plan", work):
        write(d / "scientific-plan.json", plan)
    req = {"test": "failure-only manufactured request, never a QM result"}
    write(work / "requests/fail.json", req)
    write(
        work / "jobs.json",
        {
            "jobs": {
                "fail": {
                    "memory_GiB": 4,
                    "wall_seconds": 1,
                    "request_identity": digest(req),
                }
            }
        },
    )

    def fail(*args, **kwargs):
        raise OSError("manufactured executable absence")

    monkeypatch.setattr(subprocess, "Popen", fail)
    assert run_job(REPO, root, work, Path("missing-interpreter"), "fail") == 127
    result = read(work / "processes/fail/result.json")
    assert result["status"] == "unexpected_failure"
    assert (
        "manufactured executable absence"
        in (work / "processes/fail/stderr.log").read_text()
    )
    assert (work / "processes/fail/command.json").exists()
    assert len(read(work / "attempts.json")["calls"]) == 1


def test_process_cap_prevents_launch(tmp_path, monkeypatch):
    root = tmp_path / "evidence"
    work = tmp_path / "work"
    work.mkdir()
    plan = {
        "execution_sources": {},
        "limits": {"maximum_new_gradients": 32, "campaign_wall_seconds": 36000},
    }
    for d in (root / "plan", work):
        write(d / "scientific-plan.json", plan)
    write(work / "requests/new.json", {})
    write(work / "jobs.json", {"jobs": {"new": {"request_identity": digest({})}}})
    write(
        work / "attempts.json",
        {"started_unix": time.time(), "calls": {str(i): {} for i in range(32)}},
    )

    def forbidden(*args, **kwargs):
        raise AssertionError("launched after cap")

    monkeypatch.setattr(subprocess, "Popen", forbidden)
    with pytest.raises(CompletionError, match="budget_exhausted"):
        run_job(REPO, root, work, Path("unused"), "new")
    assert not (work / "processes/new").exists()


def test_bound_name_mapping_mutation_rejected(manufactured, plan):
    bad = copy.deepcopy(plan["parameter_specifications"])
    bad[0], bad[1] = bad[1], bad[0]
    with pytest.raises(CompletionError, match="fit_bound_mapping"):
        infer(manufactured[0], bad, plan["fit"])


def test_absolute_prediction_uncertainty_for_near_zero_coefficient(manufactured):
    data, rows, truth = manufactured
    nominal = truth.copy()
    nominal[0] = 0
    plus, minus = nominal.copy(), nominal.copy()
    plus[0] = 0.2
    minus[0] = -0.2
    result = prediction_uncertainty(
        data, [plus, minus], nominal, [r["id"] for r in rows]
    )
    assert np.isfinite(result["energy_RMS95"]) and result["energy_RMS95"] > 0
    assert result["force_RMS95"] > 0
    assert not result["includes_model_discrepancy_or_basis_uncertainty"]
    with pytest.raises(CompletionError, match="two_bootstrap"):
        prediction_uncertainty(data, [plus], nominal, [r["id"] for r in rows])


def test_native_gradient_hash_mutation_rejected(tmp_path):
    work = tmp_path / "native"
    unpack(REPO / "docs/evidence/phase_4j64/native-outputs.tar.xz", work)
    file = work / "outputs/exact/center/gradient.npy"
    g = np.load(file, allow_pickle=False)
    g[0, 0] += 0.01
    np.save(file, g, allow_pickle=False)
    with pytest.raises(CompletionError):
        sample(work, "exact/center")


def test_actual_native_reconstruction_and_fresh_relocated_process(tmp_path):
    if not (ROOT / "native-outputs.tar.xz").exists():
        pytest.skip("Executed J65 campaign archive not yet retained")
    relocated = tmp_path / "relocated-repo"
    for name in [
        "docs/evidence/phase_4j54/audit/hmds/report.json",
        "docs/evidence/phase_4j60/hamiltonian-definition.json",
        "docs/evidence/phase_4j64/native-outputs.tar.xz",
        "docs/evidence/phase_4j64/baseline-calibration.json",
    ]:
        destination = relocated / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, destination)
    target = relocated / "docs/evidence/phase_4j65"
    shutil.copytree(ROOT, target)
    command = [
        str(Path(os.sys.executable)),
        str(REPO / "scripts/replay_phase_4j65.py"),
        "replay",
        "--repo",
        str(relocated),
        "--evidence",
        str(target),
    ]
    env = {
        **os.environ,
        "PYTHONPATH": f"{REPO}:{REPO / 'src'}",
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
    }
    result = subprocess.run(
        command,
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (
        "'QM_invoked': False" in result.stdout
        and "'exact_replay': True" in result.stdout
    )
