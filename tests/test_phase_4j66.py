"""Synthetic conic controls are separate from retained actual HMDS evidence."""

import copy
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, write
from island.pcff_qm_pilot.hmds_constrained import (
    canonical,
    classify,
    dual_lower_bound,
    physical_problem,
    solve,
    verify,
)
from island.pcff_qm_pilot.hmds_design import metrics
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j66 import reconstruct
from scripts.fit_phase_4j66 import attempt
from scripts.phase_4j66 import PARENTS
from scripts.replay_phase_4j66 import compare_system

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "docs/evidence/phase_4j66"


def synthetic(n=2, atoms=27):
    rng = np.random.default_rng(66)
    return {
        "names": ["manufactured"],
        "energy": rng.normal(size=(n, 1)),
        "force": rng.normal(size=(n, atoms, 3, 1)),
        "fixed_energy": np.zeros(n),
        "fixed_force": np.zeros((n, atoms, 3)),
        "QM_energy": np.zeros(n),
        "QM_force": np.zeros((n, atoms, 3)),
    }


@pytest.fixture
def backend():
    return pytest.importorskip("clarabel")


def test_manufactured_feasible_recovery(backend):
    data = synthetic()
    data["QM_energy"] = data["energy"][:, 0] * 1.3 + 7.4
    data["QM_force"] = data["force"][..., 0] * 1.3
    p = physical_problem(data)
    lo = np.array([-5.0])
    hi = -lo
    result = solve(canonical(p, lo, hi))
    check = verify(p, result["values"], lo, hi)
    assert check["original_limits_satisfied"] and check["required_rho"] < 1e-7
    assert result["values"] == pytest.approx([1.3, 7.4], abs=1e-7)


def test_LS_local_failure_can_be_constrained_feasible(backend, tmp_path):
    data = synthetic(n=1)
    data["energy"][:] = 0
    data["force"][:] = 0
    data["force"][0, 0, 0, 0] = 1
    data["force"][0, 1:11, 0, 0] = 0.4
    data["QM_force"][0, 0, 0] = 4.5
    p = physical_problem(data)
    lo = np.array([-10.0])
    hi = -lo
    least = np.linalg.lstsq(p["w"], p["y"], rcond=1e-10)[0]
    assert verify(p, least, lo, hi)["maximum_local_force_error"] > 2
    witness = solve(canonical(p, lo, hi))
    assert verify(p, witness["values"], lo, hi)["original_limits_satisfied"]
    constrained = attempt(
        p,
        lo,
        hi,
        tmp_path,
        "synthetic-LS",
        kind="least_squares",
        witness=witness["values"],
    )
    assert constrained["original_unit_verification"]["original_limits_satisfied"]
    assert constrained["solution"]["values"][0] == pytest.approx(2.5, abs=1e-5)


def test_manufactured_incompatible_checked_lower_bound(backend):
    data = synthetic()
    data["force"][:] = 0
    data["QM_force"][:, 0, 0] = 3
    p = physical_problem(data)
    lo = np.array([-2.0])
    hi = -lo
    system = canonical(p, lo, hi)
    result = solve(system)
    check = verify(p, result["values"], lo, hi)
    bound = dual_lower_bound(system, result, p, lo, hi)
    assert bound["lower_bound"] > 1.49
    assert (
        classify(result, check, bound["lower_bound"]) == "checked_lower_bound_above_one"
    )
    assert not check["original_limits_satisfied"]


@pytest.mark.parametrize(
    "status",
    ["MaxIterations", "NumericalError", "AlmostSolved", "infeasible_inaccurate"],
)
def test_inconclusive_status_is_not_infeasibility(status):
    assert (
        classify({"status": status, "rho": 1.5}, None) == "solver_failure_inconclusive"
    )
    assert (
        classify({"status": status}, {"original_limits_satisfied": True})
        == "verified_feasible_witness"
    )


def test_vector_norm_and_global_component_normalization():
    data = synthetic(n=2, atoms=3)
    data["force"][:] = 0
    data["QM_force"][:] = 0
    data["QM_force"][0, 0] = [1.3, 1.3, 1.3]
    p = physical_problem(data)
    v = verify(p, [0, 0], np.array([-1.0]), np.array([1.0]))
    assert v["maximum_local_force_error"] == pytest.approx(np.sqrt(3) * 1.3)
    assert v["force_RMS"] == pytest.approx(np.sqrt(3 * 1.3**2 / 18))
    assert not v["original_limits_satisfied"]  # components individually below2


def test_energy_intercept_affects_no_forces_and_no_per_config_refit():
    d = synthetic()
    p = physical_problem(d)
    assert np.all(p["f"][:, -1] == 0) and np.all(p["e"][:, -1] == 1)
    a = verify(p, [0, 0], np.array([-2.0]), np.array([2.0]))
    b = verify(p, [0, 1], np.array([-2.0]), np.array([2.0]))
    assert b["energy_RMSE"] == 1 and a["force_RMS"] == b["force_RMS"] == 0


def test_scaling_and_original_units(backend):
    d = synthetic()
    d["energy"] *= 1e5
    d["force"] *= 1e5
    d["QM_energy"] = d["energy"][:, 0] * 2e-5 + 3
    d["QM_force"] = d["force"][..., 0] * 2e-5
    p = physical_problem(d)
    lo = np.array([-0.001])
    hi = -lo
    s = canonical(p, lo, hi)
    assert s["scale"][0] > 1e5
    result = solve(s)
    v = verify(p, result["values"], lo, hi)
    assert v["original_limits_satisfied"] and result["values"] == pytest.approx(
        [2e-5, 3], abs=1e-6
    )


def test_bounds_and_interior_feasibility(backend):
    d = synthetic()
    d["QM_energy"] = d["energy"][:, 0]
    d["QM_force"] = d["force"][..., 0]
    p = physical_problem(d)
    lo = np.array([-1.0])
    hi = np.array([1.0])
    v = verify(p, [1, 0], lo, hi)
    assert v["active_bound_indices"] == [0]
    r = solve(canonical(p, lo, hi, interior=1e-4))
    assert r["values"][0] <= 0.99980001
    assert verify(p, r["values"], lo, hi)["active_bound_indices"] == []


def test_rank_is_distinct_from_conditioning():
    x = np.arange(30.0).reshape(10, 3)
    assert metrics(x)["rank"] < 3
    a = np.eye(3)
    a[:, 1] = a[:, 0] + 1e-5 * a[:, 1]
    assert metrics(a)["rank"] == 3 and metrics(a)["scaled_condition"] > 10000


def test_strict_physical_thresholds_not_widened():
    d = synthetic()
    d["force"][:] = 0
    d["QM_force"][:] = 0
    d["QM_force"][:, 0, 0] = 2 + 1e-11
    assert not verify(physical_problem(d), [0, 0], np.array([-1.0]), np.array([1.0]))[
        "original_limits_satisfied"
    ]


def test_actual_64_failures_preserved_and_uncertainty_unqualified():
    old = read(REPO / "docs/evidence/phase_4j65/block-uncertainty.json")
    new = read(ROOT / "block-uncertainty.json")
    assert len(new["attempts"]) == 64 and new["replacement_resamples"] == 0
    assert [r["sampled_blocks"] for r in new["attempts"]] == [
        r["sampled_blocks"] for r in old["attempts"]
    ]
    assert new["eligible_replicates"] < 0.8 * 64 and not new["qualified"]
    assert (
        sum(not r["eligible"] for r in new["attempts"])
        == 64 - new["eligible_replicates"]
    )
    assert "survivor-conditional" in new["scope"]


@pytest.mark.parametrize(
    "mutation",
    [
        "graph",
        "protocol",
        "bounds",
        "units",
        "parameter_order",
        "configuration_order",
        "atom_order",
    ],
)
def test_contract_mutations_rejected(tmp_path, mutation):
    plan = copy.deepcopy(read(ROOT / "plan/scientific-plan.json"))
    if mutation == "graph":
        plan["graph_identity"] = "bad"
    elif mutation == "protocol":
        plan["protocol"]["grid"] = [150, 974]
    elif mutation in ("bounds", "units"):
        plan["parameter_specifications"][0][
            "lower" if mutation == "bounds" else "units"
        ] = 0
    elif mutation == "parameter_order":
        plan["parameter_specifications"].reverse()
    elif mutation == "configuration_order":
        plan["configuration_order"].reverse()
    else:
        plan["atom_order"].reverse()
    write(tmp_path / "plan/scientific-plan.json", plan)
    with pytest.raises(CompletionError):
        reconstruct(REPO, tmp_path)


def test_matrix_mutation_rejected(tmp_path):
    with np.load(ROOT / "nominal-rho-canonical.npz") as v:
        stored = {k: v[k] for k in v.files}
    stored["A_data"] = stored["A_data"].copy()
    stored["A_data"][0] += 0.01
    np.savez_compressed(tmp_path / "mutated.npz", **stored)
    with np.load(ROOT / "physical-arrays.npz") as v:
        data = {k: v[k] for k in v.files}
    data["names"] = data["names"].tolist()
    plan = read(ROOT / "plan/scientific-plan.json")
    lo = np.array([s["lower"] for s in plan["parameter_specifications"]])
    hi = np.array([s["upper"] for s in plan["parameter_specifications"]])
    with pytest.raises(CompletionError):
        compare_system(
            tmp_path / "mutated.npz", canonical(physical_problem(data), lo, hi)
        )


def test_incomplete_vectors_and_nonfinite_arrays_rejected():
    d = synthetic()
    p = physical_problem(d)
    with pytest.raises(CompletionError):
        verify(p, [0], np.array([-1.0]), np.array([1.0]))
    d["QM_force"][0, 0, 0] = np.nan
    with pytest.raises(CompletionError):
        physical_problem(d)


def test_reserved_validation_not_activated():
    d = read(ROOT / "decision.json")
    assert d["validation_QM"] == d["new_QM"] == d["MD_steps"] == 0
    assert not d["automatic_assignment"] and not d["accepted_candidates"]
    assert not (ROOT / "untouched-validation.json").exists()


def test_actual_fresh_process_relocation_and_offline_replay(tmp_path):
    relocated = tmp_path / "repo"
    for name in (
        "docs/evidence/phase_4j54/audit/hmds/report.json",
        "docs/evidence/phase_4j60/hamiltonian-definition.json",
        "docs/evidence/phase_4j64/baseline-calibration.json",
    ):
        target = relocated / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, target)
    parent = relocated / "docs/evidence/phase_4j65"
    parent.mkdir(parents=True)
    for name in PARENTS:
        target = parent / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / "docs/evidence/phase_4j65" / name, target)
    target = relocated / "docs/evidence/phase_4j66"
    shutil.copytree(ROOT, target)
    r = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/replay_phase_4j66.py"),
            "--repo",
            str(relocated),
            "--evidence",
            str(target),
            "--skip-native",
        ],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHONPATH": f"{REPO}:{REPO / 'src'}",
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
        },
        check=False,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    assert "'offline': True" in r.stdout and "'new_QM': 0" in r.stdout
