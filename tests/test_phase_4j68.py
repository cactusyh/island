"""Synthetic conic controls are separate from retained real HMDS evidence."""

import copy
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, digest, pin, signed, write
from island.pcff_qm_pilot.hmds_constrained import physical_problem, solve
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_local_baseline import (
    correction,
    derivative_check,
    extend,
    fixed_families,
    references,
    span_sensitivity,
)
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.hmds_span import (
    certificate,
    check,
    classification,
    projection,
    radius,
    relaxation_system,
    system,
)
from scripts.replay_phase_4j68 import BASE, BRANCH, parent_semantics
from scripts.validate_phase_4j68_receipt import (
    DECLARATIONS,
    decision_semantics,
    semantics,
)
from tests.test_phase_4j66 import synthetic

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "docs/evidence/phase_4j68"
PARENT = REPO / "docs/evidence/phase_4j67"


def reseal(record):
    return signed({k: v for k, v in record.items() if k != "identity"})


def changed(value):
    return (
        not value
        if type(value) is bool
        else value + 1
        if type(value) is int
        else "forged"
    )


@pytest.fixture(scope="module")
def model():
    return load_model(REPO)


@pytest.fixture(scope="module")
def rows():
    return read(PARENT / "fitting-dataset.json")["observations"]


def test_manufactured_bounds_only_obstruction():
    pytest.importorskip("clarabel")
    d = synthetic()
    d["QM_energy"] = 4 * d["energy"][:, 0] + 23
    d["QM_force"] = 4 * d["force"][..., 0]
    p = physical_problem(d)
    lo, hi = np.array([-0.1]), np.array([0.1])
    bounded = solve(system(p, lo, hi, bounded=True))
    free = solve(system(p, lo, hi))
    a = check(p, None, bounded["values"], lo, hi)
    b = check(p, None, free["values"], lo, hi)
    assert a["required_rho"] > 1
    assert b["performance_limits_satisfied"] and not b["original_bounds_satisfied"]
    assert len(b["outside_original_bounds"]) == 1
    assert (
        classification(b, {}) == "independently_verified_performance_feasible_witness"
    )
    assert free["values"] == pytest.approx([4, 23], abs=1e-6)


def test_manufactured_minimum_normalized_bound_relaxation():
    pytest.importorskip("clarabel")
    d = synthetic()
    d["energy"][:] = np.array([[-0.01], [0.01]])
    d["force"][:] = 1
    d["QM_energy"] = 2.5 * d["energy"][:, 0]
    d["QM_force"][:] = 2.5
    p = physical_problem(d)
    lo, hi = np.array([-1.0]), np.array([1.0])
    s = relaxation_system(p, lo, hi)
    r = solve(s)
    assert r["rho"] == pytest.approx(0.5, abs=1e-6)  # auxiliary eta
    assert r["values"][0] == pytest.approx(2.0, abs=1e-6)
    assert check(p, None, r["values"], lo, hi)["required_rho"] <= 1 + 1e-7
    assert abs(r["values"][0]) <= 1 + 2 * r["rho"] + 1e-7


def test_manufactured_affine_span_obstruction_checked_certificate():
    pytest.importorskip("clarabel")
    d = synthetic(n=3)
    d["force"][:] = 0
    d["QM_force"][:, 0, 0] = 3
    p = physical_problem(d)
    lo, hi = np.array([-1.0]), np.array([1.0])
    s = system(p, lo, hi)
    r = solve(s)
    c = certificate(s, r, p, None, lo, hi)
    assert c["certified"] and c["lower_bound"] > 1.499
    assert (
        classification(check(p, None, r["values"], lo, hi), c)
        == "checked_lower_bound_above_one"
    )
    assert c["region"]["includes_energy_only_intercept"]


def test_certificate_includes_all_unbounded_stationarity_residuals():
    pytest.importorskip("clarabel")
    d = synthetic(n=3)
    d["force"][:] = 0
    d["QM_force"][:, 0, 0] = 3
    p = physical_problem(d)
    lo, hi = np.array([-1.0]), np.array([1.0])
    s = system(p, lo, hi)
    raw = solve(s)
    good = certificate(s, raw, p, None, lo, hi)
    bad = copy.deepcopy(raw)
    bad["dual"][1] += 1  # energy individual constraint gives intercept stationarity
    penalized = certificate(s, bad, p, None, lo, hi)
    assert penalized["ball_support"] < good["ball_support"] - 1
    assert penalized["lower_bound"] <= raw["rho"]
    assert abs(penalized["stationarity_residual"][-2]) > 0.01


def test_checked_radius_and_unbounded_energy_intercept():
    d = synthetic()
    p = physical_problem(d)
    scale = np.linalg.norm(p["w"], axis=0)
    c = radius(p, scale, 2.0)
    singular = np.linalg.svd(p["w"] / scale, compute_uv=False)
    assert c["certified"] and 0 < c["sigma_lower"] <= singular[-1]
    assert c["left_inverse_defect_F"] < 1e-10
    assert c["scaled_coefficient_radius"] > 0
    x = np.array([0.01, 1e5])
    d["QM_energy"] = d["energy"][:, 0] * x[0] + x[1]
    d["QM_force"] = d["force"][..., 0] * x[0]
    p = physical_problem(d)
    result = projection(p)
    assert result["values"] == pytest.approx(x, abs=1e-8)
    assert np.all(p["f"][:, -1] == 0) and np.all(p["e"][:, -1] == 1)
    c = radius(p, np.linalg.norm(p["w"], axis=0), 1.0)
    assert (
        np.linalg.norm(x * np.linalg.norm(p["w"], axis=0))
        <= c["scaled_coefficient_radius"]
    )


def test_rank_deficient_radius_never_ignores_unbounded_nullspace():
    d = synthetic()
    d["energy"][:] = 0
    d["force"][:] = 0
    d["QM_force"][:] = 3
    p = physical_problem(d)
    c = radius(p, np.maximum(np.linalg.norm(p["w"], axis=0), 1e-30), 10)
    assert not c["certified"] and c["rank"] == 1
    assert (
        classification(
            check(p, None, [0, 0], np.array([-1.0]), np.array([1.0])),
            {"certified": False},
        )
        == "numerically_inconclusive"
    )


@pytest.mark.parametrize("rho", [1 - 1e-8, 1, 1 + 1e-8, 1.1])
def test_no_certificate_is_not_a_status_based_infeasibility_claim(rho):
    c = {"performance_limits_satisfied": rho <= 1, "required_rho": rho}
    outcome = classification(c, {"certified": False})
    assert outcome == (
        "independently_verified_performance_feasible_witness"
        if rho <= 1
        else "boundary_or_accuracy_limited"
        if rho < 1.00001
        else "numerically_inconclusive"
    )


def test_protected_subset_cannot_be_diluted_and_upper_uses_both():
    pytest.importorskip("clarabel")
    d = synthetic(n=26)
    d["energy"][:, 0] = np.linspace(-0.1, 0.1, 26)
    d["force"][:] = 0
    d["QM_force"][:22] = 0.51
    p, old = physical_problem(d), physical_problem(d, np.arange(22))
    lo, hi = np.array([-1.0]), np.array([1.0])
    full = solve(system(p, lo, hi))
    guarded = solve(system(p, lo, hi, old))
    assert full["rho"] < 1 and guarded["rho"] == pytest.approx(1.02, abs=1e-6)
    c = certificate(system(p, lo, hi, old), guarded, p, old, lo, hi)
    assert c["lower_bound"] > 1.0199
    assert c["primal_upper_bound"] >= max(
        v["required_rho"] for v in c["upper_bound_checks"]["subsets"].values()
    )


def test_local_vector_norm_and_global_component_rms():
    d = synthetic(n=2, atoms=3)
    d["force"][:] = 0
    d["QM_force"][0, 0] = [1.3, 1.3, 1.3]
    c = check(physical_problem(d), None, [0, 0], np.array([-1.0]), np.array([1.0]))[
        "subsets"
    ]["full26"]
    assert c["maximum_local_force_error"] == pytest.approx(np.sqrt(3) * 1.3)
    assert c["force_RMS"] == pytest.approx(np.sqrt(3 * 1.3**2 / 18))
    assert not c["performance_limits_satisfied"]


@pytest.mark.parametrize("factor", [1e-6, 1e6])
def test_physical_parameter_scaling(factor):
    pytest.importorskip("clarabel")
    d = synthetic()
    d["energy"] *= factor
    d["force"] *= factor
    d["QM_energy"] = d["energy"][:, 0] * (2 / factor) + 7
    d["QM_force"] = d["force"][..., 0] * (2 / factor)
    p = physical_problem(d)
    lo, hi = np.array([-1.0]), np.array([1.0])
    r = solve(system(p, lo, hi))
    assert r["values"] == pytest.approx([2 / factor, 7], rel=1e-6, abs=1e-7)
    assert check(p, None, r["values"], lo, hi)["performance_limits_satisfied"]


@pytest.mark.parametrize("variant", ["bridge_angle", "SiO_stretch"])
def test_scalar_correction_forces_units_and_sharing(model, rows, variant):
    reference = references(model, variant)
    assert len(reference) == (1 if variant == "bridge_angle" else 2)
    checks = derivative_check(model, rows[-1]["coordinates"], variant)
    assert checks["force_finite_difference_error"] < 1e-7
    assert checks["independent_Jet_force_error"] < 1e-12
    e, f = correction(model, rows[-1]["coordinates"], variant)
    assert e.shape == (2,) and f.shape == (27, 3, 2)
    assert np.max(abs(np.sum(f, axis=0))) < 1e-12
    h = 1e-5
    q = np.array(rows[-1]["coordinates"])
    perturbed = q.copy()
    perturbed[4, 0] += h
    ep, _ = correction(model, perturbed, variant)
    assert (ep - e) / h == pytest.approx(-f[4, 0], abs=5e-5)


def test_actual_fixed_family_decomposition_closes_and_reveals_cancellation(model, rows):
    d = fixed_families(model, rows[-1]["coordinates"])
    assert d["energy_sum_error"] < 1e-10 and d["force_sum_error"] < 1e-10
    assert {
        "quartic_bond:SiO",
        "quartic_angle:SiOSi",
        "LJ96",
        "Coulomb",
        "torsion_3",
        "bond-angle",
        "bond-bond",
        "angle-angle",
    } <= set(d["families"])
    vectors = np.array([r["forces"][4] for r in d["families"].values()])
    assert sum(np.linalg.norm(vectors, axis=1)) > np.linalg.norm(sum(vectors))


def test_duplicate_or_degenerate_added_columns_are_not_independent():
    d = synthetic()
    old = physical_problem(d)
    amended = {
        **d,
        "names": ["a", "duplicate", "zero"],
        "energy": np.c_[d["energy"], d["energy"], np.zeros(2)],
        "force": np.concatenate(
            [d["force"], d["force"], np.zeros_like(d["force"])], axis=-1
        ),
    }
    c = span_sensitivity(old, physical_problem(amended))
    assert not c["independent"] and c["reason"] == "zero added column"
    amended["force"][..., -1] = d["force"][..., 0]
    amended["energy"][:, -1] = d["energy"][:, 0]
    assert not span_sensitivity(old, physical_problem(amended))["independent"]


def test_real_hypotheses_add_independent_columns(model, rows):
    from island.pcff_qm_pilot.hmds_fit import features

    d = features(model, rows)
    p = physical_problem(d)
    for variant in ("bridge_angle", "SiO_stretch"):
        s = span_sensitivity(p, physical_problem(extend(model, rows, d, variant)))
        assert s["independent"] and min(s["outside_old_span_norms"]) > 0.04


@pytest.mark.parametrize("key", list(DECLARATIONS))
def test_resigned_parent_and_new_receipt_mutations(key):
    parent, manifest, decision = (
        read(PARENT / n) for n in ("receipt.json", "manifest.json", "decision.json")
    )
    parent_semantics(parent, manifest, decision)
    if key in parent:
        parent[key] = changed(parent[key])
        with pytest.raises(CompletionError):
            parent_semantics(reseal(parent), manifest, decision)
    # Construct a semantically valid J68 receipt without relying on a later seal.
    names = (
        "bounded-Clarabel",
        "bounded-SCS",
        "bound-free-Clarabel",
        "bound-free-SCS",
        "bridge_angle-clarabel",
        "bridge_angle-scs",
        "SiO_stretch-clarabel",
        "SiO_stretch-scs",
    )
    receipt = {
        **DECLARATIONS,
        "base_commit": BASE,
        "required_J67": BASE,
        "branch": BRANCH,
        "manifest_sha256": "placeholder",
        "historical_files_preserved": 1,
        "software_counts": {},
        "terminal_unique_tests": 1,
        "focused_J68": 1,
        "historical_J54_J67": 0,
        "decision_identity": "placeholder",
        "numerical_reconstruction": {
            "observations": 26,
            "protected": 22,
            "physical_columns": 58,
            "intercept_columns": 1,
            "research_solver_calls": 8,
            "failed_processes_retained": 1,
            "new_QM": 0,
            "MD_steps": 0,
            "validation_untouched": True,
            "accepted_coefficients": 0,
            "offline": True,
            "classifications": dict.fromkeys(names, "checked_lower_bound_above_one"),
        },
    }
    m = {
        "schema": "island.numerical-evidence-manifest.v1",
        "base_commit": BASE,
        "historical_files": {"x": "y"},
        "implementation_files": {"x": "y"},
        "artifacts": {"x": "y"},
    }
    semantics(signed(receipt), signed(m))
    receipt[key] = changed(receipt[key])
    altered = signed(receipt)
    assert altered["identity"] == digest(receipt)
    with pytest.raises(CompletionError):
        semantics(altered, signed(m))


@pytest.mark.parametrize("key", [k for k in DECLARATIONS if k != "schema"])
def test_resigned_decision_linked_through_manifest_and_receipt_rejected(key, tmp_path):
    d = read(ROOT / "decision.json")
    decision_semantics(d)
    d[key] = changed(d[key])
    d = reseal(d)
    write(tmp_path / "decision.json", d)
    manifest = signed(
        {
            "schema": "island.numerical-evidence-manifest.v1",
            "artifacts": {"decision.json": pin(tmp_path / "decision.json")},
        }
    )
    write(tmp_path / "manifest.json", manifest)
    receipt = signed({"manifest_sha256": pin(tmp_path / "manifest.json")})
    assert receipt["manifest_sha256"] == pin(tmp_path / "manifest.json")
    assert manifest["artifacts"]["decision.json"] == pin(tmp_path / "decision.json")
    with pytest.raises(CompletionError):
        decision_semantics(d)


@pytest.mark.parametrize(
    "mutation", ["schema", "extra_claim", "nested_count", "parent_decision"]
)
def test_resigned_parent_manifest_and_nested_semantics(mutation):
    r, m, d = (
        read(PARENT / n) for n in ("receipt.json", "manifest.json", "decision.json")
    )
    if mutation == "schema":
        m["schema"] = "forged"
    elif mutation == "extra_claim":
        m["scientific_validation"] = "passed"
    elif mutation == "nested_count":
        r["numerical_reconstruction"]["new_gradients"] = 0
    else:
        d["training_feasible"] = True
    with pytest.raises(CompletionError):
        parent_semantics(reseal(r), reseal(m), reseal(d))


def test_preserved_observation_producers_and_reserved_validation():
    dataset = read(PARENT / "fitting-dataset.json")["observations"]
    old = read(REPO / "docs/evidence/phase_4j65/fitting-dataset.json")["observations"]
    assert dataset[:22] == old and len(dataset) == 26
    assert all(r["producer"]["kind"] == "new_J67_execution" for r in dataset[22:])
    d = read(ROOT / "decision.json")
    assert not d["validation_activated"] and not d["validation_consumed"]
    assert (
        d["new_QM"]
        == d["MD_steps"]
        == d["validation_QM"]
        == d["accepted_coefficients"]
        == 0
    )
    assert not (ROOT / "validation-consumption.json").exists()


def test_failure_and_raw_solver_state_retained():
    from xml.etree import ElementTree as ET

    from scripts.validate_phase_4j68_receipt import selected_cases

    original = ET.fromstring(
        '<testsuite><testcase classname="tests.old" name="retained"/><testcase classname="tests.test_phase_4j68" name="corrected"><failure/></testcase></testsuite>'
    )
    selected = ET.fromstring(
        '<testsuite><testcase classname="tests.old" name="retained"/></testsuite>'
    )
    replacement = ET.fromstring(
        '<testsuite><testcase classname="tests.test_phase_4j68" name="corrected"/></testsuite>'
    )
    assert len(selected_cases(original, selected, replacement)) == 1
    bad = copy.deepcopy(selected)
    bad.find("testcase").set("name", "rewritten")
    with pytest.raises(CompletionError):
        selected_cases(original, bad, replacement)
    replacement.find("testcase").append(ET.Element("failure"))
    with pytest.raises(CompletionError):
        selected_cases(original, selected, replacement)
    from island.pcff_qm_pilot.hmds_span import authorize_attempt

    rows = [
        {"name": "trial", "backend": "clarabel", "bounded": False, "status": "started"}
    ]
    assert authorize_attempt(rows, "trial", "clarabel", False)
    with pytest.raises(CompletionError):
        authorize_attempt(rows, "unrecorded", "clarabel", False)
    with pytest.raises(CompletionError):
        authorize_attempt(rows * 13, "trial", "clarabel", False)
    with pytest.raises(CompletionError):
        authorize_attempt(rows, "trial", "scs", False)
    calls = read(ROOT / "solver-attempts.json")["attempts"]
    assert len(calls) == 8
    failed = [r for r in calls if r["status"] == "failed"]
    assert [r["name"] for r in failed] == ["bridge_angle-clarabel"]
    assert not any(r["name"].endswith("retry") for r in calls)
    assert read(ROOT / "execution-environment.json")["actual_solver_successes"] == 8


def test_fresh_process_relocated_evidence_offline_replay(tmp_path):
    relocated = tmp_path / "relocated"
    shutil.copytree(
        ROOT,
        relocated,
        ignore=shutil.ignore_patterns(
            "validation", "manifest.json", "receipt.json", "solves"
        ),
    )
    command = """import importlib.abc,runpy,sys
class Deny(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname.split('.')[0] in {'psi4','clarabel','scs'}: raise AssertionError('offline replay imported native/solver')
sys.meta_path.insert(0,Deny())
repo,evidence=sys.argv[1:3]
sys.argv=['replay','--repo',repo,'--evidence',evidence]
runpy.run_path(repo+'/scripts/replay_phase_4j68.py',run_name='__main__')
"""
    result = subprocess.run(
        [sys.executable, "-c", command, str(REPO), str(relocated)],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHONPATH": str(REPO) + ":" + str(REPO / "src"),
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "mutation",
    [
        "graph",
        "atom_order",
        "protocol",
        "producer",
        "matrix",
        "units",
        "bounds",
        "parameter_order",
    ],
)
def test_semantic_problem_mutations_with_valid_new_outer_identity(
    model, rows, mutation
):
    from island.pcff_qm_pilot.hmds_span import intake

    with np.load(PARENT / "physical-arrays.npz", allow_pickle=False) as saved:
        original = {
            k: (saved[k].tolist() if k == "names" else saved[k].copy())
            for k in saved.files
        }
    spec = read(PARENT / "plan/scientific-plan.json")["parameter_specifications"]
    r, d, s = copy.deepcopy(rows), copy.deepcopy(original), copy.deepcopy(spec)
    if mutation == "graph":
        r[0]["original_request"]["graph_identity"] = "0" * 64
    elif mutation == "atom_order":
        r[0]["original_request"]["atom_order"][:2] = [2, 1]
    elif mutation == "protocol":
        r[0]["original_request"]["grid"] = [150, 974]
    elif mutation == "producer":
        r[0]["producer"]["kind"] = "fabricated_new_execution"
    elif mutation == "matrix":
        d["force"][0, 4, 0, 0] += 0.01
    elif mutation == "units":
        s[0]["units"] = "kJ/mol"
    elif mutation == "bounds":
        s[0]["lower"] -= 1
    else:
        d["names"][0:2] = d["names"][1::-1]
    new_outer = signed(
        {
            "observations": r,
            "specifications": s,
            "matrix_digest": digest(
                {
                    k: v.tolist() if isinstance(v, np.ndarray) else v
                    for k, v in d.items()
                }
            ),
        }
    )
    assert new_outer["identity"] == digest(
        {k: v for k, v in new_outer.items() if k != "identity"}
    )
    with pytest.raises(CompletionError):
        intake(model, r, d, s, rows, original, spec)
