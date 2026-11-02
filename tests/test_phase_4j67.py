"""Synthetic solver controls and additive retained-evidence mutation tests.

Synthetic data never substitute for the four actual native QM observations.
"""

import copy
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, digest, pin, signed, write
from island.pcff_qm_pilot.hmds_augmentation import (
    PARENT_DECLARATIONS,
    SELECTED,
    bootstrap_draws,
    candidate_gates,
    exact_intake,
    exact_request,
    parent_semantics,
    physically_feasible,
    preserved_producers,
    protected_system,
    solver_qualified,
    validation_activation,
    verified_subsets,
)
from island.pcff_qm_pilot.hmds_constrained import physical_problem, solve, verify
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_fit import features
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_paths import leakage
from island.pcff_qm_pilot.hmds_precision import read
from tests.test_phase_4j66 import synthetic

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "docs/evidence/phase_4j67"
PARENT = REPO / "docs/evidence/phase_4j66"


@pytest.fixture(scope="module")
def model():
    return load_model(REPO)


@pytest.fixture(scope="module")
def plan():
    return read(ROOT / "plan/scientific-plan.json")


@pytest.mark.parametrize("key", list(PARENT_DECLARATIONS))
def test_resigned_parent_scientific_declarations_rejected(key):
    receipt = read(PARENT / "receipt.json")
    value = receipt[key]
    receipt[key] = (
        not value
        if type(value) is bool
        else value + 1
        if type(value) is int
        else "forged"
    )
    receipt = signed({k: v for k, v in receipt.items() if k != "identity"})
    assert receipt["identity"] == digest(
        {k: v for k, v in receipt.items() if k != "identity"}
    )
    with pytest.raises(CompletionError):
        parent_semantics(
            receipt, read(PARENT / "manifest.json"), read(PARENT / "decision.json")
        )


@pytest.mark.parametrize(
    "kind", ["schema", "extra_semantic_claim", "decision", "nested_reconstruction"]
)
def test_resigned_manifest_and_decision_semantics(kind, tmp_path):
    r, m, d = (
        read(PARENT / n) for n in ("receipt.json", "manifest.json", "decision.json")
    )
    if kind == "schema":
        m["schema"] = "forged"
    elif kind == "extra_semantic_claim":
        m["scientific_validation"] = "passed"
    elif kind == "decision":
        d["training_feasible"] = False
    else:
        r["numerical_reconstruction"]["new_QM"] = 4
    r, m, d = (
        signed({k: v for k, v in a.items() if k != "identity"}) for a in (r, m, d)
    )
    write(tmp_path / "decision.json", d)
    m["artifacts"]["decision.json"] = pin(tmp_path / "decision.json")
    m = signed({k: v for k, v in m.items() if k != "identity"})
    write(tmp_path / "manifest.json", m)
    r["manifest_sha256"] = pin(tmp_path / "manifest.json")
    r = signed({k: v for k, v in r.items() if k != "identity"})
    assert r["manifest_sha256"] == pin(tmp_path / "manifest.json")
    assert m["artifacts"]["decision.json"] == pin(tmp_path / "decision.json")
    with pytest.raises(CompletionError):
        parent_semantics(r, m, d)


def test_parent_schema_and_semantics_pass():
    assert parent_semantics(
        *(read(PARENT / n) for n in ("receipt.json", "manifest.json", "decision.json"))
    )


def test_exact_frozen_geometry_intake(model, plan):
    expected = read(PARENT / "geometry-only-next-design.json")["selected"]
    old = read(REPO / "docs/evidence/phase_4j65/fitting-dataset.json")["observations"]
    assert tuple(r["id"] for r in expected) == SELECTED
    assert all(r["geometry_identity"] == digest(r["coordinates"]) for r in expected)
    assert exact_intake(model, expected, expected, old, plan["validation_coordinates"])
    changed = copy.deepcopy(expected)
    changed[0]["coordinates"][0][0] += 1e-12
    with pytest.raises(CompletionError):
        exact_intake(model, changed, expected, old, plan["validation_coordinates"])


@pytest.mark.parametrize("change", ["duplicate", "rigid", "equivalent_H", "whole_end"])
def test_graph_symmetry_duplicate_rejection(model, plan, change):
    rows = copy.deepcopy(read(PARENT / "geometry-only-next-design.json")["selected"])
    x = np.array(rows[0]["coordinates"])
    if change == "rigid":
        x = x @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + 3
    elif change == "equivalent_H":
        x[[9, 10]] = x[[10, 9]]
    elif change == "whole_end":
        permutation = [7, 6, 8, 9, 5, 2, 1, 3, 4, *range(19, 28), *range(10, 19)]
        x = x[np.array(permutation) - 1]

    rows[1]["coordinates"] = x.tolist()
    with pytest.raises(CompletionError):
        exact_intake(model, rows, rows, [], plan["validation_coordinates"])


def test_validation_leakage_under_rigid_transform(model, plan):
    reserved = np.array(plan["validation_coordinates"][0])
    with pytest.raises(CompletionError):
        leakage(model, [reserved + 1], plan["validation_coordinates"])


def test_new_blocks_genuinely_omitted_and_primary_steps_together(plan):
    old = read(REPO / "docs/evidence/phase_4j65/fitting-dataset.json")["observations"]
    rows = old + plan["training_requests"]
    draws = bootstrap_draws(rows)
    assert draws == read(ROOT / "plan/expanded-bootstrap-draws.json")["draws"]
    assert len(draws) == 64 and all(len(r["sampled_blocks"]) == 14 for r in draws)
    for block in SELECTED:
        assert any(block not in r["sampled_blocks"] for r in draws)
        assert any(block in r["sampled_blocks"] for r in draws)
    for draw in draws:
        indices = [
            i
            for b in draw["sampled_blocks"]
            for i, r in enumerate(rows)
            if r["block"] == b
        ]
        for b in ("CH", "SiC", "methyl", "siloxane"):
            assert len([i for i in indices if rows[i]["block"] == b]) == 4 * draw[
                "sampled_blocks"
            ].count(b)


def test_protected_subset_prevents_error_dilution():
    data = synthetic(n=26)
    data["energy"][:] = 0
    data["force"][:] = 0
    data["QM_energy"][:22] = 0.105
    p = physical_problem(data)
    old = physical_problem(data, np.arange(22))
    lo = np.array([-1.0])
    hi = -lo
    checks = verified_subsets(p, old, [0, 0], lo, hi)
    assert checks["full"]["original_limits_satisfied"]
    assert not checks["protected22"]["original_limits_satisfied"]
    assert not physically_feasible(checks)


def test_protected_conic_scaling_and_single_intercept():
    pytest.importorskip("clarabel")
    data = synthetic(n=26)
    data["QM_energy"] = 1.2 * data["energy"][:, 0] + 7
    data["QM_force"] = 1.2 * data["force"][..., 0]
    p = physical_problem(data)
    old = physical_problem(data, np.arange(22))
    lo = np.array([-2.0])
    hi = -lo
    s = protected_system(p, lo, hi, protected=old)
    result = solve(s)
    assert result["values"] == pytest.approx([1.2, 7], abs=1e-6)
    assert physically_feasible(verified_subsets(p, old, result["values"], lo, hi))
    assert np.all(p["f"][:, -1] == 0) and np.all(p["e"][:, -1] == 1)
    assert s["A"].shape[1] == 3


def test_norms_units_and_energy_offset():
    data = synthetic(n=1, atoms=27)
    data["energy"][:] = 0
    data["force"][:] = 0
    data["QM_force"][0, 0, :] = [1.2, 1.2, 1.2]
    p = physical_problem(data)
    lo = np.array([-2.0])
    hi = -lo
    a = verify(p, [0, 0], lo, hi)
    b = verify(p, [0, 0.2], lo, hi)
    assert a["maximum_local_force_error"] == pytest.approx(np.sqrt(3) * 1.2)
    assert a["force_RMS"] == pytest.approx(np.sqrt(3 * 1.2**2 / 81))
    assert not a["original_limits_satisfied"] and a["force_RMS"] == b["force_RMS"]
    assert b["signed_energy_errors"] == [0.2]


@pytest.mark.parametrize(
    "status", ["AlmostSolved", "MaxIterations", "PrimalInfeasible", "NumericalError"]
)
def test_inaccurate_status_not_qualified(status):
    assert not solver_qualified(
        {
            "status": status,
            "primal_residual": 0.0,
            "dual_residual": 0.0,
            "absolute_gap": 0.0,
        }
    )


@pytest.mark.parametrize("field", ["primal_residual", "dual_residual", "absolute_gap"])
def test_solved_string_insufficient(field):
    r = {
        "status": "Solved",
        "primal_residual": 0.0,
        "dual_residual": 0.0,
        "absolute_gap": 0.0,
    }
    r[field] = 1e-4
    assert not solver_qualified(r)


@pytest.mark.parametrize(
    "gate", ["training", "numerical", "bounds", "identifiability", "uncertainty"]
)
def test_validation_nonactivation(gate):
    gates = dict.fromkeys(
        ["training", "numerical", "bounds", "identifiability", "uncertainty"], True
    )
    gates[gate] = False
    assert not validation_activation(gates)


def test_validation_consumption_is_one_time():
    with pytest.raises(CompletionError):
        validation_activation({"all": True}, True)


def test_bound_activity_not_removed_by_status():
    trial = {
        "checks": {
            "full": {"original_limits_satisfied": True, "active_bound_indices": [3]}
        },
        "solution": {
            "status": "Solved",
            "primal_residual": 0.0,
            "dual_residual": 0.0,
            "absolute_gap": 0.0,
        },
    }
    gates = candidate_gates(
        trial, {"qualified": True}, {"rank": 59, "scaled_condition": 100}
    )
    assert not gates["bounds"] and not validation_activation(gates)


@pytest.mark.parametrize(
    "key",
    [
        "graph_identity",
        "atom_order",
        "elements",
        "basis",
        "grid",
        "thresholds",
        "settings",
    ],
)
def test_request_graph_order_protocol_mutations(plan, key):
    row = plan["training_requests"][0]
    original = request(plan, row["coordinates"])
    mutated = copy.deepcopy(original)
    if isinstance(mutated[key], list):
        mutated[key] = mutated[key][::-1]
    elif isinstance(mutated[key], dict):
        mutated[key]["ints_tolerance"] = 1e-6
    else:
        mutated[key] = "forged"
    with pytest.raises(CompletionError):
        exact_request(mutated, request(plan, row["coordinates"]))


def test_parameter_mapping_bounds_and_units_frozen(plan):
    old = read(REPO / "docs/evidence/phase_4j65/plan/scientific-plan.json")
    assert plan["parameter_specifications"] == old["parameter_specifications"]
    assert len(plan["parameter_specifications"]) == 58


def test_native_dataset_original_producers_preserved():
    data = read(ROOT / "fitting-dataset.json")
    original = read(REPO / "docs/evidence/phase_4j65/fitting-dataset.json")
    assert data["observations"][:22] == original["observations"]
    assert all(
        r["producer"]["kind"] == "new_J67_execution" for r in data["observations"][22:]
    )
    assert data["energy_offset_count"] == 1


def test_all_bootstrap_failures_retained_without_resampling():
    summary = read(ROOT / "block-uncertainty.json")
    draws = read(ROOT / "plan/expanded-bootstrap-draws.json")["draws"]
    assert len(summary["attempts"]) == 64
    assert summary["replacement_resamples"] == 0
    for old, new in zip(draws, summary["attempts"], strict=True):
        assert old["sampled_blocks"] == new["sampled_blocks"]
        assert new["eligible"] == (not new["failure_categories"])
        assert set(new["indices"]) == {
            i
            for b in old["sampled_blocks"]
            for i, r in enumerate(read(ROOT / "fitting-dataset.json")["observations"])
            if r["block"] == b
        }


def test_complete_mapping_and_offline_relocated_fresh_process(tmp_path):
    env = {
        **os.environ,
        "PYTHONPATH": str(REPO) + ":" + str(REPO / "src"),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
    }
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/replay_phase_4j67.py"),
            "--repo",
            str(REPO),
            "--evidence",
            str(ROOT),
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_physical_columns_reconstruct_and_force_shapes(model):
    rows = read(ROOT / "fitting-dataset.json")["observations"]
    actual = features(model, rows)
    with np.load(ROOT / "physical-arrays.npz", allow_pickle=False) as stored:
        assert stored["names"].tolist() == actual["names"]
        for key in actual:
            if key != "names":
                assert np.array_equal(actual[key], stored[key])


@pytest.mark.parametrize(
    "mutation", ["graph", "order", "units", "protocol", "producer"]
)
def test_semantic_intake_rejects_resigned_observation(mutation, plan):
    expected = request(plan, plan["training_requests"][0]["coordinates"])
    changed = copy.deepcopy(expected)
    if mutation == "graph":
        changed["graph_identity"] = "0" * 64
    elif mutation == "order":
        changed["atom_order"][0:2] = [2, 1]
    elif mutation == "units":
        changed["units"] = {"gradient": "kcal/mol/angstrom"}
    elif mutation == "protocol":
        changed["settings"]["dft_radial_points"] = 150
    else:
        original = [
            {
                "id": "center",
                "producer": {"kind": "reused_original_execution", "sha": "original"},
            }
        ]
        replaced = copy.deepcopy(original)
        replaced[0]["producer"]["kind"] = "new_J67_execution"
        with pytest.raises(CompletionError):
            preserved_producers(replaced, original)
        return
    assert signed(changed)["identity"] == digest(changed)
    with pytest.raises(CompletionError):
        exact_request(changed, expected)


@pytest.mark.parametrize(
    "key",
    [
        "schema",
        "candidate_qualified",
        "parameter_completeness",
        "scientific_validation",
        "chemical_transfer",
        "automatic_assignment",
        "new_QM",
        "MD_steps",
        "training_feasible",
        "simulation_readiness",
        "original_author_artifact_status",
    ],
)
def test_resigned_j67_receipt_scientific_mutations(key):
    from scripts.validate_phase_4j67_receipt import semantics

    receipt, manifest = (read(ROOT / n) for n in ("receipt.json", "manifest.json"))
    expected_feasibility = receipt["training_feasible"]
    semantics(receipt, manifest, training_feasible=expected_feasibility)
    old = receipt[key]
    receipt[key] = (
        not old if type(old) is bool else old + 1 if type(old) is int else "forged"
    )
    changed = signed({k: v for k, v in receipt.items() if k != "identity"})
    with pytest.raises(CompletionError):
        semantics(changed, manifest, training_feasible=expected_feasibility)


@pytest.mark.parametrize("phase", ["parent", "J67"])
def test_resigned_receipt_unknown_scientific_claim_rejected(phase):
    from scripts.validate_phase_4j67_receipt import semantics

    root = PARENT if phase == "parent" else ROOT
    receipt, manifest = (read(root / n) for n in ("receipt.json", "manifest.json"))
    if phase == "parent":
        decision = read(PARENT / "decision.json")
        assert parent_semantics(receipt, manifest, decision)
    else:
        semantics(receipt, manifest, training_feasible=receipt["training_feasible"])
    receipt["scientific_validation_override"] = "passed"
    receipt = signed({k: v for k, v in receipt.items() if k != "identity"})
    with pytest.raises(CompletionError):
        if phase == "parent":
            parent_semantics(receipt, manifest, decision)
        else:
            semantics(receipt, manifest, training_feasible=False)


def test_conic_protected_force_limit_cannot_be_diluted():
    pytest.importorskip("clarabel")
    data = synthetic(n=26)
    # Keep the energy column identifiable so this test isolates protected force
    # constraints. A zero column invokes the historical 1e-30 scaling floor,
    # whose conservative arithmetic allowance makes its certificate uninformative.
    data["energy"][:, 0] = np.linspace(-0.1, 0.1, 26)
    data["force"][:] = 0
    data["QM_force"][:22] = 0.51
    full = physical_problem(data)
    old = physical_problem(data, np.arange(22))
    lo = np.array([-1.0])
    hi = -lo
    unconstrained_subset = solve(protected_system(full, lo, hi))
    guarded = solve(protected_system(full, lo, hi, protected=old))
    assert unconstrained_subset["rho"] < 1
    assert guarded["rho"] == pytest.approx(1.02, abs=1e-6)
    from island.pcff_qm_pilot.hmds_constrained import dual_lower_bound

    certificate = dual_lower_bound(
        protected_system(full, lo, hi, protected=old), guarded, old, lo, hi
    )
    assert certificate["lower_bound"] > 1.0199

    assert not physically_feasible(
        verified_subsets(full, old, guarded["values"], lo, hi)
    )


def test_native_validation_nonactivation_reaches_no_worker(tmp_path, monkeypatch):
    import scripts.phase_4j67 as phase

    root = tmp_path / "evidence"
    work = tmp_path / "work"
    write(
        root / "validation-activation.json",
        signed({"activated": False, "gates": {"uncertainty": False}}),
    )
    write(work / "attempts.json", {"calls": {}})

    def forbidden(*args, **kwargs):
        raise AssertionError("No validation worker may launch")

    monkeypatch.setattr(phase, "run_job", forbidden)
    assert phase.conditional_validation(REPO, root, work, Path(sys.executable)) == 0
    assert not (root / "validation-consumption.json").exists()


def test_reconstructed_weak_sensitivities_mutation_rejected(model):
    from island.pcff_qm_pilot.hmds_augmentation import proposal_sensitivities

    proposed = read(PARENT / "geometry-only-next-design.json")
    diagnosis = read(PARENT / "block-identifiability-diagnosis.json")
    rows = proposal_sensitivities(model, proposed, diagnosis)
    assert (
        len(rows) == 4
        and max(r["maximum_absolute_sensitivity_error"] for r in rows) < 1e-12
    )
    changed = copy.deepcopy(proposed)
    changed["selected"][0]["target_weak_direction_norms"][0] += 0.001
    changed = signed({k: v for k, v in changed.items() if k != "identity"})
    with pytest.raises(CompletionError):
        proposal_sensitivities(model, changed, diagnosis)


def test_raw_solver_state_independent_closure_and_mapping():
    pytest.importorskip("clarabel")
    from island.pcff_qm_pilot.hmds_augmentation import solver_state_diagnostics

    p = physical_problem(synthetic())
    lo = np.array([-2.0])
    hi = -lo
    system = protected_system(p, lo, hi)
    result = solve(system)
    actual = solver_state_diagnostics(system, result)
    assert actual["primal_normalized_L2"] < 1e-8
    assert actual["dual_normalized_L2"] < 1e-8
    changed = copy.deepcopy(result)
    changed["values"][0] += 0.1
    with pytest.raises(CompletionError):
        solver_state_diagnostics(system, changed)
    changed = copy.deepcopy(result)
    changed["slack"][0] += 0.1
    assert (
        solver_state_diagnostics(system, changed)["primal_maximum_absolute_closure"]
        > 0.09
    )


def test_tiny_inward_selection_cannot_erase_raw_bound_activity():
    trial = {
        "checks": {
            "full": {"original_limits_satisfied": True, "active_bound_indices": []}
        },
        "raw_active_original_bounds": [3],
        "solution": {
            "status": "Solved",
            "primal_residual": 0.0,
            "dual_residual": 0.0,
            "absolute_gap": 0.0,
        },
    }
    gates = candidate_gates(
        trial, {"qualified": True}, {"rank": 59, "scaled_condition": 100}
    )
    assert not gates["bounds"] and not validation_activation(gates)


def test_interrupted_solve_artifacts_cannot_be_overwritten(tmp_path):
    from scripts.analyze_phase_4j67 import attempt

    raw = tmp_path / "failed-solver.json"
    raw.write_bytes(b"preserved interrupted solver state\n")
    with pytest.raises(CompletionError):
        attempt(
            physical_problem(synthetic()),
            np.array([-1.0]),
            np.array([1.0]),
            tmp_path,
            "failed",
            {},
        )
    assert raw.read_bytes() == b"preserved interrupted solver state\n"
