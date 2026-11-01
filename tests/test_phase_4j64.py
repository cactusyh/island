"""Numerical controls for the two-variable HMDS research extension, never QM mocks."""

import copy
from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError
from island.pcff_qm_pilot.hmds_baseline import (
    NAMES,
    calibrate,
    corrections,
    evaluate,
    extended_parts,
    objective_status,
)
from island.pcff_qm_pilot.hmds_hamiltonian import load_model, parts
from island.pcff_qm_pilot.hmds_independent import scalar
from island.pcff_qm_pilot.hmds_paths import path
from island.pcff_qm_pilot.hmds_precision import (
    check_budget,
    path_response,
    protocol_contrast,
    read,
    rounding_audit,
)
from scripts.analyze_phase_4j63 import gate
from scripts.analyze_phase_4j64 import extended_matrix
from scripts.curvature_phase_4j63 import unknown_null_check

REPO = Path(__file__).resolve().parents[1]
OLD = REPO / "docs/evidence/phase_4j63"


@pytest.fixture(scope="module")
def model():
    return load_model(REPO)


@pytest.fixture(scope="module")
def center():
    return np.array(read(OLD / "stationary-center.json")["center"]["coordinates"])


def independent_correction(model, x, values):
    total = 0.0
    for row in model["interaction_inventory"]:
        if row["family"] != "quartic_bond":
            continue
        i, j = row["sites"]
        types = {model["typing"][str(i)], model["typing"][str(j)]}
        key = (
            NAMES[0]
            if types == {"sio", "c3"}
            else NAMES[1]
            if types == {"c3", "hc"}
            else None
        )
        if key:
            total += (
                values[key]
                * (np.linalg.norm(x[i - 1] - x[j - 1]) - row["assignment"]["values"][0])
                ** 2
            )
    return total


@pytest.mark.parametrize("step", [2e-5, 1e-5])
def test_nonzero_correction_scalar_force_and_parameter_derivatives(model, center, step):
    x = center + np.random.default_rng(4064).normal(0, 0.01, center.shape)
    values = {NAMES[0]: 13.7, NAMES[1]: -4.2}  # manufactured, not physical
    c = corrections(model, x)
    vector = np.array(list(values.values()))
    assert c["energy_jacobian"] @ vector == pytest.approx(
        independent_correction(model, x, values), abs=1e-12
    )
    f = c["force_jacobian"] @ vector
    fd = np.zeros_like(x)
    for i, j in np.ndindex(x.shape):
        xp, xm = x.copy(), x.copy()
        xp[i, j] += step
        xm[i, j] -= step
        fd[i, j] = -(
            independent_correction(model, xp, values)
            - independent_correction(model, xm, values)
        ) / (2 * step)
    assert np.max(abs(f - fd)) < 2e-7
    for k, name in enumerate(NAMES):
        vp, vm = dict(values), dict(values)
        vp[name] += step
        vm[name] -= step
        de = (
            independent_correction(model, x, vp) - independent_correction(model, x, vm)
        ) / (2 * step)
        assert de == pytest.approx(c["energy_jacobian"][k], abs=1e-9)


def test_complete_manufactured_energy_and_force_require_all58(model, center):
    p = extended_parts(model, center)
    values = {n: float(0.3 + np.sin(i)) for i, n in enumerate(p["columns"])}
    result = evaluate(model, center, values)
    independent = scalar(model, center)
    expected = independent["fixed_energy"] + sum(
        independent["columns"][n] * values[n] for n in p["columns"][:-2]
    )
    expected += independent_correction(model, center, values)
    assert result["energy"] == pytest.approx(expected, abs=2e-10)
    for missing in [
        None,
        {},
        {n: values[n] for n in p["columns"][:-2]},
        {n: values[n] for n in NAMES},
        {**values, "other": 1},
    ]:
        with pytest.raises(CompletionError, match="explicit_all58"):
            evaluate(model, center, missing)
    with pytest.raises(CompletionError, match="nonfinite"):
        evaluate(model, center, {**values, NAMES[0]: float("nan")})
    assert np.max(abs(result["forces"].sum(axis=0))) < 1e-8


@pytest.mark.parametrize("name,index", [("SiC", 0), ("CH", 1)])
@pytest.mark.parametrize("h", [0.006, 0.003])
def test_unknown_null_and_correction_sensitivity(model, center, name, index, h):
    spec = next(
        p for p in read(OLD / "plan/scientific-plan.json")["paths"] if p["name"] == name
    )
    null = unknown_null_check(model, center, spec, h)
    assert len(null["unknown_columns"]) == 56
    points = [path(model, center, spec, s * h) for s in (0, 1, -1)]
    c = [corrections(model, p["coordinates"]) for p in points]
    scalar_sensitivity = (
        c[1]["energy_jacobian"] - 2 * c[0]["energy_jacobian"] + c[2]["energy_jacobian"]
    ) / h**2
    derivative = [
        -np.einsum("iap,ia->p", a["force_jacobian"], p["tangent"])
        for a, p in zip(c, points, strict=True)
    ]
    force_sensitivity = (derivative[1] - derivative[2]) / (2 * h)
    expected = np.zeros(2)
    expected[index] = 2
    assert np.max(abs(scalar_sensitivity - expected)) < 1e-8
    assert np.max(abs(force_sensitivity - expected)) < 1e-8
    # This uses every fixed source and pair contribution, not the bare bond K2.
    source = next(
        r
        for r in model["interaction_inventory"]
        if r["family"] == "quartic_bond" and set(r["sites"]) == set(spec["bond"])
    )
    assert (
        abs(
            null["verified_fixed_response"]["curvature_from_path_gradient"]
            - 2 * source["assignment"]["values"][1] / 4.184
        )
        > 0.01
    )


def test_quadratic_changes_no_source_equilibria_or_cross_kernels(model, center):
    snapshot = copy.deepcopy(model)
    original = parts(model, center)
    extended = extended_parts(model, center)
    assert np.array_equal(original["energy_jacobian"], extended["energy_jacobian"][:-2])
    assert np.array_equal(
        original["force_jacobian"], extended["force_jacobian"][:, :, :-2]
    )
    assert model == snapshot


@pytest.mark.parametrize(
    "mutation",
    ["coordinates", "atom_order", "charge", "multiplicity", "method", "settings"],
)
def test_grid_only_contrast_rejects_confounds(mutation):
    r = read(REPO / "docs/evidence/phase_4j64/plan/request-inventory.json")[
        "dense/center"
    ]
    other = copy.deepcopy(r)
    if mutation == "coordinates":
        other[mutation][0][0] += 0.01
    elif mutation == "atom_order":
        other[mutation][0:2] = other[mutation][1::-1]
    elif mutation in ("charge", "multiplicity"):
        other[mutation] += 1
    elif mutation == "method":
        other[mutation] = "b3lyp"
    else:
        other[mutation]["ints_tolerance"] = 1e-8
    with pytest.raises(CompletionError, match="confounded"):
        protocol_contrast(r, other, ["grid"])


def test_rounding_bound_does_not_hide_soft_response_discrepancy():
    center = Decimal("-893.3818442288735469")
    delta = Decimal("0.00000045")
    values = [
        (center, float(center)),
        (center + delta, float(center + delta)),
        (center + delta, float(center + delta)),
    ]
    r = rounding_audit(values, 0.02)
    assert (
        abs(
            r["decimal_energy_curvature"]
            - r["binary_subtract_before_conversion_curvature"]
        )
        < r["binary_serialization_arithmetic_bound"]
    )
    assert r["printed_decimal_rounding_bound"] < 1e-8
    assert r["binary_serialization_arithmetic_bound"] < 1e-6
    assert r["binary_serialization_arithmetic_bound"] < 0.105 / 10000


def test_curved_path_chain_rule_includes_acceleration(model, center):
    spec = next(
        s
        for s in read(OLD / "plan/scientific-plan.json")["paths"]
        if s["name"] == "siloxane"
    )
    h = 0.0001
    geos = [path(model, center, spec, v) for v in (0, h, -h)]
    # Manufactured linear Cartesian energy has zero H, nonzero g.a on a rotation.
    g = np.random.default_rng(64).normal(size=center.shape)
    samples = [
        {"energy": float(np.sum(g * p["coordinates"])), "forces": -g} for p in geos
    ]
    energies = [(Decimal(str(s["energy"])), s["energy"]) for s in samples]
    r = path_response(samples, geos, energies, h)
    expected = float(np.sum(g * geos[0]["acceleration"]))
    assert r["t_H_t_finite_step"] == 0
    assert r["curvature_from_path_gradient"] == pytest.approx(expected, abs=1e-7)
    assert r["chain_rule_sum"] == pytest.approx(expected)
    assert abs(expected) > 0.1
    samples[1]["forces"] = [[0, 0, 0]] * 26
    with pytest.raises(CompletionError, match="full_HMDS_force"):
        path_response(samples, geos, energies, h)


def test_manufactured_curvature_recovery_and_rank_failure():
    fixed = np.array([375.0, 721.0])
    v = np.array([25.3, -8.7])
    assert np.allclose(calibrate(fixed + 2 * v, fixed, 2 * np.eye(2)), v)
    for sensitivity in (np.ones((2, 2)), np.zeros((2, 2))):
        with pytest.raises(CompletionError, match="not_identifiable"):
            calibrate([1, 2], [0, 0], sensitivity)
    with pytest.raises(CompletionError, match="finite"):
        calibrate([np.nan, 1], fixed, 2 * np.eye(2))


def test_soft_torsion_failure_does_not_block_independent_curvature_calibration():
    status = objective_status(True, True, False, True)
    assert status["baseline_calibration_completed"]
    assert not status["soft_torsion_reference_precision"]
    assert not status["complete_force_fit_authorized"]
    assert status["candidate_validation"] == "not_established"
    assert status["parameter_completeness"] == "incomplete"
    assert not objective_status(True, False, True, True)[
        "baseline_calibration_completed"
    ]
    assert not gate(0.105 + 1.42, 1.42, 0.02, 0.02)["passed"]


@pytest.mark.parametrize(
    "count,elapsed,passes", [(23, 1, True), (24, 1, False), (1, 36000, False)]
)
def test_cap_counts_failed_and_interrupted_attempts(count, elapsed, passes):
    plan = {"limits": {"maximum_new_gradients": 24, "campaign_wall_seconds": 36000}}
    attempts = {str(i): {"status": "failed"} for i in range(count)}
    if passes:
        check_budget(plan, attempts, elapsed)
    else:
        with pytest.raises(CompletionError, match="budget"):
            check_budget(plan, attempts, elapsed)


def test_complete58_manufactured_joint_recovery(model, center):
    old = read(OLD / "joint-design.json")
    reference = read(OLD / "normalized-reference.json")["samples"]
    geometries = [
        center,
        *[
            s["coordinates"]
            for k, s in sorted(reference.items())
            if k.startswith("paths/")
        ],
        *[s["coordinates"] for s in old["selected_additional_configurations"]],
    ]
    a, names = extended_matrix(model, geometries, center)
    assert len(names) == 58
    scales = np.linalg.norm(a, axis=0)
    actual = np.sin(np.arange(58) + 0.3)
    observed = a @ actual
    recovered, _, rank, _ = np.linalg.lstsq(a / scales, observed, rcond=1e-10)
    assert rank == 58
    assert np.max(abs(recovered / scales - actual)) < 1e-7


def test_failed_launch_retains_diagnostics_and_consumes_budget(tmp_path):
    from island.pcff_completion import digest, write
    from scripts.phase_4j64 import run_job

    root, work = tmp_path / "evidence", tmp_path / "work"
    request = {"manufactured": "no QM executable invoked"}
    plan = {
        "execution_sources": {},
        "limits": {"maximum_new_gradients": 24, "campaign_wall_seconds": 36000},
    }
    for base in (root / "plan", work):
        write(base / "scientific-plan.json", plan)
    write(
        work / "jobs.json",
        {
            "jobs": {
                "manufactured/startup": {
                    "request_identity": digest(request),
                    "wall_seconds": 2,
                }
            }
        },
    )
    write(work / "requests/manufactured/startup.json", request)
    assert (
        run_job(
            REPO, root, work, tmp_path / "missing-interpreter", "manufactured/startup"
        )
        == 127
    )
    result = read(work / "processes/manufactured/startup/result.json")
    assert result["status"] == "unexpected_failure"
    assert (
        read(work / "attempts.json")["calls"]["manufactured/startup"]["returncode"]
        == 127
    )
    assert (
        "FileNotFoundError"
        in (work / "processes/manufactured/startup/stderr.log").read_text()
    )
    assert not (work / "outputs").exists()
    write(root / "plan/startup-maintenance.json", {"retry_job": "different"})
    with pytest.raises(CompletionError, match="additive_retry"):
        run_job(
            REPO, root, work, tmp_path / "missing-interpreter", "manufactured/startup"
        )


@pytest.mark.parametrize("name", ["SiC", "CH", "methyl", "siloxane"])
def test_independent_endpoint_tangents_and_accelerations(model, center, name):
    from scripts.independent_phase_4j64 import independent_path

    spec = next(
        p for p in read(OLD / "plan/scientific-plan.json")["paths"] if p["name"] == name
    )
    s, h = 0.015, 1e-4
    geo = path(model, center, spec, s)
    t, a = independent_path(model, center, geo["coordinates"], spec)
    pp, pm = [path(model, center, spec, s + sign * h) for sign in (1, -1)]
    assert np.max(abs(t - (pp["coordinates"] - pm["coordinates"]) / (2 * h))) < 2e-8
    assert np.max(abs(a - (pp["tangent"] - pm["tangent"]) / (2 * h))) < 2e-8


def test_quadratic_columns_rigid_invariance_and_torsion_independence(model, center):
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    c = corrections(model, center)
    rotated = corrections(model, center @ rotation.T + [2, 4, 6])
    assert np.allclose(rotated["energy_jacobian"], c["energy_jacobian"], atol=1e-14)
    assert np.allclose(
        rotated["force_jacobian"],
        np.einsum("ij,njp->nip", rotation, c["force_jacobian"]),
        atol=1e-14,
    )
    for spec in read(OLD / "plan/scientific-plan.json")["paths"]:
        if spec["kind"] == "rotation":
            moved = path(model, center, spec, 0.43)["coordinates"]
            assert (
                np.max(
                    abs(
                        corrections(model, moved)["energy_jacobian"]
                        - c["energy_jacobian"]
                    )
                )
                < 1e-14
            )


def test_baseline_rejects_graph_or_source_mutation(model, center):
    wrong = copy.deepcopy(model)
    wrong["graph_identity"] = "unsupported"
    with pytest.raises(CompletionError, match="mutation|unsupported"):
        corrections(wrong, center)
    wrong = copy.deepcopy(model)
    row = next(r for r in wrong["interaction_inventory"] if r["family"] == "quartic_bond")
    row["assignment"]["values"][0] += 0.01
    with pytest.raises(CompletionError, match="mutation|unsupported"):
        corrections(wrong, center)


def test_exact_center_correction_keeps_bitwise_protocol_separation():
    from island.pcff_qm_pilot.hmds_engine import request

    root = REPO / "docs/evidence/phase_4j64"
    plan = read(root / "plan/scientific-plan.json")
    amendment = read(root / "plan/exact-center-correction.json")
    exact = amendment["request"]
    retained = request(plan, plan["center_coordinates"], thresholds=[1e-12, 1e-10])
    assert exact["coordinates"] == plan["center_coordinates"]
    assert protocol_contrast(retained, exact, ["grid"])["changed_fields"] == ["grid"]
    original = read(root / "plan/request-inventory.json")["dense/center"]
    assert original["coordinates"] != exact["coordinates"]
    assert (
        np.max(abs(np.array(original["coordinates"]) - exact["coordinates"]))
        == amendment["maximum_original_difference_A"]
    )
    assert amendment["maximum_process_attempts_including_startup_failure"] <= 24


@pytest.fixture(scope="module")
def retained_native(tmp_path_factory):
    from island.pcff_qm_pilot.bundle import unpack

    native = tmp_path_factory.mktemp("j64-retained") / "native"
    unpack(REPO / "docs/evidence/phase_4j64/native-outputs.tar.xz", native)
    return native


def test_real_dense_native_reconstruction_before_energy_subtraction(retained_native):
    from island.pcff_qm.reference import CONVERSION
    from island.pcff_qm_pilot.hmds_precision import decimal_energy
    from scripts.analyze_phase_4j63 import independent_native, sample
    from scripts.independent_phase_4j64 import ascii_gradient, native_decimal_total

    root = REPO / "docs/evidence/phase_4j64"
    for key in ("exact/center", "dense/plus", "dense/minus", "qz/center"):
        actual = sample(retained_native, key)
        independent_native(retained_native, key, actual)
        gradient = ascii_gradient(retained_native / "outputs" / key)
        assert np.max(
            abs(
                -gradient * CONVERSION["energy_factor"] / CONVERSION["distance_factor"]
                - actual["forces"]
            )
        ) < 1e-10
        decimal, binary = decimal_energy(retained_native / "outputs" / key)
        assert decimal == native_decimal_total(retained_native / "outputs" / key)
        assert abs(float(decimal) - binary) < 2e-10
        assert actual == read(root / "normalized-reference.json")["samples"][key]
    checks = read(root / "independent-path-and-baseline-checks.json")["paths"]
    assert set(checks) == {"SiC", "CH", "siloxane"}
    assert all(
        abs(r["energy_curvature_difference_from_production"]) < 1e-9
        for r in checks.values()
    )


@pytest.mark.parametrize(
    "mutation,reason",
    [
        ("units", "unit_or_axis"),
        ("axis_order", "unit_or_axis"),
        ("scf_converged", "convergence"),
        ("multiplicity", "state"),
        ("basis", "protocol"),
        ("geometry_bohr", "geometry"),
        ("energy_hartree", "energy"),
        ("gradient", "gradient"),
        ("unverified_native", "hash"),
    ],
)
def test_real_importer_rejects_mutations(retained_native, tmp_path, mutation, reason):
    import shutil

    from island.pcff_completion import pin, write
    from island.pcff_qm_pilot.molecular_native import normalized

    key = "exact/center"
    directory = tmp_path / "mutated"
    shutil.copytree(retained_native / "outputs" / key, directory)
    request = read(retained_native / "requests" / (key + ".json"))
    result = read(directory / "result.json")
    if mutation == "units":
        result["units"]["gradient"] = "kcal/mol/angstrom"
    elif mutation == "axis_order":
        result[mutation] = "Cartesian-major"
    elif mutation == "scf_converged":
        result[mutation] = False
    elif mutation == "multiplicity":
        result[mutation] = 3
    elif mutation == "basis":
        result[mutation] = "def2-svp"
    elif mutation == "geometry_bohr":
        result[mutation][0][0] += 0.001
    elif mutation == "energy_hartree":
        result[mutation] += 0.001
    elif mutation == "gradient":
        gradient = np.load(directory / "gradient.npy")
        np.save(directory / "gradient.npy", gradient[:-1], allow_pickle=False)
        # Rebind the file to exercise dimensional validation beyond file hashing.
        result["native_files"]["gradient.npy"] = pin(directory / "gradient.npy")
    else:
        with (directory / "psi4.out").open("a") as stream:
            stream.write("\nmutation\n")
    write(directory / "result.json", result)
    with pytest.raises(CompletionError, match=reason):
        normalized(directory, request)


def test_real_ascii_missing_component_is_rejected(retained_native, tmp_path):
    import shutil

    from scripts.independent_phase_4j64 import ascii_gradient

    shutil.copyfile(
        retained_native / "outputs/exact/center/gradient.dat", tmp_path / "gradient.dat"
    )
    target = tmp_path / "gradient.dat"
    target.write_text("\n".join(target.read_text().splitlines()[:-1]) + "\n")
    with pytest.raises(CompletionError, match="coverage"):
        ascii_gradient(tmp_path)


def test_fresh_process_relocation_reconstructs_real_precision_reports(tmp_path):
    import json
    import os
    import shutil
    import subprocess
    import sys

    relocated = tmp_path / "relocated-evidence"
    shutil.copytree(REPO / "docs/evidence/phase_4j64", relocated)
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/replay_phase_4j64.py"),
            "replay",
            "--repo",
            str(REPO),
            "--evidence",
            str(relocated),
        ],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHONPATH": str(REPO) + os.pathsep + str(REPO / "src"),
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
        },
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    value = json.loads(result.stdout)
    assert value["relocated"] and value["fresh_work_directory"]
    assert not value["QM_invoked"]
    assert value["native_gradients"] == read(relocated / "compute-accounting.json")[
        "successful_gradients"
    ]
