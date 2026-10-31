"""J61 numerical baseline propagation, joint derivatives and rigid projection."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError
from island.pcff_qm_pilot.hamiltonian import evaluate, load_model
from island.pcff_qm_pilot.stretch_baseline import (
    Configuration,
    bond_kind,
    directional_curvature,
    source_baseline,
    stretching_projection,
)
from island.pcff_qm_pilot.stretch_fit import Problem, aligned_rms, leakage
from scripts.independent_phase_4j60 import scalar

REPO = Path(__file__).resolve().parents[1]
OLD = REPO / "docs/evidence/phase_4j60"
BONDS = [[1, 2], [2, 3], [1, 4], [1, 5], [1, 6], [3, 7], [3, 8], [3, 9]]


def read(path):
    return json.loads(Path(path).read_text())


@pytest.fixture(scope="module")
def model():
    return load_model(REPO)


@pytest.fixture(scope="module")
def center():
    return np.asarray(read(OLD / "plan/reference-plan.json")["geometries"]["center"])


def independent(model, xyz, p):
    changed = copy.deepcopy(model)
    for row in changed["interaction_inventory"]:
        if row["family"] == "quartic_bond":
            k = bond_kind(model, row["sites"])
            row["assignment"]["values"][0] = p[k]
            row["assignment"]["values"][1] = p[k + 2] * 4.184
    e, columns = scalar(changed, xyz)
    return e + sum(v * columns[k] for k, v in zip(sorted(columns), p[4:], strict=True))


def parameters(model):
    p = np.r_[source_baseline(model), np.linspace(-2.4, 3.1, 20)]
    p[:2] += [0.011, -0.009]
    p[2:4] *= [1.08, 0.92]
    return p


def test_original_baseline_reconstructs_J60(model, center):
    configuration = Configuration(model, center)
    p = np.r_[source_baseline(model), np.linspace(1, 20, 20)]
    new = configuration.evaluate(p)
    old = evaluate(model, center, dict(zip(configuration.columns, p[4:], strict=True)))
    assert new["energy"] == pytest.approx(old["energy"], abs=1e-11)
    assert new["forces"] == pytest.approx(old["forces"], abs=1e-10)


@pytest.mark.parametrize("h", [1e-4, 5e-5])
def test_independent_energy_Cartesian_forces(model, center, h):
    p = parameters(model)
    r = Configuration(model, center).evaluate(p)
    assert r["energy"] == pytest.approx(independent(model, center, p), abs=1e-10)
    f = np.zeros((9, 3))
    for i in range(9):
        for j in range(3):
            a, b = center.copy(), center.copy()
            a[i, j] += h
            b[i, j] -= h
            f[i, j] = -(independent(model, a, p) - independent(model, b, p)) / (2 * h)
    assert r["forces"] == pytest.approx(f, abs=2e-5)


@pytest.mark.parametrize("index", range(24))
def test_all_energy_and_mixed_force_parameter_derivatives(model, center, index):
    c = Configuration(model, center)
    p = parameters(model)
    r = c.evaluate(p)
    h = 1e-6 if index < 2 else 1e-3
    a, b = p.copy(), p.copy()
    a[index] += h
    b[index] -= h
    ra, rb = c.evaluate(a), c.evaluate(b)
    assert (ra["energy"] - rb["energy"]) / (2 * h) == pytest.approx(
        r["energy_parameter_jacobian"][index], abs=2e-6
    )
    assert (ra["forces"] - rb["forces"]) / (2 * h) == pytest.approx(
        r["force_parameter_jacobian"][:, :, index], abs=2e-6
    )


@pytest.mark.parametrize(
    "family", ["bond-bond_1_3", "end_bond-torsion_3", "middle_bond-torsion_3"]
)
def test_r0_propagation_through_fitted_cross_terms(model, center, family):
    c = Configuration(model, center)
    p = np.r_[source_baseline(model), np.zeros(20)]  # synthetic isolation control
    column = next(k for k in c.columns if k.startswith(family + ":"))
    p[4 + c.columns.index(column)] = 2.7
    shifted = p.copy()
    shifted[:2] += [0.017, -0.013]
    for q in (p, shifted):
        assert c.evaluate(q)["energy"] == pytest.approx(
            independent(model, center, q), abs=1e-10
        )
    delta = c.evaluate(shifted)["energy"] - c.evaluate(p)["energy"]
    plain, plain_shift = p.copy(), shifted.copy()
    plain[4:] = plain_shift[4:] = 0
    coupling_delta = (
        delta - c.evaluate(plain_shift)["energy"] + c.evaluate(plain)["energy"]
    )
    assert abs(coupling_delta) > 1e-5


def test_projection_rigid_modes_and_rotation_covariance(center):
    p = stretching_projection(center, BONDS)
    q = p["Q"]
    assert q.T @ q == pytest.approx(np.eye(8), abs=1e-12)
    assert p["rigid"].T @ q == pytest.approx(np.zeros((6, 8)), abs=1e-12)
    rotation, _ = np.linalg.qr(np.random.default_rng(4061).normal(size=(3, 3)))
    new = stretching_projection(center @ rotation + [2.3, -4.1, 8.2], BONDS)
    rotate = np.kron(np.eye(9), rotation.T)
    assert new["Q"] @ new["Q"].T == pytest.approx(
        rotate @ q @ q.T @ rotate.T, abs=1e-11
    )


def test_projected_Hessian_sign_units_and_step_convergence(model, center):
    q = stretching_projection(center, BONDS)["Q"]
    p = parameters(model)
    hs = []
    for step in (0.004, 0.002):
        a = [
            Configuration(model, center + (step * q[:, i]).reshape(9, 3)).evaluate(p)[
                "forces"
            ]
            for i in range(8)
        ]
        b = [
            Configuration(model, center - (step * q[:, i]).reshape(9, 3)).evaluate(p)[
                "forces"
            ]
            for i in range(8)
        ]
        response, raw, hessian = directional_curvature(a, b, q, step)
        assert q.T @ response == pytest.approx(raw, abs=1e-12)
        assert hessian == pytest.approx(hessian.T, abs=1e-12)
        hs.append(hessian)
    assert np.linalg.norm(hs[0] - hs[1]) / np.linalg.norm(hs[1]) < 1e-3
    for i in (0, 2):
        d = (0.002 * q[:, i]).reshape(9, 3)
        diagonal = (
            independent(model, center + d, p)
            - 2 * independent(model, center, p)
            + independent(model, center - d, p)
        ) / 0.002**2
        assert diagonal == pytest.approx(hs[1][i, i], rel=2e-4)


@pytest.fixture(scope="module")
def manufactured(model, center):
    values = parameters(model)
    design = read(OLD / "joint-target-design.json")
    records = []
    for k, s in enumerate(design["selected_training"][:10]):
        r = Configuration(model, s["coordinates"]).evaluate(values)
        records.append(
            {
                "key": f"manufactured-{k}",
                "coordinates": s["coordinates"],
                "energy": r["energy"] + 29.3,
                "forces": r["forces"].tolist(),
                "graph_identity": model["graph_identity"],
                "atom_order": model["atom_order"],
            }
        )
    r = Configuration(model, center).evaluate(values)
    stationary = {
        "coordinates": center.tolist(),
        "energy": r["energy"] + 29.3,
        "forces": r["forces"].tolist(),
        "graph_identity": model["graph_identity"],
        "atom_order": model["atom_order"],
    }
    q = stretching_projection(center, BONDS)["Q"]
    problem = Problem(model, records, stationary, q, np.eye(8))
    h, _ = problem.mm_curvature(values)
    problem.curvature = h
    return problem, values


def test_manufactured_full24_nonzero_recovery(manufactured):
    problem, values = manufactured
    initial = values.copy()
    initial[:2] += [-0.007, 0.006]
    initial[2:4] *= [1.08, 0.92]
    initial[4:] += 0.5
    result = problem.fit(initial)
    recovered = np.array(list(result["parameters"].values()))
    assert recovered == pytest.approx(values, abs=1e-6)
    assert result["rank"] == 24
    assert result["weighted_residual_SS"] < 1e-14
    assert result["training_reference_offset"] == pytest.approx(29.3, abs=1e-8)


def test_curvature_parameter_jacobian(manufactured):
    problem, values = manufactured
    _, j = problem.mm_curvature(values)
    for i in (0, 1, 2, 3, 14, 23):
        step = 1e-5 if i < 2 else 1e-3
        a, b = values.copy(), values.copy()
        a[i] += step
        b[i] -= step
        ha, _ = problem.mm_curvature(a)
        hb, _ = problem.mm_curvature(b)
        assert (ha - hb) / (2 * step) == pytest.approx(j[:, :, i], abs=5e-5)


def test_bounds_reject_illegal_initial(manufactured):
    problem, values = manufactured
    p = values.copy()
    p[0] += 1
    with pytest.raises(CompletionError, match="initial_parameter_bounds"):
        problem.fit(p)


def test_nonidentifiable_system_rejects(manufactured, monkeypatch):
    problem, values = manufactured
    monkeypatch.setattr(
        problem, "system", lambda x: (np.zeros(100), np.zeros((100, 24)))
    )
    with pytest.raises(CompletionError, match="unidentifiable"):
        problem.fit(values)


def test_symmetry_and_rigid_validation_leakage(center):
    changed = center[[2, 1, 0, 6, 8, 7, 3, 5, 4]] + [2, 3, -5]
    assert aligned_rms(center, changed) < 1e-10
    with pytest.raises(CompletionError, match="validation_symmetry_leakage"):
        leakage(
            [{"coordinates": center.tolist()}],
            [{"coordinates": changed.tolist()}, {"coordinates": center.tolist()}],
        )


def test_mutating_cached_source_rejected(model, center):
    changed = copy.deepcopy(model)
    c = Configuration(changed, center)
    changed["interaction_inventory"][0]["id"] = "mutated"
    with pytest.raises(CompletionError, match="residual_model_mutation"):
        c.evaluate(parameters(model))


def test_repeated_geometry_with_curvature_is_ill_conditioned(model, center):
    """A formally full-rank curvature block does not rescue uninformative data."""
    values = parameters(model)
    result = Configuration(model, center).evaluate(values)
    sample = {
        "coordinates": center.tolist(),
        "energy": result["energy"],
        "forces": result["forces"].tolist(),
        "graph_identity": model["graph_identity"],
        "atom_order": model["atom_order"],
    }
    q = stretching_projection(center, BONDS)["Q"]
    problem = Problem(model, [sample] * 4, sample, q, np.eye(8))
    problem.curvature = problem.mm_curvature(values)[0]
    _, jacobian = problem.system(values)
    scaled = jacobian / np.linalg.norm(jacobian, axis=0)
    singular = np.linalg.svd(scaled, compute_uv=False)
    assert sum(singular > singular[0] * 1e-10) == 24
    assert singular[0] / singular[-1] > 1000
    with pytest.raises(CompletionError, match="baseline_joint_unidentifiable"):
        problem.fit(values)


@pytest.fixture(scope="module")
def native_bundle(tmp_path_factory):
    from island.pcff_qm_pilot.bundle import unpack

    root = tmp_path_factory.mktemp("J61-native") / "native"
    unpack(REPO / "docs/evidence/phase_4j61/native-outputs.tar.xz", root)
    return root


def test_real_stationary_native_units_and_curvature(native_bundle, tmp_path):
    import shutil

    from scripts.analyze_phase_4j61 import curvature_report
    from scripts.independent_phase_4j59 import native_values
    from scripts.phase_4j61 import imported

    evidence = REPO / "docs/evidence/phase_4j61"
    shutil.copytree(evidence / "plan", tmp_path / "plan")
    for name in ("stationary-center.json", "stretching-subspace.json"):
        shutil.copyfile(evidence / name, tmp_path / name)
    record = imported(native_bundle, "optimization/step-012")
    e, forces = native_values(native_bundle / "outputs/optimization/step-012")
    assert record["energy"] == pytest.approx(e, abs=2e-7)
    assert record["forces"] == pytest.approx(forces, abs=1e-10)
    report, _ = curvature_report(native_bundle, tmp_path)
    assert report == read(evidence / "projected-Hessian.json")
    assert report["full_molecular_minimum_verified"] is False


@pytest.mark.parametrize("mutation", ["geometry", "units", "mapping", "convergence"])
def test_real_native_mapping_and_units_reject(native_bundle, tmp_path, mutation):
    import shutil

    from island.pcff_completion import write
    from island.pcff_qm_pilot.convergence_importer import native

    key = "curvature/h0.002-d0-s+1"
    source = native_bundle / "outputs" / key
    target = tmp_path / "changed"
    shutil.copytree(source, target)
    request = read(native_bundle / "requests" / (key + ".json"))
    result = read(target / "result.json")
    if mutation == "geometry":
        result["geometry_bohr"][0][0] += 0.001
    elif mutation == "units":
        result["units"]["gradient"] = "kJ/mol/nm"
    elif mutation == "mapping":
        result["elements"][0], result["elements"][1] = (
            result["elements"][1],
            result["elements"][0],
        )
    else:
        result["scf_converged"] = False
    write(target / "result.json", result)
    with pytest.raises(CompletionError):
        native(target, request)


def test_real_dataset_fresh_process_relocation(tmp_path):
    import os
    import shutil
    import subprocess
    import sys

    root = tmp_path / "relocated evidence with spaces"
    shutil.copytree(REPO / "docs/evidence/phase_4j61", root)
    env = {
        **os.environ,
        "PYTHONPATH": str(REPO) + ":" + str(REPO / "src"),
        "OPENBLAS_NUM_THREADS": "1",
    }
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/replay_phase_4j61.py"),
            "replay",
            "--repo",
            str(REPO),
            "--evidence",
            str(root),
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(result.stdout)
    assert report["native_samples_reconstructed"] == 76
    assert report["fit_reconstructed"] and not report["QM_invoked"]


def test_reference_protocol_requires_explicit_consistent_SCF_evidence():
    from scripts.replay_phase_4j61 import bound_protocol

    old = read(OLD / "stages/tight-d-joint-fit/normalized-reference.json")
    sample = next(iter(old["samples"].values()))
    repeated = copy.deepcopy(sample)
    repeated["protocol"]["scf_type"] = sample["protocol"]["settings"]["scf_type"]
    assert bound_protocol(sample) == bound_protocol(repeated)
    repeated["protocol"]["scf_type"] = "pk"
    with pytest.raises(CompletionError, match="conflicting_reference_SCF_type"):
        bound_protocol(repeated)
    unbound = copy.deepcopy(sample)
    del unbound["protocol"]["settings"]["scf_type"]
    with pytest.raises(CompletionError, match="unbound_reference_SCF_type"):
        bound_protocol(unbound)


def test_manufactured_outside_K2_reports_active_bound(manufactured):
    template, values = manufactured
    true = values.copy()
    true[2] = 1100  # synthetic data deliberately outside the declared fit domain
    records = []
    for original, configuration in zip(
        template.training, template.configurations, strict=True
    ):
        result = configuration.evaluate(true)
        records.append(
            {
                **original,
                "energy": result["energy"],
                "forces": result["forces"].tolist(),
            }
        )
    stationary_result = template.stationary.evaluate(true)
    center = {
        **template.center,
        "energy": stationary_result["energy"],
        "forces": stationary_result["forces"].tolist(),
    }
    curvature, _ = template.mm_curvature(true)
    problem = Problem(template.model, records, center, template.q, curvature)
    initial = true.copy()
    initial[2] = 990
    result = problem.fit(initial)
    assert "K2_SiO_kcal_mol_A2" in result["active_bounds"]
    assert result["parameters"]["K2_SiO_kcal_mol_A2"] == pytest.approx(1000, abs=1e-5)
    assert result["weighted_residual_SS"] > 1


@pytest.mark.parametrize("transform", ["proper_rotation", "end_reversal"])
def test_nonzero_baseline_energy_force_covariance(model, center, transform):
    values = parameters(model)
    original = Configuration(model, center).evaluate(values)
    if transform == "proper_rotation":
        rotation, _ = np.linalg.qr(np.random.default_rng(61).normal(size=(3, 3)))
        assert np.linalg.det(rotation) > 0
        x = center @ rotation + [3.8, -1.7, 6.2]
        expected_forces = original["forces"] @ rotation
        expected_mixed = np.einsum(
            "ijp,jk->ikp", original["force_parameter_jacobian"], rotation
        )
    else:
        order = [2, 1, 0, 6, 8, 7, 3, 5, 4]
        x = center[order]
        expected_forces = original["forces"][order]
        expected_mixed = original["force_parameter_jacobian"][order]
    changed = Configuration(model, x).evaluate(values)
    assert changed["energy"] == pytest.approx(original["energy"], abs=1e-10)
    assert changed["forces"] == pytest.approx(expected_forces, abs=1e-9)
    assert changed["energy_parameter_jacobian"] == pytest.approx(
        original["energy_parameter_jacobian"], abs=1e-9
    )
    assert changed["force_parameter_jacobian"] == pytest.approx(
        expected_mixed, abs=1e-9
    )
