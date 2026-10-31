"""Numerical controls are manufactured; physical references remain native evidence."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError
from island.pcff_qm_pilot.angular_baseline import Configuration
from island.pcff_qm_pilot.angular_fit import Problem
from island.pcff_qm_pilot.molecular_control import (
    cartesian,
    nonstretch_projection,
    response,
    rigid_projection,
    torsion,
)
from scripts.audit_phase_4j62 import build
from scripts.independent_phase_4j60 import scalar

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "docs/evidence/phase_4j62"
OLD = REPO / "docs/evidence/phase_4j61"


def read(p):
    return json.loads(Path(p).read_text())


@pytest.fixture(scope="module")
def numerical():
    base, _ = build(REPO)
    model = base.model
    x = np.array(base.center["coordinates"])
    c = Configuration(model, x)
    source = read(OLD / "trials/baseline-joint.json")["parameters"]
    p = np.r_[[source[k] for k in c.names[:24]], [0.07, 0.009, -0.006, 1.3]]
    p[4:24] = np.linspace(-1.7, 2.3, 20)
    return base, model, x, c, p


def changed_model(model, p):
    from island.pcff_qm_pilot.stretch_baseline import bond_kind

    m = copy.deepcopy(model)
    for r in m["interaction_inventory"]:
        a = r["assignment"]
        f = r["family"]
        if f == "quartic_bond":
            k = bond_kind(m, r["sites"])
            a["values"][0] = p[k]
            a["values"][1] = p[k + 2] * 4.184
        elif f == "quartic_angle":
            labels = tuple(m["typing"][str(i)] for i in r["sites"])
            labels = min(labels, labels[::-1])
            if labels == ("hsi", "sio", "hsi"):
                a["values"][0] += p[25]
            elif labels == ("hsi", "sio", "osi"):
                a["values"][0] += p[26]
            elif labels == ("sio", "osi", "sio"):
                a["values"][1] += p[27] * 4.184
        elif f == "torsion_3":
            a["values"][4] += p[24] * 4.184
    return m


def independent(model, x, p):
    e, a = scalar(changed_model(model, p), x)
    return e + sum(p[4 + i] * a[k] for i, k in enumerate(sorted(a)))


@pytest.mark.parametrize("index", range(28))
def test_exact_energy_force_parameter_derivatives(numerical, index):
    _, model, x, c, p = numerical
    r = c.evaluate(p)
    h = 1e-6 if index in (0, 1, 25, 26) else 1e-3
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
    assert r["energy"] == pytest.approx(independent(model, x, p), abs=1e-10)


@pytest.mark.parametrize("step", [1e-4, 5e-5])
def test_independent_Cartesian_force_sign_units_and_propagation(numerical, step):
    _, model, x, c, p = numerical
    r = c.evaluate(p)
    f = np.zeros((9, 3))
    for i in range(9):
        for j in range(3):
            a, b = x.copy(), x.copy()
            a[i, j] += step
            b[i, j] -= step
            f[i, j] = -(independent(model, a, p) - independent(model, b, p)) / (
                2 * step
            )
    assert r["forces"] == pytest.approx(f, abs=3e-5)


@pytest.mark.parametrize(
    "family", ["angle-angle", "angle-torsion_3", "angle-angle-torsion_1"]
)
def test_angular_equilibria_propagate_fitted_cross_families(numerical, family):
    _, model, x, c, p = numerical
    p = p.copy()
    p[4:24] = 0
    a = p.copy()
    b = p.copy()
    a[25:27] = 0
    indexes = [i for i, k in enumerate(c.names) if k.startswith(family + ":")]
    assert indexes
    for i in indexes:
        p[i] = a[i] = b[i] = 2.7
    for v in (a, b):
        assert c.evaluate(v)["energy"] == pytest.approx(
            independent(model, x, v), abs=1e-10
        )
    assert (
        max(
            abs(
                c.evaluate(a)["energy_parameter_jacobian"][indexes]
                - c.evaluate(b)["energy_parameter_jacobian"][indexes]
            )
        )
        > 1e-6
    )


def test_fixed_bond_angle_and_quartic_angle_propagate(numerical):
    from scripts.independent_phase_4j59 import source_energy

    _, model, x, c, p = numerical
    a = p.copy()
    a[25:27] = 0
    before = source_energy(changed_model(model, a)["interaction_inventory"], x)
    after = source_energy(changed_model(model, p)["interaction_inventory"], x)
    for family in ("bond-angle", "quartic_angle"):
        ids = [r["id"] for r in model["interaction_inventory"] if r["family"] == family]
        assert abs(sum(before[k] - after[k] for k in ids)) > 1e-5
    assert c.evaluate(p)["energy"] == pytest.approx(independent(model, x, p), abs=1e-10)


def test_nonstretch_projection_rank_rigid_and_stretch(numerical):
    _, _, x, _, _ = numerical
    qs = read(OLD / "stretching-subspace.json")["Q"]
    seeds = read(ROOT / "plan/scientific-plan.json")["control"]["seeds"]
    r = nonstretch_projection(x, qs, seeds)
    q = r["Q"]
    assert q.T @ q == pytest.approx(np.eye(4), abs=1e-12)
    assert np.array(qs).T @ q == pytest.approx(np.zeros((8, 4)), abs=1e-12)
    assert r["rigid_basis"].T @ q == pytest.approx(np.zeros((6, 4)), abs=1e-12)
    with pytest.raises(CompletionError, match="rank_deficiency"):
        nonstretch_projection(x, qs, [seeds[0], seeds[0]])


def test_general_dimension_torsion_gradient_and_rigid_covariance():
    x = np.array(read(ROOT / "HMDS-domain-audit.json")["coordinates"])
    sites = [3, 2, 1, 10]
    t = torsion(cartesian(x), sites)
    fd = np.zeros(x.size)
    for i in range(x.size):
        a, b = x.copy(), x.copy()
        a.flat[i] += 1e-6
        b.flat[i] -= 1e-6
        delta = torsion(cartesian(a), sites).v - torsion(cartesian(b), sites).v
        fd[i] = np.arctan2(np.sin(delta), np.cos(delta)) / 2e-6
    assert t.d == pytest.approx(fd, abs=1e-8)
    u, p = rigid_projection(x)
    assert np.max(abs(u.T @ t.d)) < 1e-10
    rot, _ = np.linalg.qr(np.random.default_rng(62).normal(size=(3, 3)))
    assert np.linalg.det(rot) > 0
    new = torsion(cartesian(x @ rot + [2, -1, 3]), sites)
    assert new.v == pytest.approx(t.v, abs=1e-12)
    assert new.d.reshape(27, 3) == pytest.approx(t.d.reshape(27, 3) @ rot, abs=1e-10)
    assert np.linalg.matrix_rank(p, tol=1e-10) == 75


def test_response_sign_and_mixed_block_symmetry(numerical):
    _, _, x, _, _ = numerical
    q = np.column_stack(
        [
            read(OLD / "stretching-subspace.json")["Q"],
            read(ROOT / "nonstretch-projection.json")["Q"],
        ]
    )
    rng = np.random.default_rng(62)
    h = rng.normal(size=(27, 27))
    h = h + h.T
    responses = []
    for step in (0.004, 0.002):
        plus = [(-h @ (x.ravel() + step * q[:, i])).reshape(9, 3) for i in range(12)]
        minus = [(-h @ (x.ravel() - step * q[:, i])).reshape(9, 3) for i in range(12)]
        r = response(plus, minus, step)
        assert r == pytest.approx(h @ q, abs=2e-10)
        responses.append(r)
    block = q.T @ responses[-1]
    assert block[:8, 8:] == pytest.approx(block[8:, :8].T, abs=1e-10)


@pytest.fixture(scope="module")
def manufactured(numerical):
    base, model, _, _, truth = numerical
    base = copy.copy(base)
    base.training = copy.deepcopy(base.training[:10])
    for s in base.training:
        r = Configuration(model, s["coordinates"]).evaluate(truth)
        s["energy"] = r["energy"] + 31.2
        s["forces"] = r["forces"].tolist()
    base.center = copy.deepcopy(base.center)
    r = Configuration(model, base.center["coordinates"]).evaluate(truth)
    base.center["forces"] = r["forces"].tolist()
    base.center["energy"] = r["energy"] + 31.2
    q = np.column_stack(
        [
            read(OLD / "stretching-subspace.json")["Q"],
            read(ROOT / "nonstretch-projection.json")["Q"],
        ]
    )
    p = Problem(base, q, np.eye(12), np.zeros((27, 4)))
    h, _, response_ = p.curvature(truth)
    p.target = h
    p.full_response = response_[:, 8:]
    return p, truth


def test_manufactured_nonzero_28_parameter_recovery(manufactured):
    problem, truth = manufactured
    initial = truth.copy()
    initial[:2] += [0.003, -0.002]
    initial[2:4] *= 1.03
    initial[4:24] += 0.2
    initial[24:] += [0.03, -0.002, 0.001, 0.2]
    r = problem.fit(initial, 28)
    recovered = np.array([r["parameters"][k] for k in problem.names])
    assert recovered == pytest.approx(truth, abs=2e-5)
    assert r["rank"] == 28
    assert r["weighted_SS"] < 1e-12


def test_combined_curvature_parameter_derivatives(manufactured):
    p, v = manufactured
    _, j, _ = p.curvature(v)
    for k in (0, 24, 25, 26, 27):
        h = 1e-5 if k in (0, 25, 26) else 1e-3
        a, b = v.copy(), v.copy()
        a[k] += h
        b[k] -= h
        ha = p.curvature(a)[0]
        hb = p.curvature(b)[0]
        assert (ha - hb) / 2 / h == pytest.approx(j[:, :, k], abs=1e-4)


def test_nested_null_extension_and_bounds(numerical, manufactured):
    _, _, _, c, p = numerical
    v = p.copy()
    v[24:] = 0
    new = c.evaluate(v)
    old = c.base.evaluate(v[:24])
    assert new["energy"] == pytest.approx(old["energy"], abs=1e-12)
    assert new["forces"] == pytest.approx(old["forces"], abs=1e-12)
    problem, truth = manufactured
    bad = truth.copy()
    bad[25] = 0.2
    with pytest.raises(CompletionError, match="bounds"):
        problem.fit(bad, 28)
    with pytest.raises(CompletionError, match="unregistered"):
        problem.fit(truth, 26)


def test_HMDS_graph_mapping_and_no_SiH_transfer():
    r = read(ROOT / "HMDS-domain-audit.json")
    assert r["atom_count"] == 27
    assert r["atom_order"] == list(range(1, 28))
    assert (
        r["graph_identity"]
        == "d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd"
    )
    assert r["formal_charge"] == 0
    assert len(r["bonds"]) == 26
    assert r["parameter_completeness"] == "incomplete"
    assert not r["SiH_parameters_transferred"]
    assert abs(r["charge_sum"]) < 1e-12


def test_failed_subprocess_retains_command_streams_and_artifacts(tmp_path):
    import sys

    from island.pcff_completion import write
    from scripts.phase_4j62 import run_process

    work = tmp_path / "work"
    (work / "requests").mkdir(parents=True)
    write(work / "requests/test.json", {"mock": "test"})
    fake = tmp_path / "fake-python"
    fake.write_text(
        "#!"
        + sys.executable
        + '\nimport sys\nfrom pathlib import Path\nPath("native-failure.out").write_text("SCF failure evidence")\nprint("native stdout")\nprint("native stderr",file=sys.stderr)\nsys.exit(7)\n'
    )
    fake.chmod(0o755)
    with pytest.raises(CompletionError, match="native_subprocess_failure_retained"):
        run_process(work, REPO, fake, "test")
    root = work / "processes/test"
    assert read(root / "result.json")["returncode"] == 7
    assert (root / "command.json").exists()
    assert "native stdout" in (root / "stdout.log").read_text()
    assert "native stderr" in (root / "stderr.log").read_text()
    assert (
        work / "scratch/test/native-failure.out"
    ).read_text() == "SCF failure evidence"


def test_actual_combined_rank_loss_rejects_acceptance(manufactured):
    problem, values = manufactured
    reduced = copy.copy(problem)
    reduced._last = None
    original = reduced.system

    def only_stationary(p):
        residual, jac = original(p)
        return residual[280:307], jac[280:307]

    reduced.system = only_stationary
    trial = reduced.fit(values, 28, maximum_evaluations=2)
    assert trial["rank"] < 28
    assert not trial["training_gates"]["full_rank_and_condition"]
    assert not trial["training_qualifies"]


def test_reserved_fresh_selection_has_no_symmetry_leakage():
    from island.pcff_qm_pilot.stretch_fit import leakage

    fresh = read(ROOT / "plan/conditional-fresh-validation.json")
    old = read(REPO / "docs/evidence/phase_4j60/joint-target-design.json")
    history = (
        old["selected_training"]
        + old["held_out"]
        + read(OLD / "plan/fresh-configurations.json")["selections"]
    )
    assert leakage(history, fresh["selections"]) == fresh["separation"]
    duplicate = copy.deepcopy(fresh["selections"])
    duplicate[1]["coordinates"] = np.array(history[0]["coordinates"])[
        [2, 1, 0, 6, 7, 8, 3, 4, 5]
    ].tolist()
    with pytest.raises(CompletionError, match="symmetry_leakage"):
        leakage(history, duplicate)


@pytest.fixture(scope="module")
def retained_native(tmp_path_factory):
    from island.pcff_qm_pilot.bundle import unpack

    work = tmp_path_factory.mktemp("retained-native") / "native"
    unpack(ROOT / "native-outputs.tar.xz", work)
    return work


@pytest.mark.parametrize(
    "key,n", [("control/h0.002-d0-s+1", 9), ("hmds-opt/step-000", 27)]
)
def test_real_native_dimension_conversion_and_independent_reconstruction(
    retained_native, key, n
):
    from island.pcff_qm_pilot.molecular_native import normalized
    from scripts.analyze_phase_4j62 import independent

    d = retained_native / "outputs" / key
    r = read(retained_native / "requests" / (key + ".json"))
    s = normalized(d, r)
    energy, force = independent(d, n)
    assert s["energy"] == pytest.approx(energy, abs=2e-7)
    assert np.array(s["forces"]).shape == (n, 3)
    assert s["forces"] == pytest.approx(force, abs=1e-10)


@pytest.mark.parametrize(
    "mutation", ["units", "order", "graph", "charge", "scf", "gradient", "output"]
)
def test_real_native_mutations_fail_closed(retained_native, tmp_path, mutation):
    import shutil

    from island.pcff_completion import pin, write
    from island.pcff_qm_pilot.molecular_native import normalized

    key = "hmds-opt/step-000"
    directory = tmp_path / "native"
    shutil.copytree(retained_native / "outputs" / key, directory)
    request = read(retained_native / "requests" / (key + ".json"))
    result = read(directory / "result.json")
    if mutation == "units":
        result["units"]["gradient"] = "kcal/mol/angstrom"
    elif mutation == "order":
        request["atom_order"][0:2] = [2, 1]
    elif mutation == "graph":
        request["graph_identity"] = "0" * 64
    elif mutation == "charge":
        result["charge"] = 1
    elif mutation == "scf":
        result["scf_converged"] = False
    elif mutation == "gradient":
        g = np.load(directory / "gradient.npy")
        g[0, 0] += 0.1
        np.save(directory / "gradient.npy", g, allow_pickle=False)
        result["native_files"]["gradient.npy"] = pin(directory / "gradient.npy")
    else:
        (directory / "psi4.out").write_text("altered native output")
    write(directory / "result.json", result)
    with pytest.raises(CompletionError):
        normalized(directory, request)


def test_relocation_and_fresh_process_reconstruction(tmp_path):
    import os
    import shutil
    import subprocess
    import sys

    evidence = tmp_path / "relocated-evidence"
    shutil.copytree(ROOT, evidence)
    child = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.replay_phase_4j62",
            "replay",
            "--repo",
            str(REPO),
            "--evidence",
            str(evidence),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        env={
            **os.environ,
            "PYTHONPATH": str(REPO) + ":" + str(REPO / "src"),
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
        },
    )
    assert child.returncode == 0, child.stdout + child.stderr
    report = json.loads(child.stdout)
    assert report["native_samples_reconstructed"] == 45
    assert not report["QM_invoked"]
    assert not report["refitting_invoked"]


def test_angular_extension_rigid_energy_force_and_parameter_covariance(numerical):
    _, model, x, c, p = numerical
    rotation, _ = np.linalg.qr(np.random.default_rng(4062).normal(size=(3, 3)))
    moved = Configuration(model, x @ rotation + [3, -2, 1]).evaluate(p)
    original = c.evaluate(p)
    assert moved["energy"] == pytest.approx(original["energy"], abs=1e-10)
    assert moved["forces"] == pytest.approx(original["forces"] @ rotation, abs=1e-9)
    assert moved["energy_parameter_jacobian"] == pytest.approx(
        original["energy_parameter_jacobian"], abs=1e-9
    )
    expected = np.einsum("iap,ab->ibp", original["force_parameter_jacobian"], rotation)
    assert moved["force_parameter_jacobian"] == pytest.approx(expected, abs=1e-8)


def test_compute_ledger_counts_failed_and_interrupted_calls(tmp_path):
    import time

    from island.pcff_completion import write
    from island.pcff_qm_pilot.control_engine import reserve

    ledger = tmp_path / "ledger.json"
    write(
        ledger,
        {
            "started_unix": time.time(),
            "calls": {
                f"job{i}": {
                    "status": ("failed" if i % 2 else "reserved_or_interrupted")
                }
                for i in range(48)
            },
        },
    )
    before = ledger.read_bytes()
    with pytest.raises(CompletionError, match="campaign_budget_exhausted"):
        reserve(tmp_path, "new-job", {"manufactured": "budget-test"})
    assert ledger.read_bytes() == before


def test_development_validation_values_do_not_enter_training_selection(
    numerical, monkeypatch
):
    base, _, _, _, _ = numerical
    _, parameters = build(REPO)
    expected = base.system(parameters)[0].copy()
    original = json.loads
    touched = []

    def changed(text, *args, **kwargs):
        record = original(text, *args, **kwargs)
        if isinstance(record, dict) and isinstance(record.get("samples"), dict):
            for key, sample in record["samples"].items():
                if key.startswith("fresh/validation-"):
                    sample["energy"] += 1e7
                    sample["forces"] = (np.array(sample["forces"]) + 1e6).tolist()
                    touched.append(key)
        return record

    monkeypatch.setattr(json, "loads", changed)
    rebuilt, values = build(REPO)
    assert len(touched) == 4
    assert values == pytest.approx(parameters, abs=0)
    assert rebuilt.system(values)[0] == pytest.approx(expected, abs=0)
