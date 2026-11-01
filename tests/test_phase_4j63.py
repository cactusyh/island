"""Numerical HMDS kernels, graph paths, and fail-visible execution contracts."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, signed
from island.pcff_qm_pilot.hmds_engine import optimization_exit, reserve
from island.pcff_qm_pilot.hmds_hamiltonian import (
    _kernels,
    evaluate,
    load_model,
    parts,
)
from island.pcff_qm_pilot.hmds_independent import energy, scalar
from island.pcff_qm_pilot.hmds_paths import (
    affected,
    automorphism_count,
    equivalent,
    leakage,
    path,
    response,
)

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "docs/evidence/phase_4j63"


def read(p):
    return json.loads(p.read_text())


@pytest.fixture(scope="module")
def model():
    return load_model(REPO)


@pytest.fixture(scope="module")
def xyz():
    return np.array(
        read(REPO / "docs/evidence/phase_4j62/HMDS-center.json")["center"][
            "coordinates"
        ]
    )


@pytest.fixture(scope="module")
def manufactured(model, xyz):
    p = parts(model, xyz)
    rng = np.random.default_rng(4063)
    return dict(
        zip(p["columns"], rng.uniform(0.02, 0.5, len(p["columns"])), strict=True)
    )


def test_independent_scalar_and_source_symmetry_reduction(model, xyz):
    p = parts(model, xyz)
    q = scalar(model, xyz)
    assert len(p["columns"]) == 56
    assert len({r["directed_query"] for r in p["query_to_columns"]}) == 26
    assert set(q["columns"]) == set(p["columns"])
    assert q["fixed_energy"] == pytest.approx(p["fixed_energy"], abs=1e-10)
    assert p["energy_jacobian"] == pytest.approx(
        [q["columns"][n] for n in p["columns"]], abs=1e-11
    )
    for k, v in p["components"].items():
        assert v == pytest.approx(q["components"][k], abs=1e-10)


@pytest.mark.parametrize("h", [2e-5, 1e-5])
def test_independent_full_Cartesian_force_and_parameter_jacobian(
    model, xyz, manufactured, h
):
    p = parts(model, xyz)
    out = evaluate(model, xyz, manufactured)
    fd = np.zeros((27, 3))
    jac = np.zeros((27, 3, 56))
    for i, j in np.ndindex(27, 3):
        a, b = xyz.copy(), xyz.copy()
        a[i, j] += h
        b[i, j] -= h
        sa, sb = scalar(model, a), scalar(model, b)
        fd[i, j] = -(
            energy(model, a, manufactured) - energy(model, b, manufactured)
        ) / (2 * h)
        jac[i, j] = [
            -(sa["columns"][n] - sb["columns"][n]) / (2 * h) for n in p["columns"]
        ]
    assert out["forces"] == pytest.approx(fd, abs=3e-5)
    assert p["force_jacobian"] == pytest.approx(jac, abs=3e-7)


@pytest.mark.parametrize("index", range(56))
def test_energy_force_parameter_derivatives(model, xyz, manufactured, index):
    p = parts(model, xyz)
    name = p["columns"][index]
    a, b = manufactured.copy(), manufactured.copy()
    h = 1e-4
    a[name] += h
    b[name] -= h
    ea, eb = evaluate(model, xyz, a), evaluate(model, xyz, b)
    assert (ea["energy"] - eb["energy"]) / (2 * h) == pytest.approx(
        p["energy_jacobian"][index], abs=1e-8
    )
    assert (ea["forces"] - eb["forces"]) / (2 * h) == pytest.approx(
        p["force_jacobian"][:, :, index], abs=1e-8
    )


@pytest.mark.parametrize("spec", read(ROOT / "plan/scientific-plan.json")["paths"])
def test_graph_fragment_path_and_tangent_acceleration(model, xyz, spec):
    h = 0.003 if spec["kind"] == "stretch" else 0.02
    p0 = path(model, xyz, spec, 0)
    p = path(model, xyz, spec, h)
    m = path(model, xyz, spec, -h)
    eps = 1e-5
    a = path(model, xyz, spec, h + eps)
    b = path(model, xyz, spec, h - eps)
    assert p["tangent"] == pytest.approx(
        (a["coordinates"] - b["coordinates"]) / (2 * eps), abs=1e-9
    )
    assert p["acceleration"] == pytest.approx(
        (a["tangent"] - b["tangent"]) / (2 * eps), abs=1e-9
    )
    changes = affected(model, xyz, p["coordinates"])
    if spec["kind"] == "rotation":
        assert all(
            r["family"] == "torsion_3" and set(r["sites"][1:3]) == set(spec["bond"])
            for r in changes
        )
        assert changes
        for r in changes:
            assert abs(r["change"]) == pytest.approx(h, abs=1e-10)
    else:
        assert len(changes) == 1 and changes[0]["family"] == "quartic_bond"
        assert changes[0]["change"] == pytest.approx(h, abs=1e-10)
    # Reverse rotation returns the exact starting geometry along the same fixed axis.
    if spec["kind"] == "rotation":
        assert path(model, p["coordinates"], spec, -h)["coordinates"] == pytest.approx(
            xyz, abs=1e-12
        )
    assert set(p0["moving_ids"]) == set(p["moving_ids"]) == set(m["moving_ids"])


@pytest.mark.parametrize("spec", read(ROOT / "plan/scientific-plan.json")["paths"])
def test_path_chain_rule_force_sign_units(model, xyz, manufactured, spec):
    h = 1e-4
    p0 = path(model, xyz, spec, 0)
    pp = path(model, xyz, spec, h)
    pm = path(model, xyz, spec, -h)
    out = []
    for p in (p0, pp, pm):
        out.append(evaluate(model, p["coordinates"], manufactured))
    r = response(out[0], out[1], out[2], p0, pp, pm, h)
    fd = (out[1]["energy"] - out[2]["energy"]) / (2 * h)
    assert -np.sum(out[0]["forces"] * p0["tangent"]) == pytest.approx(fd, abs=2e-5)
    assert r["curvature_from_path_gradient"] == pytest.approx(
        r["curvature_from_energies"], abs=0.002
    )
    assert r["curvature_from_path_gradient"] == pytest.approx(
        r["chain_rule_sum"], abs=0.001
    )
    if spec["kind"] == "rotation":
        assert abs(r["gradient_path_acceleration"]) > 1e-4
    else:
        assert r["gradient_path_acceleration"] == 0


def test_manufactured_recovery_rigid_invariance_and_graph_symmetry(
    model, xyz, manufactured
):
    rng = np.random.default_rng(163)
    rows = []
    target = []
    for _ in range(6):
        xx = xyz + rng.normal(0, 0.03, xyz.shape)
        p = parts(model, xx)
        e = evaluate(model, xx, manufactured)
        rows.extend([p["energy_jacobian"], *p["force_jacobian"].reshape(81, 56)])
        target.extend(
            [
                e["energy"] - p["fixed_energy"],
                *(e["forces"] - p["fixed_forces"]).ravel(),
            ]
        )
    a = np.array(rows)
    v, _, rank, _ = np.linalg.lstsq(a, np.array(target), rcond=1e-12)
    assert rank == 56
    assert v == pytest.approx(list(manufactured.values()), abs=2e-9)
    r, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(r) < 0:
        r[:, 0] *= -1
    before = evaluate(model, xyz, manufactured)
    after = evaluate(model, xyz @ r + [3, -2, 4], manufactured)
    assert after["energy"] == pytest.approx(before["energy"], abs=1e-10)
    assert after["forces"] == pytest.approx(before["forces"] @ r, abs=1e-9)


@pytest.mark.parametrize(
    "family",
    [
        "angle-angle",
        "angle-torsion_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-angle-torsion_1",
        "bond-bond_1_3",
    ],
)
def test_role_permutation_and_torsion_reversal(model, xyz, family):
    altered = copy.deepcopy(model)
    for row in altered["interaction_inventory"]:
        if row["family"] != family:
            continue
        s = row["sites"]
        row["sites"] = [s[3], s[1], s[2], s[0]] if family == "angle-angle" else s[::-1]
    a, b = _kernels(model, xyz), _kernels(altered, xyz)
    assert a["columns"] == b["columns"]
    assert a["energy_jacobian"] == pytest.approx(b["energy_jacobian"], abs=1e-12)
    assert a["force_jacobian"] == pytest.approx(b["force_jacobian"], abs=1e-11)
    assert a["fixed_energy"] == pytest.approx(b["fixed_energy"], abs=1e-11)


def test_all_equilibrium_dependencies_propagate(model, xyz):
    altered = copy.deepcopy(model)
    for row in altered["interaction_inventory"]:
        if row["family"] == "quartic_bond":
            row["assignment"]["values"][0] += 0.012
        elif row["family"] == "quartic_angle":
            row["assignment"]["values"][0] += 0.015
    a, b = _kernels(model, xyz), _kernels(altered, xyz)
    s = scalar(altered, xyz)
    assert b["fixed_energy"] == pytest.approx(s["fixed_energy"], abs=1e-10)
    assert b["energy_jacobian"] == pytest.approx(
        [s["columns"][n] for n in b["columns"]], abs=1e-11
    )
    for f in (
        "angle-angle",
        "angle-torsion_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-angle-torsion_1",
        "bond-bond_1_3",
    ):
        indices = [i for i, n in enumerate(a["columns"]) if n.split(":")[0] == f]
        assert (
            np.linalg.norm(
                a["energy_jacobian"][indices] - b["energy_jacobian"][indices]
            )
            > 1e-5
        )


@pytest.mark.parametrize(
    "mutation",
    ["source", "charge", "typing", "units", "graph", "conventions", "coefficient"],
)
def test_source_mutation_and_unsupported_domain_rejected(model, xyz, mutation):
    m = copy.deepcopy(model)
    if mutation == "source":
        m["source_record_sha256"] = "0" * 64
    elif mutation == "charge":
        m["charges"]["1"] += 0.001
    elif mutation == "typing":
        m["typing"]["1"] = "hsi"
    elif mutation == "units":
        m["units"]["length"] = "bohr"
    elif mutation == "graph":
        m["bonds"][0][2] = 2
    elif mutation == "conventions":
        m["conventions"]["pairs"] = "full all pairs"
    else:
        next(r for r in m["interaction_inventory"] if r["family"] == "quartic_bond")[
            "assignment"
        ]["values"][1] *= 1.01
    m = signed({k: v for k, v in m.items() if k != "identity"})
    with pytest.raises(CompletionError):
        parts(m, xyz)


def test_unknowns_are_required_and_pairs_are_bound(model, xyz, manufactured):
    with pytest.raises(CompletionError, match="explicit_all"):
        evaluate(model, xyz)
    missing = manufactured.copy()
    del missing[next(iter(missing))]
    with pytest.raises(CompletionError, match="explicit_all"):
        evaluate(model, xyz, missing)
    p = parts(model, xyz)
    assert len(p["pairs"]) == 351
    for r in p["pairs"]:
        assert (
            r["LJ_weight"]
            == r["Coulomb_weight"]
            == (0 if r["graph_distance"] < 3 else 1)
        )
        if r["graph_distance"] < 3:
            assert r["LJ_energy"] == r["Coulomb_energy"] == 0


def test_scalable_methyl_H_and_whole_end_symmetry(model, xyz):
    assert automorphism_count(model) == 3359232
    # Whole ends, methyl groups, and H labels may all permute, preserving graph colors.
    permutation = [
        7,
        6,
        8,
        9,
        5,
        2,
        1,
        3,
        4,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
    ]
    y = xyz[np.array(permutation) - 1]
    assert equivalent(model, xyz, y + [4, -2, 1])
    with pytest.raises(CompletionError, match="symmetry_leakage"):
        leakage(model, [xyz], [y])
    y = xyz.copy()
    y[9] += [0.01, 0.002, 0]
    assert not equivalent(model, xyz, y)
    assert leakage(model, [xyz], [y])["reserved_count"] == 1
    with pytest.raises(CompletionError, match="inconclusive"):
        equivalent(model, xyz, xyz, maximum_states=1)


def test_complete_manufactured_model_whole_end_and_methyl_H_permutations(
    model, xyz, manufactured
):
    end = (
        np.array(
            [
                7,
                6,
                8,
                9,
                5,
                2,
                1,
                3,
                4,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
            ]
        )
        - 1
    )
    hydrogen = np.arange(27)
    hydrogen[[9, 10, 11]] = [11, 9, 10]
    before = evaluate(model, xyz, manufactured)
    for permutation in (end, hydrogen, end[hydrogen]):
        after = evaluate(model, xyz[permutation], manufactured)
        assert after["energy"] == pytest.approx(before["energy"], abs=1e-10)
        assert after["forces"] == pytest.approx(
            before["forces"][permutation], abs=1e-10
        )


def test_optimizer_exit_status_and_failed_attempt_budget(tmp_path):
    assert optimization_exit("converged") == 0
    assert optimization_exit("budget_exhausted") == 3
    assert optimization_exit("unexpected_failure") == 1
    for i in range(80):
        reserve(tmp_path, str(i), {"coordinate": i})
    assert len(read(tmp_path / "ledger.json")["calls"]) == 80
    with pytest.raises(CompletionError, match="budget_exhausted"):
        reserve(tmp_path, "overflow", {})
    with pytest.raises(CompletionError, match="retry_identity"):
        reserve(tmp_path, "0", {})


def test_callback_cap_does_not_relabel_unexpected_metadata_failure():
    from island.pcff_qm_pilot.hmds_engine import optimization_budget_exception

    assert optimization_budget_exception(
        CompletionError("J63_optimization_callback_cap")
    )
    assert not optimization_budget_exception(
        OSError("manufactured history write failure")
    )
    assert not optimization_budget_exception(
        RuntimeError("manufactured native failure")
    )


@pytest.mark.parametrize(
    "iteration,status,message,expected",
    [
        (40, "FAILED", "Maximum number of steps exceeded: 40.", True),
        (40, "FAILED", "manufactured coordinate algorithm failure", False),
        (41, "CONTINUE", None, True),
        (41, "FAILED", None, False),
    ],
)
def test_native_iteration_cap_and_other_optimizer_failure_are_distinct(
    iteration, status, message, expected
):
    from types import FunctionType

    from island.pcff_qm_pilot.hmds_engine import optimization_budget_exception

    # Manufactured exception/frame control; no native engine or numerical QM data.
    exception_class = type(
        "OptimizationConvergenceError",
        (RuntimeError,),
        {"__module__": "psi4.driver.p4util.exceptions"},
    )

    def optimize(iteration, status, data, error):
        n, opt_status, opt_data = iteration, status, data
        assert n in (40, 41) and opt_status in ("FAILED", "CONTINUE")
        assert isinstance(opt_data, dict)
        raise error

    manufactured = FunctionType(optimize.__code__, {"__name__": "psi4.driver.driver"})
    with pytest.raises(exception_class) as caught:
        manufactured(
            iteration,
            status,
            {"error": {"error_message": message}},
            exception_class("manufactured optimizer status"),
        )
    assert optimization_budget_exception(caught.value) is expected


def test_torsion_phase_periodicity_and_source_executable_distinction(model, xyz):
    for phi in np.linspace(-np.pi, np.pi, 19):
        for n in (1, 2, 3):
            v = -1.3513
            assert v * (1 + np.cos(n * phi)) == pytest.approx(
                v * (1 - np.cos(n * phi - np.pi)), abs=1e-14
            )
            assert np.cos(n * (phi + 2 * np.pi)) == pytest.approx(
                np.cos(n * phi), abs=1e-14
            )
    assert "linear-Phi" in model["conventions"]["AAT"]
    assert (
        "cos(phi)"
        in __import__(
            "island.pcff_qm_pilot.hmds_hamiltonian", fromlist=["FORMS"]
        ).FORMS["angle-angle-torsion_1"]
    )
    # Reverse the ordered main torsions, without changing source phases or values.
    m = copy.deepcopy(model)
    for r in m["interaction_inventory"]:
        if r["family"] == "torsion_3":
            r["sites"].reverse()
    a, b = _kernels(model, xyz), _kernels(m, xyz)
    assert a["fixed_energy"] == pytest.approx(b["fixed_energy"], abs=1e-11)
    assert a["fixed_forces"] == pytest.approx(b["fixed_forces"], abs=1e-10)


def test_equilibrium_change_full_force_chain_rule(model, xyz, manufactured):
    m = copy.deepcopy(model)
    for r in m["interaction_inventory"]:
        if r["family"] == "quartic_bond":
            r["assignment"]["values"][0] += 0.015
        if r["family"] == "quartic_angle":
            r["assignment"]["values"][0] += 0.025
    p = _kernels(m, xyz)
    v = np.array([manufactured[n] for n in p["columns"]])
    f = p["fixed_forces"] + p["force_jacobian"] @ v
    d = np.random.default_rng(463).normal(size=xyz.shape)
    d /= np.linalg.norm(d)
    h = 1e-5
    assert -np.sum(f * d) == pytest.approx(
        (energy(m, xyz + h * d, manufactured) - energy(m, xyz - h * d, manufactured))
        / (2 * h),
        abs=1e-6,
    )


def test_unexpected_subprocess_failure_is_nonzero_and_artifacts_survive(tmp_path):
    from island.pcff_completion import write
    from scripts.phase_4j63 import run_process

    work = tmp_path / "work"
    work.mkdir()
    write(work / "scientific-plan.json", read(ROOT / "plan/scientific-plan.json"))
    executable = tmp_path / "fake-engine"
    executable.write_text(
        '#!/usr/bin/env python3\nimport sys,pathlib\nw=pathlib.Path(sys.argv[sys.argv.index("--work")+1]);k=sys.argv[sys.argv.index("--key")+1]\np=w/"outputs"/k;p.mkdir(parents=True);(p/"native-failure.out").write_text("intentional test failure; not physical QM")\nprint("native failed");print("diagnostic stderr",file=sys.stderr);sys.exit(1)\n'
    )
    executable.chmod(0o755)
    assert run_process(REPO, work, executable, "negative-control") == 1
    record = read(work / "processes/negative-control/result.json")
    assert record["returncode"] == 1 and not record["cleanup_performed"]
    assert (work / "outputs/negative-control/native-failure.out").is_file()
    assert (
        "native failed" in (work / "processes/negative-control/stdout.log").read_text()
    )
    assert (
        "diagnostic stderr"
        in (work / "processes/negative-control/stderr.log").read_text()
    )
    assert read(work / "processes/negative-control/command.json")["argv"][0] == str(
        executable
    )


def test_near_duplicate_pruning_respects_coordinate_RMS(model, xyz):
    other = xyz.copy()
    other[9, 0] += 5e-6
    assert np.sqrt(np.mean((other - xyz) ** 2)) < 1e-6
    assert equivalent(model, xyz, other, tolerance=1e-6)
    with pytest.raises(CompletionError, match="validation_symmetry_leakage"):
        leakage(model, [xyz], [other])


def test_independent_graph_role_closure_rejects_omission(model):
    from island.pcff_qm_pilot.hmds_hamiltonian import graph_inventory

    r = graph_inventory(model)
    assert r["counts"]["angle-angle"] == 96
    assert r["counts"]["torsion_3"] == 60
    missing = copy.deepcopy(model)
    missing["interaction_inventory"] = [
        r
        for r in missing["interaction_inventory"]
        if r["id"] != "bond-bond_1_3:1,2,5,6"
    ]
    with pytest.raises(CompletionError, match="graph_interaction_inventory_incomplete"):
        graph_inventory(missing)


def test_timeout_collection_never_promotes_optimizer_success(tmp_path):
    from island.pcff_completion import write
    from scripts.prepare_phase_4j63 import optimizer_outcome

    write(tmp_path / "processes/optimization/result.json", {"returncode": 124})
    write(
        tmp_path / "ledger.json",
        {"calls": {"optimization/step-000": {"status": "reserved_or_interrupted"}}},
    )
    result, process = optimizer_outcome(tmp_path)
    assert result["status"] == "unexpected_failure"
    assert process["returncode"] == 124
    assert not result["native_terminal_record_supplied"]
    write(tmp_path / "processes/optimization/result.json", {"returncode": 0})
    with pytest.raises(CompletionError, match="success_without_native_terminal_record"):
        optimizer_outcome(tmp_path)


def test_joint_design_rank_scaling_and_duplicate_rows():
    from island.pcff_qm_pilot.hmds_design import metrics

    a = np.random.default_rng(1463).normal(size=(150, 56))
    assert metrics(a)["rank"] == 56
    scaled = a * np.geomspace(1e-7, 1e7, 56)
    assert metrics(a)["scaled_condition"] == pytest.approx(
        metrics(scaled)["scaled_condition"]
    )
    a[:, -1] = a[:, 0]
    assert metrics(a)["rank"] == 55
    assert not metrics(np.vstack([a, a]))["identifiable_for_screen"]


def test_response_analysis_manufactured_protocol_control_is_path_specific(model, xyz):
    from scripts.analyze_phase_4j63 import analyze

    plan = read(ROOT / "plan/scientific-plan.json")

    def artificial(x):
        displacement = np.array(x) - xyz
        return {
            "coordinates": np.array(x).tolist(),
            "energy": float(0.5 * np.sum(displacement**2)),
            "forces": (-displacement).tolist(),
            "protocol": {"manufactured": True},
        }

    center, data = artificial(xyz), {}
    for spec in plan["paths"]:
        for h in spec["steps"]:
            for sign in (1, -1):
                data[f"paths/{spec['name']}-h{h:g}-s{sign:+d}"] = artificial(
                    path(model, xyz, spec, h * sign)["coordinates"]
                )
    for name, control in plan["controls"].items():
        data[f"controls/{name}-center"] = center
        for target in control["paths"]:
            spec = next(s for s in plan["paths"] if s["name"] == target)
            for sign in (1, -1):
                data[f"controls/{name}-{target}-s{sign:+d}"] = copy.deepcopy(
                    data[f"paths/{target}-h{min(spec['steps']):g}-s{sign:+d}"]
                )
    result = analyze(model, plan, center, data)
    assert result["bounded_precision_passed"]
    bad = copy.deepcopy(data)
    bad["controls/basis-siloxane-s+1"]["forces"][0][0] += 1
    result = analyze(model, plan, center, bad)
    assert not result["controls"]["basis"]["paths"]["siloxane"]["passed"]
    assert result["controls"]["basis"]["paths"]["CH"]["passed"]
    assert not result["bounded_precision_passed"]


def test_native_convergence_coordinates_and_nonconvergence_are_distinct():
    from scripts.prepare_phase_4j63 import native_convergence

    r = native_convergence(
        "\n  3 -893.38163750 -1.23e-04 o 1.08e-03 1.18e-04 3.80e-02 1.09e-02 ~\n"
    )
    assert not any(r["native_internal_threshold_checks"].values())
    assert r["final"]["maximum_internal_force"] == 0.00108
    assert "radian" in r["units_scope"]
    with pytest.raises(CompletionError, match="convergence_table_missing"):
        native_convergence("SCF succeeded; no optimizer convergence result")


@pytest.mark.parametrize("name", ["SiC", "CH"])
@pytest.mark.parametrize("h", [0.006, 0.003])
def test_pure_stretch_total_curvature_cannot_be_repaired_by_missing_couplings(
    model, xyz, manufactured, name, h
):
    from scripts.curvature_phase_4j63 import unknown_null_check

    spec = next(
        s
        for s in read(ROOT / "plan/scientific-plan.json")["paths"]
        if s["name"] == name
    )
    r = unknown_null_check(model, xyz, spec, h)
    assert len(r["unknown_columns"]) == 56
    assert r["independent_scalar_unknown_curvature_maximum"] < 1e-7
    geo = [path(model, xyz, spec, s * h) for s in (0, 1, -1)]
    physical_test_inputs = [
        evaluate(model, g["coordinates"], manufactured) for g in geo
    ]
    total = response(*physical_test_inputs, *geo, h)
    assert total["curvature_from_path_gradient"] == pytest.approx(
        r["verified_fixed_response"]["curvature_from_path_gradient"], abs=1e-8
    )


@pytest.fixture(scope="module")
def retained_native(tmp_path_factory):
    from island.pcff_qm_pilot.bundle import unpack

    directory = tmp_path_factory.mktemp("j63-native") / "relocated"
    unpack(ROOT / "native-outputs.tar.xz", directory)
    return directory


@pytest.mark.parametrize(
    "mutation",
    [
        "charge",
        "multiplicity",
        "basis",
        "axis_order",
        "units",
        "scf_converged",
        "geometry_bohr",
        "gradient",
        "printed_energy",
        "printed_atom_order",
    ],
)
def test_real_native_intake_mutation_rejection(retained_native, tmp_path, mutation):
    import shutil

    from island.pcff_completion import pin, write
    from island.pcff_qm_pilot.molecular_native import normalized

    key = "optimization/step-000"
    directory = tmp_path / "output"
    shutil.copytree(retained_native / "outputs" / key, directory)
    request = read(retained_native / "requests" / (key + ".json"))
    r = read(directory / "result.json")
    if mutation == "geometry_bohr":
        r[mutation][0][0] += 0.01
    elif mutation == "gradient":
        g = np.load(directory / "gradient.npy")
        g[0, 0] += 1
        np.save(directory / "gradient.npy", g)
        r["native_files"]["gradient.npy"] = pin(directory / "gradient.npy")
    elif mutation.startswith("printed_"):
        p = directory / "psi4.out"
        text = p.read_text()
        if mutation == "printed_energy":
            text = text.replace("@DF-RKS Final Energy:", "@DF-RKS Final Energy: 1 ")
        else:
            import re

            text = re.sub(
                r"^(\s+)C(\s+[-+\d.]+\s+[-+\d.]+\s+[-+\d.]+\s+[-+\d.]+\s*)$",
                r"\1H\2",
                text,
                count=1,
                flags=re.MULTILINE,
            )
        assert text != p.read_text()
        p.write_text(text)
        r["native_files"]["psi4.out"] = pin(p)
    else:
        r[mutation] = {
            "charge": 1,
            "multiplicity": 3,
            "basis": "def2-svp",
            "axis_order": "component-major",
            "units": {},
            "scf_converged": False,
        }[mutation]
    write(directory / "result.json", r)
    with pytest.raises(CompletionError):
        normalized(directory, request)


def test_real_reference_relocation_fresh_process(tmp_path):
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
            "scripts.replay_phase_4j63",
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
    result = json.loads(child.stdout)
    compute = read(ROOT / "compute-accounting.json")
    assert result["native_samples_reconstructed"] == compute["completed_gradients"]
    assert (
        result["attempted_gradients"]
        == compute["attempted_energy_gradient_evaluations"]
    )
    assert not result["QM_invoked"] and not result["fit_invoked"]
