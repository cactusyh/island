"""Actual J59 source evaluator, native importer, mapping and relocation controls."""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, signed
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.convergence_importer import native, normalize
from island.pcff_qm_pilot.residual import (
    conditional_columns,
    evaluate,
    lj96,
    reconstruct,
)
from island.pcff_qm_pilot.residual_analysis import controls
from scripts.independent_phase_4j59 import numerical_force, source_energy

REPO = Path(__file__).resolve().parents[1]
EVIDENCE = REPO / "docs/evidence/phase_4j59"
DIST = (
    REPO.parent / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
)


def read(path):
    return json.loads(Path(path).read_text())


@pytest.fixture(scope="module")
def model():
    return read(EVIDENCE / "interaction-table.json")


@pytest.fixture(scope="module")
def plan():
    return read(EVIDENCE / "plan/convergence-plan.json")


@pytest.fixture(scope="module")
def distribution():
    if not (DIST / "pcff.frc").is_file():
        pytest.skip("original pinned FRC source unavailable")
    return DIST


@pytest.fixture(scope="module")
def raw(tmp_path_factory):
    path = tmp_path_factory.mktemp("j59") / "native"
    unpack(EVIDENCE / "native-outputs.tar.xz", path)
    return path


def test_analytic_small_systems():
    assert controls()["passed"]


@pytest.mark.parametrize(
    "name",
    [
        "center",
        "angle-angle-0000",
        "angle-angle-0024",
        "bond-bond-0000",
        "bond-bond-0024",
    ],
)
def test_actual_fixed_source_evaluator(model, plan, name):
    x = plan["geometries"][name]
    r = evaluate(model, x, expected_identity=model["identity"])
    independent = source_energy(model["interaction_inventory"], x)
    assert r["energy"] == pytest.approx(sum(independent.values()), abs=1e-10)
    assert (
        np.max(
            np.abs(
                np.array(r["forces"])
                - numerical_force(model["interaction_inventory"], x, 5e-5)
            )
        )
        < 2e-5
    )


def test_all_interactions_accounted_without_zero(model):
    assert model["status_counts"] == {
        "fixed": 74,
        "unresolved": 115,
        "not_applicable": 8,
    }
    assert len([r for r in model["table"] if r["family"] == "LJ9-6+Coulomb"]) == 36
    assert len([r for r in model["table"] if r["family"] == "torsion-torsion_1"]) == 9
    assert model["author_confirmations"]["author_confirmed_audited"] == 28
    assert model["author_confirmations"]["source_artifact_status"] == "not_yet_supplied"
    assert not model["unknowns_set_to_zero"]


def test_complete_residual_rejected(model, plan):
    with pytest.raises(CompletionError, match="unbound_physical_remainder"):
        evaluate(model, plan["geometries"]["center"], require_complete=True)


def test_conditional_joint_rank_is_not_physical_identifiability(model, plan):
    assert len(conditional_columns(model, plan["geometries"]["center"])) == 31
    r = read(EVIDENCE / "residual-identifiability.json")
    assert r["energy_design"]["rank"] == 22
    assert r["energy_force_design"]["rank"] == 31
    assert r["energy_force_design"]["conditioning"] > 800
    assert not r["physical_identifiability"] and not r["fitted_candidates"]


@pytest.mark.parametrize(
    "field",
    [
        "family",
        "source_hashes",
        "charge_vector",
        "graph_identity",
        "atom_order",
        "energy_unit",
        "interaction_inventory",
        "source_rows",
    ],
)
def test_source_mapping_units_mutation(model, field, distribution):
    r = copy.deepcopy(model)
    if field == "interaction_inventory":
        r[field][0]["sites"] = [9]
    elif field == "source_rows":
        next(iter(r[field].values()))["source_sha256"] = "0" * 64
    elif field == "source_hashes":
        r[field]["pcff.frc"] = "0" * 64
    elif field == "charge_vector":
        r[field]["1"] = "0.62"
    elif field == "atom_order":
        r[field][0], r[field][1] = r[field][1], r[field][0]
    else:
        r[field] = "mutated"
    r = signed({k: v for k, v in r.items() if k != "identity"})
    with pytest.raises(CompletionError):
        reconstruct(r, REPO, distribution)


def test_pair_kernel_source_units_and_mixing():
    e, d = lj96(4, 4, 4, 0.5, 0.5)
    assert e == pytest.approx(-0.5)
    assert abs(d) < 1e-12
    h = 1e-5
    assert lj96(3.8, 4, 3.2, 0.5, 0.2)[1] == pytest.approx(
        (lj96(3.8 + h, 4, 3.2, 0.5, 0.2)[0] - lj96(3.8 - h, 4, 3.2, 0.5, 0.2)[0])
        / (2 * h),
        abs=1e-8,
    )


@pytest.mark.parametrize(
    "key",
    [
        "grid/angle-angle-0000",
        "grid/center",
        "qz/bond-bond-0024",
        "scf/center",
        "pk-control/center",
        "method/angle-angle-0024",
        "thread-control/angle-angle-0000",
    ],
)
def test_executed_native_import(raw, key):
    expected = read(raw / "requests" / (key + ".json"))
    r, g = native(raw / "outputs" / key, expected)
    assert r["scf_converged"] and g.shape == (9, 3) and np.isfinite(g).all()


@pytest.mark.parametrize(
    "field",
    [
        "units",
        "elements",
        "charge",
        "multiplicity",
        "functional",
        "basis",
        "scf_converged",
        "energy_hartree",
    ],
)
def test_native_metadata_mutations(raw, tmp_path, field):
    key = "grid/angle-angle-0000"
    dest = tmp_path / "job"
    shutil.copytree(raw / "outputs" / key, dest)
    r = read(dest / "result.json")
    r[field] = {
        "units": {},
        "elements": ["H"] * 9,
        "charge": 1,
        "multiplicity": 3,
        "functional": "b3lyp",
        "basis": "def2-svp",
        "scf_converged": False,
        "energy_hartree": 0.0,
    }[field]
    (dest / "result.json").write_text(json.dumps(r))
    with pytest.raises(CompletionError):
        native(dest, read(raw / "requests" / (key + ".json")))


@pytest.mark.parametrize("field", ["graph_identity", "atom_order", "coordinates"])
def test_native_input_mapping_mutation(raw, field):
    key = "grid/angle-angle-0000"
    r = read(raw / "requests" / (key + ".json"))
    if field == "coordinates":
        r[field][0][0] += 0.1
    elif field == "atom_order":
        r[field] = list(reversed(r[field]))
    else:
        r[field] = "0" * 64
    with pytest.raises(CompletionError):
        native(raw / "outputs" / key, r)


def test_every_real_dataset_reconstructed(raw):
    assert normalize(EVIDENCE, raw) == read(EVIDENCE / "normalized-convergence.json")


def test_convergence_and_fit_decision_separated():
    c = read(EVIDENCE / "QM-convergence.json")
    assert len(c["comparisons"]) == 6
    assert all(len(r["geometries"]) >= 3 for r in c["comparisons"])
    assert c["final_QM_adequacy"] == "not_established" and not c["fit_authorized"]
    assert read(EVIDENCE / "residual-identifiability.json")["decision"] == "no_fit"


def test_actual_center_trajectory(raw):
    from island.pcff_qm_pilot.convergence_analysis import center_report

    report = center_report(raw, REPO)
    assert report == read(EVIDENCE / "center-effects.json")
    assert report["optimization_status"] == "completed"
    assert len(report["native_completed_SCF_trajectory"]) == 20
    assert report["unfinished_geometry_count"] == 0
    assert report["new_center_force_max_kcal_mol_A"] < 0.05
    assert report["same_protocol_center_energy_lowering_kcal_mol"] < -0.5
    assert not report["minimum_stability_claimed"]


def test_relocated_fresh_process(tmp_path, distribution):
    moved = tmp_path / "evidence"
    shutil.copytree(EVIDENCE, moved)
    output = tmp_path / "replayed.json"
    env = {
        **os.environ,
        "PYTHONPATH": str(REPO) + ":" + str(REPO / "src"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j59",
            "replay",
            "--repo",
            str(REPO),
            "--evidence",
            str(moved),
            "--distribution",
            str(distribution),
            "--output",
            str(output),
        ],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    assert read(output) == read(EVIDENCE / "relocation-replay.json")


def test_historical_identities_remain_unchanged():
    from scripts.validate_phase_4j47_receipt import history, pin

    baseline = history(REPO, "d361f7e0a7c16a3a3a3513987971db6127690386")
    assert all(pin(REPO / name) == sha for name, sha in baseline.items())


def test_native_incomplete_gradient_even_with_updated_file_hash(raw, tmp_path):
    from island.pcff_completion import pin

    key = "grid/angle-angle-0000"
    dest = tmp_path / "job"
    shutil.copytree(raw / "outputs" / key, dest)
    np.save(
        dest / "gradient.npy",
        np.load(dest / "gradient.npy", allow_pickle=False)[:-1],
        allow_pickle=False,
    )
    r = read(dest / "result.json")
    r["native_files"]["gradient.npy"] = pin(dest / "gradient.npy")
    (dest / "result.json").write_text(json.dumps(r))
    with pytest.raises(CompletionError, match="native_gradient_mutation"):
        native(dest, read(raw / "requests" / (key + ".json")))


def test_source_file_replacement_rejected(tmp_path, distribution):
    from island.pcff_pdms import Class2AuditError, PCFFPDMSAssets

    shutil.copytree(distribution, tmp_path / "source")
    path = tmp_path / "source/pcff.frc"
    path.write_bytes(path.read_bytes() + b"\n! replaced source\n")
    with pytest.raises(Class2AuditError, match="source_hash_mutation"):
        PCFFPDMSAssets(tmp_path / "source").load()


@pytest.mark.parametrize(
    "field", ["coordinates", "graph_identity", "charge", "basis", "scf_type", "threads"]
)
def test_plan_rejects_self_consistent_request_mutation(raw, plan, field):
    from island.pcff_qm_pilot.convergence import validate_request

    work = read(EVIDENCE / "plan/campaign-fast.json")
    # This original v2 request is reconstructed from the exact archived input.
    job = read(raw / "requests/grid/angle-angle-0000.json")
    assert work["jobs"]["grid/angle-angle-0000"]
    if field == "coordinates":
        job[field][0][0] += 0.1
    elif field == "charge":
        job[field] = 1
    elif field == "threads":
        job[field] = 64
    else:
        job[field] = "mutated"
    with pytest.raises(CompletionError):
        validate_request("grid/angle-angle-0000", job, plan)


def test_native_actual_protocol_replacement(raw, tmp_path):
    from island.pcff_completion import pin

    key = "grid/angle-angle-0000"
    dest = tmp_path / "job"
    shutil.copytree(raw / "outputs" / key, dest)
    path = dest / "psi4.out"
    before = path.read_text()
    after = before.replace("Basis Set: DEF2-QZVP", "Basis Set: DEF2-SVP")
    assert after != before
    path.write_text(after)
    result = read(dest / "result.json")
    result["native_files"]["psi4.out"] = pin(path)
    (dest / "result.json").write_text(json.dumps(result))
    with pytest.raises(CompletionError, match="native_input_protocol_mutation"):
        native(dest, read(raw / "requests" / (key + ".json")))


@pytest.mark.parametrize("key", ["grid/angle-angle-0000", "grid/center"])
def test_fresh_execution_dispatches_exact_executed_adapter(raw, tmp_path, key):
    from island.pcff_completion import pin
    from island.pcff_qm_pilot.convergence import execution_adapter

    job = read(raw / "requests" / (key + ".json"))
    result = read(raw / "outputs" / key / "result.json")
    command = execution_adapter(EVIDENCE, tmp_path, sys.executable, job)
    assert command[0] == str(Path(sys.executable).resolve())
    if job.get("threads", 2) == 2:
        assert pin(Path(command[1])) == result["adapter_sha256"]
        assert execution_adapter(EVIDENCE, tmp_path, sys.executable, job) == command
    else:
        assert command[1:] == ["-m", "island.pcff_qm_pilot.convergence_fast_engine"]
        assert (
            pin(REPO / "src/island/pcff_qm_pilot/convergence_fast_engine.py")
            == result["adapter_sha256"]
        )


def test_fresh_execution_rejects_changed_historical_adapter(raw, tmp_path):
    import lzma

    from island.pcff_qm_pilot.convergence import execution_adapter

    moved = tmp_path / "evidence"
    shutil.copytree(EVIDENCE / "engine", moved / "engine")
    shutil.copytree(EVIDENCE / "plan", moved / "plan")
    path = moved / "engine/convergence_df_engine.py.xz"
    path.write_bytes(
        lzma.compress(lzma.decompress(path.read_bytes()) + b"\n# mutation\n")
    )
    job = read(raw / "requests/grid/angle-angle-0000.json")
    with pytest.raises(CompletionError, match="historical_execution_adapter_mutation"):
        execution_adapter(moved, tmp_path / "work", sys.executable, job)
