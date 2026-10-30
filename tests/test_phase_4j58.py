"""J58 real-data integrity/replay tests. No synthetic QM evidence or engine launch."""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import CompletionError, digest, inventory, pin, write
from island.pcff_qm.records import attestations
from island.pcff_qm.reference import intake
from island.pcff_qm_pilot.analysis import hessian_projection, spectrum, symmetry_key
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.campaign import execute, selected
from island.pcff_qm_pilot.evidence import replay
from island.pcff_qm_pilot.importer import native

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "docs/evidence/phase_4j58"
DIST = (
    ROOT.parent / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    temp = tmp_path_factory.mktemp("j58-real-data")
    unpack(EVIDENCE / "native-outputs.tar.xz", temp / "raw")
    unpack(EVIDENCE / "reference-datasets.tar.xz", temp / "normalized")
    return temp / "raw", temp / "normalized"


@pytest.fixture(scope="module")
def binding():
    if not (DIST / "pcff.frc").exists():
        pytest.skip("original pinned FRC source unavailable")
    scan = json.loads((EVIDENCE / "plan/angle-angle.json").read_text())
    s, a, ts = selected(ROOT, DIST, scan["center_coordinates"])
    return s, a, ts[0]


def test_verified_ancestry_and_small_bounded_plan():
    ancestry = json.loads((EVIDENCE / "plan/ancestry.json").read_text())
    assert (
        ancestry["origin_main"]
        == ancestry["actual_base"]
        == "f471e06894275b5b6485461d2bc1dcd594f7dfa7"
    )
    assert len(ancestry["prerequisites"]) == 3
    assert all(r["contained_in_main"] for r in ancestry["prerequisites"])
    plan = json.loads((EVIDENCE / "plan/pilot-plan.json").read_text())
    assert plan["budget"]["scan_gradient_jobs_max"] == 50
    assert len(plan["selected_targets"]) == 2
    assert all(len(t["atom_order"]) == 9 for t in plan["selected_targets"])
    assert plan["claims"]["parameter_completeness"] == "incomplete"


def test_new_centers_never_reuse_historical_target_identity():
    for name in ("angle-angle", "bond-bond"):
        scan = json.loads((EVIDENCE / "plan" / (name + ".json")).read_text())
        assert scan["target"]["identity"] != scan["parent_J57_target_identity"]
        assert len(scan["points"]) == 25 and len(scan["held_out_v2"]) == 8
        assert len({p["geometry_checksum"] for p in scan["points"]}) == 25
        assert scan["fit_authorized"] is False


def test_real_execution_and_failed_optimization_preserved(data):
    raw, _ = data
    failed = json.loads((raw / "optimize/result.json").read_text())
    successful = json.loads((raw / "optimize-restart/result.json").read_text())
    assert (
        failed["status"] == "failed"
        and failed["error_type"] == "OptimizationConvergenceError"
    )
    assert successful["status"] == "completed"
    campaign = json.loads((EVIDENCE / "plan/campaign.json").read_text())
    assert campaign["counts"] == {
        "scan_gradients": 50,
        "hessians": 2,
        "sensitivity_gradients": 6,
        "finite_difference_energies": 4,
    }
    for key in campaign["jobs"]:
        request = json.loads((raw / "requests" / (key + ".json")).read_text())
        r, _ = native(raw / "outputs" / key, request)
        assert r["scf_converged"] and r["engine"] == "Psi4"


def test_j57_strict_intake_on_real_dataset(data, binding):
    _, normalized = data
    s, a, target = binding
    request = json.loads((normalized / "angle-angle-request.json").read_text())
    r = intake(
        s,
        a,
        target,
        request,
        normalized / "angle-angle",
        verification_sha256=pin(normalized / "angle-angle/verification.json"),
    )
    assert r.payload == json.loads((normalized / "angle-angle-intake.json").read_text())
    assert (
        r.payload["sample_count"] == 25 and r.payload["status"] == "qm_data_available"
    )
    assert r.payload["parameter_completeness"] == "incomplete"


@pytest.mark.parametrize(
    "field",
    [
        "graph_identity",
        "typing_identity",
        "charge_identity",
        "atom_order",
        "charge_vector",
        "source_hashes",
        "units",
        "functional_form_identity",
        "target_identity",
        "methods",
    ],
)
def test_strict_intake_rejects_resigned_binding_mutation(field, data, binding):
    _, normalized = data
    s, a, target = binding
    request = json.loads((normalized / "angle-angle-request.json").read_text())
    request[field] = "mutated"
    with pytest.raises((CompletionError, AttributeError)):
        intake(
            s,
            a,
            target,
            request,
            normalized / "angle-angle",
            verification_sha256=pin(normalized / "angle-angle/verification.json"),
        )


@pytest.mark.parametrize(
    "field,value",
    [
        ("scf_converged", False),
        ("units", {}),
        ("axis_order", "axis-major"),
        ("elements", ["H"] * 9),
        ("charge", 1),
        ("multiplicity", 3),
        ("basis", "sto-3g"),
        ("functional", "HF"),
        ("geometry_bohr", [[0, 0, 0]] * 9),
        ("energy_hartree", -1.0),
    ],
)
def test_native_import_rejects_metadata_or_mapping_mutation(
    field, value, data, tmp_path
):
    raw, _ = data
    source = raw / "outputs/angle-angle/point-0000"
    dest = tmp_path / "job"
    shutil.copytree(source, dest)
    r = json.loads((dest / "result.json").read_text())
    r[field] = value
    write(dest / "result.json", r)
    request = json.loads((source / "input.json").read_text())
    with pytest.raises(CompletionError):
        native(dest, request)


@pytest.mark.parametrize("mutation", ["nonfinite", "incomplete", "sign", "mapping"])
def test_native_gradient_mutation_rejected_even_with_updated_hash(
    mutation, data, tmp_path
):
    raw, _ = data
    source = raw / "outputs/angle-angle/point-0000"
    dest = tmp_path / "job"
    shutil.copytree(source, dest)
    g = np.load(dest / "gradient.npy", allow_pickle=False)
    if mutation == "nonfinite":
        g[0, 0] = float("nan")
    elif mutation == "incomplete":
        g = g[:-1]
    elif mutation == "sign":
        g = -g
    else:
        g = g[::-1]
    np.save(dest / "gradient.npy", g, allow_pickle=False)
    r = json.loads((dest / "result.json").read_text())
    r["native_files"]["gradient.npy"] = pin(dest / "gradient.npy")
    write(dest / "result.json", r)
    with pytest.raises(CompletionError):
        native(dest, json.loads((source / "input.json").read_text()))


@pytest.mark.parametrize("mutation", ["coverage", "duplicate", "force", "units"])
def test_reference_coverage_and_values_cannot_be_resigned(
    mutation, data, binding, tmp_path
):
    _, normalized = data
    s, a, target = binding
    dest = tmp_path / "dataset"
    shutil.copytree(normalized / "angle-angle", dest)
    request = json.loads((normalized / "angle-angle-request.json").read_text())
    samples = json.loads((dest / "samples.json").read_text())
    if mutation == "coverage":
        samples.pop()
    elif mutation == "duplicate":
        samples[1]["coordinates"] = copy.deepcopy(samples[0]["coordinates"])
    elif mutation == "force":
        samples[0]["forces"] = samples[0]["forces"][:-1]
    else:
        request["energy_factor"] *= 4.184
    write(dest / "samples.json", samples)
    request["files"] = inventory(dest)
    with pytest.raises(CompletionError):
        intake(
            s,
            a,
            target,
            request,
            dest,
            verification_sha256=pin(dest / "verification.json"),
        )


def test_independent_real_unit_and_multistep_fd_checks():
    check = json.loads((EVIDENCE / "independent-unit-checks.json").read_text())
    assert check["production_importer_used"] is False and len(check["samples"]) == 50
    assert len(check["finite_difference"]) == 2
    assert all(r["absolute_error"] <= 0.01 for r in check["finite_difference"])


def test_hessian_rigid_modes_projected_without_full_positive_definiteness(data):
    raw, _ = data
    job = raw / "outputs/angle-angle/hessian"
    r, _ = native(job, json.loads((job / "input.json").read_text()))
    h = np.load(job / "hessian.npy", allow_pickle=False)
    result = hessian_projection(r["geometry_bohr"], h)
    assert result["rigid_mode_rank"] == 6 and result["vibrational_dimension"] == 21
    assert result["strict_full_positive_definiteness_required"] is False


def test_identifiability_does_not_equate_qm_curve_to_one_cross_term():
    r = json.loads((EVIDENCE / "identifiability.json").read_text())
    assert r["created_parameters"] == [] and not r["fit_completed"]
    assert r["other_contributions_not_zero"] and r["affected_interactions"]
    assert r["discrepancy_counterexample"]["nullity"] > 0
    assert all(x["fitted_energy_error"] is None for x in r["partial_residuals"])
    assert all(
        x["unexplained_force_rms_kcal_mol_angstrom"] > 0 for x in r["partial_residuals"]
    )
    assert r["parameter_completeness"] == "incomplete"


def test_source_role_symmetry_and_rank_diagnostics():
    assert symmetry_key("angle-angle:a,b,c,d") == symmetry_key("angle-angle:d,b,c,a")
    assert symmetry_key("bond-bond_1_3:a,b,c,d") == symmetry_key(
        "bond-bond_1_3:d,c,b,a"
    )
    assert spectrum([[1, 2], [2, 4]])["rank"] == 1
    assert spectrum([[1, 2], [2, 4]])["conditioning"] is None


def test_user_attestations_stay_distinct_from_author_artifacts():
    assert len(attestations()) == 28
    assert all(
        r.payload["source_artifact_status"] == "not_yet_supplied"
        for r in attestations()
    )
    assert all(
        r.payload["evidence_kind"] == "user_attested_author_confirmation"
        for r in attestations()
    )


def test_relocation_replay_and_completed_resume_without_engine(data):
    raw, normalized = data
    if not (DIST / "pcff.frc").exists():
        pytest.skip("original pinned FRC source unavailable")
    r = replay(ROOT, DIST, EVIDENCE, raw, normalized)
    assert r["native_completed_jobs"] == 62 and r["fit_completion"] == "not_run"
    rows = execute(ROOT, EVIDENCE, raw, sys.executable)
    assert len(rows) == 62 and all(x["status"] == "verified_resume_skip" for x in rows)


def test_fresh_process_reconstruction_under_different_hash_seeds(data):
    raw, normalized = data
    if not (DIST / "pcff.frc").exists():
        pytest.skip("original pinned FRC source unavailable")
    code = "import json,sys; from island.pcff_qm_pilot.evidence import replay; print(json.dumps(replay(*sys.argv[1:]),sort_keys=True))"
    outputs = []
    for seed in (17, 91):
        p = subprocess.run(
            [
                sys.executable,
                "-c",
                code,
                str(ROOT),
                str(DIST),
                str(EVIDENCE),
                str(raw),
                str(normalized),
            ],
            cwd=raw.parent,
            env={
                **os.environ,
                "PYTHONPATH": str(ROOT) + ":" + str(ROOT / "src"),
                "PYTHONHASHSEED": str(seed),
                "OMP_NUM_THREADS": "1",
            },
            capture_output=True,
            text=True,
            check=True,
        )
        outputs.append(json.loads(p.stdout))
    assert outputs[0] == outputs[1]


def test_archive_path_escape_rejected(tmp_path):
    import io
    import tarfile

    archive = tmp_path / "bad.tar.xz"
    with tarfile.open(archive, "w:xz") as handle:
        item = tarfile.TarInfo("../escape")
        item.size = 1
        handle.addfile(item, io.BytesIO(b"x"))
    with pytest.raises(CompletionError):
        unpack(archive, tmp_path / "unpack")


def test_historical_inputs_are_preserved():
    from scripts.validate_phase_4j47_receipt import history
    from scripts.validate_phase_4j47_receipt import pin as historical_pin

    baseline = history(ROOT, "f471e06894275b5b6485461d2bc1dcd594f7dfa7")
    assert all(historical_pin(ROOT / name) == sha for name, sha in baseline.items())


def test_hessian_metadata_recovery_is_additive_and_cannot_be_resigned(data, tmp_path):
    from island.pcff_qm_pilot.recovery import recovered

    raw, _ = data
    source = raw / "outputs/angle-angle/hessian"
    assert json.loads((source / "result.json").read_text())["status"] == "failed"
    sidecar = json.loads((source / "metadata-recovery.json").read_text())
    assert sidecar["original_failure_preserved"] and not sidecar["QM_rerun"]
    dest = tmp_path / "job"
    shutil.copytree(source, dest)
    sidecar["metadata"]["energy_hartree"] += 0.01
    sidecar["identity"] = digest({k: v for k, v in sidecar.items() if k != "identity"})
    write(dest / "metadata-recovery.json", sidecar)
    with pytest.raises(CompletionError):
        recovered(dest)
