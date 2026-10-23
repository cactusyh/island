"""J38 software negative controls and actual execution; no equilibrium evidence."""

import json
import os
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from island.analysis.npt import FLAGS, _signed
from island.exceptions import IslandError
from scripts.phase_4j38_contract import (
    BASE,
    CASES,
    GATES,
    base_row,
    gate,
    validate_matrix,
    validate_row,
)
from scripts.phase_4j38_fixture import BENCHMARKS, build
from scripts.phase_4j38_observables import unwrap
from scripts.validate_phase_4j38_receipt import validate_case_bundle


def synthetic_row(key):
    """Contract-only record; never publication or scientific acceptance evidence."""
    row = base_row(key)
    row["diagnostic"] = {"purpose": "software schema negative controls"}
    return _signed(row)


@pytest.mark.parametrize("spec", BENCHMARKS, ids=lambda s: s.key)
def test_psmiles_complete_deterministic_benchmark(spec):
    a, evidence = build(spec.key)
    b, second = build(spec.key)
    assert evidence == second and a.packed.graph.identity == b.packed.graph.identity
    assert len(a.chains) == spec.chains and evidence["psmiles"] == spec.psmiles
    assert evidence["setup"]["crosslink_count"] == (3 if spec.key == "network" else 0)
    assert evidence["setup"]["component_count"] == (
        1 if spec.key == "network" else spec.chains
    )
    assert a.packed.system.box.periodic == (True, True, True)
    assert all(s.element != "*" for s in a.packed.system.topology.sites.values())


def test_all_blocked_is_readable_but_never_acceptance():
    matrix = _signed(
        dict(
            schema="island.j38-matrix.v1",
            base_commit=BASE,
            cases={k: synthetic_row(k) for k in CASES},
            acceptance_met=False,
            software_acceptance_met=False,
            **FLAGS,
        )
    )
    assert not validate_matrix(matrix)["software_acceptance_met"]


@pytest.mark.parametrize(
    "mutation",
    [
        "success",
        "origin",
        "family",
        "benchmark",
        "readiness",
        "missing_gate",
        "equilibrium",
        "missing_checkpoint",
        "missing_reconstruction",
        "missing_source",
        "wrong_inventory",
    ],
)
def test_resigned_contract_mutations_fail(mutation):
    row = base_row("oplsaa/pe")
    for g in GATES[:4]:
        row["gates"][g] = gate("passed", "synthetic schema fixture only")
    row.update(
        execution_met=True,
        software_acceptance_met=True,
        charge_origin="OPLS-AA source",
        frame_count=321,
        density_blocks=4,
    )
    for k in (
        "split_equal",
        "checkpoint_byte_equal",
        "rng_equal",
        "relocation",
        "fresh_process",
        "box_volume_changed",
    ):
        row[k] = True
    for k in (
        "final_graph_identity",
        "packing_identity",
        "source_identity",
        "assignment_identity",
        "charge_identity",
        "typed_graph_identity",
        "workflow_identity",
        "trajectory_identity",
        "checkpoint_identity",
        "analysis_identity",
        "comparison_identity",
    ):
        row[k] = "a" * 64
    validate_row("oplsaa/pe", _signed(row))
    if mutation == "success":
        row["acceptance_met"] = True
    elif mutation == "origin":
        row["charge_origin"] = "AM1-BCC"
    elif mutation == "family":
        row["family"] = "gaff2"
    elif mutation == "benchmark":
        row["benchmark"] = "network"
    elif mutation == "readiness":
        row["production_validated"] = True
    elif mutation == "missing_gate":
        row["gates"].pop("numerical_reference")
    elif mutation == "equilibrium":
        row["gates"]["thermodynamic_equilibration"] = gate("passed", "short run")
    elif mutation == "missing_checkpoint":
        row.pop("checkpoint_identity")
    elif mutation == "missing_reconstruction":
        row["fresh_process"] = False
    elif mutation == "missing_source":
        row["gates"]["source_completeness"] = gate("blocked", "no source")
    else:
        row["case"] = "unknown/pe"
    with pytest.raises(ValueError):
        validate_row("oplsaa/pe", _signed(row))


def test_image_translation_and_cycle_controls():
    xyz = np.array([[1.0, 0, 0], [2.0, 0, 0], [3.0, 0, 0]])
    adj = {0: {1, 2}, 1: {0, 2}, 2: {0, 1}}
    expected = unwrap(xyz, np.ones(3) * 10, adj)
    translated = xyz + np.array([[0, 0, 0], [10, -20, 10], [-10, 10, 20]])
    assert np.allclose(unwrap(translated, np.ones(3) * 10, adj), expected)
    with pytest.raises(ValueError, match="cycle"):
        unwrap(np.array([[0.0, 0, 0], [4.0, 0, 0], [8.0, 0, 0]]), np.ones(3) * 10, adj)
    with pytest.raises(ValueError, match="half-box"):
        unwrap(np.array([[0.0, 0, 0], [5.0, 0, 0]]), np.ones(3) * 10, {0: {1}, 1: {0}})


@pytest.fixture(scope="module")
def real(tmp_path_factory):
    root = tmp_path_factory.mktemp("j38-real") / "pe"
    python = Path("../island-validation/phase4g1-env/bin/python")
    env = {
        **os.environ,
        "PYTHONPATH": ".:src",
        "OPENMM_CPU_THREADS": "1",
        "OMP_NUM_THREADS": "1",
    }
    command = [
        str(python),
        "-m",
        "scripts.validate_phase_4j38",
        "--output",
        str(root),
        "--case",
        "oplsaa/pe",
    ]
    child = subprocess.run(
        command, env=env, capture_output=True, text=True, check=False
    )
    assert child.returncode == 1, child.stdout + child.stderr
    row = json.loads((root / "result.json").read_text())
    assert row["execution_met"] and row["software_acceptance_met"], row
    return root, row, python, env


def test_actual_pipeline_does_not_claim_equilibrium(real):
    root, row, python, env = real
    assert (
        row["split_equal"]
        and row["rng_equal"]
        and row["checkpoint_byte_equal"]
        and row["fresh_process"]
    )
    assert not row["acceptance_met"] and row["density_blocks"] == 4
    assert json.loads((root / "physics.json").read_text())["passed"]
    # Relocate whole case, reconstruct observables and all workflow bindings in a child.
    moved = root.parent / "moved"
    shutil.copytree(root, moved)
    code = 'import json,sys; from scripts.validate_phase_4j38_receipt import validate_case_bundle; from pathlib import Path; p=Path(sys.argv[1]); validate_case_bundle(p,json.loads((p/"result.json").read_text()))'
    child = subprocess.run(
        [str(python), "-c", code, str(moved)],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stdout + child.stderr


@pytest.mark.parametrize(
    "field",
    [
        "density",
        "mapping",
        "source",
        "trajectory",
        "configuration",
        "checkpoint",
        "observables",
    ],
)
def test_actual_bundle_mutations_fail(real, tmp_path, field):
    root, row, _, _ = real
    copy = tmp_path / "copy"
    shutil.copytree(root, copy)
    if field == "density":
        payload = json.loads((copy / "observables.json").read_text())
        payload.pop("identity")
        payload["frames"][0]["density_g_cm3"] += 0.01
        payload = _signed(payload)
        (copy / "observables.json").write_text(json.dumps(payload))
        row = deepcopy(row)
        row["analysis_identity"] = payload["identity"]
    elif field == "mapping":
        payload = json.loads((copy / "fixture.json").read_text())
        payload["endpoints"][0].reverse()
        (copy / "fixture.json").write_text(json.dumps(payload))
    elif field == "source":
        row = deepcopy(row)
        row["source_identity"] = "0" * 64
    elif field == "trajectory":
        (copy / "trajectory.npz").write_bytes(b"not an MD trajectory")
    elif field == "configuration":
        payload = json.loads((copy / "uninterrupted/manifest.json").read_text())
        payload["config"]["pressure_bar"] = 2
        (copy / "uninterrupted/manifest.json").write_text(json.dumps(payload))
    elif field == "checkpoint":
        payload = json.loads((copy / "resumed/frames.json").read_text())
        payload["frames"][-1]["checkpoint_base64"] = "AAAA"
        (copy / "resumed/frames.json").write_text(json.dumps(payload))
    else:
        payload = json.loads((copy / "observables.json").read_text())
        payload.pop("identity")
        payload["frames"][0]["chain_rg_angstrom"][0] += 1
        payload = _signed(payload)
        (copy / "observables.json").write_text(json.dumps(payload))
        row = deepcopy(row)
        row["analysis_identity"] = payload["identity"]
    with pytest.raises((ValueError, OSError, KeyError, IslandError)):
        validate_case_bundle(copy, row)


def test_missing_new_graph_resp_stays_nonzero(tmp_path):
    python = "../island-validation/phase4g1-env/bin/python"
    root = tmp_path / "unverified-resp"
    child = subprocess.run(
        [
            python,
            "-m",
            "scripts.validate_phase_4j38",
            "--output",
            str(root),
            "--case",
            "gaff2/pe",
        ],
        env={**os.environ, "PYTHONPATH": ".:src"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 1 and (root / "result.json").is_file()
    row = json.loads((root / "result.json").read_text())
    assert (
        not row["execution_met"] and row["diagnostic"]["stage"] == "charge_provenance"
    )
    assert row["gates"]["source_completeness"]["status"] == "blocked"


@pytest.mark.parametrize(
    "mutation",
    [
        "all_blocked_success",
        "missing_case",
        "key_identity",
        "receipt_manifest",
        "test_counts",
        "readiness",
        "pcff_success",
        "am1bcc_success",
        "specification",
    ],
)
def test_consistently_resigned_receipt_manifest_mutations(mutation):
    from dataclasses import asdict

    from scripts.validate_phase_4j38_receipt import validate_structure

    matrix = dict(
        schema="island.j38-matrix.v1",
        base_commit=BASE,
        cases={k: synthetic_row(k) for k in CASES},
        acceptance_met=False,
        software_acceptance_met=False,
        **FLAGS,
    )
    matrix["pcff_control"] = {
        "acceptance_met": True,
        "dynamics_executed": False,
        "diagnostic": {"classification": "expected_typing_rejection"},
    }
    matrix["native_am1bcc"] = {"status": "blocked", "acceptance_met": False}
    manifest = {
        "schema": "island.j38-reproducibility.v1",
        "base_commit": BASE,
        "benchmarks": [asdict(b) for b in BENCHMARKS],
    }
    receipt = dict(
        schema="island.j38-receipt.v1",
        base_commit=BASE,
        cases=deepcopy(matrix["cases"]),
        acceptance_met=False,
        software_acceptance_met=False,
        tests={
            k: {"passed": 1, "failed": 0, "errors": 0}
            for k in ("focused", "j22-j37", "full_supported")
        },
        **FLAGS,
    )
    if mutation == "all_blocked_success":
        receipt["acceptance_met"] = matrix["acceptance_met"] = True
    elif mutation == "missing_case":
        matrix["cases"].pop("gaff2/network")
        receipt["cases"] = deepcopy(matrix["cases"])
    elif mutation == "key_identity":
        row = dict(matrix["cases"]["oplsaa/pe"])
        row.pop("identity")
        row["benchmark"] = "peo"
        matrix["cases"]["oplsaa/pe"] = _signed(row)
        receipt["cases"] = deepcopy(matrix["cases"])
    elif mutation == "test_counts":
        receipt["tests"]["focused"]["failed"] = 1
    elif mutation == "readiness":
        receipt["production_validated"] = True
    elif mutation == "pcff_success":
        matrix["pcff_control"]["dynamics_executed"] = True
    elif mutation == "am1bcc_success":
        matrix["native_am1bcc"]["acceptance_met"] = True
    elif mutation == "specification":
        manifest["benchmarks"][0]["dp"] = 500
    matrix = _signed(matrix)
    manifest["matrix_identity"] = matrix["identity"]
    manifest = _signed(manifest)
    receipt.update(
        matrix_identity=matrix["identity"], manifest_identity=manifest["identity"]
    )
    if mutation == "receipt_manifest":
        receipt["manifest_identity"] = "0" * 64
    receipt = _signed(receipt)
    with pytest.raises(ValueError):
        validate_structure(receipt, manifest, matrix)
