"""Fail-closed corrective protocol, mutation and real retained-run controls."""

import copy
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

from island.analysis.equilibration_protocol import (
    CorrectiveProtocol,
    assess,
    prerequisites,
    run_until_stable,
    validate_assessment,
)
from island.analysis.npt_continuation import load_continuation
from island.exceptions import AnalysisError


@pytest.mark.parametrize(
    "change",
    [
        {"replicas": 2},
        {"frame_interval": 4},
        {"target_temperature_kelvin": 100},
        {"maximum_windows": 0},
        {"burn_in": "adaptive_best_fit"},
    ],
)
def test_versioned_protocol_immutable(change):
    with pytest.raises(AnalysisError):
        CorrectiveProtocol(**change)


def test_no_dynamics_before_reference_and_pressure():
    calls = []
    with pytest.raises(AnalysisError, match="pressure_unavailable"):
        run_until_stable(
            lambda *args: calls.append(args),
            reference={"passed": False},
            finite_difference={"passed": False},
            pressure={"target_pressure_bar": 1},
        )
    assert calls == []


@pytest.mark.parametrize(
    "pressure",
    [
        None,
        {"status": "ok"},
        {"target_pressure_bar": 1},
        {"status": "measured_and_independently_validated"},
    ],
)
def test_target_pressure_not_measurement(pressure):
    assert "pressure_unavailable" in prerequisites(
        {"passed": True}, {"passed": True}, pressure
    )


@pytest.fixture(scope="module")
def historical(tmp_path_factory):
    if os.environ.get("ISLAND_J40_RUN"):
        return Path(os.environ["ISLAND_J40_RUN"])
    root = tmp_path_factory.mktemp("j40-c1-history")
    with tarfile.open("docs/evidence/phase_4j40/artifacts.tar.xz") as archive:
        archive.extractall(root, filter="data")
    return root / "run"


@pytest.fixture(scope="module")
def assessment_inputs(historical):
    case = historical / "opls-network"
    workflows = [
        load_continuation(case / f"replica-{i}/uninterrupted")[0] for i in range(2)
    ]
    return workflows, {
        "reference": json.loads((case / "cross-engine.json").read_text()),
        "finite_difference": {"passed": False},
        "pressure": {"status": "pressure_unavailable"},
    }


def test_all_five_failures_reproduced(assessment_inputs):
    workflows, inputs = assessment_inputs
    record = assess(workflows, **inputs)
    failures = record["failed_stages"]
    assert "pressure_unavailable" in failures
    assert "cross_engine_discrepancy" in failures
    assert "finite_difference_interval_not_converged" in failures
    assert any("temperature_target_unmet" in x for x in failures)
    assert any("insufficient_sampling" in x for x in failures)
    assert record["scientific_acceptance"] is False


@pytest.mark.parametrize(
    "mutation",
    [
        "pressure",
        "temperature",
        "source",
        "force",
        "finite_difference",
        "acceptance",
        "protocol",
    ],
)
def test_bound_assessment_mutations(assessment_inputs, mutation):
    workflows, inputs = assessment_inputs
    record = assess(workflows, **inputs)
    changed = copy.deepcopy(record)
    if mutation == "pressure":
        changed["pressure"] = {"target_pressure_bar": 1}
    elif mutation == "temperature":
        changed["analysis"]["replicas"][0]["metrics"]["temperature"]["mean"] = 300
    elif mutation == "source":
        changed["reference_identity"] = "0" * 64
    elif mutation == "force":
        changed["failed_stages"] = []
    elif mutation == "finite_difference":
        changed["finite_difference_identity"] = "1" * 64
    elif mutation == "acceptance":
        changed["scientific_acceptance"] = True
    else:
        changed["protocol"]["statistical_policy"]["maximum_temperature_error"] = 1
    with pytest.raises(AnalysisError):
        validate_assessment(changed, workflows, **inputs)


def test_persisted_diagnostics_fresh_interpreter(
    historical, tmp_path, assessment_inputs
):
    workflows, inputs = assessment_inputs
    expected = assess(workflows, **inputs)
    target = tmp_path / "assessment.json"
    target.write_text(json.dumps(expected))
    code = """import json,sys
from pathlib import Path
from island.analysis.equilibration_protocol import validate_assessment
from island.analysis.npt_continuation import load_continuation
p=Path(sys.argv[1]);w=[load_continuation(p/f"replica-{i}/uninterrupted")[0] for i in range(2)]
validate_assessment(json.load(open(sys.argv[2])),w,reference=json.load(open(p/"cross-engine.json")),finite_difference={"passed":False},pressure={"status":"pressure_unavailable"})
"""
    proc = subprocess.run(
        [sys.executable, "-c", code, str(historical / "opls-network"), str(target)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr


def test_force_audit_requires_exact_atom_mapping():
    from scripts.phase_4j40_c1_finite_difference import audit

    with pytest.raises(AnalysisError, match="ordering"):
        audit(
            {"numerical_data": {"sites": [{"id": 1}, {"id": 2}]}},
            [[0, 0, 0], [1, 1, 1]],
            [2, 1],
        )


@pytest.mark.parametrize(
    "mutation",
    [
        "mass",
        "charge",
        "box",
        "cutoff",
        "mixing",
        "scaling",
        "ewald",
        "coordinates",
        "atom_order",
    ],
)
def test_independent_input_parity_mutations(historical, tmp_path, mutation):
    from scripts.phase_4j38_fixture import build
    from scripts.phase_4j40_c1_parity import audit

    fixture, _ = build("network")
    source = historical / "opls-network"
    payload = json.loads((source / "backend.json").read_text())
    data = (source / "cross-engine/system.data").read_text()
    script = (source / "cross-engine/in.periodic").read_text()
    if mutation == "mass":
        data = data.replace("1 15.999000000000001", "1 16.999000000000001")
    elif mutation == "charge":
        lines = data.splitlines()
        i = next(i for i, line in enumerate(lines) if line.startswith("Atoms")) + 2
        row = lines[i].split()
        row[3] = str(float(row[3]) + 0.01)
        lines[i] = " ".join(row)
        data = "\n".join(lines) + "\n"
    elif mutation == "box":
        data = data.replace("xlo xhi", "xlo zhi")
    elif mutation == "cutoff":
        script = script.replace("lj/cut/coul/long 5", "lj/cut/coul/long 6")
    elif mutation == "mixing":
        script = script.replace("mix geometric", "mix arithmetic")
    elif mutation == "scaling":
        script = script.replace("coul 0 0 0.5", "coul 0 0 1")
    elif mutation == "ewald":
        script = script.replace("kspace_style ewald", "kspace_style pppm")
    elif mutation == "coordinates":
        lines = data.splitlines()
        i = next(i for i, line in enumerate(lines) if line.startswith("Atoms")) + 2
        row = lines[i].split()
        row[4] = str(float(row[4]) + 0.1)
        lines[i] = " ".join(row)
        data = "\n".join(lines) + "\n"
    else:
        lines = data.splitlines()
        i = next(i for i, line in enumerate(lines) if line.startswith("Atoms")) + 2
        lines[i], lines[i + 1] = lines[i + 1], lines[i]
        data = "\n".join(lines) + "\n"
    (tmp_path / "system.data").write_text(data)
    (tmp_path / "in.periodic").write_text(script)
    with pytest.raises(AnalysisError):
        audit(
            payload,
            {
                i: fixture.packed.system.coordinates.get(i)
                for i in sorted(fixture.packed.system.topology.sites)
            },
            tmp_path,
        )


@pytest.fixture(scope="module")
def corrective(tmp_path_factory):
    root = tmp_path_factory.mktemp("c1-corrective")
    with tarfile.open("docs/evidence/phase_4j40_c1/artifacts.tar.xz") as archive:
        archive.extractall(root, filter="data")
    return root / "run"


def test_corrective_bundle_offline(corrective, historical):
    from scripts.validate_phase_4j40_c1_receipt import validate_run

    assert validate_run(corrective, historical)["scientific_acceptance"] is False


@pytest.mark.parametrize(
    "mutation",
    ["graph", "origin", "case", "acceptance", "geometry", "atom_order", "parity"],
)
def test_resigned_corrective_mutation_rejected(
    corrective, historical, tmp_path, mutation
):
    import shutil

    from island.analysis.npt import _body, _signed
    from scripts.validate_phase_4j40_c1_receipt import validate_run

    shutil.copytree(corrective, tmp_path / "run")
    root = tmp_path / "run"
    case = root / "opls-network"
    matrix = json.loads((root / "matrix.json").read_text())
    if mutation in ("graph", "origin", "case"):
        cell = json.loads((case / "result.json").read_text())
        field = {"graph": "graph_identity", "origin": "charge_origin", "case": "case"}[
            mutation
        ]
        cell[field] = "incorrect"
        cell = _signed(_body(cell))
        (case / "result.json").write_text(json.dumps(cell))
        matrix["cases"]["opls-network"] = cell
    elif mutation == "acceptance":
        matrix["scientific_acceptance"] = True
    else:
        target = case / (
            "input-parity.json"
            if mutation == "parity"
            else "packed-finite-difference.json"
        )
        record = json.loads(target.read_text())
        if mutation == "geometry":
            record["coordinates_nm"][0][0] += 0.01
        elif mutation == "atom_order":
            record["atom_ids"].reverse()
        else:
            record["mass"] = "unverified"
        target.write_text(json.dumps(_signed(_body(record))))
    (root / "matrix.json").write_text(json.dumps(_signed(_body(matrix))))
    with pytest.raises((AnalysisError, ValueError)):
        validate_run(root, historical)
