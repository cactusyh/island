"""Software acceptance controls using retained diagnostics, no chemistry execution."""

import importlib
import json
import sys
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from island.exceptions import PCFFError
from island.forcefields import PreparedForceFieldError
from island.workflows import storage


@pytest.fixture
def harness(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend("scripts")
    driver = importlib.import_module("inspect_pcff_tranche_controls")
    retained = json.loads(Path("docs/evidence/phase_4j13_gaps.json").read_text())
    observations = {
        c["case"]: deepcopy(c)
        for c in retained["controls"] + retained["closed_shell_charged_control"]
    }
    declarations = {
        name: Path("docs/evidence/" + file)
        for name, file in [
            ("main", "phase_4j13_declaration.json"),
            ("charged", "phase_4j13_charge_declaration.json"),
        ]
    }
    graphs = {}
    for path in declarations.values():
        d = json.loads(path.read_text())
        graphs.update(
            {c.get("smiles", c.get("psmiles")): c["name"] for c in d["control_graphs"]}
        )
    source = SimpleNamespace(
        identity={"sha256": d["source"]["sha256"]}, raw=b"SYNTHETIC DRIVER STUB"
    )
    monkeypatch.setattr(driver, "load_pcff_source", lambda path: source)
    monkeypatch.setattr(driver, "read_source", lambda raw: None)

    def system(smiles, **kwargs):
        return SimpleNamespace(
            name=graphs[smiles], to_dict=lambda: {"software_fixture": graphs[smiles]}
        )

    monkeypatch.setattr(driver, "build_linear_polymer", system)
    monkeypatch.setattr(driver, "from_smiles", system)

    def typing(system, source, **kwargs):
        r = observations[system.name]
        return SimpleNamespace(
            complete=r["typing_complete"],
            assignments={},
            payload={"assignments": {}, "diagnostics": r["typing_diagnostics"]},
        )

    monkeypatch.setattr(driver, "type_pcff_atoms", typing)
    monkeypatch.setattr(
        driver,
        "inspect_pcff_operational_support",
        lambda s, *a, **k: observations[s.name],
    )
    monkeypatch.setattr(
        driver,
        "charge_raw",
        lambda raw, s, t: observations[s.name]["independent_charge"],
    )

    def rejected(system, *args):
        kind, message = observations[system.name]["public_rejection"].split(": ", 1)
        raise {
            "PCFFError": PCFFError,
            "PreparedForceFieldError": PreparedForceFieldError,
        }[kind](message)

    monkeypatch.setattr(driver, "prepare_forcefield", rejected)

    def run(*, set_name="main", defect=None):
        declaration = declarations[set_name]
        if defect:
            d = json.loads(declaration.read_text())
            if defect == "empty":
                d["control_graphs"] = []
            elif defect == "missing":
                d["control_graphs"] = d["control_graphs"][:-1]
            elif defect == "duplicate":
                d["control_graphs"].append(deepcopy(d["control_graphs"][0]))
            elif defect == "wrong-set":
                d["control_graphs"] = json.loads(declarations["charged"].read_text())[
                    "control_graphs"
                ]
            declaration = tmp_path / "changed.json"
            declaration.write_text(json.dumps(d))
        output = tmp_path / "out"
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "controls",
                "--source",
                "synthetic.frc",
                "--declaration",
                str(declaration),
                "--output",
                str(output),
            ],
        )
        code = driver.main() or 0
        return code, storage.decode(storage.read_json(output / "report.json"))

    return driver, run, observations


@pytest.mark.parametrize("defect", ["success", "runtime", "wrong-reason"])
def test_rejects_false_public_rejections(harness, monkeypatch, defect):
    driver, run, _ = harness

    def bad(*args):
        if defect == "success":
            return object()
        if defect == "runtime":
            raise RuntimeError("backend unavailable")
        raise PCFFError("wrong source hash")

    monkeypatch.setattr(driver, "prepare_forcefield", bad)
    code, report = run()
    assert code != 0
    assert report["acceptance_passed"] is False


@pytest.mark.parametrize("defect", ["empty", "missing", "duplicate", "wrong-set"])
def test_required_controls_fail_closed(harness, defect):
    _, run, _ = harness
    code, report = run(defect=defect)
    assert code != 0
    assert report["acceptance_passed"] is False


@pytest.mark.parametrize("set_name", ["main", "charged"])
def test_expected_chemical_rejection_passes(harness, set_name):
    _, run, _ = harness
    code, _ = run(set_name=set_name)
    assert code == 0


@pytest.mark.parametrize(
    "defect",
    [
        "typing-reason",
        "charge-stage",
        "model-stage",
        "public-model-reason",
        "increment-path",
    ],
)
def test_correct_exception_does_not_hide_wrong_native_stage(harness, defect):
    _, run, observations = harness
    set_name = "main"
    if defect == "typing-reason":
        observations["protonated_amine"]["typing_diagnostics"][0]["reason"] = (
            "missing_source"
        )
    elif defect == "charge-stage":
        observations["peroxide"]["charges_complete"] = False
        observations["peroxide"]["independent_charge"]["complete"] = False
    elif defect == "model-stage":
        observations["peroxide"]["model_complete"] = True
    elif defect == "public-model-reason":
        observations["peroxide"]["public_rejection"] = (
            "PreparedForceFieldError: Incomplete PCFF model: unrelated failure"
        )
    else:
        set_name = "charged"
        observations["protonated_sidechain"]["independent_charge"]["missing"][0][
            "searches"
        ][0]["types"] = ["na", "hn2"]
    code, report = run(set_name=set_name)
    assert code == 1 and report["acceptance_passed"] is False
    assert any(c["acceptance_failures"] for c in report["cases"])


def test_inspection_runtime_error_is_not_a_chemical_rejection(harness, monkeypatch):
    driver, run, _ = harness

    def failed(*args, **kwargs):
        raise RuntimeError("broken inspection")

    monkeypatch.setattr(driver, "inspect_pcff_operational_support", failed)
    code, report = run()
    assert code == 1
    peroxide = next(c for c in report["cases"] if c["case"] == "peroxide")
    assert peroxide["execution_error"]["stage"] == "charge_model_inspection"
    assert "public_rejection" not in peroxide


@pytest.mark.parametrize("name", ["main", "charged-sidechain"])
@pytest.mark.parametrize("defect", ["empty", "duplicate", "missing"])
def test_rebinding_contract_cannot_drop_required_cases(
    monkeypatch, tmp_path, name, defect
):
    monkeypatch.syspath_prepend("scripts")
    gates = importlib.import_module("pcff_tranche_expectations")
    contract = json.loads(gates.CONTRACT_PATH.read_text())
    selected = contract["control_sets"][name]
    d = json.loads(Path(selected["declaration"]).read_text())
    if defect == "duplicate":
        d["control_graphs"].append(deepcopy(d["control_graphs"][0]))
    elif defect == "empty":
        d["control_graphs"] = []
    else:
        d["control_graphs"] = d["control_graphs"][:-1]
    raw = storage.json_bytes(d)
    selected["declaration_sha256"] = storage.checksum(raw)
    path = tmp_path / "expectations.json"
    path.write_text(json.dumps(contract))
    monkeypatch.setattr(gates, "CONTRACT_PATH", path)
    with pytest.raises(ValueError, match="required controls"):
        gates.contract_for(d, raw, {"sha256": contract["source_sha256"]})


def test_cli_failure_exit_preserves_report_without_scientific_execution(tmp_path):
    import subprocess

    output = tmp_path / "run"
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/inspect_pcff_tranche_controls.py",
            "--source",
            str(tmp_path / "absent.frc"),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert completed.returncode == 1
    report = storage.decode(storage.read_json(output / "report.json"))
    assert report["acceptance_passed"] is False
    assert report["failures"] and report["cases"] == []
