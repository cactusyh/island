"""Synthetic evidence tests; no QM or claims of experimental accuracy."""

import subprocess
import sys
from copy import deepcopy

import numpy as np
import pytest
from test_charge_conservation import synthetic_observation
from test_oligomer_charge_references import molecule

from island.charge_references import (
    PEChargePrediction,
    PETemplateModel,
    fit_pe_template,
    load_pe_record,
    pe_target_correspondence,
    predict_pe_charges,
    project_charge_observation,
    save_record,
    validate_pe_templates,
)
from island.charge_references.pe_template import TRAINING, _solve
from island.charge_references.projection import POLICY
from island.charge_references.records import pack
from island.exceptions import ChargeReferenceError
from island.workflows.bundle import system_data


@pytest.fixture(scope="module")
def models():
    training = [
        project_charge_observation(synthetic_observation(dp=d, seed=s), policy=POLICY)
        for d, s in TRAINING
    ]
    return [fit_pe_template(training, mode=m) for m in ("baseline", "conserving")]


def test_independent_constrained_solution_and_rank():
    x, info = _solve(np.eye(2), [1.0, 3.0], [[1.0, 1.0]])
    assert x == pytest.approx([-1.0, 1.0], abs=1e-14)
    assert info["objective_e2"] == pytest.approx(8.0)
    with pytest.raises(ChargeReferenceError, match="Underdetermined"):
        _solve([[1.0, 0.0]], [1.0], [])
    with pytest.raises(ChargeReferenceError, match="Dependent"):
        _solve(np.eye(2), [1.0, 3.0], [[1.0, 1.0], [2.0, 2.0]])


@pytest.mark.parametrize("dp", [3, 4, 5, 20, 50, 100])
def test_conservation_order_and_ownership(models, dp):
    system = molecule(dp)
    before = system_data(system)
    result = predict_pe_charges(models[1], system)
    q = result.payload["prediction"]
    assert set(q["assignments_e"]) == set(system.topology.sites)
    assert abs(q["total_e"]) < 1e-12
    assert (
        result.payload["prediction"]
        == predict_pe_charges(models[1], molecule(dp, reverse=True)).payload[
            "prediction"
        ]
    )
    system.coordinates.translate((10.0, 2.0, 5.0))
    shifted = predict_pe_charges(models[1], system).payload["prediction"]
    assert shifted == q
    assert before["metadata"] == system.metadata
    q["assignments_e"].clear()
    assert result.payload["prediction"]["assignments_e"]


def test_reversal_and_hydrogen_weights(models):
    s = molecule(5)
    original = predict_pe_charges(models[0], s).payload["prediction"]
    for atom in s.topology.sites.values():
        atom.metadata["repeat_unit_index"] = 4 - atom.metadata["repeat_unit_index"]
        if atom.element == "C":
            atom.metadata["source_repeat_atom_index"] = (
                3 - atom.metadata["source_repeat_atom_index"]
            )
    p = s.metadata["polymer"]
    p["head_site_id"], p["tail_site_id"] = p["tail_site_id"], p["head_site_id"]
    assert (
        predict_pe_charges(models[0], s).payload["prediction"]["assignments_e"]
        == original["assignments_e"]
    )
    fit = models[0].payload["fit"]
    # Independent weighted means over each class, individual H multiplicity.
    sums = np.zeros(6)
    weights = np.zeros(6)
    for target in models[0].payload["training_projections"]:
        source = target["observation"]["data"]["system"]
        from island.workflows.bundle import system_from

        pred = predict_pe_charges(models[0], system_from(source)).payload["prediction"]
        q = target["projection"]["projected_charges"]
        for site, cls in pred["class_indices"].items():
            sums[cls] += q[site] / len(q) / 4
            weights[cls] += 1 / len(q) / 4
    assert fit["coefficients_e"] == pytest.approx(sums / weights, abs=1e-14)


@pytest.mark.parametrize(
    "change", ["peo", "end", "mass", "stereo", "missing_h", "definition", "oversize"]
)
def test_scope(change):
    s = molecule(oxygen=change == "peo", dp=167 if change == "oversize" else 4)
    if change == "end":
        s.metadata["polymer"]["head_site_id"] = 999
    if change == "mass":
        s.topology.sites[118].mass = 13.0
    if change == "stereo":
        s.topology.sites[118].metadata["cip_label"] = "R"
    if change == "definition":
        s.metadata["polymer"]["source_psmiles"] = "[*]CC[*]"
    if change == "missing_h":
        h = next(i for i, a in s.topology.sites.items() if a.element == "H")
        del s.topology.sites[h]
    with pytest.raises(ChargeReferenceError):
        pe_target_correspondence(s)


def test_tamper_and_exclusive_load(models, tmp_path):
    p = models[1].payload
    p["fit"]["coefficients_e"][0] += 0.01
    with pytest.raises(ChargeReferenceError):
        PETemplateModel(pack(p)).validate_integrity()
    pred = predict_pe_charges(models[1], molecule(4))
    p = pred.payload
    p["prediction"]["assignments_e"][118] += 0.1
    with pytest.raises(ChargeReferenceError):
        PEChargePrediction(pack(p)).validate_integrity()
    save_record(pred, tmp_path / "prediction.json")
    assert load_pe_record(tmp_path / "prediction.json").identity == pred.identity
    with pytest.raises(ChargeReferenceError):
        save_record(pred, tmp_path / "prediction.json")


def test_frozen_split_validation(models):
    target = project_charge_observation(
        synthetic_observation(dp=5, seed=314159), policy=POLICY
    )
    identities = [m.identity for m in models]
    result = validate_pe_templates(models, [target], frozen_model_identities=identities)
    assert not result.payload["validation"]["complete_held_out_targets"]
    assert [m.identity for m in models] == identities
    with pytest.raises(ChargeReferenceError):
        validate_pe_templates(models, [target], frozen_model_identities=["a" * 64] * 2)
    training = deepcopy(models[0].payload["training_projections"])
    training[0] = target.payload
    from island.charge_references import ChargeProjection

    with pytest.raises(ChargeReferenceError, match="split"):
        fit_pe_template([ChargeProjection(pack(p)) for p in training], mode="baseline")


def test_offline_optional_isolation(models, tmp_path):
    save_record(models[1], tmp_path / "m.json")
    script = """
import sys
for name in ('rdkit', 'parmed', 'openmm', 'scipy'): sys.modules[name] = None
from island.charge_references import load_pe_record
load_pe_record(sys.argv[1]).validate_integrity()
"""
    subprocess.run([sys.executable, "-c", script, str(tmp_path / "m.json")], check=True)


def test_prediction_integrity_source_and_model_classes(models):
    model = models[1].payload
    model["fit"]["hydrogen_multiplicities"][0] = 2
    with pytest.raises(ChargeReferenceError):
        PETemplateModel(pack(model)).validate_integrity()
    pred = predict_pe_charges(models[1], molecule(4)).payload
    pred["source_system"]["sites"][0]["mass"] = 13.0
    with pytest.raises(ChargeReferenceError):
        PEChargePrediction(pack(pred)).validate_integrity()


def test_cli_failure_and_explicit_gate(monkeypatch, tmp_path):
    from scripts import validate_pe_charge_transfer as cli

    # Injected orchestration outcomes only, never counted as live calculations.
    monkeypatch.setattr(cli, "execute", lambda *args: False)
    assert cli.main(["--execute", "--output", str(tmp_path)]) == 1
    monkeypatch.setattr(cli, "execute", lambda *args: True)
    assert cli.main(["--execute", "--output", str(tmp_path)]) == 0

    def failure(*args):
        raise ChargeReferenceError("injected unavailable backend")

    monkeypatch.setattr(cli, "execute", failure)
    assert cli.main(["--execute", "--output", str(tmp_path)]) == 1
    assert "injected unavailable" in (tmp_path / "failure.json").read_text()


def test_missing_training_never_publishes(tmp_path):
    from scripts import validate_pe_charge_transfer as cli

    destination = tmp_path / "new"
    assert cli.main(["--freeze", "--output", str(destination)]) == 1
    assert not destination.exists()
