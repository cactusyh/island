"""Bounded software controls, not numerical or real-source acceptance."""

import importlib
from copy import deepcopy

import pytest


@pytest.fixture
def gate(monkeypatch):
    monkeypatch.syspath_prepend("scripts")
    return importlib.import_module("validate_pcff_urethanes")


def declared_records(gate):
    expected = gate.contract()["cases"]
    records = []
    for name, e in expected.items():
        if e["role"] == "negative_unchanged_v5":
            r = deepcopy(e["expected"])
        else:
            r = {
                "typing_complete": True,
                "charges_complete": True,
                "model_complete": e["model_complete"],
                "labels": e["expected_motif_counts"].copy(),
                "model_diagnostics": deepcopy(e["model_diagnostics"]),
                "public_exception": deepcopy(e["public_exception"]),
                "public_preparation_complete": e["model_complete"],
                "independent_types_and_charges": True,
            }
        r["case"] = name
        records.append(r)
    return records, expected


def test_declared_controls_and_reordering(gate):
    r, e = declared_records(gate)
    assert gate.assess_set(r, e) == []
    assert gate.assess_set(r[::-1], e) == []


@pytest.mark.parametrize(
    "change",
    [
        {"public_preparation_complete": True},
        {"execution_error": {"type": "RuntimeError", "message": "unrelated"}},
        {"public_exception": {"type": "PCFFError", "message": "wrong source"}},
        {"charges_complete": True},
        {"independent_charge": {}},
    ],
)
def test_negative_stage_reason_and_failure(gate, change):
    r, e = declared_records(gate)
    c = next(c for c in r if c["case"] == "guanidinium")
    c.update(change)
    assert gate.assess_set(r, e)


@pytest.mark.parametrize(
    "change",
    [
        {"public_preparation_complete": True},
        {"model_complete": True},
        {"model_diagnostics": []},
        {"independent_types_and_charges": False},
        {"labels": {}},
    ],
)
def test_expected_missing_polymer_couplings_are_not_success(gate, change):
    r, e = declared_records(gate)
    c = next(c for c in r if c["case"] == "urethane_DP3")
    c.update(change)
    assert gate.assess_set(r, e)


def test_required_set_cannot_be_empty_missing_duplicate_or_unknown(gate):
    r, e = declared_records(gate)
    for bad in ([], r[:-1], r + [r[0]], r + [{"case": "unknown"}]):
        assert gate.assess_set(bad, e)


def test_contract_tamper_rejected(gate, monkeypatch, tmp_path):
    p = tmp_path / "contract.json"
    p.write_text("{}")
    monkeypatch.setattr(gate, "CONTRACT", p)
    with pytest.raises(ValueError, match="contract changed"):
        gate.contract()
