"""Software gate/reference tests; invented values are not scientific acceptance."""

import importlib
from copy import deepcopy

import pytest


@pytest.fixture
def gate(monkeypatch):
    monkeypatch.syspath_prepend("scripts")
    return importlib.import_module("validate_pcff_amines")


def positive():
    return {
        "case": "primary",
        "typing_complete": True,
        "charges_complete": True,
        "model_complete": True,
        "public_preparation_complete": True,
        "independent_types_and_charges": True,
    }


def negative():
    return {
        "case": "guanidinium",
        "typing_complete": True,
        "charges_complete": False,
        "model_complete": False,
        "public_exception": {
            "type": "PCFFError",
            "message": "Complete typing and charges required",
        },
        "independent_charge": {
            "complete": False,
            "missing": [],
            "components": [
                {
                    "formal": 1,
                    "known_contribution_total": "0.9999",
                    "residual": "-0.0001",
                }
            ],
        },
        "model_diagnostics": [
            {
                "reason": "native_charge_incomplete; parameterization and model construction not attempted"
            }
        ],
    }


def test_positive_and_expected_chemical_failure(gate):
    assert gate.assess(positive()) == []
    assert gate.assess(negative()) == []


@pytest.mark.parametrize(
    "change",
    [
        {"public_preparation_complete": True},
        {
            "execution_error": {
                "stage": "preparation",
                "type": "RuntimeError",
                "message": "backend",
            }
        },
        {"public_exception": {"type": "PCFFError", "message": "wrong source"}},
        {"charges_complete": True},
        {"independent_charge": {"complete": False, "missing": [], "components": []}},
    ],
)
def test_negative_gate_rejects_wrong_success_exception_stage_or_reason(gate, change):
    r = negative()
    r.update(change)
    assert gate.assess(r)


@pytest.mark.parametrize("names", [[], ["primary"], ["primary", "primary"]])
def test_gate_requires_exact_case_set(gate, names):
    assert gate.assess_set([dict(positive(), case=name) for name in names])


def test_aa_reference_rejects_center_moving_reverse(gate, tmp_path):
    ref = importlib.import_module("pcff_j14_reference")
    raw = tmp_path / "synthetic.frc"
    raw.write_text(
        "#angle-angle cff91\n1 1 h+ c n+ c 99\n1 1 c n+ c h+ 2\n1 1 c n+ h+ c 3\n"
    )
    sections = {
        "Angles": [
            ["1", "1", "1", "2", "3"],
            ["2", "2", "1", "2", "4"],
            ["3", "2", "3", "2", "4"],
        ],
        "Angle Coeffs": [["1", "109"], ["2", "110"]],
        "AngleAngle Coeffs": [["1", "99", "99", "3", "109", "110", "110"]],
        "Impropers": [["1", "1", "1", "2", "3", "4"]],
    }
    labels = {101: "c", 208: "n+", 511: "c", 702: "h+"}
    commands, evidence = ref.aa_raw_reference_coefficients(sections, labels, raw)
    assert evidence[0]["corrected"] == [2, 2, 3, 109, 110, 110]
    assert evidence[0]["raw_converter"][:2] == [99, 99]
    assert commands == ["improper_coeff 1 aa 2 2 3 109 110 110"]
    bad = deepcopy(labels)
    bad[208] = "na"
    with pytest.raises(ValueError, match="coverage"):
        ref.aa_raw_reference_coefficients(sections, bad, raw)
