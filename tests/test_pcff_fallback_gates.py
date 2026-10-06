"""Software-only acceptance controls using labelled retained numerical receipts.

No negative case launches scientific execution. Artifact reconstruction is tested
separately with real retained bundles; these tests isolate aggregate decisions.
"""

import json
import sys
from copy import deepcopy
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pcff_fallback_gates as gates


@pytest.fixture
def report(monkeypatch):
    evidence = json.loads(Path("docs/evidence/phase_4j2.json").read_text())
    case = next(x for x in evidence["cases"] if x["case"] == "dichlorine")
    # Explicit software stub, not a native validation or new LAMMPS result.
    monkeypatch.setattr(
        gates,
        "verify_artifacts",
        lambda *args: {
            "parameter_fingerprint": case["model_identity"],
            "model_fingerprint": "f" * 64,
        },
    )
    return {
        "cases": [deepcopy(case)],
        "term_checks": deepcopy(evidence["term_checks"]),
        "source_unchanged": True,
        "inventory": evidence["source"],
    }


def check(report):
    return gates.assignment_and_term_gate(report, Path("unused"), Path("unused.frc"))


def test_valid_retained_software_control(report):
    assert check(report)["passed"]


@pytest.mark.parametrize(
    "defect",
    [
        "absent",
        "empty",
        "duplicate",
        "untyped",
        "charge_failure",
        "model_failure",
        "preparation_error",
        "no_numerical",
        "one_numerical",
        "numerical_failure",
        "nonfinite",
        "duplicate_frame",
        "late_error",
        "no_bundle",
        "bad_identity",
        "missing_term",
        "duplicate_term",
        "failed_term",
        "term_failure",
        "bad_source",
        "missing_source",
        "wrong_smiles",
        "failed_fd",
    ],
)
def test_fail_closed(report, defect):
    c = report["cases"][0]
    if defect == "absent":
        c.update(case="water", smiles="O")
    elif defect == "empty":
        report["cases"] = []
    elif defect == "duplicate":
        report["cases"].append(deepcopy(c))
    elif defect in ("untyped", "charge_failure", "model_failure"):
        c[
            {
                "untyped": "typed",
                "charge_failure": "charges_complete",
                "model_failure": "model_complete",
            }[defect]
        ] = False
    elif defect in ("preparation_error", "late_error"):
        c["error"] = "injected failure, possibly after both comparisons"
    elif defect == "no_numerical":
        c.pop("numerical")
    elif defect == "one_numerical":
        c["numerical"].pop()
    elif defect == "numerical_failure":
        c["numerical"][1]["passed"] = False
    elif defect == "nonfinite":
        c["numerical"][1]["force_max_error"] = float("nan")
    elif defect == "duplicate_frame":
        for x in c["numerical"]:
            x["frame"] = 0
    elif defect == "no_bundle":
        c.pop("bundle_hashes")
    elif defect == "bad_identity":
        c["bundle_identity"] = None
    elif defect == "missing_term":
        report["term_checks"].pop()
    elif defect == "duplicate_term":
        report["term_checks"][0] = deepcopy(report["term_checks"][1])
    elif defect == "failed_term":
        report["term_checks"][0]["energy_error"] = 10
    elif defect == "term_failure":
        report["term_failure"] = "injected late error"
    elif defect == "bad_source":
        report["source_unchanged"] = False
    elif defect == "missing_source":
        report.pop("source_unchanged")
    elif defect == "wrong_smiles":
        c["smiles"] = "O=O"
    elif defect == "failed_fd":
        report["term_checks"][0]["finite_difference_errors"][-1] = 1
    result = check(report)
    assert not result["passed"] and result["diagnostics"]


def test_artifact_validation_failure_cannot_be_ignored(report, monkeypatch):
    def fail(*args):
        raise ValueError("wrong source/model or missing/changed bundle file")

    monkeypatch.setattr(gates, "verify_artifacts", fail)
    assert not check(report)["passed"]


def test_named_correction_order_and_binding(report):
    baseline = deepcopy(report)
    other = {"case": "water", "smiles": "O"}
    report["cases"].insert(0, other)
    assert gates.corrected_case(baseline, report)["case"] == "dichlorine"
    for key in ("case", "smiles", "model_identity"):
        changed = deepcopy(report)
        changed["cases"][1][key] = "wrong"
        with pytest.raises(ValueError):
            gates.corrected_case(baseline, changed)
    changed = deepcopy(report)
    changed["inventory"] = {}
    with pytest.raises(ValueError):
        gates.corrected_case(baseline, changed)
    for rows in ([], [other], [report["cases"][1]] * 2):
        with pytest.raises(ValueError):
            gates.corrected_case(baseline, {**report, "cases": rows})


def test_operational_summary_recomputes_prerequisites(report):
    gate = check(report)
    workflow = {
        **gate["verified"],
        "passed": True,
        "inputs_unchanged": True,
        "rng_equal": True,
    }
    assert gates.operational_gate(gate, workflow, {"unchanged": True})
    assert not gates.operational_gate({"passed": False}, workflow, {"unchanged": True})
    for key in (
        "parameter_fingerprint",
        "model_fingerprint",
        "passed",
        "inputs_unchanged",
        "rng_equal",
    ):
        assert not gates.operational_gate(
            gate, {**workflow, key: False}, {"unchanged": True}
        )
    assert not gates.operational_gate(
        gate, {**workflow, "error": "failed"}, {"unchanged": True}
    )
