"""Injected software outcomes test CLI gates, not scientific integration."""

import copy

import pytest

from island.workflows import storage
from scripts import validate_long_chain_provided as cli


@pytest.mark.parametrize(
    "passed,workflow,required,expected",
    [
        (0, False, False, 1),
        (2, False, False, 1),
        (3, False, False, 0),
        (3, False, True, 1),
        (3, True, True, 0),
    ],
)
def test_requested_cli_gate(
    tmp_path, monkeypatch, passed, workflow, required, expected
):
    rows = [
        {
            "case": case[0],
            "status": (
                "parameterization_and_numerical_acceptance_passed"
                if i < passed
                else "failed"
            ),
        }
        for i, case in enumerate(cli.CASES)
    ]
    rows[-1].update(
        workflow={
            "status": "completed" if workflow else "stage_failed",
            "accepted_step": 4 if workflow else 0,
        },
        analyzed_samples=5 if workflow else 0,
        analysis_verified=workflow,
    )
    original = copy.deepcopy(rows)

    def injected(root, amberhome):
        root.mkdir()
        return rows

    monkeypatch.setattr(cli, "run", injected)
    output = tmp_path / "out"
    args = ["--output", str(output), "--amberhome", "/missing/ambertools"]
    if required:
        args.append("--require-workflow")
    assert cli.main(args) == expected
    report = storage.read_json(output / "report.json")
    assert report["requested_gate_passed"] is (expected == 0)
    assert report["cases"] == original == rows
    assert report["gates"]["parameterization_numerical"] is (passed == 3)


def test_missing_cases_cannot_pass():
    assert not any(cli.acceptance_gates([]).values())
