"""Explicitly synthetic aggregate controls; no scientific calculations."""

import sys
from copy import deepcopy
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
from summarize_pcff_organic_domains import MANDATORY, numerical_gate


def software_report():
    case = {
        "passed": True,
        "bundle_hashes": {"manifest.json": "a" * 64},
        "prepared_identity": "b" * 64,
        "model_identity": "c" * 64,
        "parameter_inventory_verified": True,
        "offline_loading": True,
        "child_pid": 2,
        "parent_pid": 1,
        "reconstruction_max_error": 0.0,
        "frames": [
            {
                "energy_error": 0.0,
                "force_max_error": 0.0,
                "finite_difference_errors": [1e-6, 1e-8, 1e-7],
            }
            for _ in range(2)
        ],
    }
    return {
        "source_unchanged": True,
        "cases": [dict(deepcopy(case), case=name) for name in MANDATORY],
    }


def test_positive_and_order_independence():
    r = software_report()
    r["cases"].reverse()
    assert numerical_gate(r, MANDATORY)


@pytest.mark.parametrize(
    "kind",
    [
        "empty",
        "missing",
        "duplicate",
        "failed",
        "late_error",
        "bundle",
        "source",
        "frame",
        "fd",
        "same_process",
        "delta",
        "names",
    ],
)
def test_fail_closed(kind):
    r = software_report()
    c = r["cases"][0]
    names = MANDATORY
    if kind == "empty":
        r["cases"] = []
    elif kind == "missing":
        r["cases"].pop()
    elif kind == "duplicate":
        r["cases"][1] = deepcopy(c)
    elif kind == "failed":
        c["passed"] = False
    elif kind == "late_error":
        c["error"] = "bundle publication failed"
    elif kind == "bundle":
        c["bundle_hashes"] = {}
    elif kind == "source":
        r["source_unchanged"] = False
    elif kind == "frame":
        c["frames"].pop()
    elif kind == "fd":
        c["frames"][0]["finite_difference_errors"][-1] = float("nan")
    elif kind == "same_process":
        c["child_pid"] = c["parent_pid"]
    elif kind == "delta":
        c["reconstruction_max_error"] = 0.001
    elif kind == "names":
        names = []
        r["cases"] = []
    assert not numerical_gate(r, names)
