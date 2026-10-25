"""C13 refuses unsupported qualifications and evidence mutation."""

import copy
import json
import os
from pathlib import Path

import pytest

from scripts.phase_4j40_c11_numeric import signed
from scripts.phase_4j40_c13 import POLICY, characterize, validate_summary


@pytest.fixture
def summary():
    path = (
        Path(os.environ.get("ISLAND_J40_C13_RUN", "/tmp/j40-c13-validation/study"))
        / "c13-summary.json"
    )
    if not path.exists():
        path = Path("docs/evidence/phase_4j40_c13/summary.json")
    if not path.exists():
        pytest.skip("Executed C13 evidence unavailable")
    return json.loads(path.read_text())


def test_executed_matrix(summary):
    assert validate_summary(summary)
    assert len(summary["rows"]) >= 90
    assert len(summary["pressure"]) == 12
    assert any(not r["existing_energy_force_bounds_passed"] for r in summary["rows"])


@pytest.mark.parametrize(
    "key,value",
    [
        ("energy_absolute_kj_mol", 1.0),
        ("force_absolute_kj_mol_nm", 1.0),
        ("pressure_bound_bar", 100),
        ("status", "qualified"),
    ],
)
def test_policy_mutation(summary, key, value):
    mutated = copy.deepcopy(summary)
    mutated["policy"][key] = value
    mutated["policy"] = signed(mutated["policy"])
    with pytest.raises(ValueError, match="policy"):
        validate_summary(signed(mutated))


@pytest.mark.parametrize(
    "key,value",
    [
        ("long_npt", "executed"),
        ("pressure_status", "validated"),
        ("ewald_status", "validated"),
        ("production_validated", True),
        ("scientific_acceptance", True),
    ],
)
def test_no_gate_promotion(summary, key, value):
    with pytest.raises(ValueError, match="readiness"):
        validate_summary(signed(dict(summary, **{key: value})))


@pytest.mark.parametrize(
    "field", ["atoms", "frame_identity", "box_lengths_nm", "metrics", "settings"]
)
def test_row_mutation(summary, field):
    mutated = copy.deepcopy(summary)
    mutated["rows"][0][field] = None
    with pytest.raises(ValueError, match="characterization"):
        validate_summary(signed(mutated))


def test_bounds_unchanged():
    assert POLICY["energy_absolute_kj_mol"] == 1e-5
    assert POLICY["force_absolute_kj_mol_nm"] == 1e-4
    assert POLICY["pressure_bound_bar"] is None


def test_frame_mutation_rejected():
    with pytest.raises(ValueError, match="Mutated frame"):
        characterize(
            {"frame_identity": "forged"}, {"identity": "forged", "atom_ids": [2, 1]}
        )
