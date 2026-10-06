"""Synthetic source controls for the independent completion-obligation audit."""

import importlib
import json
from collections import Counter
from pathlib import Path

import pytest


@pytest.fixture
def audit(monkeypatch):
    monkeypatch.syspath_prepend("scripts")
    return importlib.import_module("audit_pcff_completion")


def test_increment_paths_preserve_orientation_and_explicit_zero(audit):
    # Labels and numbers are synthetic, not PCFF chemical/reference data.
    rows = audit.read_source(b"""#equivalence cff91
1 1 X X x X X X
1 1 Y Y y Y Y Y
#auto_equivalence cff91_auto
1 1 X X xa X X X X X X X
1 1 Y Y ya Y Y Y Y Y Y Y
#bond_increments cff91_auto
1 1 y x 0.1 -0.2
1 1 xa ya 0 0
""")
    forward = audit.increment_searches(rows, ["X", "Y"])
    reverse = audit.increment_searches(rows, ["Y", "X"])
    assert forward[0]["candidates"] == []
    assert forward[1]["candidates"][0]["endpoint_values"] == ["-0.2", "0.1"]
    assert reverse[1]["candidates"][0]["endpoint_values"] == ["0.1", "-0.2"]
    assert forward[2]["candidates"][0]["endpoint_values"] == ["0", "0"]
    assert forward[1]["equivalence_rows"]
    assert not any(s["candidates"] for s in audit.increment_searches(rows, ["X", "Z"]))


def test_unpinned_source_cannot_claim_completion_obligations(audit):
    with pytest.raises(ValueError, match="exact declared PCFF source"):
        audit.source_obligations(b"synthetic unrelated source", [])


@pytest.mark.parametrize("label", ["Br", "Cl"])
def test_uppercase_halogen_taxonomy_preserves_source_charge_description(audit, label):
    receipt = json.loads(Path("docs/evidence/phase_4j11_obligations.json").read_text())
    row = next(r for r in receipt["unresolved_labels"] if r["label"] == label)
    assert row["description"] == {"Br": "bromine ion", "Cl": "chlorine ion"}[label]
    assert audit.source_domain(label, row["description"]) == "halogen"
    assert not row["implemented_and_independently_verified"]


def test_charge_taxonomy_uses_source_description_not_element_or_label_sign(audit):
    receipt = json.loads(Path("docs/evidence/phase_4j11_obligations.json").read_text())
    rows = {r["label"]: r for r in receipt["unresolved_labels"]}
    assert audit.source_domain("ca+", rows["ca+"]["description"]) == "ionic"
    for label in ("ci", "h+", "hi", "n+", "n1", "ni"):
        assert (
            audit.source_domain(label, rows[label]["description"])
            == "charged_environment"
        )
    for label in ("s-", "n2"):
        assert (
            audit.source_domain(label, rows[label]["description"])
            == "organic_alias_or_specialized_environment"
        )
    assert audit.source_domain("Na", rows["Na"]["description"]) == "metal"
    # Synthetic description controls: labels/element names alone do not imply ions.
    assert (
        audit.source_domain("ca+", "synthetic neutral environment")
        == "organic_alias_or_specialized_environment"
    )
    assert audit.source_domain("synthetic", "explicit ion") == "ionic"


def test_corrected_receipt_counts_and_unresolved_status(audit):
    receipt = json.loads(Path("docs/evidence/phase_4j11_obligations.json").read_text())
    rows = receipt["unresolved_labels"]
    assert len(rows) == 47
    assert all(not r["implemented_and_independently_verified"] for r in rows)
    assert not receipt["full_source_acceptance"]
    assert receipt["production_validated"] is False
    assert receipt["simulation_readiness"] == "not_established"
    expected = {
        "halogen": 2,
        "ionic": 1,
        "charged_environment": 6,
        "metal": 16,
        "zeolite_surface": 10,
        "organic_alias_or_specialized_environment": 12,
    }
    assert receipt["unresolved_domain_counts"] == expected
    assert Counter(r["domain"] for r in rows) == expected
    assert all(
        r["domain"] == audit.source_domain(r["label"], r["description"]) for r in rows
    )
