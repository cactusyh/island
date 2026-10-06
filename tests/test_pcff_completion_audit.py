"""Synthetic source controls for the independent completion-obligation audit."""

import importlib

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
