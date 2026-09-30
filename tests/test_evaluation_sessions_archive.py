"""Opt-in session/fresh/original-prmtop acceptance, including 1000 fixed steps."""

import json
import os
from pathlib import Path

import pytest


@pytest.mark.archived_amber_reference
def test_archived_sessions(monkeypatch):
    supplied = os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_MANIFEST")
    if not supplied:
        pytest.skip(
            "set ISLAND_AMBERTOOLS_REFERENCE_MANIFEST for session archive acceptance"
        )
    for dependency in ("openmm", "parmed", "rdkit", "scipy"):
        pytest.importorskip(dependency)
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    from validate_evaluation_sessions import acceptance
    from validate_singlepoint_references import load_archived_preparation

    manifest = Path(supplied)
    cases = json.loads(manifest.read_text())["cases"]
    assert len(cases) == 5
    for row in cases:
        system, preparation, directory = load_archived_preparation(manifest.parent, row)
        report, _, _ = acceptance(system, preparation, directory, row["case"])
        assert report["minimization"]["session_contexts"] == 2
        for run in report["trajectories"]:
            assert run["termination_reason"] == "completed"
            assert run["session_contexts"] == 2
            assert run["evaluations"] == run["steps"] + 2
        if row["case"] == "phenol_gaff2_am1bcc":
            assert report["trajectories"][1]["steps"] == 1000
