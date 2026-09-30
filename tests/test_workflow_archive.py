"""Opt-in real archived-parameter workflow; no live parameterization."""

import os
import sys
from pathlib import Path

import pytest


def test_archived_pe_workflow(tmp_path):
    manifest = os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_MANIFEST")
    if not manifest:
        pytest.skip(
            "set ISLAND_AMBERTOOLS_REFERENCE_MANIFEST for archived workflow acceptance"
        )
    for name in ("openmm", "parmed", "rdkit", "scipy"):
        pytest.importorskip(name)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from validate_workflow import run

    report = run(tmp_path / "workflow", Path(manifest))
    assert report["acceptance"] == "passed"
    assert report["rng_state_equal"] and report["steps"] == 60
