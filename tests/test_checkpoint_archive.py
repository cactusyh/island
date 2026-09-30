"""Opt-in retained capped PE and assigned Cl/Br continuation."""

import os
import sys
from pathlib import Path

import pytest


def test_archived_continuation():
    manifest = os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_MANIFEST")
    if not manifest:
        pytest.skip(
            "set ISLAND_AMBERTOOLS_REFERENCE_MANIFEST for archived checkpoint continuation"
        )
    for dependency in ("openmm", "parmed", "rdkit", "scipy"):
        pytest.importorskip(dependency)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from validate_checkpoint_references import run

    report = run(Path(manifest))
    assert all(r["completed"] and r["rng_state_equal"] for r in report["cases"])
