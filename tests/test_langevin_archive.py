"""Opt-in checksummed five-molecule finite-temperature acceptance."""

import os
import sys
from pathlib import Path

import pytest


def test_archived_thermal_cases():
    manifest = os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_MANIFEST")
    if not manifest:
        pytest.skip(
            "set ISLAND_AMBERTOOLS_REFERENCE_MANIFEST for archived Langevin acceptance"
        )
    for dependency in ("openmm", "parmed", "rdkit", "scipy"):
        pytest.importorskip(dependency)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from validate_langevin_references import run

    report = run(Path(manifest))
    assert len(report["cases"]) == 5
    assert all(
        row["completed"]
        and row["final_verified"]
        and row["contexts"] == 2
        and row["evaluations"] == 202
        for row in report["cases"]
    )
