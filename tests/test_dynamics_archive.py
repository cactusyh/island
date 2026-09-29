"""Opt-in NVE acceptance using external checksum-verified Amber references."""

import os
from pathlib import Path

import pytest


@pytest.mark.archived_amber_reference
def test_archived_nve(monkeypatch):
    manifest = os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_MANIFEST")
    if not manifest:
        pytest.skip(
            "set ISLAND_AMBERTOOLS_REFERENCE_MANIFEST for archived NVE acceptance"
        )
    for dependency in ("openmm", "scipy", "parmed", "rdkit"):
        pytest.importorskip(dependency)
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    from validate_dynamics_references import OPTIONS, run

    report = run(Path(manifest))
    assert {row["case"] for row in report["cases"]} == {
        "phenol_gaff_am1bcc",
        "phenol_gaff2_am1bcc",
        "phenol_gaff2_provided",
        "pe_dp3_gaff2_provided",
        "halomethane_gaff2_provided",
    }
    for row in report["cases"]:
        assert row["minimization_converged"] and row["completed"], row
        assert row["completed_steps"] == OPTIONS.steps
        assert row["evaluations"] == OPTIONS.steps + 2
        assert row["max_abs_energy_deviation_kj_mol"] <= OPTIONS.max_energy_deviation
