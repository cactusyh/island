"""Opt-in archived real references, distinct from offline numerics and regeneration."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.archived_amber_reference


def test_archived_real_singlepoint_references(tmp_path):
    configured = os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_MANIFEST")
    if not configured:
        pytest.skip("set ISLAND_AMBERTOOLS_REFERENCE_MANIFEST to the retained manifest")
    manifest = Path(configured)
    if not manifest.is_file():
        pytest.skip(f"retained reference manifest unavailable: {manifest}")
    pytest.importorskip("openmm")
    pytest.importorskip("parmed")
    output = tmp_path / "acceptance.json"
    script = Path(__file__).parents[1] / "scripts/validate_singlepoint_references.py"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--preparation-manifest",
            str(manifest),
            "--output",
            str(output),
        ],
        check=True,
        timeout=300,
    )
    records = json.loads(output.read_text())["cases"]
    assert len(records) == 15
    assert {row["case"] for row in records} == {
        "phenol_gaff_am1bcc",
        "phenol_gaff2_am1bcc",
        "phenol_gaff2_provided",
        "pe_dp3_gaff2_provided",
        "halomethane_gaff2_provided",
    }
    for row in records:
        assert row["historical_signatures_verified"] is True
        report = row["charge_report"]
        assert report["charges_modified"] is False
        if row["case"].endswith("am1bcc"):
            assert report["charge_tolerance_e"] == 0.002
            assert abs(report["final_charge_residual_e"]) <= 0.002
