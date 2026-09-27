"""Opt-in real AmberTools integration; normal CI never mocks these references."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.ambertools_integration


def test_real_ambertools_reference_generation(tmp_path):
    if os.environ.get("ISLAND_RUN_AMBERTOOLS_INTEGRATION") != "1":
        pytest.skip("set ISLAND_RUN_AMBERTOOLS_INTEGRATION=1 for real QM/tool execution")
    required = ("antechamber", "parmchk2", "tleap")
    missing = [name for name in required if shutil.which(name) is None]
    if missing:
        pytest.skip(f"AmberTools unavailable: {', '.join(missing)}")
    pytest.importorskip("parmed")
    pytest.importorskip("rdkit")
    script = Path(__file__).parents[1] / "scripts/generate_ambertools_references.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path)],
        check=True, timeout=3600,
    )
    manifest = json.loads((tmp_path / "references.json").read_text())
    assert manifest["status"] == "generated_with_actual_ambertools"
    assert {row["case"] for row in manifest["cases"]} == {
        "phenol_gaff_am1bcc", "phenol_gaff2_am1bcc",
        "phenol_gaff2_provided", "pe_dp3_gaff2_provided",
    }
    for row in manifest["cases"]:
        assert row["record"]["artifact_sha256"]["result.prmtop"]
        assert Path(row["record"]["artifact_dir"]).is_dir()
        checks = row["independent_conversion_checks"]
        assert all(family in checks for family in (
            "bond", "angle", "proper", "improper", "lj", "coulomb",
        ))
        if row["case"].startswith("phenol"):
            assert row["improper_count"] > 0
            assert "converted_kj_mol" in checks["improper"]
