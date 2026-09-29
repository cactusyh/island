"""Opt-in real AmberTools integration; normal CI never mocks these references."""

import json
import os
import shutil
import subprocess
import sys
from hashlib import sha256
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
    output = Path(os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_OUTPUT", str(tmp_path)))
    subprocess.run(
        [sys.executable, str(script), "--output", str(output)],
        check=True, timeout=3600,
    )
    manifest = json.loads((output / "references.json").read_text())
    assert manifest["status"] == "generated_with_actual_ambertools"
    assert {row["case"] for row in manifest["cases"]} == {
        "phenol_gaff_am1bcc", "phenol_gaff2_am1bcc",
        "phenol_gaff2_provided", "pe_dp3_gaff2_provided",
        "halomethane_gaff2_provided",
    }
    for row in manifest["cases"]:
        assert row["record"]["artifact_sha256"]["result.prmtop"]
        assert Path(row["record"]["artifact_dir"]).is_dir()
        assert row["source_exclusion_count"] > 0
        assert row["source_14_pair_count"] >= 0
        directory = Path(row["record"]["artifact_dir"])
        assert directory == output.resolve() / row["artifact_subdirectory"]
        for name, expected in row["retained_file_sha256"].items():
            assert sha256((directory / name).read_bytes()).hexdigest() == expected
        assert (directory / "input_system.json").is_file()
        input_lines = (directory / "input.mol2").read_text().splitlines()
        typed_lines = (directory / "typed.mol2").read_text().splitlines()

        def mol2_names(lines):
            begin = lines.index("@<TRIPOS>ATOM") + 1
            end = lines.index("@<TRIPOS>BOND")
            return {line.split()[1] for line in lines[begin:end]}

        expected_names = set(row["record"]["lineage"]["input_name_to_site_id"])
        assert mol2_names(input_lines) == mol2_names(typed_lines) == expected_names
        if row["case"].endswith("provided"):
            assert not (directory / "sqm.out").exists()
        else:
            assert (directory / "sqm.out").is_file()
        checks = row["independent_conversion_checks"]
        assert all(family in checks for family in (
            "bond", "angle", "proper", "improper", "lj", "coulomb",
        ))
        if row["case"].startswith("phenol"):
            assert row["improper_count"] > 0
            assert "converted_kj_mol" in checks["improper"]
        if row["case"].startswith("halomethane"):
            assert any(name.startswith("Cl") for name in expected_names)
            assert any(name.startswith("Br") for name in expected_names)
            assert len(row["expected_cip_by_site"]) == 1
