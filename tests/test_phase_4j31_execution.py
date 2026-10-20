"""Opt-in-by-availability real native assignments; no fabricated parameters."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.phase_4j31_contract import require_case
from scripts.phase_4j31_network_fixture import FRC, ROOT, XML


@pytest.mark.parametrize("family", ("pcff", "oplsaa", "gaff", "gaff2"))
@pytest.mark.parametrize("kind", ("branched", "crosslinked"))
def test_real_final_graph_source_completeness_and_exact_execution(
    family, kind, tmp_path
):
    pytest.importorskip("openmm")
    executable = Path(sys.executable)
    amberhome = Path(
        os.environ.get("ISLAND_J31_AMBERHOME", ROOT / "miniforge3/envs/md_rg")
    )
    if family == "pcff" and not FRC.is_file():
        pytest.skip("Pinned PCFF source is unavailable")
    if family == "oplsaa":
        executable = Path(
            os.environ.get(
                "ISLAND_J31_OPLS_PYTHON",
                ROOT / "island-validation/phase4g1-env/bin/python",
            )
        )
        if not executable.is_file() or not XML.is_file():
            pytest.skip("Pinned Foyer environment/source is unavailable")
    if family in {"gaff", "gaff2"} and not (amberhome / "bin/antechamber").is_file():
        pytest.skip("Real AmberTools installation is unavailable")
    output = tmp_path / "acceptance"
    environment = {
        **os.environ,
        "OPENMM_CPU_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "PATH": str(amberhome / "bin") + os.pathsep + os.environ["PATH"],
    }
    command = [
        str(executable),
        "-m",
        "scripts.validate_phase_4j31_network",
        "--output",
        str(output),
        "--family",
        family,
        "--kind",
        kind,
        "--amberhome",
        str(amberhome),
    ]
    process = subprocess.run(
        command, env=environment, capture_output=True, text=True, check=False
    )
    assert process.returncode == 0, process.stdout + process.stderr
    result = json.loads((output / "result.json").read_text())
    require_case(f"{family}/{kind}", result)
    assert result["acceptance_met"] is True
    assert result["production_validated"] is False
    assert result["simulation_readiness"] == "not_established"
    if family == "pcff":
        assert result["status"] == "blocked"
        assert result["dynamics_executed"] is False
        assert result["diagnostic"]["classification"] == "expected_typing_rejection"
    else:
        assert result["status"] == "complete"
        assert result["dynamics_executed"] is True
        assert result["execution"]["comparison"]["passed"]
        assert result["execution"]["frame_count"] == 5
        assert result["execution"]["separate_process_reconstruction"] == "PASS"
        assert result["independent_assignment"]["bonds"] > 0
        assert result["independent_assignment"]["angles"] > 0
        assert result["independent_assignment"]["proper_paths"] > 0
