"""Opt-in archive acceptance; ordinary tests require no external AmberTools data."""

import importlib.util
import os
from pathlib import Path

import pytest


@pytest.mark.archived_amber_reference
def test_archived_minimization(monkeypatch):
    manifest = os.environ.get("ISLAND_AMBERTOOLS_REFERENCE_MANIFEST")
    if not manifest:
        pytest.skip("set ISLAND_AMBERTOOLS_REFERENCE_MANIFEST for retained real cases")
    pytest.importorskip("scipy")
    pytest.importorskip("openmm")
    pytest.importorskip("parmed")
    pytest.importorskip("rdkit")
    scripts = Path(__file__).resolve().parents[1] / "scripts"
    monkeypatch.syspath_prepend(str(scripts))
    spec = importlib.util.spec_from_file_location(
        "minimization_acceptance", scripts / "validate_minimization_references.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = module.run(Path(manifest))
    required = {
        "phenol_gaff2_am1bcc",
        "pe_dp3_gaff2_provided",
        "halomethane_gaff2_provided",
    }
    assert required <= {row["case"] for row in report["cases"]}
    for row in report["cases"]:
        assert row["converged"], row
        assert row["final_fmax_kj_mol_angstrom"] <= row["options"]["force_tolerance"]
        assert row["initial_fmax_kj_mol_angstrom"] > 1
        assert row["final_energy_kj_mol"] < row["initial_energy_kj_mol"]
