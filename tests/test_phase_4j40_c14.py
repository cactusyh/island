"""Backend scopes, fail-closed permission and unchanged C13 evidence."""

import copy
import json
from pathlib import Path

import pytest

from island.analysis.backend_readiness import (
    BINDINGS,
    OPENMM,
    SCOPES,
    assess_periodic_backend_readiness,
    require_openmm_periodic_workflow,
    validate_backend_readiness,
)
from island.exceptions import AnalysisError
from scripts.phase_4j40_c11_numeric import signed
from scripts.phase_4j40_c14 import build_report, load_source
from scripts.validate_phase_4j40_c14_receipt import validate_report

SOURCE = Path("docs/evidence/phase_4j40_c13")


def inputs():
    bindings = {key: signed({"unit_test_binding": key})["identity"] for key in BINDINGS}
    checks = {
        key: {
            "name": key,
            "scope": scope,
            "bindings": bindings.copy(),
            "status": "passed",
            "reason": "Synthetic contract test only, not scientific evidence",
            "evidence_sha256": ["a" * 64],
        }
        for key, scope in SCOPES.items()
    }
    return bindings, checks


def test_openmm_scope_does_not_require_lammps_equivalence():
    bindings, checks = inputs()
    for key in checks:
        if not key.startswith("openmm."):
            checks[key]["status"] = "blocked"
    r = assess_periodic_backend_readiness(bindings=bindings, checks=checks)
    require_openmm_periodic_workflow(r, expected_bindings=bindings)
    assert r["cross_engine_equivalence_status"] == "ewald_unvalidated"
    assert r["production_validated"] is False
    assert r["scientific_npt_acceptance"] is False
    assert r["scientific_equilibrium_status"] == "not_established"


@pytest.mark.parametrize("missing", OPENMM)
def test_each_openmm_requirement_is_mandatory(missing):
    bindings, checks = inputs()
    del checks["openmm." + missing]
    r = assess_periodic_backend_readiness(bindings=bindings, checks=checks)
    with pytest.raises(AnalysisError, match="workflow blocked"):
        require_openmm_periodic_workflow(r, expected_bindings=bindings)


def test_export_and_runtime_cannot_claim_equivalence():
    bindings, checks = inputs()
    only = {k: v for k, v in checks.items() if k.startswith("lammps")}
    r = assess_periodic_backend_readiness(bindings=bindings, checks=only)
    assert r["lammps_periodic_export_status"] == "exact_frame_exported_only"
    assert r["lammps_kspace_runtime_status"] == "native_kspace_executed_only"
    assert r["openmm_periodic_backend_status"] == "software_tests_only"
    assert r["cross_engine_equivalence_status"] == "ewald_unvalidated"
    assert r["pressure_validation_status"] == "pressure_unavailable"


@pytest.mark.parametrize("key", BINDINGS)
def test_stale_binding_rejected(key):
    bindings, checks = inputs()
    checks["openmm.exact_graph_parameters"]["bindings"][key] = "b" * 64
    with pytest.raises(AnalysisError, match="Stale"):
        assess_periodic_backend_readiness(bindings=bindings, checks=checks)


@pytest.mark.parametrize(
    "field,value",
    [
        ("scope", "lammps_export"),
        ("status", True),
        ("evidence_sha256", []),
        ("name", "lammps_export.exact_frame"),
    ],
)
def test_wrong_evidence_scope_and_invalid_witness(field, value):
    bindings, checks = inputs()
    checks["openmm.exact_graph_parameters"][field] = value
    with pytest.raises(AnalysisError):
        assess_periodic_backend_readiness(bindings=bindings, checks=checks)


@pytest.mark.parametrize(
    "field",
    [
        "openmm_periodic_backend_status",
        "cross_engine_equivalence_status",
        "pressure_validation_status",
        "scientific_equilibrium_status",
        "scientific_npt_acceptance",
        "production_validated",
        "simulation_readiness",
    ],
)
def test_persistence_rejects_even_resigned_decision_mutation(field):
    bindings, _ = inputs()
    r = assess_periodic_backend_readiness(bindings=bindings, checks={})
    r[field] = "forged"
    with pytest.raises(AnalysisError, match="Changed"):
        validate_backend_readiness(signed(r), expected_bindings=bindings)


def test_caller_mutation_does_not_change_record():
    bindings, checks = inputs()
    r = assess_periodic_backend_readiness(bindings=bindings, checks=checks)
    checks["openmm.internal_energy_force"]["status"] = "blocked"
    assert validate_backend_readiness(r, expected_bindings=bindings)


@pytest.fixture(scope="module")
def report():
    return build_report(SOURCE)


def test_retained_source_counts_and_residuals(report):
    summary = json.loads((SOURCE / "summary.json").read_text())
    assert report["comparisons"] == summary["rows"]
    assert report["native_evaluations"] == 300
    assert report["comparison_count"] == 108
    assert report["passing_comparisons"] == 13
    assert report["new_qualified_equivalence_policy"] is None
    for cell in report["cells"]:
        assert not cell["decision"]["openmm_periodic_workflow_permitted"]
        assert cell["decision"]["pressure_validation_status"] == "pressure_unavailable"
        assert (
            cell["decision"]["lammps_kspace_runtime_status"]
            == "native_kspace_executed_only"
        )


@pytest.mark.parametrize(
    "field,value",
    [
        ("comparison_count", 109),
        ("passing_comparisons", 108),
        ("native_evaluations", 301),
        ("cross_engine_equivalence_status", "validated"),
        ("pressure_validation_status", "validated"),
        ("new_qualified_equivalence_policy", {}),
        ("scientific_npt_acceptance", True),
    ],
)
def test_report_mutation_cannot_promote_evidence(report, field, value):
    mutated = copy.deepcopy(report)
    mutated[field] = value
    with pytest.raises(ValueError, match="Changed backend report"):
        validate_report(signed(mutated), SOURCE)


def test_c13_receipt_replacement_rejected(tmp_path):
    r = json.loads((SOURCE / "receipt.json").read_text())
    r["individual_existing_bound_passes"] = 108
    (tmp_path / "receipt.json").write_text(json.dumps(signed(r)))
    with pytest.raises(ValueError, match="Changed C13 receipt"):
        load_source(tmp_path)


def test_cross_engine_requires_its_own_complete_evidence():
    bindings, checks = inputs()
    checks["cross_engine.numerical_policy"]["status"] = "blocked"
    r = assess_periodic_backend_readiness(bindings=bindings, checks=checks)
    assert r["openmm_periodic_workflow_permitted"]
    assert r["cross_engine_equivalence_status"] == "ewald_unvalidated"


def test_boolean_type_mutation_cannot_bypass_identity():
    bindings, _ = inputs()
    r = assess_periodic_backend_readiness(bindings=bindings, checks={})
    r["production_validated"] = 0
    with pytest.raises(AnalysisError, match="Changed"):
        validate_backend_readiness(r, expected_bindings=bindings)


@pytest.mark.parametrize("section", ["comparisons", "source_hashes", "cells"])
def test_residuals_sources_and_system_identity_are_immutable(report, section):
    mutated = copy.deepcopy(report)
    if section == "comparisons":
        mutated[section][0]["metrics"]["absolute_energy_error_kj_mol"] = 0.0
    elif section == "source_hashes":
        mutated[section]["artifacts.tar.xz"] = "a" * 64
    else:
        mutated[section][0]["decision"]["bindings"]["charges"] = "b" * 64
    with pytest.raises(ValueError, match="Changed backend report"):
        validate_report(signed(mutated), SOURCE)
