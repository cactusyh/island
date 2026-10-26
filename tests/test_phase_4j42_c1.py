"""Strict blocked-record validation; no external parameters or execution assets."""

import json
import shutil
import subprocess
import sys
from dataclasses import replace

import pytest

from island.charge_references.records import pack
from island.forcefields import (
    ExternalPDMSBlockedError,
    ExternalPDMSPreparation,
    create_evaluator,
)
from island.forcefields.external_pdms import bundle_template
from scripts.phase_4j42_c1 import CASES, execute, fixture, prepare, reconstruct


@pytest.fixture
def outcomes(tmp_path):
    execute(tmp_path / "study")
    return tmp_path / "study"


@pytest.mark.parametrize("case", CASES)
def test_public_boundaries_preserve_full_result_and_diagnostics(case, tmp_path):
    execute(tmp_path / "study")
    folder = tmp_path / "study" / case
    result = ExternalPDMSPreparation((folder / "preparation.json").read_text())
    system = fixture(case)
    result.validate_integrity(system)
    with pytest.raises(ExternalPDMSBlockedError) as caught:
        create_evaluator(system, result)
    error = caught.value
    assert error.result is result
    assert error.status == result.payload["status"]
    assert error.diagnostics == result.payload["diagnostics"]
    assert error.validation_diagnostic is None
    error.diagnostics.clear()
    assert result.payload["diagnostics"]


@pytest.mark.parametrize("case", CASES[:5])
def test_typing_blocked_graph_is_not_revalidated(case, tmp_path, monkeypatch):
    import island.forcefields.external_pdms as schema

    system = fixture(case)
    result = prepare(system, tmp_path, {"file": None, "pin": None})
    assert result.payload["status"] == "typing_blocked"
    assert result.payload["requirements"] is None

    def forbidden(*args):
        pytest.fail("Invalid graph was validated a second time")

    monkeypatch.setattr(schema, "requirements", forbidden)
    result.validate_integrity(system)
    with pytest.raises(ExternalPDMSBlockedError) as caught:
        create_evaluator(system, result)
    assert caught.value.status == "typing_blocked"
    assert caught.value.diagnostics[0]["stage"] == "final_graph_invalid"


@pytest.mark.parametrize(
    "status", ["prepared", "asset_unavailable", "parameter_incomplete"]
)
def test_resigned_typing_status_mutation_is_rejected(status, tmp_path):
    result = prepare(fixture("invalid_valence"), tmp_path, {"file": None, "pin": None})
    p = result.payload
    p["status"] = status
    with pytest.raises(ValueError):
        replace(result, json_text=pack(p))


@pytest.mark.parametrize(
    "mutation", ["inventory", "stage", "context", "missing", "dependencies", "import"]
)
def test_typing_blocked_record_cannot_acquire_fake_evidence(mutation, tmp_path):
    result = prepare(fixture("invalid_valence"), tmp_path, {"file": None, "pin": None})
    p = result.payload
    if mutation == "inventory":
        p["requirements"] = bundle_template(fixture("asset_unavailable"))["coverage"]
    elif mutation == "stage":
        p["diagnostics"][0]["stage"] = "asset_unavailable"
    elif mutation == "context":
        p["diagnostics"][0]["context"] = ""
    elif mutation == "missing":
        p["missing_interactions"] = {"fake": {}}
    elif mutation == "dependencies":
        p["unresolved_dependencies"] = ["fake"]
    else:
        p["import_validation"] = {"status": "parameter_bundle_validated"}
    with pytest.raises(ValueError, match="typing_blocked_record_changed"):
        replace(result, json_text=pack(p))


@pytest.mark.parametrize(
    "case", ["asset_unavailable", "missing_hash", "malformed_bundle"]
)
def test_valid_graph_status_rechecks_inventory_and_rejects_stale_graph(case, outcomes):
    folder = outcomes / case
    result = ExternalPDMSPreparation((folder / "preparation.json").read_text())
    system = fixture(case)
    system.topology.sites[min(system.topology.sites)].metadata["chain_id"] = "stale"
    with pytest.raises(ValueError, match="stale_final_graph"):
        result.validate_integrity(system)
    with pytest.raises(ExternalPDMSBlockedError) as caught:
        create_evaluator(system, result)
    assert caught.value.result is result
    assert caught.value.validation_diagnostic["context"] == "stale_final_graph"


def test_validation_errors_remain_classified_after_constructor_bypass(tmp_path):
    result = prepare(fixture("invalid_valence"), tmp_path, {"file": None, "pin": None})
    p = result.payload
    p["status"] = "prepared"
    object.__setattr__(result, "json_text", pack(p))
    with pytest.raises(ExternalPDMSBlockedError) as caught:
        create_evaluator(fixture("invalid_valence"), result)
    assert (
        caught.value.result is result and caught.value.status == "parameter_incomplete"
    )
    assert caught.value.validation_diagnostic["stage"] == "blocked_record_invalid"
    object.__setattr__(result, "json_text", "{")
    with pytest.raises(ExternalPDMSBlockedError):
        create_evaluator(fixture("invalid_valence"), result)


def test_malformed_system_attribute_is_typing_blocked(tmp_path):
    result = prepare(object(), tmp_path, {"file": None, "pin": None})
    assert result.payload["status"] == "typing_blocked"
    with pytest.raises(ExternalPDMSBlockedError):
        create_evaluator(object(), result)


@pytest.mark.parametrize(
    "case", ["asset_unavailable", "missing_hash", "invalid_valence"]
)
def test_workflow_and_minimization_never_enter_numerical_backends(
    case, outcomes, tmp_path, monkeypatch
):
    from island import evaluation
    from island.forcefields import PreparedForceField, save_prepared_forcefield
    from island.forcefields.bundles import PreparedBundleError
    from island.workflows import prepared as workflow

    def forbidden(*a, **kw):
        pytest.fail("Blocked external preparation entered numerical execution")

    monkeypatch.setattr(evaluation, "PCFFSinglePointEvaluator", forbidden)
    monkeypatch.setattr(evaluation, "OpenMMSinglePointEvaluator", forbidden)
    for name in (
        "minimize_geometry",
        "initialize_velocities",
        "run_dynamics_segment",
        "create_dynamics_checkpoint",
    ):
        monkeypatch.setattr(workflow, name, forbidden)
    monkeypatch.setattr(subprocess, "run", forbidden)
    folder = outcomes / case
    result = ExternalPDMSPreparation((folder / "preparation.json").read_text())
    system = fixture(case)
    assert not isinstance(result, PreparedForceField)
    with pytest.raises(ExternalPDMSBlockedError):
        workflow.minimize_geometry(system, workflow.create_evaluator(system, result))
    target = tmp_path / "prepared"
    with pytest.raises(
        PreparedBundleError, match="Validated PreparedForceField required"
    ):
        save_prepared_forcefield(system, result, target)
    assert not target.exists()
    assert not list(tmp_path.rglob("*checkpoint*"))
    assert not list(tmp_path.rglob("*trajectory*"))


def test_all_blocked_statuses_reconstruct_after_fresh_process_relocation(
    outcomes, tmp_path
):
    target = tmp_path / "relocated"
    shutil.copytree(outcomes, target)
    child = subprocess.run(
        [sys.executable, "scripts/phase_4j42_c1.py", "--reconstruct", str(target)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr
    assert json.loads(child.stdout) == reconstruct(outcomes)
    for p in outcomes.rglob("*"):
        if p.is_file():
            assert p.read_bytes() == (target / p.relative_to(outcomes)).read_bytes()
