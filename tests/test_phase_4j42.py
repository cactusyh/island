"""J42 integration controls: no complete supplied asset, no numerical authority."""

import json
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace
from hashlib import sha256

import pytest

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import IslandError
from island.forcefields import (
    ExternalPDMSBlockedError,
    ExternalPDMSOptions,
    ExternalPDMSPreparation,
    ForceFieldRequest,
    ForceFieldRequestError,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.external_pdms import LIMITS, bundle_template, requirements
from island.forcefields.external_pdms_adapter import FAMILY
from scripts.phase_4j41_c5 import CASES, MUTATIONS, fixture, mutated_record
from scripts.phase_4j42 import prepare


@pytest.fixture(scope="module")
def system():
    return fixture(CASES[1])


@pytest.mark.parametrize("case", CASES, ids=[c[0] for c in CASES])
def test_public_preparation_preserves_full_graph_and_dependency_inventory(case):
    s = fixture(case)
    p = prepare(s).payload
    r = requirements(s)
    assert p["status"] == "asset_unavailable" and p["requirements"] == r
    assert p["missing_interactions"] == r["interactions"]
    expected = [
        {"consumer": k, "dependency": d}
        for k, v in r["interactions"].items()
        for d in v["dependencies"]
    ]
    assert p["unresolved_dependencies"] == expected
    assert p["limitations"] == LIMITS
    assert p["parameter_assignment"] is p["execution_profile"] is None
    assert p["backend_executions"] == [] and not p["execution_authorized"]
    assert p["native_pcff_acceptance_gate"] == p["acceptance_gate"] == 1


def test_no_pin_means_no_candidate_read(system, tmp_path, monkeypatch):
    import island.forcefields.external_pdms as schema

    def forbidden(*a, **kw):
        pytest.fail("Candidate opened without an independent pin")

    monkeypatch.setattr(schema, "load_external_pdms_bundle", forbidden)
    path = tmp_path / "candidate.json"
    for pin in (None, "", "bad", "0" * 63):
        p = prepare(system, path, pin).payload
        assert p["status"] == "parameter_incomplete"
        assert p["diagnostics"] == [{"stage": "independent_bundle_sha256_required"}]


def test_no_implicit_discovery_or_source_fallback(system, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "pdms.json").write_text("do not discover")
    (tmp_path / "pcff.frc").write_text("do not discover")
    assert prepare(system).payload["status"] == "asset_unavailable"


@pytest.mark.parametrize("mutation", MUTATIONS)
def test_c5_mutations_rejected_through_public_dispatch(system, tmp_path, mutation):
    raw, stage = mutated_record(bundle_template(system), mutation)
    path = tmp_path / "bundle.json"
    path.write_text(raw)
    p = prepare(system, path, sha256(raw.encode()).hexdigest()).payload
    assert p["status"] == "parameter_incomplete"
    assert stage in {d["stage"] for d in p["diagnostics"]}
    assert p["backend_executions"] == [] and p["parameter_assignment"] is None


@pytest.mark.parametrize("raw", [b"{", b"\xff", b"null", b"[]"])
def test_malformed_external_files_return_classified_results(system, tmp_path, raw):
    path = tmp_path / "invalid.json"
    path.write_bytes(raw)
    p = prepare(system, path, sha256(raw).hexdigest()).payload
    assert p["status"] == "parameter_incomplete"
    assert len(p["missing_interactions"]) == len(requirements(system)["interactions"])


def test_wrong_pin_cannot_be_repaired_by_resigning_the_candidate(system, tmp_path):
    path = tmp_path / "bundle.json"
    path.write_text(pack(bundle_template(system)))
    p = prepare(system, path, "0" * 64).payload
    assert "bundle_pin_mismatch" in {d["stage"] for d in p["diagnostics"]}
    assert p["execution_profile"] is None


@pytest.mark.parametrize("text", ["[SiH2]O", "C[Si](C)(C)[O-]", "COP(=O)(OC)OC"])
def test_wrong_final_chemistry_returns_typing_blocked(text):
    p = prepare(from_smiles(text, random_seed=2026)).payload
    assert p["status"] == "typing_blocked"
    assert p["requirements"] is None
    assert p["diagnostics"][0]["stage"] == "final_graph_invalid"


def test_blocked_result_cannot_reach_evaluation_or_minimization(system, monkeypatch):
    from island import evaluation
    from island.minimization import minimize_geometry

    def forbidden(*a, **kw):
        pytest.fail("Numerical engine invoked for blocked external preparation")

    monkeypatch.setattr(evaluation, "PCFFSinglePointEvaluator", forbidden)
    result = prepare(system)
    with pytest.raises(ExternalPDMSBlockedError, match="asset_unavailable") as error:
        create_evaluator(system, result)
    assert error.value.result is result
    # The existing workflow creates its evaluator before minimization or MD.
    with pytest.raises(ExternalPDMSBlockedError):
        minimize_geometry(system, create_evaluator(system, result))


def test_owned_report_and_rejected_execution_authority(system):
    result = prepare(system)
    p = result.payload
    p["execution_authorized"] = True
    with pytest.raises(ValueError, match="cannot_grant_execution_authority"):
        replace(result, json_text=pack(p))
    assert not result.payload["execution_authorized"]
    with pytest.raises(FrozenInstanceError):
        result.json_text = "changed"
    raw = json.loads(result.json_text)
    raw["payload"]["status"] = "prepared"
    with pytest.raises(IslandError):
        ExternalPDMSPreparation(json.dumps(raw))


def test_stale_component_and_graph_rejected_before_backend(system):
    from island.workflows.bundle import system_data, system_from

    changed = system_from(system_data(system))
    changed.topology.sites[min(changed.topology.sites)].metadata["chain_id"] = "stale"
    with pytest.raises(ExternalPDMSBlockedError) as error:
        create_evaluator(changed, prepare(system))
    assert error.value.validation_diagnostic["context"] == "stale_final_graph"


def test_external_options_cannot_impersonate_native_family():
    with pytest.raises(ForceFieldRequestError, match="Options do not belong"):
        ForceFieldRequest("pcff", ExternalPDMSOptions())
    with pytest.raises(ForceFieldRequestError, match="native records"):
        ForceFieldRequest(FAMILY, ExternalPDMSOptions(), typed_graph=object())


def test_request_revalidates_tampered_options(system):
    request = ForceFieldRequest(FAMILY, ExternalPDMSOptions())
    object.__setattr__(request.options, "bundle_path", 17)
    with pytest.raises(ForceFieldRequestError, match="Invalid backend options"):
        prepare_forcefield(system, request)


def test_adapter_selection_stays_lazy():
    code = """import sys
from island.forcefields import ExternalPDMSOptions, ForceFieldRequest
ForceFieldRequest('external_pdms_class2', ExternalPDMSOptions())
assert not set(('rdkit','foyer','parmed','openmm','scipy')) & sys.modules.keys()
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_relocation_and_fresh_process_reconstruct_blocked_boundary(tmp_path):
    import shutil

    from scripts.phase_4j42 import execute, reconstruct

    source, target = tmp_path / "original", tmp_path / "relocated"
    report = execute(source)
    assert [c["atoms"] for c in report["cases"]] == [9, 27, 16, 26, 46]
    shutil.copytree(source, target)
    expected = reconstruct(source)
    child = subprocess.run(
        [sys.executable, "scripts/phase_4j42.py", "--reconstruct", str(target)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr
    assert json.loads(child.stdout) == expected
    for path in source.rglob("*"):
        if path.is_file():
            assert path.read_bytes() == (target / path.relative_to(source)).read_bytes()


def test_full_source_acceptance_stays_nonzero_without_exception(tmp_path):
    r = subprocess.run(
        [sys.executable, "scripts/phase_4j42.py", "--output", str(tmp_path / "study")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 1 and "Traceback" not in r.stderr
    p = json.loads(r.stdout)
    assert p["execution_profiles"] == [] and p["native_pcff_acceptance_gate"] == 1
    assert all(c["status"] == "asset_unavailable" for c in p["cases"])


@pytest.mark.parametrize("target", ["receipt", "manifest", "archive"])
def test_tampered_j42_publication_rejected(tmp_path, target):
    from scripts.phase_4j42 import BASE, BRANCH
    from scripts.validate_phase_4j41_receipt import sha, signed
    from scripts.validate_phase_4j42_receipt import validate

    manifest = signed(
        {"archives": {"execution.tar.xz": "0" * 64, "validation.tar.xz": "0" * 64}}
    )
    if target == "manifest":
        manifest["forged"] = True
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    receipt = signed(
        {
            "schema": "island.j42-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "manifest_sha256": sha(tmp_path / "manifest.json"),
            "acceptance_gate": 1,
            "native_pcff_acceptance_gate": 1,
            "status": "asset_unavailable",
            "execution_profiles": [],
            "backend_executions": [],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    if target == "receipt":
        receipt["forged"] = True
    (tmp_path / "receipt.json").write_text(json.dumps(receipt))
    (tmp_path / "execution.tar.xz").write_bytes(b"tampered")
    with pytest.raises(
        ValueError,
        match={
            "receipt": "Receipt checksum changed",
            "manifest": "Manifest checksum changed",
            "archive": "Archive changed",
        }[target],
    ):
        validate(tmp_path)
