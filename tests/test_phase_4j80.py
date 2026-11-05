"""Synthetic gate controls never establish source authority or scientific data."""

import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from island.pcff_pdms.pilot_intake import gates
from scripts import phase_4j80_pilot as j80
from scripts import validate_phase_4j80_receipt as validator

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / j80.ROOT


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j80.scope()


@pytest.fixture(scope="module")
def docs():
    return j80.reconstruct(ROOT)


def fixture_contract():
    """Manufactured byte/audit controls; no paper or force-field facts."""
    raw = b"%PDF-1.4\nsynthetic software fixture, not a scientific source\n"
    artifact = {
        "id": "synthetic-artifact",
        "sha256": hashlib.sha256(raw).hexdigest(),
        "kind": "primary_PDF",
        "authority": "original_publisher",
        "authority_audit_id": "independently-pinned-synthetic-audit",
        "immutable_filename": "synthetic.pdf",
        "source_version": "fixture-v1",
        "publication_year": 1997,
        "license_or_access_status": "synthetic_fixture",
        "citation": {
            "title": "fixture",
            "authors": ["fixture"],
            "identifier": "synthetic",
        },
    }
    row = {
        "candidate_id": "synthetic-candidate",
        "artifact_id": artifact["id"],
        "source_sha256": artifact["sha256"],
        "location": "synthetic/table-1/row-1",
        "reporting_mode": "directly_reported",
        "source_family": "published_PDMS",
        "raw_values_and_units": {"value": "synthetic", "unit": "synthetic"},
    }
    fields = {
        k: {
            "value": dict.fromkeys(names, "synthetic_documented_value"),
            "evidence_refs": ["synthetic/table-1"],
        }
        for k, names in gates.CONVENTIONS.items()
    }
    fields["vdw_function"]["value"] = {
        "source_family": "lj_9_6",
        "target_family": "lj_9_6",
        "source_equation": gates.VDW_9_6,
        "target_equation": gates.VDW_9_6,
    }
    fields["sigma_epsilon_minimum"]["value"]["sigma_meaning"] = (
        "minimum_energy_distance"
    )
    fields["charges_and_neutrality"]["value"].update(
        repeat_charge_e=0,
        capped_chain_charge_e=0,
        atom_type_charge_mapping={"synthetic_si_type": "0.0"},
        type_counts_by_scope={
            name: {"synthetic_si_type": n}
            for name, n in (("repeat", 1), ("hmds", 2), ("dp2", 3), ("dp4", 5))
        },
    )
    fields["energy_and_length_units"]["value"] = {
        "energy": "kJ/mol",
        "length": "angstrom",
        "angle_deviation": "radian",
    }
    fields["cross_terms"]["value"]["complete_required_coverage"] = True
    fields["internal_and_end_group_typing"]["value"].update(
        source_type_ids=["synthetic_si_type"], cap_transfer_supported=True
    )
    fields["finite_and_periodic_scope"]["value"].update(
        finite_capped=True, periodic_support=False
    )
    manifest = {"artifacts": [artifact], "row_evidence": [row]}
    audit = {
        artifact["authority_audit_id"]: {
            "artifact_metadata": copy.deepcopy(artifact),
            "verified_rows": {row["candidate_id"]: copy.deepcopy(row)},
            "convention_assertions": copy.deepcopy(fields),
            "evidence_kind": "synthetic_contract_control",
        }
    }
    return manifest, {artifact["id"]: raw}, audit, {"fields": fields}


def run_fixture(manifest, contents, audits, contract):
    a = gates.source_gate(manifest, contents, audits, ["synthetic-candidate"])
    return a, gates.convention_gate(contract, a)


def test_real_run_stops_at_source_gate(docs):
    report = docs["gate-report.json"]
    assert report["A"]["passed"] is False
    assert "primary_artifact_unavailable" in report["A"]["blocking_reasons"]
    assert report["B"]["status"] == report["C"]["status"] == "not_evaluated"
    assert report["scientific_work_stopped_at"] == "Gate_A"
    assert docs["result.json"]["decision"] == "qm_pilot_blocked"


def test_complete_manifest_preserves_unverified_claims(docs):
    m = docs["source-provenance-manifest.json"]
    assert m["artifacts"] == []
    assert m["authority_audits_available"] == 0
    assert len(m["context_files"]) == 27
    assert len(m["row_evidence"]) == 20
    assert all(len(f["sha256"]) == 64 for f in m["context_files"])
    assert all(
        r["source_sha256"] is None and r["location"] is None for r in m["row_evidence"]
    )
    assert all(
        r["reporting_mode"] == "unverified_user_reported" for r in m["row_evidence"]
    )
    assert m["publisher_access"]["full_text_recovered"] is False


def test_synthetic_complete_gates_do_not_launch_or_accept():
    a, b = run_fixture(*fixture_contract())
    assert a["passed"] and b["passed"]
    c = gates.pilot_gate(a, b)
    assert not c["launch_allowed"]
    assert "explicit_QM_execution_profile_missing" in c["blocking_reasons"]
    assert not gates.pilot_gate(a, b, {"claimed_authorization": True})["launch_allowed"]


@pytest.mark.parametrize(
    "field",
    (
        "sha256",
        "citation",
        "publication_year",
        "source_version",
        "immutable_filename",
        "license_or_access_status",
    ),
)
def test_missing_source_provenance_rejected(field):
    m, contents, audits, contract = fixture_contract()
    del m["artifacts"][0][field]
    assert not run_fixture(m, contents, audits, contract)[0]["passed"]


def test_mutated_primary_bytes_rejected():
    m, contents, audits, contract = fixture_contract()
    contents["synthetic-artifact"] += b"mutation"
    a, _ = run_fixture(m, contents, audits, contract)
    assert not a["passed"] and any("hash" in r for r in a["blocking_reasons"])


def test_self_attestation_and_publisher_metadata_cannot_supply_authority():
    m, contents, _, contract = fixture_contract()
    m["artifacts"][0]["author_confirmed_reported"] = True
    assert not run_fixture(m, contents, {}, contract)[0]["passed"]
    m, contents, audits, contract = fixture_contract()
    m["artifacts"][0]["kind"] = "publisher_metadata"
    audits["independently-pinned-synthetic-audit"]["artifact_metadata"] = copy.deepcopy(
        m["artifacts"][0]
    )
    assert not run_fixture(m, contents, audits, contract)[0]["passed"]


@pytest.mark.parametrize(
    "mutation", ("location", "raw_values_and_units", "source_sha256", "reporting_mode")
)
def test_resigned_row_transcription_and_location_rejected(mutation):
    m, contents, audits, contract = fixture_contract()
    m["row_evidence"][0][mutation] = "fabricated"
    assert not run_fixture(m, contents, audits, contract)[0]["passed"]


@pytest.mark.parametrize("mode", ("inferred", "reconstructed", "analog"))
def test_non_reported_values_blocked_even_with_bytes(mode):
    m, contents, audits, contract = fixture_contract()
    m["row_evidence"][0]["reporting_mode"] = mode
    audits["independently-pinned-synthetic-audit"]["verified_rows"][
        "synthetic-candidate"
    ] = copy.deepcopy(m["row_evidence"][0])
    a, b = run_fixture(m, contents, audits, contract)
    assert not a["passed"]
    assert (
        gates.classify(m["row_evidence"][0], a, b) == "inferred_or_analog_value_blocked"
    )


@pytest.mark.parametrize(
    "mutation",
    ("row_duplicate", "row_missing", "artifact_duplicate", "byte_identity_duplicate"),
)
def test_duplicate_and_missing_identity_coverage_rejected(mutation):
    m, contents, audits, contract = fixture_contract()
    if mutation == "row_duplicate":
        m["row_evidence"].append(copy.deepcopy(m["row_evidence"][0]))
    elif mutation == "row_missing":
        m["row_evidence"].clear()
    else:
        art = copy.deepcopy(m["artifacts"][0])
        if mutation == "byte_identity_duplicate":
            art["id"] = "copy"
            contents["copy"] = contents["synthetic-artifact"]
        m["artifacts"].append(art)
    assert not run_fixture(m, contents, audits, contract)[0]["passed"]


@pytest.mark.parametrize("field", tuple(gates.CONVENTIONS))
def test_every_missing_convention_blocks(field):
    m, contents, audits, contract = fixture_contract()
    contract["fields"][field]["value"] = None
    a, b = run_fixture(m, contents, audits, contract)
    assert a["passed"] and not b["passed"]
    assert "unresolved_convention:" + field in b["blocking_reasons"]


@pytest.mark.parametrize(
    "field", ("energy_and_length_units", "torsion", "combination_rules", "cross_terms")
)
def test_convention_mutation_cannot_be_self_resigned(field):
    m, contents, audits, contract = fixture_contract()
    first = next(iter(contract["fields"][field]["value"]))
    contract["fields"][field]["value"][first] = "new_unverified_convention"
    assert not run_fixture(m, contents, audits, contract)[1]["passed"]


@pytest.mark.parametrize(
    "mutation",
    (
        "12_6",
        "sigma",
        "cap_charge",
        "charge_arithmetic",
        "element_types",
        "cross_terms",
        "unit",
    ),
)
def test_even_audited_incompatible_conventions_fail(mutation):
    m, contents, audits, contract = fixture_contract()
    f = contract["fields"]
    if mutation == "12_6":
        f["vdw_function"]["value"]["target_family"] = "lj_12_6"
    elif mutation == "sigma":
        f["sigma_epsilon_minimum"]["value"]["sigma_meaning"] = "zero_crossing_distance"
    elif mutation == "cap_charge":
        f["charges_and_neutrality"]["value"]["capped_chain_charge_e"] = 0.240
    elif mutation == "charge_arithmetic":
        f["charges_and_neutrality"]["value"]["atom_type_charge_mapping"][
            "synthetic_si_type"
        ] = ".240"
    elif mutation == "unit":
        f["energy_and_length_units"]["value"]["energy"] = "guessed_energy_unit"
    elif mutation == "element_types":
        f["internal_and_end_group_typing"]["value"]["source_type_ids"] = ["Si", "O"]
    else:
        f["cross_terms"]["value"]["complete_required_coverage"] = False
    audits["independently-pinned-synthetic-audit"]["convention_assertions"] = (
        copy.deepcopy(f)
    )
    assert not run_fixture(m, contents, audits, contract)[1]["passed"]


def test_literature_verified_is_not_native_and_caps_need_separate_authority():
    m, contents, audits, contract = fixture_contract()
    a, b = run_fixture(m, contents, audits, contract)
    row = m["row_evidence"][0]
    assert gates.classify(row, a, b) == "literature_candidate_verified"
    row["end_group"] = True
    row["audited_end_group_scope"] = True  # untrusted claim cannot grant scope
    assert gates.classify(row, a, b) == "end_group_scope_blocked"


@pytest.mark.parametrize(
    "family", ("GAFF", "OPLS", "CVFF", "COMPASS", "CFF91", "analog")
)
def test_cross_family_promotion_rejected(family):
    m, contents, audits, contract = fixture_contract()
    a, b = run_fixture(m, contents, audits, contract)
    row = {**m["row_evidence"][0], "source_family": family}
    assert gates.classify(row, a, b) == "inferred_or_analog_value_blocked"


def test_native_identity_cannot_be_invented_after_source_audit():
    m, contents, audits, contract = fixture_contract()
    a, b = run_fixture(m, contents, audits, contract)
    row = {
        **m["row_evidence"][0],
        "source_family": "native_PCFF",
        "native_identity": {"source_row_id": "invented"},
    }
    assert gates.classify(row, a, b) == "inferred_or_analog_value_blocked"


def test_synthetic_native_identity_needs_every_source_and_convention_gate():
    m, contents, audits, contract = fixture_contract()
    row = m["row_evidence"][0]
    row["source_family"] = "native_PCFF"
    row["native_identity"] = {
        "upstream_release": "synthetic-release",
        "source_row_id": "synthetic-row",
        "ordered_source_types": ["synthetic_si_type", "synthetic_o_type"],
        "exact_environment": {"graph_identity": "synthetic-graph"},
        "equation": "synthetic-source-equation",
        "units_and_phase": {"value": "synthetic-source-contract"},
        "ordering": "synthetic-source-order",
        "exclusion_and_mixing_policy": {"value": "synthetic-source-policy"},
        "provenance": {"artifact": "synthetic-artifact"},
    }
    audits["independently-pinned-synthetic-audit"]["verified_rows"][
        "synthetic-candidate"
    ] = copy.deepcopy(row)
    a, b = run_fixture(m, contents, audits, contract)
    assert gates.classify(row, a, b) == "exact_native_row"
    assert (
        gates.classify(row, a, {"passed": False})
        == "exact_native_row_convention_blocked"
    )
    assert j80.SCOPE["runtime_registration"] is False
    assert j80.SCOPE["native_exact_rows_resolved"] == 0


def split_samples():
    return [
        {
            "id": "train",
            "role": "training",
            "geometry_identity": "g1",
            "conformer_group": "c1",
            "geometry_family": "bond",
            "graph_identity": "hmds",
        },
        {
            "id": "valid",
            "role": "validation",
            "geometry_identity": "g2",
            "conformer_group": "c2",
            "geometry_family": "coupled",
            "graph_identity": "dp2",
        },
        {
            "id": "test",
            "role": "test",
            "geometry_identity": "g3",
            "conformer_group": "c3",
            "geometry_family": "equilibrium",
            "graph_identity": "dp4",
        },
    ]


@pytest.mark.parametrize(
    "mutation",
    ("none", "id", "geometry_identity", "conformer_group", "family", "held_out"),
)
def test_training_validation_and_topology_leakage(mutation):
    samples = split_samples()
    if mutation == "none":
        assert gates.validate_split(samples)
        return
    if mutation in ("id", "geometry_identity", "conformer_group"):
        samples[1][mutation] = samples[0][mutation]
    elif mutation == "family":
        samples[1]["geometry_family"] = samples[0]["geometry_family"]
    else:
        samples[2]["graph_identity"] = samples[0]["graph_identity"]
    with pytest.raises(ValueError):
        gates.validate_split(samples)


@pytest.mark.parametrize(
    "kind",
    (
        "identified",
        "rank_deficient",
        "under_determined",
        "ill_conditioned",
        "low_signal",
        "bounds",
    ),
)
def test_identifiability_before_any_fitting(kind):
    import numpy as np

    a = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [-1.0, 1.0]])
    bounds = [{"parameter": p, "lower": -1.0, "upper": 1.0} for p in ("a", "b")]
    if kind == "rank_deficient":
        a[:, 1] = a[:, 0]
    elif kind == "under_determined":
        a = a[:1]
    elif kind == "ill_conditioned":
        a[:, 1] = a[:, 0] + 1e-8 * a[:, 1]
    elif kind == "low_signal":
        a *= 1e-16
    elif kind == "bounds":
        bounds[0]["lower"] = None
        with pytest.raises(ValueError):
            gates.synthetic_diagnostics(a, ["a", "b"], bounds)
        return
    r = gates.synthetic_diagnostics(a, ["a", "b"], bounds)
    assert r["identifiable"] == (kind == "identified")
    assert r["null_space_dimension"] == (
        1 if kind in ("rank_deficient", "under_determined") else 0
    )
    assert r["accepted_parameters"] == r["fitted_parameters"] == []
    assert r["scientific_QM_data"] is False


def test_no_QM_requests_or_scientific_design_generated(docs):
    q = docs["qm-pilot-request-manifest.json"]
    assert q["requests"] == q["completed_requests"] == []
    assert q["QM_input_decks_generated"] == 0
    assert q["execution_enabled"] is False
    assert [m["name"] for m in q["prospective_contract_only"]["molecules"]] == [
        "hmds",
        "dp2",
        "dp4",
    ]
    assert q["prospective_contract_only"]["held_out_topology"] == "dp4"
    assert q["prospective_contract_only"]["full_423_campaign_allowed"] is False
    c = docs["identifiability-and-conditioning.json"]
    assert not c["computed"]
    assert len(c["blocks"]) == 6
    assert all(
        b["rank"] is None and b["null_space_dimension"] is None for b in c["blocks"]
    )


def test_complete_candidate_and_target_classification(docs):
    r = docs["candidate-row-classification.json"]
    assert len(r["candidate_rows"]) == 20
    assert r["candidate_class_counts"]["literature_candidate_convention_blocked"] == 20
    assert len(r["extension_targets"]) == 129
    assert r["unresolved_cross_terms"] == 87
    assert r["unresolved_end_group_terms"] == 69
    assert {t["classification"] for t in r["extension_targets"]} <= set(
        gates.ROW_CLASSES
    )
    assert all(not t["resolved"] for t in r["extension_targets"])


@pytest.mark.parametrize(
    "name",
    (
        "source-provenance-manifest.json",
        "convention-contract.json",
        "candidate-row-classification.json",
        "qm-pilot-request-manifest.json",
        "identifiability-and-conditioning.json",
        "result.json",
    ),
)
def test_resigned_semantic_report_mutation_rejected(tmp_path, name):
    root = tmp_path / "evidence"
    shutil.copytree(ROOT, root)
    d = j80.read(root / name)
    d["fabricated_scientific_pass"] = True
    d.pop("identity")
    j80.write(root / name, j80.signed(d))
    with pytest.raises(ValueError, match="semantic_mutation"):
        j80.replay(root)


@pytest.mark.parametrize(
    "name",
    (
        "source-catalog.json",
        "trusted-authority-audits.json",
        "source-access.json",
        "parent-inputs.json",
    ),
)
def test_frozen_source_and_authority_mutation_rejected(tmp_path, name):
    root = tmp_path / "evidence"
    shutil.copytree(ROOT, root)
    p = root / name
    p.write_bytes(p.read_bytes() + b" ")
    with pytest.raises(ValueError, match="source_hash_mutation"):
        j80.reconstruct(root)


@pytest.mark.parametrize(
    "field",
    (
        "completed_QM_requests",
        "accepted_parameters",
        "runtime_registration",
        "automatic_assignment_enabled",
        "production_validated",
        "J80_QM_attempts",
    ),
)
def test_resigned_receipt_activation_rejected(field):
    r = {
        "schema": "island.j80-receipt.v1",
        "base": j80.BASE,
        "branch": j80.BRANCH,
        "decision": "qm_pilot_blocked",
        "scope": dict(j80.SCOPE),
    }
    m = {
        "schema": "island.j80-manifest.v1",
        "base": j80.BASE,
        "branch": j80.BRANCH,
        "freeze_sha256": j80.FREEZE_SHA,
        "scope": dict(j80.SCOPE),
    }
    r["scope"][field] = True if type(r["scope"][field]) is bool else 1
    with pytest.raises(ValueError):
        validator.semantic(j80.signed(r), j80.signed(m))


@pytest.mark.parametrize(
    "name", ("sqm", "psi4", "lammps", "msi2lmp", "synthetic_fit_driver")
)
def test_forbidden_fake_executables_never_launch(tmp_path, name):
    marker = tmp_path / "marker"
    executable = tmp_path / name
    executable.write_text(
        f"#!{sys.executable}\nfrom pathlib import Path\nPath({str(marker)!r}).touch()\n"
    )
    executable.chmod(0o755)
    with pytest.raises((PermissionError, RuntimeError)):
        subprocess.run([str(executable)], check=True)
    assert not marker.exists()


def test_child_inherits_fence_and_cannot_import_backend(tmp_path):
    marker = tmp_path / "imported"
    (tmp_path / "psi4.py").write_text(
        f"from pathlib import Path\nPath({str(marker)!r}).touch()\n"
    )
    code = (
        f"import sys; sys.path.insert(0,{str(tmp_path)!r}); "
        "from scripts.phase_4j76_source_resolution import scope; scope(); "
        "import importlib.util; assert importlib.util.find_spec('psi4'); import psi4"
    )
    p = subprocess.run(
        [sys.executable, "-c", code],
        env=os.environ.copy(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode != 0
    assert "scope" in p.stderr.lower()
    assert not marker.exists()


def test_relocated_fresh_process_replay(tmp_path):
    root = tmp_path / "moved"
    shutil.copytree(ROOT, root)
    p = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j80_pilot",
            "replay",
            "--directory",
            str(root),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode == 0, p.stderr + p.stdout
    assert json.loads(p.stdout) == j80.replay(ROOT)


def test_J79_and_all_history_preserved():
    assert validator.history(REPO, ROOT) > 10322


@pytest.mark.parametrize(
    "control",
    (
        "QM",
        "MD",
        "converter",
        "physical_fitting",
        "evaluator",
        "historical_reserved_validation",
    ),
)
def test_scientific_controls_not_launched(control):
    pytest.skip("J80 Gate A blocked: " + control + " is not launched")
