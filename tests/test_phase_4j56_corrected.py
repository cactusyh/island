"""Corrected J56 software contracts. Sentinel arrays are not QM calculations."""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import (
    CompletionError,
    classify,
    digest,
    inventory,
    load_assessment,
    pin,
    signed,
    write,
)
from island.pcff_completion.confirmation import (
    STATES,
    reported_confirmation,
    validate_primary,
    validate_reported,
)
from island.pcff_completion.qm import ADAPTERS, asset_status, deform, templates
from island.pcff_completion.readiness import readiness
from island.pcff_completion.references import reference_intake
from island.pcff_completion.targets import SETTINGS, QMTarget, configure, targets

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/evidence/phase_4j54/audit"
DIST = (
    ROOT.parent / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
)


def load(name="hmds"):
    if not (DIST / "pcff.frc").is_file():
        pytest.skip("pinned source unavailable")
    return load_assessment(AUDIT / name, pin(AUDIT / name / "manifest.json"), DIST)


@pytest.fixture(scope="module")
def hmds():
    return load()


def reference(root, case):
    """Files mimic a producer attestation solely to test integrity validation."""
    s, a = case
    p = a.payload
    t = targets(s, a, "hmds")[0]
    xyz = [s.coordinates.get(i).tolist() for i in p["atom_order"]]
    samples = []
    for i in range(3):
        x = copy.deepcopy(xyz)
        x[0][0] += i * 0.01
        samples.append(
            {
                "id": str(i),
                "graph_identity": p["graph_identity"],
                "bond_order_identity": p["graph_identity"],
                "atom_order": p["atom_order"],
                "coordinates": x,
                "energy": float(i),
                "charge_state": 0,
                "method_id": "sentinel",
                "split": "held_out" if i == 2 else "train",
            }
        )
    q = {
        "schema": "island.pcff-qm-reference.v1",
        "family": "pcff",
        "target_identity": t.payload["identity"],
        **{
            k: p[k]
            for k in (
                "graph_identity",
                "typing_identity",
                "charge_identity",
                "atom_order",
                "charge_vector",
                "source_hashes",
                "domain",
            )
        },
        "units": t.payload["units"],
        "samples_file": "samples.json",
        "license_file": "LICENSE",
        "permission": "licensed",
        "methods": {
            "sentinel": {
                "method": "software-only",
                "basis": "none",
                "dispersion": "none",
                "software": "sentinel",
                "version": "test",
                "convergence": "test",
                "multiplicity": 1,
            }
        },
        "output_verification_file": "verification.json",
    }
    (root / "LICENSE").write_text("Software test only; no scientific reference data")
    (root / "importer.py").write_text("# software sentinel only\n")
    records = []
    for sample in samples:
        name = sample["id"] + ".out"
        (root / name).write_text("Not a QM output: software integrity test")
        records.append(
            {
                "sample_id": sample["id"],
                "output_file": name,
                "output_sha256": pin(root / name),
                "normalized_sample_checksum": digest(sample),
                "status": "verified",
            }
        )
    write(root / "samples.json", samples)
    write(
        root / "verification.json",
        {
            "importer_file": "importer.py",
            "importer_sha256": pin(root / "importer.py"),
            "implementation_version": "test",
            "records": records,
        },
    )
    q["files"] = inventory(root)
    return t, q, pin(root / "verification.json")


def test_corrected_28_and_19_no_author_confirmation_blockers():
    queries, plans = {}, set()
    for name in ("disiloxane", "hmds", "dp2", "dp4", "dp6"):
        s, a = load(name)
        c = classify(s, a)
        for row in c["terms"]:
            assert row["status"] in STATES
            queries[row["query"]] = row["status"]
        assert all(
            r["status"] == "qm_fit_required" for r in c["missing_parameter_blockers"]
        )
        assert all(
            r["status"] == "author_confirmed" for r in c["effective_execution_contract"]
        )
        for t in targets(s, a, name):
            assert t.validate(s, a) == t.payload
            plans.add(t.payload["query"])
    assert list(queries.values()).count("author_confirmed") == 28
    assert list(queries.values()).count("qm_fit_required") == 19
    assert len(plans) == 19


@pytest.mark.parametrize(
    "field",
    [
        "artifact",
        "artifact_sha256",
        "recorded_date",
        "scope",
        "affected_interaction_families",
        "status",
        "numerical_values",
        "execution_authority",
        "author_confirmation_date",
    ],
)
def test_reported_confirmation_tampering_rejects_even_resigned(field):
    r = reported_confirmation()
    r[field] = "changed"
    r = signed({k: v for k, v in r.items() if k != "identity"})
    with pytest.raises(CompletionError, match="author_confirmation_mutation"):
        validate_reported(r)


def test_no_primary_author_artifact_or_numeric_confirmation_invented():
    r = validate_reported(reported_confirmation())
    assert r["status"] == "author_confirmed"
    assert r["artifact_kind"] == "user_reported_author_confirmation"
    assert r["primary_artifact_status"] == "not_supplied"
    assert r["author_confirmation_date"] is None
    assert r["numerical_values"] is None and not r["execution_authority"]


@pytest.mark.parametrize(
    "field",
    [
        "graph_identity",
        "typing_identity",
        "charge_identity",
        "source_hashes",
        "query",
        "functional_form",
        "atom_order",
        "charge_vector",
        "domain",
    ],
)
def test_target_mutation_cannot_be_hidden_by_checksum(hmds, field):
    s, a = hmds
    p = targets(s, a, "hmds")[0].payload
    p[field] = "changed"
    p = signed({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(CompletionError, match="qm_target_binding_mutation"):
        QMTarget(json.dumps(p)).validate(s, a)


def test_no_unrequested_QM_defaults(hmds):
    s, a = hmds
    t = targets(s, a, "hmds")[0]
    assert all(t.payload[k] is None for k in SETTINGS)
    with pytest.raises(CompletionError, match="explicit_qm_settings_required"):
        configure(t, {})


@pytest.mark.parametrize(
    "mutation",
    [
        "graph",
        "type",
        "charge",
        "source",
        "family",
        "order",
        "units",
        "license",
        "hash",
        "output",
        "missing_graph",
        "mixed_charge",
        "mixed_method",
        "forces",
        "nan",
        "sample_order",
        "duplicate",
        "unverified",
        "raw_replace",
        "importer_replace",
        "charge_state",
        "bond_order",
        "no_method",
    ],
)
def test_reference_negative_controls(hmds, tmp_path, mutation):
    s, a = hmds
    t, q, vsha = reference(tmp_path, hmds)
    samples = json.loads((tmp_path / "samples.json").read_text())
    if mutation in ("graph", "type", "charge", "source", "family", "order", "units"):
        key = {
            "graph": "graph_identity",
            "type": "typing_identity",
            "charge": "charge_vector",
            "source": "source_hashes",
            "family": "family",
            "order": "atom_order",
            "units": "units",
        }[mutation]
        q[key] = "changed"
    elif mutation == "license":
        q["permission"] = "unknown"
    elif mutation == "hash":
        q["files"]["samples.json"] = "0" * 64
    elif mutation == "output":
        q["output_verification_file"] = None
    elif mutation == "raw_replace":
        (tmp_path / "0.out").write_text("replacement")
    elif mutation == "importer_replace":
        (tmp_path / "importer.py").write_text("replacement")
    elif mutation == "mixed_method":
        q["methods"]["other"] = q["methods"]["sentinel"].copy()
    elif mutation == "no_method":
        q["methods"]["sentinel"]["basis"] = ""
    elif mutation == "unverified":
        vsha = "0" * 64
    else:
        if mutation == "missing_graph":
            samples[0].pop("graph_identity")
        elif mutation == "mixed_charge":
            samples[0]["charge_state"] = 1
        elif mutation == "charge_state":
            for row in samples:
                row["charge_state"] = 1
        elif mutation == "forces":
            samples[0]["forces"] = [[1, 2, 3]]
        elif mutation == "nan":
            # write() deliberately rejects NaN; the runtime must reject foreign JSON too.
            samples[0]["energy"] = float("nan")
        elif mutation == "sample_order":
            samples[0]["atom_order"] = samples[0]["atom_order"][::-1]
        elif mutation == "duplicate":
            samples[1]["coordinates"] = samples[0]["coordinates"]
        elif mutation == "bond_order":
            samples[0]["bond_order_identity"] = "wrong"
        (tmp_path / "samples.json").write_text(json.dumps(samples))
        q["files"] = inventory(tmp_path)
    with pytest.raises(CompletionError):
        reference_intake(s, a, t, q, tmp_path, verification_sha256=vsha)


def test_reference_integrity_is_not_a_fitted_term(hmds, tmp_path):
    s, a = hmds
    t, q, sha = reference(tmp_path, hmds)
    r = reference_intake(s, a, t, q, tmp_path, verification_sha256=sha)
    assert r["status"] == "reference_integrity_passed"
    assert r["parameter_completeness"] == "incomplete" and not r["execution_authority"]
    moved = tmp_path.parent / (tmp_path.name + "-moved")
    shutil.copytree(tmp_path, moved)
    assert reference_intake(s, a, t, q, moved, verification_sha256=sha) == r


@pytest.mark.parametrize("engine", ADAPTERS)
def test_adapter_missing_and_templates_never_claim_QM(hmds, tmp_path, engine):
    s, a = hmds
    t = targets(s, a, "hmds")[0]
    assert (
        asset_status(engine, "/no/such/quantum-executable")["status"]
        == "asset_unavailable"
    )
    settings = {
        "qm_method": "unexecuted test",
        "basis_set": "unexecuted test",
        "dispersion_treatment": "explicit test choice",
        "convergence_criteria": {"test": True},
        "deformation_range": {"units": "angstrom", "values": [-0.01, 0.01]},
        "expected_uncertainty": {"status": "unknown"},
        "held_out_thresholds": {"test": True},
    }
    bond = next(iter(s.topology.bonds.values()))
    scan = {
        "kind": "bond",
        "sites": [bond.site1, bond.site2],
        "units": "angstrom",
        "displacements": [-0.01, 0.01],
    }
    r = templates(
        s,
        a,
        t,
        settings,
        scan,
        tmp_path / "inputs",
        engine=engine,
        executable="/no/such/quantum-executable",
    )
    assert not r["QM_executed"] and r["reference_data_status"] == "not_supplied"
    assert r["adapter"]["status"] == "asset_unavailable"
    assert len(list((tmp_path / "inputs").glob("*.xyz"))) == 2
    assert "<" in next((tmp_path / "inputs").glob("*.template")).read_text()


@pytest.mark.parametrize("kind", ["bond", "angle", "torsion"])
def test_graph_deformations_preserve_order_and_finite_coordinates(hmds, kind):
    s, a = hmds
    row = next(
        r
        for r in a.payload["interaction_inventory"]
        if r["family"]
        == {"bond": "quartic_bond", "angle": "quartic_angle", "torsion": "torsion_3"}[
            kind
        ]
    )
    scan = {
        "kind": kind,
        "sites": row["sites"],
        "units": "angstrom" if kind == "bond" else "radian",
        "displacements": [-0.02, 0.0, 0.02],
    }
    geometries = deform(s, a, scan)
    assert all(g["atom_order"] == a.payload["atom_order"] for g in geometries)
    original = np.array([s.coordinates.get(i) for i in a.payload["atom_order"]])
    assert np.allclose(geometries[1]["coordinates"], original, rtol=0, atol=1e-14)
    assert np.isfinite(np.array([g["coordinates"] for g in geometries])).all()
    assert geometries[0]["geometry_checksum"] != geometries[2]["geometry_checksum"]


def test_readiness_relocation_fresh_process(hmds, tmp_path):
    s, a = hmds
    shutil.copytree(AUDIT / "hmds", tmp_path / "bundle")
    shutil.copytree(DIST, tmp_path / "sources")
    expected = readiness(s, a, "hmds")
    for seed in (1, 99397):
        out = tmp_path / f"{seed}.json"
        run = subprocess.run(
            [
                sys.executable,
                "-m",
                "island.pcff_completion.readiness",
                "readiness",
                "--bundle",
                str(tmp_path / "bundle"),
                "--manifest-sha256",
                pin(tmp_path / "bundle/manifest.json"),
                "--distribution",
                str(tmp_path / "sources"),
                "--fixture",
                "hmds",
                "--json",
                str(out),
            ],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            capture_output=True,
            text=True,
            check=False,
        )
        assert run.returncode == 2, run.stdout + run.stderr
        assert json.loads(out.read_text()) == expected
    assert expected["cross_engine_equivalence"] == "not_claimed"
    assert (
        expected["parameter_completeness"] == "incomplete"
        and not expected["MD_authorized"]
    )


def test_primary_confirmation_evidence_mutation(hmds, tmp_path):
    s, a = hmds
    queries = sorted(
        {r["query"] for r in classify(s, a)["effective_execution_contract"]}
    )
    (tmp_path / "author.txt").write_text(
        "Software test confirmation only; not an author message"
    )
    entries = []
    for q in queries:
        entry = {
            "query": q,
            "status": "author_confirmed",
            "interpreted_value": "test only",
            "units": "test",
            "location": "line 1",
            "raw_evidence": "Software test confirmation only",
            "source_file": "author.txt",
            "source_sha256": pin(tmp_path / "author.txt"),
        }
        entries.append({**entry, "evidence_checksum": digest(entry)})
    r = {
        "family": "pcff",
        "author": "software sentinel",
        "scope": "exact test query set",
        "date": "2026-10-30",
        "source_hashes": a.payload["source_hashes"],
        "files": inventory(tmp_path),
        "entries": entries,
    }
    assert (
        validate_primary(r, tmp_path, digest(r), queries)["status"]
        == "author_confirmed"
    )
    (tmp_path / "author.txt").write_text("replacement")
    with pytest.raises(CompletionError, match="confirmation_inventory_mutation"):
        validate_primary(r, tmp_path, digest(r), queries)


def test_fixture_name_cannot_change_graph_domain(hmds):
    s, a = hmds
    with pytest.raises(CompletionError, match="fixture_graph_mismatch"):
        targets(s, a, "dp6")


def test_reference_fresh_process_reconstruction(hmds, tmp_path):
    s, a = hmds
    data = tmp_path / "payload"
    data.mkdir()
    target, request, sha = reference(data, hmds)
    expected = reference_intake(s, a, target, request, data, verification_sha256=sha)
    moved = tmp_path / "moved"
    shutil.copytree(data, moved)
    request_path = tmp_path / "request.json"
    write(request_path, request)
    for seed in (1, 99397):
        out = tmp_path / f"reference-{seed}.json"
        run = subprocess.run(
            [
                sys.executable,
                "-m",
                "island.pcff_completion.readiness",
                "reference-intake",
                "--bundle",
                str(AUDIT / "hmds"),
                "--manifest-sha256",
                pin(AUDIT / "hmds/manifest.json"),
                "--distribution",
                str(DIST),
                "--fixture",
                "hmds",
                "--query",
                target.payload["query"],
                "--request",
                str(request_path),
                "--request-sha256",
                pin(request_path),
                "--data-root",
                str(moved),
                "--verification-sha256",
                sha,
                "--json",
                str(out),
            ],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
        )
        assert run.returncode == 0, run.stdout + run.stderr
        assert json.loads(out.read_text()) == expected


def test_family_only_corrected_planning_api(hmds):
    from island.pcff_completion.readiness import prepare_polymer_forcefield
    from island.pcff_pdms import use_pcff_pdms_assets

    s, a = hmds
    with use_pcff_pdms_assets(a.assets):
        result = prepare_polymer_forcefield(s, forcefield="pcff")
    p = result.validate_integrity(s)
    assert p["parameter_completeness"] == "incomplete"
    assert all(
        r["status"] == "qm_fit_required" for r in p["missing_parameter_blockers"]
    )
    assert all(
        r["status"] == "author_confirmed" for r in p["effective_execution_contract"]
    )
    assert p["cross_engine_equivalence"] == "not_claimed"
    with pytest.raises(CompletionError, match="wrong_profile_family"):
        prepare_polymer_forcefield(s, forcefield="opls")


def test_explicit_equivalence_claim_is_independent(hmds, tmp_path):
    out = tmp_path / "claims.json"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "island.pcff_completion.readiness",
            "readiness",
            "--bundle",
            str(AUDIT / "hmds"),
            "--manifest-sha256",
            pin(AUDIT / "hmds/manifest.json"),
            "--distribution",
            str(DIST),
            "--fixture",
            "hmds",
            "--require-cross-engine-equivalence",
            "--json",
            str(out),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 2, run.stdout + run.stderr
    p = json.loads(out.read_text())
    assert p["claim_failures"] == [
        {
            "claim": "cross_engine_equivalence",
            "field": "cross_engine_equivalence",
            "required": "passed",
            "observed": "not_claimed",
        }
    ]
    assert p["parameter_completeness"] == "incomplete"
