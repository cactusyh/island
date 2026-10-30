"""J56 contracts and fail-closed controls; no reference coefficients or MD."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from island.forcefields.selection.inventory import canonical
from island.pcff_completion import (
    CLASSES,
    SUPPORTED,
    CompletionError,
    CompletionRecord,
    classify,
    code_identity,
    derive,
    equation,
    fit,
    intake,
    inventory,
    load_assessment,
    pin,
    request_template,
    signed,
    write,
)
from island.pcff_pdms import REGISTERED_PROFILES

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/evidence/phase_4j54/audit"
DIST = (
    ROOT.parent / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
)


def load(name="hmds", distribution=DIST, root=AUDIT):
    if not (distribution / "pcff.frc").is_file():
        pytest.skip("pinned source distribution unavailable")
    return load_assessment(
        root / name, pin(root / name / "manifest.json"), distribution
    )


@pytest.fixture(scope="module")
def hmds():
    return load()


def dummy_contract(tmp_path, hmds, family="angle-angle"):
    """Intake-only software sentinel data; never a physical fitting dataset."""
    s, assessment = hmds
    p = assessment.payload
    row = next(r for r in p["missing_terms"] if r["family"] == family)
    q = request_template(s, assessment, row["id"])
    xyz = [s.coordinates.get(i).tolist() for i in p["atom_order"]]
    data = [
        {
            "id": str(i),
            "atom_order": p["atom_order"],
            "coordinates": xyz,
            "energy": float(i + 1),
            "split": "train" if i < 8 else "held_out",
        }
        for i in range(10)
    ]
    write(tmp_path / "samples.json", data)
    (tmp_path / "LICENSE.txt").write_text(
        "Software contract sentinel only, no scientific reference data"
    )
    q.update(samples_file="samples.json", files=inventory(tmp_path))
    q["provenance"] = {
        "kind": "reference_engine",
        "citation": "J56 software contract only",
        "license": (tmp_path / "LICENSE.txt").read_text(),
        "license_file": "LICENSE.txt",
        "permission_status": "licensed",
        "data_origin": "nonphysical test sentinel",
        "reference_version": "software-only-v1",
        "method": "not a QM calculation",
        "basis": "not a basis",
        "software": "test data producer",
    }
    n = len(q["fit"]["fields"])
    q["fit"].update(
        initial_values=[0.125] * n,
        lower_bounds=[-100.0] * n,
        upper_bounds=[100.0] * n,
        weights={s["id"]: 1.0 for s in data},
    )
    return q


def test_all_47_queries_and_2916_sites_classified_without_derivations():
    expected = json.loads((AUDIT / "missing-rows.json").read_text())
    queries = {}
    total = 0
    for name in ("disiloxane", "hmds", "dp2", "dp4", "dp6"):
        s, a = load(name)
        p = classify(s, a)
        total += len(p["terms"])
        assert p["parameter_completeness"] == "incomplete"
        assert p["created_terms"] == p["registered_profiles"] == []
        assert not p["MD_authorized"] and not p["production_validated"]
        assert p["cross_engine_equivalence"] == "not_claimed"
        for r in p["terms"]:
            assert r["classification"] in CLASSES
            assert r["source_status"] == r["status"]
            assert r["historical_source_status"] == "source_missing"
            assert not r["analytical_derivation_authorized"]
            queries.setdefault(r["query"], set()).add(r["classification"])
    assert total == 2916 and len(queries) == 47
    assert set(queries) == {q["query"] for q in expected["queries"]}
    assert sum(v == {"qm_fit_required"} for v in queries.values()) == 19
    assert sum(v == {"author_confirmed"} for v in queries.values()) == 28
    assert not REGISTERED_PROFILES


@pytest.mark.parametrize(
    "family",
    [
        "angle-angle",
        "angle-angle-torsion_1",
        "bond-bond_1_3",
        "angle-torsion_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
    ],
)
def test_no_independent_coefficient_is_derived_from_equilibria(hmds, family):
    s, a = hmds
    target = next(r["id"] for r in a.payload["missing_terms"] if r["family"] == family)
    p = derive(s, a, target)
    assert p["status"] == "blocked" and p["result"] is None
    assert p["origin_on_success"] == "derived_external"
    assert p["formula"] is None and p["uncertainty"]["value"] is None
    assert p["source_equation"]["coefficient_value_equation"] is None
    assert CompletionRecord(canonical(p)).validate_integrity(s, a) == p


def test_template_and_no_data_reject_fit(hmds, tmp_path):
    s, a = hmds
    target = a.payload["missing_terms"][0]["id"]
    q = request_template(s, a, target)
    assert all(x is None for x in q["fit"]["initial_values"])
    with pytest.raises(CompletionError, match="reference_data_unavailable"):
        fit(s, a, q, tmp_path)


@pytest.mark.parametrize(
    "family",
    [
        "angle-angle",
        "bond-bond_1_3",
        "angle-angle-torsion_1",
        "angle-torsion_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
    ],
)
def test_all_missing_families_accept_data_contract_but_block_unisolated_fit(
    hmds, tmp_path, family
):
    s, a = hmds
    q = dummy_contract(tmp_path, hmds, family)
    p = intake(s, a, q, tmp_path)
    assert p["status"] == "blocked"
    assert (
        "total_reference_requires_verified_joint_residual_decomposition"
        in p["blockers"]
    )
    assert p["created_terms"] == []
    if "torsion" in family:
        assert "author_confirmed_query_is_not_one_of_19_qm_targets" in p["blockers"]
    with pytest.raises(CompletionError, match="fitting_blocked"):
        fit(s, a, q, tmp_path)
    assert (
        CompletionRecord(canonical(p)).validate_integrity(s, a, data_root=tmp_path) == p
    )


def test_all_required_class2_record_families_have_explicit_equations(hmds):
    _, a = hmds
    for family in SUPPORTED:
        e = equation(a, family)
        assert e["source_sha256"] == a.payload["source_hashes"]["pcff.frc"]
        assert "E =" in e["raw_text"] and e["fields"] and e["source_units"]


@pytest.mark.parametrize(
    "mutation",
    [
        "family",
        "profile",
        "source",
        "graph",
        "typing",
        "charge_identity",
        "charges",
        "atom_order",
        "domain",
        "query",
        "form",
        "units",
        "license",
        "permission",
        "hash",
        "bounds",
        "weights",
        "optimizer",
        "method",
        "reference_version",
        "geometry",
        "force",
        "hessian",
        "split",
        "traversal",
    ],
)
def test_mutations_reject_before_fitting(hmds, tmp_path, mutation):
    s, a = hmds
    q = dummy_contract(tmp_path, hmds)
    if mutation == "family":
        q["family"] = "opls"
    elif mutation == "profile":
        q["candidate_profile"] = "island_external_opls_pdms_xiang2025_dft_local_v1"
    elif mutation == "source":
        q["source_hashes"]["pcff.frc"] = "0" * 64
    elif mutation in ("graph", "typing", "charge_identity"):
        q[
            {
                "graph": "graph_identity",
                "typing": "typing_identity",
                "charge_identity": "charge_identity",
            }[mutation]
        ] = "0" * 64
    elif mutation == "charges":
        q["charge_vector"]["1"] += 0.1
    elif mutation == "atom_order":
        q["atom_order"] = q["atom_order"][::-1]
    elif mutation == "domain":
        q["domain"]["kind"] = "Si-O-C"
    elif mutation == "query":
        q["query"] += ",wildcard"
    elif mutation == "form":
        q["equation_identity"] = "0" * 64
    elif mutation == "units":
        q["units"]["angles"] = "degree"
    elif mutation == "license":
        q["provenance"]["license"] = "invented grant"
    elif mutation == "permission":
        q["provenance"]["permission_status"] = "unknown"
    elif mutation == "hash":
        q["files"]["samples.json"] = "0" * 64
    elif mutation == "bounds":
        q["fit"]["lower_bounds"] = [999]
    elif mutation == "weights":
        q["fit"]["weights"]["0"] = -1
    elif mutation == "optimizer":
        q["fit"]["optimizer"] = "guess"
    elif mutation == "method":
        q["provenance"]["software"] = None
    elif mutation == "reference_version":
        q["provenance"]["reference_version"] = ""
    elif mutation == "traversal":
        q["samples_file"] = "../samples.json"
    else:
        data = json.loads((tmp_path / "samples.json").read_text())
        if mutation == "geometry":
            data[0]["coordinates"].pop()
        elif mutation == "force":
            data[0]["forces"] = [[1, 2, 3]]
        elif mutation == "hessian":
            data[0]["hessian"] = [[1]]
        elif mutation == "split":
            data[0]["split"] = "secret"
        write(tmp_path / "samples.json", data)
        q["files"] = inventory(tmp_path)
    with pytest.raises(CompletionError):
        fit(s, a, q, tmp_path)


def test_raw_force_hessian_vibration_experiment_channels_are_intake_only(
    hmds, tmp_path
):
    s, a = hmds
    q = dummy_contract(tmp_path, hmds)
    data = json.loads((tmp_path / "samples.json").read_text())
    n = len(q["atom_order"])
    data[0].update(
        forces=[[1.0, 2.0, 3.0]] * n,
        hessian=[[1.0] * (3 * n)] * (3 * n),
        vibrations={"units": "cm^-1", "frequencies": [25.0]},
        observables=[
            {
                "name": "mode",
                "units": "cm^-1",
                "method": "software sentinel",
                "value": 25.0,
                "uncertainty": 1.0,
            }
        ],
    )
    write(tmp_path / "samples.json", data)
    q["files"] = inventory(tmp_path)
    p = intake(s, a, q, tmp_path)
    assert p["status"] == "blocked" and not p["created_terms"]


def test_rank_deficient_component_reference_does_not_create_coefficient(hmds, tmp_path):
    s, a = hmds
    q = dummy_contract(tmp_path, hmds, "bond-bond_1_3")
    target = next(r for r in a.payload["missing_terms"] if r["id"] == q["target_id"])
    data = json.loads((tmp_path / "samples.json").read_text())
    i, j, _, _ = target["sites"]
    dep = next(
        r
        for r in a.payload["interaction_inventory"]
        if r["id"] == target["dependencies"][0]
    )
    r0 = dep["assignment"]["values"][0]
    for row in data:
        row["coordinates"][q["atom_order"].index(i)] = [0.0, 0.0, 0.0]
        row["coordinates"][q["atom_order"].index(j)] = [r0, 0.0, 0.0]
    write(tmp_path / "samples.json", data)
    write(
        tmp_path / "raw.json",
        {
            "query": q["query"],
            "target_id": q["target_id"],
            "graph_identity": q["graph_identity"],
            "samples": data,
        },
    )
    (tmp_path / "producer.txt").write_text("Software negative control only; no QM data")
    write(
        tmp_path / "isolation.json",
        {
            "query": q["query"],
            "target_id": q["target_id"],
            "graph_identity": q["graph_identity"],
            "method": "negative test",
            "code_file": "producer.txt",
            "raw_reference_file": "raw.json",
            "equation_identity": q["equation_identity"],
            "reference_command": "not executed",
        },
    )
    q.update(
        files=inventory(tmp_path),
        energy_scope="isolated_target_component",
        isolation_evidence_file="isolation.json",
    )
    assert intake(s, a, q, tmp_path)["blockers"] == []
    with pytest.raises(CompletionError, match="not_identifiable"):
        fit(s, a, q, tmp_path)
    raw = json.loads((tmp_path / "raw.json").read_text())
    raw["target_id"] = "another_site_with_the_same_type_query"
    write(tmp_path / "raw.json", raw)
    q["files"] = inventory(tmp_path)
    with pytest.raises(CompletionError, match="raw_reference_identity_mutation"):
        fit(s, a, q, tmp_path)


@pytest.mark.parametrize("mutation", ["result", "origin", "uncertainty", "code"])
def test_resigned_completion_mutation_rejected(hmds, mutation):
    s, a = hmds
    p = derive(s, a, a.payload["missing_terms"][0]["id"])
    if mutation == "result":
        p["result"] = {"value": 0}
    elif mutation == "origin":
        p["origin_on_success"] = "native_pcff"
    elif mutation == "uncertainty":
        p["uncertainty"] = {"status": "certain", "value": 0}
    else:
        p["code_identity"] = "0" * 64
    p = signed({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(CompletionError, match="completion_evidence_mutation"):
        CompletionRecord(canonical(p)).validate_integrity(s, a)


def test_live_graph_charge_source_mutation_rejected(hmds, tmp_path):
    s, a = hmds
    altered = s.copy()
    position = altered.coordinates.get(1).copy()
    position[0] += 1
    altered.coordinates.set(1, position)
    with pytest.raises(Exception, match="mutation"):
        classify(altered, a)
    shutil.copytree(DIST, tmp_path / "sources")
    (tmp_path / "sources/pcff.frc").write_text("replacement")
    with pytest.raises(Exception, match="source_hash_mutation"):
        load(distribution=tmp_path / "sources")


def test_recursive_inventory_and_symlinks(hmds, tmp_path):
    nested = tmp_path / "a"
    nested.mkdir()
    (nested / "data.txt").write_text("reference")
    assert inventory(tmp_path) == {"a/data.txt": pin(nested / "data.txt")}
    (tmp_path / "link").symlink_to(nested / "data.txt")
    with pytest.raises(CompletionError, match="symlink"):
        inventory(tmp_path)


def test_relocation_and_two_fresh_processes(hmds, tmp_path):
    s, a = hmds
    shutil.copytree(AUDIT / "hmds", tmp_path / "bundle")
    shutil.copytree(DIST, tmp_path / "sources")
    expected = classify(s, a)
    for seed in (1, 99397):
        output = tmp_path / f"classification-{seed}.json"
        cmd = [
            sys.executable,
            "-m",
            "island.pcff_completion",
            "classify",
            "--bundle",
            str(tmp_path / "bundle"),
            "--manifest-sha256",
            pin(tmp_path / "bundle/manifest.json"),
            "--distribution",
            str(tmp_path / "sources"),
            "--json",
            str(output),
        ]
        r = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            check=False,
        )
        assert r.returncode == 0, r.stdout + r.stderr
        assert json.loads(output.read_text()) == expected
    assert code_identity() == expected["code_identity"]


def test_cli_without_data_fails_closed(hmds, tmp_path):
    output = tmp_path / "result.json"
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "island.pcff_completion",
            "fit",
            "--bundle",
            str(AUDIT / "hmds"),
            "--manifest-sha256",
            pin(AUDIT / "hmds/manifest.json"),
            "--distribution",
            str(DIST),
            "--json",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 2 and "reference_data_unavailable" in r.stdout
    assert not output.exists()


def test_schema_and_template_do_not_grant_derivation_authority(hmds):
    from island.pcff_completion.schema import SCHEMA

    assert SCHEMA["properties"]["family"] == {"const": "pcff"}
    assert SCHEMA["$defs"]["derived_record"]["properties"]["origin"] == {
        "const": "derived_external"
    }
    assert SCHEMA["$defs"]["fitted_record"]["properties"]["origin"] == {
        "const": "fitted_external"
    }
    assert "uncertainty" in SCHEMA["$defs"]["fitted_record"]["required"]
    s, a = hmds
    with pytest.raises(CompletionError, match="verified_J54_assessment_required"):
        intake(s, a.payload, {}, None)
    with pytest.raises(CompletionError, match="exact_missing_target_required"):
        derive(s, a, "made-up:1,2,3")


@pytest.mark.parametrize(
    "defect", ["held_out", "samples", "force_only", "initial", "omitted_file"]
)
def test_insufficient_reference_cannot_create_term(hmds, tmp_path, defect):
    s, a = hmds
    q = dummy_contract(tmp_path, hmds)
    data = json.loads((tmp_path / "samples.json").read_text())
    if defect == "held_out":
        for row in data:
            row["split"] = "train"
    elif defect == "samples":
        data = data[:1]
    elif defect == "force_only":
        for row in data:
            row.pop("energy")
            row["forces"] = [[1.0, 2.0, 3.0]] * len(q["atom_order"])
    elif defect == "initial":
        q["fit"]["initial_values"] = [None]
    elif defect == "omitted_file":
        (tmp_path / "unlisted.txt").write_text("replacement")
    if defect in ("held_out", "samples", "force_only"):
        write(tmp_path / "samples.json", data)
        q["files"] = inventory(tmp_path)
        q["fit"]["weights"] = {s["id"]: 1.0 for s in data}
    with pytest.raises(CompletionError):
        fit(s, a, q, tmp_path)


def test_intake_record_relocation_preserves_provenance_and_blockers(hmds, tmp_path):
    s, a = hmds
    payload = tmp_path / "payload"
    payload.mkdir()
    q = dummy_contract(payload, hmds)
    p = intake(s, a, q, payload)
    relocated = tmp_path / "relocated-payload"
    shutil.copytree(payload, relocated)
    assert (
        CompletionRecord(canonical(p)).validate_integrity(s, a, data_root=relocated)
        == p
    )
    record, output = tmp_path / "record.json", tmp_path / "reconstructed.json"
    write(record, p)
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "island.pcff_completion",
            "validate",
            "--bundle",
            str(AUDIT / "hmds"),
            "--manifest-sha256",
            pin(AUDIT / "hmds/manifest.json"),
            "--distribution",
            str(DIST),
            "--record",
            str(record),
            "--record-sha256",
            pin(record),
            "--data-root",
            str(relocated),
            "--json",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    # Reconstruction preserves the blocked status, hence exit 2.
    assert run.returncode == 2 and output.is_file(), run.stdout + run.stderr
    assert json.loads(output.read_text()) == p
    (relocated / "LICENSE.txt").write_text("replacement")
    with pytest.raises(CompletionError, match="inventory_or_hash_mutation"):
        CompletionRecord(canonical(p)).validate_integrity(s, a, data_root=relocated)
