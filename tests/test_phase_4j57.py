"""J57 contracts and real geometry checks. Sentinel outputs are software tests only."""

import copy
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from island.pcff_completion import (
    CompletionError,
    canonical,
    digest,
    inventory,
    pin,
    signed,
    write,
)
from island.pcff_completion.qm import ADAPTERS, asset_status
from island.pcff_qm import prepare_polymer_forcefield
from island.pcff_qm.__main__ import load_cases, main
from island.pcff_qm.geometry import geometries, measure
from island.pcff_qm.inputs import deck, generate, validate_bundle
from island.pcff_qm.readiness import closure
from island.pcff_qm.records import (
    AUDIT_SCOPE,
    READINESS_STATES,
    Record,
    attestations,
    confirmation_inventory,
    original_artifact,
    record,
    validate_attestation,
)
from island.pcff_qm.reference import CONVERSION, fit_readiness, intake
from island.pcff_qm.supporting import supporting_scans
from island.pcff_qm.targets import BINDINGS, FIXTURES, targets, validate_target
from island.pcff_qm.targets import inventory as target_inventory

ROOT = Path(__file__).resolve().parents[1]
DIST = (
    ROOT.parent / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files"
)


@pytest.fixture(scope="module")
def cases():
    if not (DIST / "pcff.frc").is_file():
        pytest.skip("pinned PCFF source unavailable")
    return load_cases(ROOT, DIST)


@pytest.fixture(scope="module")
def chosen(cases):
    return target_inventory(cases)


def restore(value):
    return Record(canonical(value))


def test_28_audited_attestations_preserve_evidence_distinctions():
    rs = attestations()
    assert len(rs) == len({r.payload["query"] for r in rs}) == 28
    for r in rs:
        p = validate_attestation(r.payload)
        assert p["status"] == "author_confirmed_audited"
        assert p["authority"] == "user_reported_author_confirmation"
        assert p["evidence_kind"] == "user_attested_author_confirmation"
        assert p["source_artifact_status"] == "not_yet_supplied"
        assert p["audit_scope"] == AUDIT_SCOPE
        for invented in (
            "author_email",
            "author_date",
            "filename",
            "sha256",
            "numerical_values",
            "exact_author_text",
        ):
            assert invented not in p
    c = confirmation_inventory()
    assert (
        c["original_author_artifacts"]
        == c["independently_hashed_author_evidence"]
        == []
    )


@pytest.mark.parametrize(
    "field",
    [
        "status",
        "authority",
        "source_artifact_status",
        "audit_scope",
        "query",
        "source_hashes",
        "evidence_checksum",
    ],
)
def test_attestation_mutation_even_resigned(field):
    p = attestations()[0].payload
    p[field] = "mutated"
    p = signed({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(CompletionError):
        validate_attestation(p)


def test_additive_original_artifact_requires_supplied_evidence(tmp_path):
    # Deliberately synthetic document, never recorded in the phase evidence.
    a = attestations()[0].payload
    (tmp_path / "test.txt").write_text(
        "SYNTHETIC software fixture, no author confirmation"
    )
    req = {
        "attestation_identity": a["identity"],
        "files": inventory(tmp_path),
        "filename": "test.txt",
        "sha256": pin(tmp_path / "test.txt"),
        "date": "2026-01-01",
        "author_source_identity": "synthetic test",
        "location": "line 1",
        "exact_text": "SYNTHETIC software fixture",
        "permission_status": "licensed",
        "independent_reconstruction": {
            "status": "passed",
            "method": "software fixture",
            "implementation_version": "test",
            "result": "fixture",
        },
    }
    before = a.copy()
    added = original_artifact(a, req, tmp_path, expected_sha256=digest(req)).payload
    assert added["parent_attestation_identity"] == a["identity"]
    assert a == before and a["source_artifact_status"] == "not_yet_supplied"
    req["sha256"] = "0" * 64
    with pytest.raises(CompletionError, match="artifact_hash"):
        original_artifact(a, req, tmp_path, expected_sha256=digest(req))


def test_exact_19_target_inventory_and_family_counts(chosen):
    p = chosen.payload
    assert len(p["targets"]) == 19
    assert Counter(r["target"]["interaction_type"] for r in p["targets"]) == {
        "angle-angle": 12,
        "bond-bond_1_3": 7,
    }
    assert p["joint_fits_authorized"] == []
    for r in p["targets"]:
        t = r["target"]
        assert t["status"] == "qm_fit_required" and t["readiness"] == "qm_data_missing"
        assert t["qm_protocol"]["method_adequacy"].startswith("not_established")
        assert len(t["atom_order"]) == len(t["atom_labels"])
        assert t["scan"]["joint_fit"] is False
        assert (
            t["hessian_frequency_requirement"]["frequencies_alone_sufficient"] is False
        )


@pytest.mark.parametrize("fixture", FIXTURES)
def test_confirmation_contract_and_independent_statuses(cases, fixture):
    s, a = cases[fixture]
    p = closure(s, a, fixture).payload
    assert all(
        r["status"] == "author_confirmed_audited"
        for r in p["effective_execution_contract"]
    )
    assert all(
        r["status"] == "qm_fit_required" for r in p["missing_parameter_blockers"]
    )
    assert not p["original_artifact_absence_blocks_conventions"]
    assert p["parameter_completeness"] == "incomplete"
    assert p["cross_engine_equivalence"] == "not_claimed"
    assert p["typing_status"] == p["charge_status"] == "complete"
    assert p["registered_profiles"] == p["fitted_parameters"] == []
    assert p["production_validated"] is False and p["MD_authorized"] is False


def test_family_only_api_and_no_profile_registration(cases):
    s, a = cases["hmds"]
    r = prepare_polymer_forcefield(s, forcefield="pcff", assets=a.assets)
    p = r.validate_integrity(s)
    assert p == closure(s, a, "hmds").payload
    with pytest.raises(CompletionError, match="wrong_profile_family"):
        prepare_polymer_forcefield(s, forcefield="opls", assets=a.assets)


@pytest.mark.parametrize("index", range(19))
def test_all_target_geometries_are_measured_not_assumed(cases, chosen, index):
    t = chosen.payload["targets"][index]["target"]
    s, a = cases[t["molecular_fixture"]]
    before = [s.coordinates.get(i).tolist() for i in t["atom_order"]]
    gs = geometries(s, a, t["scan"])
    assert len(gs) == 25 and len({g["geometry_checksum"] for g in gs}) == 25
    assert {g["split"] for g in gs} == {"train", "held_out"}
    for g in gs:
        measured = [
            measure(g["coordinates"], t["atom_order"], c)
            for c in t["scan"]["coordinates"]
        ]
        assert (
            np.max(
                np.abs(
                    np.array(measured)
                    - np.array(g["baseline_coordinates"])
                    - g["offsets"]
                )
            )
            < 1e-10
        )
        assert "energy" not in g and "forces" not in g
    assert before == [s.coordinates.get(i).tolist() for i in t["atom_order"]]


@pytest.mark.parametrize("fixture", FIXTURES)
def test_supporting_scan_families_and_joint_coordinates(cases, fixture):
    s, a = cases[fixture]
    scans = supporting_scans(s, a, fixture)
    labels = {r.payload["label"] for r in scans}
    assert {
        "Si-O-bond",
        "Si-O-Si-angle",
        "joint-bond-bond",
        "joint-bond-angle",
        "joint-angle-angle",
    } <= labels
    if fixture != "disiloxane":
        assert {
            "Si-C-bond",
            "joint-angle-torsion",
            "joint-angle-angle-torsion",
        } <= labels
    for r in scans:
        g = geometries(s, a, r.payload["scan"])
        assert all(max(abs(v) for v in x["coordinate_residuals"]) < 1e-10 for x in g)
        assert r.payload["creates_missing_target"] is False


@pytest.mark.parametrize("engine", ADAPTERS)
def test_adapter_templates_and_missing_executable(cases, tmp_path, engine):
    s, a = cases["hmds"]
    t = targets(s, a, "hmds")[0]
    out = tmp_path / engine
    m = generate(s, a, t, out, engines=(engine,))
    assert validate_bundle(out, pin(out / "manifest.json")) == m
    assert len(list(out.glob(f"point-*.{engine}.*"))) == 25
    assert m["QM_executed"] is False
    assert (
        asset_status(engine, "/nonexistent/j57-engine")["status"] == "asset_unavailable"
    )
    central = deck(engine, "H 0 0 0", hessian=True)
    assert "0 1" in central or "spin=0" in central or "MULTIPLICITY 1" in central
    (out / "point-0000.xyz").write_text("mutated")
    with pytest.raises(CompletionError, match="inventory_mutation"):
        validate_bundle(out, pin(out / "manifest.json"))


@pytest.mark.parametrize(
    "field",
    [
        "graph_identity",
        "atom_order",
        "charge_vector",
        "source_hashes",
        "scan",
        "functional_form",
        "spin_state",
        "qm_protocol",
        "convention_contract",
        "input_geometry_checksum",
    ],
)
def test_target_mutations_resigned_fail_reconstruction(cases, field):
    s, a = cases["hmds"]
    p = targets(s, a, "hmds")[0].payload
    p[field] = "mutated"
    with pytest.raises(CompletionError, match="target_binding_mutation"):
        validate_target(s, a, record({k: v for k, v in p.items() if k != "identity"}))


def reference_fixture(root, case):
    """Sufficient shape/coverage sentinel, not QM or scientific reference evidence."""
    s, a = case
    t = targets(s, a, "hmds")[0]
    p = t.payload
    n = len(p["atom_order"])
    samples = []
    proofs = []
    for i, g in enumerate(geometries(s, a, p["scan"])):
        sample = {
            k: g[k] for k in ("id", "atom_order", "coordinates", "offsets", "split")
        }
        sample.update(
            {
                "graph_identity": p["graph_identity"],
                "graph": p["graph"],
                "energy": float(i),
                "forces": np.ones((n, 3)).tolist(),
                "charge_state": 0,
                "spin_state": p["spin_state"],
                "method_id": "test",
            }
        )
        if all(v == 0 for v in g["offsets"]):
            sample["hessian"] = np.eye(3 * n).tolist()
        samples.append(sample)
        name = g["id"] + ".out"
        (root / name).write_text("SOFTWARE SENTINEL ONLY; not an actual QM result\n")
        proofs.append(
            {
                "sample_id": g["id"],
                "output_files": {name: pin(root / name)},
                "normalized_sample_checksum": digest(sample),
                "status": "verified_by_supplied_importer",
                "raw_observables": {
                    "units": {
                        "energy": "hartree",
                        "gradient": "hartree/bohr",
                        "hessian": "hartree/bohr^2",
                    },
                    "energy": sample["energy"] / CONVERSION["energy_factor"],
                    "gradient": (
                        -np.array(sample["forces"])
                        * CONVERSION["distance_factor"]
                        / CONVERSION["energy_factor"]
                    ).tolist(),
                    **(
                        {
                            "hessian": (
                                np.array(sample["hessian"])
                                * CONVERSION["distance_factor"] ** 2
                                / CONVERSION["energy_factor"]
                            ).tolist()
                        }
                        if "hessian" in sample
                        else {}
                    ),
                },
            }
        )
    (root / "LICENSE").write_text("Software test fixture only; no scientific data")
    (root / "importer.py").write_text("# software verification sentinel\n")
    write(root / "samples.json", samples)
    request = {
        "schema": "island.pcff-qm-reference.v2",
        "family": "pcff",
        "target_identity": p["identity"],
        **{k: p[k] for k in BINDINGS},
        "graph": p["graph"],
        "functional_form_identity": p["functional_form"]["identity"],
        "units": p["units"],
        "samples_file": "samples.json",
        "verification_file": "verification.json",
        "license_file": "LICENSE",
        "permission_status": "licensed",
        "methods": {
            "test": {
                **{
                    k: p["qm_protocol"][k]
                    for k in ("method", "basis_set", "dispersion_treatment")
                },
                "software": "software sentinel",
                "version": "test",
                "calculation_settings": p["qm_protocol"],
            }
        },
        "conversion_source": CONVERSION["source"],
        "energy_factor": CONVERSION["energy_factor"],
        "distance_factor": CONVERSION["distance_factor"],
    }
    write(
        root / "verification.json",
        {
            "importer_file": "importer.py",
            "importer_sha256": pin(root / "importer.py"),
            "implementation_version": "test",
            "records": proofs,
            "normalization": {
                "energy": "hartree to kcal/mol by declared factor",
                "forces": "minus gradient; hartree/bohr to kcal/mol/angstrom",
                **{
                    k: request[k]
                    for k in ("conversion_source", "energy_factor", "distance_factor")
                },
            },
        },
    )
    request["files"] = inventory(root)
    return t, request, pin(root / "verification.json")


def test_strict_reference_integrity_not_fit_or_qm_claim(cases, tmp_path):
    t, r, sha = reference_fixture(tmp_path, cases["hmds"])
    d = intake(*cases["hmds"], t, r, tmp_path, verification_sha256=sha)
    assert d.payload["status"] == "qm_data_available"
    assert d.payload["parameter_completeness"] == "incomplete"
    assert fit_readiness(t).payload["status"] == "qm_data_missing"
    assert fit_readiness(t, d).payload["status"] == "qm_data_available"
    assert set(READINESS_STATES) == {
        "qm_data_missing",
        "qm_data_available",
        "fit_ready",
        "fit_completed",
        "fit_validated",
        "fit_rejected",
    }


@pytest.mark.parametrize(
    "field",
    [
        "graph_identity",
        "graph",
        "atom_order",
        "charge_vector",
        "charge_identity",
        "source_hashes",
        "functional_form_identity",
        "units",
        "target_identity",
        "permission_status",
        "family",
    ],
)
def test_reference_binding_mutations(cases, tmp_path, field):
    t, r, sha = reference_fixture(tmp_path, cases["hmds"])
    r[field] = "mutated"
    with pytest.raises(CompletionError):
        intake(*cases["hmds"], t, r, tmp_path, verification_sha256=sha)


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_graph",
        "atom_order",
        "charge",
        "spin",
        "force",
        "nonfinite",
        "energy",
        "hessian",
        "duplicate",
        "coverage",
        "units",
        "method",
        "output",
        "license",
        "settings",
        "form",
    ],
)
def test_reference_sample_and_raw_output_rejections(cases, tmp_path, mutation):
    t, r, sha = reference_fixture(tmp_path, cases["hmds"])
    samples = json.loads((tmp_path / "samples.json").read_text())
    if mutation == "missing_graph":
        samples[0].pop("graph_identity")
    if mutation == "atom_order":
        samples[0]["atom_order"].reverse()
    if mutation == "charge":
        samples[0]["charge_state"] = 1
    if mutation == "spin":
        samples[0]["spin_state"] = {"multiplicity": 3, "two_S": 2}
    if mutation == "force":
        samples[0]["forces"].pop()
    if mutation == "nonfinite":
        samples[0]["energy"] = float("nan")
    if mutation == "energy":
        samples[0].pop("energy")
    if mutation == "hessian":
        next(s for s in samples if "hessian" in s).pop("hessian")
    if mutation == "duplicate":
        samples[1]["coordinates"] = samples[0]["coordinates"]
    if mutation == "coverage":
        samples.pop()
    if mutation == "units":
        r.pop("units")
    if mutation == "method":
        r["methods"]["second"] = copy.deepcopy(r["methods"]["test"])
    if mutation == "output":
        (tmp_path / "point-0000.out").write_text("replaced")
    if mutation == "license":
        r["permission_status"] = "unknown"
    if mutation == "settings":
        r["methods"]["test"]["calculation_settings"] = {}
    if mutation == "form":
        r["functional_form_identity"] = "wrong"
    # Resign the inventory/proofs so each control exercises semantics, not just
    # the outer checksum. No untrusted producer field bypasses scientific checks.
    (tmp_path / "samples.json").write_text(json.dumps(samples))
    proof = json.loads((tmp_path / "verification.json").read_text())
    for s in samples:
        for row in proof["records"]:
            if row["sample_id"] == s["id"]:
                row["normalized_sample_checksum"] = (
                    digest(s) if mutation != "nonfinite" else "nonfinite"
                )
    write(tmp_path / "verification.json", proof)
    r["files"] = inventory(tmp_path)
    sha = pin(tmp_path / "verification.json")
    with pytest.raises(CompletionError):
        intake(*cases["hmds"], t, r, tmp_path, verification_sha256=sha)


def test_cli_independent_claim_requirement_and_attestations(tmp_path):
    out = tmp_path / "report.json"
    assert main(["confirmations", "--json", str(out)]) == 0
    assert len(json.loads(out.read_text())["confirmations"]) == 28
    assert (
        main(
            ["confirmations", "--require-cross-engine-equivalence", "--json", str(out)]
        )
        == 2
    )
    assert main(["confirmations", "--forcefield", "opls", "--json", str(out)]) == 2


def test_relocation_fresh_process_and_hash_seed(cases, tmp_path):
    # Reconstruct using a relocated repository audit and source directory.
    repo = tmp_path / "relocated"
    dest = repo / "docs/evidence/phase_4j54/audit"
    for name in FIXTURES:
        shutil.copytree(ROOT / "docs/evidence/phase_4j54/audit" / name, dest / name)
    shutil.copytree(DIST, tmp_path / "sources")
    from scripts.phase_4j57 import reconstruction

    expected = reconstruction(cases)
    script = "from island.pcff_qm.__main__ import load_cases; from scripts.phase_4j57 import reconstruction; from pathlib import Path; import sys,json; print(json.dumps(reconstruction(load_cases(Path(sys.argv[1]),Path(sys.argv[2])))))"
    for seed in (1, 99397):
        run = subprocess.run(
            [sys.executable, "-c", script, str(repo), str(tmp_path / "sources")],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            capture_output=True,
            text=True,
            check=False,
        )
        assert run.returncode == 0, run.stderr
        assert json.loads(run.stdout) == expected


def test_bad_fit_controls_cannot_assert_ready(cases, tmp_path):
    t, request, sha = reference_fixture(tmp_path, cases["hmds"])
    data = intake(*cases["hmds"], t, request, tmp_path, verification_sha256=sha)
    controls = dict.fromkeys(
        (
            "objective",
            "weights",
            "constraints",
            "bounds",
            "initial_values",
            "optimizer",
            "identifiable_source_bound_residual",
            "uncertainty_plan",
            "held_out_thresholds",
            "independent_reconstruction_plan",
        ),
        "asserted",
    )
    with pytest.raises(CompletionError, match="dataset_replay_required"):
        fit_readiness(t, data, controls)
    with pytest.raises(CompletionError, match="fit_objective_not_supported"):
        fit_readiness(
            t,
            data,
            controls,
            system=cases["hmds"][0],
            assessment=cases["hmds"][1],
            root=tmp_path,
            verification_sha256=sha,
        )


@pytest.mark.parametrize("mutation", ["unit_factor", "sign", "raw_energy", "raw_units"])
def test_independent_raw_conversion_reconstruction(cases, tmp_path, mutation):
    t, request, _sha = reference_fixture(tmp_path, cases["hmds"])
    proof = json.loads((tmp_path / "verification.json").read_text())
    if mutation == "unit_factor":
        request["energy_factor"] *= 4.184
        proof["normalization"]["energy_factor"] = request["energy_factor"]
    if mutation == "sign":
        raw = proof["records"][0]["raw_observables"]
        raw["gradient"] = (-np.array(raw["gradient"])).tolist()
    if mutation == "raw_energy":
        proof["records"][0]["raw_observables"]["energy"] += 1
    if mutation == "raw_units":
        proof["records"][0]["raw_observables"]["units"] = {}
    write(tmp_path / "verification.json", proof)
    request["files"] = inventory(tmp_path)
    with pytest.raises(CompletionError):
        intake(
            *cases["hmds"],
            t,
            request,
            tmp_path,
            verification_sha256=pin(tmp_path / "verification.json"),
        )
