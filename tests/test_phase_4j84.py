"""Manufactured source bundles test software only, never scientific qualification."""

import json
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from island.forcefields.assignment_core import normalize, rebuild, type_graph
from island.forcefields.assignment_core.schema import CROSS_TERMS, FAMILIES, data
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.bundle_adapters import adapt, capture, ingest, proposals
from island.forcefields.bundle_adapters.schema import (
    ADAPTER,
    AUDIT,
    BUNDLE,
    SOURCE,
    manifest_identity,
    sha,
    source_digest,
)
from scripts import phase_4j84_bundles as j84
from tests.test_phase_4j82 import control, graph

REPO = Path(__file__).resolve().parents[1]
WHEN = "2026-11-06T00:00:00Z"
MOTIFS = j84.read(REPO / j84.FIXTURE)["motifs"]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j84.scope()


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")
    return sha(path.read_bytes())


def bundle(
    tmp_path, family="OPLS-AA", g=None, change=None, real_branch=False, split=False
):
    """Explicit fake values from the existing J82 software controls, not stock data."""
    g, p = control(g, family)
    records = []

    def row(kind, name, value):
        records.append(
            {
                "record_id": kind + ":" + name,
                "kind": kind,
                "location": "Synthetic test table; row " + name,
                "origin": "directly_reported",
                "value": value,
            }
        )

    for r in p.rules:
        row("atom_type", r.id, {k: v for k, v in data(r).items() if k != "evidence_id"})
    for f in p.forms:
        row("form", f.family, data(f))
    for c in p.conventions:
        row("convention", c.name, {"name": c.name, "value": c.value})
    for r in p.parameters:
        row(
            "parameter",
            r.id,
            {
                k: v
                for k, v in data(r).items()
                if k
                not in (
                    "profile_family",
                    "evidence_id",
                    "convention_identity",
                    "source_location",
                )
            },
        )
    row("graph_scope", "graph", {"graph_identity": g.identity})
    conventions = {c.name: c.value for c in p.conventions}
    m = {
        "schema": BUNDLE,
        "bundle_id": "software_" + family,
        "profile_id": "audited_control_" + family,
        "source_uri": "urn:island:synthetic-software-control:" + family,
        "source_sha256": "",
        "family": family,
        "version": "synthetic-v1",
        "files": [],
        "unit_convention": {
            n: conventions[n] for n in ("energy_unit", "length_unit", "angular_unit")
        },
        "functional_form_convention": data(p.forms),
        "mixing_rule": conventions["mixing_rule"],
        "policy_14": {
            "lj": conventions["scaling_14_lj"],
            "coulomb": conventions["scaling_14_coulomb"],
        },
        "charge_convention": conventions["charge_origin"],
        "provenance": {
            "citation": "Manufactured software fixture; not scientific data",
            "authors": "ISLAND test authors",
            "year": 2026,
            "persistent_id": "urn:synthetic:software-only",
            "license": "CC0 synthetic test data",
            "redistribution": "permitted",
            "acquisition_method": "generated software control",
            "authorization": "authorized_for_local_readonly_ingestion",
            "origin": "authorized_table" if real_branch else "synthetic_control",
            "upstream_identity": "manufactured test release",
        },
        "ingestion_timestamp": WHEN,
        "adapter_version": ADAPTER,
        "synthetic": not real_branch,
        "evidence_status": "verified_stock" if real_branch else "derived_provisional",
        "topology_support": list(p.topology_support),
        "identity": "",
    }
    if change:
        change(m, records)
    root = tmp_path / "bundle"
    root.mkdir(parents=True)
    parts = (
        (records[: len(records) // 2], records[len(records) // 2 :])
        if split
        else (records,)
    )
    for i, items in enumerate(parts):
        relative = f"sources/records_{i}.json"
        path = root / relative
        digest = put(
            path,
            {
                "schema": SOURCE,
                "family": family,
                "version": m["version"],
                "records": items,
            },
        )
        m["files"].append(
            {
                "path": relative,
                "sha256": digest,
                "size": path.stat().st_size,
                "role": "record_source",
                "family": family,
                "version": m["version"],
            }
        )
    primary = root / "sources/original.txt"
    primary.write_text("Synthetic original artifact; no physical parameters.\n")
    m["files"].append(
        {
            "path": "sources/original.txt",
            "sha256": sha(primary.read_bytes()),
            "size": primary.stat().st_size,
            "role": "primary_artifact",
            "family": family,
            "version": m["version"],
        }
    )
    m["source_sha256"] = source_digest(m["files"])
    m["identity"] = manifest_identity(m)
    manifest = root / "manifest.json"
    manifest_sha = put(manifest, m)
    cap = capture(manifest, manifest_sha, capture_timestamp=WHEN)
    _, proposed = proposals(cap)
    # This deliberately manufactured audit exercises the acceptance code path.
    # No production function creates, signs or approves independent audits.
    audit = {
        "schema": AUDIT,
        "manifest_identity": m["identity"],
        "reviewer": "synthetic test auditor",
        "review_timestamp": WHEN,
        "authority_basis": "Manufactured control only; not a scientific review",
        "policy": "independent_source_review",
        "synthetic": m["synthetic"],
        "source_reviews": [
            {
                "path": f["path"],
                "sha256": f["sha256"],
                "provenance_verified": True,
                "license_verified": True,
                "claims": [
                    {
                        "record_id": r["record"]["record_id"],
                        "claim_sha256": r["claim_sha256"],
                    }
                    for r in proposed
                    if r["file"] == f["path"]
                ],
            }
            for f in m["files"]
        ],
    }
    audit_path = tmp_path / "independent-audit.json"
    audit_sha = put(audit_path, audit)
    return g, (manifest, manifest_sha, audit_path, audit_sha)


def load(args):
    return ingest(*args, capture_timestamp=WHEN)


def test_missing_source_never_creates_scientific_data():
    result = ingest(capture_timestamp=WHEN)
    assert result.status == "source_bundle_missing"
    assert (
        result.profile is None
        and not result.trusted_evidence
        and not result.source_hashes
    )
    assert result.accepted_coefficients == 0 and not result.automatic_assignment_enabled
    assert result.registry().coverage(graph()).ranked_profiles == ()


@pytest.mark.parametrize("family", FAMILIES)
def test_all_family_adapters_map_complete_controls(tmp_path, family):
    g, args = bundle(tmp_path, family, split=True)
    result = load(args)
    assert result.status == "audited_bundle_candidate" and not result.blocking_reasons
    p = result.profile
    assert p.family == family and p.functional_class == (
        "II" if family in ("PCFF", "COMPASS") else "I"
    )
    report = result.registry().coverage(g)
    assert all(
        c.coverage_decision == "provisional" for c in report.capability_matrix.cells
    )
    assert not report.profile_reports[0].missing_terms
    assert set(result.source_hashes) == set(report.source_hashes)
    assert len(result.source_hashes) == 3
    evidence = json.loads(result.evidence_graph_json)
    assert all(
        r["pointer"].startswith("/records/")
        and r["record"]["location"]
        and r["record"]["value"]
        for r in evidence["records"]
    )
    assert not evidence["convention_conversions"]
    assert not report.automatic_assignment_enabled and not report.production_validated
    assert report.simulation_readiness == "not_established"


@pytest.mark.parametrize("motif", MOTIFS, ids=lambda m: m["id"])
def test_each_organic_motif_can_use_source_bound_adapter(tmp_path, motif):
    g = normalize(j84.spec_from_data(motif))
    _, args = bundle(tmp_path, g=g)
    report = load(args).registry().coverage(g)
    assert report.profile_reports[0].status == "synthetic_complete"
    assert not report.profile_reports[0].missing_terms
    assert report.source_hashes


def test_manufactured_audit_branch_does_not_enable_runtime(tmp_path):
    g, args = bundle(tmp_path, real_branch=True)
    report = load(args).registry().coverage(g)
    assert all(
        c.coverage_decision == "accepted" for c in report.capability_matrix.cells
    )
    assert report.accepted_coefficients == 0 and not report.runtime_registration
    assert report.profile_reports[0].status == "coverage_complete_not_activated"


def test_no_audit_means_no_profile_even_with_verified_label(tmp_path):
    _, args = bundle(tmp_path, real_branch=True)
    result = ingest(args[0], args[1], capture_timestamp=WHEN)
    assert result.profile is None
    assert "independent_source_audit_missing" in result.blocking_reasons
    assert result.source_hashes and json.loads(result.evidence_graph_json)["records"]


def test_source_records_and_locations_are_preserved_exactly(tmp_path):
    _, args = bundle(tmp_path)
    cap = capture(*args, capture_timestamp=WHEN)
    result = adapt(cap)
    evidence = json.loads(result.evidence_graph_json)
    for file in cap["files"]:
        if file["raw_text"] is not None:
            assert (
                file["raw_text"].encode()
                == (args[0].parent / file["path"]).read_bytes()
            )
    for row in evidence["records"]:
        assert row["source_sha256"] in result.source_hashes
        if row["record"]["kind"] == "parameter":
            source = row["record"]["value"]
            mapped = next(p for p in result.profile.parameters if p.id == source["id"])
            assert list(mapped.values) == source["values"]
            assert (
                mapped.source_location
                == row["file"] + "#" + row["pointer"] + "; " + row["record"]["location"]
            )


@pytest.mark.parametrize("kind", ("manifest", "audit"))
def test_unpinned_bundle_or_audit_is_rejected(tmp_path, kind):
    _, args = bundle(tmp_path)
    values = list(args)
    values[1 if kind == "manifest" else 3] = None
    with pytest.raises(ValueError, match="external_hash_required"):
        load(values)


def test_missing_source_file_never_uses_a_checkout_copy(tmp_path):
    _, args = bundle(tmp_path)
    m = j84.read(args[0])
    (args[0].parent / m["files"][0]["path"]).unlink()
    with pytest.raises(ValueError, match="unavailable_or_symlink"):
        load(args)


@pytest.mark.parametrize("which", ("manifest", "source", "audit"))
def test_external_SHA_mutations_rejected(tmp_path, which):
    _, args = bundle(tmp_path)
    m = j84.read(args[0])
    path = (
        args[0]
        if which == "manifest"
        else args[2]
        if which == "audit"
        else args[0].parent / m["files"][0]["path"]
    )
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="SHA256_mutation"):
        load(args)


@pytest.mark.parametrize("which", ("manifest", "audit", "record"))
def test_captured_resigned_fabrication_rejected(tmp_path, which):
    _, args = bundle(tmp_path)
    cap = capture(*args, capture_timestamp=WHEN)
    if which == "manifest":
        cap["manifest"]["version"] = "fabricated"
    elif which == "audit":
        cap["audit"]["source_reviews"][0]["claims"] = []
    else:
        next(f for f in cap["files"] if f["raw_text"] is not None)["raw_text"] += " "
    with pytest.raises(ValueError, match="mutated"):
        adapt(cap)


def test_resigned_audit_cannot_hide_missing_claims(tmp_path):
    _, args = bundle(tmp_path)
    audit = j84.read(args[2])
    audit["source_reviews"][0]["claims"] = []
    new_pin = put(args[2], audit)
    with pytest.raises(ValueError, match="audit_claim_set_mismatch"):
        load((*args[:3], new_pin))


@pytest.mark.parametrize("field", ("provenance_verified", "license_verified"))
def test_source_authority_and_license_must_be_audited(tmp_path, field):
    _, args = bundle(tmp_path)
    audit = j84.read(args[2])
    audit["source_reviews"][0][field] = False
    new_pin = put(args[2], audit)
    with pytest.raises(ValueError, match="not_audited"):
        load((*args[:3], new_pin))


@pytest.mark.parametrize(
    "mutation",
    (
        "mixing",
        "policy_14",
        "charge",
        "units",
        "form",
        "missing_convention",
        "analog",
        "polyjarvis",
        "license",
        "authorization",
        "status",
    ),
)
def test_convention_and_authority_fail_closed_despite_audit(tmp_path, mutation):
    def change(m, rows):
        if mutation == "mixing":
            m["mixing_rule"] = "different_rule"
        elif mutation == "policy_14":
            m["policy_14"]["lj"] = "different_scaling"
        elif mutation == "charge":
            m["charge_convention"] = "unverified_charge_model"
        elif mutation == "units":
            m["unit_convention"]["length_unit"] = "changed_unit"
        elif mutation == "form":
            m["functional_form_convention"][0]["equation"] = "9_6_cannot_be_12_6"
        elif mutation == "missing_convention":
            rows[:] = [
                r
                for r in rows
                if not (
                    r["kind"] == "convention" and r["value"]["name"] == "exclusions_13"
                )
            ]
        elif mutation == "analog":
            m["provenance"]["origin"] = "analog"
        elif mutation == "polyjarvis":
            m["provenance"]["upstream_identity"] = "PolyJarvis numeric patch"
        elif mutation == "license":
            m["provenance"]["license"] = "unknown"
        elif mutation == "authorization":
            m["provenance"]["authorization"] = "unconfirmed"
        else:
            m["evidence_status"] = "local_patch"

    _, args = bundle(tmp_path, change=change)
    result = load(args)
    assert result.profile is None and result.blocking_reasons
    assert result.accepted_coefficients == 0


@pytest.mark.parametrize(
    "family",
    ("bond", "angle", "torsion", "improper", "nonbonded", "charge", *CROSS_TERMS),
)
def test_missing_term_never_zero_filled(tmp_path, family):
    def change(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (r["kind"] == "parameter" and r["value"]["family"] == family)
        ]

    g, args = bundle(tmp_path, "PCFF", graph("*CCC*"), change)
    result = load(args)
    assert result.profile is not None
    assert not any(p.family == family for p in result.profile.parameters)
    report = result.registry().coverage(g)
    assert any(m.family == family for m in report.profile_reports[0].missing_terms)
    assert all(c.coverage_decision == "blocked" for c in report.capability_matrix.cells)


def test_incomplete_numeric_record_not_filled(tmp_path):
    def change(_, rows):
        next(r for r in rows if r["kind"] == "parameter")["value"]["values"] = []

    with pytest.raises(ValueError, match="no_zero_fill"):
        bundle(tmp_path, change=change)


def test_cross_term_dependency_failure_remains_explicit(tmp_path):
    def change(_, rows):
        next(
            r
            for r in rows
            if r["kind"] == "parameter" and r["value"]["family"] == "bond_angle"
        )["value"]["equilibrium_references"] = []

    g, args = bundle(tmp_path, "PCFF", change=change)
    report = load(args).registry().coverage(g)
    assert any(
        m.reason == "equilibrium_reference_value_mismatch"
        for m in report.profile_reports[0].missing_terms
    )


@pytest.mark.parametrize("compatible", (True, False))
def test_9_6_nonbonded_is_preserved_and_never_coerced_to_12_6(tmp_path, compatible):
    equation = "epsilon * (2*(sigma/r)^9 - 3*(sigma/r)^6)"

    def change(m, rows):
        for f in m["functional_form_convention"]:
            if f["family"] == "nonbonded":
                f["equation"] = (
                    equation if compatible else "4*epsilon*((sigma/r)^12-(sigma/r)^6)"
                )
        for r in rows:
            if r["kind"] == "form" and r["value"]["family"] == "nonbonded":
                r["value"]["equation"] = equation
            elif r["kind"] == "parameter" and r["value"]["family"] == "nonbonded":
                r["value"]["form"]["equation"] = equation
                r["value"]["values"] = ["1.00000", "2.00000"]
            elif r["kind"] == "convention" and r["value"]["name"] == "nonbonded_form":
                r["value"]["value"] = equation

    _, args = bundle(tmp_path, "PCFF", change=change)
    result = load(args)
    if compatible:
        rows = [r for r in result.profile.parameters if r.family == "nonbonded"]
        assert rows and all(
            r.form.equation == equation and r.values == ("1.00000", "2.00000")
            for r in rows
        )
        assert not json.loads(result.evidence_graph_json)["convention_conversions"]
    else:
        assert result.profile is None
        assert "functional_form_units_or_ordering_mismatch" in result.blocking_reasons


def test_branch_and_crosslink_edits_require_retyping_and_audited_scope(tmp_path):
    from island.forcefields.assignment_core.schema import Branch, Crosslink

    g, args = bundle(tmp_path, g=graph("*CCC*", 2))
    result = load(args)
    typed = type_graph(g, result.profile, result.trusted_evidence)
    for changed in (
        rebuild(g, branches=(Branch(2, "*C"),)),
        rebuild(g, crosslinks=(Crosslink(2, 5),)),
    ):
        with pytest.raises(ValueError, match="topology_changed_retyping_required"):
            typed.validate(changed, result.profile)
        report = result.registry().coverage(changed)
        assert (
            "exact_graph_scope_not_evidenced"
            in report.profile_reports[0].blocking_reasons
        )
        assert report.profile_reports[0].missing_terms


def test_same_priority_atom_type_candidates_block_and_replay(tmp_path):
    def change(_, rows):
        old = next(
            r
            for r in rows
            if r["kind"] == "atom_type" and r["value"]["atom_type"] == "control_C"
        )
        new = json.loads(json.dumps(old))
        new["record_id"] = "atom_type:competing_C"
        new["value"]["id"] = "competing_C"
        new["value"]["atom_type"] = "other_C"
        rows.append(new)

    g, args = bundle(tmp_path, change=change)
    result = load(args)
    a, b = result.registry().coverage(g), load(args).registry().coverage(g)
    assert a == b
    assert any("ambiguous_site" in r for r in a.profile_reports[0].blocking_reasons)
    assert all(c.coverage_decision == "blocked" for c in a.capability_matrix.cells)


def test_explicit_typing_priority_is_deterministic(tmp_path):
    def change(_, rows):
        old = next(
            r
            for r in rows
            if r["kind"] == "atom_type" and r["value"]["atom_type"] == "control_C"
        )
        new = json.loads(json.dumps(old))
        new["record_id"] = "atom_type:low_priority_C"
        new["value"].update(id="low_priority_C", atom_type="other_C", priority=-1)
        rows.insert(0, new)

    g, args = bundle(tmp_path, change=change)
    outcome = load(args)
    typed = type_graph(g, outcome.profile, outcome.trusted_evidence)
    assert "other_C" not in dict(typed.assignments).values()
    assert (
        outcome.registry().coverage(g).profile_reports[0].status == "synthetic_complete"
    )


def test_repeat_cut_end_groups_and_copolymer_require_new_scope_and_typing(tmp_path):
    from rdkit import Chem

    original = graph("*CO*", 2)
    _, args = bundle(tmp_path, g=original)
    outcome = load(args)
    typed = type_graph(original, outcome.profile, outcome.trusted_evidence)
    reverse = graph("*OC*", 2)
    assert Chem.MolToSmiles(molecule(original)) == Chem.MolToSmiles(molecule(reverse))
    variants = [
        reverse,
        rebuild(original, end_groups=("*C", "*O")),
        rebuild(
            original, repeat_units=(("A", "*CO*"), ("B", "*CC*")), sequence=("A", "B")
        ),
    ]
    for g in variants:
        with pytest.raises(ValueError, match="topology_changed_retyping_required"):
            typed.validate(g, outcome.profile)
        fresh = type_graph(g, outcome.profile, outcome.trusted_evidence)
        assert fresh.graph_identity == g.identity
        report = outcome.registry().coverage(g)
        assert (
            "exact_graph_scope_not_evidenced"
            in report.profile_reports[0].blocking_reasons
        )
        assert all(
            c.coverage_decision == "blocked" for c in report.capability_matrix.cells
        )


def test_source_family_mixing_rejected(tmp_path):
    _, args = bundle(tmp_path, "GAFF2")
    m = j84.read(args[0])
    m["files"][0]["family"] = "GAFF"
    m["identity"] = manifest_identity(m)
    pin = put(args[0], m)
    with pytest.raises(ValueError, match="cross_family_or_version_overlay"):
        load((args[0], pin, *args[2:]))


@pytest.mark.parametrize(
    "mutation",
    (
        "license",
        "duplicate_file",
        "duplicate_hash",
        "path_escape",
        "symlink",
        "record_identity",
    ),
)
def test_source_contract_integrity(tmp_path, mutation):
    if mutation == "record_identity":

        def change(_, rows):
            rows.append(rows[0])

        with pytest.raises(ValueError, match="duplicate_source_record_identity"):
            bundle(tmp_path, change=change)
        return
    _, args = bundle(tmp_path)
    m = j84.read(args[0])
    if mutation == "license":
        del m["provenance"]["license"]
    elif mutation == "duplicate_file":
        m["files"].append(m["files"][0])
    elif mutation == "duplicate_hash":
        copy = dict(m["files"][0])
        copy["path"] = "sources/duplicate.json"
        m["files"].append(copy)
    elif mutation == "path_escape":
        m["files"][0]["path"] = "../secret.json"
    else:
        path = args[0].parent / m["files"][0]["path"]
        moved = tmp_path / "outside.json"
        path.rename(moved)
        path.symlink_to(moved)
    m["source_sha256"] = source_digest(m["files"])
    m["identity"] = manifest_identity(m)
    pin = put(args[0], m)
    with pytest.raises(ValueError):
        load((args[0], pin, *args[2:]))


def test_bundle_reingestion_after_directory_relocation(tmp_path):
    g, args = bundle(tmp_path / "original")
    first = load(args)
    root = tmp_path / "relocated"
    shutil.copytree(args[0].parent, root)
    audit = root.parent / "relocated-audit.json"
    shutil.copyfile(args[2], audit)
    second = load((root / args[0].name, args[1], audit, args[3]))
    assert data(first) == data(second)
    assert (
        first.registry().coverage(g).replay_identity
        == second.registry().coverage(g).replay_identity
    )
    later = ingest(
        root / args[0].name,
        args[1],
        audit,
        args[3],
        capture_timestamp="2026-11-06T01:00:00Z",
    )
    assert (
        first.bundle_identity == later.bundle_identity
        and first.profile.identity == later.profile.identity
    )


@pytest.mark.parametrize(
    "field,value",
    (
        ("automatic_assignment_enabled", True),
        ("production_validated", True),
        ("accepted_coefficients", 1),
        ("runtime_registration", True),
        ("simulation_readiness", "ready"),
    ),
)
def test_readiness_cannot_be_enabled(field, value):
    with pytest.raises(ValueError, match="init=False"):
        replace(ingest(capture_timestamp=WHEN), **{field: value})


@pytest.mark.parametrize(
    "name",
    (
        "sqm",
        "lammps",
        "msi2lmp",
        "fit_hmds",
        "reserved_validation",
        "generate_parameters",
    ),
)
def test_no_forbidden_child_execution(tmp_path, name):
    marker = tmp_path / "executed"
    fake = tmp_path / name
    fake.write_text("#!/bin/sh\ntouch " + str(marker) + "\n")
    fake.chmod(0o700)
    code = (
        "from scripts.phase_4j76_source_resolution import scope; scope(); import subprocess; subprocess.run(["
        + repr(str(fake))
        + "],check=True)"
    )
    child = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert child.returncode and not marker.exists()


@pytest.fixture(scope="module")
def replay_artifacts(tmp_path_factory):
    root = tmp_path_factory.mktemp("j84_replay") / "evidence"
    j84.emit(REPO, root)
    assert (
        j84.replay(root)["decision"] == "audited_bundle_ingestion_ready_runtime_blocked"
    )
    return root


def test_benchmark_deterministic_missing_bundles(replay_artifacts):
    result = j84.replay(replay_artifacts)
    assert result["motifs"] == 8 and result["candidate_decisions"] == 40
    assert result["class_i_reports"] == 24 and result["class_ii_reports"] == 16
    assert (
        result["authorized_bundles_supplied"] == 0
        and result["accepted_coefficients"] == 0
    )
    reports = j84.read(replay_artifacts / "organic-coverage-report.json")
    assert all(
        r["missing_terms"]
        for m in reports["motifs"]
        for r in m["report"]["profile_reports"]
    )


def test_cli_ingest_and_fresh_replay_with_pinned_capture(tmp_path):
    _, args = bundle(tmp_path)
    out = tmp_path / "ingested"
    command = [
        sys.executable,
        "-m",
        "scripts.phase_4j84_bundles",
        "ingest",
        "--directory",
        str(out),
        "--manifest",
        str(args[0]),
        "--manifest-sha256",
        args[1],
        "--audit",
        str(args[2]),
        "--audit-sha256",
        args[3],
        "--capture-timestamp",
        WHEN,
    ]
    run = subprocess.run(command, capture_output=True, text=True, check=False)
    assert run.returncode == 0, run.stderr + run.stdout
    pin = json.loads(run.stdout)["snapshot_sha256"]
    child = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j84_bundles",
            "replay",
            "--directory",
            str(out),
            "--snapshot-sha256",
            pin,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr + child.stdout
    assert json.loads(child.stdout)["audited_bundle_candidates"] == 1
    with pytest.raises(ValueError, match="pin_mismatch"):
        j84.replay(out)


@pytest.mark.parametrize(
    "mutation",
    ("scope", "decision", "counts", "missing", "source_hash", "evidence", "class"),
)
def test_strict_receipt_rejects_resigned_semantic_mutations(
    tmp_path, monkeypatch, replay_artifacts, mutation
):
    from scripts import validate_phase_4j84_receipt as validator

    target = tmp_path / "evidence"
    shutil.copytree(replay_artifacts, target)
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"synthetic_control_only": True}
    )
    validator.seal(REPO, target)
    receipt = j84.read(target / "receipt.json")
    if mutation == "scope":
        receipt["scope"]["automatic_assignment_enabled"] = True
    elif mutation == "decision":
        receipt["decision"] = "assigned"
    elif mutation == "counts":
        receipt["result"]["authorized_bundles_supplied"] = 1
    else:
        name = "organic-coverage-report.json"
        obj = j84.read(target / name)
        report = obj["motifs"][0]["report"]
        if mutation == "missing":
            report["unresolved_terms"] = []
        elif mutation == "source_hash":
            report["source_hashes"] = ["0" * 64]
        elif mutation == "evidence":
            report["evidence_records"] = []
        else:
            report["profile_reports"][0]["functional_class"] = "I"
        j84.write(
            target / name, j84.signed({k: v for k, v in obj.items() if k != "identity"})
        )
        manifest = j84.read(target / "manifest.json")
        manifest["artifacts"][name] = j84.pin(target / name)
        j84.write(
            target / "manifest.json",
            j84.signed({k: v for k, v in manifest.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j84.pin(target / "manifest.json")
    j84.write(
        target / "receipt.json",
        j84.signed({k: v for k, v in receipt.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="mutation"):
        validator.validate(REPO, target)
