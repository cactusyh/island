"""Manufactured software controls only; no scientific force-field qualification."""

import json
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from island.forcefields.assignment_core import normalize, rebuild, type_graph
from island.forcefields.assignment_core.schema import FAMILIES
from island.forcefields.bundle_adapters import adapt, capture
from island.forcefields.bundle_adapters.qualification import (
    qualify,
    rank_qualifications,
)
from scripts import phase_4j85_qualification as j85
from tests.test_phase_4j82 import graph
from tests.test_phase_4j84 import WHEN, bundle

REPO = Path(__file__).resolve().parents[1]
MOTIFS = j85.read(REPO / j85.FIXTURE)["motifs"]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j85.scope()


def supplied(tmp_path, family="OPLS-AA", g=None, change=None, authorized=False):
    def declaration(m, records):
        records.append(
            {
                "record_id": "convention:energy_convention",
                "kind": "convention",
                "location": "Manufactured test table; explicit energy convention row",
                "origin": "directly_reported",
                "value": {
                    "name": "energy_convention",
                    "value": "Synthetic energy expression exactly as declared; no implicit prefactors",
                },
            }
        )
        if change:
            change(m, records)

    g, args = bundle(tmp_path, family, g=g, change=declaration, real_branch=authorized)
    return g, capture(*args, capture_timestamp=WHEN), args


def reports(result):
    return json.loads(result.graph_reports_json)


@pytest.mark.parametrize("family", FAMILIES)
def test_missing_bundle_explicitly_requires_source(family):
    result = qualify(capture(capture_timestamp=WHEN), (graph(),), family=family)
    assert result.decision == "source_bundle_required"
    assert not result.source_hashes and result.accepted_coefficients == 0
    assert result.source_status == "source_bundle_missing"
    assert not result.source_authorized
    assert reports(result)[0]["missing_terms"]
    assert reports(result)[0]["coverage"]["interaction_14"]["covered"] == 0


@pytest.mark.parametrize("family", FAMILIES)
def test_authorized_audit_branch_qualifies_only_declared_graph_software_control(
    tmp_path, family
):
    # Fabricated audit metadata intentionally exercises the trusted-input branch.
    # This fixture is NOT an authorized real force-field bundle or scientific result.
    g, cap, _ = supplied(tmp_path, family, authorized=True)
    result = qualify(cap, (g,), family=family)
    assert result.decision == "authorized_bundle_qualified_runtime_blocked"
    assert result.source_authorized
    assert (
        not result.blocking_reasons
        and result.source_status == "audited_bundle_candidate"
    )
    assert all(
        c["covered"] == c["required"] for c in reports(result)[0]["coverage"].values()
    )
    assert reports(result)[0]["coverage"]["interaction_14"]["required"] == 2
    cross = reports(result)[0]["coverage"]["class_ii_cross_terms"]
    assert cross["status"] == (
        "complete"
        if result.functional_class == "II"
        else "not_required_by_class_i_contract"
    )
    assert (
        cross["required"] > 0
        if result.functional_class == "II"
        else cross["required"] == 0
    )
    assert result.accepted_coefficients == 0 and not result.runtime_registration
    assert not result.automatic_assignment_enabled and not result.production_validated
    assert result.simulation_readiness == "not_established"
    evidence = json.loads(result.evidence_graph_json)
    assert set(evidence["source_hashes"]) == set(result.source_hashes)
    assert {f["sha256"] for f in evidence["adapter"]["source_files"]} == set(
        result.source_hashes
    )


@pytest.mark.parametrize("motif", MOTIFS, ids=lambda m: m["id"])
def test_each_organic_graph_software_coverage_with_audited_mapping(tmp_path, motif):
    g = normalize(j85.spec_from_data(motif))
    _, cap, _ = supplied(tmp_path, g=g, authorized=True)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert result.decision == "authorized_bundle_qualified_runtime_blocked"
    assert not reports(result)[0]["missing_terms"]


def test_synthetic_control_cannot_qualify_scientifically(tmp_path):
    g, cap, _ = supplied(tmp_path)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert result.decision == "source_bundle_required" and not result.source_authorized
    assert (
        "synthetic_software_control_not_scientific_evidence" in result.blocking_reasons
    )
    assert not reports(result)[0]["missing_terms"]


@pytest.mark.parametrize("empty", ((), iter(())))
def test_empty_graph_input_never_qualifies(tmp_path, empty):
    _, cap, _ = supplied(tmp_path, authorized=True)
    with pytest.raises(ValueError, match="qualification_requires_declared_graphs"):
        qualify(cap, empty, family="OPLS-AA")


def test_one_shot_graph_iterator_preserves_required_checks(tmp_path):
    g, cap, _ = supplied(tmp_path, authorized=True)
    first = qualify(cap, (g,), family="OPLS-AA")
    assert first == qualify(cap, iter((g,)), family="OPLS-AA")
    assert reports(first) and reports(first)[0]["graph_identity"] == g.identity


@pytest.mark.parametrize(
    "kind",
    ("no_audit", "unauthorized", "local_patch", "fallback", "polyjarvis", "license"),
)
def test_unauthorized_sources_never_qualify(tmp_path, kind):
    def change(m, _):
        if kind == "unauthorized":
            m["provenance"]["authorization"] = "not_authorized"
        elif kind in ("local_patch", "fallback"):
            m["evidence_status"] = kind
        elif kind == "polyjarvis":
            m["source_uri"] = "urn:polyjarvis:unverified-patch"
        elif kind == "license":
            m["provenance"]["license"] = "unknown"

    g, cap, args = supplied(tmp_path, change=change, authorized=True)
    if kind == "no_audit":
        cap = capture(*args[:2], capture_timestamp=WHEN)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert result.decision == "source_bundle_required" and result.blocking_reasons
    assert not result.source_authorized
    assert result.source_status == "source_parsed_audit_or_convention_blocked"
    assert result.source_hashes and result.accepted_coefficients == 0


@pytest.mark.parametrize("which", ("source", "manifest", "audit"))
def test_source_hash_mismatch_rejected_before_qualification(tmp_path, which):
    _, _, args = supplied(tmp_path)
    target = (
        args[0]
        if which == "manifest"
        else args[2]
        if which == "audit"
        else args[0].parent / "sources/records_0.json"
    )
    target.write_bytes(target.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA256_mutation"):
        capture(*args, capture_timestamp=WHEN)


@pytest.mark.parametrize(
    "mutation",
    ("unit", "parameter_unit", "energy_convention", "undefined_type", "dependency"),
)
def test_source_semantics_block_even_with_independent_audit(tmp_path, mutation):
    def change(m, rows):
        if mutation == "unit":
            m["unit_convention"]["energy_unit"] = "different_unit"
        elif mutation == "energy_convention":
            rows[:] = [
                r for r in rows if r["record_id"] != "convention:energy_convention"
            ]
        else:
            r = next(
                r["value"]
                for r in rows
                if r["kind"] == "parameter" and r["value"]["family"] == "bond"
            )
            if mutation == "parameter_unit":
                r["form"]["units"][0] = "incorrect_unit"
            elif mutation == "undefined_type":
                r["atom_types"][0] = "undeclared_type"
            else:
                r["dependency_keys"] = ["0" * 64]

    g, cap, _ = supplied(tmp_path, change=change, authorized=True)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert result.decision == "parameter_bundle_incomplete" and result.blocking_reasons
    assert result.source_authorized
    assert result.accepted_coefficients == 0
    if mutation == "undefined_type":
        assert "undefined_atom_type:undeclared_type" in result.blocking_reasons


@pytest.mark.parametrize(
    "family,section",
    (("OPLS-AA", "torsion"), ("PCFF", "bond_angle"), ("COMPASS", "angle_angle")),
)
def test_missing_parameter_sections_no_zero_substitution(tmp_path, family, section):
    def change(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (r["kind"] == "parameter" and r["value"]["family"] == section)
        ]

    g, cap, _ = supplied(tmp_path, family, change=change, authorized=True)
    result = qualify(cap, (g,), family=family)
    report = reports(result)[0]
    assert result.decision == "parameter_bundle_incomplete"
    assert report["sections"][section]["source_records"] == 0
    assert report["sections"][section]["covered_instances"] == 0
    assert report["sections"][section]["required_instances"] > 0
    assert any(
        m["family"] == section and "no_fallback" in m["reason"]
        for m in report["missing_terms"]
    )


def test_missing_terminal_terms_stay_missing(tmp_path):
    g = graph("*CC*", 2, end_groups=("*F", "*O"))

    def change(_, rows):
        rows[:] = [
            r
            for r in rows
            if not (
                r["kind"] == "parameter" and "control_F" in r["value"]["atom_types"]
            )
        ]

    _, cap, _ = supplied(tmp_path, g=g, change=change, authorized=True)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert result.decision == "parameter_bundle_incomplete"
    assert reports(result)[0]["coverage"]["atom_typing"]["status"] == "complete"
    assert reports(result)[0]["coverage"]["nonbonded"]["status"] == "incomplete"
    assert reports(result)[0]["missing_terms"]


def test_unused_source_row_with_undefined_type_blocks_qualification(tmp_path):
    def change(_, rows):
        old = next(
            r
            for r in rows
            if r["kind"] == "parameter" and r["value"]["family"] == "bond"
        )
        new = json.loads(json.dumps(old))
        new["record_id"] = "unused_bond_record"
        new["value"].update(
            id="unused_bond",
            environment="f" * 64,
            atom_types=["unused_undeclared", "control_C"],
        )
        rows.append(new)

    g, cap, _ = supplied(tmp_path, change=change, authorized=True)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert not reports(result)[0]["missing_terms"]
    assert "undefined_atom_type:unused_undeclared" in result.blocking_reasons
    assert result.decision == "parameter_bundle_incomplete"


def test_ambiguous_atom_types_are_not_tie_broken_by_guessing(tmp_path):
    def change(_, rows):
        old = next(
            r
            for r in rows
            if r["kind"] == "atom_type" and r["value"]["atom_type"] == "control_C"
        )
        new = json.loads(json.dumps(old))
        new["record_id"] = "atom_type:competing"
        new["value"].update(id="competing", atom_type="competing_C")
        rows.append(new)

    g, cap, _ = supplied(tmp_path, change=change, authorized=True)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert result == qualify(cap, (g,), family="OPLS-AA")
    assert result.decision == "parameter_bundle_incomplete"
    assert any("ambiguous_site" in r for r in result.blocking_reasons)


def test_unsupported_motif_blocked_despite_complete_software_rows(tmp_path):
    g = graph("*CB(O)C*")
    _, cap, _ = supplied(tmp_path, g=g, authorized=True)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert "unsupported_motif_not_in_catalog" in result.blocking_reasons
    assert result.decision == "parameter_bundle_incomplete"


@pytest.mark.parametrize("psmiles", ("*CCO*", "*CC(=O)O*"))
def test_unclassified_terminal_chemistry_is_not_hidden_by_internal_motif(
    tmp_path, psmiles
):
    g = graph(psmiles, 2)
    _, cap, _ = supplied(tmp_path, g=g, authorized=True)
    result = qualify(cap, (g,), family="OPLS-AA")
    assert result.decision == "parameter_bundle_incomplete"
    assert "unsupported_motif_not_in_catalog" in result.blocking_reasons
    assert not reports(result)[0]["missing_terms"]


@pytest.mark.parametrize("change", ("cut", "ends", "junction", "branch", "crosslink"))
def test_topology_changes_require_retyping_and_new_evidence(tmp_path, change):
    from island.forcefields.assignment_core.schema import Branch, Crosslink

    g = graph("*CO*", 2) if change in ("cut", "ends", "junction") else graph("*CCC*", 2)
    _, cap, _ = supplied(tmp_path, g=g, authorized=True)
    p = adapt(cap)
    typed = type_graph(g, p.profile, p.trusted_evidence)
    if change == "cut":
        changed = graph("*OC*", 2)
    elif change == "ends":
        changed = rebuild(g, end_groups=("*C", "*O"))
    elif change == "junction":
        changed = rebuild(
            g, repeat_units=(("A", "*CO*"), ("B", "*CC*")), sequence=("A", "B")
        )
    elif change == "branch":
        changed = rebuild(g, branches=(Branch(2, "*C"),))
    else:
        changed = rebuild(g, crosslinks=(Crosslink(2, 5),))
    with pytest.raises(ValueError, match="topology_changed_retyping_required"):
        typed.validate(changed, p.profile)
    assert (
        type_graph(changed, p.profile, p.trusted_evidence).graph_identity
        == changed.identity
    )
    result = qualify(cap, (changed,), family="OPLS-AA")
    assert "exact_graph_scope_not_evidenced" in result.blocking_reasons
    assert reports(result)[0]["coverage"]["topology"]["status"] == "incomplete"


def test_ranking_and_ties_are_deterministic_and_never_mix(tmp_path):
    outcomes = []
    for name in ("z_profile", "a_profile"):
        g, cap, _ = supplied(
            tmp_path / name,
            change=lambda m, _, profile_id=name: m.update(profile_id=profile_id),
            authorized=True,
        )
        outcomes.append(qualify(cap, (g,), family="OPLS-AA"))
    a = rank_qualifications(outcomes)
    assert a == rank_qualifications(outcomes[::-1])
    assert a == rank_qualifications(iter(outcomes))
    assert [r.profile_id for r in a] == ["a_profile", "z_profile"]
    with pytest.raises(ValueError, match="duplicate_candidate_bundle"):
        rank_qualifications((outcomes[0], outcomes[0]))
    with pytest.raises(ValueError, match="cross_family_qualification_forbidden"):
        qualify(cap, (g,), family="GAFF2")


def test_relocated_reingestion_preserves_qualification(tmp_path):
    g, cap, args = supplied(tmp_path / "original", authorized=True)
    root = tmp_path / "relocated"
    shutil.copytree(args[0].parent, root)
    audit = tmp_path / "audit-copy.json"
    shutil.copyfile(args[2], audit)
    fresh = capture(
        root / "manifest.json", args[1], audit, args[3], capture_timestamp=WHEN
    )
    first = qualify(cap, (g,), family="OPLS-AA")
    assert first == qualify(fresh, (g,), family="OPLS-AA")


@pytest.mark.parametrize(
    "field,value",
    (
        ("automatic_assignment_enabled", True),
        ("runtime_registration", True),
        ("production_validated", True),
        ("accepted_coefficients", 1),
    ),
)
def test_readiness_cannot_be_enabled(field, value):
    r = qualify(capture(capture_timestamp=WHEN), (graph(),), family="GAFF")
    with pytest.raises(ValueError, match="init=False"):
        replace(r, **{field: value})


@pytest.mark.parametrize(
    "name",
    (
        "sqm",
        "lammps",
        "msi2lmp",
        "fit_parameters",
        "generate_parameters",
        "reserved_validation",
    ),
)
def test_prohibited_child_cannot_execute(tmp_path, name):
    marker = tmp_path / "ran"
    fake = tmp_path / name
    fake.write_text("#!/bin/sh\ntouch " + str(marker) + "\n")
    fake.chmod(0o700)
    child = subprocess.run(
        [
            sys.executable,
            "-c",
            "from scripts.phase_4j76_source_resolution import scope; scope(); import subprocess; subprocess.run(["
            + repr(str(fake))
            + "])",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode and not marker.exists()


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    target = tmp_path_factory.mktemp("j85") / "evidence"
    j85.emit(REPO, target)
    return target


def test_default_benchmark_missing_bundles_has_all_dimensions(artifacts):
    result = j85.replay(artifacts)
    assert result["decision"] == "source_bundle_required"
    assert result["motifs"] == 8 and result["candidate_decisions"] == 40
    assert result["class_i_reports"] == 24 and result["class_ii_reports"] == 16
    assert result["qualified_graph_reports"] == result["accepted_coefficients"] == 0
    rows = j85.read(artifacts / "organic-coverage-report.json")["reports"]
    assert all(
        len(r["coverage"]) == 10 and r["missing_terms"] and r["blocking_reasons"]
        for r in rows
    )
    assert j85.read(artifacts / "source-evidence-graph.json")["source_hashes"] == []


def test_cli_supplied_bundle_and_fresh_process_replay(tmp_path):
    _, _, args = supplied(tmp_path)
    target = tmp_path / "intake"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j85_qualification",
            "ingest",
            "--directory",
            str(target),
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
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr + run.stdout
    output = json.loads(run.stdout)
    assert output["result"]["decision"] == "source_bundle_required"
    assert output["result"]["authorized_bundles_supplied"] == 0
    assert output["result"]["synthetic_bundles_supplied"] == 1
    relocated = tmp_path / "relocated-evidence"
    shutil.copytree(target, relocated)
    child = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j85_qualification",
            "replay",
            "--directory",
            str(relocated),
            "--snapshot-sha256",
            output["snapshot_sha256"],
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr + child.stdout
    assert json.loads(child.stdout) == output["result"]


@pytest.mark.parametrize(
    "mutation", ("scope", "decision", "source", "missing", "14", "class", "receipt")
)
def test_strict_receipt_rejects_resigned_mutations(
    tmp_path, monkeypatch, artifacts, mutation
):
    from scripts import validate_phase_4j85_receipt as validator

    target = tmp_path / "evidence"
    shutil.copytree(artifacts, target)
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"software_control_only": True}
    )
    validator.seal(REPO, target)
    receipt = j85.read(target / "receipt.json")
    if mutation in ("scope", "decision"):
        if mutation == "scope":
            receipt["scope"]["automatic_assignment_enabled"] = True
        else:
            receipt["decision"] = "authorized_bundle_qualified_runtime_blocked"
    else:
        name = (
            "qualification-receipt.json"
            if mutation == "receipt"
            else "organic-coverage-report.json"
        )
        obj = j85.read(target / name)
        if mutation == "receipt":
            obj["decision"] = "authorized_bundle_qualified_runtime_blocked"
        elif mutation == "source":
            obj["reports"][0]["source_hashes"] = ["0" * 64]
        elif mutation == "missing":
            obj["reports"][0]["missing_terms"] = []
        elif mutation == "14":
            obj["reports"][0]["coverage"]["interaction_14"]["covered"] = 2
        else:
            obj["reports"][0]["functional_class"] = "III"
        j85.write(
            target / name, j85.signed({k: v for k, v in obj.items() if k != "identity"})
        )
        manifest = j85.read(target / "manifest.json")
        manifest["artifacts"][name] = j85.pin(target / name)
        j85.write(
            target / "manifest.json",
            j85.signed({k: v for k, v in manifest.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j85.pin(target / "manifest.json")
    j85.write(
        target / "receipt.json",
        j85.signed({k: v for k, v in receipt.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="mutation"):
        validator.validate(REPO, target)
