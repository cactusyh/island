"""Manufactured software controls only; no scientific force-field validation."""

import json
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import pytest

from island.forcefields.assignment_core import normalize, rebuild, type_graph
from island.forcefields.assignment_core.inventory import requests
from island.forcefields.assignment_core.schema import (
    CROSS_TERMS,
    EVIDENCE_STATUSES,
    READINESS,
    claim_identity,
    data,
    identity,
)
from island.forcefields.coverage_registry import (
    ChemicalMotifCatalog,
    ForceFieldProfileRegistry,
)
from scripts import phase_4j83_registry as j83
from tests.test_phase_4j82 import control, graph

REPO = Path(__file__).resolve().parents[1]
FIXTURES = j83.read(REPO / j83.FIXTURE)["motifs"]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j83.scope()


@pytest.mark.parametrize("fixture", FIXTURES, ids=lambda f: f["id"])
def test_motif_replay_and_no_hidden_atoms(fixture):
    g = normalize(j83.spec_from_data(fixture))
    c = ChemicalMotifCatalog()
    a, b = c.detect(g), c.detect(g)
    assert a == b and a
    assert {s for m in a for s in m.sites} == set(g.atom_order)
    assert all(len(m.sites) == len(set(m.sites)) for m in a)
    assert identity(a) == identity(b)


@pytest.mark.parametrize(
    "fixture,motif",
    (
        ("hydrocarbon", "aliphatic_ch"),
        ("aromatic", "aromatic"),
        ("ether", "ether"),
        ("ester", "ester"),
        ("carbonyl", "carbonyl"),
        ("amide", "amide"),
        ("nitrile", "nitrile"),
        ("halogenated", "halogenated_carbon"),
        ("sulfur", "sulfur"),
        ("phosphorus", "phosphorus"),
        ("silicon", "silicon"),
        ("end_groups", "end_groups"),
        ("crosslink", "crosslink_junctions"),
        ("branch", "crosslink_junctions"),
        ("unclassified_boron", "unclassified"),
    ),
)
def test_representative_motifs(fixture, motif):
    spec = next(f for f in FIXTURES if f["id"] == fixture)
    matches = ChemicalMotifCatalog().detect(normalize(j83.spec_from_data(spec)))
    assert motif in {m.motif for m in matches}


def test_siloxane_and_overlapping_carbonyl_are_descriptive():
    c = ChemicalMotifCatalog()
    assert "siloxane" in {m.motif for m in c.detect(graph("*O[Si](C)(C)*", 2))}
    m = c.detect(graph("*CC(=O)N*"))
    carbonyl = next(v for v in m if v.motif == "carbonyl")
    amide = next(v for v in m if v.motif == "amide")
    assert set(carbonyl.sites) <= set(amide.sites)


def test_end_group_builder_provenance_and_junction_boundary():
    from island.forcefields.assignment_core.schema import Crosslink

    g = graph("*CCC*", 2, end_groups=("*F", "*O"), crosslinks=(Crosslink(2, 5),))
    c = ChemicalMotifCatalog().detect(g)
    from island.forcefields.assignment_core.topology import molecule

    mol = molecule(g)
    caps = {s for m in c if m.motif == "end_groups" for s in m.sites}
    assert {
        a.GetIdx() + 1 for a in mol.GetAtoms() if a.GetSymbol() in ("F", "O")
    } <= caps
    report = ForceFieldProfileRegistry.contracts().coverage(g)
    cell = report.capability_matrix.lookup(
        "j82_contract_pcff", "crosslink_junctions", "finite_crosslinked"
    )
    assert cell is not None and {2, 5} <= set(cell.sites)
    assert "bond:2,5" in cell.request_ids
    assert any("topology_architecture_unsupported" in r for r in cell.reasons)


def test_empty_profiles_never_imply_scientific_coverage():
    r = ForceFieldProfileRegistry.contracts()
    assert {p.family for p in r.profiles} == {
        "OPLS-AA",
        "GAFF",
        "GAFF2",
        "PCFF",
        "COMPASS",
    }
    assert all(not p.parameters and not p.rules for p in r.profiles)
    result = r.coverage(graph())
    assert all(c.coverage_decision == "blocked" for c in result.capability_matrix.cells)
    assert all(
        c.reasons and c.unresolved_term_ids for c in result.capability_matrix.cells
    )
    assert all(data(result)[k] == v for k, v in READINESS.items())
    assert not result.runtime_registration and not result.source_hashes


@pytest.mark.parametrize("family", ("OPLS-AA", "GAFF", "GAFF2", "PCFF", "COMPASS"))
def test_synthetic_completeness_provisional_not_scientific(family):
    g, p = control(family=family)
    r = ForceFieldProfileRegistry((p,)).coverage(g)
    assert all(c.coverage_decision == "provisional" for c in r.capability_matrix.cells)
    assert all(not c.unresolved_term_ids for c in r.capability_matrix.cells)
    for c in r.capability_matrix.cells:
        cross = dict(c.coverage)["cross_terms"]
        assert cross.status == (
            "complete"
            if p.functional_class == "II"
            else "not_required_by_class_i_contract"
        )
        assert cross.required > 0 if p.functional_class == "II" else cross.required == 0
        assert "derived_provisional" in c.evidence_statuses
    assert r.source_hashes == (p.evidence[0].source_sha256,)
    assert dict(r.evidence_records)[p.id] == p.evidence
    assert not r.automatic_assignment_enabled and r.accepted_coefficients == 0


@pytest.mark.parametrize("status", EVIDENCE_STATUSES)
def test_status_preserved_without_promotion(status):
    g, p = control()
    e = replace(p.evidence[0], status=status, synthetic=False)
    p = replace(p, purpose="evidence_candidate", evidence=(e,))
    r = ForceFieldProfileRegistry((p,)).coverage(g)
    assert all(c.coverage_decision == "blocked" for c in r.capability_matrix.cells)
    assert all(status in c.evidence_statuses for c in r.capability_matrix.cells)
    assert r.evidence_records[0][1] == (e,)
    assert r.unresolved_terms[0][1]


def audited_control():
    """Manufactured audit, still only a test of the coverage acceptance branch."""
    g, p = control()
    e = replace(p.evidence[0], synthetic=False, status="verified_stock")
    p = replace(p, purpose="evidence_candidate", evidence=(e,))
    claims = [
        claim_identity(k, v)
        for k, values in (
            ("typing_rule", p.rules),
            ("parameter", p.parameters),
            ("convention", p.conventions),
            ("form", p.forms),
        )
        for v in values
    ]
    claims.append(claim_identity("graph_scope", g.identity))
    e = replace(e, audited_claims=tuple(sorted(set(claims))))
    return g, replace(p, evidence=(e,)), e


def test_audited_exact_graph_coverage_acceptance_keeps_runtime_blocked():
    g, p, e = audited_control()
    r = ForceFieldProfileRegistry((p,), (e,)).coverage(g)
    assert all(c.coverage_decision == "accepted" for c in r.capability_matrix.cells)
    assert not any(v for _, v in r.unresolved_terms)
    assert not r.automatic_assignment_enabled and not r.production_validated
    assert r.simulation_readiness == "not_established"
    assert all(
        c.reasons == ("audited_exact_graph_coverage_only_runtime_disabled",)
        for c in r.capability_matrix.cells
    )


@pytest.mark.parametrize("mutation", ("hash", "parameter", "graph_scope", "convention"))
def test_separate_audit_rejects_fabrication(mutation):
    g, p, e = audited_control()
    if mutation == "hash":
        p = replace(p, evidence=(replace(e, source_sha256="0" * 64),))
    elif mutation == "parameter":
        row = next(r for r in p.parameters if r.family == "bond")
        p = replace(
            p,
            parameters=tuple(
                replace(r, values=("3", "4")) if r.id == row.id else r
                for r in p.parameters
            ),
        )
    elif mutation == "graph_scope":
        g = rebuild(g, degree_of_polymerization=2)
    else:
        p = replace(
            p,
            conventions=tuple(
                replace(c, value="fabricated") if c.name == "mixing_rule" else c
                for c in p.conventions
            ),
        )
    result = ForceFieldProfileRegistry((p,), (e,)).coverage(g)
    assert all(c.coverage_decision == "blocked" for c in result.capability_matrix.cells)


@pytest.mark.parametrize(
    "family",
    (*CROSS_TERMS, "improper", "charge", "nonbonded", "torsion", "bond", "angle"),
)
def test_missing_terms_explicit_no_zero_fill(family):
    g, p = control(graph("*CCC*"), family="PCFF")
    p = replace(p, parameters=tuple(r for r in p.parameters if r.family != family))
    result = ForceFieldProfileRegistry((p,)).coverage(g)
    missing = dict(result.unresolved_terms)[p.id]
    assert any(m.family == family for m in missing)
    assert all(c.coverage_decision == "blocked" for c in result.capability_matrix.cells)
    assert (
        sum(
            cov.required - cov.covered
            for c in result.capability_matrix.cells
            for _, cov in c.coverage
        )
        > 0
    )


def test_cross_term_dependencies_preserved_in_motif_matrix():
    g, p = control(graph("*CCOCC*"), "PCFF")
    r = ForceFieldProfileRegistry((p,)).coverage(g)
    cell = r.capability_matrix.lookup(p.id, "ether", g.architecture)
    reqs = {v.id: v for v in requests(g, type_graph(g, p), p)}
    assert all(
        set(reqs[rid].dependencies) <= set(cell.request_ids) for rid in cell.request_ids
    )
    # Coverage cells overlap; they are never summed as unique global coverage.
    assert len(
        {rid for c in r.capability_matrix.cells for rid in c.request_ids}
    ) == len(reqs)
    assert sum(len(c.request_ids) for c in r.capability_matrix.cells) > len(reqs)


def test_unknown_chemistry_and_periodic_topology_block_even_complete_controls():
    for g in (graph("*CB(O)C*"), graph(periodic_topology=(("boundary", "periodic"),))):
        _, p = control(g)
        r = ForceFieldProfileRegistry((p,)).coverage(g)
        assert any(c.coverage_decision == "blocked" for c in r.capability_matrix.cells)
        assert dict(r.unresolved_motifs)[p.id]


def test_registry_and_catalog_identity_invariants():
    r = ForceFieldProfileRegistry.contracts()
    assert replace(r, profiles=tuple(reversed(r.profiles))).identity == r.identity
    with pytest.raises(ValueError, match="duplicate_profile"):
        replace(r, profiles=(r.profiles[0], r.profiles[0]))
    with pytest.raises(ValueError, match="duplicate_motif"):
        replace(r.catalog, motifs=(r.catalog.motifs[0], r.catalog.motifs[0]))
    with pytest.raises(FrozenInstanceError):
        r.profiles = ()
    assert ForceFieldProfileRegistry().coverage(graph()).ranked_profiles == ()
    assert (
        r.coverage(graph()).capability_matrix.lookup("absent", "ether", "finite_linear")
        is None
    )


def test_forged_graph_metadata_rejected():
    g = graph()
    with pytest.raises(ValueError, match="graph_spec_or_normalizer_mutation"):
        ChemicalMotifCatalog().detect(replace(g, architecture="finite_crosslinked"))


@pytest.mark.parametrize(
    "mutation", ("unit", "equation", "phase", "dependency", "conflict")
)
def test_convention_conflicts_and_dependencies_remain_explained(mutation):
    g, p = control(family="PCFF")
    row = next(r for r in p.parameters if r.family == "bond_angle")
    if mutation == "unit":
        bad = replace(row, form=replace(row.form, units=("unmatched", "units")))
    elif mutation == "equation":
        bad = replace(row, form=replace(row.form, equation="different_equation"))
    elif mutation == "phase":
        p = replace(
            p,
            conventions=tuple(
                replace(c, value="changed_phase") if c.name == "torsion_phase" else c
                for c in p.conventions
            ),
        )
        bad = row
    elif mutation == "dependency":
        bad = replace(row, equilibrium_references=())
    else:
        p = replace(
            p,
            parameters=(*p.parameters, replace(row, id="conflict", values=("3", "4"))),
        )
        bad = row
    p = replace(p, parameters=tuple(bad if r.id == row.id else r for r in p.parameters))
    result = ForceFieldProfileRegistry((p,)).coverage(g)
    assert all(c.coverage_decision == "blocked" for c in result.capability_matrix.cells)
    reasons = [m.reason for m in dict(result.unresolved_terms)[p.id]]
    expected = {
        "unit": "equation_units_or_ordering_mismatch",
        "equation": "equation_units_or_ordering_mismatch",
        "phase": "convention_identity_mismatch",
        "dependency": "equilibrium_reference_value_mismatch",
        "conflict": "conflicting_exact_parameters_no_selection",
    }[mutation]
    assert any(expected in r for r in reasons)


@pytest.mark.parametrize("other", ("GAFF", "GAFF2", "COMPASS", "PCFF"))
def test_no_cross_family_overlay(other):
    _, p = control()
    with pytest.raises(ValueError, match="cross_family_evidence"):
        replace(p, evidence=(replace(p.evidence[0], source_family=other),))


@pytest.mark.parametrize(
    "name",
    (
        "sqm",
        "lammps",
        "msi2lmp",
        "run_qm",
        "fit_hmds",
        "reserved_validation",
        "generate_parameters",
    ),
)
def test_forbidden_child_never_runs(tmp_path, name):
    marker = tmp_path / "executed"
    fake = tmp_path / name
    fake.write_text("#!/bin/sh\ntouch " + str(marker) + "\n")
    fake.chmod(0o700)
    command = (
        "from scripts.phase_4j76_source_resolution import scope; scope(); import subprocess; subprocess.run(["
        + repr(str(fake))
        + "],check=True)"
    )
    child = subprocess.run(
        [sys.executable, "-c", command], capture_output=True, text=True, check=False
    )
    assert child.returncode != 0 and not marker.exists()
    assert (
        "reject" in child.stderr.lower()
        or "forbidden" in child.stderr.lower()
        or "denied" in child.stderr.lower()
    )


@pytest.fixture(scope="module")
def replay_artifacts(tmp_path_factory):
    root = tmp_path_factory.mktemp("j83_replay") / "evidence"
    j83.emit(REPO, root)
    assert (
        j83.replay(root)["decision"] == "coverage_registry_implemented_runtime_blocked"
    )
    return root


def test_fresh_process_and_relocated_evidence_replay(tmp_path, replay_artifacts):
    dest = tmp_path / "evidence"
    shutil.copytree(replay_artifacts, dest)
    cmd = [
        sys.executable,
        "-m",
        "scripts.phase_4j83_registry",
        "replay",
        "--directory",
        str(dest),
    ]
    child = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert child.returncode == 0, child.stderr
    assert json.loads(child.stdout) == j83.replay(replay_artifacts)


@pytest.mark.parametrize(
    "mutation",
    (
        "scope",
        "readiness",
        "decision",
        "count",
        "missing",
        "hash",
        "evidence",
        "class",
        "matrix",
    ),
)
def test_resigned_receipt_semantic_mutations(
    tmp_path, monkeypatch, replay_artifacts, mutation
):
    from scripts import validate_phase_4j83_receipt as validator

    target = tmp_path / "evidence"
    shutil.copytree(replay_artifacts, target)
    # Only this manufactured receipt control replaces terminal/Git audit inputs.
    # The actual delivery uses the unmodified strict validator.
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"synthetic_control_only": True}
    )
    validator.seal(REPO, target)
    receipt = j83.read(target / "receipt.json")
    if mutation == "scope":
        receipt["scope"]["J83_QM_attempts"] = 1
    elif mutation == "readiness":
        receipt["scope"]["automatic_assignment_enabled"] = True
    elif mutation == "decision":
        receipt["decision"] = "assigned"
    elif mutation == "count":
        receipt["result"]["capability_cells"] = 0
    else:
        name = (
            "capability-matrix.json" if mutation == "matrix" else "coverage-report.json"
        )
        obj = j83.read(target / name)
        if mutation == "matrix":
            obj["matrices"][0]["matrix"]["cells"][0]["coverage_decision"] = "accepted"
        else:
            r = obj["motifs"][0]["report"]
            if mutation == "missing":
                r["unresolved_terms"] = []
            elif mutation == "hash":
                r["source_hashes"] = ["0" * 64]
            elif mutation == "evidence":
                r["evidence_records"] = []
            else:
                r["profile_reports"][0]["functional_class"] = "I"
        j83.write(
            target / name, j83.signed({k: v for k, v in obj.items() if k != "identity"})
        )
        manifest = j83.read(target / "manifest.json")
        manifest["artifacts"][name] = j83.pin(target / name)
        j83.write(
            target / "manifest.json",
            j83.signed({k: v for k, v in manifest.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j83.pin(target / "manifest.json")
    j83.write(
        target / "receipt.json",
        j83.signed({k: v for k, v in receipt.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="mutation"):
        validator.validate(REPO, target)


def test_captured_source_pin_rejects_resigned_input(tmp_path, replay_artifacts):
    dest = tmp_path / "evidence"
    shutil.copytree(replay_artifacts, dest)
    obj = j83.read(dest / "captured-input.json")
    obj["motifs"][0]["degree_of_polymerization"] = 5
    j83.write(dest / "captured-input.json", obj)
    with pytest.raises(ValueError, match="pin_mismatch"):
        j83.replay(dest)
