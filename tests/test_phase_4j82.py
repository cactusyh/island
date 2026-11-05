"""Synthetic software controls only. None qualifies a scientific force field."""

import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from island.forcefields.assignment_core import (
    PolymerSpec,
    builtin_contracts,
    normalize,
    rebuild,
    resolve,
    type_graph,
    validate_class_i,
    validate_class_ii,
)
from island.forcefields.assignment_core.inventory import requests
from island.forcefields.assignment_core.profiles import (
    required_conventions,
    required_families,
)
from island.forcefields.assignment_core.schema import (
    CROSS_TERMS,
    EVIDENCE_STATUSES,
    Branch,
    Convention,
    Crosslink,
    EvidenceRecord,
    ForceFieldProfile,
    ParameterRecord,
    TermForm,
    TypingRule,
    claim_identity,
    data,
    identity,
)
from island.forcefields.assignment_core.topology import molecule
from scripts import phase_4j82_coverage as j82

REPO = Path(__file__).resolve().parents[1]
MOTIFS = j82.read(REPO / j82.FIXTURE)["motifs"]


@pytest.fixture(scope="module", autouse=True)
def fenced():
    j82.scope()


def graph(psmiles="*CC*", dp=1, **kw):
    return normalize(PolymerSpec((("A", psmiles),), dp, **kw))


def control(g=None, family="OPLS-AA"):
    """Explicit synthetic numbers exercise completeness; never a scientific profile."""
    g = g or graph()
    base = next(p for p in builtin_contracts() if p.family == family)
    e = EvidenceRecord(
        "synthetic-source",
        "derived_provisional",
        identity(("synthetic", family)),
        "Synthetic software control",
        "generated test rows",
        family,
        True,
    )
    elements = sorted(
        {a.GetSymbol(): a.GetAtomicNum() for a in molecule(g).GetAtoms()}.items()
    )
    rules = tuple(
        TypingRule("rule_" + symbol, "control_" + symbol, f"[#{number}:1]", e.id)
        for symbol, number in elements
    )
    forms = tuple(
        TermForm(
            f,
            "synthetic_" + f,
            ("q",) if f == "charge" else ("reference", "k"),
            ("e",) if f == "charge" else ("synthetic_unit", "synthetic_unit"),
            "explicit_order_no_symmetry",
            ("reference",) if f in ("bond", "angle") else (),
        )
        for f in required_families(base)
    )
    conventions = tuple(
        Convention(
            n,
            "0"
            if n == "charge_tolerance"
            else "synthetic_nonbonded"
            if n == "nonbonded_form"
            else "explicit_synthetic_" + n,
            e.id,
        )
        for n in required_conventions(base)
    )
    p = ForceFieldProfile(
        "control_" + family,
        family,
        base.functional_class,
        rules,
        forms,
        conventions,
        (e,),
        (),
        (g.identity,),
        ("finite_linear", "finite_branched", "finite_crosslinked"),
        "synthetic_control",
    )
    typed = type_graph(g, p)
    rows = requests(g, typed, p)
    by_id = {r.id: r for r in rows}
    parameters = {}
    for r in rows:
        f = next(f for f in forms if f.family == r.family)
        parameters[r.parameter_key] = ParameterRecord(
            "parameter_" + r.parameter_key,
            family,
            r.family,
            r.atom_types,
            r.environment,
            f,
            ("0",) if r.family == "charge" else ("1", "2"),
            e.id,
            p.convention_identity,
            tuple(sorted({by_id[d].parameter_key for d in r.dependencies})),
            source_location="synthetic_generated_row:" + r.parameter_key,
        )
    refs = {}
    for k, r in parameters.items():
        refs[k] = tuple(
            (f, dict(zip(r.form.fields, r.values, strict=True))[f])
            for f in r.form.equilibrium_fields
        )
    parameters = {
        k: replace(
            r, equilibrium_references=tuple((d, refs[d]) for d in r.dependency_keys)
        )
        for k, r in parameters.items()
    }
    return g, replace(p, parameters=tuple(parameters[k] for k in sorted(parameters)))


def report(g, p):
    typed = type_graph(g, p)
    return (validate_class_ii if p.functional_class == "II" else validate_class_i)(
        g, typed, p
    )


@pytest.mark.parametrize("motif", MOTIFS, ids=lambda m: m["id"])
def test_deterministic_motifs_final_graph(motif):
    spec = j82.spec_from_data(motif)
    a, b = normalize(spec), normalize(spec)
    assert data(a) == data(b) and a.identity == b.identity
    assert all(
        a.GetAtomicNum() and not a.GetNumImplicitHs() and not a.GetNumExplicitHs()
        for a in molecule(a).GetAtoms()
    )
    assert a.blockers == (
        ("periodic_topology_requires_explicit_image_graph_adapter",)
        if spec.periodic_topology
        else ()
    )


@pytest.mark.parametrize("family", ("OPLS-AA", "GAFF", "GAFF2", "PCFF", "COMPASS"))
def test_complete_synthetic_profile_never_runtime_assigns(family):
    g, p = control(family=family)
    r = report(g, p)
    assert r.status == "synthetic_complete" and not r.missing_terms
    assert not r.blocking_reasons
    d = resolve(g, (p,))
    assert d.status == "blocked" and not d.automatic_assignment_enabled
    assert d.accepted_coefficients == 0 and not d.production_validated
    assert d.simulation_readiness == "not_established"


def test_empty_scientific_contracts_no_unverified_numbers():
    g = graph()
    d = resolve(g, builtin_contracts())
    assert len(d.reports) == 5
    assert all(p.parameters == () and p.rules == () for p in builtin_contracts())
    for r in d.reports:
        assert r.status == "incomplete" and r.missing_terms
        assert r.source_hashes == ()
        assert "no_evidenced_typing_rules" in r.blocking_reasons
        assert any(b.startswith("convention_missing:") for b in r.blocking_reasons)


def test_class_i_class_ii_reports_distinct():
    a, pa = control()
    b, pb = control(family="PCFF")
    ra, rb = report(a, pa), report(b, pb)
    assert ra.functional_class == "I" and rb.functional_class == "II"
    assert dict(ra.coverage)["cross_terms"].status == "not_required_by_class_i_contract"
    assert dict(rb.coverage)["cross_terms"].required > 0
    with pytest.raises(ValueError, match="class_ii_validator"):
        validate_class_ii(a, type_graph(a, pa), pa)
    with pytest.raises(ValueError, match="class_i_validator"):
        validate_class_i(b, type_graph(b, pb), pb)


@pytest.mark.parametrize("family", CROSS_TERMS)
def test_every_required_class_ii_coupling_checked(family):
    g, p = control(graph("*CCC*", dp=2), "PCFF")
    assert any(r.family == family for r in p.parameters)
    p = replace(p, parameters=tuple(r for r in p.parameters if r.family != family))
    result = report(g, p)
    assert any(r.family == family for r in result.missing_terms)
    assert result.status == "incomplete"


def test_missing_dependency_blocks_coupled_terms():
    g, p = control(family="PCFF")
    p = replace(p, parameters=tuple(r for r in p.parameters if r.family != "bond"))
    r = report(g, p)
    assert any(
        m.reason == "required_coupling_dependency_unresolved" for m in r.missing_terms
    )


def test_wrong_dependency_parameter_identity():
    g, p = control(family="PCFF")
    rows = tuple(
        replace(r, dependency_keys=()) if r.family == "bond_angle" else r
        for r in p.parameters
    )
    r = report(g, replace(p, parameters=rows))
    assert any(
        m.reason == "class_ii_equilibrium_dependency_mismatch" for m in r.missing_terms
    )


@pytest.mark.parametrize(
    "operation", ("end", "branch", "crosslink", "dp", "stereo", "periodic")
)
def test_topology_change_requires_retyping(operation):
    g, p = control(graph("*CCC*", dp=2))
    typed = type_graph(g, p)
    changes = {
        "end": {"end_groups": ("*C", None)},
        "branch": {"branches": (Branch(2, "*C"),)},
        "crosslink": {"crosslinks": (Crosslink(2, 5),)},
        "dp": {"degree_of_polymerization": 3},
        "stereo": {"stereo_seed": 42},
        "periodic": {"periodic_topology": (("boundary", "periodic"),)},
    }
    altered = rebuild(g, **changes[operation])
    with pytest.raises(ValueError, match="retyping_required"):
        validate_class_i(altered, typed, p)


def test_typing_priority_and_tie_conflict():
    g, p = control()
    carbon = next(r for r in p.rules if r.atom_type == "control_C")
    special = replace(carbon, id="competing", atom_type="special_C")
    conflict = replace(p, rules=(*p.rules, special))
    typed = type_graph(g, conflict)
    assert any("ambiguous_site" in d for d in typed.diagnostics)
    priority = replace(p, rules=(*p.rules, replace(special, priority=1)))
    assert {t for _, t in type_graph(g, priority).assignments} == {
        "special_C",
        "control_H",
    }
    explicit = replace(p, rules=(*p.rules, replace(special, overrides=(carbon.id,))))
    assert type_graph(g, explicit).assignments == type_graph(g, priority).assignments


def test_override_cycle_rejected():
    g, p = control()
    c = next(r for r in p.rules if r.atom_type == "control_C")
    special = replace(c, id="cycle", overrides=(c.id,))
    rules = tuple(
        replace(r, overrides=("cycle",)) if r.id == c.id else r for r in p.rules
    )
    with pytest.raises(ValueError, match="cycle"):
        type_graph(g, replace(p, rules=(*rules, special)))


def test_stereospecific_smarts_not_silently_ignored():
    g, p = control(graph("*C[C@H](C)*"))
    p = replace(
        p,
        rules=(
            *p.rules,
            TypingRule("chiral", "chiral_type", "[C@:1]", p.evidence[0].id),
        ),
    )
    assert any("chiral_engine" in b for b in report(g, p).blocking_reasons)


def test_no_element_sequence_substitution_for_missing_type():
    g, p = control()
    p = replace(p, rules=tuple(r for r in p.rules if r.atom_type != "control_H"))
    r = report(g, p)
    assert dict(r.coverage)["atom_typing"].covered == 2
    assert any(m.reason == "untyped_or_ambiguous_site" for m in r.missing_terms)


@pytest.mark.parametrize("field", ("equation", "units", "ordering"))
def test_equation_units_order_mismatch_blocks(field):
    g, p = control()
    row = next(r for r in p.parameters if r.family == "bond")
    value = (
        ("wrong_unit",) * len(row.form.units)
        if field == "units"
        else "different_convention"
    )
    bad = replace(row, form=replace(row.form, **{field: value}))
    r = report(
        g,
        replace(
            p, parameters=tuple(bad if t.id == row.id else t for t in p.parameters)
        ),
    )
    assert any("equation_units_or_ordering" in m.reason for m in r.missing_terms)


def test_nonbonded_9_6_vs_12_6_never_coerced():
    g, p = control()
    conventions = tuple(
        replace(c, value="9-6") if c.name == "nonbonded_form" else c
        for c in p.conventions
    )
    assert (
        "nonbonded_9_6_vs_12_6_or_equation_mismatch"
        in report(g, replace(p, conventions=conventions)).blocking_reasons
    )


@pytest.mark.parametrize("status", EVIDENCE_STATUSES)
def test_evidence_labels_alone_cannot_promote(status):
    g, p = control()
    e = replace(p.evidence[0], status=status, synthetic=False)
    p = replace(p, purpose="evidence_candidate", evidence=(e,))
    r = report(g, p)
    assert r.status == "incomplete"
    assert any(
        "unqualified_evidence" in b or "independent_evidence_audit_missing" in b
        for b in r.blocking_reasons
    )
    if status in ("verified_stock", "verified_extension"):
        claims = [claim_identity("typing_rule", r) for r in p.rules]
        claims += [claim_identity("parameter", r) for r in p.parameters]
        claims += [claim_identity("convention", c) for c in p.conventions]
        claims += [claim_identity("form", f) for f in p.forms]
        claims += [claim_identity("graph_scope", g.identity)]
        e = replace(e, audited_claims=tuple(sorted(set(claims))))
        p = replace(p, evidence=(e,))
        t = type_graph(g, p, (e,))
        r = validate_class_i(g, t, p, (e,))
        assert r.status == "coverage_complete_not_activated"
        assert resolve(g, (p,), (e,)).automatic_assignment_enabled is False


def test_source_hash_mutation_invalidates_independent_audit():
    g, p = control()
    e = replace(p.evidence[0], status="verified_stock", synthetic=False)
    bad = replace(e, source_sha256="0" * 64)
    p = replace(p, purpose="evidence_candidate", evidence=(bad,))
    r = validate_class_i(g, type_graph(g, p, (e,)), p, (e,))
    assert any("independent_evidence_audit_missing" in b for b in r.blocking_reasons)


@pytest.mark.parametrize("family", ("GAFF", "GAFF2", "PCFF", "COMPASS"))
def test_cross_family_overlays_rejected(family):
    _, p = control()
    bad = replace(p.parameters[0], profile_family=family)
    with pytest.raises(ValueError, match="cross_family_parameter"):
        replace(p, parameters=(bad,))
    with pytest.raises(ValueError, match="cross_family_evidence"):
        replace(p, evidence=(replace(p.evidence[0], source_family=family),))


def test_compass_pcff_and_gaff_versions_distinct():
    contracts = {p.family: p for p in builtin_contracts()}
    assert contracts["PCFF"].identity != contracts["COMPASS"].identity
    assert contracts["GAFF"].identity != contracts["GAFF2"].identity


def test_conflicting_rows_preserved_and_blocked():
    g, p = control()
    a = next(r for r in p.parameters if r.family == "bond")
    conflict = replace(a, id="unequal", values=("1", "3"))
    r = report(g, replace(p, parameters=(*p.parameters, conflict)))
    m = next(m for m in r.missing_terms if m.family == "bond")
    assert m.reason == "conflicting_exact_parameters_no_selection"
    assert set(m.candidate_ids) == {a.id, "unequal"}


def test_equal_duplicate_keeps_both_identities():
    g, p = control()
    a = next(r for r in p.parameters if r.family == "bond")
    r = report(g, replace(p, parameters=(*p.parameters, replace(a, id="equal"))))
    assert r.status == "synthetic_complete"
    assert any(set(ids) == {a.id, "equal"} for _, ids in r.matched_records)


def test_duplicate_identity_rejected():
    _, p = control()
    with pytest.raises(ValueError, match="duplicate_parameter_identity"):
        replace(p, parameters=(*p.parameters, p.parameters[0]))


def test_no_missing_value_zero_fill_or_wildcard():
    _, p = control()
    a = p.parameters[0]
    with pytest.raises(ValueError, match="no_zero_fill"):
        replace(a, values=())
    with pytest.raises(ValueError, match="wildcard"):
        replace(a, atom_types=("*",) * len(a.atom_types))


def test_explicit_evidenced_nonapplicability_is_not_zero_fill():
    g, p = control(family="PCFF")
    a = next(r for r in p.parameters if r.family == "improper")
    absent = replace(
        a,
        values=(),
        applicability="not_applicable",
        absence_reason="synthetic omission control",
    )
    assert (
        report(
            g,
            replace(
                p, parameters=tuple(absent if r.id == a.id else r for r in p.parameters)
            ),
        ).status
        == "synthetic_complete"
    )
    with pytest.raises(ValueError, match="unevidenced_term_omission"):
        replace(absent, absence_reason="")


def test_charges_not_neutralized_implicitly():
    g, p = control()
    rows = tuple(
        replace(r, values=("0.5",)) if r.family == "charge" else r for r in p.parameters
    )
    assert any(
        "component_charge_not_conserved" in b
        for b in report(g, replace(p, parameters=rows)).blocking_reasons
    )


def test_typing_snapshot_and_order_mutation():
    g, p = control()
    typed = type_graph(g, p)
    with pytest.raises(ValueError, match="typing_snapshot_mutation"):
        validate_class_i(g, replace(typed, assignments=typed.assignments[::-1]), p)
    bad = replace(g, atom_order=g.atom_order[::-1])
    with pytest.raises(ValueError, match="atom_order"):
        molecule(bad)


@pytest.mark.parametrize(
    "kwargs",
    (
        {"degree_of_polymerization": 0},
        {"degree_of_polymerization": True},
        {"sequence": ("undefined",)},
        {"connection_points": (("A", 9, 10),)},
        {"crosslinks": (Crosslink(1, 1),)},
        {"branches": (Branch(999, "*C"),)},
    ),
)
def test_invalid_topology_fails_explicitly(kwargs):
    base = {"repeat_units": (("A", "*CC*"),), "degree_of_polymerization": 1}
    base.update(kwargs)
    with pytest.raises(ValueError):
        normalize(PolymerSpec(**base))


def test_candidate_ranking_independent_of_input_order():
    g, complete = control()
    empty = builtin_contracts()
    a = resolve(g, (*empty, complete))
    b = resolve(g, (complete, *reversed(empty)))
    assert data(a) == data(b) and a.ranked_profiles[0] == complete.id
    assert resolve(g, ()).blocking_reasons == (
        "J82_runtime_assignment_disabled",
        "no_candidate_profiles",
    )


def test_unique_parameter_identity_not_instance_alias():
    g, p = control(graph("*CC*", 3))
    rows = requests(g, type_graph(g, p), p)
    assert len({r.parameter_key for r in rows}) < len(rows)
    assert len({r.id for r in rows}) == len(rows)


def test_exact_environment_and_type_mutation_no_fallback():
    g, p = control()
    rows = tuple(
        replace(r, environment="0" * 64) if r.family == "bond" else r
        for r in p.parameters
    )
    assert any(
        m.family == "bond" for m in report(g, replace(p, parameters=rows)).missing_terms
    )
    rows = tuple(
        replace(r, atom_types=("another_type",) * len(r.atom_types))
        for r in p.parameters
    )
    assert report(g, replace(p, parameters=rows)).matched_records == ()


def test_no_runtime_module_mutation():
    import island.polymer as public
    from island.forcefields.pcff import registry

    before = public.prepare_polymer
    registry_bytes = Path(registry.__file__).read_bytes()
    g, p = control()
    resolve(g, (p,))
    assert public.prepare_polymer is before
    assert Path(registry.__file__).read_bytes() == registry_bytes
    assert not hasattr(
        __import__("island.forcefields.assignment_core", fromlist=[""]), "register"
    )


@pytest.mark.parametrize(
    "name",
    (
        "sqm",
        "lammps",
        "msi2lmp",
        "run_qm",
        "run_md",
        "fit_hmds",
        "evaluate_energy",
        "reserved_validation",
    ),
)
def test_prohibited_synthetic_child_cannot_execute(tmp_path, name):
    fake = tmp_path / name
    sentinel = tmp_path / "executed"
    fake.write_text("#!/bin/sh\nexit 0\n")
    fake.chmod(0o700)
    code = (
        "import subprocess; from pathlib import Path; "
        f"subprocess.run([{str(fake)!r}],check=True); Path({str(sentinel)!r}).touch()"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert result.returncode != 0 and not sentinel.exists()
    assert "forbidden" in result.stderr.lower() or "scope" in result.stderr.lower()


def test_fence_required_for_cli(tmp_path, monkeypatch):
    def unavailable():
        raise RuntimeError("fence_unavailable")

    monkeypatch.setattr(j82, "scope", unavailable)
    with pytest.raises(RuntimeError, match="fence_unavailable"):
        j82.emit(REPO, tmp_path / "no_fence")
    assert not (tmp_path / "no_fence").exists()


def test_replay_and_resigned_semantic_mutation(tmp_path):
    target = tmp_path / "evidence"
    j82.emit(REPO, target)
    expected = j82.replay(target)
    assert expected["scope"]["accepted_coefficients"] == 0
    result = j82.read(target / "result.json")
    result["scope"]["automatic_assignment_enabled"] = True
    j82.write(
        target / "result.json",
        j82.signed({k: v for k, v in result.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="resigned_semantic_mutation"):
        j82.replay(target)


def test_relocated_replay(tmp_path):
    target = tmp_path / "relocated"
    j82.emit(REPO, target)
    # Replay needs no original coordinates or force-field assets.
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.phase_4j82_coverage",
            "replay",
            "--directory",
            str(target),
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert (
        json.loads(result.stdout)["decision"]
        == "coverage_core_implemented_runtime_blocked"
    )


def test_historical_receipt_strict_validation_and_byte_preservation():
    from scripts.validate_phase_4j81_receipt import validate

    result = validate(REPO, REPO / "docs/evidence/phase_4j81")
    assert result["status"] == "passed"
    assert len(j82.historical(REPO)) > 10619


def test_equilibrium_reference_mutation_blocks_cross_term():
    g, p = control(family="PCFF")
    row = next(r for r in p.parameters if r.family == "bond_angle")
    bad = replace(row, equilibrium_references=())
    r = report(
        g,
        replace(
            p, parameters=tuple(bad if a.id == row.id else a for a in p.parameters)
        ),
    )
    assert any(
        m.reason == "equilibrium_reference_value_mismatch" for m in r.missing_terms
    )


@pytest.mark.parametrize("resolution", ("united_atom", "coarse_grained"))
def test_resolution_overlay_rejected(resolution):
    _, p = control()
    with pytest.raises(ValueError, match="resolution_overlay"):
        replace(p, evidence=(replace(p.evidence[0], resolution=resolution),))


def test_independently_audited_source_does_not_authenticate_fabricated_rows():
    g, p = control()
    e = replace(p.evidence[0], synthetic=False, status="verified_stock")
    claims = [claim_identity("typing_rule", r) for r in p.rules]
    claims += [claim_identity("parameter", r) for r in p.parameters]
    claims += [claim_identity("convention", c) for c in p.conventions]
    claims += [claim_identity("form", f) for f in p.forms]
    claims += [claim_identity("graph_scope", g.identity)]
    e = replace(e, audited_claims=tuple(sorted(set(claims))))
    p = replace(p, purpose="evidence_candidate", evidence=(e,))
    row = next(r for r in p.parameters if r.family == "bond")
    p = replace(
        p,
        parameters=tuple(
            replace(r, values=("1", "999")) if r.id == row.id else r
            for r in p.parameters
        ),
    )
    r = validate_class_i(g, type_graph(g, p, (e,)), p, (e,))
    assert any("claim_not_independently_audited" in m.reason for m in r.missing_terms)


def test_schema_does_not_retain_mutable_input_aliases():
    values = [["A", "*CC*"]]
    spec = PolymerSpec(values, 1)
    values[0][1] = "*CO*"
    assert spec.repeat_units == (("A", "*CC*"),)
    rule = TypingRule("r", "t", "[C:1]", "e", overrides=["base"])
    assert rule.overrides == ("base",)


@pytest.fixture(scope="module")
def replay_artifacts(tmp_path_factory):
    root = tmp_path_factory.mktemp("j82_replay_artifacts") / "evidence"
    j82.emit(REPO, root)
    return root


@pytest.mark.parametrize(
    "mutation",
    ("scope", "readiness", "decision", "counts", "missing_terms", "source_hash"),
)
def test_strict_validator_rejects_resigned_semantic_mutations(
    tmp_path, monkeypatch, replay_artifacts, mutation
):
    import shutil

    from scripts import validate_phase_4j82_receipt as validator

    target = tmp_path / "evidence"
    shutil.copytree(replay_artifacts, target)
    # Synthetic terminal-audit stand-ins isolate receipt semantics in this control.
    # The actual sealed delivery uses the unpatched validator and real Git/fence data.
    monkeypatch.setattr(validator, "historical", lambda _: {})
    monkeypatch.setattr(
        validator, "verification", lambda *_: {"synthetic_control_only": True}
    )
    validator.seal(REPO, target)
    receipt = j82.read(target / "receipt.json")
    if mutation == "scope":
        receipt["scope"]["J82_QM_attempts"] = 1
    elif mutation == "readiness":
        receipt["scope"]["automatic_assignment_enabled"] = True
    elif mutation == "decision":
        receipt["decision"] = "assigned"
    elif mutation == "counts":
        receipt["result"]["candidate_decisions"] = 0
    else:
        name = "coverage-report.json"
        artifact = j82.read(target / name)
        report0 = artifact["motifs"][0]["decision"]["reports"][0]
        if mutation == "missing_terms":
            report0["missing_terms"] = []
        else:
            report0["source_hashes"] = ["0" * 64]
        j82.write(
            target / name,
            j82.signed({k: v for k, v in artifact.items() if k != "identity"}),
        )
        manifest = j82.read(target / "manifest.json")
        manifest["artifacts"][name] = j82.pin(target / name)
        j82.write(
            target / "manifest.json",
            j82.signed({k: v for k, v in manifest.items() if k != "identity"}),
        )
        receipt["manifest_sha256"] = j82.pin(target / "manifest.json")
    j82.write(
        target / "receipt.json",
        j82.signed({k: v for k, v in receipt.items() if k != "identity"}),
    )
    with pytest.raises(ValueError, match="mutation"):
        validator.validate(REPO, target)
