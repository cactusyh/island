"""Manufactured software controls only. No scientific coefficients or coverage."""

import copy
import json
import shutil
from pathlib import Path

import pytest

from island.forcefields.assignment_core.schema import data, identity
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.opls_providers import (
    FLAGS,
    OPLSProviderRegistry,
    benchmarks,
    finite_graph,
)
from island.forcefields.opls_providers.export import export_parameters, read_export
from island.forcefields.opls_providers.registry import obligations
from island.forcefields.opls_providers.schema import CONVENTIONS, FORMS, PROVIDERS
from scripts.phase_4j73_source_intake import pin, read, signed, write
from scripts.phase_4j90_providers import emit, replay


def settings():
    return copy.deepcopy(benchmarks()["PE"])


def conventions():
    return {
        "energy_unit": "kJ/mol",
        "length_unit": "nm",
        "charge_unit": "e",
        "angle_equilibrium_unit": "degree",
        "angle_force_unit": "radian",
        "nonbonded": "lj12_6_sigma_zero",
        "mixing_rule": "geometric_sigma_epsilon",
        "exclusions": "graph_distance_1_2",
        "scale14_lj": "0.5",
        "scale14_coulomb": "0.5",
        "long_range": "finite_NoCutoff",
        "cutoff": "none_finite",
        "tail": "none",
        "torsion_phase": "cos(n*phi-phase)",
        "torsion_normalization": "form_explicit",
        "improper_ordering": "source_ordered_center_second",
        "atom_typing": "explicit_provider_site_assignment",
        "charge_origin": "manufactured_software_control_not_scientific_data",
    }


def control(
    directory, *, provider="ligpargen", spec=None, change=None, native=False, db=False
):
    """Never persists into the checkout; generated literals test schema only."""
    directory.mkdir(parents=True, exist_ok=True)
    spec = spec or settings()
    graph = finite_graph(spec)
    mol = molecule(graph)
    needed, exclusions = obligations(graph)
    atoms = [
        {
            "site": a.GetIdx() + 1,
            "atom_type": "control_" + str(a.GetIdx() + 1),
            "element": a.GetSymbol(),
            "mass": str(a.GetMass()),
            "location": "table:software_atoms;row:" + str(a.GetIdx() + 1),
        }
        for a in mol.GetAtoms()
    ]
    records = []
    forms = {
        "bond": "bond_harmonic",
        "angle": "angle_harmonic",
        "torsion": "torsion_periodic",
        "charge": "charge",
        "nonbonded": "lj12_6",
    }
    literals = {
        "bond": [["0.1", "1"]],
        "angle": [["100", "1"]],
        "torsion": [["180", "1", "2"]],
        "charge": [["0"]],
        "nonbonded": [["0.1", "1"]],
    }
    for family, sites in needed.items():
        if family not in (*forms, "improper"):
            continue
        for i, s in enumerate(sites):
            records.append(
                {
                    "id": family + str(i),
                    "family": family,
                    "sites": list(s),
                    "form": forms.get(family),
                    "values": copy.deepcopy(literals.get(family, [])),
                    "applicability": "not_applicable"
                    if family == "improper"
                    else "required",
                    "reason": "explicit software-only omission"
                    if family == "improper"
                    else "",
                    "location": "table:software;row:" + family + str(i),
                }
            )
    source = {
        "schema": "island.opls-provider-records.v1",
        "family": "OPLS-AA",
        "graph_identity": graph.identity,
        "forms": {name: copy.deepcopy(FORMS[name]) for name in forms.values()},
        "atoms": atoms,
        "records": records,
        "pairs14": [list(p) for p in needed["interaction_14"]],
        "exclusions": [list(p) for p in exclusions],
    }
    manifest = {
        "schema": "island.opls-provider-bundle.v1",
        "provider": provider,
        "provider_version": "software-test",
        "source_version": "software-test",
        "profile_id": "not_a_scientific_profile",
        "family": "OPLS-AA",
        "graph_identity": graph.identity,
        "files": [],
        "license": {
            "identifier": "synthetic",
            "authorization": "authorized_local_use",
            "redistribution": "not_a_scientific_source",
        },
        "provenance": {
            "citation": "manufactured software controls only",
            "acquisition": "test generated",
            "source_uri": "test://never-scientific",
        },
        "conventions": conventions(),
        "synthetic": True,
    }
    if change:
        change(manifest, source)
    write(directory / "graph.json", data(graph))
    if native:
        # Create genuine GROMACS section syntax from software literals, not a
        # canonical JSON payload disguised as a native provider output.
        lines = ["[ defaults ]", "1 3 yes 0.5 0.5", "[ atomtypes ]"]
        lines += [f"{a['atom_type']} {a['mass']} 0 A 0.1 1" for a in source["atoms"]]
        lines += ["[ moleculetype ]", "Control 3", "[ atoms ]"]
        lines += [
            f"{a['site']} {a['atom_type']} 1 CTL {a['element']} {a['site']} 0 {a['mass']}"
            for a in source["atoms"]
        ]
        for family, section, func in (
            ("bond", "bonds", "1"),
            ("angle", "angles", "1"),
            ("torsion", "dihedrals", "9"),
        ):
            lines.append("[ " + section + " ]")
            for r in source["records"]:
                if r["family"] == family:
                    lines += [
                        " ".join([*map(str, r["sites"]), func, *v]) for v in r["values"]
                    ]
        lines += ["[ pairs ]", *[f"{a} {b} 1" for a, b in source["pairs14"]]]
        (directory / "parameters.itp").write_text("\n".join(lines) + "\n")
        write(
            directory / "omissions.json",
            [
                {"sites": r["sites"], "reason": r["reason"], "location": r["location"]}
                for r in source["records"]
                if r["family"] == "improper"
            ],
        )
        record_files = [
            ("parameters.itp", "parameters", "gromacs_itp"),
            ("omissions.json", "omissions", "improper_omissions_json"),
        ]
    else:
        if db:
            source = {"schema": "island.oplsaa-db-snapshot.v1", "entries": [source]}
        write(directory / "parameters.json", source)
        record_files = [
            (
                "parameters.json",
                "parameters",
                "oplsaa_db_snapshot" if db else "canonical_json",
            )
        ]
    for name, role, fmt in [
        ("graph.json", "graph", "chemical_graph_json"),
        *record_files,
    ]:
        p = directory / name
        manifest["files"].append(
            {
                "path": name,
                "sha256": pin(p),
                "size": p.stat().st_size,
                "role": role,
                "format": fmt,
            }
        )
    write(directory / "manifest.json", manifest)
    return {
        "provider": provider,
        "manifest": directory / "manifest.json",
        "manifest_sha256": pin(directory / "manifest.json"),
    }, spec


def audited(inputs, spec):
    registry = OPLSProviderRegistry()
    provisional = {**inputs, "provider": "ligpargen"}
    # Generate audit proposals by direct source parsing, not by assuming audit.
    m = read(inputs["manifest"])
    m["provider"] = "ligpargen"
    write(inputs["manifest"], m)
    provisional["manifest_sha256"] = pin(inputs["manifest"])
    result = registry.assign(spec, **provisional)
    m["provider"] = inputs["provider"]
    write(inputs["manifest"], m)
    inputs["manifest_sha256"] = pin(inputs["manifest"])
    report = result.report
    audit = {
        "schema": "island.opls-provider-audit.v1",
        "reviewer": "independent software reviewer control",
        "authority_basis": "caller pinned synthetic review; not scientific authority",
        "independent": True,
        "manifest_sha256": inputs["manifest_sha256"],
        "graph_identity": report["graph_identity"],
        "convention_identity": identity(m["conventions"]),
        "claim_hashes": [r["claim_sha256"] for r in report["record_claims"]],
        "atom_claim_hashes": [r["claim_sha256"] for r in report["atom_claims"]],
    }
    path = inputs["manifest"].parent / "audit.json"
    write(path, audit)
    return {**inputs, "audit": path, "audit_sha256": pin(path)}


@pytest.mark.parametrize("provider", PROVIDERS)
def test_provider_interfaces_and_complete_schema_software_only(tmp_path, provider):
    inputs, spec = control(tmp_path, provider=provider, db=provider == "oplsaa_db")
    if provider == "audited_bundle":
        inputs = audited(inputs, spec)
    result = OPLSProviderRegistry().assign(spec, **inputs)
    assert result.report["decision"] == (
        "audited_exact_graph"
        if provider == "audited_bundle"
        else "provisional_parameterization"
    ), result.report["blocking_reasons"]
    assert result.report["scientific_graphs_covered"] == 0
    assert all(r["status"] == "complete" for r in result.report["coverage"].values())
    assert {k: result.report[k] for k in FLAGS} == FLAGS
    assert result.canonical.payload["graph_identity"] == finite_graph(spec).identity


@pytest.mark.parametrize("provider", ["polypargen", "ligpargen"])
def test_native_provider_gromacs_sections_and_line_pointers(tmp_path, provider):
    inputs, spec = control(tmp_path, provider=provider, native=True)
    result = OPLSProviderRegistry().assign(spec, **inputs)
    assert result.report["decision"] == "provisional_parameterization", result.report[
        "blocking_reasons"
    ]
    assert all(
        "line:" in r["source"]["pointer"] or r["applicability"] == "not_applicable"
        for r in result.report["record_claims"]
    )


@pytest.mark.parametrize(
    "family", ["bond", "angle", "torsion", "improper", "charge", "nonbonded"]
)
def test_missing_terms_remain_explicit_no_zero_or_fallback(tmp_path, family):
    def change(m, source):
        source["records"] = [r for r in source["records"] if r["family"] != family]

    inputs, spec = control(tmp_path, change=change)
    result = OPLSProviderRegistry().assign(spec, **inputs)
    assert result.report["decision"] == "parameter_bundle_incomplete"
    assert any(r["family"] == family for r in result.report["missing_terms"])
    assert not any(r["family"] == family for r in result.canonical.payload["records"])


@pytest.mark.parametrize("target", ["parameters.json", "graph.json", "manifest.json"])
def test_source_hash_mutation_rejected(tmp_path, target):
    inputs, spec = control(tmp_path)
    path = tmp_path / target
    path.write_bytes(path.read_bytes() + b" ")
    result = OPLSProviderRegistry().assign(spec, **inputs)
    assert result.report["decision"] == "parameter_bundle_incomplete"
    assert "SHA256_mutation" in str(result.report["blocking_reasons"])


def test_missing_pin_never_opens_source(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("external IO attempted")

    monkeypatch.setattr("island.forcefields.opls_providers.registry.pinned", forbidden)
    r = OPLSProviderRegistry().assign(
        settings(), provider="ligpargen", manifest=tmp_path / "untrusted.json"
    )
    assert r.report["decision"] == "source_bundle_required"
    assert r.report["source_files_inspected"] == 0


def test_authorization_checked_before_parameter_sources(tmp_path, monkeypatch):
    inputs, spec = control(
        tmp_path,
        change=lambda m, s: m["license"].update(authorization="not_authorized"),
    )
    (tmp_path / "parameters.json").unlink()
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert r.report["decision"] == "source_bundle_required"
    assert r.report["source_files_inspected"] == 0


@pytest.mark.parametrize("provider", ["GAFF", "PCFF", ["ligpargen", "polypargen"]])
def test_one_provider_only_no_family_fallback(provider):
    with pytest.raises(ValueError, match="one_explicit_OPLS_provider"):
        OPLSProviderRegistry().assign(settings(), provider=provider)


def test_conflicting_provider_manifest_rejected(tmp_path):
    inputs, spec = control(tmp_path)
    r = OPLSProviderRegistry().assign(spec, **{**inputs, "provider": "polypargen"})
    assert "conflicting_provider" in str(r.report["blocking_reasons"])


def test_duplicate_and_conflicting_rows_no_tie_break(tmp_path):
    def change(m, s):
        r = copy.deepcopy(s["records"][0])
        r["id"] += "duplicate"
        r["values"][0][-1] = "2"
        s["records"].append(r)

    inputs, spec = control(tmp_path, change=change)
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "conflicting_or_duplicate_parameter" in str(r.report["blocking_reasons"])


@pytest.mark.parametrize("field", CONVENTIONS)
def test_complete_convention_contract_required(tmp_path, field):
    inputs, spec = control(tmp_path, change=lambda m, s: m["conventions"].pop(field))
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert r.report["decision"] == "parameter_bundle_incomplete"
    assert r.canonical is None


@pytest.mark.parametrize(
    "field,value",
    [
        ("energy_unit", "kcal/mol"),
        ("length_unit", "angstrom"),
        ("nonbonded", "lj9_6"),
        ("mixing_rule", "arithmetic"),
        ("torsion_phase", "unknown"),
        ("improper_ordering", "unknown"),
    ],
)
def test_no_convention_coercion(tmp_path, field, value):
    inputs, spec = control(
        tmp_path, change=lambda m, s: m["conventions"].update({field: value})
    )
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "unsupported_convention:" + field in r.report["blocking_reasons"]
    assert r.report["decision"] == "parameter_bundle_incomplete"


def test_source_form_units_cannot_be_reinterpreted(tmp_path):
    inputs, spec = control(
        tmp_path,
        change=lambda m, s: s["forms"]["bond_harmonic"]["units"].__setitem__(
            0, "angstrom"
        ),
    )
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "equation_unit_ordering_mismatch" in str(r.report["blocking_reasons"])


def test_charge_conservation_without_renormalization(tmp_path):
    def change(m, s):
        next(r for r in s["records"] if r["family"] == "charge")["values"] = [["0.1"]]

    inputs, spec = control(tmp_path, change=change)
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "charge_conservation_failed" in r.report["blocking_reasons"]
    assert r.report["decision"] == "parameter_bundle_incomplete"


@pytest.mark.parametrize("name", list(benchmarks()))
def test_final_capped_graphs_deterministic_not_scientific_coverage(name):
    s = benchmarks()[name]
    graph = finite_graph(s)
    assert graph == finite_graph(copy.deepcopy(s))
    assert not graph.blockers
    assert all(a.GetAtomicNum() > 0 for a in molecule(graph).GetAtoms())
    r = OPLSProviderRegistry().assign(s, provider="ligpargen")
    assert r.report["decision"] == "source_bundle_required"
    assert r.report["scientific_graphs_covered"] == 0
    assert all(c["covered"] == 0 for c in r.report["coverage"].values())


@pytest.mark.parametrize("field", list(settings()))
def test_finite_settings_never_implicit(field):
    s = settings()
    s.pop(field)
    with pytest.raises(ValueError, match="explicit_finite_chain_settings"):
        finite_graph(s)


def test_end_group_and_repeat_count_changes_require_reparameterization(tmp_path):
    inputs, spec = control(tmp_path)
    for change in ({"end_groups": ["*O", "*C"]}, {"repeat_count": 3}):
        changed = {**spec, **change}
        assert finite_graph(changed).identity != finite_graph(spec).identity
        r = OPLSProviderRegistry().assign(changed, **inputs)
        assert "topology_changed_provider_reassignment" in str(
            r.report["blocking_reasons"]
        )


def test_stereochemistry_identity_and_unassigned_rejection():
    s = benchmarks()["PLA"]
    other = {**s, "psmiles": s["psmiles"].replace("@", "@@")}
    assert finite_graph(s).identity != finite_graph(other).identity
    with pytest.raises(ValueError, match="achiral_claim"):
        finite_graph({**s, "tacticity": "not_applicable", "stereochemistry": "achiral"})


def test_atom_order_mutation_rejected_even_with_resigned_file_pin(tmp_path):
    inputs, spec = control(tmp_path)
    p = tmp_path / "graph.json"
    g = read(p)
    g["atom_order"].reverse()
    write(p, g)
    m = read(inputs["manifest"])
    m["files"][0].update(sha256=pin(p), size=p.stat().st_size)
    write(inputs["manifest"], m)
    inputs["manifest_sha256"] = pin(inputs["manifest"])
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "graph_bytes_atom_order_stereo_or_caps" in str(r.report["blocking_reasons"])


@pytest.mark.parametrize("engine", ["gromacs", "lammps", "openmm"])
def test_metadata_preserving_engine_round_trip_no_execution(tmp_path, engine):
    s = {**settings(), "intended_export_engine": engine}
    inputs, s = control(tmp_path / "sources", spec=s)
    p = tmp_path / ("export." + engine)
    metadata = export_parameters(s, engine=engine, destination=p, **inputs)
    value = read_export(
        p,
        expected_sha256=metadata["sha256"],
        expected_source_hashes=metadata["source_hashes"],
        expected_graph_identity=metadata["graph_identity"],
    )
    assert value["receipt"]["non_production"] is True
    assert value["receipt"]["source_hashes"] == metadata["source_hashes"]
    assert all(checksum in p.read_text() for checksum in metadata["source_hashes"])
    assert value["canonical"]["graph_identity"] == finite_graph(s).identity
    assert value["canonical"]["synthetic"] is True
    if engine == "openmm":
        assert 'loadableForceField="false"' in p.read_text()
    p.write_text(p.read_text() + "altered")
    with pytest.raises(ValueError, match="export_body_or_metadata_mutation"):
        read_export(
            p,
            expected_sha256=pin(p),
            expected_source_hashes=metadata["source_hashes"],
            expected_graph_identity=metadata["graph_identity"],
        )


def test_export_rehashes_sources_and_rejects_checkout(tmp_path):
    inputs, spec = control(tmp_path)
    with pytest.raises(ValueError, match="outside_Git_checkout"):
        export_parameters(
            spec, engine="gromacs", destination=Path.cwd() / "forbidden.itp", **inputs
        )
    (tmp_path / "parameters.json").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="source_bundle_required_no_export"):
        export_parameters(
            spec, engine="gromacs", destination=tmp_path / "out.itp", **inputs
        )


def test_missing_improper_omissions_in_native_source_blocked(tmp_path):
    inputs, spec = control(tmp_path, native=True)
    m = read(inputs["manifest"])
    m["files"] = [f for f in m["files"] if f["role"] != "omissions"]
    write(inputs["manifest"], m)
    inputs["manifest_sha256"] = pin(inputs["manifest"])
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert r.report["decision"] == "parameter_bundle_incomplete"
    assert any(t["family"] == "improper" for t in r.report["missing_terms"])


def test_native_macros_and_includes_never_executed(tmp_path):
    inputs, spec = control(tmp_path, native=True)
    p = tmp_path / "parameters.itp"
    p.write_text('#include "unknown.itp"\n' + p.read_text())
    m = read(inputs["manifest"])
    f = next(f for f in m["files"] if f["role"] == "parameters")
    f.update(sha256=pin(p), size=p.stat().st_size)
    write(inputs["manifest"], m)
    inputs["manifest_sha256"] = pin(inputs["manifest"])
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "preprocessor_include_or_macro" in str(r.report["blocking_reasons"])


def test_audited_provider_requires_independent_exact_scope_audit(tmp_path):
    inputs, spec = control(tmp_path, provider="audited_bundle")
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert r.report["decision"] == "source_bundle_required"
    assert r.canonical is None
    inputs = audited(inputs, spec)
    a = read(inputs["audit"])
    a["claim_hashes"] = []
    write(inputs["audit"], a)
    inputs["audit_sha256"] = pin(inputs["audit"])
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "audit_claim_coverage_mismatch" in str(r.report["blocking_reasons"])


def test_opaque_attachment_hash_pinned_without_parsing(tmp_path):
    inputs, spec = control(tmp_path)
    p = tmp_path / "licensed.bin"
    p.write_bytes(bytes(range(256)))
    m = read(inputs["manifest"])
    m["files"].append(
        {
            "path": p.name,
            "sha256": pin(p),
            "size": p.stat().st_size,
            "role": "opaque",
            "format": "opaque",
        }
    )
    write(inputs["manifest"], m)
    inputs["manifest_sha256"] = pin(inputs["manifest"])
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert pin(p) in r.report["source_hashes"]
    p.write_bytes(b"mutation")
    assert (
        OPLSProviderRegistry().assign(spec, **inputs).report["decision"]
        == "parameter_bundle_incomplete"
    )


def test_no_source_benchmark_deterministic_replay_and_resigned_receipt_rejection(
    tmp_path,
):
    root = tmp_path / "report"
    result = emit(Path.cwd(), root)
    assert replay(root) == result["result"]
    assert result["result"]["status_counts"]["source_bundle_required"] == 20
    assert result["result"]["source_rows_inspected"] == 0
    assert result["result"]["scientific_graphs_covered"] == 0
    r = read(root / "result.json")
    r["scope"]["runtime_registration"] = True
    write(root / "result.json", signed({k: v for k, v in r.items() if k != "identity"}))
    with pytest.raises(ValueError, match="semantic_mutation"):
        replay(root)


def test_external_deterministic_and_relocated_replay_requires_originals(tmp_path):
    inputs, spec = control(tmp_path / "external")
    request = {
        "schema": "island.opls-polymer-request.v1",
        "provider": "ligpargen",
        "settings": spec,
        "graph_identity": finite_graph(spec).identity,
    }
    p = tmp_path / "external/request.json"
    write(p, request)
    external = {k: v for k, v in inputs.items() if k != "provider"} | {
        "request": p,
        "request_sha256": pin(p),
    }
    root = tmp_path / "reports"
    result = emit(Path.cwd(), root, **external)
    assert replay(root, result["snapshot_sha256"], **external) == result["result"]
    moved = tmp_path / "relocated"
    shutil.copytree(tmp_path / "external", moved)
    relocated = {
        k: (moved / v.name if isinstance(v, Path) else v) for k, v in external.items()
    }
    assert replay(root, result["snapshot_sha256"], **relocated) == result["result"]
    assert "values" not in json.dumps(read(root / "provider-results.json"))
    with pytest.raises(ValueError, match="originals_and_pins"):
        replay(root, result["snapshot_sha256"])
    (moved / "parameters.json").write_bytes(b"mutated")
    with pytest.raises(ValueError, match="semantic_mutation|unavailable_or_mutated"):
        replay(root, result["snapshot_sha256"], **relocated)


def test_no_runtime_registry_changes():
    from island.forcefields.opls_registry.registry import canonical_registry

    before = canonical_registry()
    OPLSProviderRegistry().assign(settings(), provider="ligpargen")
    assert canonical_registry() == before


def test_real_provider_bundle_not_supplied():
    pytest.skip(
        "No real user-supplied authorized provider bundle/pins; synthetic controls are not scientific data."
    )


@pytest.mark.parametrize(
    "kind", ["atom", "pair14", "location", "empty_values", "improper_order"]
)
def test_explicit_coverage_and_location_failures(tmp_path, kind):
    def change(m, source):
        if kind == "atom":
            source["atoms"].pop()
        elif kind == "pair14":
            source["pairs14"].pop()
        elif kind == "location":
            source["records"][0]["location"] = ""
        elif kind == "empty_values":
            source["records"][0]["values"] = []
        else:
            row = next(r for r in source["records"] if r["family"] == "improper")
            row["sites"][0], row["sites"][1] = row["sites"][1], row["sites"][0]

    inputs, spec = control(tmp_path, change=change)
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert r.report["decision"] == "parameter_bundle_incomplete"
    assert r.report["blocking_reasons"]
    if kind == "atom":
        assert r.report["coverage"]["atom_typing"]["status"] == "incomplete"
    if kind == "pair14":
        assert r.report["coverage"]["interaction_14"]["status"] == "incomplete"


def test_provenance_only_without_parameters(tmp_path):
    inputs, spec = control(tmp_path)
    m = read(inputs["manifest"])
    m["files"] = [f for f in m["files"] if f["role"] == "graph"]
    write(inputs["manifest"], m)
    inputs["manifest_sha256"] = pin(inputs["manifest"])
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert r.report["decision"] == "provenance_only"
    assert r.canonical is None


def test_symlink_and_path_substitution_rejected(tmp_path):
    inputs, spec = control(tmp_path)
    original = tmp_path / "parameters.json"
    other = tmp_path / "original.json"
    original.rename(other)
    original.symlink_to(other)
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert "source_symlink_or_path_substitution" in str(r.report["blocking_reasons"])


def test_strict_receipt_reconstructs_semantics_even_after_resigning(
    tmp_path, monkeypatch
):
    from scripts import validate_phase_4j90_receipt as validator

    root = tmp_path / "evidence"
    emit(Path.cwd(), root)
    # Software unit control for the outer seal. Actual terminal verification and
    # complete historical preservation are independently run for the delivery.
    monkeypatch.setattr(
        validator, "verification", lambda repo, root: {"software_control": True}
    )
    monkeypatch.setattr(
        validator, "historical", lambda repo: {"software-control": "unchanged"}
    )
    validator.seal(Path.cwd(), root)
    assert validator.validate(Path.cwd(), root)["status"] == "passed"
    r = read(root / "receipt.json")
    r["scope"]["accepted_coefficients"] = 1
    write(
        root / "receipt.json", signed({k: v for k, v in r.items() if k != "identity"})
    )
    with pytest.raises(ValueError, match="scope_or_readiness"):
        validator.validate(Path.cwd(), root)


@pytest.mark.parametrize(
    "form,values", [("torsion_rb", [["1"] * 6]), ("torsion_opls", [["1"] * 4])]
)
def test_torsion_forms_retained_without_reinterpretation(tmp_path, form, values):
    def change(m, s):
        s["forms"][form] = copy.deepcopy(FORMS[form])
        for row in s["records"]:
            if row["family"] == "torsion":
                row.update(form=form, values=copy.deepcopy(values))

    inputs, spec = control(tmp_path, change=change)
    r = OPLSProviderRegistry().assign(spec, **inputs)
    assert r.report["decision"] == "provisional_parameterization"
    p = tmp_path / "export.itp"
    export_parameters(spec, engine="gromacs", destination=p, **inputs)
    if form == "torsion_rb":
        spec = {**spec, "intended_export_engine": "lammps"}
        # Settings change must rebind the graph manifest before another export.
        inputs, spec = control(tmp_path / "lammps", spec=spec, change=change)
        with pytest.raises(ValueError, match="no_torsion_reinterpretation"):
            export_parameters(
                spec, engine="lammps", destination=tmp_path / "coeff", **inputs
            )
