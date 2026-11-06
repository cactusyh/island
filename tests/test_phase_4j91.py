"""J91 software-only controls. The optional real integration never fabricates data."""

import copy
import os
import shutil
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from island.forcefields.assignment_core.schema import identity
from island.forcefields.opls_providers.export import read_export
from island.forcefields.opls_providers.intake91 import (
    CONVERSION,
    DEFINITIONS,
    export,
    intake,
    pair_claims,
)
from island.forcefields.opls_providers.schema import FLAGS, FORMS
from scripts.phase_4j73_source_intake import pin, read, signed, write
from scripts.phase_4j91_intake import emit, historical, replay
from tests.test_phase_4j90 import audited, control, settings


def bundle(root, *, change=None, native=False, engine="gromacs", conversion=False):
    spec = settings()
    spec["intended_export_engine"] = engine
    inputs, spec = control(
        root, provider="polypargen", spec=spec, change=change, native=native
    )
    m = read(inputs["manifest"])
    if native:
        from island.forcefields.opls_providers import finite_graph
        from island.forcefields.opls_providers.registry import obligations

        _, exclusions = obligations(finite_graph(spec))
        path = root / "parameters.itp"
        with path.open("a") as f:
            f.write(
                "[ exclusions ]\n" + "\n".join(f"{a} {b}" for a, b in exclusions) + "\n"
            )
        for entry in m["files"]:
            if entry["path"] == path.name:
                entry.update(sha256=pin(path), size=path.stat().st_size)
        write(inputs["manifest"], m)
        inputs["manifest_sha256"] = pin(inputs["manifest"])
    inputs = audited(inputs, spec)
    source_audit = read(root / "audit.json")
    m = read(inputs["manifest"])
    q = {
        "schema": "island.j91.request.v1",
        "provider": "polypargen",
        "settings": spec,
        "graph_identity": m["graph_identity"],
        "source_authority": "synthetic software control only",
        "self_contained_finite_molecule": True,
        "expected_conventions": m["conventions"],
        "expected_forms": {
            family: copy.deepcopy(FORMS[form])
            for family, form in {
                "bond": "bond_harmonic",
                "angle": "angle_harmonic",
                "torsion": "torsion_periodic",
                "improper": "improper_harmonic",
                "charge": "charge",
                "nonbonded": "lj12_6",
            }.items()
        },
        "nonbonded_definitions": DEFINITIONS,
    }
    write(root / "request.json", q)
    f = next(f for f in m["files"] if f["role"] == "parameters")
    if native:
        from island.forcefields.opls_providers.parsers import gromacs

        _, _, pairs, exclusions = gromacs(
            (root / f["path"]).read_bytes(), f, q["graph_identity"], m["conventions"]
        )
    else:
        s = read(root / f["path"])
        pairs, exclusions = s["pairs14"], s["exclusions"]
    pairs_meta, exclusions_meta = pair_claims(
        (root / f["path"]).read_bytes(), f, pairs, exclusions
    )
    write(
        root / "audit.json",
        {
            "schema": "island.j91.audit.v1",
            "request_sha256": pin(root / "request.json"),
            "source_audit": source_audit,
            "pair_claim_hashes": [identity(r) for r in pairs_meta],
            "exclusion_claim_hashes": [identity(r) for r in exclusions_meta],
            "export_conversions": [CONVERSION] if conversion else [],
        },
    )
    return inputs_for(root)


def inputs_for(root):
    return {
        "source_dir": root,
        **{
            n + "_sha256": pin(root / (n + ".json"))
            for n in ("request", "manifest", "audit")
        },
    }


def edit_control(inputs, name, mutate):
    root = inputs["source_dir"]
    obj = read(root / (name + ".json"))
    mutate(obj)
    write(root / (name + ".json"), obj)
    if name != "audit":
        a = read(root / "audit.json")
        if name == "request":
            a["request_sha256"] = pin(root / "request.json")
        if name == "manifest":
            a["source_audit"]["manifest_sha256"] = pin(root / "manifest.json")
        write(root / "audit.json", a)
    return inputs_for(root)


@pytest.mark.parametrize("native", [False, True])
def test_complete_source_contract_software_only(tmp_path, native):
    inputs = bundle(tmp_path / "external", native=native)
    result = intake(**inputs)
    r = result.report
    assert r["decision"] == "audited_exact_graph", r["blocking_reasons"]
    assert r["scientific_graphs_covered"] == 0 and r["qualified_graph"] is None
    assert r["synthetic"] and r["source_files_inspected"] == (3 if native else 2)
    assert r["coverage"]["exclusion"]["status"] == "complete"
    assert not r["missing_terms"]
    assert all(
        t["status"] == "source_record_present" for t in r["interaction_instances"]
    )
    assert r["source_rows_inspected"] > 0
    assert all(r[k] == v for k, v in FLAGS.items())


def test_missing_source_does_not_construct_fake_graph_or_inspect(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("source inspected")

    monkeypatch.setattr("island.forcefields.opls_providers.intake91.pinned", forbidden)
    r = intake().report
    assert r["decision"] == "source_bundle_required"
    assert r["requested_graph"] is r["provider"] is r["qualified_graph"] is None
    assert r["coverage"] == {} and r["missing_terms"] == []
    assert r["source_files_inspected"] == r["source_rows_inspected"] == 0


@pytest.mark.parametrize(
    "missing", ["source_dir", "request_sha256", "manifest_sha256", "audit_sha256"]
)
def test_missing_input_pin_never_opens_controls(tmp_path, monkeypatch, missing):
    inputs = bundle(tmp_path / "external")
    del inputs[missing]

    def forbidden(*args):
        raise AssertionError("read before complete pins")

    monkeypatch.setattr("island.forcefields.opls_providers.intake91.pinned", forbidden)
    assert intake(**inputs).report["decision"] == "source_bundle_required"


@pytest.mark.parametrize("which", ["request", "manifest", "audit"])
def test_missing_control_never_inspects_sources(tmp_path, which):
    inputs = bundle(tmp_path / "external")
    (inputs["source_dir"] / (which + ".json")).unlink()
    r = intake(**inputs).report
    assert (
        r["decision"] == "source_bundle_required" and r["source_files_inspected"] == 0
    )


@pytest.mark.parametrize(
    "name,key",
    [
        ("request", "source_authority"),
        ("request", "settings"),
        ("manifest", "provider_version"),
        ("manifest", "license"),
        ("audit", "source_audit"),
    ],
)
def test_incomplete_metadata_no_parameter_inspection(tmp_path, monkeypatch, name, key):
    inputs = bundle(tmp_path / "external")
    inputs = edit_control(inputs, name, lambda v: v.pop(key))
    from island.forcefields.opls_providers.intake91 import pinned

    def guard(path, sha):
        assert Path(path).name in ("request.json", "manifest.json", "audit.json")
        return pinned(path, sha)

    monkeypatch.setattr("island.forcefields.opls_providers.intake91.pinned", guard)
    r = intake(**inputs).report
    assert (
        r["decision"] == "source_bundle_required" and r["source_files_inspected"] == 0
    )


@pytest.mark.parametrize(
    "path",
    ["parameters.json", "graph.json", "request.json", "manifest.json", "audit.json"],
)
def test_byte_mutation_and_replay_tamper(tmp_path, path):
    inputs = bundle(tmp_path / "external")
    out = tmp_path / "receipt"
    original = emit(Path.cwd(), out, **inputs)
    assert replay(out, original["snapshot_sha256"], **inputs) == original["result"]
    with (inputs["source_dir"] / path).open("a") as f:
        f.write(" ")
    assert intake(**inputs).report["decision"] not in (
        "audited_exact_graph",
        "provisional_parameterization",
    )
    with pytest.raises(ValueError, match="semantic_mutation"):
        replay(out, original["snapshot_sha256"], **inputs)


@pytest.mark.parametrize(
    "change", ["unlisted", "missing", "symlink", "opaque_mutation"]
)
def test_complete_directory_census_and_opaque_pins(tmp_path, change):
    inputs = bundle(tmp_path / "external")
    root = inputs["source_dir"]
    if change == "unlisted":
        (root / "secret.bin").write_bytes(b"do not inspect")
    elif change == "missing":
        (root / "parameters.json").unlink()
    elif change == "symlink":
        (root / "parameters.json").rename(tmp_path / "outside.json")
        (root / "parameters.json").symlink_to(tmp_path / "outside.json")
    else:
        f = root / "opaque.bin"
        f.write_bytes(b"private original")
        inputs = edit_control(
            inputs,
            "manifest",
            lambda m: m["files"].append(
                {
                    "path": f.name,
                    "sha256": pin(f),
                    "size": f.stat().st_size,
                    "role": "opaque",
                    "format": "opaque",
                }
            ),
        )
        f.write_bytes(b"private modified")
    r = intake(**inputs).report
    assert r["decision"] == "parameter_bundle_incomplete"
    if change != "opaque_mutation":
        assert r["source_files_inspected"] == 0


@pytest.mark.parametrize(
    "key,value",
    [
        ("energy_unit", "kcal/mol"),
        ("length_unit", "Angstrom"),
        ("mixing_rule", "arithmetic_sigma_geometric_epsilon"),
        ("scale14_lj", "0.25"),
        ("scale14_coulomb", "0.25"),
        ("torsion_normalization", "missing_half_factor"),
        ("improper_ordering", "center_first"),
        ("nonbonded", "lj9_6"),
    ],
)
def test_requested_convention_mismatch(tmp_path, key, value):
    inputs = bundle(tmp_path / "external")
    inputs = edit_control(
        inputs, "request", lambda q: q["expected_conventions"].update({key: value})
    )
    r = intake(**inputs).report
    assert r["decision"] == "parameter_bundle_incomplete"
    assert "requested_convention_mismatch:" + key in r["convention_blockers"]
    assert r["source_files_inspected"] == 0


@pytest.mark.parametrize(
    "forms",
    [
        None,
        [],
        {},
        "bad",
        {"bond": {}},
        {
            k: None
            for k in ("bond", "angle", "torsion", "improper", "charge", "nonbonded")
        },
    ],
)
def test_malformed_forms_structured_blocker(tmp_path, forms):
    inputs = bundle(tmp_path / "external")
    inputs = edit_control(inputs, "request", lambda q: q.update(expected_forms=forms))
    r = intake(**inputs).report
    assert (
        r["decision"] == "source_bundle_required" and r["source_files_inspected"] == 0
    )


@pytest.mark.parametrize(
    "family",
    [
        "bond",
        "angle",
        "torsion",
        "improper",
        "charge",
        "nonbonded",
        "exclusion",
        "interaction_14",
    ],
)
def test_missing_interaction_no_zero_or_derived_exclusion(tmp_path, family):
    def change(m, source):
        if family in ("exclusion", "interaction_14"):
            source["exclusions" if family == "exclusion" else "pairs14"].pop()
        else:
            source["records"].remove(
                next(r for r in source["records"] if r["family"] == family)
            )

    inputs = bundle(tmp_path / "external", change=change)
    r = intake(**inputs).report
    assert r["decision"] == "parameter_bundle_incomplete", r["blocking_reasons"]
    assert any(t["family"] == family for t in r["missing_terms"])
    assert r["coverage"][family]["covered"] < r["coverage"][family]["required"]


@pytest.mark.parametrize(
    "mutation",
    [
        "identity",
        "caps",
        "repeat_count",
        "protonation",
        "charge",
        "stereo",
        "tacticity",
        "multiple_components",
        "radical",
    ],
)
def test_exact_finite_graph_scope_mismatch(tmp_path, mutation):
    inputs = bundle(tmp_path / "external")

    def change(q):
        if mutation == "identity":
            q["graph_identity"] = "0" * 64
        elif mutation == "caps":
            q["settings"]["end_groups"] = ["*O", "*C"]
        else:
            field, val = {
                "repeat_count": ("repeat_count", 3),
                "protonation": ("protonation_state", "guess"),
                "charge": ("charge_state", 1),
                "stereo": ("stereochemistry", "explicit"),
                "tacticity": ("tacticity", "isotactic"),
                "multiple_components": ("psmiles", "*CC*.C"),
                "radical": ("psmiles", "*[CH]C*"),
            }[mutation]
            q["settings"][field] = val

    inputs = edit_control(inputs, "request", change)
    r = intake(**inputs).report
    assert (
        r["decision"] == "source_bundle_required" and r["source_files_inspected"] == 0
    )


def test_source_graph_order_mutation(tmp_path):
    inputs = bundle(tmp_path / "external")
    root = inputs["source_dir"]
    graph = read(root / "graph.json")
    graph["atom_order"] = list(reversed(graph["atom_order"]))
    write(root / "graph.json", graph)

    def modify(m):
        f = next(f for f in m["files"] if f["role"] == "graph")
        f.update(
            sha256=pin(root / "graph.json"), size=(root / "graph.json").stat().st_size
        )

    inputs = edit_control(inputs, "manifest", modify)
    r = intake(**inputs).report
    assert (
        "source_graph_identity_order_stereo_caps_mismatch" in r["blocking_reasons"][0]
    )


def test_charge_disagreement(tmp_path):
    def change(m, s):
        next(r for r in s["records"] if r["family"] == "charge")["values"] = [["1"]]

    r = intake(**bundle(tmp_path / "external", change=change)).report
    assert "charge_conservation_failed" in r["blocking_reasons"]


@pytest.mark.parametrize(
    "kind",
    [
        "provider",
        "independent",
        "authorization",
        "audit_scope",
        "pair_claim",
        "exclusion_claim",
    ],
)
def test_authority_and_single_provider_guards(tmp_path, kind):
    inputs = bundle(tmp_path / "external")
    if kind == "provider":
        inputs = edit_control(
            inputs, "request", lambda q: q.update(provider=["polypargen", "ligpargen"])
        )
    elif kind == "authorization":
        inputs = edit_control(
            inputs, "manifest", lambda m: m["license"].update(authorization="unknown")
        )
    else:

        def change(a):
            if kind == "independent":
                a["source_audit"]["independent"] = False
            elif kind == "audit_scope":
                a["source_audit"]["graph_identity"] = "0" * 64
            else:
                a[
                    "pair_claim_hashes"
                    if kind == "pair_claim"
                    else "exclusion_claim_hashes"
                ].pop()

        inputs = edit_control(inputs, "audit", change)
    r = intake(**inputs).report
    assert r["decision"] in ("source_bundle_required", "parameter_bundle_incomplete")
    assert r["scientific_graphs_covered"] == 0


@pytest.mark.parametrize("engine", ["gromacs", "lammps", "openmm"])
def test_fragment_exports_roundtrip_no_engine_execution(tmp_path, engine):
    inputs = bundle(tmp_path / "external", engine=engine, conversion=engine == "lammps")
    destination = tmp_path / (
        "export." + {"gromacs": "itp", "lammps": "coeff", "openmm": "xml"}[engine]
    )
    info = export(destination=destination, engine=engine, **inputs)
    payload = read_export(
        destination,
        expected_sha256=info["sha256"],
        expected_source_hashes=info["source_hashes"],
        expected_graph_identity=info["graph_identity"],
    )
    assert payload["receipt"]["non_production"]
    assert payload["canonical"]["synthetic"]
    if engine == "openmm":
        assert ET.parse(destination).getroot().attrib["loadableForceField"] == "false"
    if engine == "lammps":
        assert payload["report"]["export_conversions"] == [CONVERSION]


def test_lammps_conversion_requires_independent_audit(tmp_path):
    inputs = bundle(tmp_path / "external", engine="lammps")
    with pytest.raises(ValueError, match="conversion_audit_required"):
        export(destination=tmp_path / "output", engine="lammps", **inputs)
    assert not (tmp_path / "output").exists()


def test_exports_rehash_and_remain_outside_worktree(tmp_path):
    inputs = bundle(tmp_path / "external")
    with pytest.raises(ValueError, match="outside_Git_worktree"):
        export(destination=Path.cwd() / "forbidden.itp", engine="gromacs", **inputs)
    (inputs["source_dir"] / "parameters.json").write_bytes(b"changed")
    with pytest.raises(ValueError, match="incomplete_result_export_forbidden"):
        export(destination=tmp_path / "output", engine="gromacs", **inputs)


def test_deterministic_relocated_external_replay_and_missing_originals(tmp_path):
    inputs = bundle(tmp_path / "external")
    out = tmp_path / "reports"
    emitted = emit(Path.cwd(), out, **inputs)
    assert emitted["result"] == replay(out, emitted["snapshot_sha256"], **inputs)
    relocated = tmp_path / "relocated"
    shutil.copytree(inputs["source_dir"], relocated)
    shutil.rmtree(inputs["source_dir"])
    with pytest.raises(ValueError, match="semantic_mutation"):
        replay(out, emitted["snapshot_sha256"], **inputs)
    assert (
        replay(out, emitted["snapshot_sha256"], **{**inputs, "source_dir": relocated})
        == emitted["result"]
    )
    with pytest.raises(ValueError, match="originals_and_pins"):
        replay(out, emitted["snapshot_sha256"])
    serialized = "".join(p.read_text() for p in out.iterdir())
    assert '"values"' not in serialized and '"mass"' not in serialized


def test_missing_source_replay_and_resigned_semantic_rejection(tmp_path):
    out = tmp_path / "receipt"
    emitted = emit(Path.cwd(), out)
    assert replay(out) == emitted["result"]
    result = read(out / "result.json")
    result["accepted_coefficients"] = 1
    from scripts.phase_4j73_source_intake import signed

    write(
        out / "result.json", signed({k: v for k, v in result.items() if k != "_seal"})
    )
    with pytest.raises(ValueError, match="semantic_mutation"):
        replay(out)


def test_historical_files_and_main_untouched():
    assert len(historical(Path.cwd())) >= 12832


def test_real_external_intake_only_when_explicitly_supplied(tmp_path):
    root = os.environ.get("ISLAND_OPLS_SOURCE_DIR")
    if root is None:
        pytest.skip(
            "No explicitly supplied external OPLS source directory; no real scientific coverage."
        )
    pins = {
        n + "_sha256": os.environ.get("ISLAND_OPLS_" + n.upper() + "_SHA256")
        for n in ("request", "manifest", "audit")
    }
    result = intake(root, **pins).report
    assert not result["synthetic"], (
        "Synthetic controls cannot exercise the real-source path"
    )
    assert all(result[k] == v for k, v in FLAGS.items())
    out = tmp_path / "real-report"
    e = emit(Path.cwd(), out, source_dir=root, **pins)
    assert replay(out, e["snapshot_sha256"], source_dir=root, **pins) == e["result"]
    # Missing pins remain blocked without opening originals; no scientific success asserted.
    assert result["decision"] in (
        "source_bundle_required",
        "provenance_only",
        "parameter_bundle_incomplete",
        "provisional_parameterization",
        "audited_exact_graph",
    )


def test_strict_receipt_reconstructs_semantics_even_after_resigning(
    tmp_path, monkeypatch
):
    from scripts import validate_phase_4j91_receipt as validator

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
