"""Independent parameter claims, immutable source identity and backend observations."""

import json
import os
import shutil
import subprocess
import sys
from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import pytest

from island.chemistry import from_smiles
from island.forcefield_policy import (
    ClaimRequirements,
    OPLSAssets,
    build_openmm,
    check_openmm,
    export_lammps,
    prepare_polymer_forcefield,
)
from island.forcefield_policy.claims import CONVENTIONS, assess
from island.forcefield_policy.evidence import (
    reconstruct,
    record_static_comparison,
    save,
)
from island.forcefields.opls_registry import ForceFieldRequest, OPLSResolutionError
from island.forcefields.opls_registry.registry import ARCHIVE_SHA256, NATIVE
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import canonical, requirements
from island.opls_resolution import ConventionSources
from scripts.phase_4j50 import fixture

ROOT = Path(__file__).resolve().parents[1]
XML = (
    ROOT.parent
    / "island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml"
)
ARCHIVE = ROOT / "island_ffref/pcff/PDMS_benchmark-main.zip"
PAPER = ROOT / "island_ffref/pcff/jp4c08471.pdf"
SUPPLEMENT = ROOT / "island_ffref/pcff/jp4c08471_si_001.pdf"


@pytest.fixture(scope="module")
def assets():
    assert XML.is_file() and ARCHIVE.is_file(), (
        "J55 validation requires the installed pinned assets"
    )
    return OPLSAssets(
        ARCHIVE, ARCHIVE_SHA256, XML, sha256(XML.read_bytes()).hexdigest()
    )


@pytest.fixture(scope="module")
def sources():
    return ConventionSources(PAPER, SUPPLEMENT)


@pytest.fixture(scope="module")
def native(assets):
    system = from_smiles("CC", random_seed=2026)
    return system, prepare_polymer_forcefield(system, forcefield="opls", assets=assets)


@pytest.fixture(scope="module", params=[0, 2, 4, 6])
def pdms(request, assets, sources):
    s = fixture(request.param)
    return s, prepare_polymer_forcefield(
        s, forcefield="opls", assets=assets, sources=sources
    )


def diagnostic(path, system, prepared, *, failed=True):
    # Controlled external observations test only the import/status policy.
    # They are never represented as a real engine run or equivalence proof.
    value = {
        "assignment_identity": prepared.payload["assignment_identity"],
        "input_identity": prepared._system_identity,
        "reference_engine": "OpenMM",
        "observed_engine": "LAMMPS",
        "scope": "static_components",
        "test_control": True,
        "measurements": [
            {
                "observable": "total",
                "kind": "energy",
                "reference": -1.0,
                "observed": -0.99 if failed else -1.0,
                "bound": 3e-4,
            }
        ],
    }
    path.write_text(canonical(value))
    return sha256(path.read_bytes()).hexdigest()


def test_complete_without_equivalence(native):
    _, p = native
    r = p.payload
    assert r["parameter_completeness"] == "complete" and r["parameter_complete"]
    assert r["typing_status"] == r["charge_status"] == "complete"
    assert r["cross_engine_equivalence"] == "not_claimed"
    assert r["lammps_runtime_status"] == r["openmm_build_status"] == "not_run"
    assert r["scientific_validation"] == r["production_readiness"] == "not_established"
    assert r["selected_profile"] == NATIVE
    assert not r["missing_terms"] and not r["missing_conventions"]
    p.require(ClaimRequirements(parameter_complete=True))


@pytest.mark.parametrize("failed", [True, False])
def test_optional_diagnostic_never_invalidates_complete_assignment(
    native, tmp_path, failed
):
    s, p = native
    file = tmp_path / "comparison.json"
    q = record_static_comparison(
        s, p, file, expected_sha256=diagnostic(file, s, p, failed=failed)
    )
    assert q.payload["parameter_completeness"] == "complete"
    assert q.payload["source_assignment"] == p.payload["source_assignment"]
    assert q.payload["cross_engine_equivalence"] == (
        "failed" if failed else "not_claimed"
    )
    q.require(ClaimRequirements(parameter_complete=True))
    with pytest.raises(OPLSResolutionError, match="cross_engine_equivalence"):
        q.require(ClaimRequirements(cross_engine_equivalence=True))
    assert not q.payload["production_validated"]


def test_pdms_source_gaps_are_not_numerical_gaps(pdms, tmp_path):
    s, p = pdms
    r = p.payload
    assert (
        r["selected_profile"] == PROFILE and r["source_kind"] == "external_literature"
    )
    assert r["typing_status"] == r["charge_status"] == "complete"
    assert not r["missing_terms"]
    assert r["parameter_completeness"] == "incomplete"
    assert r["missing_conventions"] == ["dielectric", "periodic_image_policy"]
    assert len(r["source_convention_gaps"]) == 9
    assert r["cross_engine_equivalence"] == "not_claimed"
    assert r["charge_record"]["normalization"] is False
    assert r["charge_record"]["redistribution"] is False
    q = export_lammps(s, p, tmp_path / "a")
    q = export_lammps(s, q, tmp_path / "b")
    assert q.payload["lammps_export_status"] == "passed"
    assert q.payload["lammps_runtime_status"] == "not_run"
    assert all(
        (tmp_path / "a" / f).read_bytes() == (tmp_path / "b" / f).read_bytes()
        for f in q.payload["lammps_export"]["files"]
    )
    q, evaluator = build_openmm(s, q)
    try:
        assert q.payload["openmm_build_status"] == "passed"
        assert q.payload["parameter_completeness"] == "incomplete"
        assert q.payload["cross_engine_equivalence"] == "not_claimed"
    finally:
        evaluator.close()


def test_failed_comparison_does_not_block_pdms_preparation(pdms, tmp_path):
    s, p = pdms
    file = tmp_path / "external.json"
    q = record_static_comparison(s, p, file, expected_sha256=diagnostic(file, s, p))
    q = export_lammps(s, q, tmp_path / "export")
    q, model = build_openmm(s, q)
    try:
        assert q.payload["cross_engine_equivalence"] == "failed"
        assert (
            q.payload["lammps_export_status"]
            == q.payload["openmm_build_status"]
            == "passed"
        )
        assert q.payload["missing_conventions"] == p.payload["missing_conventions"]
    finally:
        model.close()


@pytest.mark.parametrize(
    "claim",
    [
        "parameter_complete",
        "lammps_runtime",
        "openmm_runtime",
        "cross_engine_equivalence",
        "scientific_validation",
    ],
)
def test_claims_are_independent(native, claim):
    r = native[1].payload
    result = ClaimRequirements(**{claim: True}).failures(r)
    assert bool(result) == (claim != "parameter_complete")
    if result:
        assert [v["claim"] for v in result] == [claim]


def test_build_and_runtime_are_independent(native, tmp_path):
    s, p = native
    q, model = build_openmm(s, p)
    assert model is not None and q.payload["openmm_build_status"] == "passed"
    q.require(ClaimRequirements(parameter_complete=True))
    assert q.payload["openmm_runtime_status"] == "not_run"
    q = check_openmm(s, q)
    assert q.payload["openmm_runtime_status"] == "passed"
    q.require(ClaimRequirements(openmm_runtime=True, parameter_complete=True))
    q = export_lammps(s, q, tmp_path / "unsupported")
    assert q.payload["lammps_export_status"] == "blocked"
    assert q.payload["parameter_completeness"] == "complete"
    assert not (tmp_path / "unsupported").exists()


@pytest.mark.parametrize("field", CONVENTIONS)
def test_missing_units_scaling_and_conventions_reject(native, field):
    p = native[1].payload
    c = deepcopy(p["conventions"])
    c.pop(field)
    r = assess(
        requirements(native[0], "opls"),
        [r["assignment"] for r in p["interaction_inventory"]],
        c,
        typing_complete=True,
        charge_complete=True,
        source_valid=True,
    )
    assert r["parameter_completeness"] == "incomplete"
    assert field in r["missing_conventions"]


@pytest.mark.parametrize(
    "family",
    [
        "bond-bond",
        "bond-angle",
        "bond-bond_1_3",
        "angle-angle",
        "angle-torsion_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-angle-torsion_1",
    ],
)
def test_missing_class2_couplings_never_zero_filled(native, family):
    s, p = native
    req = requirements(s, "pcff")
    assert any(r["family"] == family for r in req)
    assigned = [
        {**r, "status": "assigned", "value": "bounded audit control"}
        for r in req
        if r["family"] != family
    ]
    result = assess(
        req,
        assigned,
        p.payload["conventions"],
        typing_complete=True,
        charge_complete=True,
        source_valid=True,
    )
    assert result["parameter_completeness"] == "incomplete"
    assert any(
        r["family"] == family and r["assignment"] is None
        for r in result["missing_terms"]
    )


def test_dependency_chain_rejected(native):
    req = [
        {"id": "coupling", "dependencies": ["equilibrium"]},
        {"id": "equilibrium", "dependencies": []},
    ]
    result = assess(
        req,
        [{"id": "coupling", "status": "assigned", "value": 3}],
        native[1].payload["conventions"],
        typing_complete=True,
        charge_complete=True,
        source_valid=True,
    )
    assert {r["id"] for r in result["missing_terms"]} == {"coupling", "equilibrium"}


@pytest.mark.parametrize(
    "field",
    [
        "parameter_completeness",
        "source_identity",
        "charge_record",
        "selected_profile",
        "cross_engine_equivalence",
    ],
)
def test_report_tampering_rejected(native, field):
    s, p = native
    raw = p.payload
    raw[field] = "tampered"
    with pytest.raises(OPLSResolutionError):
        replace(p, json_text=canonical(raw)).validate_integrity(s)


def test_source_replacement_invalidates_complete(native, assets, tmp_path):
    path = tmp_path / "oplsaa.xml"
    shutil.copyfile(XML, path)
    s = from_smiles("CC", random_seed=2026)
    p = prepare_polymer_forcefield(s, assets=replace(assets, native_source_path=path))
    assert p.payload["parameter_complete"]
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(OPLSResolutionError):
        p.validate_integrity(s)


@pytest.mark.parametrize("mutation", ["charge", "profile", "family"])
def test_canonical_constraints_invalidate_completeness(native, assets, mutation):
    s, p = native
    kwargs = (
        {"graph_charges": {**p.payload["charge_record"], "charge_origin": "mutated"}}
        if mutation == "charge"
        else {"profile" if mutation == "profile" else "family": "pcff"}
    )
    q = prepare_polymer_forcefield(
        s, assets=assets, request=ForceFieldRequest(**kwargs)
    )
    assert q.payload["parameter_completeness"] == "incomplete"
    assert q.payload["source_integrity"] == "blocked"


def test_graph_mutation_rejected(native):
    s, p = native
    changed = deepcopy(s)
    changed.topology.sites[next(iter(changed.topology.sites))].metadata[
        "radical_electrons"
    ] = 1
    with pytest.raises(OPLSResolutionError):
        p.validate_integrity(changed)


def test_diagnostic_mutation(native, tmp_path):
    s, p = native
    path = tmp_path / "comparison.json"
    q = record_static_comparison(s, p, path, expected_sha256=diagnostic(path, s, p))
    path.write_text("{}")
    with pytest.raises(OPLSResolutionError):
        q.validate_integrity(s)


def test_relocation_and_fresh_process(native, assets, tmp_path):
    s, p = native
    path = tmp_path / "comparison.json"
    p = record_static_comparison(s, p, path, expected_sha256=diagnostic(path, s, p))
    p = check_openmm(s, p)
    pin = save(s, p, tmp_path / "bundle")
    moved = tmp_path / "moved"
    shutil.copytree(tmp_path / "bundle", moved)
    shutil.copyfile(XML, tmp_path / "source.xml")
    local_assets = replace(
        assets,
        archive_path=None,
        archive_sha256=None,
        native_source_path=tmp_path / "source.xml",
    )
    restored_s, restored_p = reconstruct(
        moved, manifest_sha256=pin, assets=local_assets
    )
    assert restored_p.payload == p.payload
    restored_p.validate_integrity(restored_s)
    code = """
import json,sys
from island.forcefield_policy import OPLSAssets
from island.forcefield_policy.evidence import reconstruct
s,p=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],assets=OPLSAssets(native_source_path=sys.argv[3],native_source_sha256=sys.argv[4]))
p.validate_integrity(s)
print(p.json_text)
"""
    for seed in ["1", "99397"]:
        run = subprocess.run(
            [
                sys.executable,
                "-c",
                code,
                str(moved),
                pin,
                str(tmp_path / "source.xml"),
                local_assets.native_source_sha256,
            ],
            env={**os.environ, "PYTHONHASHSEED": seed},
            check=False,
            capture_output=True,
            text=True,
        )
        assert run.returncode == 0, run.stderr
        assert json.loads(run.stdout) == p.payload


@pytest.mark.parametrize(
    "flag,code",
    [
        ("--require-parameter-complete", 0),
        ("--require-lammps-runtime", 1),
        ("--require-openmm-runtime", 0),
        ("--require-cross-engine-equivalence", 1),
        ("--require-scientific-validation", 1),
    ],
)
def test_cli_flags_independent(assets, tmp_path, flag, code):
    config = tmp_path / "assets.json"
    config.write_text(
        json.dumps(
            {
                "native_source_path": str(XML),
                "native_source_sha256": assets.native_source_sha256,
            }
        )
    )
    report = tmp_path / "report.json"
    cmd = [
        sys.executable,
        "-m",
        "island.inspect_parameters",
        "--smiles",
        "CC",
        "--forcefield",
        "opls",
        "--assets",
        str(config),
        "--assets-sha256",
        sha256(config.read_bytes()).hexdigest(),
        "--json",
        str(report),
        flag,
    ]
    run = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert run.returncode == code, run.stderr
    result = json.loads(report.read_text())
    if flag != "--require-lammps-runtime":
        assert result["parameter_completeness"] == "complete"
        assert result["cross_engine_equivalence"] == "not_claimed"


def test_all_base_files_unchanged():
    from scripts.phase_4j55 import BASE
    from scripts.validate_phase_4j47_receipt import history, pin

    old = history(ROOT, BASE)
    assert all(pin(ROOT / n) == value for n, value in old.items())


def test_complete_external_software_domain_not_gated_by_diagnostic(tmp_path):
    from test_extensible_opls import analytical_xml, manifest_for

    from island.forcefields.opls_external import ExternalOPLSOptions
    from island.forcefields.selection.registry import USER_OPLS

    # Existing labelled analytical fixture, not a new physical parameter profile.
    xml = analytical_xml.__wrapped__()
    source, _ = manifest_for(tmp_path / "source", xml)
    s = from_smiles("CC", random_seed=2026)
    p = prepare_polymer_forcefield(
        s,
        forcefield="opls",
        source_profile=USER_OPLS,
        charge_method="provided",
        options=ExternalOPLSOptions(source, sha256(source.read_bytes()).hexdigest()),
    )
    assert p.payload["parameter_completeness"] == "complete", p.payload["diagnostics"]
    evidence = tmp_path / "comparison.json"
    p = record_static_comparison(
        s, p, evidence, expected_sha256=diagnostic(evidence, s, p)
    )
    assert (
        p.payload["parameter_complete"]
        and p.payload["cross_engine_equivalence"] == "failed"
    )
    p.require(ClaimRequirements(parameter_complete=True))


def test_pcff_real_missing_terms(assets):
    from test_phase_4j46 import native_options

    s = from_smiles("CC", random_seed=2026)
    p = prepare_polymer_forcefield(
        s, forcefield="pcff", charge_method="native", options=native_options("pcff")
    )
    assert p.payload["parameter_completeness"] == "incomplete"
    assert p.payload["missing_terms"]
    assert any(
        "torsion" in r["family"] or "bond" in r["family"]
        for r in p.payload["missing_terms"]
    )
    assert p.payload["cross_engine_equivalence"] == "not_claimed"


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
@pytest.mark.parametrize("method", ["am1bcc", "provided", "resp_psi4"])
def test_amber_charge_provenance_retained_with_missing_asset(family, method):
    from island.forcefields.ambertools import AmberToolsOptions

    s = from_smiles("CC", random_seed=2026)
    options = AmberToolsOptions(
        force_field=family,
        charge_method="provided" if method != "am1bcc" else "am1bcc",
        charge_source="external declared charge" if method != "am1bcc" else None,
        provided_charges={i: 0.0 for i in s.topology.sites}
        if method != "am1bcc"
        else None,
    )
    p = prepare_polymer_forcefield(
        s, forcefield=family, charge_method=method, options=options
    )
    assert p.payload["parameter_completeness"] == "incomplete"
    assert p.payload["source_assignment"]["charge_method"] == method
    assert p.payload["cross_engine_equivalence"] == "not_claimed"
    assert not p.payload["production_validated"]


@pytest.mark.parametrize("kind", ["source_archive", "paper", "charge_vector"])
def test_pdms_source_and_charge_mutation(assets, sources, tmp_path, kind):
    s = fixture(0)
    if kind == "source_archive":
        path = tmp_path / "archive.zip"
        shutil.copyfile(ARCHIVE, path)
        p = prepare_polymer_forcefield(
            s, assets=replace(assets, archive_path=path), sources=sources
        )
        path.write_bytes(path.read_bytes() + b"mutated")
        with pytest.raises(OPLSResolutionError):
            p.validate_integrity(s)
    elif kind == "paper":
        path = tmp_path / "paper.pdf"
        shutil.copyfile(PAPER, path)
        p = prepare_polymer_forcefield(
            s, assets=assets, sources=replace(sources, paper=path)
        )
        path.write_bytes(path.read_bytes() + b"mutated")
        with pytest.raises(OPLSResolutionError):
            p.validate_integrity(s)
    else:
        p = prepare_polymer_forcefield(s, assets=assets, sources=sources)
        charges = deepcopy(p.payload["charge_record"])
        charges["vector"][0][1] = "0"
        q = prepare_polymer_forcefield(
            s,
            assets=assets,
            sources=sources,
            request=ForceFieldRequest(graph_charges=charges),
        )
        assert q.payload["parameter_completeness"] == "incomplete"


def test_public_psmiles_path_only_family(assets, sources):
    from island.core import SimulationBox
    from island.forcefield_policy import prepare_polymer

    p = prepare_polymer(
        psmiles="[*:1]O[Si](C)(C)[*:2]",
        forcefield="opls",
        dp=2,
        head_end_group="[*][Si](C)(C)C",
        tail_end_group="[*]C",
        box=SimulationBox(40, 40, 40),
        assets=assets,
        sources=sources,
    )
    assert p.payload["selected_profile"] == PROFILE
    assert p.payload["charge_record"]["charge_generation_method"] == "author_dft_npa"
    p.validate_integrity()


@pytest.mark.parametrize(
    "field",
    ["charge_record", "source_identity", "selected_profile", "parameter_complete"],
)
def test_resigned_mutations_independently_reconstructed(native, field):
    from island.forcefields.selection.inventory import digest

    s, p = native
    raw = p.payload
    raw.pop("identity")
    raw[field] = "tampered"
    text = canonical({**raw, "identity": digest(raw)})
    with pytest.raises(OPLSResolutionError):
        replace(p, json_text=text, _original=text).validate_integrity(s)
