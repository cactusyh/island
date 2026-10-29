"""Unified selection, real bounded source assignments and data-only completion."""

import json
import os
import shutil
import subprocess
import sys
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from island import SimulationBox
from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields.preparation import OPLSOptions, PCFFOptions
from island.forcefields.selection import (
    CompletionRecord,
    SourceSelection,
    UnifiedForceFieldError,
    capability_registry,
    create_evaluator,
    derive_geometric_pair,
    fit_harmonic_bond,
    prepare_polymer_forcefield,
)
from island.forcefields.selection.inventory import canonical, digest
from island.forcefields.selection.registry import OPLS, PCFF, PDMS, USER_OPLS

ROOT = Path(__file__).resolve().parents[1]
PCFF_FILE = Path(
    "/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/tools/msi2lmp/frc_files/pcff.frc"
)
OPLS_FILE = (
    ROOT.parent
    / "island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml"
)
ZIP = ROOT / "island_ffref/pcff/PDMS_benchmark-main.zip"
ZIP_HASH = "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4"


def pin(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def native_options(family):
    path = PCFF_FILE if family == "pcff" else OPLS_FILE
    if not path.exists():
        pytest.skip("Explicit authorized local source absent")
    return (
        PCFFOptions(path, lj=(0, 0, 1), coulomb=(0, 0, 1))
        if family == "pcff"
        else OPLSOptions(path)
    )


def prepare(system, family="opls", **kwargs):
    return prepare_polymer_forcefield(
        system, family, charge_method="native", options=native_options(family), **kwargs
    )


def test_registry_owned_and_stage_separation():
    registry = capability_registry()
    assert len({r["profile_identity"] for r in registry}) == len(registry)
    assert {r["family"] for r in registry} == {"opls", "pcff", "gaff", "gaff2"}
    for r in registry:
        assert r["stages"]["scientifically_validated"] is False
        assert not r["production_validated"]
        assert r["simulation_readiness"] == "not_established"
        assert all(
            field in r
            for field in (
                "supported_elements",
                "supported_atom_environments",
                "bond_terms",
                "angle_terms",
                "torsions",
                "impropers",
                "class2_cross_terms",
                "nonbonded_terms",
                "charge_methods",
                "scaling_14",
                "mixing_rule",
                "exclusions",
                "long_range_electrostatics",
                "backend_support",
                "periodic_backend_support",
            )
        )
    assert (
        next(r for r in registry if r["profile_identity"] == PDMS)["source_kind"]
        == "external_literature"
    )
    registry.clear()
    assert capability_registry()


@pytest.mark.parametrize(
    "family,smiles",
    [("pcff", "C"), ("opls", "CC"), ("opls", "c1ccccc1"), ("opls", "COCCOC")],
)
def test_real_assignment_no_source_fallback(family, smiles):
    s = from_smiles(smiles, random_seed=46)
    r = prepare(s, family)
    p = r.payload
    assert p["status"] == "prepared", p["diagnostics"]
    assert all(
        p["stages"][k]
        for k in (
            "typing_complete",
            "charge_complete",
            "parameter_complete",
            "executable",
        )
    )
    assert (
        not p["stages"]["numerically_verified"]
        and not p["stages"]["scientifically_validated"]
    )
    assert (
        p["typed_graph_identity"] and p["charge_identity"] and p["parameter_identity"]
    )
    assert not p["blocked_terms"]
    assert p["assigned_terms"] == sorted(p["assigned_terms"])
    assert r.require_executable(s).identity
    s.coordinates.set(next(iter(s.topology.sites)), [10, 20, 30])
    s.metadata["coordinate_history"] = ["changed"]
    r.validate_integrity(s)


def test_pcff_strict_missing_coupling_and_dependency_chain():
    s = from_smiles("CC", random_seed=46)
    r = prepare(s, "pcff")
    p = r.payload
    assert p["status"] == "parameter_incomplete"
    missing = [a for a in p["requirements"] if a["id"] in p["blocked_terms"]]
    assert missing and all(a["family"] == "bond-bond_1_3" for a in missing)
    assert all(a["dependencies"] and a.get("supplied_types") for a in missing)
    assert not p["stages"]["executable"]
    with pytest.raises(UnifiedForceFieldError) as error:
        create_evaluator(s, r)
    assert error.value.result is r and error.value.status == "parameter_incomplete"


@pytest.mark.parametrize("smiles", ["C[Si](C)(C)O[Si](C)(C)C", "C[SiH3]", "OP(=O)(O)O"])
@pytest.mark.parametrize("family", ["pcff", "opls"])
def test_native_si_p_never_fallback(smiles, family):
    p = prepare(from_smiles(smiles, random_seed=46), family).payload
    assert p["status"] == "typing_ambiguous"
    assert p["blocked_terms"] and not p["stages"]["executable"]
    assert p["completion_records"] == []


@pytest.mark.parametrize(
    "psmiles", ["[*:1]CC[*:2]", "[*:1]COC[*:2]", "[*:1]c1ccc([*:2])cc1"]
)
def test_psmiles_linear_chain(psmiles):
    s = build_linear_polymer(
        psmiles,
        dp=2,
        head_end_group="[*]C",
        tail_end_group="[*]C",
        add_hydrogens=True,
        random_seed=46,
    )
    r = prepare(s)
    assert r.payload["status"] in (
        "prepared",
        "parameter_incomplete",
        "charge_incomplete",
    ), r.payload["diagnostics"]
    assert r.payload["stages"]["typing_complete"]
    assert r.payload["final_graph"]["system_metadata"]["polymer"]


def test_multichain_and_crosslinked_final_graphs():
    # Two components and a chemically valid bridge/ring network are inspected
    # as final graphs; no reuse of pre-crosslink repeat-unit atom labels.
    for smiles, components in [("CC.CCO", 2), ("C1CC2CCC1C2", 1)]:
        s = from_smiles(smiles, random_seed=46)
        p = prepare(s).payload
        assert len(p["final_graph"]["components"]) == components
        assert p["status"] in ("prepared", "parameter_incomplete")
        assert p["stages"]["typing_complete"]
        assert (
            len([r for r in p["requirements"] if r["family"] == "bond"])
            == s.number_of_bonds
        )


@pytest.mark.parametrize("family", ["pcff", "opls", "gaff", "gaff2"])
def test_unavailable_explicit_assets(family, monkeypatch):
    def forbidden(*a, **kw):
        raise AssertionError("Engine constructed for blocked request")

    monkeypatch.setattr("island.forcefields.research.create_evaluator", forbidden)
    s = from_smiles("CC", random_seed=46)
    r = prepare_polymer_forcefield(
        s, family, charge_method="am1bcc" if family.startswith("gaff") else "native"
    )
    assert r.payload["status"] == "asset_unavailable"
    with pytest.raises(UnifiedForceFieldError) as error:
        create_evaluator(s, r)
    assert error.value.result is r


@pytest.mark.parametrize(
    "defect", ["family", "policy", "profile", "charge", "source_profile"]
)
def test_request_mismatch_classified(defect):
    kwargs = {
        "forcefield": "opls",
        "charge_method": "native",
        "options": native_options("opls"),
    }
    if defect == "family":
        kwargs["forcefield"] = "unknown"
    elif defect == "policy":
        kwargs["resolution_policy"] = "zero_fill"
    elif defect == "profile":
        kwargs["source_profile"] = PCFF
    elif defect == "charge":
        kwargs["charge_method"] = "resp_psi4"
    else:
        kwargs["source_profile"] = "absent"
    r = prepare_polymer_forcefield(from_smiles("CC", random_seed=46), **kwargs)
    assert r.payload["status"] in ("source_incompatible", "charge_incomplete")
    assert not r.payload["stages"]["executable"]


@pytest.mark.parametrize("defect", ["bond", "charge", "order", "component", "mass"])
def test_stale_graph_rejected(defect):
    s = from_smiles("CC", random_seed=46)
    r = prepare(s)
    i = next(iter(s.topology.sites))
    if defect == "bond":
        s.topology.bonds.pop(next(iter(s.topology.bonds)))
    elif defect == "charge":
        s.topology.sites[i].formal_charge = 1
    elif defect == "order":
        s.topology.sites = dict(reversed(list(s.topology.sites.items())))
    elif defect == "component":
        s.topology.sites[i].metadata["chain_id"] = 99
    else:
        s.topology.sites[i].mass += 1
    with pytest.raises(UnifiedForceFieldError):
        r.require_executable(s)


@pytest.mark.parametrize(
    "mutation", ["status", "charges", "requirements", "source", "profile"]
)
def test_resigned_report_mutation_rejected(mutation):
    s = from_smiles("CC", random_seed=46)
    r = prepare(s)
    p = r.payload
    if mutation == "status":
        p["status"] = "asset_unavailable"
    elif mutation == "charges":
        p["charge_identity"] = "0" * 64
    elif mutation == "requirements":
        p["requirements"].clear()
    elif mutation == "source":
        p["source_hashes"] = ["0" * 64]
    else:
        p["profile_identity"] = PDMS
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    bad = replace(r, json_text=canonical(p))
    with pytest.raises(UnifiedForceFieldError):
        bad.require_executable(s)


def test_source_mutation_and_relocation(tmp_path):
    path = tmp_path / "opls.xml"
    shutil.copyfile(native_options("opls").source_path, path)
    s = from_smiles("CC", random_seed=46)
    r = prepare_polymer_forcefield(
        s, "opls", charge_method="native", options=OPLSOptions(path)
    )
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(UnifiedForceFieldError, match="source_mutation"):
        r.validate_integrity(s)


@pytest.mark.parametrize(
    "policy", ["derived_from_declared_formula", "fit_from_external_reference"]
)
def test_completion_never_guesses_missing_class2(policy):
    r = prepare(from_smiles("CC", random_seed=46), "pcff", resolution_policy=policy)
    assert r.payload["status"] == "parameter_incomplete"
    assert r.payload["diagnostics"][0]["stage"] == "completion_evidence_unavailable"
    assert r.payload["blocked_terms"]


def lj_inputs(root):
    raw = "a 3.0 0.2\nb 4.0 0.8\n"
    (root / "rows.txt").write_text(raw)
    return [
        {
            "file": "rows.txt",
            "sha256": pin(root / "rows.txt"),
            "line": n,
            "raw": line,
            "units": ["angstrom", "kcal/mol"],
            "columns": [1, 2],
            "values": line.split()[1:],
            "label": line.split()[0],
            "citation": "Analytical software fixture only",
            "source_identity": OPLS,
        }
        for n, line in enumerate(raw.splitlines(), 1)
    ]


def test_derived_evidence_independent_reconstruction_and_policy(tmp_path):
    a, b = lj_inputs(tmp_path)
    r = derive_geometric_pair(tmp_path, OPLS, a, b)
    p = r.validate_integrity(tmp_path)
    assert p["kind"] == "derived" and not p["execution_authority"]
    assert float(p["values"][0]) == pytest.approx(np.sqrt(12))
    assert float(p["values"][1]) == pytest.approx(0.4)
    s = from_smiles("CC", random_seed=46)
    selection = SourceSelection(OPLS, native_options("opls"))
    out = prepare_polymer_forcefield(
        s,
        "opls",
        source_profile=selection,
        charge_method="native",
        resolution_policy="derived_from_declared_formula",
        completion_records=[r],
        completion_root=tmp_path,
    )
    assert (
        out.payload["status"] == "prepared"
        and out.payload["completion_records"][0]["kind"] == "derived"
    )
    blocked = prepare(s, completion_records=[r], completion_root=tmp_path)
    assert blocked.payload["status"] == "source_incompatible"
    with pytest.raises(ValueError, match="No declared geometric"):
        derive_geometric_pair(tmp_path, PCFF, a, b)


@pytest.mark.parametrize(
    "defect", ["hash", "raw", "units", "value", "source", "path", "label"]
)
def test_derived_bad_inputs(tmp_path, defect):
    a, b = lj_inputs(tmp_path)
    if defect == "hash":
        a["sha256"] = "0" * 64
    elif defect == "raw":
        a["raw"] = "a 9 9"
    elif defect == "units":
        a["units"] = ["nm", "kJ/mol"]
    elif defect == "value":
        a["values"] = ["9", "9"]
    elif defect == "source":
        a["source_identity"] = PCFF
    elif defect == "label":
        a["label"] = "*"
    else:
        a["file"] = "../rows.txt"
    with pytest.raises(ValueError):
        derive_geometric_pair(tmp_path, OPLS, a, b)


def scan(root):
    r = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    data = {
        "schema": "island_isolated_bond_scan_v1",
        "source_identity": OPLS,
        "reference_kind": "reference_engine",
        "reference_provenance": "Analytical test only, NOT physical force-field authority",
        "citation": "J46 software test",
        "functional_form": "K*(r-r0)^2+offset",
        "units": {"r": "angstrom", "energy": "kcal/mol"},
        "isolated_bond_contribution": True,
        "labels": ["a", "b"],
        "r": r,
        "energy": [100 * (x - 1.25) ** 2 + 2 for x in r],
        "weights": [1] * len(r),
    }
    (root / "scan.json").write_text(json.dumps(data))
    return data


def test_fit_external_provenance_and_reconstruction(tmp_path):
    scan(tmp_path)
    r = fit_harmonic_bond(
        tmp_path,
        OPLS,
        reference_file="scan.json",
        expected_sha256=pin(tmp_path / "scan.json"),
    )
    p = r.validate_integrity(tmp_path)
    assert p["kind"] == "fitted_external" and not p["execution_authority"]
    assert float(p["values"]["K"]) == pytest.approx(100)
    assert float(p["values"]["r0"]) == pytest.approx(1.25)
    assert p["uncertainty"]["degrees_of_freedom"] == 3
    out = prepare(
        from_smiles("CC", random_seed=46),
        resolution_policy="fit_from_external_reference",
        completion_records=[r],
        completion_root=tmp_path,
    )
    assert out.payload["status"] == "prepared"
    p["values"]["K"] = "0"
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(ValueError, match="Contradictory"):
        CompletionRecord(canonical(p)).validate_integrity(tmp_path)


@pytest.mark.parametrize(
    "defect", ["units", "reference", "isolation", "rank", "weights", "family"]
)
def test_fit_bad_reference(tmp_path, defect):
    data = scan(tmp_path)
    if defect == "units":
        data["units"]["r"] = "unknown"
    elif defect == "reference":
        data["reference_provenance"] = ""
    elif defect == "isolation":
        data["isolated_bond_contribution"] = False
    elif defect == "rank":
        data["r"] = [1.0] * 6
    elif defect == "weights":
        data["weights"][0] = -1
    else:
        data["source_identity"] = PCFF
    (tmp_path / "scan.json").write_text(json.dumps(data))
    with pytest.raises(ValueError):
        fit_harmonic_bond(
            tmp_path,
            OPLS,
            reference_file="scan.json",
            expected_sha256=pin(tmp_path / "scan.json"),
        )


@pytest.mark.parametrize("dp", [0, 2, 4])
def test_external_pdms_existing_research_adapter(dp):
    if not ZIP.exists():
        pytest.skip("User author archive absent")
    from island.forcefields.pdms_opls import PDMSOPLSOptions

    if dp:
        s = build_linear_polymer(
            "[*:1]O[Si](C)(C)[*:2]",
            dp=dp,
            head_end_group="[*][Si](C)(C)C",
            tail_end_group="[*]C",
            add_hydrogens=True,
            random_seed=46,
        )
        anchors = (s.metadata["polymer"]["head_site_id"],)
    else:
        s = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=46)
        anchors = ()
    s.box = SimulationBox(40, 40, 40)
    r = prepare_polymer_forcefield(
        s,
        "opls",
        source_profile=PDMS,
        charge_method="author_dft_npa",
        options=PDMSOPLSOptions(ZIP, ZIP_HASH, e_oxygen_site_ids=anchors),
    )
    p = r.payload
    assert p["status"] == "prepared", p["diagnostics"]
    assert p["source_kind"] == "external_literature" and not p["blocked_terms"]
    assert r.require_executable(s).payload["model"]["paper_Qeq_model"] is False
    assert not p["stages"]["scientifically_validated"]


def test_user_manifest_adapter_projection(tmp_path):
    # Reuse existing explicitly analytical test XML, never literature data.
    from test_extensible_opls import analytical_xml, manifest_for

    from island.forcefields.opls_external import ExternalOPLSOptions

    xml = analytical_xml.__wrapped__()
    path, _document = manifest_for(tmp_path / "source", xml)
    r = prepare_polymer_forcefield(
        from_smiles("CC", random_seed=46),
        "opls",
        source_profile=USER_OPLS,
        charge_method="provided",
        options=ExternalOPLSOptions(path, pin(path)),
    )
    assert r.payload["status"] == "prepared", r.payload["diagnostics"]
    assert r.payload["source_kind"] == "user_provided"
    assert (
        r.payload["external_assignment"]["profile_identity"]
        == "island_external_opls_analytical_test_v1"
    )
    (path.parent / "LICENSE").write_text("mutated")
    with pytest.raises(UnifiedForceFieldError):
        r.validate_integrity()


@pytest.mark.parametrize("family", ["gaff", "gaff2"])
@pytest.mark.parametrize(
    "method", ["am1bcc", "provided", "external_provided", "resp_psi4"]
)
def test_amber_charge_origin_boundary(family, method):
    from island.forcefields.ambertools import AmberToolsOptions

    s = from_smiles("CC", random_seed=46)
    kwargs = {
        "force_field": family,
        "charge_method": "am1bcc" if method == "am1bcc" else "provided",
    }
    if method != "am1bcc":
        kwargs.update(
            provided_charges={i: 0.0 for i in s.topology.sites},
            charge_source="explicit external software test vector",
        )
    r = prepare_polymer_forcefield(
        s, family, charge_method=method, options=AmberToolsOptions(**kwargs)
    )
    assert r.payload["status"] in ("asset_unavailable", "charge_incomplete")
    if method == "resp_psi4":
        assert (
            r.payload["diagnostics"][0]["stage"] == "independent_resp_evidence_required"
        )
    else:
        assert r.payload["charge_origin"]["method"] in ("AM1-BCC", method)


def test_fresh_process_relocation_and_cli(tmp_path):
    source = tmp_path / "first.xml"
    shutil.copyfile(native_options("opls").source_path, source)
    env = {
        **os.environ,
        "PYTHONPATH": str(ROOT / "src"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }

    def run(source, output, smiles="CC"):
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "island.inspect_forcefield",
                "--smiles",
                smiles,
                "--forcefield",
                "opls",
                "--charge-method",
                "native",
                "--source",
                str(source),
                "--json",
                str(output),
            ],
            cwd=tmp_path,
            env=env,
            check=False,
            capture_output=True,
            text=True,
        )

    first = tmp_path / "first.json"
    assert run(source, first).returncode == 0
    other = tmp_path / "relocated"
    other.mkdir()
    moved = other / "second.xml"
    shutil.move(source, moved)
    second = other / "second.json"
    assert run(moved, second).returncode == 0
    assert json.loads(first.read_text()) == json.loads(second.read_text())
    blocked = other / "blocked.json"
    assert run(moved, blocked, "C[SiH3]").returncode == 1
    assert json.loads(blocked.read_text())["status"] == "typing_ambiguous"


def test_invalid_graph_is_classified_and_no_engine():
    s = from_smiles("CC", random_seed=46)
    s.topology.sites[1].atomic_number = 14
    r = prepare(s)
    assert r.payload["status"] == "typing_ambiguous"
    with pytest.raises(UnifiedForceFieldError):
        create_evaluator(s, r)


@pytest.mark.parametrize("family,smiles", [("pcff", "C"), ("opls", "CC")])
def test_unified_evaluator_delegates_exact_native_model(family, smiles):
    from island.forcefields.preparation import create_evaluator as existing

    s = from_smiles(smiles, random_seed=46)
    r = prepare(s, family)
    a = create_evaluator(s, r).evaluate_fresh()
    b = existing(s, r.require_executable(s)).evaluate_fresh()
    assert a == b
    assert np.isfinite(a.potential_energy)
    # Software equality does not qualify periodic science or a whole profile.
    assert not r.payload["stages"]["scientifically_validated"]


def test_charge_failure_classified_from_backend_cause():
    from island.exceptions import ChargeAssignmentError
    from island.forcefields.selection.api import _failure_status

    outer = ValueError("converter/import boundary")
    outer.__cause__ = ChargeAssignmentError("source charge vector incomplete")
    assert _failure_status(outer, "source") == "charge_incomplete"


def test_source_checksum_is_roundtrip_stable_for_double_digit_ids():
    assert digest({2: "a", 11: "b"}) == digest({"2": "a", "11": "b"})
    with pytest.raises(ValueError):
        canonical({"invalid": float("nan")})


def test_blocked_invalid_graph_preserves_classification():
    s = from_smiles("CC", random_seed=46)
    s.topology.sites[1].atomic_number = 14
    r = prepare(s)
    with pytest.raises(UnifiedForceFieldError) as error:
        create_evaluator(s, r)
    assert error.value.status == "typing_ambiguous"
    assert error.value.result is r


def test_receipt_identity_parser_rejects_tamper(tmp_path):
    from scripts.validate_phase_4j46_receipt import read

    path = tmp_path / "record.json"
    value = {"schema": "analytical_test"}
    path.write_text(json.dumps({**value, "identity": digest(value)}))
    assert read(path)["schema"] == "analytical_test"
    path.write_text(path.read_text().replace("analytical_test", "tampered_test"))
    with pytest.raises(ValueError, match="Changed identity"):
        read(path)
