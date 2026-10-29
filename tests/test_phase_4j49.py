"""Raw electrostatic audits cannot silently promote the source profile to MD."""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from island.electrostatic_audit import (
    COMPONENTS,
    InputBinding,
    bind_inputs,
    finite_difference,
)
from island.forcefields.opls_registry import (
    ForceFieldRequest,
    OPLSAssets,
    OPLSResolutionError,
    Preparation,
    prepare_forcefield,
    prepare_polymer_forcefield,
    use_opls_assets,
)
from island.forcefields.opls_registry.registry import ARCHIVE_SHA256
from island.forcefields.pdms_opls import PROFILE, PDMSOPLSError
from island.forcefields.selection.inventory import canonical, digest
from island.opls_lammps import run_md
from scripts.phase_4j49 import (
    MODES,
    STEPS,
    ChargedReference,
    component_forces,
    packed_components,
)
from tests.test_phase_4j48 import ARCHIVE, fixture


@pytest.fixture(scope="module")
def assets():
    if not ARCHIVE.is_file():
        pytest.skip("Local source ZIP unavailable")
    return OPLSAssets(ARCHIVE, ARCHIVE_SHA256)


def polynomial(x):
    bonded = float(np.sum(x * x))
    lj = float(np.sum(0.3 * x**4))
    coulomb = float(np.sum(2 * x))
    return {
        "bonded": bonded,
        "lj": lj,
        "coulomb": coulomb,
        "total": bonded + lj + coulomb,
    }


def polynomial_forces(x):
    bonded, lj, coulomb = -2 * x, -1.2 * x**3, np.full_like(x, -2)
    return {
        "bonded": bonded,
        "lj": lj,
        "coulomb": coulomb,
        "total": bonded + lj + coulomb,
    }


def test_independent_fd_preserves_every_raw_value():
    xyz = np.array([[0.5, 0.7, 0.9], [1.1, 1.3, 1.5]])
    original = xyz.copy()
    record = finite_difference(polynomial, xyz, polynomial_forces(xyz))
    assert np.array_equal(xyz, original)
    assert record["passed"] and record["bound"] == 1e-4
    assert record["comparisons_per_component"] == 18
    assert record["steps_angstrom"] == list(STEPS)
    for row in record["raw_records"]:
        n, a, h = row["atom_index"], row["axis"], row["h_angstrom"]
        plus, minus = xyz.copy(), xyz.copy()
        plus[n, a] += h
        minus[n, a] -= h
        for component, values in row["components"].items():
            assert values["positive_energy_kj_mol"] == polynomial(plus)[component]
            assert values["negative_energy_kj_mol"] == polynomial(minus)[component]
            assert values["central_difference_force"] == -(
                polynomial(plus)[component] - polynomial(minus)[component]
            ) / (2 * h)


@pytest.mark.parametrize("scale", [-1, 10, 0.1, 4.184, 1 / 4.184])
def test_detect_sign_and_unit_errors_without_correction(scale):
    xyz = np.array([[0.7, 0.9, 1.1]])
    forces = {k: f * scale for k, f in polynomial_forces(xyz).items()}
    r = finite_difference(polynomial, xyz, forces)
    assert not r["passed"]
    assert not any(r["passed_by_component"].values())
    for row in r["raw_records"]:
        coul = row["components"]["coulomb"]
        assert coul["analytical_force"] == -2 * scale
        if scale == -1:
            assert coul["sign_reversal_error"] < 1e-9
        elif scale == 10:
            assert coul["nm_to_angstrom_error"] < 1e-9
        elif scale == 1 / 4.184:
            assert coul["kcal_to_kj_error"] < 1e-9


@pytest.mark.parametrize(
    "bad", ["bound", "steps", "shape", "nan", "forces", "energies", "empty"]
)
def test_fd_malformed_input_classified(bad):
    xyz = np.ones((1, 3))
    kwargs = {}
    force = polynomial_forces(xyz)
    energy = polynomial
    if bad == "bound":
        kwargs["bound"] = 1
    elif bad == "steps":
        kwargs["steps"] = (1e-4,)
    elif bad == "shape":
        xyz = np.ones(3)
    elif bad == "nan":
        xyz[0, 0] = np.nan
    elif bad == "forces":
        force.pop("coulomb")
    elif bad == "energies":
        energy = lambda x: {"total": 1}
    elif bad == "empty":
        kwargs["probes"] = []
    with pytest.raises(OPLSResolutionError):
        finite_difference(energy, xyz, force, **kwargs)


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_canonical_binding_preserves_provenance_and_source_numbers(assets, dp):
    system = fixture(dp)
    with use_opls_assets(assets):
        prepared = prepare_polymer_forcefield(system, forcefield="opls")
    binding = bind_inputs(system, prepared)
    assert binding.validate(system, prepared) == binding.payload
    p = binding.payload
    assert p["family"] == "opls" and p["selected_profile"] == PROFILE
    assert p["charge_origin"] == prepared.payload["charge_origin"]
    assert p["domain"] == prepared.payload["domain"]
    assert p["numerical_contract"]["LAMMPS"]["tolerance"] == 1e-4
    assert p["numerical_contract"]["OpenMM"]["tolerance"] == 1e-7
    assert p["source_identity"] and p["graph_identity"] and p["atoms"]
    assert binding.identity == digest(binding.payload)


@pytest.mark.parametrize(
    "mutation",
    [
        "charge",
        "order",
        "exclusions",
        "coul14",
        "lj14",
        "mixing",
        "box",
        "cutoff",
        "reciprocal",
        "source",
        "profile",
        "units",
    ],
)
def test_resigned_binding_mutation_rejected(assets, mutation):
    system = fixture()
    prepared = prepare_forcefield(system, assets=assets)
    b = bind_inputs(system, prepared)
    p = b.payload
    if mutation == "charge":
        p["atoms"][0]["charge_e"] = "0.001"
    elif mutation == "order":
        p["atom_order"].reverse()
    elif mutation == "exclusions":
        p["special_pairs"].pop()
    elif mutation == "box":
        p["box_angstrom"][0] = 41
    elif mutation == "source":
        p["source_identity"] = "0" * 64
    elif mutation == "profile":
        p["selected_profile"] = "pcff"
    elif mutation == "reciprocal":
        p["numerical_contract"]["LAMMPS"]["tolerance"] = 1e-10
    else:
        p["numerical_contract"][
            {
                "coul14": "coulomb_14",
                "lj14": "lj_14",
                "mixing": "lj_mixing",
                "cutoff": "coulomb_cutoff_angstrom",
                "units": "units",
            }[mutation]
        ] = "changed"
    with pytest.raises(
        OPLSResolutionError, match="electrostatic_input_binding_mutation"
    ):
        InputBinding(canonical(p)).validate(system, prepared)


@pytest.mark.parametrize("mutation", ["charge", "source", "profile", "family"])
def test_resigned_canonical_assignment_rejected_before_audit(assets, mutation):
    system = fixture()
    prepared = prepare_forcefield(system, assets=assets)
    p = prepared.payload
    if mutation == "charge":
        p["graph_charges"]["vector"][0][1] = "0.1"
    else:
        p[
            {
                "source": "source_identity",
                "profile": "selected_profile",
                "family": "family",
            }[mutation]
        ] = "pcff"
    p["identity"] = digest({k: v for k, v in p.items() if k != "identity"})
    forged = Preparation(canonical(p), prepared.request, assets)
    with pytest.raises(OPLSResolutionError):
        bind_inputs(system, forged)


def test_incomplete_archive_and_wrong_source_hash(assets, tmp_path):
    file = tmp_path / "incomplete.zip"
    file.write_bytes(ARCHIVE.read_bytes()[:1024])
    for options in (
        replace(assets, archive_path=file),
        replace(assets, archive_sha256="0" * 64),
    ):
        p = prepare_forcefield(fixture(), assets=options)
        with pytest.raises(OPLSResolutionError):
            bind_inputs(fixture(), p)


@pytest.mark.parametrize(
    "smiles", ["[SiH3]O[SiH3]", "CO[Si](C)(C)C", "C[Si](O)(C)O[Si](C)(C)C"]
)
def test_unsupported_terminals_fail_closed(assets, smiles):
    from island.chemistry import from_smiles

    system = from_smiles(smiles, random_seed=2026)
    p = prepare_forcefield(system, assets=assets)
    with pytest.raises(OPLSResolutionError):
        bind_inputs(system, p)


def test_wrong_family_rejected(assets):
    s = fixture()
    p = prepare_forcefield(s, ForceFieldRequest("pcff", profile=PROFILE), assets=assets)
    with pytest.raises(OPLSResolutionError):
        bind_inputs(s, p)


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_no_engine_or_checkpoint_before_source_gates(assets, tmp_path, dp):
    s = fixture(dp)
    p = prepare_forcefield(s, assets=assets)
    with pytest.raises(OPLSResolutionError, match="scientifically_unqualified"):
        run_md(s, p, tmp_path / "blocked", engine="unavailable", steps=10)
    assert not (tmp_path / "blocked").exists()
    assert not p.payload["production_validated"]
    assert p.payload["simulation_readiness"] == "not_established"


@pytest.mark.parametrize("periodic", [False, True])
def test_generic_charged_backend_fd(periodic):
    pytest.importorskip("openmm")
    xyz = np.array(
        [
            [0.41 if periodic else 14.0, 19.0, 20.0],
            [38.17 if periodic else 17.17, 20.23, 21.09],
        ]
    )
    reference = ChargedReference(xyz, (1.0, -1.0))
    try:
        force = component_forces(reference)

        def evaluate(x):
            reference.set(x)
            return packed_components(reference)

        fd = finite_difference(evaluate, xyz, force)
        assert fd["passed"]
        assert fd["max_error_by_component"]["coulomb"] < 1e-4
    finally:
        reference.close()


def test_diagnostic_controls_do_not_mutate_source_contract():
    assert MODES["author"] == ""
    assert "1e-10" in MODES["converged-pppm"]
    assert "ewald" in MODES["matched-ewald"]
    from island.forcefields.opls_registry import canonical_registry

    p = next(
        e for e in canonical_registry()["entries"] if e["profile_identity"] == PROFILE
    )
    assert not p["backend_support"]["MD_LAMMPS"]
    assert "PPPM 1e-4" in p["long_range_convention"]["LAMMPS"]


def test_kernel_patch_is_explicit_and_source_bound(tmp_path):
    from scripts.build_phase_4j49_kernel import patch
    from tests.test_phase_4j48 import ROOT

    source = ROOT.parent / "lammps/lammps/src"
    if not source.exists():
        pytest.skip("Local LAMMPS source unavailable")
    import shutil

    for name in ("pair_lj_charmm_coul_long.cpp", "pair_coul_long.cpp"):
        target = tmp_path / "src/KSPACE" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / "KSPACE" / name, target)
    r = patch(tmp_path)
    assert [v["occurrences"] for v in r.values()] == [3, 2]
    assert all(v["original_sha256"] != v["modified_sha256"] for v in r.values())
    with pytest.raises(ValueError):
        patch(tmp_path)


def test_committed_raw_evidence_if_available():
    root = Path(__file__).resolve().parents[1] / "docs/evidence/phase_4j49/audit"
    if not (root / "summary.json").exists():
        pytest.skip("J49 receipt not generated yet")
    p = json.loads((root / "summary.json").read_text())
    assert p["end_to_end_acceptance_gate"] == 1 and not p["author_profile_MD_executed"]
    for dp in (0, 2, 4, 6):
        modes = json.loads((root / f"components-dp{dp}.json").read_text())
        assert not modes["author"]["finite_difference"]["passed"]
        assert modes["exact-erfc-ewald"]["finite_difference"]["passed"]
        assert all(
            v <= 3e-4 for v in modes["exact-erfc-ewald"]["energy_error_kj_mol"].values()
        )


@pytest.mark.parametrize("blocked", ["energy", "force", "FD"])
def test_diagnostic_workflow_rejects_failed_static_before_engine(
    assets, tmp_path, blocked
):
    from scripts.phase_4j49_workflow import run_control

    s = fixture()
    with use_opls_assets(assets):
        p = prepare_polymer_forcefield(s, forcefield="opls")
    record = {
        "energy_error_kj_mol": {k: 0.0 for k in COMPONENTS},
        "force_error_kj_mol_angstrom": {k: 0.0 for k in COMPONENTS},
        "finite_difference": {"passed": True},
    }
    if blocked == "FD":
        record["finite_difference"]["passed"] = False
    else:
        record[
            {"energy": "energy_error_kj_mol", "force": "force_error_kj_mol_angstrom"}[
                blocked
            ]
        ]["coulomb"] = 0.01
    with pytest.raises(PDMSOPLSError, match="diagnostic_static_gate_failed") as error:
        run_control(
            s, p, tmp_path / "blocked", engine=None, static_control=record, solver=""
        )
    assert isinstance(error.value, PDMSOPLSError)
    assert not (tmp_path / "blocked").exists()


def test_committed_workflow_keeps_author_scope_separate():
    root = Path(__file__).resolve().parents[1] / "docs/evidence/phase_4j49/audit"
    if not (root / "summary.json").exists():
        pytest.skip("J49 evidence not generated yet")
    summary = json.loads((root / "summary.json").read_text())
    assert not summary["author_profile_MD_executed"]
    assert summary["diagnostic_MD_executed"]
    assert summary["bounded_diagnostic_acceptance_gate"] == 0
    for dp in (0, 2, 4, 6):
        directory = root / f"workflow-dp{dp}"
        w = summary["diagnostic_workflows"][str(dp)]
        assert w["NVT"]["steps"] == 40
        assert w["checkpoint_bytes_equal"] and w["complete_RNG_state_equal"]
        assert (
            w["RNG_reconstructed_in_fresh_process"]
            and w["future_random_stream_probe_equal"]
        )
        rng = json.loads((directory / "rng-continuation.json").read_text())
        restored = np.random.Generator(np.random.PCG64(0))
        restored.bit_generator.state = rng["after_initial_velocity_generation"]
        assert (
            restored.bit_generator.random_raw(16).tolist()
            == rng["future_random_raw_16_probe"]
        )
        assert json.loads((directory / "static-verification.json").read_text())[
            "finite_difference"
        ]["passed"]
        assert (directory / "whole.restart").read_bytes() == (
            directory / "resumed.restart"
        ).read_bytes()
        assert (directory / "rng-initial.json").read_bytes() == (
            directory / "rng-continuation.json"
        ).read_bytes()
        assert not w["canonical_author_MD_authority"]


def test_missing_parameter_row_cannot_reach_numerical_audit(assets, monkeypatch):
    from island.forcefields import pdms_opls

    original = pdms_opls.load_source

    def incomplete(options):
        source = original(options)
        source["rows"] = [r for r in source["rows"] if r["family"] != "angle"]
        return source

    monkeypatch.setattr(pdms_opls, "load_source", incomplete)
    system = fixture()
    with use_opls_assets(assets):
        p = prepare_polymer_forcefield(system, forcefield="opls")
    assert p.payload["status"] == "parameter_incomplete"
    assert any(r["family"] == "angle" for r in p.payload["missing_terms"])
    with pytest.raises(OPLSResolutionError):
        bind_inputs(system, p)
