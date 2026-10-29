"""Source authority is independent of backend choice, numerical success and opt-in."""

import os
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace

import pytest

from island.forcefields.opls_registry import OPLSAssets, OPLSResolutionError
from island.forcefields.opls_registry.registry import ARCHIVE_SHA256
from island.forcefields.selection.inventory import canonical, digest
from island.opls_conventions import (
    CLASSIFICATIONS,
    LAYERS,
    MODES,
    PROFILE,
    REQUIRED_SETTINGS,
    signed,
)
from island.opls_execution import (
    LAMMPSOptions,
    generate_lammps,
    input_text,
    prepare_polymer_forcefield,
    reconstruct,
    run_md,
    run_static,
)
from tests.test_phase_4j48 import ARCHIVE, LMP, LMP_PIN, fixture


@pytest.fixture(scope="module")
def assets():
    if not ARCHIVE.is_file():
        pytest.skip("Explicit local research archive unavailable")
    return OPLSAssets(ARCHIVE, ARCHIVE_SHA256)


@pytest.fixture(scope="module")
def prepared(assets):
    system = fixture(2)
    return system, prepare_polymer_forcefield(system, forcefield="opls", assets=assets)


def test_canonical_call_and_four_separate_immutable_records(prepared):
    system, p = prepared
    assert p.validate_integrity(system)
    report = p.payload
    assert report["family"] == "opls" and report["selected_profile"] == PROFILE
    assert report["charge_origin"] == p.assignment.payload["charge_origin"]
    assert report["domain"] == p.assignment.payload["domain"]
    assert report["assignment_identity"] == p.assignment.identity
    assert report["identity"] == digest(
        {k: v for k, v in report.items() if k != "identity"}
    )
    assert set(p.records()) == set(LAYERS)
    for name, record in p.records().items():
        data = record.payload
        assert data["kind"] == name
        assert record.identity == digest(
            {k: v for k, v in data.items() if k != "identity"}
        )
        with pytest.raises(FrozenInstanceError):
            record.json_text = "{}"
        data["kind"] = "changed"
        assert record.payload["kind"] == name
    assert not report["default_MD_enabled"]
    assert not report[LAYERS[3]]["current_graph_numerically_verified_by_preparation"]
    assert report[LAYERS[3]]["prior_verified_cases_only"] == [
        "HMDS",
        "capped_DP2",
        "capped_DP4",
        "capped_DP6",
    ]
    assert not report["production_validated"]
    assert report["simulation_readiness"] == "not_established"


def test_source_defaults_are_not_promoted_or_inferred(prepared):
    _, p = prepared
    source = p.payload[LAYERS[0]]
    backend = p.payload[LAYERS[2]]
    assert set(source["settings"]) == set(REQUIRED_SETTINGS)
    for name, r in source["settings"].items():
        assert r["classification"] in CLASSIFICATIONS
        if r["value"] is None:
            assert r["classification"] == "unresolved"
            assert r["missing_evidence"] and r["absence_evidence_sha256"]
            assert (
                backend["settings"][name]["classification"]
                == "selected_by_backend_policy"
            )
            assert not backend["settings"][name]["source_authority"]
    assert source["settings"]["dielectric"]["value"] is None
    boundary = source["settings"]["boundary_conditions"]
    assert boundary["value"] == "p p p"
    assert boundary["classification"] == "explicitly_stated_by_source"
    assert boundary["evidence"][0]["source_medium"] == "PDF"
    assert boundary["evidence"][0]["page"] == 3
    assert (
        source["settings"]["special_bonds"]["evidence"][0]["raw_text"]
        == "special_bonds lj/coul 0.0 0.0 0.5"
    )
    assert source["settings"]["reciprocal_tolerance"]["value"] == "0.0001"
    assert backend["settings"]["dielectric"]["value"] == "1"


@pytest.mark.parametrize("layer", LAYERS)
def test_resigning_any_layer_cannot_authorize_execution(prepared, layer):
    system, p = prepared
    value = p.payload
    value[layer]["authority"] = "native_authorized"
    value[layer] = signed({k: v for k, v in value[layer].items() if k != "identity"})
    value = signed({k: v for k, v in value.items() if k != "identity"})
    with pytest.raises(OPLSResolutionError, match="stale_or_mutated"):
        replace(p, json_text=canonical(value)).validate_integrity(system)


@pytest.mark.parametrize(
    "setting",
    [
        "lammps_units",
        "special_bonds",
        "pair_cutoff_angstrom",
        "reciprocal_tolerance",
        "dielectric",
        "lj_mixing",
        "torsion_convention",
        "pressure_convention",
    ],
)
@pytest.mark.parametrize("layer", [LAYERS[0], LAYERS[2]])
def test_setting_mutations_even_resigned_fail_before_output(
    prepared, setting, layer, tmp_path
):
    system, p = prepared
    value = p.payload
    value[layer]["settings"][setting]["value"] = "tampered"
    value[layer] = signed({k: v for k, v in value[layer].items() if k != "identity"})
    value = signed({k: v for k, v in value.items() if k != "identity"})
    with pytest.raises(OPLSResolutionError, match="stale_or_mutated"):
        generate_lammps(
            system, replace(p, json_text=canonical(value)), tmp_path / "out"
        )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("action", ["minimize", "nvt", "npt"])
def test_modes_never_bypass_source_or_scientific_gates(
    prepared, mode, action, tmp_path
):
    system, p = prepared
    # Changing a requested mode is legitimate, but never changing its authority.
    p = prepare_polymer_forcefield(
        system, assets=p.assignment.assets, execution_mode=mode
    )

    class NoEngine:
        def validate(self):
            pytest.fail("blocked preparation reached an engine")

    with pytest.raises(OPLSResolutionError) as exc:
        run_md(
            system,
            p,
            tmp_path / "blocked",
            engine=NoEngine(),
            research_opt_in=True,
            action=action,
        )
    diagnostics = exc.value.diagnostics[0]
    if mode == "research_opt_in":
        assert diagnostics["stage"] == "source_conventions_unresolved"
        assert "dielectric" in {r["setting"] for r in diagnostics["missing"]}
        assert all(r["missing_evidence"] for r in diagnostics["missing"])
    elif mode == "diagnostic":
        assert diagnostics["stage"] == "diagnostic_mode_does_not_authorize_MD"
    else:
        assert diagnostics["stage"] == "production_candidate_unavailable"
    assert not (tmp_path / "blocked").exists()


def test_research_requires_separate_explicit_opt_in(prepared, tmp_path):
    system, p = prepared
    p = prepare_polymer_forcefield(
        system, assets=p.assignment.assets, execution_mode="research_opt_in"
    )
    with pytest.raises(OPLSResolutionError, match="explicit_research_opt_in_required"):
        generate_lammps(system, p, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_unresolved_inventory_cannot_be_deleted_to_enable_MD(prepared, tmp_path):
    system, p = prepared
    value = p.payload
    value["unresolved_settings"] = []
    value["research_MD_available"] = True
    value[LAYERS[3]]["research_MD_authorized"] = True
    value = signed({k: v for k, v in value.items() if k != "identity"})
    with pytest.raises(OPLSResolutionError):
        generate_lammps(
            system, replace(p, json_text=canonical(value)), tmp_path / "out"
        )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "mutation", ["source_archive", "graph", "charge", "order", "box", "profile"]
)
def test_chemical_or_source_mutation_cannot_reuse_policy(
    prepared, assets, mutation, tmp_path
):
    import copy

    original, p = prepared
    system = copy.deepcopy(original)
    if mutation == "source_archive":
        archive = tmp_path / "changed.zip"
        archive.write_bytes(ARCHIVE.read_bytes() + b"changed")
        p = replace(
            p,
            assignment=replace(
                p.assignment, assets=replace(assets, archive_path=archive)
            ),
        )
    elif mutation == "graph":
        system.topology.bonds.pop(next(iter(system.topology.bonds)))
    elif mutation == "order":
        system.topology.sites = dict(reversed(list(system.topology.sites.items())))
    elif mutation == "charge":
        value = p.assignment.payload
        value["graph_charges"]["vector"][0][1] = "0.001"
        value = signed({k: v for k, v in value.items() if k != "identity"})
        p = replace(p, assignment=replace(p.assignment, json_text=canonical(value)))
    elif mutation == "profile":
        p = replace(
            p,
            assignment=replace(
                p.assignment, request=replace(p.assignment.request, family="pcff")
            ),
        )
    else:
        from island.core import SimulationBox

        system.box = SimulationBox(41, 41, 41)
    with pytest.raises(OPLSResolutionError):
        generate_lammps(system, p, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("family", ["pcff", "compass", "gaff", "gaff2"])
def test_wrong_family_has_no_policy_adapter(family):
    with pytest.raises(OPLSResolutionError, match="OPLS_family_required"):
        prepare_polymer_forcefield(None, forcefield=family)


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_actual_static_deterministic_export_relocated_fresh_process(
    assets, dp, tmp_path
):
    if not LMP.is_file():
        pytest.skip("Explicit pinned KSPACE executable unavailable")
    system = fixture(dp)
    p = prepare_polymer_forcefield(system, forcefield="opls", assets=assets)
    one, two = tmp_path / "one", tmp_path / "two"
    assert generate_lammps(system, p, one) == generate_lammps(system, p, two)
    text = input_text(system, p)
    for literal in [
        "dielectric 1",
        "special_bonds lj/coul 0 0 0.5",
        "pair_modify mix geometric",
        "kspace_style pppm 0.0001",
        "lj/charmm/coul/long 9 11",
        "dihedral_style opls",
    ]:
        assert literal in text
    result = run_static(system, p, tmp_path / "run", engine=LAMMPSOptions(LMP, LMP_PIN))
    assert not result["MD_executed"] and result["observed"]["forces_kj_mol_angstrom"]
    assert not result["source_conventions_resolved"]
    moved = tmp_path / "relocated"
    shutil.copytree(tmp_path / "run", moved)
    s, q = reconstruct(moved, manifest_sha256=result["manifest_sha256"], assets=assets)
    assert q.identity == p.identity and input_text(s, q) == text
    script = """import sys
from island.opls_execution import reconstruct
from island.forcefields.opls_registry import OPLSAssets
s,p=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],assets=OPLSAssets(sys.argv[3],sys.argv[4]))
print(p.identity)
"""
    archive = tmp_path / "moved.zip"
    shutil.copyfile(ARCHIVE, archive)
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            script,
            str(moved),
            result["manifest_sha256"],
            str(archive),
            ARCHIVE_SHA256,
        ],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONHASHSEED": "99397"},
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == p.identity
    (moved / "lammps/system.in").write_text(
        text.replace("dielectric 1", "dielectric 2")
    )
    with pytest.raises(OPLSResolutionError, match="workflow_file_mutation"):
        reconstruct(moved, manifest_sha256=result["manifest_sha256"], assets=assets)


def test_backend_policy_is_not_a_new_physical_profile(prepared):
    system, original = prepared
    p = prepare_polymer_forcefield(
        system,
        assets=original.assignment.assets,
        backend_policy="exact_erfc_diagnostic_v1",
    )
    assert p.assignment.json_text == original.assignment.json_text
    assert p.payload[LAYERS[0]] == original.payload[LAYERS[0]]
    assert p.payload[LAYERS[1]] == original.payload[LAYERS[1]]
    assert p.identity != original.identity
    b = p.payload[LAYERS[2]]
    assert (
        b["settings"]["kspace_style"]["classification"] == "selected_by_backend_policy"
    )
    assert b["settings"]["reciprocal_tolerance"]["source_value"] == "0.0001"
    assert "table 0" in input_text(system, p)
    assert "kspace_style ewald 0.0000001" in input_text(system, p)
    assert not p.payload["research_MD_available"]


def test_wrong_pinned_engine_rejected_before_output(prepared, tmp_path):
    system, p = prepared
    with pytest.raises(OPLSResolutionError, match="backend_policy_executable_mismatch"):
        run_static(
            system, p, tmp_path / "out", engine=LAMMPSOptions("absent", "0" * 64)
        )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "setting",
    [
        "dielectric",
        "special_bonds",
        "lj_mixing",
        "reciprocal_tolerance",
        "pair_cutoff_angstrom",
        "torsion_convention",
        "pressure_convention",
    ],
)
def test_unresolved_selected_backend_setting_has_no_implicit_default(
    prepared, setting, monkeypatch, tmp_path
):
    import island.opls_conventions as conventions

    system, p = prepared
    original = conventions._backend

    def incomplete_backend(source, policy):
        value = original(source, policy)
        value["settings"][setting] = {
            "value": None,
            "classification": "unresolved",
            "evidence": [],
            "missing_evidence": "test of incomplete numerical policy",
        }
        return signed({k: v for k, v in value.items() if k != "identity"})

    monkeypatch.setattr(conventions, "_backend", incomplete_backend)
    incomplete = prepare_polymer_forcefield(system, assets=p.assignment.assets)
    with pytest.raises(OPLSResolutionError) as error:
        generate_lammps(system, incomplete, tmp_path / "out")
    assert error.value.diagnostics[0]["stage"] == "unresolved_backend_settings"
    assert error.value.diagnostics[0]["missing"][0]["setting"] == setting
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "options",
    [
        {"execution_mode": "auto_MD"},
        {"backend_policy": "implicit_defaults"},
        {"execution_mode": None},
    ],
)
def test_unknown_mode_or_policy_rejected(options):
    with pytest.raises(OPLSResolutionError, match="unknown_execution_mode_or_policy"):
        prepare_polymer_forcefield(None, **options)


@pytest.mark.parametrize("present", [False, True])
def test_malformed_reconstruction_is_classified(tmp_path, present):
    from island.pdms_lammps import pin

    expected = "0" * 64
    if present:
        (tmp_path / "manifest.json").write_text("{}")
        expected = pin(tmp_path / "manifest.json")
    with pytest.raises(OPLSResolutionError):
        reconstruct(tmp_path, manifest_sha256=expected, assets=OPLSAssets())
