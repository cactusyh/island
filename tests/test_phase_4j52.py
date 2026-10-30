"""Permanent absence is scoped to pinned evidence and cannot authorize MD."""

import copy
import json
import os
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace

import pytest

from island.forcefields.opls_registry import OPLSAssets, OPLSResolutionError
from island.forcefields.selection.inventory import canonical
from island.opls_conventions import LAYERS, MODES, POLICIES, REQUIRED_SETTINGS, signed
from island.opls_final_execution import (
    LAMMPSOptions,
    generate_lammps,
    prepare_polymer_forcefield,
    reconstruct,
    run_md,
    run_static,
)
from island.opls_final_resolution import (
    NINE_CONVENTIONS,
    RESOLUTIONS,
    ConventionSources,
    corpus_audit,
    use_convention_sources,
)
from island.polymer import use_opls_assets
from tests.test_phase_4j50 import ARCHIVE, ARCHIVE_SHA256, LMP, LMP_PIN, fixture
from tests.test_phase_4j51 import PAPER, SUPPLEMENT


@pytest.fixture(scope="module")
def prepared():
    if not all(p.is_file() for p in (ARCHIVE, PAPER, SUPPLEMENT)):
        pytest.skip("Explicit pinned source corpus unavailable")
    system = fixture(2)
    with (
        use_opls_assets(OPLSAssets(ARCHIVE, ARCHIVE_SHA256)),
        use_convention_sources(ConventionSources(PAPER, SUPPLEMENT)),
    ):
        preparation = prepare_polymer_forcefield(system, forcefield="opls")
    return system, preparation


def with_assignment(p, assignment):
    return replace(
        p, prior=replace(p.prior, prior=replace(p.prior.prior, assignment=assignment))
    )


def test_canonical_family_only_and_historical_identity(prepared):
    s, p = prepared
    p.validate_integrity(s)
    assert p.payload["family"] == "opls"
    assert p.payload["selected_profile"] == p.assignment.payload["selected_profile"]
    assert p.payload["charge_origin"] == p.assignment.payload["charge_origin"]
    assert p.payload["J51_resolution_identity"] == p.prior.identity
    assert all(p.payload[layer] == p.prior.payload[layer] for layer in LAYERS)
    assert p.payload["production_validated"] is False
    assert p.payload["simulation_readiness"] == "not_established"
    with pytest.raises(FrozenInstanceError):
        p.json_text = "{}"


def test_complete_corpus_with_absent_and_out_of_scope_matches():
    audit = corpus_audit()
    assert audit["no_network_discovery"]
    assert "new evidence requires a new audited version" in audit["closure_scope"]
    assert len(audit["assets"]) == 64
    assert sum(a["kind"] == "archive_member" for a in audit["assets"]) == 32
    for asset in audit["assets"]:
        assert asset["lines_searched"] and asset["source_sha256"]
        assert set(asset["matches"]) == set(NINE_CONVENTIONS)
    assert audit["prior_comparison_records_read"]
    # A different family's explicit dielectric is retained, never borrowed.
    dreiding = next(
        a
        for a in audit["assets"]
        if a["source_file"] == "param/Dreiding/dr_system.in.init"
    )
    assert dreiding["matches"]["dielectric"]
    assert dreiding["scope"] == "other_family_or_readme_not_OPLS_authority"
    opls = next(
        a
        for a in audit["assets"]
        if a["source_file"] == "param/OPLS/opls_system.in.init"
    )
    assert not opls["matches"]["dielectric"]
    audit["assets"].clear()
    assert corpus_audit()["assets"]


@pytest.mark.parametrize("setting", NINE_CONVENTIONS)
def test_exact_permanent_missing_artifact_and_separate_backend(prepared, setting):
    _, p = prepared
    rows = p.payload["final_resolutions"]
    assert [r["setting"] for r in rows] == list(NINE_CONVENTIONS)
    row = next(r for r in rows if r["setting"] == setting)
    assert row["resolution"] == "unresolved_permanently" in RESOLUTIONS
    assert row["value"] is None and not row["source_authoritative"]
    assert row["exact_missing_artifact"] and row["default_scientifically_unjustified"]
    assert row["affected_observables"] and row["affected_execution_modes"]
    assert len(row["all_searched_evidence"]) == 64
    assert row["evidence_checksum"] and row["independent_reconstruction"]
    assert row["backend"]["resolution"] == "backend_resolved"
    assert not row["backend"]["source_authoritative"]


@pytest.mark.parametrize("setting", NINE_CONVENTIONS)
def test_final_row_cannot_be_resigned_into_authority(prepared, setting, tmp_path):
    s, p = prepared
    value = p.payload
    row = next(r for r in value["final_resolutions"] if r["setting"] == setting)
    row["resolution"] = "source_resolved"
    row["value"] = row["backend"]["value"]
    value["permanent_blockers"] = value["unresolved_conventions"] = []
    value = signed({k: v for k, v in value.items() if k != "identity"})
    with pytest.raises(
        OPLSResolutionError, match="mutated_final_convention_resolution"
    ):
        generate_lammps(s, replace(p, json_text=canonical(value)), tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("setting", REQUIRED_SETTINGS)
def test_all_conventions_including_units_special_torsion_cutoff_kspace_mutation(
    prepared, setting, tmp_path
):
    s, p = prepared
    value = p.payload
    value[LAYERS[0]]["settings"][setting]["value"] = "mutated"
    value = signed({k: v for k, v in value.items() if k != "identity"})
    with pytest.raises(OPLSResolutionError):
        generate_lammps(s, replace(p, json_text=canonical(value)), tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("source", ["paper", "supplement", "archive"])
def test_original_source_or_author_archive_replacement(prepared, source, tmp_path):
    s, p = prepared
    path = tmp_path / "replacement"
    path.write_bytes(b"different source bytes")
    if source == "archive":
        p = with_assignment(
            p,
            replace(
                p.assignment, assets=replace(p.assignment.assets, archive_path=path)
            ),
        )
    else:
        p = replace(
            p, prior=replace(p.prior, sources=replace(p.sources, **{source: path}))
        )
    with pytest.raises(OPLSResolutionError):
        generate_lammps(s, p, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mutation", ["graph", "charge", "order", "box", "family"])
def test_graph_charge_order_box_and_family_bound(prepared, mutation, tmp_path):
    original, p = prepared
    system = copy.deepcopy(original)
    if mutation == "graph":
        system.topology.bonds.pop(next(iter(system.topology.bonds)))
    elif mutation == "order":
        system.topology.sites = dict(reversed(list(system.topology.sites.items())))
    elif mutation == "box":
        from island.core import SimulationBox

        system.box = SimulationBox(41, 40, 40)
    elif mutation == "family":
        p = with_assignment(
            p,
            replace(p.assignment, request=replace(p.assignment.request, family="pcff")),
        )
    else:
        value = p.assignment.payload
        value["graph_charges"]["vector"][0][1] = "0.001"
        p = with_assignment(
            p,
            replace(
                p.assignment,
                json_text=canonical(
                    signed({k: v for k, v in value.items() if k != "identity"})
                ),
            ),
        )
    with pytest.raises(OPLSResolutionError):
        generate_lammps(system, p, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("action", ["static", "minimize", "nvt", "npt"])
def test_modes_reject_before_engine_or_outputs(prepared, mode, action, tmp_path):
    s, original = prepared
    p = prepare_polymer_forcefield(
        s,
        forcefield="opls",
        assets=original.assignment.assets,
        sources=original.sources,
        execution_mode=mode,
    )
    if mode == "diagnostic" and action == "static":
        assert p.require_execution(s)
        return
    with pytest.raises(OPLSResolutionError) as exc:
        if action == "static":
            run_static(s, p, tmp_path / "out", engine=object(), research_opt_in=True)
        else:
            run_md(
                s,
                p,
                tmp_path / "out",
                engine=object(),
                action=action,
                research_opt_in=True,
            )
    if mode == "research_opt_in":
        assert (
            exc.value.diagnostics[0]["stage"]
            == "source_conventions_unresolved_permanently"
        )
        assert [r["setting"] for r in exc.value.diagnostics[0]["missing"]] == list(
            NINE_CONVENTIONS
        )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("policy", POLICIES)
def test_backend_experiment_does_not_close_source(prepared, policy):
    s, original = prepared
    p = prepare_polymer_forcefield(
        s,
        forcefield="opls",
        assets=original.assignment.assets,
        sources=original.sources,
        backend_policy=policy,
    )
    assert p.assignment.json_text == original.assignment.json_text
    assert len(p.payload["permanent_blockers"]) == 9


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_actual_static_deterministic_relocation_and_fresh_process(
    prepared, dp, tmp_path
):
    if not LMP.is_file():
        pytest.skip("Pinned KSPACE executable unavailable")
    _, original = prepared
    s = fixture(dp)
    p = prepare_polymer_forcefield(
        s,
        forcefield="opls",
        assets=original.assignment.assets,
        sources=original.sources,
    )
    assert generate_lammps(s, p, tmp_path / "a") == generate_lammps(
        s, p, tmp_path / "b"
    )
    result = run_static(s, p, tmp_path / "run", engine=LAMMPSOptions(LMP, LMP_PIN))
    assert result["final_resolution_identity"] == p.identity
    assert not result["MD_authorized"] and result["observed"]["forces_kj_mol_angstrom"]
    moved = tmp_path / "moved"
    shutil.copytree(tmp_path / "run", moved)
    for path, name in (
        (ARCHIVE, "source.zip"),
        (PAPER, "paper.pdf"),
        (SUPPLEMENT, "si.pdf"),
    ):
        shutil.copyfile(path, tmp_path / name)
    command = [
        sys.executable,
        "-m",
        "scripts.phase_4j52",
        "reconstruct",
        "--directory",
        str(moved),
        "--manifest-sha256",
        result["manifest_sha256"],
        "--archive",
        str(tmp_path / "source.zip"),
        "--archive-sha256",
        ARCHIVE_SHA256,
        "--paper",
        str(tmp_path / "paper.pdf"),
        "--supplement",
        str(tmp_path / "si.pdf"),
    ]
    for seed in (1, 99397):
        proc = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
        )
        assert proc.returncode == 0, proc.stderr
        assert json.loads(proc.stdout)["identity"] == p.identity
    (moved / "final-resolution.json").write_text("{}")
    with pytest.raises(OPLSResolutionError, match="workflow_file_mutation"):
        reconstruct(
            moved,
            manifest_sha256=result["manifest_sha256"],
            assets=original.assignment.assets,
            sources=original.sources,
        )
