"""Source closure is scoped: an engine definition cannot authorize author MD."""

import copy
import json
import os
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import pytest

from island.forcefields.opls_registry import OPLSAssets, OPLSResolutionError
from island.forcefields.selection.inventory import canonical, digest
from island.opls_conventions import LAYERS, MODES, POLICIES, REQUIRED_SETTINGS, signed
from island.opls_resolution import (
    EVIDENCE_CLASSES,
    REMAINING_SOURCE_GAPS,
    RESOLUTIONS,
    TEN_CONVENTIONS,
    ConventionSources,
    prepare_polymer_forcefield,
    reconcile,
    reject_conflicts,
    use_convention_sources,
)
from island.opls_resolved_execution import (
    LAMMPSOptions,
    generate_lammps,
    reconstruct,
    run_md,
    run_static,
)
from tests.test_phase_4j50 import ARCHIVE, ARCHIVE_SHA256, LMP, LMP_PIN, fixture

PAPER = Path("island_ffref/pcff/jp4c08471.pdf")
SUPPLEMENT = Path("island_ffref/pcff/jp4c08471_si_001.pdf")


@pytest.fixture(scope="module")
def prepared():
    if not all(p.is_file() for p in (ARCHIVE, PAPER, SUPPLEMENT)):
        pytest.skip("Explicit local source assets unavailable")
    assets = OPLSAssets(ARCHIVE, ARCHIVE_SHA256)
    sources = ConventionSources(PAPER, SUPPLEMENT)
    system = fixture(2)
    with use_convention_sources(sources):
        p = prepare_polymer_forcefield(system, forcefield="opls", assets=assets)
    return system, p


def test_canonical_selection_and_four_unchanged_layers(prepared):
    s, p = prepared
    p.validate_integrity(s)
    value = p.payload
    assert value["family"] == "opls"
    assert value["selected_profile"] == p.assignment.payload["selected_profile"]
    assert value["charge_origin"] == p.assignment.payload["charge_origin"]
    for layer in LAYERS:
        assert value[layer] == p.prior.payload[layer]
    assert not value["research_MD_available"]
    assert not value["production_validated"]
    assert value["simulation_readiness"] == "not_established"
    with pytest.raises(FrozenInstanceError):
        p.json_text = "{}"
    value["family"] = "other"
    assert p.payload["family"] == "opls"


def test_ten_explicit_gaps_and_separate_backend_resolution(prepared):
    _, p = prepared
    report = p.payload
    assert report["ten_conventions"] == list(TEN_CONVENTIONS)
    assert {r["setting"] for r in report["unresolved_conventions"]} == set(
        REMAINING_SOURCE_GAPS
    )
    assert {r["setting"] for r in report["rows"]} == set(REQUIRED_SETTINGS)
    for r in report["rows"]:
        assert r["resolution"] in RESOLUTIONS
        assert r["backend"]["resolution"] == "backend_resolved"
        assert r["evidence_checksum"]
        assert r["units"] and r["conversion_formula"] and r["conversion_implementation"]
        if r["setting"] in REMAINING_SOURCE_GAPS:
            assert r["value"] is None and r["resolution"] == "unresolved"
            assert not r["source_authoritative"]
            assert r["missing_evidence"] and r["absence_inventory_sha256"]
            assert r["backend"]["value"] is not None
            assert not r["backend"]["source_authoritative"]
            assert all(
                c["authority"] == "backend_only" for c in r["backend"]["evidence"]
            )
    assert report["source_scope_review"]["status"] == "unresolved"
    assert (
        report["source_scope_review"]["archive_pair_style"]
        != report["source_scope_review"]["SI_Table_S4_pair_style"]
    )


@pytest.mark.parametrize("setting", REQUIRED_SETTINGS)
def test_every_convention_mutation_rejected_even_after_resigning(
    prepared, setting, tmp_path
):
    system, p = prepared
    report = p.payload
    r = next(r for r in report["rows"] if r["setting"] == setting)
    r["value"] = "mutated"
    r["resolution"] = "source_resolved"
    r["identity"] = digest({k: v for k, v in r.items() if k != "identity"})
    report = signed({k: v for k, v in report.items() if k != "identity"})
    with pytest.raises(OPLSResolutionError, match="mutated_convention_resolution"):
        generate_lammps(
            system, replace(p, json_text=canonical(report)), tmp_path / "out"
        )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("setting", TEN_CONVENTIONS)
def test_backend_resolution_cannot_be_promoted_to_source(prepared, setting, tmp_path):
    system, p = prepared
    report = p.payload
    r = next(r for r in report["rows"] if r["setting"] == setting)
    r["backend"]["value"] = "changed"
    r["source_status"] = "source_resolved"
    report["unresolved_conventions"] = []
    report = signed({k: v for k, v in report.items() if k != "identity"})
    with pytest.raises(OPLSResolutionError):
        run_md(
            system,
            replace(p, json_text=canonical(report)),
            tmp_path / "out",
            research_opt_in=True,
        )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "kind,status",
    list(
        zip(
            EVIDENCE_CLASSES,
            (
                "source_resolved",
                "source_resolved",
                "source_resolved",
                "converter_resolved",
                "backend_resolved",
                "derived_resolved",
                "unresolved",
            ),
        )
    ),
)
def test_explicit_resolution_result(kind, status):
    claim = {
        "value": None if kind == "unresolved" else "test_only",
        "evidence_class": kind,
        "evidence": [{"test_only": True}],
    }
    r = reconcile("test_only", [claim])
    assert r["resolution"] == status


@pytest.mark.parametrize("setting", REQUIRED_SETTINGS)
def test_source_converter_conflict_preserves_both_and_rejects(setting):
    claims = [
        {
            "value": "source value",
            "evidence_class": "source_code_or_input",
            "evidence": [{"file": "test-only-source", "raw_text": "source value"}],
        },
        {
            "value": "different converter value",
            "evidence_class": "converter_defined",
            "evidence": [
                {"file": "test-only-converter", "raw_text": "different value"}
            ],
        },
    ]
    r = reconcile(setting, claims)
    assert r == reconcile(setting, list(reversed(claims)))
    assert r["resolution"] == "conflicting" and r["value"] is None
    with pytest.raises(OPLSResolutionError, match="conflicting_conventions") as exc:
        reject_conflicts([r])
    assert exc.value.diagnostics[0]["conflicts"][0]["candidates"] == r["candidates"]


@pytest.mark.parametrize("source", ["paper", "supplement", "archive"])
def test_source_file_replacement_rejected(prepared, source, tmp_path):
    system, p = prepared
    if source == "archive":
        path = tmp_path / "source.zip"
        path.write_bytes(ARCHIVE.read_bytes() + b"mutation")
        p = replace(
            p,
            prior=replace(
                p.prior,
                assignment=replace(
                    p.assignment, assets=replace(p.assignment.assets, archive_path=path)
                ),
            ),
        )
    else:
        path = tmp_path / "source.pdf"
        path.write_bytes(b"replacement")
        p = replace(p, sources=replace(p.sources, **{source: path}))
    with pytest.raises(OPLSResolutionError):
        generate_lammps(system, p, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mutation", ["graph", "charge", "order", "box", "family"])
def test_exact_input_identities(prepared, mutation, tmp_path):
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
        a = replace(p.assignment, request=replace(p.assignment.request, family="pcff"))
        p = replace(p, prior=replace(p.prior, assignment=a))
    else:
        value = p.assignment.payload
        value["graph_charges"]["vector"][0][1] = "0.001"
        a = replace(
            p.assignment,
            json_text=canonical(
                signed({k: v for k, v in value.items() if k != "identity"})
            ),
        )
        p = replace(p, prior=replace(p.prior, assignment=a))
    with pytest.raises(OPLSResolutionError):
        generate_lammps(system, p, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("action", ["static", "minimize", "nvt", "npt"])
def test_execution_policy_before_engine_or_outputs(prepared, mode, action, tmp_path):
    system, original = prepared
    p = prepare_polymer_forcefield(
        system,
        forcefield="opls",
        assets=original.assignment.assets,
        sources=original.sources,
        execution_mode=mode,
    )
    if mode == "diagnostic" and action == "static":
        assert p.require_execution(system, action)
        return
    with pytest.raises(OPLSResolutionError) as exc:
        # A sentinel engine would fail if any code tried to inspect it.
        if action == "static":
            run_static(
                system, p, tmp_path / "out", engine=object(), research_opt_in=True
            )
        else:
            run_md(
                system,
                p,
                tmp_path / "out",
                engine=object(),
                action=action,
                research_opt_in=True,
            )
    if mode == "research_opt_in":
        assert exc.value.diagnostics[0]["stage"] == "source_conventions_unresolved"
        assert {r["setting"] for r in exc.value.diagnostics[0]["missing"]} == set(
            REMAINING_SOURCE_GAPS
        )
    assert not (tmp_path / "out").exists()


def test_opt_in_is_separate_and_missing_assets_are_classified(prepared, tmp_path):
    s, p = prepared
    q = prepare_polymer_forcefield(
        s,
        forcefield="opls",
        assets=p.assignment.assets,
        sources=p.sources,
        execution_mode="research_opt_in",
    )
    with pytest.raises(OPLSResolutionError, match="explicit_research_opt_in_required"):
        run_md(s, q, tmp_path / "out")
    with pytest.raises(
        OPLSResolutionError, match="explicit_convention_sources_required"
    ):
        prepare_polymer_forcefield(s)
    with pytest.raises(OPLSResolutionError, match="convention_source_unavailable"):
        replace(p.sources, paper=tmp_path / "absent.pdf").validate()


@pytest.mark.parametrize("dp", [0, 2, 4, 6])
def test_deterministic_static_relocation_separate_process(prepared, dp, tmp_path):
    if not LMP.is_file():
        pytest.skip("Pinned KSPACE engine unavailable")
    _, original = prepared
    s = fixture(dp)
    p = prepare_polymer_forcefield(
        s,
        forcefield="opls",
        assets=original.assignment.assets,
        sources=original.sources,
    )
    assert generate_lammps(s, p, tmp_path / "one") == generate_lammps(
        s, p, tmp_path / "two"
    )
    result = run_static(s, p, tmp_path / "run", engine=LAMMPSOptions(LMP, LMP_PIN))
    assert result["resolution_identity"] == p.identity
    assert not result["MD_authorized"] and result["observed"]["forces_kj_mol_angstrom"]
    moved = tmp_path / "moved"
    shutil.copytree(tmp_path / "run", moved)
    for path, name in (
        (ARCHIVE, "source.zip"),
        (PAPER, "paper.pdf"),
        (SUPPLEMENT, "si.pdf"),
    ):
        shutil.copyfile(path, tmp_path / name)
    script = """import sys
from island.opls_resolution import ConventionSources
from island.opls_resolved_execution import reconstruct
from island.forcefields.opls_registry import OPLSAssets
_,p=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],assets=OPLSAssets(sys.argv[3],sys.argv[4]),sources=ConventionSources(sys.argv[5],sys.argv[6]))
print(p.identity)
"""
    for seed in (1, 99397):
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                script,
                str(moved),
                result["manifest_sha256"],
                str(tmp_path / "source.zip"),
                ARCHIVE_SHA256,
                str(tmp_path / "paper.pdf"),
                str(tmp_path / "si.pdf"),
            ],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
        )
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.strip() == p.identity
    (moved / "resolution.json").write_text("{}")
    with pytest.raises(OPLSResolutionError, match="workflow_file_mutation"):
        reconstruct(
            moved,
            manifest_sha256=result["manifest_sha256"],
            assets=original.assignment.assets,
            sources=original.sources,
        )


@pytest.mark.parametrize("policy", POLICIES)
def test_backend_selection_does_not_resolve_source_gap(prepared, policy):
    s, p = prepared
    q = prepare_polymer_forcefield(
        s,
        forcefield="opls",
        assets=p.assignment.assets,
        sources=p.sources,
        backend_policy=policy,
    )
    assert p.assignment.json_text == q.assignment.json_text
    assert len(q.payload["unresolved_conventions"]) == 9


def test_malformed_reconstruction_classified(tmp_path):
    with pytest.raises(OPLSResolutionError):
        reconstruct(
            tmp_path, manifest_sha256="0" * 64, assets=OPLSAssets(), sources=None
        )


def test_no_uncited_values_or_unknown_classes():
    with pytest.raises(OPLSResolutionError, match="uncited_convention_claim"):
        reconcile(
            "test",
            [{"value": 1, "evidence_class": "source_code_or_input", "evidence": []}],
        )
    with pytest.raises(OPLSResolutionError, match="unknown_convention_evidence_class"):
        reconcile("test", [{"value": 1, "evidence_class": "common_default"}])


def test_tail_is_proven_by_style_restriction_not_default(prepared):
    _, p = prepared
    tail = next(r for r in p.payload["rows"] if r["setting"] == "lj_tail_correction")
    assert tail["resolution"] == "converter_resolved" and tail["value"] == "no"
    assert not tail["source_authoritative"]
    assert "cannot use" in tail["conversion_formula"]
    assert any(
        "None of the" in c["raw_text"] for c in tail["candidates"][0]["evidence"]
    )


def test_report_json_order_independent(prepared):
    _, p = prepared
    data = p.payload
    assert json.loads(canonical(data)) == json.loads(p.json_text)
