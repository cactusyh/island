"""Preflight never promises coverage; caches must retain chemistry and source identity."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from rdkit import Chem

from island.exceptions import ForceFieldError
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_diagnostics import (
    assignment_failure,
    preflight_assignment,
)


@pytest.fixture
def mock_source(monkeypatch, tmp_path):
    from island.forcefields import assignment_diagnostics, emc_graph

    monkeypatch.setattr(
        assignment_diagnostics,
        "_field_sources",
        lambda *a: (tmp_path, tmp_path, tmp_path, {"mock": "test-only"}),
    )
    monkeypatch.setattr(
        emc_graph,
        "_field_elements",
        lambda *a: {
            "carbon": "C",
            "hydrogen": "H",
            "oxygen": "O",
            "nitrogen": "N",
            "chloride": "Cl",
        },
    )


def test_preflight_does_not_embed_or_claim_assignment(mock_source, monkeypatch):
    from island.forcefields import _emc_graph_input

    monkeypatch.setattr(
        _emc_graph_input,
        "_embed",
        lambda *a: pytest.fail("preflight attempted embedding"),
    )
    report = preflight_assignment(normalize(PolymerSpec((("A", "*CCO*"),), 3)), "opls")
    assert report["ready_for_assignment"]
    assert report["parameter_completeness"] == report["md"] == "unmeasured"
    assert report["atoms"] == 23
    assert len(report["environments"]) == report["atoms"]


def test_preflight_missing_element_and_ua_isotope(mock_source):
    si = preflight_assignment(Chem.MolFromSmiles("C[Si](C)(C)C"), "pcff")
    assert not si["ready_for_assignment"]
    assert si["blockers"] == [{"code": "element_absent_from_source", "sites": [1]}]
    ua = preflight_assignment(Chem.MolFromSmiles("[13CH3]C"), "opls/2024/opls-ua")
    assert ua["blockers"] == [{"code": "united_atom_isotope_unsupported", "sites": [0]}]


def test_ion_pair_warns_without_inventing_parameters(mock_source):
    report = preflight_assignment(Chem.MolFromSmiles("C[N+](C)(C)C.[Cl-]"), "opls")
    assert report["formal_charge"] == 0
    assert report["component_charges"] == [1, -1]
    assert report["warnings"] and report["ready_for_assignment"]
    assert report["parameter_completeness"] == "unmeasured"


def test_failure_localization_is_backend_evidence_only():
    graph = normalize(PolymerSpec((("A", "*CCO*"),), 2))
    exc = ForceFieldError(
        "increment pair [a:b] not found; no angle coefficients found for [a:b:c]"
    )
    row = assignment_failure(exc, graph)
    assert set(row["categories"]) == {"charge_increment", "angle"}
    assert not row["source_sites"] and not row["environments"]
    exc.source_sites = [3]
    row = assignment_failure(exc, graph)
    assert row["environments"][0]["element"] == "O"
    assert row["environments"][0]["site_id"] == 3
    assert row["localization"] == "backend_reported_sites"


def test_closed_shell_preflight_rejects_radicals(mock_source):
    with pytest.raises(ForceFieldError, match="closed-shell"):
        preflight_assignment(Chem.MolFromSmiles("[CH3]"), "opls")


def test_checksum_stable_across_integer_key_json_roundtrip():
    from island.forcefields._assignment_cache import digest

    data = {"site_groups": {i: [i] for i in range(1, 30)}}
    assert digest(data) == digest(json.loads(json.dumps(data)))


@pytest.fixture
def installed_field(tmp_path):
    pytest.importorskip("pyemc")
    from island.forcefields.emc_profiles import install

    field = tmp_path / "field"
    install("opls-polymer-v1", field)
    return field


def test_cache_hit_tamper_and_source_mutation(installed_field, tmp_path, monkeypatch):
    from island.builders import build_linear_polymer
    from island.forcefields import emc_graph

    graph = build_linear_polymer("*CC*", dp=2, generate_3d=False)
    cache = tmp_path / "cache"
    original = emc_graph._build_script
    calls = []

    def tracked(*args, **kw):
        calls.append(1)
        return original(*args, **kw)

    monkeypatch.setattr(emc_graph, "_build_script", tracked)
    first = emc_graph.assign_emc_graph(
        graph, "opls", field_dir=installed_field, cache_dir=cache
    )
    second = emc_graph.assign_emc_graph(
        graph, "opls", field_dir=installed_field, cache_dir=cache
    )
    assert first == second and len(calls) == 1
    source = installed_field / "opls/2024/opls-aa.top"
    source.write_bytes(source.read_bytes() + b"\n")
    third = emc_graph.assign_emc_graph(
        graph, "opls", field_dir=installed_field, cache_dir=cache
    )
    assert len(calls) == 2
    assert third.field_sha256 != first.field_sha256
    assert len(list(cache.glob("*.json"))) == 2
    for path in cache.glob("*.json"):
        raw = json.loads(path.read_text())
        raw["result"]["parameters"] += "\ninvalid"
        path.write_text(json.dumps(raw))
    with pytest.raises(ForceFieldError, match="Invalid assignment cache"):
        emc_graph.assign_emc_graph(
            graph, "opls", field_dir=installed_field, cache_dir=cache
        )
    assert len(calls) == 2


def test_cache_key_includes_ids_geometry_stereo_and_seed():
    import numpy as np

    from island.forcefields._assignment_cache import digest, request

    def source(smiles):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        return SimpleNamespace(
            mol=mol,
            sites=list(range(mol.GetNumAtoms())),
            generated=(),
            xyz=np.zeros((mol.GetNumAtoms(), 3)),
            bounds=np.array([[0, 20]] * 3),
            coordinate_source="test",
        )

    first = source("C/C=C/C")
    a = request(first, "opls", 4, {"source": "a"})
    variants = [
        request(first, "opls", 5, {"source": "a"}),
        request(first, "opls", 4, {"source": "b"}),
        request(source(r"C/C=C\C"), "opls", 4, {"source": "a"}),
    ]
    first.xyz[0, 0] = 1
    variants.append(request(first, "opls", 4, {"source": "a"}))
    first.sites = list(reversed(first.sites))
    variants.append(request(first, "opls", 4, {"source": "a"}))
    assert len({digest(x) for x in [a, *variants]}) == 6


def test_failure_is_not_cached(installed_field, tmp_path, monkeypatch):
    from island.forcefields import emc_graph

    def failed(*a, **kw):
        raise ForceFieldError("missing parameter")

    monkeypatch.setattr(emc_graph, "_build_script", failed)
    cache = tmp_path / "cache"
    with pytest.raises(ForceFieldError, match="missing parameter"):
        emc_graph.assign_emc_graph(
            Chem.MolFromSmiles("CC"), "opls", field_dir=installed_field, cache_dir=cache
        )
    assert not list(cache.glob("*.json"))


def test_cache_cannot_modify_source_tree(installed_field):
    from island.forcefields import emc_graph

    with pytest.raises(ValueError, match="outside parameter"):
        emc_graph.assign_emc_graph(
            Chem.MolFromSmiles("CC"),
            "opls",
            field_dir=installed_field,
            cache_dir=installed_field / "cache",
        )


def test_source_change_during_assignment_is_not_cached(
    installed_field, tmp_path, monkeypatch
):
    from island.forcefields import emc_graph

    original = emc_graph._build_script

    def mutate(*a, **kw):
        result = original(*a, **kw)
        p = installed_field / "opls/2024/opls-aa.top"
        p.write_bytes(p.read_bytes() + b"\n")
        return result

    monkeypatch.setattr(emc_graph, "_build_script", mutate)
    cache = tmp_path / "cache"
    with pytest.raises(ForceFieldError, match="changed during assignment"):
        emc_graph.assign_emc_graph(
            Chem.MolFromSmiles("CC"), "opls", field_dir=installed_field, cache_dir=cache
        )
    assert not list(cache.glob("*.json"))


def test_capability_probe_uses_graph_default_and_retains_preflight(
    installed_field, tmp_path
):
    from island.forcefields.capabilities import probe

    report = probe(
        "*CCO*", "opls", dp=2, field_dir=installed_field, cache_dir=tmp_path / "cache"
    )
    assert report["assignment"] == "passed"
    assert report["preflight"]["parameter_completeness"] == "unmeasured"
    assert report["dynamics"] == "unmeasured"
    rejected = probe(
        "*CC*", "opls", dp=2, field_dir=installed_field, stereo_policy="polyjarvis"
    )
    assert rejected["assignment"] == "failed"
    assert "select emc_mode" in rejected["diagnostic"]
    legacy = probe(
        "*CC*",
        "opls",
        dp=2,
        field_dir=installed_field,
        emc_mode="repeat",
        stereo_policy="polyjarvis",
    )
    assert legacy["assignment"] == "passed"


def test_summary_retains_timeouts_and_unmatched_pairs():
    from scripts.benchmark_final_graphs import summarize

    rows = [
        {
            "id": "a",
            "panel": "fixed982",
            "field": "opls",
            "mode": "repeat",
            "assigned": True,
            "strict": False,
        },
        {
            "id": "a",
            "panel": "fixed982",
            "field": "opls",
            "mode": "graph",
            "assigned": False,
            "strict": False,
            "failure": {"categories": ["timeout"]},
        },
        {
            "id": "b",
            "panel": "fixed982",
            "field": "opls",
            "mode": "graph",
            "assigned": True,
            "strict": True,
        },
    ]
    summary = summarize(rows, rows + [{"unstarted": True, "panel": "fixed982"}])
    assert not summary["complete"]
    assert summary["paired"]["opls"] == {"repeat_only": 1}
    assert summary["groups"]["fixed982:opls:graph"]["attempted"] == 2
    assert summary["groups"]["fixed982:opls:graph"]["failure_categories"] == {
        "timeout": 1
    }


def test_fixed_panel_and_challenge_denominators():
    from scripts.benchmark_final_graphs import cases, challenges

    repo = Path(__file__).resolve().parents[1]
    assert len(cases(repo / "docs/evidence/reference_assignment.csv.gz")) == 982
    extra = challenges()
    assert len(extra) == 34 and len({row["id"] for row in extra}) == 34
    assert len([row for row in extra if "smiles" in row]) == 3


@pytest.mark.parametrize(
    "residue,smiles",
    [
        ("DMES", "CSC"),
        ("EMES", "CSCC"),
        ("MSME", "CS(C)(=O)=O"),
        ("ETBR", "CCBr"),
        ("ETCL", "CCCl"),
    ],
)
def test_extended_native_residue_types_and_charges(installed_field, residue, smiles):
    from island.forcefields.emc_graph import _field_elements, assign_emc_graph
    from scripts.check_opls_polymer_profile import source_atoms, typed_atoms

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    expected = source_atoms(
        (installed_field / "opls/2024/opls-aa.define").read_text(),
        residue,
        mol,
        _field_elements(installed_field, "opls/2024/opls-aa"),
    )
    actual = typed_atoms(assign_emc_graph(mol, "opls", field_dir=installed_field))
    assert actual.keys() == expected.keys()
    for i, (kind, charge) in expected.items():
        assert actual[i][0] == kind
        assert actual[i][1] == pytest.approx(charge, abs=1e-6)


def test_missing_source_torsion_remains_a_failure(installed_field):
    from island.forcefields.emc_graph import assign_emc_graph

    with pytest.raises(ForceFieldError, match="torsion coefficients") as caught:
        assign_emc_graph(
            Chem.MolFromSmiles("CS(C)=O"), "opls", field_dir=installed_field
        )
    failure = assignment_failure(caught.value)
    assert "torsion" in failure["categories"]
    assert any("c4s:s3:c4s:h1" in t for t in failure["parameter_tuples"])


def test_timeout_keeps_completed_preflight_and_denominator(tmp_path, monkeypatch):
    import signal
    import subprocess

    from scripts import benchmark_final_graphs as benchmark

    job = {"id": "slow", "panel": "challenge", "field": "opls", "mode": "graph"}
    killed = []

    class Process:
        pid = 123456789

        def __init__(self, argv, **kwargs):
            self.waited = False
            request = Path(argv[-1])
            benchmark.write(
                request.with_name(
                    request.name.replace(".input.json", ".progress.json")
                ),
                {**job, "preflight": {"atoms": 10}, "graph_valid": True},
            )

        def wait(self, timeout=None):
            if not self.waited:
                self.waited = True
                raise subprocess.TimeoutExpired("test-only-worker", timeout)
            return -9

    monkeypatch.setattr(benchmark.subprocess, "Popen", Process)
    monkeypatch.setattr(
        benchmark.os, "killpg", lambda pid, sig: killed.append((pid, sig))
    )
    result = benchmark.execute(job, tmp_path, 0.01)
    assert killed == [(Process.pid, signal.SIGKILL)]
    assert not result["assigned"] and not result["strict"]
    assert result["preflight"] == {"atoms": 10}
    assert result["failure"]["categories"] == ["timeout"]


def test_complete_summary_rejects_duplicate_rows():
    from scripts.benchmark_final_graphs import summarize

    a = {
        "id": "a",
        "panel": "fixed982",
        "field": "opls",
        "mode": "graph",
        "assigned": True,
        "strict": True,
    }
    b = {**a, "id": "b"}
    assert not summarize([a, a], [a, b])["complete"]


@pytest.mark.parametrize("dp", [1, 2, 4])
@pytest.mark.parametrize(
    "repeat", ["*n1cccc1CC*", "*c1ccc(N2C(=O)c3cc4C(=O)N(*)C(=O)c4cc3C2=O)cc1"]
)
def test_aromatic_n_termination_matches_explicit_hydrogen(dp, repeat):
    from island.builders import build_linear_polymer
    from island.chemistry import to_rdkit

    default = to_rdkit(build_linear_polymer(repeat, dp=dp, generate_3d=False)).mol
    capped = to_rdkit(
        build_linear_polymer(
            repeat,
            dp=dp,
            generate_3d=False,
            head_end_group="*[H]",
            tail_end_group="*[H]",
        )
    ).mol
    assert Chem.MolToSmiles(Chem.RemoveHs(default)) == Chem.MolToSmiles(
        Chem.RemoveHs(capped)
    )
    assert not any(a.GetNumRadicalElectrons() for a in default.GetAtoms())
    assert any(
        a.GetIsAromatic()
        and a.GetAtomicNum() == 7
        and any(n.GetAtomicNum() == 1 for n in a.GetNeighbors())
        for a in default.GetAtoms()
    )


def test_legacy_benchmark_is_not_blocked_by_island_graph_failure(monkeypatch):
    from island.exceptions import PolymerBuildError
    from island.forcefields import assignment, assignment_core
    from scripts.benchmark_final_graphs import run_case

    called = []

    def fail(*a, **kw):
        raise PolymerBuildError("test-only graph failure")

    def native(*a, **kw):
        called.append(kw["emc_mode"])
        return SimpleNamespace(
            bond_stereo_omitted=False,
            field_sha256={},
            emc_version="test-only",
            data="\nAtoms\n\n1 1 1 0 0 0 0 # c4\n\n",
            parameters="test-only",
            emc_smiles="*CC*",
        )

    monkeypatch.setattr(assignment_core, "normalize", fail)
    monkeypatch.setattr(assignment, "assign_polymer_parameters", native)
    result = run_case(
        {
            "id": "test",
            "psmiles": "*CC*",
            "dp": 4,
            "end_groups": [None, None],
            "field": "opls",
            "field_dir": "unused",
            "seed": 1,
            "mode": "repeat",
        }
    )
    assert result["assigned"] and called == ["repeat"]
    assert not result["graph_valid"] and result["preflight_failure"]


def test_ring_scan_retries_only_explicit_depth_failure(installed_field, monkeypatch):
    from island.forcefields import emc_graph

    original = emc_graph._build_script
    depths = []

    def scan(*args, **kw):
        depths.append(args[-1])
        if len(depths) == 1:
            raise ForceFieldError("increase scanning depth?")
        return original(*args, **kw)

    monkeypatch.setattr(emc_graph, "_build_script", scan)
    result = emc_graph.assign_emc_graph(
        Chem.MolFromSmiles("CC"), "opls", field_dir=installed_field
    )
    assert depths == [False, True] and result.automatic_ring_depth


def test_normal_scan_matches_automatic_scan_coefficients(installed_field, monkeypatch):
    from island.epoxy.chemistry import monomer
    from island.forcefields import emc_graph

    mol = monomer("COc1ccccc1", 173)
    bounded = emc_graph.assign_emc_graph(mol, "opls", field_dir=installed_field)
    original = emc_graph._build_script

    def auto(*args, **kw):
        args = list(args)
        args[-1] = True
        return original(*args, **kw)

    monkeypatch.setattr(emc_graph, "_build_script", auto)
    reference = emc_graph.assign_emc_graph(mol, "opls", field_dir=installed_field)
    assert bounded.parameters == reference.parameters
    for section in ("Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"):
        assert emc_graph._section(bounded.data, section) == emc_graph._section(
            reference.data, section
        )
    assert not bounded.automatic_ring_depth
