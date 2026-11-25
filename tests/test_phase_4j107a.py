"""Equivalent backend representations preserve graph identity and bounded failures."""

import os
import sys

import pytest
from rdkit import Chem

from island.exceptions import ForceFieldError
from island.forcefields._emc_graph_input import emc_groups
from island.forcefields.emc import _run
from island.forcefields.emc_graph import assign_emc_graph


@pytest.mark.parametrize(
    "traversal", ["input", "canonical", "reverse", "permuted-1", "permuted-2"]
)
def test_traversal_retains_source_bonds_and_components(traversal):
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)c1ccc2ccccc2c1.CCO"))
    smiles, orders = emc_groups(mol, False, False, traversal=traversal)
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    observed = set()
    sites = []
    for name, text in smiles.items():
        parsed = Chem.MolFromSmiles(text, parser)
        order = orders[name]
        sites.extend(order)
        assert parsed.GetNumAtoms() == len(order)
        for a, original in zip(parsed.GetAtoms(), order, strict=True):
            assert a.GetAtomicNum() == mol.GetAtomWithIdx(original).GetAtomicNum()
        for b in parsed.GetBonds():
            observed.add(
                (
                    tuple(
                        sorted((order[b.GetBeginAtomIdx()], order[b.GetEndAtomIdx()]))
                    ),
                    b.GetBondType(),
                )
            )
    expected = {
        (tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))), b.GetBondType())
        for b in mol.GetBonds()
    }
    assert observed == expected
    assert sorted(sites) == list(range(mol.GetNumAtoms()))


def test_backend_timeout_is_explicit(tmp_path):
    with pytest.raises(ForceFieldError, match="backend timed out"):
        _run(
            [sys.executable, "-c", "import time; time.sleep(10)"],
            tmp_path,
            os.environ,
            timeout=0.05,
        )


@pytest.mark.parametrize("timeout", [True, -1, 0, float("nan"), float("inf")])
def test_timeout_validation_before_backend(timeout):
    with pytest.raises(ValueError, match="backend_timeout"):
        assign_emc_graph(Chem.MolFromSmiles("CC"), "pcff", backend_timeout=timeout)


def test_geometry_failure_retains_stage(monkeypatch):
    from island.forcefields import _emc_graph_input
    from island.forcefields.assignment_diagnostics import assignment_failure

    def fail(*args):
        raise ForceFieldError("ETKDG test failure")

    monkeypatch.setattr(_emc_graph_input, "graph_input", fail)
    with pytest.raises(ForceFieldError) as caught:
        assign_emc_graph(Chem.MolFromSmiles("CC"), "pcff")
    report = assignment_failure(caught.value)
    assert report["stage"] == "input_geometry" and report["elapsed_seconds"] >= 0


@pytest.fixture
def field(tmp_path):
    pytest.importorskip("pyemc")
    from island.forcefields.emc_profiles import install

    path = tmp_path / "field"
    install("opls-polymer-v1", path)
    return path


def test_retry_retains_failed_attempt_and_exact_output(field, monkeypatch):
    from island.forcefields import emc_graph

    mol = Chem.MolFromSmiles("COc1ccccc1")
    reference = assign_emc_graph(
        mol, "opls", field_dir=field, serialization="canonical"
    )
    original = emc_graph._build_script
    calls = []

    def fail_first(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise ForceFieldError("Inconsistent connectivity in rings")
        return original(*args, **kwargs)

    monkeypatch.setattr(emc_graph, "_build_script", fail_first)
    result = assign_emc_graph(mol, "opls", field_dir=field)
    assert result.serialization == "canonical"
    assert result.data == reference.data and result.parameters == reference.parameters
    assert len(result.assignment_attempts) == 2
    assert not result.assignment_attempts[0]["passed"]
    assert result.assignment_attempts[-1]["passed"]
    assert result.timings_seconds["input_geometry"] >= 0


def test_runtime_failure_is_not_a_typing_retry(field, monkeypatch):
    from island.forcefields import emc_graph

    calls = []

    def fail(*args, **kwargs):
        calls.append(1)
        raise ForceFieldError("EMC backend timed out")

    monkeypatch.setattr(emc_graph, "_build_script", fail)
    with pytest.raises(ForceFieldError) as caught:
        assign_emc_graph(Chem.MolFromSmiles("c1ccccc1"), "opls", field_dir=field)
    assert len(calls) == 1
    assert caught.value.assignment_stage == "backend"


def test_explicit_serialization_participates_in_cache_identity(field, tmp_path):
    cache = tmp_path / "cache"
    mol = Chem.MolFromSmiles("COc1ccccc1")
    for traversal in ["input", "canonical"]:
        a = assign_emc_graph(
            mol, "opls", field_dir=field, serialization=traversal, cache_dir=cache
        )
        assert a.serialization == traversal
    assert len(list(cache.glob("*.json"))) == 2


def test_macrocycle_uses_graph_ring_size(field):
    mol = Chem.MolFromSmiles("C1CCCCCCCCCCC1")
    a = assign_emc_graph(mol, "opls", field_dir=field, backend_timeout=10)
    assert a.ring_scan_depth == 12 and not a.automatic_ring_depth
    assert len(a.site_to_atom_id) == Chem.AddHs(mol).GetNumAtoms()


def test_explicit_attempt_timeout_can_retry_equivalent_graph(field, monkeypatch):
    from island.forcefields import emc_graph

    original = emc_graph._build_script
    calls = []

    def first_timeout(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise ForceFieldError("EMC backend timed out after 1 seconds")
        return original(*args, **kwargs)

    monkeypatch.setattr(emc_graph, "_build_script", first_timeout)
    a = assign_emc_graph(
        Chem.MolFromSmiles("COc1ccccc1"), "opls", field_dir=field, backend_timeout=10
    )
    assert len(calls) == 2 and a.serialization == "canonical"
    assert len(a.assignment_attempts) == 2
