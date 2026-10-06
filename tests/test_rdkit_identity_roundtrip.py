"""Chemical identity round trips; synthetic PCFF source used only in software tests."""

from copy import deepcopy

import pytest

pytest.importorskip("rdkit")
from test_pcff_automatic import synthetic  # noqa: F401
from test_pcff_domains import domains_source  # noqa: F401

from island.chemistry import from_rdkit, from_smiles, to_rdkit
from island.exceptions import RDKitConversionError


@pytest.mark.parametrize("smiles", ["[2H]O[2H]", "[CH3]", "CCO"])
def test_isotope_radical_ordinary_roundtrip(smiles):
    system = from_smiles(smiles, random_seed=2026)
    before = deepcopy(system)
    graph = to_rdkit(system)
    restored = from_rdkit(
        graph.mol,
        site_ids=[
            graph.rdkit_index_to_site_id[i] for i in range(graph.mol.GetNumAtoms())
        ],
    ).system
    for i, atom in system.topology.sites.items():
        other = restored.topology.sites[i]
        assert other.metadata.get("isotope", 0) == atom.metadata.get("isotope", 0)
        assert other.metadata.get("radical_electrons", 0) == atom.metadata.get(
            "radical_electrons", 0
        )
        assert other.mass == atom.mass
    assert system.topology.sites == before.topology.sites
    assert system.topology.bonds == before.topology.bonds


def test_pcff_deuterium_after_roundtrip(domains_source):  # noqa: F811
    from island.forcefields.pcff import type_pcff_atoms

    system = from_smiles("[2H]O[2H]", random_seed=2026)
    restored = from_rdkit(to_rdkit(system).mol).system
    before = type_pcff_atoms(
        system, domains_source, profile="island_pcff_source_graph_v3"
    )
    after = type_pcff_atoms(
        restored, domains_source, profile="island_pcff_source_graph_v3"
    )
    assert after.complete
    assert after.assignments == before.assignments
    assert sorted(after.assignments.values()) == ["dw", "dw", "o*"]
    assert after.identity == before.identity


@pytest.mark.parametrize("key", ["isotope", "radical_electrons"])
@pytest.mark.parametrize(
    "value", [True, False, -1, 1.5, "2", None, float("nan"), 2**40]
)
def test_malformed_identity_metadata_rejected(key, value):
    s = from_smiles("CC", random_seed=2026)
    s.topology.sites[1].metadata[key] = value
    with pytest.raises(RDKitConversionError, match=key):
        to_rdkit(s)


def test_no_isotope_inference_or_input_mass_repair():
    s = from_smiles("O", random_seed=2026)
    hydrogen = next(a for a in s.topology.sites.values() if a.element == "H")
    hydrogen.mass = 2.014
    assert to_rdkit(s).mol.GetAtomWithIdx(hydrogen.id - 1).GetIsotope() == 0
    assert hydrogen.mass == 2.014
