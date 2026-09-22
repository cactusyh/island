from copy import deepcopy
from itertools import combinations

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from island.builders import build_linear_polymer, build_polymer_from_sequence
from island.chemistry import to_rdkit
from island.conformations import (
    ConnectionBondLengthPolicy,
    LocalTemplateConformationGenerator,
)
from island.conformations.sterics import build_excluded_pairs
from island.exceptions import (
    ConformationGenerationError,
    PolymerBuildError,
    UnsupportedConformationError,
)


def chemistry_snapshot(system):
    result = deepcopy(system.to_dict())
    result.pop("coordinates")
    polymer = result["metadata"]["polymer"]
    polymer.pop("coordinates", None)
    polymer.pop("coordinate_generation", None)
    return result


def coordinate_array(system):
    return np.array(
        [system.coordinates.get(site_id) for site_id in system.topology.sites]
    )


def geometry_cip(system):
    converted = to_rdkit(system)
    Chem.RemoveStereochemistry(converted.mol)
    Chem.AssignStereochemistryFrom3D(converted.mol, replaceExistingTags=True)
    values = []
    for site_id, atom_index in converted.site_id_to_rdkit_index.items():
        site = system.topology.get_site(site_id)
        if not site.metadata.get("controllable_stereocenter"):
            continue
        atom = converted.mol.GetAtomWithIdx(atom_index)
        values.append(
            (
                site.metadata["repeat_unit_index"],
                atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
            )
        )
    return [state for _, state in sorted(values)]


@pytest.mark.parametrize("dp", [50, 100])
def test_pe_long_chain_never_embeds_the_full_polymer(monkeypatch, dp):
    embedded_sizes = []
    original = AllChem.EmbedMolecule

    def instrumented(mol, *args, **kwargs):
        embedded_sizes.append(mol.GetNumAtoms())
        return original(mol, *args, **kwargs)

    monkeypatch.setattr(AllChem, "EmbedMolecule", instrumented)
    system = build_linear_polymer(
        "[*]CC[*]",
        dp=dp,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    generation = system.metadata["polymer"]["coordinate_generation"]
    assert embedded_sizes
    assert max(embedded_sizes) == generation["maximum_embedded_atom_count"]
    assert max(embedded_sizes) < system.number_of_sites
    assert generation["full_polymer_embedded"] is False
    assert generation["maximum_embedded_atom_count"] <= 14
    assert np.all(np.isfinite(coordinate_array(system)))


def test_local_templates_preserve_authoritative_chemistry_ids_and_provenance():
    reference = build_linear_polymer("[*]CC[*]", dp=12, generate_3d=False)
    generated = build_linear_polymer(
        "[*]CC[*]", dp=12, coordinate_method="local_templates"
    )
    assert chemistry_snapshot(generated) == chemistry_snapshot(reference)
    assert list(generated.topology.sites) == list(reference.topology.sites)
    assert generated.topology.bonds == reference.topology.bonds


def test_local_template_geometry_and_sterics_are_valid():
    system = build_linear_polymer(
        "[*]CC[*]",
        dp=20,
        coordinate_method="local_templates",
        template_seed=7,
        assembly_seed=8,
    )
    for bond in system.topology.bonds.values():
        distance = np.linalg.norm(
            system.coordinates.get(bond.site1) - system.coordinates.get(bond.site2)
        )
        assert 0.7 < distance < 2.1
    for center in system.topology.sites:
        for left, right in combinations(sorted(system.topology.neighbors(center)), 2):
            first = system.coordinates.get(left) - system.coordinates.get(center)
            second = system.coordinates.get(right) - system.coordinates.get(center)
            cosine = np.dot(first, second) / (
                np.linalg.norm(first) * np.linalg.norm(second)
            )
            angle = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
            assert 10.0 < angle < 179.9
    excluded = build_excluded_pairs(system)
    minimum = float("inf")
    site_ids = list(system.topology.sites)
    for index, site1 in enumerate(site_ids):
        for site2 in site_ids[index + 1 :]:
            if tuple(sorted((site1, site2))) in excluded:
                continue
            minimum = min(
                minimum,
                np.linalg.norm(
                    system.coordinates.get(site1) - system.coordinates.get(site2)
                ),
            )
    assert minimum >= 0.8


def test_heteroatom_and_sequence_defined_copolymer_are_supported():
    hetero = build_linear_polymer("[*]CO[*]", dp=8, coordinate_method="local_templates")
    assert any(site.element == "O" for site in hetero.topology.sites.values())
    copolymer = build_polymer_from_sequence(
        {"A": "[*:1]CC[*:2]", "B": "[*:2]CO[*:1]"},
        ["A", "B", "A", "B", "B", "A"],
        coordinate_method="local_templates",
        template_seed=31,
        assembly_seed=44,
    )
    assert copolymer.metadata["polymer"]["sequence"] == [
        "A",
        "B",
        "A",
        "B",
        "B",
        "A",
    ]
    assert np.all(np.isfinite(coordinate_array(copolymer)))


@pytest.mark.parametrize(
    ("tacticity", "expected"),
    [
        ("isotactic", ["R"] * 6),
        ("syndiotactic", ["R", "S", "R", "S", "R", "S"]),
    ],
)
def test_controlled_stereochemistry_is_validated_from_coordinates(tacticity, expected):
    system = build_linear_polymer(
        "[*:1]N[C@H](F)C[*:2]",
        dp=6,
        tacticity=tacticity,
        coordinate_method="local_templates",
        template_seed=10,
        assembly_seed=20,
    )
    assert system.metadata["polymer"]["stereochemical_sequence"] == expected
    assert geometry_cip(system) == expected


def test_template_and_assembly_seeds_are_reproducible_and_independent():
    kwargs = {
        "psmiles": "[*]CC[*]",
        "dp": 15,
        "coordinate_method": "local_templates",
        "template_seed": 101,
    }
    first = build_linear_polymer(**kwargs, assembly_seed=202)
    same = build_linear_polymer(**kwargs, assembly_seed=202)
    different = build_linear_polymer(**kwargs, assembly_seed=203)
    assert np.allclose(coordinate_array(first), coordinate_array(same))
    assert not np.allclose(coordinate_array(first), coordinate_array(different))
    assert chemistry_snapshot(first) == chemistry_snapshot(different)


def test_local_template_option_conflicts_and_failures_are_clear(monkeypatch):
    with pytest.raises(PolymerBuildError, match="conflicts"):
        build_linear_polymer(
            "[*]CC[*]",
            dp=3,
            generate_3d=False,
            coordinate_method="local_templates",
        )

    source = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=False)
    before = source.to_dict()
    monkeypatch.setattr(AllChem, "EmbedMolecule", lambda *args, **kwargs: -1)
    with pytest.raises(ConformationGenerationError, match="template embedding"):
        LocalTemplateConformationGenerator().generate(source)
    assert source.to_dict() == before


def test_invalid_junction_length_policy_is_rejected_without_mutation():
    class InvalidLength(ConnectionBondLengthPolicy):
        def length(self, system, site1, site2, bond_order):
            return float("nan")

    source = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=False)
    before = source.to_dict()
    with pytest.raises(UnsupportedConformationError, match="Invalid connection"):
        LocalTemplateConformationGenerator(bond_length_policy=InvalidLength()).generate(
            source
        )
    assert source.to_dict() == before


def test_aromatic_ring_geometry_is_rigid_with_local_templates():
    system = build_linear_polymer(
        "[*]c1ccccc1[*]",
        dp=4,
        coordinate_method="local_templates",
        template_seed=5,
        assembly_seed=9,
    )
    ring_ids = [
        site.id
        for site in system.topology.sites.values()
        if site.metadata["repeat_unit_index"] == 2 and site.element == "C"
    ]
    points = np.array([system.coordinates.get(site_id) for site_id in ring_ids])
    centered = points - points.mean(axis=0)
    assert np.linalg.svd(centered)[1][-1] < 1e-5
