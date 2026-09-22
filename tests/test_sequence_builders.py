import pytest
from rdkit import Chem

from island.builders import (
    PolymerSequence,
    build_alternating_copolymer,
    build_block_copolymer,
    build_polymer_from_sequence,
    build_random_copolymer,
)
from island.chemistry import to_rdkit
from island.exceptions import PolymerBuildError

A = "[*:1]CC[*:2]"
B = "[*:1]CO[*:2]"


def assert_linear_chemical_system(system: object, number_of_units: int) -> None:
    assert len(system.topology.connected_components()) == 1
    assert all(site.atomic_number != 0 for site in system.topology.sites.values())
    cross_repeat_bonds = [
        bond
        for bond in system.topology.bonds.values()
        if system.topology.get_site(bond.site1).metadata["repeat_unit_index"]
        != system.topology.get_site(bond.site2).metadata["repeat_unit_index"]
    ]
    assert len(cross_repeat_bonds) == number_of_units - 1
    mol = to_rdkit(system).mol
    assert not any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
    Chem.SanitizeMol(mol)


def test_polymer_sequence_preserves_order_count_composition_and_metadata() -> None:
    sequence = PolymerSequence(
        ("A", "A", "C", "B", "A"), {"strategy": "explicit", "user": True}
    )
    assert sequence.identities == ("A", "A", "C", "B", "A")
    assert sequence.number_of_repeat_units == 5
    assert sequence.composition_counts == {"A": 3, "C": 1, "B": 1}
    assert sequence.generation_metadata == {"strategy": "explicit", "user": True}


def test_generic_a_b_a_sequence_is_reproduced_in_provenance() -> None:
    system = build_polymer_from_sequence(
        {"A": A, "B": B}, ["A", "B", "A"], generate_3d=False
    )
    polymer = system.metadata["polymer"]
    assert polymer["polymer_type"] == "sequence_defined"
    assert polymer["sequence"] == ["A", "B", "A"]
    assert polymer["composition_counts"] == {"A": 2, "B": 1}
    assert polymer["repeat_unit_definitions"] == {"A": A, "B": B}
    for repeat_index, expected_type in enumerate(("A", "B", "A")):
        observed = {
            site.metadata["repeat_unit_type"]
            for site in system.topology.sites.values()
            if site.metadata["repeat_unit_index"] == repeat_index
        }
        assert observed == {expected_type}
    assert_linear_chemical_system(system, 3)


def test_alternating_dp_is_total_repeat_unit_count() -> None:
    system = build_alternating_copolymer(A, B, dp=5, generate_3d=False)
    polymer = system.metadata["polymer"]
    assert polymer["polymer_type"] == "alternating_copolymer"
    assert polymer["sequence"] == ["A", "B", "A", "B", "A"]
    assert polymer["number_of_repeat_units"] == 5
    assert polymer["composition_counts"] == {"A": 3, "B": 2}
    assert_linear_chemical_system(system, 5)


def test_multiblock_a3_b2_a1_sequence() -> None:
    system = build_block_copolymer([(A, 3), (B, 2), (A, 1)], generate_3d=False)
    polymer = system.metadata["polymer"]
    assert polymer["polymer_type"] == "block_copolymer"
    assert polymer["sequence"] == ["A", "A", "A", "B", "B", "A"]
    assert polymer["composition_counts"] == {"A": 4, "B": 2}
    assert polymer["sequence_generation"]["blocks"] == [
        {"repeat_unit_type": "A", "count": 3},
        {"repeat_unit_type": "B", "count": 2},
        {"repeat_unit_type": "A", "count": 1},
    ]
    assert_linear_chemical_system(system, 6)


def test_random_sequence_and_composition_are_reproducible() -> None:
    kwargs = {
        "repeat_units": {"A": A, "B": B},
        "dp": 7,
        "fractions": {"A": 0.6, "B": 0.4},
        "sequence_seed": 41,
        "generate_3d": False,
    }
    first = build_random_copolymer(**kwargs, random_seed=1)
    second = build_random_copolymer(**kwargs, random_seed=999)
    first_metadata = first.metadata["polymer"]
    second_metadata = second.metadata["polymer"]
    assert first_metadata["sequence"] == second_metadata["sequence"]
    assert first_metadata["composition_counts"] == {"A": 4, "B": 3}
    assert first_metadata["sequence_generation"]["allocation_method"] == (
        "largest_remainder"
    )
    assert len(first_metadata["sequence"]) == 7
    assert_linear_chemical_system(first, 7)


def test_different_random_seeds_can_produce_different_sequences() -> None:
    kwargs = {
        "repeat_units": {"A": A, "B": B},
        "dp": 20,
        "fractions": {"A": 0.5, "B": 0.5},
        "generate_3d": False,
    }
    first = build_random_copolymer(**kwargs, sequence_seed=1)
    second = build_random_copolymer(**kwargs, sequence_seed=2)
    assert (
        first.metadata["polymer"]["sequence"] != second.metadata["polymer"]["sequence"]
    )


def test_each_repeat_units_labeled_orientation_is_respected() -> None:
    repeat_units = {
        "A": "[*:1]CN[*:2]",
        "B": "[*:2]CO[*:1]",
    }
    system = build_polymer_from_sequence(
        repeat_units, ["A", "B", "A"], generate_3d=False
    )
    polymer = system.metadata["polymer"]
    assert system.topology.get_site(polymer["head_site_id"]).element == "C"
    assert system.topology.get_site(polymer["tail_site_id"]).element == "N"

    connections: list[tuple[str, str]] = []
    for bond in system.topology.bonds.values():
        first = system.topology.get_site(bond.site1)
        second = system.topology.get_site(bond.site2)
        first_index = first.metadata["repeat_unit_index"]
        second_index = second.metadata["repeat_unit_index"]
        if first_index == second_index:
            continue
        if first_index < second_index:
            connections.append((first.element, second.element))
        else:
            connections.append((second.element, first.element))
    assert connections == [("N", "O"), ("C", "C")]


def test_every_atom_and_generated_hydrogen_has_repeat_unit_type() -> None:
    system = build_polymer_from_sequence(
        {"A": A, "B": B}, ["A", "B", "A"], generate_3d=False
    )
    for site in system.topology.sites.values():
        repeat_type = site.metadata["repeat_unit_type"]
        assert repeat_type in {"A", "B"}
        if site.element != "H":
            continue
        assert site.metadata["generated_hydrogen"] is True
        parent_id = next(iter(system.topology.neighbors(site.id)))
        parent = system.topology.get_site(parent_id)
        assert parent.metadata["repeat_unit_type"] == repeat_type
        assert (
            parent.metadata["repeat_unit_index"] == site.metadata["repeat_unit_index"]
        )


def test_explicit_sequence_rejects_empty_or_undefined_labels() -> None:
    with pytest.raises(PolymerBuildError, match="cannot be empty"):
        build_polymer_from_sequence({"A": A}, [], generate_3d=False)
    with pytest.raises(PolymerBuildError, match="undefined"):
        build_polymer_from_sequence({"A": A}, ["A", "B"], generate_3d=False)


@pytest.mark.parametrize(
    "fractions",
    [
        {"A": 0.7, "B": 0.2},
        {"A": 1.1, "B": -0.1},
        {"A": 1.0},
    ],
)
def test_random_fraction_validation(fractions: dict[str, float]) -> None:
    with pytest.raises(PolymerBuildError, match="Fraction|fractions"):
        build_random_copolymer(
            {"A": A, "B": B},
            dp=5,
            fractions=fractions,
            generate_3d=False,
        )


def test_largest_remainder_tie_uses_definition_order() -> None:
    system = build_random_copolymer(
        {"A": A, "B": B},
        dp=10,
        fractions={"A": 0.55, "B": 0.45},
        sequence_seed=7,
        generate_3d=False,
    )
    assert system.metadata["polymer"]["composition_counts"] == {"A": 6, "B": 4}
