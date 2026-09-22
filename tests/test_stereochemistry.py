import pytest
from rdkit import Chem

from island.builders import build_linear_polymer, build_polymer_from_sequence
from island.chemistry import (
    generate_stereochemical_sequence,
    inspect_repeat_unit_stereochemistry,
    invert_repeat_unit_stereochemistry,
    parse_psmiles,
    to_rdkit,
)
from island.exceptions import (
    PolymerBuildError,
    TacticityError,
    UnsupportedStereochemistryError,
)

STEREOGENIC_REPEAT = "[*:1]N[C@H](F)C[*:2]"
MULTICENTER_REPEAT = "[*:1]N[C@H](F)[C@H](Cl)C[*:2]"


def actual_controlled_cip_sequence(system: object) -> list[str]:
    converted = to_rdkit(system)
    mol = converted.mol
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centers: list[tuple[int, str]] = []
    for site_id, rdkit_index in converted.site_id_to_rdkit_index.items():
        site = system.topology.get_site(site_id)
        if not site.metadata.get("controllable_stereocenter"):
            continue
        atom = mol.GetAtomWithIdx(rdkit_index)
        assert atom.HasProp("_CIPCode")
        centers.append((site.metadata["repeat_unit_index"], atom.GetProp("_CIPCode")))
    return [state for _, state in sorted(centers)]


def test_bare_string_polymer_sequence_is_rejected() -> None:
    with pytest.raises(PolymerBuildError, match="bare string"):
        build_polymer_from_sequence(
            {"A": "[*]CC[*]", "B": "[*]CO[*]"},
            sequence="ABA",
            generate_3d=False,
        )


def test_single_backbone_stereocenter_is_detected() -> None:
    repeat = parse_psmiles(STEREOGENIC_REPEAT)
    inspection = inspect_repeat_unit_stereochemistry(repeat)
    assert len(inspection.relevant_centers) == 1
    center = inspection.controllable_center
    assert center is not None
    assert center.cip_label == "R"
    assert center.atom_index in inspection.backbone_path
    assert repeat.mol.GetAtomWithIdx(center.atom_index).GetAtomicNum() != 0


def test_repeat_unit_inversion_flips_graph_stereochemistry_and_preserves_attachments() -> (
    None
):
    repeat = parse_psmiles(STEREOGENIC_REPEAT)
    mirrored = invert_repeat_unit_stereochemistry(repeat)
    original_center = inspect_repeat_unit_stereochemistry(repeat).controllable_center
    mirrored_center = inspect_repeat_unit_stereochemistry(mirrored).controllable_center
    assert original_center is not None and mirrored_center is not None
    assert original_center.cip_label == "R"
    assert mirrored_center.cip_label == "S"
    assert mirrored.head == repeat.head
    assert mirrored.tail == repeat.tail
    assert [bond.GetBondType() for bond in mirrored.mol.GetBonds()] == [
        bond.GetBondType() for bond in repeat.mol.GetBonds()
    ]
    assert [atom.GetFormalCharge() for atom in mirrored.mol.GetAtoms()] == [
        atom.GetFormalCharge() for atom in repeat.mol.GetAtoms()
    ]


def test_achiral_repeat_unit_rejects_explicit_tacticity() -> None:
    with pytest.raises(UnsupportedStereochemistryError, match="has none"):
        build_linear_polymer("[*]CC[*]", dp=3, tacticity="isotactic", generate_3d=False)


def test_multiple_backbone_stereocenters_are_unsupported() -> None:
    inspection = inspect_repeat_unit_stereochemistry(parse_psmiles(MULTICENTER_REPEAT))
    assert len(inspection.relevant_centers) == 2
    with pytest.raises(UnsupportedStereochemistryError, match="found 2"):
        build_linear_polymer(
            MULTICENTER_REPEAT, dp=3, tacticity="isotactic", generate_3d=False
        )


def test_unassigned_backbone_stereocenter_is_unsupported() -> None:
    with pytest.raises(UnsupportedStereochemistryError, match="unassigned"):
        build_linear_polymer(
            "[*:1]NC(F)C[*:2]",
            dp=3,
            tacticity="isotactic",
            generate_3d=False,
        )


@pytest.mark.parametrize(
    ("tacticity", "expected"),
    [
        ("isotactic", ["R", "R", "R", "R", "R", "R"]),
        ("syndiotactic", ["R", "S", "R", "S", "R", "S"]),
    ],
)
def test_regular_tacticity_matches_actual_rdkit_cip_sequence(
    tacticity: str, expected: list[str]
) -> None:
    system = build_linear_polymer(
        STEREOGENIC_REPEAT, dp=6, tacticity=tacticity, generate_3d=False
    )
    assert system.metadata["polymer"]["stereochemical_sequence"] == expected
    assert actual_controlled_cip_sequence(system) == expected


def test_atactic_sequence_is_reproducible_and_seeded_independently() -> None:
    first = build_linear_polymer(
        STEREOGENIC_REPEAT,
        dp=20,
        tacticity="atactic",
        stereo_seed=11,
        random_seed=1,
        generate_3d=False,
    )
    second = build_linear_polymer(
        STEREOGENIC_REPEAT,
        dp=20,
        tacticity="atactic",
        stereo_seed=11,
        random_seed=999,
        generate_3d=False,
    )
    different = build_linear_polymer(
        STEREOGENIC_REPEAT,
        dp=20,
        tacticity="atactic",
        stereo_seed=12,
        random_seed=1,
        generate_3d=False,
    )
    first_sequence = first.metadata["polymer"]["stereochemical_sequence"]
    assert first_sequence == second.metadata["polymer"]["stereochemical_sequence"]
    assert first_sequence != different.metadata["polymer"]["stereochemical_sequence"]
    assert first_sequence.count("S") == 10
    assert actual_controlled_cip_sequence(first) == first_sequence


def test_stereo_seed_is_independent_of_chemical_sequence_generation_metadata() -> None:
    system = build_polymer_from_sequence(
        {"A": STEREOGENIC_REPEAT},
        ["A"] * 8,
        tacticity="atactic",
        stereo_seed=44,
        random_seed=99,
        generate_3d=False,
    )
    polymer = system.metadata["polymer"]
    assert polymer["sequence"] == ["A"] * 8
    assert polymer["stereo_seed"] == 44
    assert polymer["coordinates"] == "rdkit_2d"


def test_stereochemical_provenance_is_stored_per_repeat_and_atom() -> None:
    system = build_linear_polymer(
        STEREOGENIC_REPEAT, dp=5, tacticity="syndiotactic", generate_3d=False
    )
    expected = ["R", "S", "R", "S", "R"]
    polymer = system.metadata["polymer"]
    assert polymer["tacticity"] == "syndiotactic"
    assert polymer["stereochemical_sequence"] == expected
    assert polymer["stereochemistry_convention"] == "final_graph_absolute_cip"
    for site in system.topology.sites.values():
        repeat_index = site.metadata["repeat_unit_index"]
        assert (
            site.metadata["repeat_unit_stereochemical_state"] == expected[repeat_index]
        )
    controlled = [
        site
        for site in system.topology.sites.values()
        if site.metadata.get("controllable_stereocenter")
    ]
    assert len(controlled) == 5
    assert [site.metadata["cip_label"] for site in controlled] == expected


def test_mixed_repeat_units_reject_tacticity_control() -> None:
    with pytest.raises(
        UnsupportedStereochemistryError, match="single repeat-unit type"
    ):
        build_polymer_from_sequence(
            {"A": STEREOGENIC_REPEAT, "B": "[*]CC[*]"},
            ["A", "B", "A"],
            tacticity="isotactic",
            generate_3d=False,
        )


def test_stereochemical_sequence_generation_validation_and_fraction() -> None:
    sequence = generate_stereochemical_sequence(
        5,
        "atactic",
        reference_state="S",
        stereo_seed=9,
        atactic_fraction=0.4,
    )
    assert sequence.states.count("R") == 2
    assert sequence.states.count("S") == 3
    with pytest.raises(TacticityError, match="tacticity"):
        generate_stereochemical_sequence(3, "invalid", reference_state="R")
    with pytest.raises(TacticityError, match="between 0 and 1"):
        generate_stereochemical_sequence(
            3, "atactic", reference_state="R", atactic_fraction=1.1
        )
