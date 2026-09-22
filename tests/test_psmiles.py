import pytest
from rdkit import Chem

from island.chemistry import RepeatUnit, parse_psmiles
from island.exceptions import InvalidRepeatUnitError, PSMILESError


def test_unlabeled_attachments_follow_original_parse_order() -> None:
    repeat = parse_psmiles("[*]CO[*]")
    assert isinstance(repeat, RepeatUnit)
    assert repeat.original_psmiles == "[*]CO[*]"
    assert repeat.head.role == "head"
    assert repeat.tail.role == "tail"
    assert repeat.head.dummy_atom_index < repeat.tail.dummy_atom_index
    assert repeat.mol.GetAtomWithIdx(repeat.head.neighbor_atom_index).GetSymbol() == "C"
    assert repeat.mol.GetAtomWithIdx(repeat.tail.neighbor_atom_index).GetSymbol() == "O"


def test_explicit_labels_determine_asymmetric_orientation() -> None:
    repeat = parse_psmiles("[*:2]CO[*:1]")
    assert repeat.head.atom_map_number == 1
    assert repeat.tail.atom_map_number == 2
    assert repeat.mol.GetAtomWithIdx(repeat.head.neighbor_atom_index).GetSymbol() == "O"
    assert repeat.mol.GetAtomWithIdx(repeat.tail.neighbor_atom_index).GetSymbol() == "C"


@pytest.mark.parametrize(
    ("psmiles", "count"),
    [("CC", 0), ("[*]CC", 1), ("[*]C([*])[*]", 3)],
)
def test_repeat_unit_requires_exactly_two_dummies(psmiles: str, count: int) -> None:
    with pytest.raises(InvalidRepeatUnitError, match=rf"found {count}"):
        parse_psmiles(psmiles)


@pytest.mark.parametrize("psmiles", ["[*:1]CC[*]", "[*:1]CC[*:1]", "[*:3]CC[*:4]"])
def test_ambiguous_or_conflicting_labels_are_rejected(psmiles: str) -> None:
    with pytest.raises(InvalidRepeatUnitError, match="labels are ambiguous"):
        parse_psmiles(psmiles)


def test_dummy_to_dummy_structure_is_rejected() -> None:
    with pytest.raises(InvalidRepeatUnitError, match="connect to real"):
        parse_psmiles("[*][*]")


def test_dummy_with_multiple_neighbors_is_rejected() -> None:
    mol = Chem.MolFromSmiles("C[*](C)C[*]", sanitize=False)
    assert mol is not None
    with pytest.raises(InvalidRepeatUnitError, match="exactly one neighbor"):
        parse_psmiles("C[*](C)C[*]")


def test_chemically_invalid_psmiles_is_rejected() -> None:
    with pytest.raises(PSMILESError, match="invalid PSMILES"):
        parse_psmiles("not a smiles")


def test_attachment_bond_information_is_preserved() -> None:
    repeat = parse_psmiles("[*]=CC=[*]")
    assert repeat.head.bond_type == Chem.BondType.DOUBLE
    assert repeat.tail.bond_order == 2.0
