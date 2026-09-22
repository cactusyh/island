"""Internal final-graph tetrahedral stereochemistry assignment."""

from collections.abc import Sequence

from rdkit import Chem

from island.exceptions import UnsupportedStereochemistryError


def assign_cip_sequence(
    mol: Chem.Mol, atom_indices: Sequence[int], states: Sequence[str]
) -> None:
    """Assign and verify requested absolute CIP states on a complete graph."""
    if len(atom_indices) != len(states):
        raise UnsupportedStereochemistryError(
            "Stereocenter indices and stereochemical states must have equal length"
        )
    for atom_index, target in zip(atom_indices, states, strict=True):
        atom = mol.GetAtomWithIdx(atom_index)
        atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if _cip_label(atom) == target:
            continue
        atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if _cip_label(atom) != target:
            raise UnsupportedStereochemistryError(
                f"Final polymer atom {atom_index} cannot realize requested "
                f"{target} stereochemistry"
            )
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    observed = tuple(_cip_label(mol.GetAtomWithIdx(index)) for index in atom_indices)
    if observed != tuple(states):
        raise UnsupportedStereochemistryError(
            "Final polymer stereochemistry does not match the requested sequence: "
            f"requested {tuple(states)}, observed {observed}"
        )


def verify_cip_sequence(
    mol: Chem.Mol, atom_indices: Sequence[int], states: Sequence[str]
) -> None:
    """Verify an already assigned CIP sequence, for example after adding H atoms."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    observed = tuple(_cip_label(mol.GetAtomWithIdx(index)) for index in atom_indices)
    if observed != tuple(states):
        raise UnsupportedStereochemistryError(
            "Polymer stereochemistry changed during graph processing: "
            f"requested {tuple(states)}, observed {observed}"
        )


def _cip_label(atom: Chem.Atom) -> str | None:
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
