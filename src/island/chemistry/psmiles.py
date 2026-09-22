"""PSMILES parsing for linear, bifunctional repeat units."""

from dataclasses import dataclass
from typing import Literal

from rdkit import Chem

from island.exceptions import InvalidRepeatUnitError, PSMILESError


@dataclass(frozen=True)
class AttachmentPoint:
    """A dummy atom and the real atom to which it is attached."""

    dummy_atom_index: int
    neighbor_atom_index: int
    bond_type: Chem.BondType
    bond_order: float
    role: Literal["head", "tail"]
    atom_map_number: int = 0


@dataclass(frozen=True)
class RepeatUnit:
    """A parsed, oriented, linear PSMILES repeat unit."""

    original_psmiles: str
    mol: Chem.Mol
    head: AttachmentPoint
    tail: AttachmentPoint

    @property
    def attachment_points(self) -> tuple[AttachmentPoint, AttachmentPoint]:
        return (self.head, self.tail)

    @classmethod
    def from_psmiles(cls, psmiles: str) -> "RepeatUnit":
        return parse_psmiles(psmiles)


def parse_psmiles(psmiles: str) -> RepeatUnit:
    """Parse a bifunctional linear repeat unit without reordering its atoms.

    Labels ``[*:1]`` and ``[*:2]`` mean head and tail, respectively. If neither
    dummy is labeled, their RDKit parse order determines head then tail.
    """
    if not isinstance(psmiles, str) or not psmiles.strip():
        raise PSMILESError("PSMILES must be a non-empty string")
    mol = Chem.MolFromSmiles(psmiles, sanitize=False)
    if mol is None:
        raise PSMILESError(f"Chemically invalid PSMILES: {psmiles!r}")
    try:
        Chem.SanitizeMol(mol)
    except Exception as error:
        raise PSMILESError(f"Chemically invalid PSMILES: {psmiles!r}") from error

    dummies = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(dummies) != 2:
        raise InvalidRepeatUnitError(
            "A linear repeat unit requires exactly two attachment dummy atoms; "
            f"found {len(dummies)}"
        )
    for dummy in dummies:
        if dummy.GetDegree() != 1:
            raise InvalidRepeatUnitError(
                f"Attachment dummy atom {dummy.GetIdx()} must have exactly one neighbor"
            )
        if dummy.GetNeighbors()[0].GetAtomicNum() == 0:
            raise InvalidRepeatUnitError(
                "Attachment dummy atoms must connect to real repeat-unit atoms"
            )

    ordered = sorted(dummies, key=lambda atom: atom.GetIdx())
    labels = [atom.GetAtomMapNum() for atom in ordered]
    if labels == [0, 0]:
        head_dummy, tail_dummy = ordered
    elif set(labels) == {1, 2}:
        head_dummy = next(atom for atom in ordered if atom.GetAtomMapNum() == 1)
        tail_dummy = next(atom for atom in ordered if atom.GetAtomMapNum() == 2)
    else:
        raise InvalidRepeatUnitError(
            "Attachment labels are ambiguous: use no labels, or exactly [*:1] "
            "for head and [*:2] for tail"
        )

    return RepeatUnit(
        original_psmiles=psmiles,
        mol=mol,
        head=_attachment_point(mol, head_dummy, "head"),
        tail=_attachment_point(mol, tail_dummy, "tail"),
    )


def _attachment_point(
    mol: Chem.Mol, dummy: Chem.Atom, role: Literal["head", "tail"]
) -> AttachmentPoint:
    neighbor = dummy.GetNeighbors()[0]
    bond = mol.GetBondBetweenAtoms(dummy.GetIdx(), neighbor.GetIdx())
    if bond is None:
        raise InvalidRepeatUnitError("Attachment dummy has no attachment bond")
    return AttachmentPoint(
        dummy_atom_index=dummy.GetIdx(),
        neighbor_atom_index=neighbor.GetIdx(),
        bond_type=bond.GetBondType(),
        bond_order=float(bond.GetBondTypeAsDouble()),
        role=role,
        atom_map_number=dummy.GetAtomMapNum(),
    )
