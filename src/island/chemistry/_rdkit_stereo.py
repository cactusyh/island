"""Internal helpers for restoring graph stereochemistry in RDKit adapters."""

from rdkit import Chem

from island.core import MolecularSystem
from island.exceptions import RDKitConversionError


def enforce_stored_cip_labels(
    mol: Chem.Mol,
    system: MolecularSystem,
    site_id_to_rdkit_index: dict[int, int],
) -> None:
    """Set tetrahedral tags so stored R/S labels survive graph reconstruction."""
    targets = {
        site_id_to_rdkit_index[site_id]: site.metadata.get("cip_label")
        for site_id, site in system.topology.sites.items()
        if site.metadata.get("cip_label") in {"R", "S"}
    }
    for atom_index, target in targets.items():
        atom = mol.GetAtomWithIdx(atom_index)
        atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if _cip_label(atom) == target:
            continue
        atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if _cip_label(atom) != target:
            raise RDKitConversionError(
                f"Cannot restore stored {target} stereochemistry for site mapped "
                f"to RDKit atom {atom_index}"
            )
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    incorrect = {
        atom_index: (target, _cip_label(mol.GetAtomWithIdx(atom_index)))
        for atom_index, target in targets.items()
        if _cip_label(mol.GetAtomWithIdx(atom_index)) != target
    }
    if incorrect:
        raise RDKitConversionError(
            f"Stored stereochemistry changed during RDKit conversion: {incorrect}"
        )


def _cip_label(atom: Chem.Atom) -> str | None:
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
