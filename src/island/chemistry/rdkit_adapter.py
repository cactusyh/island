"""Conversions between RDKit molecules and ISLAND's core data model."""

from collections.abc import Iterable
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem

from island.chemistry._rdkit_stereo import enforce_stored_cip_labels
from island.core import AtomSite, BeadSite, Coordinates, MolecularSystem, Topology
from island.exceptions import (
    EmbeddingError,
    MissingConformerError,
    RDKitConversionError,
    UnsupportedRepresentationError,
)


@dataclass(frozen=True)
class RDKitToSystemResult:
    system: MolecularSystem
    rdkit_index_to_site_id: dict[int, int]
    site_id_to_rdkit_index: dict[int, int]


@dataclass(frozen=True)
class SystemToRDKitResult:
    mol: Chem.Mol
    site_id_to_rdkit_index: dict[int, int]
    rdkit_index_to_site_id: dict[int, int]


def from_rdkit(
    mol: Chem.Mol, *, site_ids: Iterable[int] | None = None, conformer_id: int = -1
) -> RDKitToSystemResult:
    """Convert a conformer-bearing RDKit molecule to an atomistic system."""
    if mol is None:
        raise RDKitConversionError("Cannot convert a null RDKit molecule")
    if mol.GetNumConformers() == 0:
        raise MissingConformerError(
            "RDKit molecule has no conformer; generate coordinates explicitly"
        )
    ids = list(range(1, mol.GetNumAtoms() + 1)) if site_ids is None else list(site_ids)
    _validate_site_ids(ids, mol.GetNumAtoms())
    rdkit_to_site = dict(enumerate(ids))
    site_to_rdkit = {site_id: index for index, site_id in rdkit_to_site.items()}
    try:
        conformer = mol.GetConformer(conformer_id)
    except ValueError as error:
        raise MissingConformerError(
            f"RDKit conformer {conformer_id} does not exist"
        ) from error

    topology = Topology()
    coordinates = Coordinates()
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        site_id = rdkit_to_site[index]
        topology.add_site(
            AtomSite(
                id=site_id,
                name=f"{atom.GetSymbol()}{index + 1}",
                mass=float(atom.GetMass()),
                element=atom.GetSymbol(),
                atomic_number=atom.GetAtomicNum(),
                formal_charge=atom.GetFormalCharge(),
                metadata={
                    "aromatic": atom.GetIsAromatic(),
                    "hybridization": str(atom.GetHybridization()),
                    "chiral_tag": str(atom.GetChiralTag()),
                    "cip_label": (
                        atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
                    ),
                },
            )
        )
        point = conformer.GetAtomPosition(index)
        coordinates.set(site_id, [point.x, point.y, point.z])
    for bond in mol.GetBonds():
        topology.add_bond(
            rdkit_to_site[bond.GetBeginAtomIdx()],
            rdkit_to_site[bond.GetEndAtomIdx()],
            order=float(bond.GetBondTypeAsDouble()),
            aromatic=bond.GetIsAromatic(),
        )
    topology.rebuild_derived_interactions()
    system = MolecularSystem(topology, coordinates, representation="atomistic")
    system.validate()
    return RDKitToSystemResult(system, rdkit_to_site, site_to_rdkit)


def to_rdkit(system: MolecularSystem) -> SystemToRDKitResult:
    """Convert an atomistic ISLAND system to an RDKit molecule."""
    if system.representation != "atomistic":
        raise UnsupportedRepresentationError(
            f"RDKit conversion requires an atomistic system, got {system.representation!r}"
        )
    if any(
        isinstance(site, BeadSite) or not isinstance(site, AtomSite)
        for site in system.topology.sites.values()
    ):
        raise UnsupportedRepresentationError(
            "RDKit conversion supports AtomSite objects only"
        )
    system.validate()
    editable = Chem.RWMol()
    site_to_rdkit: dict[int, int] = {}
    for site_id, site in system.topology.sites.items():
        assert isinstance(site, AtomSite)
        site_to_rdkit[site_id] = editable.AddAtom(_atom_from_site(site))
    for bond in system.topology.bonds.values():
        editable.AddBond(
            site_to_rdkit[bond.site1],
            site_to_rdkit[bond.site2],
            _rdkit_bond_type(bond.order, aromatic=bond.aromatic),
        )
        if bond.aromatic:
            rd_bond = editable.GetBondBetweenAtoms(
                site_to_rdkit[bond.site1], site_to_rdkit[bond.site2]
            )
            rd_bond.SetIsAromatic(True)
            editable.GetAtomWithIdx(site_to_rdkit[bond.site1]).SetIsAromatic(True)
            editable.GetAtomWithIdx(site_to_rdkit[bond.site2]).SetIsAromatic(True)
    converted = editable.GetMol()
    try:
        Chem.SanitizeMol(converted)
    except Exception as error:
        raise RDKitConversionError(
            "RDKit could not sanitize the converted graph"
        ) from error
    enforce_stored_cip_labels(converted, system, site_to_rdkit)
    conformer = Chem.Conformer(converted.GetNumAtoms())
    conformer.Set3D(True)
    for site_id, rdkit_index in site_to_rdkit.items():
        x, y, z = system.coordinates.get(site_id)
        conformer.SetAtomPosition(rdkit_index, (float(x), float(y), float(z)))
    converted.AddConformer(conformer, assignId=True)
    rdkit_to_site = {index: site_id for site_id, index in site_to_rdkit.items()}
    return SystemToRDKitResult(converted, site_to_rdkit, rdkit_to_site)


def from_smiles(
    smiles: str,
    *,
    add_hydrogens: bool = True,
    generate_3d: bool = True,
    random_seed: int = 0xF00D,
    site_ids: Iterable[int] | None = None,
) -> MolecularSystem:
    """Create an atomistic system from SMILES with deterministic embedding."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RDKitConversionError(f"Invalid SMILES: {smiles!r}")
    if add_hydrogens:
        mol = Chem.AddHs(mol)
    if not generate_3d:
        raise MissingConformerError(
            "from_smiles requires generate_3d=True because coordinates cannot be omitted"
        )
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise EmbeddingError(f"RDKit failed to embed SMILES: {smiles!r}")
    return from_rdkit(mol, site_ids=site_ids).system


def _validate_site_ids(site_ids: list[int], atom_count: int) -> None:
    if len(site_ids) != atom_count:
        raise RDKitConversionError(
            f"Expected {atom_count} site IDs, received {len(site_ids)}"
        )
    if any(
        not isinstance(site_id, int) or isinstance(site_id, bool)
        for site_id in site_ids
    ):
        raise RDKitConversionError("All site IDs must be integers")
    if len(set(site_ids)) != len(site_ids):
        raise RDKitConversionError("Site IDs must be unique")


def _atom_from_site(site: AtomSite) -> Chem.Atom:
    try:
        atom = Chem.Atom(
            site.atomic_number if site.atomic_number is not None else site.element
        )
    except Exception as error:
        raise RDKitConversionError(
            f"Cannot create an RDKit atom from site {site.id}"
        ) from error
    atom.SetFormalCharge(site.formal_charge)
    chiral_tag = site.metadata.get("chiral_tag")
    chiral_tags = {
        "CHI_UNSPECIFIED": Chem.ChiralType.CHI_UNSPECIFIED,
        "CHI_TETRAHEDRAL_CW": Chem.ChiralType.CHI_TETRAHEDRAL_CW,
        "CHI_TETRAHEDRAL_CCW": Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
    }
    if chiral_tag in chiral_tags:
        atom.SetChiralTag(chiral_tags[chiral_tag])
    elif chiral_tag is not None:
        raise RDKitConversionError(
            f"Site {site.id} has unsupported chiral tag {chiral_tag!r}"
        )
    atom.SetIsAromatic(bool(site.metadata.get("aromatic", False)))
    return atom


def _rdkit_bond_type(order: float | None, *, aromatic: bool) -> Chem.BondType:
    if aromatic:
        return Chem.BondType.AROMATIC
    mapping = {
        1.0: Chem.BondType.SINGLE,
        2.0: Chem.BondType.DOUBLE,
        3.0: Chem.BondType.TRIPLE,
    }
    try:
        return mapping[float(order)]
    except (KeyError, TypeError, ValueError) as error:
        raise RDKitConversionError(
            f"Unsupported chemical bond order: {order!r}"
        ) from error
