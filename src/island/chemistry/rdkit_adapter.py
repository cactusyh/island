"""Coordinate-aware conversions between RDKit and ISLAND's core data model."""

from collections.abc import Iterable
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem

from island.chemistry._rdkit_stereo import enforce_stored_cip_labels
from island.chemistry.rdkit_graph import system_to_rdkit_graph
from island.core import AtomSite, Coordinates, MolecularSystem, Topology
from island.exceptions import (
    EmbeddingError,
    MissingConformerError,
    RDKitConversionError,
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
                    "no_implicit_hydrogens": atom.GetNoImplicit(),
                    "explicit_hydrogen_count": atom.GetNumExplicitHs(),
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
    """Convert an ISLAND system graph and coordinates to an RDKit molecule."""
    system.coordinates.validate(system.topology)
    graph = system_to_rdkit_graph(system)
    converted = graph.mol
    enforce_stored_cip_labels(converted, system, graph.site_id_to_rdkit_index)
    conformer = Chem.Conformer(converted.GetNumAtoms())
    conformer.Set3D(True)
    for site_id, rdkit_index in graph.site_id_to_rdkit_index.items():
        x, y, z = system.coordinates.get(site_id)
        conformer.SetAtomPosition(rdkit_index, (float(x), float(y), float(z)))
    converted.AddConformer(conformer, assignId=True)
    return SystemToRDKitResult(
        converted,
        graph.site_id_to_rdkit_index,
        graph.rdkit_index_to_site_id,
    )


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
