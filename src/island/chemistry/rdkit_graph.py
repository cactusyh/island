"""Graph-only conversion between ISLAND chemistry and RDKit."""

from dataclasses import dataclass

from rdkit import Chem

from island.core import AtomSite, BeadSite, MolecularSystem, Topology
from island.exceptions import RDKitConversionError, UnsupportedRepresentationError


@dataclass(frozen=True)
class GraphToRDKitResult:
    """A coordinate-free RDKit graph and explicit stable-ID mappings."""

    mol: Chem.Mol
    site_id_to_rdkit_index: dict[int, int]
    rdkit_index_to_site_id: dict[int, int]


def topology_to_rdkit_graph(
    topology: Topology, *, representation: str = "atomistic"
) -> GraphToRDKitResult:
    """Convert chemical topology without requiring or generating coordinates."""
    if representation != "atomistic":
        raise UnsupportedRepresentationError(
            f"RDKit graph conversion requires atomistic representation, got "
            f"{representation!r}"
        )
    topology.validate_bond_graph()
    if any(
        isinstance(site, BeadSite) or not isinstance(site, AtomSite)
        for site in topology.sites.values()
    ):
        raise UnsupportedRepresentationError(
            "RDKit graph conversion supports AtomSite objects only"
        )
    editable = Chem.RWMol()
    site_to_rdkit: dict[int, int] = {}
    for site_id in sorted(topology.sites):
        site = topology.sites[site_id]
        assert isinstance(site, AtomSite)
        site_to_rdkit[site_id] = editable.AddAtom(atom_site_to_rdkit(site))
    for bond_key in sorted(topology.bonds):
        bond = topology.bonds[bond_key]
        editable.AddBond(
            site_to_rdkit[bond.site1],
            site_to_rdkit[bond.site2],
            rdkit_bond_type(bond.order, aromatic=bond.aromatic),
        )
        if bond.aromatic:
            rdkit_bond = editable.GetBondBetweenAtoms(
                site_to_rdkit[bond.site1], site_to_rdkit[bond.site2]
            )
            rdkit_bond.SetIsAromatic(True)
            editable.GetAtomWithIdx(site_to_rdkit[bond.site1]).SetIsAromatic(True)
            editable.GetAtomWithIdx(site_to_rdkit[bond.site2]).SetIsAromatic(True)
    converted = editable.GetMol()
    try:
        Chem.SanitizeMol(converted)
    except Exception as error:
        raise RDKitConversionError(
            "RDKit could not sanitize the converted molecular graph"
        ) from error
    rdkit_to_site = {index: site_id for site_id, index in site_to_rdkit.items()}
    return GraphToRDKitResult(converted, site_to_rdkit, rdkit_to_site)


def system_to_rdkit_graph(system: MolecularSystem) -> GraphToRDKitResult:
    """Convert only a system's authoritative graph, ignoring coordinates."""
    return topology_to_rdkit_graph(
        system.topology, representation=system.representation
    )


def atom_site_to_rdkit(site: AtomSite) -> Chem.Atom:
    """Create an RDKit atom from supported chemical identity fields."""
    try:
        atom = Chem.Atom(
            site.atomic_number if site.atomic_number is not None else site.element
        )
    except Exception as error:
        raise RDKitConversionError(
            f"Cannot create an RDKit atom from site {site.id}"
        ) from error
    atom.SetFormalCharge(site.formal_charge)
    atom.SetNoImplicit(bool(site.metadata.get("no_implicit_hydrogens", False)))
    explicit_hydrogens = site.metadata.get("explicit_hydrogen_count", 0)
    if not isinstance(explicit_hydrogens, int) or explicit_hydrogens < 0:
        raise RDKitConversionError(
            f"Site {site.id} has invalid explicit hydrogen count {explicit_hydrogens!r}"
        )
    atom.SetNumExplicitHs(explicit_hydrogens)
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


def rdkit_bond_type(order: float | None, *, aromatic: bool) -> Chem.BondType:
    """Map supported chemical bond orders to RDKit bond types."""
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
