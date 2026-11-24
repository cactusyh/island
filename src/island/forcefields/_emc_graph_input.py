"""Convert finite molecular graphs to explicit sites, coordinates and EMC groups."""

import json
from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from island.chemistry import system_to_rdkit_graph, to_rdkit
from island.chemistry._rdkit_stereo import enforce_stored_cip_labels
from island.core import AtomSite, Coordinates, MolecularSystem, SimulationBox, Topology
from island.core.coordinate_provenance import previous_coordinate_source
from island.exceptions import ForceFieldError
from island.forcefields.assignment_core.schema import ChemicalGraph
from island.forcefields.assignment_core.topology import molecule
from island.graph.final import FinalChemicalGraph


@dataclass
class GraphInput:
    mol: Chem.Mol
    sites: list[int]
    generated: tuple[int, ...]
    xyz: np.ndarray
    bounds: np.ndarray
    coordinate_source: str


def _from_final(graph):
    p = graph.payload
    top = Topology()
    for row in p["sites"]:
        top.add_site(
            AtomSite(
                id=row["id"],
                name=row["element"],
                element=row["element"],
                atomic_number=row["atomic_number"],
                mass=row["mass"],
                formal_charge=row["formal_charge"],
                metadata=row["metadata"],
            )
        )
    for bond in p["bonds"]:
        top.add_bond(*bond["sites"], order=bond["order"], aromatic=bond["aromatic"])
    cell = p["periodic"]
    box = (
        SimulationBox(*cell["lengths"], periodic=tuple(cell["periodic"]))
        if cell
        else None
    )
    return MolecularSystem(top, Coordinates(), box=box)


def _unwrap(mol, xyz, lengths):
    """Reconstruct a finite molecule from orthorhombic minimum-image bond vectors."""
    out = xyz.copy()
    for fragment in Chem.GetMolFrags(mol):
        visited = {fragment[0]}
        stack = [fragment[0]]
        while stack:
            i = stack.pop()
            for atom in mol.GetAtomWithIdx(i).GetNeighbors():
                j = atom.GetIdx()
                delta = xyz[j] - xyz[i]
                candidate = out[i] + delta - lengths * np.rint(delta / lengths)
                if j not in visited:
                    visited.add(j)
                    stack.append(j)
                    out[j] = candidate
                elif not np.allclose(out[j], candidate, rtol=0, atol=1e-6):
                    raise ForceFieldError(
                        "A bonded cycle winds through the box; a finite graph is required"
                    )
    return out


def _embed(mol, seed):
    indices = []
    fragments = Chem.GetMolFrags(mol, asMols=True, fragsMolAtomMapping=indices)
    xyz = np.empty((mol.GetNumAtoms(), 3))
    end = 0.0
    for part, ids in zip(fragments, indices, strict=True):
        options = AllChem.ETKDGv3()
        options.randomSeed = seed
        if AllChem.EmbedMolecule(part, options) != 0:
            options.useRandomCoords = True
            if AllChem.EmbedMolecule(part, options) != 0:
                raise ForceFieldError("ETKDG failed; supply a 3D conformer")
        coords = part.GetConformer().GetPositions()
        coords[:, 0] += end - coords[:, 0].min()
        end = coords[:, 0].max() + 8.0
        xyz[list(ids)] = coords
    return xyz


def _check_stereo(mol, xyz):
    specified = Chem.StereoSpecified.Specified
    expected = {
        (s.type, s.centeredOn): s.descriptor
        for s in Chem.FindPotentialStereo(mol)
        if s.specified == specified
    }
    if not expected:
        return
    check = Chem.Mol(mol)
    check.RemoveAllConformers()
    conf = Chem.Conformer(len(xyz))
    for i, pos in enumerate(xyz):
        conf.SetAtomPosition(i, pos)
    check.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(check)
    actual = {
        (s.type, s.centeredOn): s.descriptor
        for s in Chem.FindPotentialStereo(check)
        if s.specified == specified
    }
    if any(actual.get(k) != v for k, v in expected.items()):
        raise ForceFieldError("3D coordinates disagree with declared stereochemistry")


def graph_input(graph, seed):
    box = None
    if isinstance(graph, FinalChemicalGraph):
        graph = _from_final(graph)
    if isinstance(graph, MolecularSystem):
        box = graph.box
        if len(graph.coordinates):
            converted = to_rdkit(graph)
            if previous_coordinate_source(graph.metadata) == "rdkit_2d":
                converted.mol.RemoveAllConformers()
        else:
            converted = system_to_rdkit_graph(graph)
            enforce_stored_cip_labels(
                converted.mol, graph, converted.site_id_to_rdkit_index
            )
        mol = converted.mol
        sites = [converted.rdkit_index_to_site_id[i] for i in range(mol.GetNumAtoms())]
    elif isinstance(graph, ChemicalGraph):
        if graph.blockers:
            raise ForceFieldError(f"Unresolved input graph: {graph.blockers}")
        mol = molecule(graph)
        sites = list(graph.atom_order)
    elif isinstance(graph, Chem.Mol):
        mol = Chem.Mol(graph)
        sites = list(range(mol.GetNumAtoms()))
    else:
        raise TypeError(
            "Expected a MolecularSystem, ChemicalGraph, FinalChemicalGraph or RDKit Mol"
        )

    Chem.SanitizeMol(mol)
    if not sites or any(
        a.GetAtomicNum() == 0 or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        raise ForceFieldError("Supply a terminated, closed-shell molecular graph")
    has_3d = mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()
    if box:
        lengths = np.asarray(box.lengths)
        if not all(box.periodic) or not np.isfinite(lengths).all() or min(lengths) <= 0:
            raise ForceFieldError(
                "Graph MD requires a finite, fully periodic orthorhombic cell"
            )
        if has_3d:
            # Place any missing H using local, contiguous bond vectors.
            conf = mol.GetConformer()
            for i, pos in enumerate(_unwrap(mol, conf.GetPositions(), lengths)):
                conf.SetAtomPosition(i, pos)
        elif len(Chem.GetMolFrags(mol)) > 1:
            raise ForceFieldError(
                "Supply packed 3D positions for multiple molecules in a specified box"
            )
    mol = Chem.AddHs(mol, addCoords=has_3d)
    generated = tuple(
        range(max(sites) + 1, max(sites) + 1 + mol.GetNumAtoms() - len(sites))
    )
    sites.extend(generated)
    xyz = mol.GetConformer().GetPositions() if has_3d else _embed(mol, seed)
    if not np.isfinite(xyz).all():
        raise ForceFieldError("Non-finite starting coordinates")
    if box:
        bounds = np.column_stack((np.zeros(3), lengths))
    else:
        bounds = np.column_stack((xyz.min(axis=0) - 12, xyz.max(axis=0) + 12))
    _check_stereo(mol, xyz)
    return GraphInput(
        mol, sites, generated, xyz, bounds, "input_3d" if has_3d else "etkdg"
    )


def emc_groups(mol, ua, kekule):
    """Retain the writer's actual SMILES traversal, including every H in all-atom mode."""
    indices = []
    parts = Chem.GetMolFrags(mol, asMols=True, fragsMolAtomMapping=indices)
    groups = {}
    orders = {}
    for n, (part, ids) in enumerate(zip(parts, indices, strict=True)):
        # PSF segment names are four characters. Different input components
        # must remain identifiable even when EMC sorts its clusters.
        name = "g" + np.base_repr(n, base=36).lower()
        if len(name) > 4:
            raise ForceFieldError(
                "Too many distinct EMC group labels for the PSF segment column"
            )
        for i, atom in enumerate(part.GetAtoms()):
            atom.SetIntProp("_island_index", ids[i])
            atom.SetAtomMapNum(0)
            if ua and atom.GetIsotope():
                raise ForceFieldError(
                    "United-atom isotope mass assignment requires an all-atom field"
                )
            atom.SetIsotope(0)
        if ua:
            part = Chem.RemoveHs(part)
        # EMC cannot read slash bonds. Geometry and declared stereo are checked
        # on the original molecule, not on EMC's disposable packing coordinates.
        for bond in part.GetBonds():
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            if bond.GetBondDir() in (
                Chem.BondDir.ENDUPRIGHT,
                Chem.BondDir.ENDDOWNRIGHT,
            ):
                bond.SetBondDir(Chem.BondDir.NONE)
        groups[name] = Chem.MolToSmiles(part, canonical=False, kekuleSmiles=kekule)
        orders[name] = [
            part.GetAtomWithIdx(i).GetIntProp("_island_index")
            for i in json.loads(part.GetProp("_smilesAtomOutputOrder"))
        ]
    return groups, orders
