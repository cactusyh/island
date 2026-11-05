"""Reuse the PSMILES builder; perform every graph edit before explicit-H typing."""

import json
from dataclasses import replace

from island.forcefields.assignment_core.schema import (
    ChemicalGraph,
    PolymerSpec,
    require,
)


def molecule(graph):
    from rdkit import Chem

    mols = Chem.JSONToMols(graph.molecule_json)
    require(len(mols) == 1, "one_final_graph_required")
    mol = mols[0]
    require(
        graph.atom_order == tuple(range(1, mol.GetNumAtoms() + 1)),
        "atom_order_or_mapping_mutation",
    )
    require(
        all(
            a.GetAtomicNum() > 0
            and not a.GetNumImplicitHs()
            and not a.GetNumExplicitHs()
            for a in mol.GetAtoms()
        ),
        "final_graph_requires_explicit_H_no_dummies",
    )
    return mol


def as_topology(graph):
    """Adapt into the existing coordinate-free typing/inventory interfaces."""
    from island.core import AtomSite, Topology

    mol = molecule(graph)
    topology = Topology()
    for a in mol.GetAtoms():
        site = a.GetIdx() + 1
        topology.add_site(
            AtomSite(
                id=site,
                name=f"{a.GetSymbol()}{site}",
                mass=float(a.GetMass()),
                element=a.GetSymbol(),
                atomic_number=a.GetAtomicNum(),
                formal_charge=a.GetFormalCharge(),
                metadata={
                    "aromatic": a.GetIsAromatic(),
                    "isotope": a.GetIsotope(),
                    "radical_electrons": a.GetNumRadicalElectrons(),
                    "no_implicit_hydrogens": True,
                    "chiral_tag": str(a.GetChiralTag()),
                },
            )
        )
    for b in mol.GetBonds():
        topology.add_bond(
            b.GetBeginAtomIdx() + 1,
            b.GetEndAtomIdx() + 1,
            order=b.GetBondTypeAsDouble(),
            aromatic=b.GetIsAromatic(),
        )
    return topology


def _consume_h(mol, index, count=1):
    a = mol.GetAtomWithIdx(index)
    require(a.GetTotalNumHs() >= count, "attachment_requires_replaceable_hydrogen")
    # Bracket-defined [C@H] carries an explicit H count, not an H site yet.
    explicit = min(count, a.GetNumExplicitHs())
    if explicit:
        a.SetNumExplicitHs(a.GetNumExplicitHs() - explicit)


def normalize(spec: PolymerSpec):
    """Deterministic graph construction only; never embed, evaluate or optimize."""
    import rdkit
    from rdkit import Chem

    from island.builders.polymer import build_polymer_from_sequence
    from island.chemistry.psmiles import parse_psmiles
    from island.chemistry.rdkit_adapter import to_rdkit

    library = dict(spec.repeat_units)
    for label, head, tail in spec.connection_points:
        require(label in library, "undefined_connection_repeat")
        repeat = parse_psmiles(library[label])
        require(
            (head, tail)
            == (repeat.head.dummy_atom_index, repeat.tail.dummy_atom_index),
            "connection_points_disagree_with_PSMILES_orientation",
        )
    system = build_polymer_from_sequence(
        library,
        spec.sequence or (spec.repeat_units[0][0],) * spec.degree_of_polymerization,
        generate_3d=False,
        add_hydrogens=False,
        tacticity=spec.tacticity,
        stereo_seed=spec.stereo_seed,
        head_end_group=spec.end_groups[0],
        tail_end_group=spec.end_groups[1],
    )
    mol = to_rdkit(system).mol
    mol.RemoveAllConformers()  # Builder's 2D depiction is never a physical geometry.
    heavy_count = mol.GetNumAtoms()
    for branch in sorted(spec.branches, key=lambda b: (b.anchor_site, b.psmiles)):
        require(
            type(branch.anchor_site) is int and 1 <= branch.anchor_site <= heavy_count,
            "branch_anchor_not_in_pre_edit_graph",
        )
        piece = Chem.MolFromSmiles(branch.psmiles)
        require(piece is not None, "invalid_branch_SMILES")
        dummies = [a for a in piece.GetAtoms() if not a.GetAtomicNum()]
        require(
            len(dummies) == 1 and dummies[0].GetDegree() == 1,
            "branch_requires_one_attachment_dummy",
        )
        dummy = dummies[0]
        neighbor = dummy.GetNeighbors()[0]
        require(
            neighbor.GetAtomicNum() > 0
            and piece.GetBondBetweenAtoms(
                dummy.GetIdx(), neighbor.GetIdx()
            ).GetBondType()
            == Chem.BondType.SINGLE,
            "branch_requires_single_attachment_bond",
        )
        _consume_h(mol, branch.anchor_site - 1)
        offset = mol.GetNumAtoms()
        edit = Chem.RWMol(Chem.CombineMols(mol, piece))
        edit.AddBond(
            branch.anchor_site - 1, offset + neighbor.GetIdx(), Chem.BondType.SINGLE
        )
        edit.RemoveAtom(offset + dummy.GetIdx())
        mol = edit.GetMol()
        Chem.SanitizeMol(mol)
    for cross in sorted(
        spec.crosslinks, key=lambda c: (min(c.site1, c.site2), max(c.site1, c.site2))
    ):
        require(
            all(
                type(i) is int and 1 <= i <= heavy_count
                for i in (cross.site1, cross.site2)
            ),
            "crosslink_sites_not_in_pre_edit_graph",
        )
        require(
            cross.site1 != cross.site2 and cross.order == 1,
            "crosslink_requires_distinct_sites_and_single_bond",
        )
        a, b = cross.site1 - 1, cross.site2 - 1
        require(mol.GetBondBetweenAtoms(a, b) is None, "crosslink_bond_already_exists")
        _consume_h(mol, a)
        _consume_h(mol, b)
        edit = Chem.RWMol(mol)
        edit.AddBond(a, b, Chem.BondType.SINGLE)
        mol = edit.GetMol()
        Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    blockers = (
        ("periodic_topology_requires_explicit_image_graph_adapter",)
        if spec.periodic_topology
        else ()
    )
    architecture = (
        "finite_crosslinked"
        if spec.crosslinks
        else "finite_branched"
        if spec.branches
        else "finite_linear"
    )
    encoded = json.dumps(
        json.loads(Chem.MolsToJSON([mol])),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    graph = ChemicalGraph(
        spec,
        encoded,
        tuple(range(1, mol.GetNumAtoms() + 1)),
        architecture,
        blockers,
        rdkit_version=rdkit.__version__,
    )
    molecule(graph)
    return graph


def rebuild(graph, **spec_changes):
    """Graph edits create a new identity. Previously typed snapshots cannot be reused."""
    return normalize(replace(graph.spec, **spec_changes))
