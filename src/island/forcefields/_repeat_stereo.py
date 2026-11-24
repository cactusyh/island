"""Transfer repeat-local stereo to the finished chain using the builder's atom map."""

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

from island.chemistry import parse_psmiles, to_rdkit
from island.exceptions import EmbeddingError


def _transfer(mol, fragment, mapping):
    for bond in fragment.GetBonds():
        controls = bond.GetStereoAtoms()
        if len(controls) != 2:
            continue
        a, b = mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()]
        dest = mol.GetBondBetweenAtoms(a, b)
        left, right = (mapping[idx] for idx in controls)
        if dest.GetBeginAtomIdx() != a:
            left, right = right, left
        dest.SetStereoAtoms(left, right)
        dest.SetStereo(bond.GetStereo())

    for atom in fragment.GetAtoms():
        if atom.GetChiralTag() not in (
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
        ):
            continue
        dest = mol.GetAtomWithIdx(mapping[atom.GetIdx()])
        source_order = [mapping[v.GetIdx()] for v in atom.GetNeighbors()]
        actual = [v.GetIdx() for v in dest.GetNeighbors()]
        # RDKit puts an implicit stereocentre H last in its parity convention.
        source_order.extend(v for v in actual if v not in source_order)
        parity = [source_order.index(v) for v in actual]
        dest.SetChiralTag(atom.GetChiralTag())
        if sum(i > j for k, i in enumerate(parity) for j in parity[k + 1 :]) % 2:
            dest.InvertChirality()


def stereo_chain(system, psmiles, dp, seed, embed, end_groups):
    graph = to_rdkit(system)
    mol = graph.mol
    unit = parse_psmiles(psmiles)
    # The two placeholders can be equivalent in isolation and distinct in
    # the finished chain. Only clean redundant stereo after substitution.
    Chem.SetBondStereoFromDirections(unit.mol)
    maps = [{} for _ in range(dp)]
    caps = {}
    cap_maps = {"head": {}, "tail": {}}
    for i, site in system.topology.sites.items():
        meta = site.metadata
        index = graph.site_id_to_rdkit_index[i]
        if "source_repeat_atom_index" in meta:
            maps[meta["repeat_unit_index"]][meta["source_repeat_atom_index"]] = index
        if "end_group_role" in meta:
            role = meta["end_group_role"]
            caps[index] = role
            if "source_end_group_atom_index" in meta:
                cap_maps[role][meta["source_end_group_atom_index"]] = index

    # A dummy's replacement can be in the neighbouring repeat or in its cap.
    for n, mapping in enumerate(maps):
        for point, other, step in (
            (unit.head, unit.tail, -1),
            (unit.tail, unit.head, 1),
        ):
            if 0 <= n + step < dp:
                idx = maps[n + step][other.neighbor_atom_index]
            else:
                parent = mol.GetAtomWithIdx(mapping[point.neighbor_atom_index])
                idx = next(
                    a.GetIdx()
                    for a in parent.GetNeighbors()
                    if caps.get(a.GetIdx()) == point.role
                )
            mapping[point.dummy_atom_index] = idx

        _transfer(mol, unit.mol, mapping)

    for role, cap_smiles in end_groups.items():
        if not any(c in cap_smiles for c in "/\\@"):
            continue
        fragment = Chem.MolFromSmiles(cap_smiles)
        mapping = cap_maps[role]
        dummy = next(a for a in fragment.GetAtoms() if a.GetAtomicNum() == 0)
        parent_key = role + "_site_id"
        mapping[dummy.GetIdx()] = graph.site_id_to_rdkit_index[
            system.metadata["polymer"][parent_key]
        ]
        _transfer(mol, fragment, mapping)

    mol.RemoveAllConformers()
    Chem.SetDoubleBondNeighborDirections(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    if embed:
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise EmbeddingError("RDKit could not embed the stereospecified chain")
    else:
        rdDepictor.Compute2DCoords(mol)
    return mol, tuple(graph.rdkit_index_to_site_id[i] for i in range(mol.GetNumAtoms()))


def radonpy_order(mol, sites, psmiles, dp):
    """Use RadonPy's atom and bond traversal, carrying over ISLAND stereo and coordinates."""
    import json

    from radonpy.core.poly import polymerize_MolFromSmiles

    from island.exceptions import ForceFieldError

    if "[*:" in psmiles:
        unit = parse_psmiles(psmiles)
        for atom in unit.mol.GetAtoms():
            atom.SetAtomMapNum(0)
        psmiles = Chem.MolToSmiles(
            unit.mol, rootedAtAtom=unit.head.dummy_atom_index, canonical=False
        )
    if dp == 1:
        # RadonPy's linear builder needs at least one inter-repeat bond.
        # Its single-molecule entry point reads a capped SMILES and adds H.
        capped = Chem.MolToSmiles(Chem.RemoveHs(mol))
        native = Chem.AddHs(Chem.MolFromSmiles(capped))
    else:
        native = polymerize_MolFromSmiles(psmiles, n=dp, terminal="[H]")
    if native is None:
        raise ForceFieldError(
            "RadonPy could not construct this H-capped repeat in native order"
        )

    def traversal(graph):
        stripped = Chem.Mol(graph)
        Chem.RemoveStereochemistry(stripped)
        for atom in stripped.GetAtoms():
            atom.SetAtomMapNum(0)
        smiles = Chem.MolToSmiles(stripped)
        return smiles, json.loads(stripped.GetProp("_smilesAtomOutputOrder"))

    source_smiles, source = traversal(mol)
    native_smiles, dest = traversal(native)
    if native_smiles != source_smiles:
        raise ForceFieldError(
            "RadonPy and ISLAND produced different chain connectivity"
        )
    mapping = dict(zip(source, dest, strict=True))
    Chem.RemoveStereochemistry(native)
    for atom in native.GetAtoms():
        atom.SetAtomMapNum(0)
    _transfer(native, mol, mapping)
    Chem.SetDoubleBondNeighborDirections(native)
    Chem.AssignStereochemistry(native, force=True, cleanIt=True)
    conf = Chem.Conformer(native.GetNumAtoms())
    original = mol.GetConformer()
    conf.Set3D(original.Is3D())
    ids = [0] * len(sites)
    for i, j in mapping.items():
        conf.SetAtomPosition(j, original.GetAtomPosition(i))
        ids[j] = sites[i]
    native.AddConformer(conf)
    return native, tuple(ids)
