"""Map completed ISLAND chains to the native RadonPy atom traversal."""

from rdkit import Chem

from island.chemistry import parse_psmiles
from island.chemistry._stereo_transfer import transfer_stereo as _transfer


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
