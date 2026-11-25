"""Transfer stereo through atom substitution, including implicit terminal H."""

from rdkit import Chem


def transfer_bond_stereo(mol, fragment, mapping):
    """Map a fragment onto a completed graph; -1 denotes a terminal implicit H."""
    for bond in fragment.GetBonds():
        controls = bond.GetStereoAtoms()
        if len(controls) != 2:
            continue
        if bond.GetBeginAtomIdx() not in mapping or bond.GetEndAtomIdx() not in mapping:
            continue
        a, b = mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()]
        dest = mol.GetBondBetweenAtoms(a, b)
        left, right = (mapping[idx] for idx in controls)
        stereo = bond.GetStereo()
        for side, center, other in ((0, a, b), (1, b, a)):
            if (left, right)[side] != -1:
                continue
            choices = [
                n.GetIdx()
                for n in mol.GetAtomWithIdx(center).GetNeighbors()
                if n.GetIdx() != other
            ]
            if not choices:  # Terminal CH2 double bond is no longer stereogenic.
                break
            if side == 0:
                left = choices[0]
            else:
                right = choices[0]
            stereo = {
                Chem.BondStereo.STEREOE: Chem.BondStereo.STEREOZ,
                Chem.BondStereo.STEREOZ: Chem.BondStereo.STEREOE,
                Chem.BondStereo.STEREOCIS: Chem.BondStereo.STEREOTRANS,
                Chem.BondStereo.STEREOTRANS: Chem.BondStereo.STEREOCIS,
            }[stereo]
        else:
            if dest.GetBeginAtomIdx() != a:
                left, right = right, left
            dest.SetStereoAtoms(left, right)
            dest.SetStereo(stereo)


def transfer_stereo(mol, fragment, mapping):
    """Preserve bond stereo and tetrahedral neighbor parity through substitution."""
    transfer_bond_stereo(mol, fragment, mapping)
    for atom in fragment.GetAtoms():
        if atom.GetChiralTag() not in (
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
        ):
            continue
        dest = mol.GetAtomWithIdx(mapping[atom.GetIdx()])
        source_order = [mapping[v.GetIdx()] for v in atom.GetNeighbors()]
        actual = [v.GetIdx() for v in dest.GetNeighbors()]
        if -1 in source_order:
            actual.append(-1)
        source_order.extend(v for v in actual if v not in source_order)
        parity = [source_order.index(v) for v in actual]
        dest.SetChiralTag(atom.GetChiralTag())
        if sum(i > j for k, i in enumerate(parity) for j in parity[k + 1 :]) % 2:
            dest.InvertChirality()


def verify_bond_geometry(mol):
    """Infer E/Z independently from the conformer, without reusing graph tags."""
    from island.exceptions import UnsupportedStereochemistryError

    expected = {
        s.centeredOn: s.descriptor
        for s in Chem.FindPotentialStereo(mol)
        if s.type == Chem.StereoType.Bond_Double
        and s.specified == Chem.StereoSpecified.Specified
    }
    if not expected:
        return
    check = Chem.Mol(mol)
    Chem.RemoveStereochemistry(check)
    Chem.AssignStereochemistryFrom3D(check, replaceExistingTags=True)
    observed = {
        s.centeredOn: s.descriptor
        for s in Chem.FindPotentialStereo(check)
        if s.type == Chem.StereoType.Bond_Double
        and s.specified == Chem.StereoSpecified.Specified
    }
    if any(observed.get(i) != state for i, state in expected.items()):
        raise UnsupportedStereochemistryError(
            "3D coordinates disagree with declared double-bond stereochemistry"
        )
