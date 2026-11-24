"""Explicit-H DGEBA/PACM chemistry. RDKit indices are stable site identities."""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem

DGEBA = "CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1"
PACM = "NC1CCC(CC1)CC1CCC(N)CC1"
EPOXIDE = Chem.MolFromSmarts("[CH2;r3]1[CH;r3][O;r3]1")


def monomer(smiles, seed):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise ValueError("Monomer embedding failed")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    return mol


def epoxides(mol):
    Chem.GetSymmSSSR(mol)
    return {c: o for c, _, o in mol.GetSubstructMatches(EPOXIDE)}


def amines(mol):
    sites = {}
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 7:
            sites[a.GetIdx()] = tuple(
                sorted(h.GetIdx() for h in a.GetNeighbors() if h.GetAtomicNum() == 1)
            )
    return sites


@dataclass(frozen=True)
class Addition:
    carbon: int
    oxygen: int
    nitrogen: int
    hydrogen: int


def react(mol, events):
    rings, nitrogens = epoxides(mol), amines(mol)
    edit = Chem.RWMol(mol)
    used = set()
    for e in events:
        touched = {e.carbon, e.oxygen, e.nitrogen, e.hydrogen}
        if (
            len(touched) != 4
            or used & touched
            or rings.get(e.carbon) != e.oxygen
            or e.hydrogen not in nitrogens.get(e.nitrogen, ())
            or mol.GetBondBetweenAtoms(e.carbon, e.nitrogen) is not None
        ):
            raise ValueError("Ineligible or overlapping epoxy–amine addition")
        used.update(touched)
        edit.RemoveBond(e.carbon, e.oxygen)
        edit.RemoveBond(e.nitrogen, e.hydrogen)
        edit.AddBond(e.carbon, e.nitrogen, Chem.BondType.SINGLE)
        edit.AddBond(e.oxygen, e.hydrogen, Chem.BondType.SINGLE)
    product = edit.GetMol()
    Chem.SanitizeMol(product)
    return product


def reference_reactions(seed=2026):
    d, p = monomer(DGEBA, seed), monomer(PACM, seed)
    n = min(amines(p))
    h1, h2 = amines(p)[n]
    c, o = next(iter(epoxides(d).items()))
    pre1 = Chem.CombineMols(p, d)
    offset = p.GetNumAtoms()
    first = Addition(offset + c, offset + o, n, h1)
    post1 = react(pre1, [first])
    pre2 = Chem.CombineMols(post1, d)
    offset = post1.GetNumAtoms()
    second = Addition(offset + c, offset + o, n, h2)
    return d, p, [(pre1, first, post1), (pre2, second, react(pre2, [second]))]
