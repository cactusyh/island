"""Mapped epoxy ring opening by amine or thiol donors. RDKit indices are stable site identities."""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem

DGEBA = "CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1"
PACM = "NC1CCC(CC1)CC1CCC(N)CC1"
PETMP = "C(COC(=O)CCS)(COC(=O)CCS)(COC(=O)CCS)COC(=O)CCS"
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


def donor_sites(mol, atomic_number):
    sites = {}
    for a in mol.GetAtoms():
        if (
            a.GetAtomicNum() == atomic_number
            and not a.GetFormalCharge()
            and not a.GetIsAromatic()
            and (
                atomic_number != 7 or a.GetHybridization() == Chem.HybridizationType.SP3
            )
            and (atomic_number != 16 or a.GetTotalValence() == 2)
        ):
            sites[a.GetIdx()] = tuple(
                sorted(h.GetIdx() for h in a.GetNeighbors() if h.GetAtomicNum() == 1)
            )
    return sites


def amines(mol):
    return donor_sites(mol, 7)


@dataclass(frozen=True)
class Addition:
    carbon: int
    oxygen: int
    nitrogen: int
    hydrogen: int

    @property
    def donor(self):
        return self.nitrogen


@dataclass(frozen=True)
class ThiolAddition:
    carbon: int
    oxygen: int
    sulfur: int
    hydrogen: int

    @property
    def donor(self):
        return self.sulfur


@dataclass(frozen=True)
class RingOpening:
    name: str
    hardener: str
    element: int
    sites_per_hardener: int
    max_additions: int
    event_type: type
    reference_names: tuple[str, ...]
    product_donor_angle_deg: float | None = None

    def donors(self, mol):
        return donor_sites(mol, self.element)

    def apply(self, edit, event):
        """Commit the mapped ring opening and donor-H transfer to a writable graph."""
        edit.RemoveBond(event.carbon, event.oxygen)
        edit.RemoveBond(event.donor, event.hydrogen)
        edit.AddBond(event.carbon, event.donor, Chem.BondType.SINGLE)
        edit.AddBond(event.oxygen, event.hydrogen, Chem.BondType.SINGLE)

    def event(self, record):
        return self.event_type(
            **{k: record[k] for k in self.event_type.__dataclass_fields__}
        )


AMINE = RingOpening(
    "pacm", PACM, 7, 2, 2, Addition, ("dgeba", "pacm", "secondary", "tertiary")
)
# Equilibrium angle of the native PCFF c,sc,c thioether angle (degrees).
THIOL = RingOpening(
    "petmp", PETMP, 16, 4, 1, ThiolAddition, ("dgeba", "petmp", "thioether"), 97.5
)
TEMPLATES = {t.name: t for t in (AMINE, THIOL)}


def react(mol, events, template=AMINE):
    rings, donors = epoxides(mol), template.donors(mol)
    edit = Chem.RWMol(mol)
    used = set()
    for e in events:
        touched = {e.carbon, e.oxygen, e.donor, e.hydrogen}
        if (
            len(touched) != 4
            or used & touched
            or rings.get(e.carbon) != e.oxygen
            or not isinstance(e, template.event_type)
            or e.hydrogen not in donors.get(e.donor, ())
            or mol.GetBondBetweenAtoms(e.carbon, e.donor) is not None
        ):
            raise ValueError("Ineligible or overlapping epoxy addition")
        used.update(touched)
        template.apply(edit, e)
    product = edit.GetMol()
    Chem.SanitizeMol(product)
    return product


def reference_reactions(seed=2026, template=AMINE):
    d, p = monomer(DGEBA, seed), monomer(template.hardener, seed)
    sites = template.donors(p)
    n = min(sites)
    if (
        len(sites) != template.sites_per_hardener
        or len(sites[n]) != template.max_additions
    ):
        raise ValueError(
            "Reference hardener functionality disagrees with reaction template"
        )
    c, o = next(iter(epoxides(d).items()))
    product = p
    reactions = []
    for h in sites[n]:
        offset = product.GetNumAtoms()
        pre = Chem.CombineMols(product, d)
        event = template.event_type(offset + c, offset + o, n, h)
        product = react(pre, [event], template)
        reactions.append((pre, event, product))
    return d, p, reactions


def reference_models(seed, template):
    d, p, rx = reference_reactions(seed, template)
    return list(
        zip(template.reference_names, [d, p, *(post for _, _, post in rx)], strict=True)
    )
