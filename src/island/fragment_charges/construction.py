"""Optional RDKit construction of the small capped computational fragment."""

from island.charge_references.records import pack
from island.workflows.bundle import system_data

from .core import CappedFragment, boundary, fragment_data, need


@boundary
def prepare_capped_fragment(psmiles, *, seed=2026):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem

    from island.chemistry import RepeatUnit, from_rdkit

    need(type(seed) is int and 0 <= seed < 2**31, "Invalid fragment embedding seed")
    repeat = RepeatUnit.from_psmiles(psmiles)
    need(
        repeat.head.atom_map_number == 1 and repeat.tail.atom_map_number == 2,
        "Use explicit [*:1] and [*:2] attachments",
    )
    need(
        all(p.bond_order == 1 for p in repeat.attachment_points),
        "Only single attachment bonds supported",
    )
    need(len(Chem.GetMolFrags(repeat.mol)) == 1, "Disconnected repeat unsupported")
    for atom in repeat.mol.GetAtoms():
        need(
            atom.GetFormalCharge() == 0
            and atom.GetNumRadicalElectrons() == 0
            and atom.GetIsotope() == 0
            and atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED,
            "Charged, radical, isotope or assigned stereochemistry unsupported",
        )
    need(
        all(
            b.GetStereo() == Chem.BondStereo.STEREONONE
            and b.GetBondDir() == Chem.BondDir.NONE
            for b in repeat.mol.GetBonds()
        ),
        "Assigned bond stereochemistry unsupported",
    )
    rw = Chem.RWMol(repeat.mol)
    caps = {}
    for point in repeat.attachment_points:
        need(
            repeat.mol.GetAtomWithIdx(point.neighbor_atom_index).GetAtomicNum() > 1,
            "Attachment parent must be heavy",
        )
        rw.ReplaceAtom(point.dummy_atom_index, Chem.Atom("H"))
        caps[point.role] = {
            "site_id": point.dummy_atom_index + 1,
            "parent_local_index": point.neighbor_atom_index,
            "dummy_source_index": point.dummy_atom_index,
            "atom_map_number": point.atom_map_number,
        }
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    need(
        AllChem.EmbedMolecule(mol, parameters) == 0,
        "Fragment embedding failed; no retry",
    )
    converted = from_rdkit(mol)
    system = converted.system
    heavy = {
        a.GetIdx(): converted.rdkit_index_to_site_id[a.GetIdx()]
        for a in repeat.mol.GetAtoms()
        if a.GetAtomicNum() > 1
    }
    cap_ids = {v["site_id"] for v in caps.values()}
    native = {
        i: sorted(
            h
            for h in system.topology.neighbors(s)
            if system.topology.sites[h].element == "H" and h not in cap_ids
        )
        for i, s in heavy.items()
    }
    for i, s in heavy.items():
        system.topology.sites[s].metadata["source_repeat_atom_index"] = i
    for role, cap in caps.items():
        system.topology.sites[cap["site_id"]].metadata.update(
            fragment_cap_role=role, dummy_source_index=cap["dummy_source_index"]
        )
    for hs in native.values():
        for h in hs:
            system.topology.sites[h].metadata["fragment_native_hydrogen"] = True
    system.metadata.update(
        coordinate_source="fragment_etkdg_v3",
        fragment_psmiles=psmiles,
        fragment_seed=seed,
    )
    p = {
        "schema": "island_capped_fragment_v1",
        "system": system_data(system),
        "definition": psmiles,
        "seed": seed,
        "construction_version": "RDKit " + rdBase.rdkitVersion,
        "mapping": {"heavy": heavy, "native_h": native, "caps": caps},
        "data": {},
    }
    p["data"] = fragment_data(p)
    result = CappedFragment(pack(p))
    result.validate_integrity()
    return result
