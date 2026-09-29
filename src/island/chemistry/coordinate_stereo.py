"""Shared assigned tetrahedral labels and independent Cartesian validation."""

import numpy as np

from island.exceptions import StereochemistryError


def assigned_cip_labels(system, converted=None):
    """Resolve explicit stored or graph tags using the existing RDKit convention."""
    if converted is None and not any(
        site.metadata.get("cip_label") is not None
        or site.metadata.get("chiral_tag") not in (None, "CHI_UNSPECIFIED")
        for site in system.topology.sites.values()
    ):
        return {}
    from rdkit import Chem

    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    converted = converted or system_to_rdkit_graph(system)
    Chem.AssignStereochemistry(converted.mol, cleanIt=True, force=True)
    expected = {}
    for site, index in converted.site_id_to_rdkit_index.items():
        atom = converted.mol.GetAtomWithIdx(index)
        stored = system.topology.sites[site].metadata.get("cip_label")
        observed = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        if stored is not None and stored not in ("R", "S"):
            raise StereochemistryError(
                f"Unsupported CIP label at site {site}: {stored}"
            )
        if stored and observed and stored != observed:
            raise StereochemistryError(
                f"Site {site} stored CIP {stored} conflicts with chemical graph {observed}"
            )
        if stored or observed:
            expected[site] = stored or observed
        elif system.topology.sites[site].metadata.get("chiral_tag") not in (
            None,
            "CHI_UNSPECIFIED",
        ):
            raise StereochemistryError(
                f"Assigned tetrahedral center {site} is indeterminate"
            )
    return expected


def validate_coordinate_stereochemistry(system, positions, expected, *, stage):
    """Erase tags and infer handedness from 3D; reject degenerate assigned centers."""
    if not expected:
        return
    from rdkit import Chem

    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    for center in expected:
        neighbors = sorted(system.topology.neighbors(center))
        if len(neighbors) not in (3, 4):
            raise StereochemistryError(
                f"{stage}: indeterminate tetrahedral site {center}"
            )
        anchor = positions[neighbors[3]] if len(neighbors) == 4 else positions[center]
        vectors = [np.array(positions[site]) - anchor for site in neighbors[:3]]
        volume6 = abs(float(np.dot(vectors[0], np.cross(vectors[1], vectors[2]))))
        if not np.isfinite(volume6) or volume6 < 1e-3:
            raise StereochemistryError(
                f"{stage}: planar/degenerate tetrahedral geometry at site {center}"
            )
    converted = system_to_rdkit_graph(system)
    mol = Chem.Mol(converted.mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for site, index in converted.site_id_to_rdkit_index.items():
        conf.SetAtomPosition(index, positions[site])
    mol.AddConformer(conf, assignId=True)
    Chem.RemoveStereochemistry(mol)
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    for site, label in expected.items():
        atom = mol.GetAtomWithIdx(converted.site_id_to_rdkit_index[site])
        actual = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        if actual != label:
            raise StereochemistryError(
                f"{stage}: site {site} stereochemistry expected {label}, observed {actual}"
            )
