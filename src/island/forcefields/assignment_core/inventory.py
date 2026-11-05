"""Stable site instances, unique environment keys and explicit coupling dependencies."""

from itertools import combinations

from island.forcefields.assignment_core.schema import TermRequest, identity
from island.forcefields.assignment_core.topology import as_topology, molecule
from island.forcefields.parameters.inventory import derive_interaction_inventory


def key(family, sites):
    sites = tuple(sites)
    if family in ("bond", "angle", "torsion"):
        sites = min(sites, sites[::-1])
    return family + ":" + ",".join(map(str, sites))


def environment(mol, sites):
    """Radius-two, stereo-aware rooted fragment; role order survives type aliases."""
    from rdkit import Chem

    fragment = Chem.Mol(mol)
    for atom in fragment.GetAtoms():
        atom.SetAtomMapNum(0)
    included = {i - 1 for i in sites}
    frontier = set(included)
    for _ in range(2):
        frontier = {
            n.GetIdx() for i in frontier for n in mol.GetAtomWithIdx(i).GetNeighbors()
        } - included
        included |= frontier
    for role, site in enumerate(sites, 1):
        fragment.GetAtomWithIdx(site - 1).SetAtomMapNum(role)
    text = Chem.MolFragmentToSmiles(
        fragment,
        atomsToUse=sorted(included),
        canonical=True,
        isomericSmiles=True,
        allHsExplicit=True,
    )
    roles = [
        (
            mol.GetAtomWithIdx(i - 1).GetDegree(),
            mol.GetAtomWithIdx(i - 1).GetTotalValence(),
        )
        for i in sites
    ]
    return identity(("radius_two_role_ordered_v1", text, roles))


def requests(graph, typed, profile):
    typed.validate(graph, profile)
    topology = as_topology(graph)
    base = derive_interaction_inventory(topology)
    mol = molecule(graph)
    types = dict(typed.assignments)
    rows = {}
    environments = {}

    def add(family, sites, dependencies=()):
        sites = tuple(sites)
        rid = key(family, sites)
        if sites not in environments:
            environments[sites] = environment(mol, sites)
        row = TermRequest(
            rid,
            family,
            sites,
            tuple(types.get(i, "") for i in sites),
            environments[sites],
            tuple(sorted(dependencies)),
        )
        if rid in rows and rows[rid] != row:
            raise ValueError("conflicting_request_identity:" + rid)
        rows[rid] = row

    for i in graph.atom_order:
        add("nonbonded", (i,))
        add("charge", (i,))
    for sites in base.bonds:
        add("bond", sites)
    for a, b, c in base.angles:
        add("angle", (a, b, c))
        if profile.functional_class == "II":
            lengths = (key("bond", (a, b)), key("bond", (b, c)))
            add("bond_bond", (a, b, c), lengths)
            add("bond_angle", (a, b, c), (*lengths, key("angle", (a, b, c))))
    for a, b, c, d in base.proper_torsions:
        add("torsion", (a, b, c, d))
        if profile.functional_class == "II":
            ends = (key("bond", (a, b)), key("bond", (c, d)))
            bends = (key("angle", (a, b, c)), key("angle", (b, c, d)))
            for family, deps in (
                ("bond_bond_13", ends),
                ("end_bond_torsion", (*ends, key("torsion", (a, b, c, d)))),
                (
                    "middle_bond_torsion",
                    (key("bond", (b, c)), key("torsion", (a, b, c, d))),
                ),
                ("angle_torsion", (*bends, key("torsion", (a, b, c, d)))),
                ("angle_angle_torsion", (*bends, key("torsion", (a, b, c, d)))),
            ):
                add(family, (a, b, c, d), deps)
    for center in graph.atom_order:
        for a, b, c in combinations(sorted(topology.neighbors(center)), 3):
            # Conservative enumeration: omission requires an exact evidenced N/A row.
            add("improper", (a, center, b, c))
            if profile.functional_class == "II":
                for shared in (a, b, c):
                    left, right = sorted({a, b, c} - {shared})
                    add(
                        "angle_angle",
                        (left, center, shared, right),
                        (
                            key("angle", (left, center, shared)),
                            key("angle", (shared, center, right)),
                        ),
                    )
    if profile.functional_class == "II":
        paths = set()
        for a, b, c, d in base.proper_torsions:
            for e in sorted(topology.neighbors(d) - {a, b, c}):
                p = (a, b, c, d, e)
                paths.add(min(p, p[::-1]))
            for e in sorted(topology.neighbors(a) - {b, c, d}):
                p = (e, a, b, c, d)
                paths.add(min(p, p[::-1]))
        for a, b, c, d, e in sorted(paths):
            add(
                "torsion_torsion",
                (a, b, c, d, e),
                (key("torsion", (a, b, c, d)), key("torsion", (b, c, d, e))),
            )
    return tuple(rows[k] for k in sorted(rows))
