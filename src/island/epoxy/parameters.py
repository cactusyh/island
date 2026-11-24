"""Transport native PCFF terms from mapped reaction products to periodic sites."""

import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem

from island.forcefields.assignment import assign_emc_graph
from island.forcefields.emc_graph import _section

from .chemistry import reference_reactions

SECTIONS = ("Bonds", "Angles", "Dihedrals", "Impropers")
KINDS = ("atom", "bond", "angle", "dihedral", "improper")


def rows(data, name):
    return [
        line.partition("#")[0].split() for line in _section(data, name).splitlines()
    ]


def neighborhoods(mol):
    """Rooted, canonical labelled two-hop graphs and their actual atom traversal."""
    adj = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
    # Include full degree/H count on boundary atoms. This distinguishes primary,
    # secondary and tertiary amines even at the outer edge of a local fragment.
    codes = [
        a.GetAtomicNum() * 128
        + a.GetDegree() * 16
        + a.GetTotalNumHs(includeNeighbors=True) * 2
        + int(a.GetIsAromatic())
        for a in mol.GetAtoms()
    ]
    if any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        raise ValueError(
            "This template set covers neutral, unlabelled epoxy and amine reactants"
        )
    for i, neighbors in enumerate(adj):
        ids = {i, *neighbors}
        for j in neighbors:
            ids.update(adj[j])
        sites = sorted(ids)
        local = Chem.RWMol()
        indices = {}
        for j in sites:
            atom = Chem.Atom(6)
            atom.SetNoImplicit(True)
            atom.SetIsotope(codes[j] + (10000 if j == i else 0))
            indices[j] = local.AddAtom(atom)
        for j in sites:
            for k in adj[j]:
                if k in ids and k > j:
                    local.AddBond(
                        indices[j],
                        indices[k],
                        mol.GetBondBetweenAtoms(j, k).GetBondType(),
                    )
        frag = local.GetMol()
        frag.UpdatePropertyCache(strict=False)
        key = Chem.MolToSmiles(frag)
        order = [sites[j] for j in json.loads(frag.GetProp("_smilesAtomOutputOrder"))]
        yield i, key, order


def compile_library(directory, *, field="pcff", seed=2026):
    """Build a versioned, portable library; fail on any missing native parameter."""
    if field != "pcff":
        raise ValueError("Periodic template transport is qualified for pcff")
    from dataclasses import asdict

    out = Path(directory)
    out.mkdir(parents=True, exist_ok=False)
    d, p, reactions = reference_reactions(seed)
    samples = [
        ("dgeba", d),
        ("pacm", p),
        ("secondary", reactions[0][2]),
        ("tertiary", reactions[1][2]),
    ]
    coeffs = {k: {} for k in KINDS}
    pair = {}
    library = {}
    conventions = None
    for name, mol in samples:
        mol = Chem.Mol(mol)
        mol.RemoveAllConformers()
        native = assign_emc_graph(mol, field, seed=seed)
        (out / (name + ".json")).write_text(json.dumps(asdict(native)))
        commands = [
            line.split("#")[0].split() for line in native.parameters.splitlines()
        ]
        styles = [
            r
            for r in commands
            if r
            and (
                r[0].endswith("_style")
                or r[0] in ("pair_modify", "special_bonds")
                or (r[0] == "variable" and r[1] in ("cutoff", "charge_cutoff"))
            )
        ]
        if conventions is not None and styles != conventions:
            raise ValueError(
                "Native force-field conventions changed between reference molecules"
            )
        conventions = styles
        mapping = {}
        for k in KINDS:
            cmd = "mass" if k == "atom" else k + "_coeff"
            lines = defaultdict(list)
            for line in native.parameters.splitlines():
                pre, _, label = line.partition("#")
                r = pre.split()
                if r and r[0] == cmd:
                    lines[int(r[1])].append((r[2:], label.strip()))
            for local, records in lines.items():
                key = json.dumps(records)
                table = coeffs[k]
                mapping[k, local] = table.setdefault(key, len(table) + 1)
        for r in commands:
            if r and r[0] == "pair_coeff":
                ij = tuple(sorted(mapping["atom", int(t)] for t in r[1:3]))
                if pair.setdefault(ij, r[3:]) != r[3:]:
                    raise ValueError("Conflicting pair coefficients")
        atom_ids = {i: group[0] for i, group in native.site_groups.items()}
        atoms = {
            atom_ids[int(r[0])]: [mapping["atom", int(r[2])], float(r[3])]
            for r in rows(native.data, "Atoms")
        }
        owned = defaultdict(list)
        for section, k in zip(SECTIONS, KINDS[1:]):
            for _, t, *rr in rows(native.data, section):
                sites = [atom_ids[int(x)] for x in rr]
                anchors = (
                    sites
                    if k == "bond"
                    else sites[1:3]
                    if k == "dihedral"
                    else sites[1:2]
                )
                for a in anchors:
                    owned[a].append((section, mapping[k, int(t)], sites))
        for i, key, order in neighborhoods(mol):
            ix = {j: k for k, j in enumerate(order)}
            record = [
                atoms[i],
                [[sec, t, [ix[s] for s in ss]] for sec, t, ss in owned[i]],
            ]
            previous = library.setdefault(key, record)
            if previous[0] != record[0]:
                raise ValueError(f"Ambiguous PCFF atom environment at {name} site {i}")
    lines = [" ".join(r) for r in conventions]
    for k, table in coeffs.items():
        cmd = "mass" if k == "atom" else k + "_coeff"
        for key, i in table.items():
            for values, label in json.loads(key):
                lines.append(f"{cmd} {i} {' '.join(values)} # {label}")
    for (i, j), val in sorted(pair.items()):
        lines.append(f"pair_coeff {i} {j} {' '.join(val)}")
    contents = {
        "stencils": library,
        "counts": {k: len(v) for k, v in coeffs.items()},
        "parameters": "\n".join(lines) + "\n",
        "field": field,
        "field_sha256": native.field_sha256,
        "emc_version": native.emc_version,
    }
    (out / "library.json").write_text(json.dumps(contents))
    return contents


@lru_cache(maxsize=4)
def _native_labels(parameters):
    # Labels are part of the native coefficient records. E.g. the carbon
    # beyond an ether may acquire a new PCFF type after epoxide opening,
    # even when a two-hop skeleton around the ether remains isomorphic.
    labels = {k: {} for k in KINDS}
    ids = {k: defaultdict(set) for k in KINDS}
    for line in parameters.splitlines():
        pre, _, comment = line.partition("#")
        row = pre.split()
        if not row:
            continue
        k = "atom" if row[0] == "mass" else row[0].removesuffix("_coeff")
        if k not in labels:
            continue
        t = int(row[1])
        names = tuple(comment.strip().split(","))
        if k == "improper":
            # EMC names Class-II improper types center-first; LAMMPS lists it second.
            names = (names[1], names[0], *names[2:])
        labels[k][t] = names
        ids[k][names].add(t)
    return labels, ids


def assign_network(mol, lib):
    """Type an explicit-H graph locally, independent of global embeddability."""
    atoms = []
    terms = {s: [] for s in SECTIONS}
    for i, key, order in neighborhoods(mol):
        if key not in lib["stencils"]:
            raise ValueError(f"Uncovered local chemistry at site {i}: {key}")
        atom, local = lib["stencils"][key]
        atoms.append(atom)
        for sec, t, indices in local:
            sites = [order[k] for k in indices]
            if sec == "Bonds" and i != min(sites):
                continue
            if sec == "Dihedrals" and i != min(sites[1:3]):
                continue
            terms[sec].append((t, sites))
    labels, ids = _native_labels(lib["parameters"])
    for section, k in zip(SECTIONS, KINDS[1:], strict=True):
        for index, (t, sites) in enumerate(terms[section]):
            names = tuple(labels["atom"][atoms[j][0]][0] for j in sites)
            if names == labels[k][t]:
                continue
            options = ids[k].get(names)
            if not options and k != "improper":
                sites = list(reversed(sites))
                options = ids[k].get(names[::-1])
            if not options or len(options) != 1:
                raise ValueError(
                    f"Ambiguous or missing {k} coefficients for sites {sites}, types {names}"
                )
            terms[section][index] = (next(iter(options)), sites)
    if (
        len(terms["Bonds"]) != mol.GetNumBonds()
        or abs(sum(q for _, q in atoms) - Chem.GetFormalCharge(mol)) > 1e-6
    ):
        raise ValueError("Mapped PCFF graph failed bond or total charge accounting")
    return atoms, terms


def data_text(atoms, terms, counts, positions, box, membership, velocities=None):
    if len(atoms) != len(positions) or not np.isfinite(positions).all():
        raise ValueError("Invalid particle positions")
    lines = ["ISLAND periodic epoxy network (real units)", "", f"{len(atoms)} atoms"]
    lines.extend(f"{len(v)} {k.lower()}" for k, v in terms.items())
    lines.extend(f"{v} {k} types" for k, v in counts.items())
    lines.extend(f"0 {v:.15g} {k}lo {k}hi" for k, v in zip("xyz", box, strict=True))
    lines.extend(["", "Atoms", ""])
    for i, ((t, q), pos, m) in enumerate(
        zip(atoms, positions, membership, strict=True), 1
    ):
        lines.append(f"{i} {m} {t} {q} " + " ".join(f"{x:.15g}" for x in pos))
    for sec, rows_ in terms.items():
        lines.extend(["", sec, ""])
        for i, (t, sites) in enumerate(rows_, 1):
            lines.append(f"{i} {t} " + " ".join(str(v + 1) for v in sites))
    if velocities is not None:
        lines.extend(["", "Velocities", ""])
        for i, vel in enumerate(velocities, 1):
            lines.append(str(i) + " " + " ".join(f"{v:.15g}" for v in vel))
    return "\n".join(lines) + "\n"
