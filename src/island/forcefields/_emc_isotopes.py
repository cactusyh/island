"""Carry explicit isotope masses onto EMC's all-atom topology without changing its potential."""

import re
from collections import defaultdict
from itertools import product

from rdkit import Chem

from island.builders import build_linear_polymer
from island.chemistry import system_to_rdkit_graph
from island.exceptions import ForceFieldError


def _section(data, name):
    header = "\n" + name + "\n\n"
    return data.split(header, 1)[1].split("\n\n", 1)[0]


def _skeleton(numbers, edges):
    graph = Chem.RWMol()
    for number in numbers:
        atom = Chem.Atom(number)
        atom.SetNoImplicit(True)
        graph.AddAtom(atom)
    for a, b in edges:
        graph.AddBond(a, b, Chem.BondType.SINGLE)
    mol = graph.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol


def apply_isotope_masses(data, params, psmiles, dp):
    system = build_linear_polymer(
        psmiles, dp=dp, head_end_group="*[H]", tail_end_group="*[H]", generate_3d=False
    )
    ref = system_to_rdkit_graph(system).mol
    numbers = [a.GetAtomicNum() for a in ref.GetAtoms()]
    skeleton = _skeleton(
        numbers, [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in ref.GetBonds()]
    )
    # LAMMPS bonds lack chemical bond orders. Check non-uniformly labelled
    # symmetry classes before trusting an arbitrary full-graph isomorphism.
    classes = defaultdict(list)
    for i, rank in enumerate(Chem.CanonicalRankAtoms(skeleton, breakTies=False)):
        classes[rank].append(i)
    mixed = [
        group
        for group in classes.values()
        if len({ref.GetAtomWithIdx(i).GetIsotope() for i in group}) > 1
    ]
    chemical_ranks = Chem.CanonicalRankAtoms(
        ref, breakTies=False, includeIsotopes=False
    )
    for group in mixed:
        if len({chemical_ranks[i] for i in group}) > 1 or "@" in psmiles:
            raise ForceFieldError("Ambiguous isotope placement in the EMC connectivity")
    table = Chem.GetPeriodicTable()
    masses = _section(data, "Masses")
    elements = {}
    for row in masses.splitlines():
        t, mass = row.split()[:2]
        matches = [
            n
            for n in set(numbers)
            if abs(table.GetAtomicWeight(n) - float(mass)) < 0.15
        ]
        if len(matches) != 1:
            raise ForceFieldError("Cannot identify an EMC element for isotope mapping")
        elements[int(t)] = matches[0]

    atoms = _section(data, "Atoms")
    rows = [row.split() for row in atoms.splitlines()]
    atom_ids = [int(row[0]) for row in rows]
    indices = {site: i for i, site in enumerate(atom_ids)}
    molecules = defaultdict(list)
    for i, row in enumerate(rows):
        molecules[int(row[1])].append(i)
    bonds = defaultdict(list)
    for row in _section(data, "Bonds").splitlines():
        _, _, a, b = row.split()[:4]
        i, j = indices[int(a)], indices[int(b)]
        bonds[int(rows[i][1])].append((i, j))

    isotope_rows = []
    for molecule, group in molecules.items():
        local = {i: n for n, i in enumerate(group)}
        observed = _skeleton(
            [elements[int(rows[i][2])] for i in group],
            [(local[a], local[b]) for a, b in bonds[molecule]],
        )
        match = observed.GetSubstructMatch(skeleton)
        if (
            len(match) != len(group)
            or len(match) != len(numbers)
            or observed.GetNumBonds() != skeleton.GetNumBonds()
        ):
            raise ForceFieldError(
                "EMC connectivity differs from the isotope-labelled chain"
            )
        for equivalent in mixed:
            candidates = [rows[group[match[i]]] for i in equivalent]
            if len({(row[2], float(row[3])) for row in candidates}) != 1:
                raise ForceFieldError(
                    "Ambiguous EMC type or charge for isotope placement"
                )
        for atom in ref.GetAtoms():
            iso = atom.GetIsotope()
            if not iso:
                continue
            row = rows[group[match[atom.GetIdx()]]]
            isotope_rows.append((row, atom.GetAtomicNum(), iso))
    return split_isotope_types(data, params, rows, isotope_rows)


def split_isotope_types(data, params, rows, isotope_rows):
    table = Chem.GetPeriodicTable()
    masses = _section(data, "Masses")
    atoms = _section(data, "Atoms")
    parent_types = [int(row.split()[0]) for row in masses.splitlines()]
    replacements = {}
    splits = defaultdict(list)
    new_masses = []
    for row, number, iso in isotope_rows:
        parent = int(row[2])
        key = parent, iso
        if key not in replacements:
            child = max(parent_types) + len(replacements) + 1
            mass = table.GetMassForIsotope(number, iso)
            if mass <= 0:
                raise ForceFieldError(
                    f"Unknown isotope mass for {table.GetElementSymbol(number)}-{iso}"
                )
            replacements[key] = child
            splits[parent].append(child)
            new_masses.append(
                f"{child} {mass:.10f} # {table.GetElementSymbol(number)}-{iso}"
            )
        row[2] = str(replacements[key])

    # Clone explicit LJ overrides as well as diagonal rows. All bonded type IDs
    # and force constants remain those emitted by the provider.
    extra_coeffs = []
    for line in params.splitlines():
        fields = line.split()
        if not fields or fields[0] != "pair_coeff":
            continue
        i, j = map(int, fields[1:3])
        for a, b in sorted(
            {tuple(sorted(pair)) for pair in product([i, *splits[i]], [j, *splits[j]])}
        ):
            if (a, b) != (i, j):
                extra_coeffs.append(
                    " ".join(["pair_coeff", str(a), str(b), *fields[3:]])
                )
    data = data.replace(
        masses + "\n\n", masses + "\n" + "\n".join(new_masses) + "\n\n", 1
    )
    data = data.replace(
        atoms + "\n\n", "\n".join(" ".join(r) for r in rows) + "\n\n", 1
    )
    data = re.sub(
        r"(?m)^(\s*)\d+\s+atom types$",
        rf"\g<1>{len(parent_types) + len(new_masses)} atom types",
        data,
        count=1,
    )
    return data, params + "\n" + "\n".join(extra_coeffs) + "\n"
