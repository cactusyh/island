"""Independent declared CHON fixture labels and raw-source verification.

No production PCFF perception, selection or numerical coefficients are used.
This bounded experiment reader is not a runtime or general chemical typer.
"""

from collections import Counter

from pcff_j13_reference import require
from pcff_j13_reference import selection_check as check_rows


def labels_for(system):
    atoms = system.topology.sites
    neighbors = {i: [] for i in atoms}
    orders = {i: [] for i in atoms}
    for b in system.topology.bonds.values():
        for a, c in ((b.site1, b.site2), (b.site2, b.site1)):
            neighbors[a].append(c)
            orders[a].append(b.order)

    def acyl(i):
        return (
            atoms[i].element == "C"
            and sorted(orders[i]) == [1, 1, 2]
            and Counter(atoms[j].element for j in neighbors[i]) == {"N": 1, "O": 2}
        )

    labels = {}
    for i, atom in atoms.items():
        require(
            atom.formal_charge == 0
            and not atom.metadata.get("isotope")
            and not atom.metadata.get("radical_electrons"),
            "Outside declared reference state",
        )
        ns = neighbors[i]
        e = Counter(atoms[j].element for j in ns)
        if atom.element == "C":
            if atom.metadata.get("aromatic"):
                require(
                    sorted(orders[i]) == [1, 1.5, 1.5],
                    "Reference aromatic carbon mismatch",
                )
                labels[i] = "cp"
            elif acyl(i):
                labels[i] = "c_2"
            else:
                require(sorted(orders[i]) == [1] * 4, "Reference carbon mismatch")
                labels[i] = "c" + str(e["H"]) if e["H"] else "c"
        elif atom.element == "N":
            require(
                sorted(orders[i]) == [1] * 3 and not atom.metadata.get("aromatic"),
                "Reference N mismatch",
            )
            labels[i] = (
                "n_2"
                if any(acyl(j) for j in ns)
                else "nn"
                if any(atoms[j].metadata.get("aromatic") for j in ns)
                else "na"
            )
        elif atom.element == "O":
            if sorted(orders[i]) == [2]:
                labels[i] = "o_1"
            else:
                require(
                    sorted(orders[i]) == [1, 1] and e == {"C": 2},
                    "Reference O mismatch",
                )
                labels[i] = "o_2" if any(acyl(j) for j in ns) else "oc"
    for i, atom in atoms.items():
        if atom.element == "H":
            require(len(neighbors[i]) == 1, "Reference H mismatch")
            label = labels[neighbors[i][0]]
            require(
                label in {"c", "c1", "c2", "c3", "cp", "na", "nn", "n_2"},
                "Reference H parent mismatch",
            )
            labels[i] = (
                "hn2" if label == "n_2" else "hn" if label in {"na", "nn"} else "hc"
            )
    require(set(labels) == set(atoms), "Incomplete independent labels")
    return labels


def selection_check(system, model, source):
    return check_rows(system, model, source, reference_labels=labels_for(system))
