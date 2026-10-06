"""J3 graph predicates. No element-only generic fallback or charge repair.

Rules supplement the frozen J1/J2 predicates. An unresolved atom blocks its
component. Source descriptions bound each rule; narrower conditions avoid
claiming all chemistry associated with a generic source label.
"""

from .automatic import _components
from .expanded import MASSES, NUMBERS, environments, recognize

PROFILE_NAME = "island_pcff_source_graph_v3"
RULES = {
    "h": "neutral explicit H-H single bond only; FRC 98,789,1741",
    "s": "neutral nonaromatic H-S-H, two single bonds; FRC 167,794,1743,2402",
    "dw": "isotope=2 H, mass 2.014..2.014102, bonded to neutral fully deuterated water; FRC 95,238,364",
    "n": "neutral nonaromatic three-single-bond N attached to exactly one ordinary amide carbonyl and C/H; exclude carbamate/urea; FRC 116,259,800",
    "npc": "neutral aromatic five/six-ring N, two aromatic bonds and one single external C, no H; FRC 142,281; pinned LUNAR 876-881",
    "nn": "neutral nonring three-single-bond N attached to aromatic C, other C/H, no carbonyl; FRC 140,279; pinned LUNAR 904-907; nb alias remains unresolved",
    "n=": "neutral acyclic C=N-H with one single H; FRC 127; pinned LUNAR 795-797",
    "n=1": "neutral acyclic CH2=N-C, single substituent C; FRC 128; pinned LUNAR 800-802",
    "n=2": "neutral acyclic substituted C=N-C, no terminal CH2; FRC 129; pinned LUNAR 805-807",
    "p": "neutral tetra-coordinate P, one terminal P=O and three single C/H/O neighbors; FRC 163,305; no inference for other oxidation/coordination states",
    "hp": "single H on the declared neutral tetra-coordinate P=O motif; FRC 111,253,792",
    "o=": "terminal O double bonded to the declared neutral tetra-coordinate P=O motif; FRC 148,290; O=O remains outside this profile",
    "o": "neutral acyclic peroxide oxygen, one O single bond and one H/C single bond; FRC 143,285,961",
    "ho2": "H on the declared peroxide/ordinary carboxylic acid OH; FRC 108,250; pinned LUNAR 513-516; alcohol ho remains unchanged",
}
EVIDENCE = {
    "rules": RULES,
    "precedence": "specific declared motif supplements unresolved J2 sites; amine hn2 retained; isotope exception only verified D2O",
    "base_charge": "zero only; source endpoint increments must reproduce component formal charge; otherwise incomplete",
    "unresolved": "charged heterocycles, guanidinium resonance, generic aliases, sulfur oxidation, metal/zeolite/base-charge conventions remain explicit",
    "isotope": "D is encoded as element H, atomic_number 1, isotope 2, mass [2.014,2.014102]; source dw element D is an explicit representation bridge, no mass mutation",
}


def recognize_domains(graph):
    answers, env, issues = recognize(
        graph, elemental_halogens=True, defer_components=True
    )
    atoms, adj, bonds, _ = environments(graph)
    additions = {}

    def el(i):
        return atoms[i]["element"]

    def order(i, j):
        return bonds[frozenset((i, j))]["order"]

    def normal(i):
        a = atoms[i]
        return (
            a["formal_charge"] == 0
            and el(i) in NUMBERS
            and a["atomic_number"] == NUMBERS[el(i)]
            and not (el(i) == "H" and env[i]["aromatic"])
            and not a["metadata"]["isotope"]
            and not a["metadata"]["radical_electrons"]
            and abs(a["mass"] - MASSES[el(i)]) <= 0.02
        )

    def carbonyl(i):
        return (
            el(i) == "C"
            and len(adj[i]) == 3
            and any(el(j) == "O" and order(i, j) == 2 for j in adj[i])
        )

    def ordinary_amide(i):
        return carbonyl(i) and sorted(el(j) for j in adj[i]) in (
            ["C", "N", "O"],
            ["H", "N", "O"],
        )

    def poxo(i):
        return (
            normal(i)
            and el(i) == "P"
            and not env[i]["aromatic"]
            and len(adj[i]) == 4
            and env[i]["bond_orders"] == [1, 1, 1, 2]
            and sum(
                el(j) == "O" and order(i, j) == 2 and len(adj[j]) == 1 for j in adj[i]
            )
            == 1
            and all(el(j) in ("C", "H", "O") and normal(j) for j in adj[i])
        )

    def peroxide(i):
        return (
            normal(i)
            and el(i) == "O"
            and not env[i]["aromatic"]
            and env[i]["smallest_ring"] is None
            and env[i]["bond_orders"] == [1, 1]
            and sum(el(j) == "O" for j in adj[i]) == 1
            and all(el(j) in ("H", "C", "O") for j in adj[i])
        )

    def acid_oxygen(i):
        return (
            normal(i)
            and el(i) == "O"
            and env[i]["bond_orders"] == [1, 1]
            and not env[i]["aromatic"]
            and env[i]["neighbor_elements"] == {"C": 1, "H": 1}
            and any(carbonyl(j) for j in adj[i] if el(j) == "C")
        )

    def deuterium(i):
        a = atoms[i]
        return (
            el(i) == "H"
            and a["atomic_number"] == 1
            and a["formal_charge"] == 0
            and a["metadata"]["isotope"] == 2
            and not env[i]["aromatic"]
            and not a["metadata"]["radical_electrons"]
            and 2.014 <= a["mass"] <= 2.014102
            and env[i]["bond_orders"] == [1]
        )

    for i in atoms:
        e = env[i]
        ns = adj[i]
        orders = e["bond_orders"]
        element = el(i)
        if deuterium(i):
            parent = next(iter(ns))
            if (
                normal(parent)
                and el(parent) == "O"
                and len(adj[parent]) == 2
                and all(deuterium(j) for j in adj[parent])
            ):
                additions[i] = "dw"
            continue
        if not normal(i):
            continue
        if (
            element == "N"
            and not e["aromatic"]
            and orders == [1, 1, 1]
            and all(el(j) in ("C", "H") for j in ns)
        ):
            acyl = [j for j in ns if carbonyl(j)]
            if len(acyl) == 1 and ordinary_amide(acyl[0]):
                additions[i] = "n"
            elif (
                not acyl
                and e["smallest_ring"] is None
                and any(env[j]["aromatic"] for j in ns)
            ):
                additions[i] = "nn"
        elif (
            element == "N"
            and e["aromatic"]
            and e["smallest_ring"] in (5, 6)
            and orders == [1, 1.5, 1.5]
            and all(el(j) == "C" for j in ns)
        ):
            additions[i] = "npc"
        elif (
            element == "N"
            and not e["aromatic"]
            and e["smallest_ring"] is None
            and orders == [1, 2]
        ):
            double = next(j for j in ns if order(i, j) == 2)
            single = next(j for j in ns if order(i, j) == 1)
            if (
                el(double) == "C"
                and normal(double)
                and env[double]["bond_orders"] == [1, 1, 2]
                and not env[double]["aromatic"]
                and env[double]["smallest_ring"] is None
                and all(el(j) in ("C", "H") for j in adj[double] - {i})
                and el(single) in ("C", "H")
            ):
                additions[i] = (
                    "n="
                    if el(single) == "H"
                    else "n=1"
                    if env[double]["hydrogens"] == 2
                    else "n=2"
                )
                additions[double] = (
                    "c="
                    if env[double]["hydrogens"] == 2
                    else "c=1"
                    if el(single) == "H"
                    else "c=2"
                )
                # C=N carbon uses the same terminal/next-to-terminal distinction.
        elif element == "P" and poxo(i):
            additions[i] = "p"
        elif (
            element == "O"
            and not e["aromatic"]
            and orders == [2]
            and poxo(next(iter(ns)))
        ):
            additions[i] = "o="
        elif peroxide(i):
            additions[i] = "o"
        elif acid_oxygen(i):
            additions[i] = "oh"
        elif (
            element == "S"
            and not e["aromatic"]
            and orders == [1, 1]
            and all(el(j) == "H" for j in ns)
        ):
            additions[i] = "s"
        elif element == "H" and orders == [1]:
            parent = next(iter(ns))
            if (
                el(parent) == "H"
                and normal(parent)
                and env[parent]["bond_orders"] == [1]
            ):
                additions[i] = "h"
            elif poxo(parent):
                additions[i] = "hp"
            elif peroxide(parent) or acid_oxygen(parent):
                additions[i] = "ho2"
            elif el(parent) == "N" and env[parent]["bond_orders"] == [1, 2]:
                additions[i] = "hn"
    # Keep all invalid graph/charge/isotope/environment diagnostics that weren't
    # resolved by a declared motif. Never resurrect a partially valid component.
    issues = [
        d for d in issues if not (len(d["sites"]) == 1 and d["sites"][0] in additions)
    ]
    answers.update(additions)
    blocked = {
        i
        for c in _components(adj)
        if any(set(d["sites"]) & set(c) for d in issues)
        for i in c
    }
    return {i: t for i, t in answers.items() if i not in blocked}, env, issues
