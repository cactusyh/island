"""Explicit J4 chemical motifs; supplements, never reinterprets, J3 records.

Source atom/equivalence/automatic row evidence is recorded in each rule. These
are deliberately bounded chemical domains, not universal meanings of labels.
"""

from .automatic import _components
from .domains import recognize_domains
from .expanded import MASSES, NUMBERS, environments

PROFILE_NAME = "island_pcff_source_graph_v4"
RULES = {
    "h": "neutral isolated H-F/Cl/Br/I single bond, standard masses; atom 98; increments 475,697,769,790; ordinary 240; automatic 366; no ionic dissociation model",
    "s'": "neutral acyclic terminal thiocarbonyl S=C(H/C)2; atom 166; ordinary 308; automatic 426; increments 580,607,635",
    "nh+": "one protonated N in isolated six-membered aromatic C5NH ring, unsubstituted; atom 134; ordinary 276; automatic 399",
    "c+": "central neutral C in explicit guanidinium C(NH2)3 resonance graph, one C=N+ and two C-N; atom 68; ordinary 210; automatic 336",
    "nr": "terminal NH2 at declared guanidinium center; one formal +1 N with C=N; atom 140; ordinary 282; automatic 405; n2 alias not selected",
    "n3n": "neutral nonaromatic three-single-bond N in three-membered ordinary lactam ring; atom 122; ordinary 264; automatic n mapping; no generic strained-amine override",
    "n4n": "neutral nonaromatic three-single-bond N in four-membered ordinary lactam ring; atom 125; ordinary 267; automatic n mapping; no generic strained-amine override",
}
EVIDENCE = {
    "rules": RULES,
    "precedence": "whole verified motifs override only their own sites; unresolved component blocks assignment; J3 remains unchanged",
    "charge_convention": "existing exact/ordinary/automatic increments and zero base charge; component formal-charge agreement required, no normalization",
    "scope": "sulfur oxidation, phosphorus charge states, imidazolium aliases and other resonance environments remain unresolved unless an older explicit rule applies",
}


def recognize_organic(graph):
    answers, env, issues = recognize_domains(graph, defer_components=True)
    atoms, adj, bonds, _ = environments(graph)
    added = {}

    def el(i):
        return atoms[i]["element"]

    def order(i, j):
        return bonds[frozenset((i, j))]["order"]

    def standard(i, charge=0):
        a = atoms[i]
        return (
            el(i) in NUMBERS
            and a["atomic_number"] == NUMBERS[el(i)]
            and a["formal_charge"] == charge
            and not a["metadata"]["isotope"]
            and not a["metadata"]["radical_electrons"]
            and abs(a["mass"] - MASSES[el(i)]) <= 0.02
        )

    def hydrogen(i):
        return (
            standard(i)
            and el(i) == "H"
            and env[i]["bond_orders"] == [1]
            and not env[i]["aromatic"]
            and env[i]["smallest_ring"] is None
        )

    for i in atoms:
        e = env[i]
        if hydrogen(i):
            j = next(iter(adj[i]))
            if (
                el(j) in ("F", "Cl", "Br", "I")
                and standard(j)
                and env[j]["bond_orders"] == [1]
                and not env[j]["aromatic"]
                and env[j]["smallest_ring"] is None
            ):
                added.update({i: "h", j: el(j).lower()})
        if not standard(i):
            continue
        if (
            el(i) == "S"
            and e["bond_orders"] == [2]
            and not e["aromatic"]
            and e["smallest_ring"] is None
        ):
            j = next(iter(adj[i]))
            others = adj[j] - {i}
            if (
                el(j) == "C"
                and standard(j)
                and env[j]["bond_orders"] == [1, 1, 2]
                and not env[j]["aromatic"]
                and env[j]["smallest_ring"] is None
                and all(standard(k) and el(k) in ("C", "H") for k in others)
            ):
                added[i] = "s'"
                added[j] = "c=" if all(el(k) == "H" for k in others) else "c=2"
        if (
            el(i) == "N"
            and e["smallest_ring"] in (3, 4)
            and not e["aromatic"]
            and e["bond_orders"] == [1, 1, 1]
            and all(standard(j) and el(j) in ("C", "H") for j in adj[i])
        ):
            carbonyls = [
                j
                for j in adj[i]
                if el(j) == "C"
                and env[j]["bond_orders"] == [1, 1, 2]
                and env[j]["neighbor_elements"] == {"C": 1, "N": 1, "O": 1}
                and env[j]["smallest_ring"] == e["smallest_ring"]
                and any(
                    el(k) == "O"
                    and standard(k)
                    and order(j, k) == 2
                    and len(adj[k]) == 1
                    for k in adj[j]
                )
            ]
            if len(carbonyls) == 1:
                added[i] = "n3n" if e["smallest_ring"] == 3 else "n4n"
                added.update({j: "hn" for j in adj[i] if hydrogen(j)})
        if (
            el(i) == "C"
            and not e["aromatic"]
            and e["smallest_ring"] is None
            and e["bond_orders"] == [1, 1, 2]
            and all(el(j) == "N" for j in adj[i])
        ):
            ns = adj[i]
            if all(
                standard(j, 1 if order(i, j) == 2 else 0)
                and not env[j]["aromatic"]
                and env[j]["smallest_ring"] is None
                and len(adj[j]) == 3
                and all(hydrogen(k) for k in adj[j] - {i})
                for j in ns
            ):
                added[i] = "c+"
                for j in ns:
                    added[j] = "nr"
                    added.update({k: "hn" for k in adj[j] - {i}})
    # Isolated pyridinium is a whole-ring rule: no fused/extra hetero rings,
    # substituents, isotope, radical, or alternate formal-charge placement.
    for i in atoms:
        if not (
            el(i) == "N"
            and standard(i, 1)
            and env[i]["aromatic"]
            and env[i]["smallest_ring"] == 6
            and env[i]["bond_orders"] == [1, 1.5, 1.5]
        ):
            continue
        ring, todo = {i}, [i]
        while todo:
            for j in adj[todo.pop()]:
                if env[j]["aromatic"] and j not in ring:
                    ring.add(j)
                    todo.append(j)
        if len(ring) != 6 or any(el(j) != "C" or not standard(j) for j in ring - {i}):
            continue
        if not all(
            env[j]["bond_orders"] == [1, 1.5, 1.5]
            and len(adj[j] & ring) == 2
            and len(adj[j] - ring) == 1
            and all(hydrogen(k) for k in adj[j] - ring)
            for j in ring
        ):
            continue
        added[i] = "nh+"
        added.update({j: "cp" for j in ring - {i}})
        for j in ring:
            added.update({k: "hn" if j == i else "hc" for k in adj[j] - ring})
    issues = [
        d for d in issues if not (len(d["sites"]) == 1 and d["sites"][0] in added)
    ]
    answers.update(added)
    blocked = {
        i
        for c in _components(adj)
        if any(set(d["sites"]) & set(c) for d in issues)
        for i in c
    }
    return {i: t for i, t in answers.items() if i not in blocked}, env, issues
