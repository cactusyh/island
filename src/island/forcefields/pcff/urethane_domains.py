"""Opt-in source-family coherence for bounded urethanes and benzenoid amines.

Independent local predicates; no external typer, converter, or repeat metadata.
No equivalence/parameter/charge policy is changed.
"""

from .amine_domains import recognize_amines
from .expanded import MASSES, NUMBERS, environments

PROFILE_NAME = "island_pcff_source_graph_v6"
RULES = {
    "hn2": "neutral standard-state acyclic R-O-C(=O)-N(H/Csp3)2; R saturated C, N three single bonds; H on n_2 only; source ref7 urethane family, not a global H alias",
    "hn": "neutral acyclic three-single-bond aromatic-amine N with one/two H, one/two C neighbors including a verified six-carbon benzenoid ring; other C sp3 or same benzenoid environment; hn->h* ref1 family; no aromatic-ring N or acyl N",
}
EVIDENCE = {
    "parent_profile": "island_pcff_source_graph_v5",
    "convention": "bounded_urethane_ref7_benzenoid_amine_ref1_v1",
    "source_revision": "e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
    "rules": RULES,
    "atom_rows": {
        "c_2": 83,
        "o_1": 149,
        "o_2": 150,
        "n_2": 129,
        "hn2": 106,
        "nn": 137,
        "hn": 105,
        "h*": 99,
    },
    "increments": {
        "c_n_2": 512,
        "c_2_n_2": 690,
        "hn2_n_2": 816,
        "cp_nn": 729,
        "c_nn": 515,
        "hstar_nn": 806,
    },
    "authority": "FRC atom/equivalence families and references 1/7/8; neutral carbamate N is distinct from ordinary amide/urea N; no charge repair or availability-based generic fallback",
    "external_evidence": {
        "lunar_revision": "67dabeda9e6bd3cc8f968aa0c6a88730886138ef",
        "urethane": "PCFF.py n_2 on C3H7NO2 urethane example; v6 replaces whole-formula recognition with a bounded local motif, not claimed exact LUNAR parity",
        "aromatic": "PCFF.py nn and typing_functions.is_aromatic_amine_nitrogen exclude ring N and require direct aromatic carbon; v6 further bounds ring/charge/isotope/valence",
        "hydrogen_disagreement": "LUNAR generic non-sp3-amine N-H uses hn including urethane; v6 selects the separately described hn2/n_2 ref7 family. Aromatic nn/hn agrees with generic N-H family, not inherited hn2 decision",
        "pysimm_revision": "fb33814128189a99d3c8d7c4eddac2dfe44262fe",
        "scope": "pysimm PCFF.py corroborates carbamate c_2, carbonyl o_1, ester o_2 and generic N-bound h*, but does not establish the specialized n_2/hn2 convention",
    },
    "exclusions": "no urea, acylated aromatic N, cyclic urethane N, fused/hetero aromatic N attachment, ionic, radical or isotopic override; all inherited failures remain diagnostic",
}


def recognize_urethanes(graph):
    answers, env, issues = recognize_amines(graph)
    atoms, adj, bonds, _ = environments(graph)

    def el(i):
        return atoms[i]["element"]

    def ordinary(i, element, orders, aromatic=False):
        a = atoms[i]
        return (
            i in answers
            and el(i) == element
            and a["atomic_number"] == NUMBERS[element]
            and a["formal_charge"] == 0
            and not a["metadata"]["isotope"]
            and not a["metadata"]["radical_electrons"]
            and abs(a["mass"] - MASSES[element]) <= 0.02
            and env[i]["aromatic"] is aromatic
            and env[i]["bond_orders"] == orders
        )

    def carbon(i):
        return ordinary(i, "C", [1, 1, 1, 1])

    def hydrogen(i):
        return ordinary(i, "H", [1])

    def benzenoid(i):
        # Require the entire aromatic component to be one six-carbon cycle.
        # This excludes fused and hetero rings rather than guessing from degree.
        seen, todo = set(), [i]
        while todo:
            j = todo.pop()
            if j in seen:
                continue
            if not ordinary(j, "C", [1, 1.5, 1.5], aromatic=True):
                return False
            seen.add(j)
            todo.extend(
                k
                for k in adj[j]
                if bonds[frozenset((j, k))]["order"] == 1.5 and k not in seen
            )
        return len(seen) == 6 and all(sum(k in seen for k in adj[j]) == 2 for j in seen)

    for n in atoms:
        if not ordinary(n, "N", [1, 1, 1]) or env[n]["smallest_ring"] is not None:
            continue
        hs = [j for j in adj[n] if hydrogen(j)]
        cs = [j for j in adj[n] if el(j) == "C"]
        if len(hs) + len(cs) != 3:
            continue
        acyl = [j for j in cs if ordinary(j, "C", [1, 1, 2])]
        if len(acyl) == 1 and all(carbon(j) for j in cs if j not in acyl):
            c = acyl[0]
            os = [j for j in adj[c] if ordinary(j, "O", [2])]
            alkoxy = [j for j in adj[c] if ordinary(j, "O", [1, 1])]
            if len(os) != 1 or len(alkoxy) != 1 or env[c]["smallest_ring"] is not None:
                continue
            o = alkoxy[0]
            other = adj[o] - {c}
            if (
                len(other) != 1
                or not carbon(next(iter(other)))
                or env[o]["smallest_ring"] is not None
            ):
                continue
            # Source-family checks precede the H specialization, without aliases.
            if [answers[j] for j in (n, c, os[0], o)] != ["n_2", "c_2", "o_1", "o_2"]:
                continue
            env[n]["nitrogen_family_rule"] = {
                "id": "urethane_ref7",
                "carbonyl": c,
                "carbonyl_oxygen": os[0],
                "alkoxy_oxygen": o,
                "hydrogens": sorted(hs),
            }
            for h in hs:
                answers[h] = "hn2"
                env[h]["nitrogen_family_rule"] = {
                    "id": "urethane_ref7_hydrogen",
                    "parent": n,
                }
        elif (
            len(hs) in (1, 2)
            and cs
            and any(benzenoid(j) for j in cs)
            and all(carbon(j) or benzenoid(j) for j in cs)
            and answers[n] == "nn"
        ):
            env[n]["nitrogen_family_rule"] = {
                "id": "benzenoid_amine_ref1",
                "carbon_neighbors": sorted(cs),
                "hydrogens": sorted(hs),
            }
            for h in hs:
                answers[h] = "hn"
                env[h]["nitrogen_family_rule"] = {
                    "id": "benzenoid_amine_ref1_hydrogen",
                    "parent": n,
                }
    return answers, env, issues
