"""Opt-in, local ref-1 aliphatic amine convention; historical rules stay intact.

These predicates are independently implemented from source descriptions and
family evidence. They do not copy the LUNAR or pysimm implementations.
"""

from .expanded import MASSES, NUMBERS, environments
from .organic_domains import recognize_organic

PROFILE_NAME = "island_pcff_source_graph_v5"
EVIDENCE = {
    "convention": "bounded_ref1_aliphatic_amines_v1",
    "source_revision": "e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
    "atom_rows": {
        "hn": 105,
        "hn2": 106,
        "h*": 99,
        "h+": 100,
        "na": 130,
        "n4": 123,
        "n+": 118,
        "n_2": 129,
    },
    "neutral": "acyclic neutral closed-shell three-single-bond N with 1..3 saturated neutral C neighbors and remaining explicit H; na/hn ref1, hn ordinary and increment equivalence h*",
    "cation": "acyclic closed-shell formal +1 four-single-bond N with 1..4 saturated neutral C neighbors and remaining explicit neutral H; n4/h+ ref1, n4 equivalence n+; H formal charge stays zero",
    "precedence": "replace only verified N/H sites in an already resolved v4 component; no amide, aromatic, resonance, strained ring, isotope or radical override",
    "increments": {"c_na": 513, "hstar_na": 804, "c_nplus": 508, "hplus_nplus": 813},
    "family_evidence": "ref1 (Biosym 25-Dec-1991) labels and equivalences; hn2 ref8 and hn2-n_2 ref7 row816 belong to a distinct refined family, not an alias",
    "independent_comparison": {
        "pysimm_revision": "fb33814128189a99d3c8d7c4eddac2dfe44262fe",
        "pysimm_file_sha256": "4ed9372aef3446537e1388c1902922d953fc9d0b32bdbe3ca1c5d2b726a459c9",
        "scope": "pysimm/forcefield/pcff.py H-on-N selects h*; corroborates hydrogen family only, not all N predicates or protonated typing",
        "lunar_revision": "67dabeda9e6bd3cc8f968aa0c6a88730886138ef",
        "disagreement": "LUNAR PCFF.py selects hn2 on sp3 amines and hn on n4; v5 intentionally selects the coherent ref1 family, not claimed LUNAR agreement",
    },
    "limitations": "a declared source-family convention, not a universal/commercial PCFF type assignment; charges and all active model terms must independently validate",
}


def recognize_amines(graph):
    answers, env, issues = recognize_organic(graph)
    atoms, adj, _, _ = environments(graph)

    def standard(i, charge=0):
        a = atoms[i]
        return (
            a["element"] in ("C", "H", "N")
            and a["atomic_number"] == NUMBERS[a["element"]]
            and a["formal_charge"] == charge
            and not a["metadata"]["isotope"]
            and not a["metadata"]["radical_electrons"]
            and abs(a["mass"] - MASSES[a["element"]]) <= 0.02
            and not env[i]["aromatic"]
        )

    for i, a in atoms.items():
        if a["element"] != "N" or i not in answers:
            continue
        charge = a["formal_charge"]
        degree = 3 if charge == 0 else 4
        if (
            charge not in (0, 1)
            or not standard(i, charge)
            or env[i]["smallest_ring"] is not None
            or env[i]["bond_orders"] != [1] * degree
        ):
            continue
        carbons = [j for j in adj[i] if atoms[j]["element"] == "C"]
        hydrogens = [j for j in adj[i] if atoms[j]["element"] == "H"]
        if (
            not carbons
            or len(carbons) + len(hydrogens) != degree
            or not all(
                standard(j) and env[j]["bond_orders"] == [1] * 4 for j in carbons
            )
            or not all(standard(j) and env[j]["bond_orders"] == [1] for j in hydrogens)
            or not all(j in answers for j in adj[i])
        ):
            continue
        answers[i] = "na" if charge == 0 else "n4"
        rule = "ref1_neutral_amine" if charge == 0 else "ref1_tetrahedral_ammonium"
        env[i]["amine_rule"] = {
            "id": rule,
            "carbon_neighbors": sorted(carbons),
            "hydrogen_neighbors": sorted(hydrogens),
            "formal_charge": charge,
        }
        for j in hydrogens:
            answers[j] = "hn" if charge == 0 else "h+"
            env[j]["amine_rule"] = {"id": rule + "_hydrogen", "parent": i}
    return answers, env, issues
