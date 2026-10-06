"""Versioned graph rules for the full-source coverage project.

This is an explicitly partial implementation of the full pinned inventory.
Every rule is a chemical predicate; unresolved chemistry has no catch-all type.
Neither graph recognition nor charge completion asserts parameter coverage.
"""

from collections import Counter, deque
from copy import deepcopy
from math import fsum

from island.charge_references.records import pack

from .automatic import _components, chemical_graph, graph_system
from .charges import bond_selection, identity
from .source import FLAGS, PIN, boundary, records, require, select

PROFILE_NAME = "island_pcff_source_graph_v1"
TYPING_SCHEMA = "island_pcff_source_typing_v1"
CHARGE_SCHEMA = "island_pcff_source_charges_v1"
PROFILE = {
    "name": PROFILE_NAME,
    "source_sha256": PIN["sha256"],
    "implementation": "verified_local_environments_shortest_cycles_v1",
    "hydroxyl_convention": "ordinary_alcohol_oh_ho_v1",
    "charge_policy": "zero_base_direct_then_bond_equivalence_v1",
    "auto_equivalence_charges": "not_enabled_reference_branch_unreachable",
    "rule_evidence": "pinned FRC atom descriptions and LUNAR PCFF.py 67dabeda9e6bd3cc8f968aa0c6a88730886138ef; independent rules, no assumed assignments",
    "explicit_aliases": {"hc": ["h"], "n4": ["n+"]},
    "scope": "partial full-source implementation; inspect per-site diagnostics and coverage ledger",
}
NUMBERS = {
    "H": 1,
    "He": 2,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "Br": 35,
    "Kr": 36,
    "I": 53,
    "Xe": 54,
}
# Standard-weight range guard also detects isotope masses when an older input
# adapter omitted isotope metadata. Never alter a supplied mass. Values are the
# audited FRC element weights; 0.02 Da accommodates its documented rounding and
# current standard-weight conventions, not an isotope prediction policy.
MASSES = {
    "H": 1.00797,
    "He": 4.003,
    "C": 12.01115,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984,
    "Ne": 20.183,
    "Si": 28.086,
    "P": 30.9738,
    "S": 32.064,
    "Cl": 35.453,
    "Ar": 39.944,
    "Br": 79.909,
    "Kr": 83.8,
    "I": 126.9044,
    "Xe": 131.3,
}
# One evidence label per predicate branch. These are not arbitrary generic fallbacks.
RULE_TYPES = {
    "c",
    "c1",
    "c2",
    "c3",
    "c3h",
    "c3m",
    "c4h",
    "c4m",
    "co",
    "coh",
    "cp",
    "c5",
    "cs",
    "c=",
    "c=1",
    "c=2",
    "ct",
    "c_0",
    "c_1",
    "c_2",
    "cz",
    "c-",
    "hc",
    "ho",
    "hn2",
    "hn",
    "hs",
    "hsi",
    "o*",
    "hw",
    "oc",
    "oh",
    "o3e",
    "o4e",
    "op",
    "o_1",
    "o_2",
    "oo",
    "oz",
    "o-",
    "na",
    "n4",
    "n_2",
    "np",
    "nh",
    "nt",
    "n3m",
    "n4m",
    "s3e",
    "s4e",
    "sp",
    "sc",
    "sh",
    "s1",
    "si",
    "sio",
    "osi",
    "f",
    "cl",
    "br",
    "i",
    "he",
    "ne",
    "ar",
    "kr",
    "xe",
}


def environments(graph):
    atoms = {a["id"]: a for a in graph["sites"]}
    neighbors = {i: set() for i in atoms}
    bonds = {}
    for b in graph["bonds"]:
        i, j = b["sites"]
        neighbors[i].add(j)
        neighbors[j].add(i)
        bonds[frozenset((i, j))] = b
    env = {}
    for i, a in atoms.items():
        # Smallest cycle through i, independent of atom numbering/cycle basis.
        shortest = None
        for start in neighbors[i]:
            queue, distances = deque([start]), {start: 0}
            while queue:
                j = queue.popleft()
                if j != start and j in neighbors[i]:
                    n = distances[j] + 2
                    shortest = n if shortest is None else min(shortest, n)
                    break
                for k in neighbors[j] - {i}:
                    if k not in distances:
                        distances[k] = distances[j] + 1
                        queue.append(k)
        env[i] = {
            "element": a["element"],
            "degree": len(neighbors[i]),
            "hydrogens": sum(atoms[j]["element"] == "H" for j in neighbors[i]),
            "neighbor_elements": dict(
                sorted(Counter(atoms[j]["element"] for j in neighbors[i]).items())
            ),
            "neighbors": sorted(neighbors[i]),
            "smallest_ring": shortest,
            "bond_orders": sorted(
                bonds[frozenset((i, j))]["order"] for j in neighbors[i]
            ),
            "aromatic": a["metadata"]["aromatic"],
        }
    return atoms, neighbors, bonds, env


def recognize(graph, *, elemental_halogens=False, defer_components=False):
    atoms, adj, bonds, env = environments(graph)
    answers, issues = {}, []
    for bond in graph["bonds"]:
        if bond["aromatic"] != (bond["order"] == 1.5) or (
            bond["aromatic"]
            and not all(atoms[i]["metadata"]["aromatic"] for i in bond["sites"])
        ):
            issues.append(
                {
                    "sites": bond["sites"],
                    "reason": "inconsistent_aromatic_bond_semantics",
                }
            )

    def element(i):
        return atoms[i]["element"]

    def order(i, j):
        return bonds[frozenset((i, j))]["order"]

    def carbonyl(i):
        return element(i) == "C" and any(
            element(j) == "O" and order(i, j) == 2 for j in adj[i]
        )

    def carboxylate(i):
        return (
            element(i) == "C"
            and len(adj[i]) == 3
            and sum(element(j) == "O" and len(adj[j]) == 1 for j in adj[i]) == 2
            and sum(atoms[j]["formal_charge"] for j in adj[i]) == -1
            and atoms[i]["formal_charge"] == 0
            and sorted(order(i, j) for j in adj[i]) == [1, 1, 2]
            and all(
                (
                    element(j) == "O"
                    and atoms[j]["formal_charge"] == (-1 if order(i, j) == 1 else 0)
                )
                or (element(j) in ("C", "H") and order(i, j) == 1)
                for j in adj[i]
            )
        )

    for i, a in atoms.items():
        e = env[i]
        el = element(i)
        d = e["degree"]
        h = e["hydrogens"]
        ring = e["smallest_ring"]
        ns = adj[i]
        elements = e["neighbor_elements"]
        orders = e["bond_orders"]
        label = None
        reason = None
        if el not in NUMBERS or a["atomic_number"] != NUMBERS.get(el):
            reason = "element_or_atomic_number_outside_implemented_rules"
        elif (
            a["metadata"]["isotope"]
            or a["metadata"]["radical_electrons"]
            or abs(a["mass"] - MASSES[el]) > 0.02
        ):
            reason = "isotope_radical_or_nonstandard_mass_rule_unresolved"
        elif a["formal_charge"] and not (
            el == "O"
            and a["formal_charge"] == -1
            and d == 1
            and carboxylate(next(iter(ns)))
            or el == "N"
            and a["formal_charge"] == 1
            and orders == [1, 1, 1, 1]
            and all(element(j) in ("C", "H") for j in ns)
        ):
            reason = "charge_state_rule_unresolved"
        elif (
            el in ("He", "Ne", "Ar", "Kr", "Xe")
            and d == 0
            or el in ("F", "Cl", "Br", "I")
            and d == 1
            and orders == [1]
            and (
                element(next(iter(ns))) == "C"
                or (elemental_halogens and element(next(iter(ns))) == el)
            )
        ):
            label = el.lower()
        elif el == "C":
            if carboxylate(i):
                label = "c-"
            elif (
                e["aromatic"]
                and d == 3
                and ring in (5, 6)
                and orders in ([1, 1.5, 1.5], [1.5, 1.5, 1.5])
            ):
                label = (
                    "cs"
                    if ring == 5 and elements.get("S")
                    else "c5"
                    if ring == 5
                    else "cp"
                )
            elif (
                not e["aromatic"]
                and orders == [1, 3]
                and all(element(j) in ("C", "H", "N") for j in ns)
            ):
                label = "ct"
            elif not e["aromatic"] and orders == [1, 1, 2]:
                doubles = [j for j in ns if order(i, j) == 2]
                other = [j for j in ns if order(i, j) == 1]
                if element(doubles[0]) == "O" and len(adj[doubles[0]]) == 1:
                    hetero = [element(j) for j in other if element(j) in ("O", "N")]
                    if all(element(j) in ("C", "H") for j in other):
                        label = "c_0"
                    elif len(hetero) == 1:
                        label = "c_1"
                    elif len(hetero) == 2:
                        label = "cz" if hetero == ["O", "O"] else "c_2"
                elif element(doubles[0]) == "C" and all(
                    element(j) in ("C", "H", "F", "Cl", "Br", "I") for j in ns
                ):
                    terminal = sum(len(adj[j]) == 1 for j in ns) == 2
                    label = (
                        "c="
                        if terminal
                        else "c=1"
                        if env[doubles[0]]["hydrogens"] == 2
                        else "c=2"
                    )
            elif not e["aromatic"] and orders == [1, 1, 1, 1]:
                if ring in (3, 4):
                    label = f"c{ring}" + ("h" if h else "m")
                elif (
                    elements.get("O", 0) == 2
                    and all(
                        len(adj[j]) == 2 and all(element(k) == "C" for k in adj[j])
                        for j in ns
                        if element(j) == "O"
                    )
                    and h in (0, 1)
                ):
                    label = "coh" if h else "co"
                elif elements.get("O", 0) > 1:
                    reason = "multiple_oxygen_carbon_rule_unresolved"
                elif elements.get("N") and any(carbonyl(j) for j in ns):
                    reason = "amino_acid_alpha_rule_ambiguous"
                elif all(
                    element(j) in ("C", "H", "O", "N", "S", "Si", "F", "Cl", "Br", "I")
                    for j in ns
                ):
                    label = f"c{h}" if h in (1, 2, 3) else "c"
        elif el == "O":
            if d == 1 and carboxylate(next(iter(ns))):
                label = "o-"
            elif orders == [2] and element(next(iter(ns))) == "C":
                parent = next(iter(ns))
                label = (
                    "oo" if env[parent]["neighbor_elements"].get("O") == 3 else "o_1"
                )
            elif e["aromatic"] and d == 2 and ring == 5 and orders == [1.5, 1.5]:
                label = "op"
            elif not e["aromatic"] and orders == [1, 1]:
                if h == 2:
                    label = "o*"
                elif ring in (3, 4) and elements.get("C") == 2:
                    label = f"o{ring}e"
                elif elements.get("Si") == 2:
                    label = "osi"
                elif (
                    h == 1
                    and elements.get("C") == 1
                    and not any(carbonyl(j) for j in ns)
                ):
                    label = "oh"
                elif elements.get("C") == 2:
                    parents = [j for j in ns if carbonyl(j)]
                    if len(parents) == 0:
                        label = "oc"
                    elif len(parents) == 1:
                        label = (
                            "oz"
                            if env[parents[0]]["neighbor_elements"].get("O") == 3
                            else "o_2"
                        )
        elif el == "N" and not e["aromatic"]:
            if a["formal_charge"] == 1 and orders == [1, 1, 1, 1]:
                label = "n4"
            elif orders == [3] and element(next(iter(ns))) == "C":
                label = "nt"
            elif orders == [1, 1, 1] and all(element(j) in ("C", "H") for j in ns):
                if ring in (3, 4):
                    label = f"n{ring}m"
                elif any(carbonyl(j) for j in ns):
                    parents = [j for j in ns if carbonyl(j)]
                    if (
                        len(parents) == 1
                        and env[parents[0]]["neighbor_elements"].get("O") == 2
                    ):
                        label = "n_2"
                    else:
                        reason = "amide_n_vs_n_2_rule_requires_independent_fixture"
                elif any(atoms[j]["metadata"]["aromatic"] for j in ns):
                    reason = "aromatic_amine_nn_nb_choice_unresolved"
                else:
                    label = "na"
        elif el == "N" and e["aromatic"] and ring in (5, 6):
            if orders == [1.5, 1.5]:
                label = "np"
            elif orders == [1, 1.5, 1.5] and h == 1:
                label = "nh"
        elif el == "S":
            if e["aromatic"] and ring == 5 and orders == [1.5, 1.5]:
                label = "sp"
            elif orders == [1, 1]:
                if ring in (3, 4):
                    label = f"s{ring}e"
                elif elements.get("S") == 1 and elements.get("C") == 1:
                    label = "s1"
                elif elements.get("C") == 2:
                    label = "sc"
                elif elements.get("C") == 1 and h == 1:
                    label = "sh"
        elif el == "Si" and orders == [1, 1, 1, 1]:
            if all(element(j) in ("C", "H", "O") for j in ns):
                label = "sio" if elements.get("O") else "si"
        elif el == "H" and orders == [1]:
            parent = next(iter(ns))
            p = element(parent)
            if p == "C":
                label = "hc"
            elif p == "O" and env[parent]["hydrogens"] == 2:
                label = "hw"
            elif (
                p == "O"
                and env[parent]["neighbor_elements"].get("C") == 1
                and not any(carbonyl(j) for j in adj[parent])
            ):
                label = "ho"
            elif p == "N" and (
                atoms[parent]["formal_charge"] == 1 or env[parent]["aromatic"]
            ):
                label = "hn"
            elif p == "N" and env[parent]["bond_orders"] == [1, 1, 1]:
                label = "hn" if any(carbonyl(j) for j in adj[parent]) else "hn2"
            elif p == "S":
                label = "hs"
            elif p == "Si":
                label = "hsi"
        if label is None:
            issues.append(
                {
                    "sites": [i],
                    "reason": reason or "chemical_rule_ambiguous",
                    "environment": e,
                }
            )
        else:
            answers[i] = label
    if defer_components:
        return answers, env, issues
    # No successful component can conceal an unresolved connected neighborhood.
    blocked = {
        i
        for component in _components(adj)
        if any(set(x["sites"]) & set(component) for x in issues)
        for i in component
    }
    return {i: t for i, t in answers.items() if i not in blocked}, env, issues


def typing_data(graph, source, supplied=None, provenance=None, *, version=1):
    require(
        type(version) is int and version in (1, 2, 3, 4),
        "Unsupported graph profile version",
    )
    profile = deepcopy(PROFILE)
    if version == 2:
        profile.update(
            name="island_pcff_source_graph_v2",
            implementation="verified_local_environments_elemental_halogens_v2",
            elemental_halogens="neutral closed-shell homonuclear single bond; source valence-one label, native increment and automatic bond rows",
        )
    if version == 3:
        from .domains import EVIDENCE
        from .domains import PROFILE_NAME as DOMAIN_PROFILE

        profile.update(
            name=DOMAIN_PROFILE,
            implementation="audited_specific_domains_v1",
            domain_evidence=EVIDENCE,
        )
    if version == 4:
        from .organic_domains import EVIDENCE as ORGANIC_EVIDENCE

        profile.update(
            name="island_pcff_source_graph_v4",
            implementation="audited_specific_domains_v2",
            domain_evidence=ORGANIC_EVIDENCE,
        )
    source.require_assignment()
    require(
        source.identity["sha256"] == PROFILE["source_sha256"],
        "Expanded profile/source mismatch",
    )
    graph_system(graph)
    require(
        graph["representation"] == "atomistic" and not graph["has_box"],
        "Finite atomistic graph required",
    )
    if version == 4:
        from .organic_domains import recognize_organic

        automatic, env, diagnostics = recognize_organic(graph)
    elif version == 3:
        from .domains import recognize_domains

        automatic, env, diagnostics = recognize_domains(graph)
    else:
        automatic, env, diagnostics = recognize(graph, elemental_halogens=version == 2)
    assignments = automatic
    if supplied is not None:
        require(
            type(supplied) is dict
            and all(type(k) is int and type(v) is str for k, v in supplied.items()),
            "Exact integer explicit type mapping required",
        )
        require(set(supplied) == set(env), "Explicit type coverage mismatch")
        require(
            type(provenance) is str and provenance.strip(),
            "Explicit typing provenance required",
        )
        require(
            not diagnostics,
            "Explicit chemistry constraints unresolved; inspect automatic diagnostic",
        )
        # Explicit labels do not bypass chemical rules. Aliases need separately audited rules.
        require(
            all(
                label == automatic[i]
                or label in PROFILE["explicit_aliases"].get(automatic[i], [])
                for i, label in supplied.items()
            ),
            "Explicit labels contradict implemented chemical rules",
        )
        assignments = dict(supplied)
    atom_rows = records(source.inventory, "atom_types")
    entries = {}
    for i, e in env.items():
        label = assignments.get(i)
        row = (
            select([r for r in atom_rows if r["data"]["type"] == label])
            if label
            else None
        )
        if label:
            require(
                row is not None
                and (
                    row["record"]["data"]["element"] == e["element"]
                    or (version in (3, 4) and label == "dw" and e["element"] == "H")
                )
                and row["record"]["data"]["connections"] == e["degree"],
                f"Source type/chemical environment conflict at {i}: {label}",
            )
        entries[i] = {
            "type": label,
            "environment": e,
            "status": "typed" if label else "unresolved",
            "rule": f"source_graph_v1:{label}" if label else None,
            "source_atom_record": row,
        }
    return {
        "schema": TYPING_SCHEMA
        if version == 1
        else f"island_pcff_source_typing_v{version}",
        "source": source.identity,
        "graph": graph,
        "graph_identity": identity(graph),
        "profile": profile,
        "profile_identity": identity(profile),
        "origin": "explicit_checked" if supplied is not None else "automatic_graph",
        "explicit_types": supplied,
        "explicit_provenance": provenance,
        "assignments": assignments,
        "entries": entries,
        "diagnostics": diagnostics,
        "coverage": {
            "total": len(env),
            "typed": len(assignments),
            "complete": not diagnostics and len(assignments) == len(env),
        },
        **FLAGS,
    }


def charge_data(typing, source, *, resolution_policy=None):
    if resolution_policy is not None:
        from .fallbacks import validate_policy

        validate_policy(resolution_policy)
    require(
        typing["schema"]
        not in (
            "island_pcff_source_typing_v2",
            "island_pcff_source_typing_v3",
            "island_pcff_source_typing_v4",
        )
        or resolution_policy is not None,
        "Graph v2 requires explicit charge resolution policy",
    )
    require(typing["coverage"]["complete"], "Incomplete source graph typing")
    from .fallbacks import DOMAIN_POLICY

    detailed = resolution_policy == DOMAIN_POLICY
    inv = source.inventory
    incs = records(inv, "bond_increments")
    eqs = records(inv, "equivalence")
    charges = {i: [] for i in typing["assignments"]}
    contributions = []
    diagnostics = []
    for bond in typing["graph"]["bonds"]:
        a, b = bond["sites"]
        labels = [typing["assignments"][s] for s in (a, b)]
        match = bond_selection(labels, incs)
        resolved = labels
        evidence = []
        path = "direct"
        searches = [{"path": path, "types": list(labels), "matched": match is not None}]
        if match is None:
            evidence = [
                select([r for r in eqs if r["data"]["type"] == t]) for t in labels
            ]
            if all(evidence):
                resolved = [r["record"]["data"]["families"]["bond"] for r in evidence]
                match = bond_selection(resolved, incs)
                path = "equivalence.bond"
            searches.append(
                {
                    "path": "equivalence.bond",
                    "types": list(resolved),
                    "matched": match is not None,
                    "equivalence_records": evidence,
                }
            )
        if match is None and resolution_policy is not None:
            evidence = [
                select(
                    [
                        r
                        for r in records(inv, "auto_equivalence")
                        if r["data"]["type"] == t
                    ]
                )
                for t in labels
            ]
            if all(evidence):
                resolved = [
                    r["record"]["data"]["families"]["bond_increment"] for r in evidence
                ]
                match = bond_selection(resolved, incs)
                path = "auto_equivalence.bond_increment"
            searches.append(
                {
                    "path": "auto_equivalence.bond_increment",
                    "types": list(resolved),
                    "matched": match is not None,
                    "equivalence_records": evidence,
                }
            )
        if match is None:
            diagnostics.append(
                {
                    "sites": [a, b],
                    "types": labels,
                    "reason": "source_parameter_missing",
                    **({"searches": searches} if detailed else {}),
                    "detail": "No direct/ordinary increment. Automatic fallback not independently established."
                    if resolution_policy is None
                    else "No direct, ordinary bond, or automatic bond_increment row",
                    **(
                        {
                            "automatic_resolved_types": resolved,
                            "automatic_evidence": evidence,
                        }
                        if resolution_policy
                        else {}
                    ),
                }
            )
            continue
        for s, q in zip((a, b), match["increments"]):
            charges[s].append(q)
        contributions.append(
            {
                "sites": [a, b],
                "supplied_types": labels,
                "resolved_types": resolved,
                "path": path,
                "equivalence_records": evidence,
                **({"searches": searches} if detailed else {}),
                **match,
            }
        )
    totals = {s: fsum(q) for s, q in charges.items()}
    atoms, adj, _, _ = environments(typing["graph"])
    components = []
    for c in _components(adj):
        formal = sum(atoms[s]["formal_charge"] for s in c)
        total = fsum(totals[s] for s in c)
        components.append({"sites": c, "formal_charge": formal, "charge": total})
        if abs(total - formal) > PIN["tolerance_e"]:
            diagnostics.append(
                {
                    "sites": c,
                    "reason": "component_charge_mismatch",
                    "charge": total,
                    "formal_charge": formal,
                }
            )
    native = {
        "schema": CHARGE_SCHEMA,
        "typing_identity": identity(typing),
        "source": source.identity,
        "policy": resolution_policy or PROFILE["charge_policy"],
        "base_charge": 0.0,
        **(
            {
                "base_charges": dict.fromkeys(atoms, 0.0),
                "site_formal_charges": {
                    i: a["formal_charge"] for i, a in atoms.items()
                },
            }
            if detailed
            else {}
        ),
        "tolerance_e": PIN["tolerance_e"],
        "contributions": contributions,
        "partial_charges": totals,
        "components": components,
        "total_charge": fsum(totals.values()),
        "formal_charge": sum(a["formal_charge"] for a in atoms.values()),
        "unit": "elementary_charge",
        "complete": not diagnostics,
        "diagnostics": diagnostics,
        **FLAGS,
    }
    return {
        **({"resolution_policy": resolution_policy} if resolution_policy else {}),
        "schema": CHARGE_SCHEMA
        if resolution_policy is None
        else "island_pcff_source_charges_v2",
        "automatic_typing": typing,
        "automatic_typing_identity": identity(typing),
        "bridge": "source_graph_typing_native_increments_v1",
        "native_charge_record": native,
        "native_charge_identity": identity(native),
        **FLAGS,
    }


@boundary
def assign_pcff_source_types(
    system, source, types, *, provenance, profile=PROFILE_NAME
):
    """Explicit labels checked by the same chemical contract as automatic typing."""
    from .automatic import PCFFAutomaticTypingResult

    require(
        profile
        in (
            PROFILE_NAME,
            "island_pcff_source_graph_v2",
            "island_pcff_source_graph_v3",
            "island_pcff_source_graph_v4",
        ),
        "Unsupported explicit profile",
    )
    result = PCFFAutomaticTypingResult(
        pack(
            typing_data(
                chemical_graph(system),
                source,
                types,
                provenance,
                version=int(profile[-1]),
            )
        ),
        source,
    )
    result.validate_integrity(system)
    return result
