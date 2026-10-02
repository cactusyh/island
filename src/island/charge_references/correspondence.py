"""Strict, dependency-free correspondence for two explicit builder definitions."""

from itertools import pairwise
from math import isclose, isfinite

from island.exceptions import ChargeReferenceError

DEFINITIONS = {"[*:1]CC[*:2]": ("C", "C"), "[*:1]CCO[*:2]": ("C", "C", "O")}


def require(condition, message):
    if not condition:
        raise ChargeReferenceError(message)


def _correspondence(system):
    system.validate()
    p = system.metadata["polymer"]
    definition = p["source_psmiles"]
    require(definition in DEFINITIONS, "Unsupported mapped repeat definition")
    elements = DEFINITIONS[definition]
    dp = p["degree_of_polymerization"]
    require(type(dp) is int and dp >= 3 and dp % 2 == 1, "Require odd DP >=3")
    require(
        system.representation == "atomistic" and system.box is None,
        "Require nonperiodic atomistic system",
    )
    require(
        p["architecture"] == "linear" and p["polymer_type"] == "homopolymer",
        "Require linear homopolymer",
    )
    require(
        p["number_of_repeat_units"] == dp
        and p["number_of_inter_repeat_unit_bonds"] == dp - 1,
        "Repeat counts disagree",
    )
    sequence = p["sequence"]
    require(len(sequence) == dp and len(set(sequence)) == 1, "Copolymers unsupported")
    identity = sequence[0]
    require(
        p["repeat_unit_definitions"] == {identity: definition}
        and p["composition_counts"] == {identity: dp},
        "Repeat definitions disagree",
    )
    require(
        not any(k in p for k in ("tacticity", "stereochemical_sequence")),
        "Stereochemical cases unsupported",
    )
    sites = system.topology.sites
    neighbors = {i: set() for i in sites}
    actual_edges = set()
    for b in system.topology.bonds.values():
        require(
            b.order == 1 and not b.aromatic, "Only single nonaromatic bonds supported"
        )
        edge = frozenset((b.site1, b.site2))
        require(edge not in actual_edges, "Duplicate bond")
        actual_edges.add(edge)
        neighbors[b.site1].add(b.site2)
        neighbors[b.site2].add(b.site1)
    heavy = {}
    for i, s in sites.items():
        m = s.metadata
        require(
            type(i) is int and i == s.id and s.formal_charge == 0,
            "Invalid IDs or nonneutral atomic formal charge",
        )
        require(s.element in {"C", "O", "H"}, "Unsupported element")
        require(
            s.atomic_number == {"C": 6, "O": 8, "H": 1}[s.element]
            and not m.get("aromatic", False),
            "Atomic identity/aromaticity disagrees",
        )
        require(
            isfinite(s.mass)
            and isclose(
                s.mass, {"C": 12.011, "O": 15.999, "H": 1.008}[s.element], abs_tol=1e-6
            ),
            "Isotopic or nonstandard masses unsupported",
        )
        require(
            not m.get("isotope")
            and not m.get("cip_label")
            and m.get("chiral_tag", "CHI_UNSPECIFIED") == "CHI_UNSPECIFIED",
            "Isotope/stereochemistry unsupported",
        )
        require(
            m["chain_id"] == p["chain_id"] and m["repeat_unit_type"] == identity,
            "Site chain/repeat identity disagrees",
        )
        r = m["repeat_unit_index"]
        require(type(r) is int and 0 <= r < dp, "Invalid repeat index")
        if s.element != "H":
            j = m["source_repeat_atom_index"]
            require(
                type(j) is int
                and 1 <= j <= len(elements)
                and s.element == elements[j - 1],
                "Invalid source repeat atom identity",
            )
            require(
                (r, j) not in heavy and not m.get("generated_hydrogen"),
                "Ambiguous heavy identity",
            )
            heavy[r, j] = i
        else:
            require(
                m.get("generated_hydrogen") is True
                and "source_repeat_atom_index" not in m,
                "Hydrogen provenance unsupported",
            )
    require(
        set(heavy) == {(r, j) for r in range(dp) for j in range(1, len(elements) + 1)},
        "Incomplete heavy inventory",
    )
    require(
        p["head_site_id"] == heavy[0, 1]
        and p["tail_site_id"] == heavy[dp - 1, len(elements)],
        "Head/tail orientation disagrees",
    )
    path = [heavy[r, j] for r in range(dp) for j in range(1, len(elements) + 1)]
    expected_edges = {frozenset((a, b)) for a, b in pairwise(path)}
    repeats = []
    for r in range(dp):
        groups = []
        for j in range(1, len(elements) + 1):
            i = heavy[r, j]
            hs = sorted(n for n in neighbors[i] if sites[n].element == "H")
            heavy_neighbors = neighbors[i] - set(hs)
            require(
                len(hs) == {"C": 4, "O": 2}[sites[i].element] - len(heavy_neighbors),
                "Hydrogen count/valence mismatch",
            )
            for h in hs:
                require(
                    neighbors[h] == {i} and sites[h].metadata["repeat_unit_index"] == r,
                    "Hydrogen parent membership disagrees",
                )
                expected_edges.add(frozenset((i, h)))
            groups.append(
                {
                    "source_repeat_atom_index": j,
                    "element": sites[i].element,
                    "site_id": i,
                    "hydrogen_ids": hs,
                    "heavy_environment": sorted(
                        (
                            sites[n].metadata["repeat_unit_index"] - r,
                            sites[n].metadata["source_repeat_atom_index"],
                            sites[n].element,
                        )
                        for n in heavy_neighbors
                    ),
                }
            )
        repeats.append(
            {
                "repeat_index": r,
                "role": "head" if r == 0 else "tail" if r == dp - 1 else "interior",
                "central": r == dp // 2,
                "groups": groups,
            }
        )
    require(
        actual_edges == expected_edges,
        "Connectivity inconsistent with oriented repeat definition",
    )
    require(
        sum(1 + len(g["hydrogen_ids"]) for r in repeats for g in r["groups"])
        == len(sites),
        "Unmapped sites",
    )
    return {
        "compatible": True,
        "definition": definition,
        "dp": dp,
        "end_groups": "hydrogen_terminated",
        "repeats": repeats,
        "diagnostics": [],
    }


def repeat_correspondence(system):
    """Return a verified mapping or structured diagnostics; never guess atom matches."""
    try:
        return _correspondence(system)
    except Exception as error:  # noqa: BLE001 -- structured mapping diagnostics
        return {
            "compatible": False,
            "diagnostics": [{"reason": str(error), "type": type(error).__name__}],
        }
