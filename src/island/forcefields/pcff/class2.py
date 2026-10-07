"""Source-native Class II assignment, not an executable potential.

This layer leaves the historical FRC parser and signed H1/H2 contracts intact.
Historical exact/ordinary selection remains the default. Explicit J2 policy adds
source-supported lower-order terms; missing couplings still block completeness.
"""

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal
from itertools import combinations, permutations
from math import isfinite, pi
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .automatic import PCFFAutomaticChargeResult, chemical_graph
from .charges import identity
from .source import FLAGS, boundary, records, require, select
from .validation_cache import record_identity, validated_record

SCHEMA = "island_pcff_class2_assignment_v1"
POLICY = "exact_then_ordinary_family_equivalence_highest_version_v1"
# arity, equivalence family, source fields and their native dimensions.
# Angular displacements are radians; only equilibrium angles/phases convert.
FAMILIES = {
    "quartic_bond": (2, "bond", ["angstrom", "E/L^2", "E/L^3", "E/L^4"]),
    "quartic_angle": (3, "angle", ["degree", "E/A^2", "E/A^3", "E/A^4"]),
    "torsion_3": (4, "torsion", ["E", "degree"] * 3),
    "wilson_out_of_plane": (4, "out_of_plane", ["E/A^2", "degree"]),
    "nonbond(9-6)": (1, "nonbond", ["angstrom", "E"]),
    "bond-bond": (3, "angle", ["E/L^2"]),
    "bond-angle": (3, "angle", ["E/L/A"] * 2),
    "bond-bond_1_3": (4, "torsion", ["E/L^2"]),
    "end_bond-torsion_3": (4, "torsion", ["E/L"] * 6),
    "middle_bond-torsion_3": (4, "torsion", ["E/L"] * 3),
    "angle-torsion_3": (4, "torsion", ["E/A"] * 6),
    "angle-angle-torsion_1": (4, "torsion", ["E/angular_power_unresolved"]),
    "angle-angle": (4, "out_of_plane", ["E/A^2"]),
}
FIELDS = {
    "quartic_bond": ["r0", "K2", "K3", "K4"],
    "quartic_angle": ["theta0", "K2", "K3", "K4"],
    "torsion_3": ["V1", "phase1", "V2", "phase2", "V3", "phase3"],
    "wilson_out_of_plane": ["Kchi", "chi0"],
    "nonbond(9-6)": ["r_min", "epsilon"],
    "bond-bond": ["Kbb"],
    "bond-angle": ["Kleft", "Kright"],
    "bond-bond_1_3": ["Kbb13"],
    "end_bond-torsion_3": [
        "Fleft1",
        "Fleft2",
        "Fleft3",
        "Fright1",
        "Fright2",
        "Fright3",
    ],
    "middle_bond-torsion_3": ["F1", "F2", "F3"],
    "angle-torsion_3": ["Fleft1", "Fleft2", "Fleft3", "Fright1", "Fright2", "Fright3"],
    "angle-angle-torsion_1": ["Kaat"],
    "angle-angle": ["Kaa"],
}

CONVENTIONS = {
    "profile": "island_pcff_class2_source_records_v1",
    "selection": POLICY,
    "energy": "kJ/mol",
    "length": "angstrom",
    "angle": "radian",
    "source_energy": "kcal/mol",
    "kcal_to_kj": 4.184,
    "polynomial_prefactor": "no one-half; source coefficients retained",
    "angular_coefficients": "radian displacements; no degree scaling",
    "distance_9_6": "minimum-energy distance r*, not LJ12-6 sigma",
    "mixing": {
        "distance": "((ri**6+rj**6)/2)**(1/6)",
        "epsilon": "2*sqrt(ei*ej)*ri**3*rj**3/(ri**6+rj**6)",
        "policy": "source sixth-power",
    },
    "angle_angle_roles": "I-J-K and K-J-L; J center, K shared arm; swap I/L only",
    "torsion_roles": "I-J-K-L consecutive bonded sites; reverse swaps end blocks",
    "wilson_applicability": "three-connected center only; none in saturated explicit-H CHO",
    "angle_angle_multiplicity": "three angle pairs per neighbor triple, including degree-four centers",
    "unresolved_physics": [
        "FRC torsion_3 plus-cos comment versus pinned LAMMPS minus-cos implementation; no evaluator conversion claimed",
        "FRC angle-angle-torsion_1 linear-Phi comment versus pinned LAMMPS cos(Phi) implementation",
        "source-specific LJ and Coulomb exclusions/scaling not declared by FRC; no inherited defaults",
        "empty torsion-torsion_1 section: applicability/absence not established",
    ],
    "reference_revision": "e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
}


def conversion(unit):
    return pi / 180 if unit == "degree" else 4.184 if unit.startswith("E") else 1.0


@boundary
def inspect_pcff_class2(source):
    """An additive semantic catalog; raw unsupported sections remain in source.inventory."""
    inventory = source.inventory
    catalog = {}
    sections = []
    for section in inventory["sections"]:
        name = section["name"]
        interpreted = name in FAMILIES and section["namespace"] == "cff91"
        sections.append(
            {
                "name": name,
                "namespace": section["namespace"],
                "line": section["line"],
                "interpretation": "class2_records"
                if interpreted
                else section["interpretation"],
            }
        )
        if not interpreted:
            continue
        arity, _, units = FAMILIES[name]
        for line in section["lines"]:
            text = line["raw"].strip()
            if not text or text.startswith(("!", ">", "@")):
                continue
            try:
                fields = text.split("!", 1)[0].split()
                version = Decimal(fields[0])
                require(
                    version.is_finite() and version >= 0 and fields[1].isdigit(),
                    "Invalid version/reference",
                )
                types = fields[2 : 2 + arity]
                values = [float(v) for v in fields[2 + arity :]]
                require(
                    len(types) == arity and all(isfinite(v) for v in values),
                    "Invalid types/numbers",
                )
                original = list(values)
                symmetric = name in (
                    "bond-angle",
                    "end_bond-torsion_3",
                    "angle-torsion_3",
                ) and len(values) * 2 == len(units)
                if symmetric:
                    values = values * 2
                require(len(values) == len(units), "Truncated/extra coefficients")
                rid = f"{name}:cff91:{line['line']}"
                catalog[rid] = {
                    "id": rid,
                    "line": line["line"],
                    "raw": line["raw"],
                    "section": name,
                    "namespace": "cff91",
                    "version": fields[0],
                    "reference": fields[1],
                    "types": types,
                    "original_values": original,
                    "source_units": list(units),
                    "fields": list(FIELDS[name]),
                    "functional_form_status": "unresolved_source_vs_reference"
                    if name in ("torsion_3", "angle-angle-torsion_1")
                    else "source_records_only",
                    "expanded_values": values,
                    "expansion": "symmetric_half_repeat" if symmetric else "none",
                    "conversion_factors": [conversion(u) for u in units],
                    "normalized_values": [
                        v * conversion(u) for v, u in zip(values, units)
                    ],
                    "normalized_units": [
                        u.replace("E", "kJ/mol")
                        .replace("L", "angstrom")
                        .replace("A", "radian")
                        .replace("degree", "radian")
                        for u in units
                    ],
                }
            except Exception as error:
                from island.exceptions import PCFFError

                raise PCFFError(f"{name} line {line['line']}: {error}") from error
    return {
        "source": source.identity,
        "conventions": deepcopy(CONVENTIONS),
        "sections": sections,
        "records": catalog,
    }


def reverse_values(family, values):
    if family == "bond-angle":
        return values[::-1]
    if family in ("end_bond-torsion_3", "angle-torsion_3"):
        return values[3:] + values[:3]
    return values[:]


def resolve(family, supplied, catalog, equivalents):
    """Compare oriented coefficients, including repeated labels; never file-order ties."""
    eq_family = FAMILIES[family][1]
    paths = [("direct", supplied, [])]
    eq_types, evidence = [], []
    for label in supplied:
        eq = select([r for r in equivalents if r["data"]["type"] == label])
        if eq is None:
            eq_types = []
            break
        eq_types.append(eq["record"]["data"]["families"][eq_family])
        evidence.append(eq)
    if eq_types:
        paths.append(("ordinary_" + eq_family, eq_types, evidence))
    for path, types, eq_evidence in paths:
        candidates = []
        for row in catalog.values():
            if row["section"] != family:
                continue
            # AA center and shared arm stay fixed; never full-reverse an improper.
            reverse = (
                [types[3], types[1], types[2], types[0]]
                if family == "angle-angle"
                else types[::-1]
            )
            orientations = [
                ("forward", types),
                ("swapped_outer" if family == "angle-angle" else "reverse", reverse),
            ]
            if family == "wilson_out_of_plane":
                orientations = [
                    (
                        "arms:" + ",".join(map(str, p)),
                        [types[p[0]], types[1], types[p[1]], types[p[2]]],
                    )
                    for p in permutations((0, 2, 3))
                ]
            for orientation, query in orientations:
                if row["types"] == query:
                    values = row["normalized_values"]
                    if orientation == "reverse":
                        values = reverse_values(family, values)
                    candidates.append(
                        {
                            "record_id": row["id"],
                            "version": row["version"],
                            "orientation": orientation,
                            "values": values,
                        }
                    )
        if not candidates:
            continue
        latest = max(Decimal(r["version"]) for r in candidates)
        chosen = [r for r in candidates if Decimal(r["version"]) == latest]
        common = {
            "supplied_types": supplied,
            "resolved_types": types,
            "path": path,
            "equivalence_evidence": eq_evidence,
            "candidates": candidates,
        }
        if any(r["values"] != chosen[0]["values"] for r in chosen):
            return {
                **common,
                "status": "ambiguous",
                "reason": "conflicting highest-version oriented coefficients",
            }
        # Equivalent numerical candidates retain all identities; deterministic representative.
        chosen.sort(key=lambda r: (r["record_id"], r["orientation"]))
        return {
            **common,
            "status": "assigned",
            "selected": chosen,
            "normalized_values": chosen[0]["values"],
        }
    return {
        "status": "missing",
        "supplied_types": supplied,
        "resolved_types": eq_types,
        "reason": "no exact or ordinary family-equivalence row; no fallback",
        "candidates": [],
    }


def canonical(sites):
    return min(tuple(sites), tuple(reversed(sites)))


def interaction_requests(graph, *, expanded):
    """Authoritative ordered requests, shared by assignment and failure inspection.

    This enumerates structure only. It cannot establish charge or parameter
    validity and must never be used to turn a failed import into a valid model.
    """
    neighbors = {s["id"]: set() for s in graph["sites"]}
    for bond in graph["bonds"]:
        a, b = bond["sites"]
        neighbors[a].add(b)
        neighbors[b].add(a)
    bonds = sorted(tuple(b["sites"]) for b in graph["bonds"])
    angles = sorted(
        (a, j, b) for j in neighbors for a, b in combinations(sorted(neighbors[j]), 2)
    )
    torsions = sorted(
        {
            canonical((a, b, c, d))
            for b, c in bonds
            for a in neighbors[b] - {c}
            for d in neighbors[c] - {b}
            if a != d
        }
    )
    requests = []

    def add(family, sites, dependencies=(), reason=None):
        requests.append((family, tuple(sites), tuple(dependencies), reason))

    for i in sorted(neighbors):
        add("nonbond(9-6)", [i])
    for bond in bonds:
        add("quartic_bond", bond)
    for angle in angles:
        add("quartic_angle", angle)
    for a, b, c in angles:
        lengths = [("quartic_bond", (a, b)), ("quartic_bond", (b, c))]
        add("bond-bond", (a, b, c), lengths)
        add("bond-angle", (a, b, c), lengths + [("quartic_angle", (a, b, c))])
    for a, b, c, d in torsions:
        ends = [("quartic_bond", (a, b)), ("quartic_bond", (c, d))]
        bends = [("quartic_angle", (a, b, c)), ("quartic_angle", (b, c, d))]
        for family, deps in (
            ("torsion_3", []),
            ("bond-bond_1_3", ends),
            ("end_bond-torsion_3", ends),
            ("middle_bond-torsion_3", [("quartic_bond", (b, c))]),
            ("angle-torsion_3", bends),
            ("angle-angle-torsion_1", bends),
        ):
            add(family, (a, b, c, d), deps)
    for j in sorted(neighbors):
        for arms in combinations(sorted(neighbors[j]), 3):
            a, k, l = arms
            add(
                "wilson_out_of_plane",
                (a, j, k, l),
                reason=(
                    None
                    if expanded and len(neighbors[j]) == 3
                    else "degree-four saturated center; not a three-connected Wilson center"
                ),
            )
            for shared in arms:
                outer = sorted(set(arms) - {shared})
                i, l = outer
                add(
                    "angle-angle",
                    (i, j, shared, l),
                    [
                        ("quartic_angle", (i, j, shared)),
                        ("quartic_angle", (shared, j, l)),
                    ],
                )
    return requests, {
        "bonds": [list(b) for b in bonds],
        "angles": [list(a) for a in angles],
        "proper_torsions": [list(t) for t in torsions],
    }


def derive(charge, source, resolution_policy=None, *, resolution_cache=None):
    from .expanded import PROFILE_NAME

    auto = charge["automatic_typing"]
    expanded = auto["profile"]["name"] in (
        PROFILE_NAME,
        "island_pcff_source_graph_v2",
        "island_pcff_source_graph_v3",
        "island_pcff_source_graph_v4",
        "island_pcff_source_graph_v5",
    )
    if resolution_policy is not None:
        from .fallbacks import validate_policy

        validate_policy(resolution_policy)
        require(expanded, "Fallback requires expanded graph profile")
    require(
        charge.get("resolution_policy") == resolution_policy,
        "Charge/parameter resolution policy mismatch",
    )
    graph = auto["graph"]
    types = auto["assignments"]
    inspected = inspect_pcff_class2(source)
    if resolution_policy:
        from .catalog import inspect_pcff_full_source

        inspected = inspect_pcff_full_source(source)
        inspected["records"] = {r["id"]: r for r in inspected["records"]}
    catalog = inspected["records"]
    inventory = source.inventory
    equivalents = records(inventory, "equivalence")
    requests, inventories = interaction_requests(graph, expanded=expanded)
    assignments = []
    base = {}
    cache = deepcopy(resolution_cache) if resolution_cache is not None else {}

    def add(family, sites, dependencies=(), reason=None):
        sites = tuple(sites)
        key = f"{family}:" + ",".join(map(str, sites))
        labels = [types[i] for i in sites]
        lookup = (family, tuple(labels))
        if reason:
            found = {
                "status": "not_applicable",
                "reason": reason,
                "supplied_types": labels,
            }
        else:
            if lookup not in cache:
                if resolution_policy:
                    from .fallbacks import resolve as supplement

                    cache[lookup] = supplement(
                        family, labels, catalog, inventory, policy=resolution_policy
                    )
                else:
                    cache[lookup] = resolve(family, labels, catalog, equivalents)
            found = cache[lookup]
        deps = []
        for dep_family, dep_sites in dependencies:
            dep = base[(dep_family, canonical(dep_sites))]
            deps.append(
                {
                    "assignment_id": dep["id"],
                    "sites": list(dep_sites),
                    "status": dep["status"],
                    "equilibrium_value": dep.get("normalized_values", [None])[0],
                    "source_rows": [r["record_id"] for r in dep.get("selected", [])],
                }
            )
        entry = {
            "id": key,
            "family": found.get("selected_family", family)
            if found["status"] == "assigned"
            else family,
            **({"requested_family": family} if resolution_policy else {}),
            "sites": list(sites),
            **found,
            "dependencies": deps,
        }
        if entry["status"] == "assigned" and any(
            d["status"] != "assigned" for d in deps
        ):
            entry = {
                **entry,
                "status": "missing",
                "reason": "unresolved equilibrium dependency",
            }
        assignments.append(entry)
        if family in ("quartic_bond", "quartic_angle"):
            base[(family, canonical(sites))] = entry

    for family, sites, dependencies, reason in requests:
        add(family, sites, dependencies, reason)
    coverage = {
        f: dict(Counter(a["status"] for a in assignments if a["family"] == f))
        for f in (
            list(FAMILIES)
            + (
                ["quadratic_bond", "quadratic_angle", "torsion_1"]
                if resolution_policy
                else []
            )
        )
    }
    return {
        "schema": "island_pcff_source_class2_assignment_v2"
        if resolution_policy
        else "island_pcff_source_class2_assignment_v1"
        if expanded
        else SCHEMA,
        **({"resolution_policy": resolution_policy} if resolution_policy else {}),
        "source": source.identity,
        "charge_record": charge,
        "charge_identity": identity(charge),
        "graph_identity": identity(graph),
        "typing_identity": charge["automatic_typing_identity"],
        "conventions": (
            {
                **deepcopy(CONVENTIONS),
                "profile": "island_pcff_expanded_source_records_v1",
                "wilson_applicability": "one unordered neighbor triple at each degree-three center; center fixed, six peripheral permutations",
            }
            if expanded
            else deepcopy(CONVENTIONS)
        ),
        "source_catalog": inspected,
        "inventories": inventories,
        "assignments": assignments,
        "coverage": coverage,
        "parameter_coverage_complete": all(
            a["status"] in ("assigned", "not_applicable") for a in assignments
        ),
        "physical_model_complete": False,
        "evaluator_implemented": False,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFClass2Result:
    """Owned data-only assignment. Completeness never implies an executable potential."""

    json_text: str
    source: object

    @boundary
    @validated_record
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        require(
            p["schema"]
            in (
                SCHEMA,
                "island_pcff_source_class2_assignment_v1",
                "island_pcff_source_class2_assignment_v2",
            ),
            "Unsupported Class II schema",
        )
        charge = PCFFAutomaticChargeResult(pack(p["charge_record"]), self.source)
        charge.validate_integrity(system)
        require(charge.complete, "Complete compatible native charges required")
        require(
            pack(p)
            == pack(derive(charge.payload, self.source, p.get("resolution_policy"))),
            "Contradictory Class II assignment",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return record_identity(self)


@boundary
def assign_pcff_parameters(system, typing, charges, *, resolution_policy=None):
    """Resolve native records; supplementation requires an explicit versioned policy."""
    typing.validate_integrity(system)
    require(
        type(charges) is PCFFAutomaticChargeResult,
        "Automatic native charge result required",
    )
    charges.validate_integrity(system)
    require(
        typing.complete and charges.complete, "Complete typing and charges required"
    )
    p = charges.payload
    require(
        p["automatic_typing_identity"] == typing.identity,
        "Typing/charge identity mismatch",
    )
    require(p["automatic_typing"]["graph"] == chemical_graph(system), "Graph mismatch")
    result = PCFFClass2Result(
        pack(derive(p, charges.source, resolution_policy)), charges.source
    )
    result.validate_integrity(system)
    return result


@boundary
def save_pcff_parameters(result, path):
    require(type(result) is PCFFClass2Result, "Expected Class II result")
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_parameters(path, source, *, system=None):
    result = PCFFClass2Result(Path(path).read_text(), source)
    result.validate_integrity(system)
    return result
