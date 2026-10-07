"""Explicit, force-field-independent crosslink planning on a final chemical graph.

The caller owns the reaction-site annotation and any preceding atom removal.
Plans select bonds only; applying one never changes atoms, charges or coordinates.
"""

import random
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass, field
from hashlib import sha256
from itertools import combinations
from math import floor, isfinite
from pathlib import Path

from island.exceptions import ValidationError
from island.graph.final import (
    FinalChemicalGraph,
    final_graph,
    pack,
    transformation,
    unpack,
)

PLAN_SCHEMA = "island_crosslink_network_plan_v1"
RULE_SCHEMA = "island_reactive_site_rule_v1"
SITE_SCHEMA = "island_reactive_site_v1"
STRATEGY = "seeded_random"
RULE_FIELDS = {
    "schema",
    "name",
    "metadata_key",
    "metadata_value",
    "allowed_elements",
    "maximum_crosslinks_per_site",
    "bond_order",
    "provenance",
    "evidence",
}
SITE_FIELDS = {
    "schema",
    "site_id",
    "role",
    "element",
    "component_identity",
    "current_bond_order_valence",
    "available_capacity",
    "rule_identity",
    "provenance",
    "evidence",
}


def _require(condition, message):
    if not condition:
        raise ValidationError(message)


def _digest(value):
    return sha256(pack(value).encode()).hexdigest()


def _text(value, label):
    _require(type(value) is str and bool(value.strip()), f"{label} required")


def _evidence(value):
    _require(
        type(value) in (tuple, list)
        and all(type(item) is str and item.strip() for item in value),
        "Evidence must be a sequence of nonempty references",
    )
    return tuple(value)


def _number(value, label, *, positive=False):
    _require(
        type(value) in (int, float) and isfinite(value),
        f"{label} must be finite and nonboolean",
    )
    if positive:
        _require(value > 0, f"{label} must be positive")


def _integer(value, label, *, minimum=0):
    _require(
        type(value) is int and value >= minimum,
        f"{label} must be an integer >= {minimum}",
    )


def _pair(a, b):
    _require(type(a) is int and type(b) is int, "Pair site IDs must be integers")
    _require(a != b, "Self-pair crosslink is forbidden")
    return tuple(sorted((a, b)))


def _rule_payload(rule):
    return {
        "schema": RULE_SCHEMA,
        "name": rule.name,
        "metadata_key": rule.metadata_key,
        "metadata_value": rule.metadata_value,
        "allowed_elements": list(rule.allowed_elements),
        "maximum_crosslinks_per_site": rule.maximum_crosslinks_per_site,
        "bond_order": rule.bond_order,
        "provenance": rule.provenance,
        "evidence": list(rule.evidence),
    }


@dataclass(frozen=True)
class ReactiveSiteRule:
    name: str
    metadata_key: str
    metadata_value: str | int | float | bool
    allowed_elements: tuple[str, ...]
    maximum_crosslinks_per_site: int
    bond_order: float
    provenance: str
    evidence: tuple[str, ...]

    def __post_init__(self):
        for value, label in (
            (self.name, "Rule name"),
            (self.metadata_key, "Metadata key"),
            (self.provenance, "Rule provenance"),
        ):
            _text(value, label)
        _require(
            not any(
                token in self.metadata_key.lower()
                for token in (
                    "force_field",
                    "forcefield",
                    "atom_type",
                    "source",
                    "charge",
                    "pcff",
                    "opls",
                    "gaff",
                    "frc",
                )
            ),
            "Reactive marker must not refer to force-field or source metadata",
        )
        _require(
            type(self.metadata_value) in (str, int, float, bool),
            "Metadata value must be a JSON scalar",
        )
        if type(self.metadata_value) is float:
            _number(self.metadata_value, "Metadata value")
        _require(
            type(self.allowed_elements) in (tuple, list),
            "Allowed elements must be a sequence",
        )
        elements = tuple(self.allowed_elements)
        _require(
            elements
            and all(type(x) is str and x for x in elements)
            and len(set(elements)) == len(elements),
            "Allowed elements must be unique nonempty symbols",
        )
        object.__setattr__(self, "allowed_elements", elements)
        _integer(self.maximum_crosslinks_per_site, "Maximum site usage", minimum=1)
        _number(self.bond_order, "Bond order", positive=True)
        _require(
            self.bond_order in (1, 2, 3),
            "Only integral single, double or triple new bonds are supported",
        )
        object.__setattr__(self, "evidence", _evidence(self.evidence))

    @property
    def identity(self):
        return _digest(_rule_payload(self))

    @property
    def payload(self):
        return deepcopy(_rule_payload(self))


def _rule_from_payload(payload):
    _require(
        type(payload) is dict and set(payload) == RULE_FIELDS,
        "Malformed reactive rule",
    )
    _require(payload["schema"] == RULE_SCHEMA, "Unsupported reactive rule schema")
    _require(type(payload["allowed_elements"]) is list, "Malformed allowed elements")
    _require(type(payload["evidence"]) is list, "Malformed rule evidence")
    return ReactiveSiteRule(
        payload["name"],
        payload["metadata_key"],
        payload["metadata_value"],
        tuple(payload["allowed_elements"]),
        payload["maximum_crosslinks_per_site"],
        payload["bond_order"],
        payload["provenance"],
        tuple(payload["evidence"]),
    )


def _site_payload(site):
    return {
        "schema": SITE_SCHEMA,
        "site_id": site.site_id,
        "role": site.role,
        "element": site.element,
        "component_identity": site.component_identity,
        "current_bond_order_valence": site.current_bond_order_valence,
        "available_capacity": site.available_capacity,
        "rule_identity": site.rule_identity,
        "provenance": site.provenance,
        "evidence": list(site.evidence),
    }


@dataclass(frozen=True)
class ReactiveSite:
    site_id: int
    role: str
    element: str
    component_identity: str
    current_bond_order_valence: float
    available_capacity: int
    rule_identity: str
    provenance: str
    evidence: tuple[str, ...]
    # Retained in memory so the suggested plan_crosslinks(graph, sites, ...)
    # signature can recover the explicit rule. It is deliberately excluded
    # from the JSON record and therefore cannot alter the site identity.
    rule: ReactiveSiteRule | None = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        _require(type(self.site_id) is int, "Site ID must be an integer")
        for value, label in (
            (self.role, "Site role"),
            (self.element, "Element"),
            (self.component_identity, "Component identity"),
            (self.rule_identity, "Rule identity"),
            (self.provenance, "Site provenance"),
        ):
            _text(value, label)
        _number(self.current_bond_order_valence, "Current valence")
        _integer(self.available_capacity, "Available capacity")
        object.__setattr__(self, "evidence", _evidence(self.evidence))

    @property
    def identity(self):
        return _digest(_site_payload(self))

    @property
    def payload(self):
        return deepcopy(_site_payload(self))


def _site_from_payload(payload):
    _require(
        type(payload) is dict and set(payload) == SITE_FIELDS,
        "Malformed reactive site",
    )
    _require(payload["schema"] == SITE_SCHEMA, "Unsupported reactive site schema")
    _require(type(payload["evidence"]) is list, "Malformed site evidence")
    return ReactiveSite(
        payload["site_id"],
        payload["role"],
        payload["element"],
        payload["component_identity"],
        payload["current_bond_order_valence"],
        payload["available_capacity"],
        payload["rule_identity"],
        payload["provenance"],
        tuple(payload["evidence"]),
    )


def _graph_sites(graph):
    _require(isinstance(graph, FinalChemicalGraph), "FinalChemicalGraph required")
    graph.validate_integrity()
    payload = graph.payload
    _require(
        payload["periodic"] is None, "Periodic-image crosslink policy is unsupported"
    )
    rows = payload["sites"]
    ids = [row["id"] for row in rows]
    _require(
        all(type(i) is int for i in ids) and ids == sorted(set(ids)),
        "Stable unique integer atom IDs required",
    )
    _require(
        all(type(row["metadata"]) is dict for row in rows), "Malformed atom metadata"
    )
    bonds = payload["bonds"]
    _require(
        all(
            type(b["order"]) in (int, float) and isfinite(b["order"]) and b["order"] > 0
            for b in bonds
        ),
        "Invalid bond order",
    )
    adjacency = {i: set() for i in ids}
    for bond in bonds:
        a, b = bond["sites"]
        _require(
            a in adjacency and b in adjacency and a != b,
            "Bond references unknown or duplicate site",
        )
        adjacency[a].add(b)
        adjacency[b].add(a)
    _require(
        payload["connectivity"] == {str(i): sorted(adjacency[i]) for i in ids},
        "Graph connectivity disagrees with bonds",
    )
    components = []
    remaining = set(ids)
    while remaining:
        pending = [min(remaining)]
        component = set()
        while pending:
            i = pending.pop()
            if i in component:
                continue
            component.add(i)
            pending.extend(adjacency[i] - component)
        remaining -= component
        components.append(sorted(component))
    _require(
        payload["components"] == components, "Graph components disagree with bonds"
    )
    _require(
        payload["component_membership"]
        == {str(i): n for n, component in enumerate(components) for i in component},
        "Graph component membership mismatch",
    )
    return payload, {row["id"]: row for row in rows}


def _maximum_valence(element):
    # This is a structural ceiling; it is not a reaction or typing rule.
    return {
        "H": 1,
        "C": 4,
        "N": 4,
        "O": 2,
        "F": 1,
        "P": 5,
        "S": 6,
        "Cl": 1,
        "Br": 1,
        "I": 1,
    }.get(element)


def identify_reactive_sites(graph, rule):
    """Select only caller-marked graph sites, sorted by stable atom ID."""
    _require(isinstance(rule, ReactiveSiteRule), "ReactiveSiteRule required")
    payload, atoms = _graph_sites(graph)
    valence = {i: 0.0 for i in atoms}
    for bond in payload["bonds"]:
        for atom_id in bond["sites"]:
            valence[atom_id] += bond["order"]
    components = {
        i: _digest(component) for component in payload["components"] for i in component
    }
    result = []
    for atom_id, atom in sorted(atoms.items()):
        metadata = atom["metadata"]
        if rule.metadata_key not in metadata:
            continue
        marker = metadata[rule.metadata_key]
        _require(
            type(marker) is type(rule.metadata_value),
            f"Malformed reactive marker at site {atom_id}",
        )
        if marker != rule.metadata_value:
            continue
        _require(
            atom["element"] in rule.allowed_elements,
            f"Reactive site {atom_id} has disallowed element {atom['element']}",
        )
        _require(
            not atom["aromatic"], f"Aromatic reactive site {atom_id} is unsupported"
        )
        ceiling = _maximum_valence(atom["element"])
        _require(ceiling is not None, f"Unsupported reactive element at site {atom_id}")
        capacity = min(
            rule.maximum_crosslinks_per_site,
            floor((ceiling - valence[atom_id]) / rule.bond_order),
        )
        _require(capacity > 0, f"Reactive site {atom_id} has no valence capacity")
        result.append(
            ReactiveSite(
                atom_id,
                rule.name,
                atom["element"],
                components[atom_id],
                valence[atom_id],
                capacity,
                rule.identity,
                rule.provenance,
                rule.evidence,
                rule,
            )
        )
    return tuple(result)


def _candidates(graph_payload, sites, *, intercomponent):
    existing = {tuple(b["sites"]) for b in graph_payload["bonds"]}
    return tuple(
        (a.site_id, b.site_id)
        for a, b in combinations(sites, 2)
        if (not intercomponent or a.component_identity != b.component_identity)
        and (a.site_id, b.site_id) not in existing
    )


def _pair_identity(graph_identity, rule_identity, pair, order):
    return _digest(
        {
            "graph_identity": graph_identity,
            "rule_identity": rule_identity,
            "sites": list(pair),
            "bond_order": order,
        }
    )


def _selection(candidates, sites, target, seed, maximum):
    if target == 0:
        return []
    shuffled = list(candidates)
    random.Random(seed).shuffle(shuffled)
    limits = {site.site_id: min(site.available_capacity, maximum) for site in sites}
    usage = Counter()
    chosen = []
    for pair in shuffled:
        if all(usage[i] < limits[i] for i in pair):
            chosen.append(pair)
            usage.update(pair)
            if len(chosen) == target:
                break
    return chosen


@dataclass(frozen=True)
class CrosslinkNetworkPlan:
    json_text: str

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]

    def validate_integrity(self, graph=None):
        p = unpack(self.json_text)
        required = {
            "schema",
            "input_graph_identity",
            "reactive_sites",
            "selected_reactive_site_identities",
            "candidate_pairs",
            "candidate_pair_identities",
            "selected_crosslink_bonds",
            "target_crosslink_count",
            "target_conversion",
            "pairing_strategy",
            "require_intercomponent",
            "maximum_crosslinks_per_site",
            "random_seed",
            "bond_order",
            "rule",
            "rule_identity",
            "provenance",
            "evidence",
            "identity",
        }
        _require(
            type(p) is dict and set(p) == required and p["schema"] == PLAN_SCHEMA,
            "Malformed crosslink plan schema or fields",
        )
        _require(
            p["identity"] == _digest({k: v for k, v in p.items() if k != "identity"}),
            "Crosslink plan identity mismatch",
        )
        _text(p["provenance"], "Plan provenance")
        _evidence(p["evidence"])
        rule = _rule_from_payload(p["rule"])
        _require(
            rule.identity == p["rule_identity"], "Crosslink rule identity mismatch"
        )
        _integer(p["target_crosslink_count"], "Target crosslink count")
        _integer(p["maximum_crosslinks_per_site"], "Maximum site usage", minimum=1)
        _integer(p["random_seed"], "Seed")
        _require(
            type(p["require_intercomponent"]) is bool,
            "Intercomponent policy must be boolean",
        )
        _require(p["pairing_strategy"] == STRATEGY, "Unsupported pairing strategy")
        _require(p["bond_order"] == rule.bond_order, "Rule bond order mismatch")
        if p["target_conversion"] is not None:
            _number(p["target_conversion"], "Target conversion")
            _require(
                0 <= p["target_conversion"] <= 1,
                "Target conversion must be between zero and one",
            )
        sites = tuple(_site_from_payload(row) for row in p["reactive_sites"])
        _require(
            sites
            and tuple(s.site_id for s in sites)
            == tuple(sorted({s.site_id for s in sites})),
            "Reactive sites must have unique sorted IDs",
        )
        _require(
            all(s.rule_identity == rule.identity for s in sites),
            "Reactive site rule mismatch",
        )
        _require(
            p["selected_reactive_site_identities"] == [s.identity for s in sites],
            "Reactive site identity mismatch",
        )
        candidates = [_pair(*pair) for pair in p["candidate_pairs"]]
        _require(
            candidates == sorted(set(candidates)),
            "Duplicate or unordered candidate pair",
        )
        _require(
            p["candidate_pair_identities"]
            == [
                _pair_identity(
                    p["input_graph_identity"], rule.identity, pair, p["bond_order"]
                )
                for pair in candidates
            ],
            "Candidate pair identity mismatch",
        )
        selected = [_pair(*pair) for pair in p["selected_crosslink_bonds"]]
        _require(
            len(selected) == len(set(selected)) == p["target_crosslink_count"]
            and set(selected) <= set(candidates),
            "Duplicate, unknown or missing selected crosslink bond",
        )
        _require(
            selected
            == _selection(
                candidates,
                sites,
                p["target_crosslink_count"],
                p["random_seed"],
                p["maximum_crosslinks_per_site"],
            ),
            "Selected bonds disagree with seeded plan",
        )
        if graph is not None:
            payload, _ = _graph_sites(graph)
            _require(
                p["input_graph_identity"] == graph.identity,
                "Stale input graph identity",
            )
            _require(
                sites == identify_reactive_sites(graph, rule),
                "Reactive site records disagree with graph",
            )
            _require(
                candidates
                == list(
                    _candidates(
                        payload, sites, intercomponent=p["require_intercomponent"]
                    )
                ),
                "Candidate pairs disagree with graph",
            )


def plan_crosslinks(
    graph,
    reactive_sites,
    *,
    target_crosslinks,
    seed=2026,
    strategy=STRATEGY,
    require_intercomponent=True,
    maximum_crosslinks_per_site=1,
    bond_order=1.0,
    provenance,
    evidence,
    rule=None,
    target_conversion=None,
):
    """Shuffle lexicographically sorted eligible pairs using a local seeded RNG."""
    payload, _ = _graph_sites(graph)
    _integer(target_crosslinks, "Target crosslink count")
    _integer(seed, "Seed")
    _integer(maximum_crosslinks_per_site, "Maximum site usage", minimum=1)
    _require(
        type(require_intercomponent) is bool, "Intercomponent policy must be boolean"
    )
    _require(strategy == STRATEGY, "Unsupported pairing strategy")
    _number(bond_order, "Bond order", positive=True)
    _text(provenance, "Plan provenance")
    evidence = _evidence(evidence)
    sites = tuple(reactive_sites)
    _require(
        sites and all(isinstance(s, ReactiveSite) for s in sites),
        "Nonempty ReactiveSite sequence required",
    )
    if rule is None:
        rule = sites[0].rule
    _require(
        isinstance(rule, ReactiveSiteRule),
        "ReactiveSiteRule must be supplied or originate from identify_reactive_sites",
    )
    _require(bond_order == rule.bond_order, "Plan bond order differs from rule")
    _require(
        sites == identify_reactive_sites(graph, rule),
        "Reactive sites must exactly match rule and graph in stable ID order",
    )
    candidates = _candidates(payload, sites, intercomponent=require_intercomponent)
    selected = _selection(
        candidates, sites, target_crosslinks, seed, maximum_crosslinks_per_site
    )
    _require(
        len(selected) == target_crosslinks,
        f"Target {target_crosslinks} exceeds feasible candidate/capacity count {len(selected)}",
    )
    if target_conversion is not None:
        _number(target_conversion, "Target conversion")
        _require(
            0 <= target_conversion <= 1,
            "Target conversion must be between zero and one",
        )
    p = {
        "schema": PLAN_SCHEMA,
        "input_graph_identity": graph.identity,
        "reactive_sites": [_site_payload(s) for s in sites],
        "selected_reactive_site_identities": [s.identity for s in sites],
        "candidate_pairs": [list(pair) for pair in candidates],
        "candidate_pair_identities": [
            _pair_identity(graph.identity, rule.identity, pair, bond_order)
            for pair in candidates
        ],
        "selected_crosslink_bonds": [list(pair) for pair in selected],
        "target_crosslink_count": target_crosslinks,
        "target_conversion": target_conversion,
        "pairing_strategy": strategy,
        "require_intercomponent": require_intercomponent,
        "maximum_crosslinks_per_site": maximum_crosslinks_per_site,
        "random_seed": seed,
        "bond_order": bond_order,
        "rule": _rule_payload(rule),
        "rule_identity": rule.identity,
        "provenance": provenance,
        "evidence": list(evidence),
    }
    p["identity"] = _digest(p)
    result = CrosslinkNetworkPlan(pack(p))
    result.validate_integrity(graph)
    return result


def apply_crosslink_plan(system, graph, plan):
    """Apply selected bonds to a copy and bind the plan into graph history."""
    _require(isinstance(plan, CrosslinkNetworkPlan), "CrosslinkNetworkPlan required")
    graph.validate_integrity(system)
    for atom_id in system.topology.sites:
        _require(
            all(isfinite(float(x)) for x in system.coordinates.get(atom_id)),
            f"Nonfinite coordinates at site {atom_id}",
        )
    plan.validate_integrity(graph)
    p = plan.payload
    prior = deepcopy(graph.payload["transformations"])
    _require(
        prior == system.metadata.get("final_graph_transformations", []),
        "Stale graph transformation history",
    )
    result = system.copy()
    for a, b in p["selected_crosslink_bonds"]:
        _require(
            (a, b) not in result.topology.bonds, "Crosslink duplicates existing bond"
        )
        result.topology.add_bond(a, b, order=p["bond_order"])
    result.topology.rebuild_derived_interactions()
    result.metadata["crosslinks"] = deepcopy(result.metadata.get("crosslinks", [])) + [
        {"sites": pair, "seed": p["random_seed"], "plan_identity": plan.identity}
        for pair in p["selected_crosslink_bonds"]
    ]
    parameters = {
        "plan_identity": plan.identity,
        "bonds": p["selected_crosslink_bonds"],
        "bond_order": p["bond_order"],
    }
    after = final_graph(result)
    tr = transformation(
        graph,
        after,
        "planned_crosslink_network",
        parameters=parameters,
        seed=p["random_seed"],
        provenance=p["provenance"],
        evidence=p["evidence"],
    )
    output = final_graph(result, transformations=prior + [tr.payload])
    tr = transformation(
        graph,
        output,
        "planned_crosslink_network",
        parameters=parameters,
        seed=p["random_seed"],
        provenance=p["provenance"],
        evidence=p["evidence"],
    )
    output = final_graph(result, transformations=prior + [tr.payload])
    result.metadata["final_graph_transformations"] = deepcopy(
        output.payload["transformations"]
    )
    output.validate_integrity(result)
    tr.validate_integrity(graph, output)
    return result, output, tr


@dataclass(frozen=True)
class CrosslinkNetworkResult:
    system: object
    graph: FinalChemicalGraph
    reactive_sites: tuple[ReactiveSite, ...]
    plan: CrosslinkNetworkPlan
    transformation: object


def generate_crosslink_network(
    system,
    graph,
    rule,
    *,
    target_crosslinks,
    seed=2026,
    strategy=STRATEGY,
    require_intercomponent=True,
    maximum_crosslinks_per_site=1,
    bond_order=1.0,
    provenance,
    evidence,
):
    graph.validate_integrity(system)
    sites = identify_reactive_sites(graph, rule)
    plan = plan_crosslinks(
        graph,
        sites,
        target_crosslinks=target_crosslinks,
        seed=seed,
        strategy=strategy,
        require_intercomponent=require_intercomponent,
        maximum_crosslinks_per_site=maximum_crosslinks_per_site,
        bond_order=bond_order,
        provenance=provenance,
        evidence=evidence,
        rule=rule,
    )
    copied, output, tr = apply_crosslink_plan(system, graph, plan)
    return CrosslinkNetworkResult(copied, output, sites, plan, tr)


def save_crosslink_plan(plan, path):
    _require(isinstance(plan, CrosslinkNetworkPlan), "CrosslinkNetworkPlan required")
    plan.validate_integrity()
    target = Path(path)
    _require(not target.exists() and not target.is_symlink(), "Plan destination exists")
    _require(target.parent.is_dir(), "Plan parent does not exist")
    target.write_text(plan.json_text, encoding="utf-8")
    return target


def load_crosslink_plan(path, graph, *, expected_identity=None):
    plan = CrosslinkNetworkPlan(Path(path).read_text(encoding="utf-8"))
    plan.validate_integrity(graph)
    if expected_identity is not None:
        _require(
            plan.identity == expected_identity,
            "Crosslink plan identity differs from trusted reference",
        )
    return plan
