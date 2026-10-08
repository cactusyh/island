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

PLAN_SCHEMA = "island_crosslink_network_plan_v2"
RULE_SCHEMA = "island_reactive_site_rule_v1"
SITE_SCHEMA = "island_reactive_site_v2"
STRUCTURAL_POLICY = "island_crosslink_cho_n_structural_policy_v1"
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
    "lifetime_capacity",
    "consumed_crosslinks",
    "remaining_capacity",
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
        "lifetime_capacity": site.lifetime_capacity,
        "consumed_crosslinks": site.consumed_crosslinks,
        "remaining_capacity": site.remaining_capacity,
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
    lifetime_capacity: int = 1
    consumed_crosslinks: int = 0
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
        _integer(self.lifetime_capacity, "Lifetime capacity", minimum=1)
        _integer(self.consumed_crosslinks, "Consumed crosslinks")
        _require(
            0
            <= self.available_capacity
            <= self.lifetime_capacity - self.consumed_crosslinks,
            "Reactive site lifetime accounting mismatch",
        )
        object.__setattr__(self, "evidence", _evidence(self.evidence))

    @property
    def remaining_capacity(self):
        """Unused lifetime allowance, before the chemical valence restriction."""
        return self.lifetime_capacity - self.consumed_crosslinks

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
    _integer(payload["remaining_capacity"], "Remaining capacity")
    _require(
        payload["remaining_capacity"]
        == payload["lifetime_capacity"] - payload["consumed_crosslinks"],
        "Remaining lifetime capacity mismatch",
    )
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
        payload["lifetime_capacity"],
        payload["consumed_crosslinks"],
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


def _structural_states(payload, selected=(), order=1.0):
    """Bounded graph checks, including the proposed output, without state repair.

    Neutral C with valence 2/3 is accepted only as a caller-declared incomplete
    topology fixture; this does not assert a closed-shell electronic structure.
    Explicit radicals are unsupported because this operation cannot update them.
    """
    incident = {atom["id"]: [] for atom in payload["sites"]}
    for bond in payload["bonds"]:
        for atom_id in bond["sites"]:
            incident[atom_id].append(bond["order"])
    for pair in selected:
        for atom_id in pair:
            _require(atom_id in incident, "Unknown selected site")
            incident[atom_id].append(order)
    for atom in payload["sites"]:
        i, element, charge = atom["id"], atom["element"], atom["formal_charge"]
        meta = atom["metadata"]
        _require(type(charge) is int, f"Invalid formal charge at site {i}")
        _require(type(atom["aromatic"]) is bool, f"Invalid aromaticity at site {i}")
        _require(
            not atom["aromatic"] and not meta.get("aromatic", False),
            f"Aromatic site {i} is unsupported by {STRUCTURAL_POLICY}",
        )
        radical = meta.get("radical_electrons", 0)
        _require(
            type(radical) is int and radical == 0,
            f"Unsupported radical state at site {i}",
        )
        _require(
            not any(
                "radical" in key.lower() and key != "radical_electrons" for key in meta
            ),
            f"Ambiguous radical metadata at site {i}",
        )
        orders = incident[i]
        _require(
            all(type(v) in (int, float) and v in (1, 2, 3) for v in orders),
            f"Unsupported bond orders at site {i}",
        )
        valence = sum(orders)
        states = {
            ("H", 0): (1,),
            ("C", 0): (2, 3, 4),
            ("O", 0): (1, 2),
            ("N", 0): (3,),
            ("N", 1): (4,),
        }
        _require(
            (element, charge) in states,
            f"Unsupported element/charge state at site {i}: {element} {charge}",
        )
        _require(
            valence in states[(element, charge)],
            f"Invalid charge-aware valence at site {i}: {element} {charge}, valence {valence}",
        )
        if element == "N" and charge == 1:
            _require(
                orders == [1] * 4,
                f"Unsupported positively charged nitrogen state at site {i}",
            )
    return incident


def _consumed_crosslinks(payload):
    """Reconcile every historical crosslink with provenance and a current bond."""
    bonds = {tuple(b["sites"]): b["order"] for b in payload["bonds"]}
    recorded = {}
    limits = {}
    for tr in payload["transformations"]:
        if tr["operation"] not in ("crosslink_bonds", "planned_crosslink_network"):
            continue
        parameters = tr["parameters"]
        planned = tr["operation"] == "planned_crosslink_network"
        _require(
            type(parameters) is dict and type(parameters.get("bonds")) is list,
            "Malformed historical crosslink parameters",
        )
        _require(
            not planned or "bond_order" in parameters, "Missing historical bond order"
        )
        _require(
            "structural_policy" not in parameters
            or parameters["structural_policy"] == STRUCTURAL_POLICY,
            "Unsupported historical structural policy",
        )
        declared_limits = parameters.get("site_lifetime_capacities", {})
        _require(
            type(declared_limits) is dict, "Malformed historical lifetime capacities"
        )
        if "structural_policy" in parameters:
            _require(
                planned
                and "site_lifetime_capacities" in parameters
                and {str(i) for pair in parameters["bonds"] for i in pair}
                <= set(declared_limits),
                "Missing historical lifetime capacity evidence",
            )
        order = parameters["bond_order"] if planned else 1
        _number(order, "Historical crosslink bond order", positive=True)
        plan_id = parameters.get("plan_identity") if planned else None
        _require(
            not planned or (type(plan_id) is str and len(plan_id) == 64),
            "Ambiguous historical plan identity",
        )
        for row in parameters.get("bonds", []):
            pair = _pair(*row)
            _require(pair not in recorded, "Duplicate historical crosslink accounting")
            _require(
                bonds.get(pair) == order,
                "Historical crosslink disagrees with actual bond",
            )
            recorded[pair] = (tr["seed"], plan_id)
        for key, value in declared_limits.items():
            _integer(value, "Historical lifetime capacity", minimum=1)
            _require(
                key in payload["molecule_membership"],
                "Unknown historical capacity site",
            )
            _require(
                key not in limits or limits[key] == value,
                "Contradictory lifetime capacities",
            )
            limits[key] = value
    provenance = payload["crosslink_provenance"]
    _require(type(provenance) is list, "Malformed crosslink provenance")
    seen = set()
    consumed = Counter()
    for row in provenance:
        _require(
            type(row) is dict
            and set(row) in ({"sites", "seed"}, {"sites", "seed", "plan_identity"}),
            "Ambiguous crosslink provenance",
        )
        pair = _pair(*row["sites"])
        _require(pair not in seen, "Duplicate crosslink provenance")
        _require(pair in recorded, "Crosslink provenance lacks transformation evidence")
        _require(
            recorded[pair] == (row["seed"], row.get("plan_identity")),
            "Crosslink provenance contradicts transformation",
        )
        seen.add(pair)
        consumed.update(pair)
    _require(seen == set(recorded), "Missing crosslink provenance")
    _require(
        all(consumed[int(i)] <= cap for i, cap in limits.items()),
        "Consumed crosslinks exceed historical lifetime capacity",
    )
    return consumed, limits


def identify_reactive_sites(graph, rule):
    """Select only caller-marked graph sites, sorted by stable atom ID."""
    _require(isinstance(rule, ReactiveSiteRule), "ReactiveSiteRule required")
    payload, atoms = _graph_sites(graph)
    incident = _structural_states(payload)
    valence = {i: sum(orders) for i, orders in incident.items()}
    consumed, previous_limits = _consumed_crosslinks(payload)
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
        _require(
            atom["element"] in {"C", "O"} and atom["formal_charge"] == 0,
            f"Unsupported reactive element/charge at site {atom_id}; only neutral C/O sites are enabled",
        )
        ceiling = {"C": 4, "O": 2}[atom["element"]]
        _require(
            str(atom_id) not in previous_limits
            or previous_limits[str(atom_id)] == rule.maximum_crosslinks_per_site,
            f"Lifetime capacity changed at site {atom_id}",
        )
        _require(
            consumed[atom_id] <= rule.maximum_crosslinks_per_site,
            f"Consumed crosslinks exceed lifetime capacity at site {atom_id}",
        )
        capacity = min(
            rule.maximum_crosslinks_per_site - consumed[atom_id],
            floor((ceiling - valence[atom_id]) / rule.bond_order),
        )
        # Exhausted sites remain observable, but have zero planning capacity.
        _require(capacity >= 0, f"Reactive site {atom_id} has invalid valence capacity")
        _require(
            capacity > 0 or consumed[atom_id] > 0,
            f"Reactive site {atom_id} has no valence capacity",
        )
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
                rule.maximum_crosslinks_per_site,
                consumed[atom_id],
                rule,
            )
        )
    return tuple(result)


def _candidates(graph_payload, sites, *, intercomponent):
    existing = {tuple(b["sites"]) for b in graph_payload["bonds"]}
    return tuple(
        (a.site_id, b.site_id)
        for a, b in combinations(sites, 2)
        if a.available_capacity > 0
        and b.available_capacity > 0
        and (not intercomponent or a.component_identity != b.component_identity)
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
    limits = {
        site.site_id: site.available_capacity
        if maximum is None
        else min(site.available_capacity, maximum)
        for site in sites
    }
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
            "input_history_identity",
            "structural_policy",
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
        _require(
            p["structural_policy"] == STRUCTURAL_POLICY, "Unsupported structural policy"
        )
        _number(p["bond_order"], "Bond order", positive=True)
        _text(p["provenance"], "Plan provenance")
        _evidence(p["evidence"])
        rule = _rule_from_payload(p["rule"])
        _require(
            rule.identity == p["rule_identity"], "Crosslink rule identity mismatch"
        )
        _integer(p["target_crosslink_count"], "Target crosslink count")
        if p["maximum_crosslinks_per_site"] is not None:
            _integer(
                p["maximum_crosslinks_per_site"], "Per-batch site limit", minimum=1
            )
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
            all(i in {s.site_id for s in sites} for pair in candidates for i in pair),
            "Candidate pair contains unknown reactive site",
        )
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
                p["input_history_identity"] == payload["history_identity"],
                "Stale input graph history",
            )
            _structural_states(payload, selected, p["bond_order"])
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
    if maximum_crosslinks_per_site is not None:
        _integer(maximum_crosslinks_per_site, "Per-batch site limit", minimum=1)
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
        "input_history_identity": payload["history_identity"],
        "structural_policy": STRUCTURAL_POLICY,
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
        "structural_policy": STRUCTURAL_POLICY,
        "site_lifetime_capacities": {
            str(s["site_id"]): s["lifetime_capacity"] for s in p["reactive_sites"]
        },
    }
    membership = {
        int(i): label for i, label in graph.payload["molecule_membership"].items()
    }
    after = final_graph(result, molecule_membership=membership)
    _structural_states(after.payload)
    tr = transformation(
        graph,
        after,
        "planned_crosslink_network",
        parameters=parameters,
        seed=p["random_seed"],
        provenance=p["provenance"],
        evidence=p["evidence"],
    )
    output = final_graph(
        result, transformations=prior + [tr.payload], molecule_membership=membership
    )
    tr = transformation(
        graph,
        output,
        "planned_crosslink_network",
        parameters=parameters,
        seed=p["random_seed"],
        provenance=p["provenance"],
        evidence=p["evidence"],
    )
    output = final_graph(
        result, transformations=prior + [tr.payload], molecule_membership=membership
    )
    result.metadata["final_graph_transformations"] = deepcopy(
        output.payload["transformations"]
    )
    output.validate_integrity(result)
    _structural_states(output.payload)
    _consumed_crosslinks(output.payload)
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
