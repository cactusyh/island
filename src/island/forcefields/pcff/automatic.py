"""Bounded ISLAND graph rules, independent of construction and coordinates.

Historical manual records and the charge resolver are reused without modification.
"""

import json
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from math import isfinite
from pathlib import Path

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .charges import PCFFChargeResult, assign_pcff_charges, assign_pcff_types, identity
from .source import FLAGS, PCFFSource, boundary, records, require, select

PROFILE = json.loads(Path(__file__).with_name("typing_profile.json").read_text())
AUTO_SCHEMA = "island_pcff_automatic_typing_v1"
BRIDGE_SCHEMA = "island_pcff_automatic_charges_v1"


@boundary
def chemical_graph(system):
    """Only authoritative chemical attributes; coordinate/repeat history is not input."""
    system.topology.validate_bond_graph()
    require(bool(system.topology.sites), "Empty graph")
    sites = []
    for sid, atom in sorted(system.topology.sites.items()):
        require(
            type(sid) is int and isinstance(atom, AtomSite) and atom.id == sid,
            "Exact integer AtomSite IDs required",
        )
        require(
            type(atom.element) is str
            and bool(atom.element)
            and type(atom.atomic_number) is int
            and 1 <= atom.atomic_number <= 118,
            "Malformed element/atomic number",
        )
        if atom.element in {"C", "H", "O"}:
            require(
                atom.atomic_number == {"C": 6, "H": 1, "O": 8}[atom.element],
                "Inconsistent element/atomic number",
            )
        require(
            type(atom.mass) in (int, float) and isfinite(atom.mass) and atom.mass > 0,
            "Invalid mass",
        )
        require(type(atom.formal_charge) is int, "Invalid formal charge")
        require(type(atom.metadata) is dict, "Malformed chemical metadata")
        meta = {
            "aromatic": atom.metadata.get("aromatic", False),
            "isotope": atom.metadata.get("isotope", 0),
            "radical_electrons": atom.metadata.get("radical_electrons", 0),
            "chiral_tag": atom.metadata.get("chiral_tag", "CHI_UNSPECIFIED"),
            "cip_label": atom.metadata.get("cip_label"),
        }
        require(
            type(meta["aromatic"]) is bool
            and all(
                type(meta[k]) is int and meta[k] >= 0
                for k in ("isotope", "radical_electrons")
            ),
            "Malformed chemical state metadata",
        )
        require(
            type(meta["chiral_tag"]) is str
            and (meta["cip_label"] is None or type(meta["cip_label"]) is str),
            "Malformed stereo metadata",
        )
        sites.append(
            {
                "id": sid,
                "element": atom.element,
                "atomic_number": atom.atomic_number,
                "mass": atom.mass,
                "formal_charge": atom.formal_charge,
                "metadata": meta,
            }
        )
    bonds = []
    for key, bond in sorted(system.topology.bonds.items()):
        require(
            key == bond.key
            and type(bond.order) in (int, float)
            and isfinite(bond.order)
            and bond.order > 0
            and type(bond.aromatic) is bool,
            "Malformed bond",
        )
        bonds.append(
            {"sites": list(key), "order": bond.order, "aromatic": bond.aromatic}
        )
    require(type(system.representation) is str, "Malformed representation")
    return {
        "sites": sites,
        "bonds": bonds,
        "representation": system.representation,
        "has_box": system.box is not None,
    }


def graph_system(graph):
    require(
        set(graph) == {"sites", "bonds", "representation", "has_box"},
        "Malformed graph record",
    )
    require(type(graph["has_box"]) is bool, "Malformed periodicity")
    topology = Topology()
    for site in graph["sites"]:
        require(
            set(site)
            == {"id", "element", "atomic_number", "mass", "formal_charge", "metadata"},
            "Malformed site record",
        )
        topology.add_site(AtomSite(name=str(site["id"]), **deepcopy(site)))
    for bond in graph["bonds"]:
        require(set(bond) == {"sites", "order", "aromatic"}, "Malformed bond record")
        a, b = bond["sites"]
        require((a, b) not in topology.bonds, "Duplicate bond")
        topology.add_bond(a, b, order=bond["order"], aromatic=bond["aromatic"])
    system = MolecularSystem(
        topology, Coordinates(), representation=graph["representation"]
    )
    # Only presence is relevant to an unsupported periodic diagnostic; never evaluate a box.
    if graph["has_box"]:
        system.box = object()
    require(pack(chemical_graph(system)) == pack(graph), "Noncanonical chemical graph")
    return system


def _components(neighbors):
    remaining, components = set(neighbors), []
    while remaining:
        queue, found = [min(remaining)], set()
        while queue:
            site = queue.pop()
            if site not in found:
                found.add(site)
                queue.extend(neighbors[site] - found)
        remaining -= found
        components.append(sorted(found))
    return components


def typing_data(graph, source, profile):
    from .expanded import PROFILE_NAME
    from .expanded import typing_data as expanded_typing

    if profile in (
        PROFILE_NAME,
        "island_pcff_source_graph_v2",
        "island_pcff_source_graph_v3",
        "island_pcff_source_graph_v4",
    ):
        return expanded_typing(graph, source, version=int(profile[-1]))
    source.require_assignment()
    require(profile == PROFILE["name"], "Unsupported PCFF automatic typing profile")
    require(
        source.identity["sha256"] == PROFILE["source_sha256"],
        "Typing profile/source mismatch",
    )
    graph_system(graph)
    sites = {s["id"]: s for s in graph["sites"]}
    neighbors = {sid: set() for sid in sites}
    for bond in graph["bonds"]:
        a, b = bond["sites"]
        neighbors[a].add(b)
        neighbors[b].add(a)
    environments = {}
    issues = []
    for sid, atom in sites.items():
        counts = Counter(sites[n]["element"] for n in neighbors[sid])
        environments[sid] = {
            "element": atom["element"],
            "degree": len(neighbors[sid]),
            "neighbors": sorted(neighbors[sid]),
            "neighbor_elements": dict(sorted(counts.items())),
            "hydrogens": counts["H"],
            "oxygen_neighbors": counts["O"],
        }
        reasons = []
        e, meta = atom["element"], atom["metadata"]
        if e not in ("C", "H", "O"):
            reasons.append("unsupported_element")
        else:
            if len(neighbors[sid]) != {"C": 4, "H": 1, "O": 2}[e]:
                reasons.append("incomplete_hydrogens_or_valence")
            if not any(
                abs(atom["mass"] - mass) <= PROFILE["mass_atol_dalton"]
                for mass in PROFILE["standard_masses"][e]
            ):
                reasons.append("nonstandard_mass_or_isotope")
        if atom["formal_charge"] != 0:
            reasons.append("charged_site")
        if meta["isotope"] or meta["radical_electrons"]:
            reasons.append("isotope_or_radical")
        if meta["aromatic"]:
            reasons.append("aromaticity")
        if e == "C" and counts["O"] > 1:
            reasons.append("multiple_oxygen_carbon")
        if e == "O" and counts["O"]:
            reasons.append("peroxide")
        if e == "O" and counts["H"] == 2:
            reasons.append("water")
        if e == "H" and not (counts["C"] == 1 or counts["O"] == 1):
            reasons.append("unsupported_hydrogen_parent")
        for reason in reasons:
            issues.append({"sites": [sid], "reason": reason})
    for bond in graph["bonds"]:
        if bond["order"] != 1 or bond["aromatic"]:
            issues.append(
                {"sites": bond["sites"], "reason": "unsaturated_or_aromatic_bond"}
            )
    if graph["representation"] != "atomistic" or graph["has_box"]:
        issues.append(
            {"sites": sorted(sites), "reason": "unsupported_representation_or_box"}
        )
    blocked = set()
    for component in _components(neighbors):
        edges = sum(len(neighbors[s]) for s in component) // 2
        if edges >= len(component):
            issues.append({"sites": component, "reason": "cyclic_component"})
        if any(set(i["sites"]).intersection(component) for i in issues):
            blocked.update(component)
    entries, assignments = {}, {}
    atom_rows = records(source.inventory, "atom_types")
    for sid, env in environments.items():
        entry = {
            "environment": env,
            "matched_rules": [],
            "overridden_rules": [],
            "selected_rule": None,
            "type": None,
            "status": "unsupported" if sid in blocked else "untyped",
        }
        if sid not in blocked:

            def matches(rule, sid=sid, env=env):
                for key, expected in rule["predicate"].items():
                    if key == "parent_element":
                        actual = sites[next(iter(neighbors[sid]))]["element"]
                    else:
                        actual = env[key]
                    if actual != expected:
                        return False
                return True

            matches_ = [r for r in PROFILE["rules"] if matches(r)]
            overridden = {i for r in matches_ for i in r["overrides"]}
            active = [r for r in matches_ if r["id"] not in overridden]
            entry["matched_rules"] = sorted(r["id"] for r in matches_)
            entry["overridden_rules"] = sorted(
                overridden.intersection(entry["matched_rules"])
            )
            if len(active) == 1:
                rule = active[0]
                selected = select(
                    [r for r in atom_rows if r["data"]["type"] == rule["type"]]
                )
                require(
                    selected is not None
                    and selected["record"]["data"]["element"] == env["element"]
                    and selected["record"]["data"]["connections"] == env["degree"],
                    "Rule/source atom-type inconsistency",
                )
                entry.update(
                    selected_rule=rule["id"],
                    type=rule["type"],
                    status="typed",
                    rule_version=rule["version"],
                    source_atom_record=selected,
                )
                assignments[sid] = rule["type"]
            else:
                entry["status"] = "ambiguous" if active else "untyped"
                issues.append(
                    {
                        "sites": [sid],
                        "reason": entry["status"],
                        "rules": entry["matched_rules"],
                    }
                )
        entries[sid] = entry
    return {
        "schema": AUTO_SCHEMA,
        "source": source.identity,
        "graph": graph,
        "graph_identity": identity(graph),
        "profile": deepcopy(PROFILE),
        "profile_identity": identity(PROFILE),
        "assignments": assignments,
        "entries": entries,
        "diagnostics": issues,
        "coverage": {
            "total": len(sites),
            "typed": len(assignments),
            "complete": not issues and len(assignments) == len(sites),
        },
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFAutomaticTypingResult:
    """Owned decisions and diagnostic coverage for one exact graph/profile/source."""

    json_text: str
    source: PCFFSource

    @boundary
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        from .expanded import TYPING_SCHEMA
        from .expanded import typing_data as expanded_typing

        if p["schema"] in (
            TYPING_SCHEMA,
            "island_pcff_source_typing_v2",
            "island_pcff_source_typing_v3",
            "island_pcff_source_typing_v4",
        ):
            expected = expanded_typing(
                p["graph"],
                self.source,
                p["explicit_types"],
                p["explicit_provenance"],
                version=int(p["schema"][-1]),
            )
        else:
            expected = typing_data(p["graph"], self.source, p["profile"]["name"])
        require(pack(p) == pack(expected), "Contradictory automatic typing record")
        if system is not None:
            require(
                pack(chemical_graph(system)) == pack(p["graph"]),
                "Automatic typing chemical graph mismatch",
            )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @property
    def complete(self):
        return self.payload["coverage"]["complete"]

    @property
    def assignments(self):
        p = self.payload
        require(
            p["coverage"]["complete"],
            "Incomplete automatic typing; inspect diagnostics",
        )
        return p["assignments"]


@boundary
def type_pcff_atoms(system, source, *, profile="island_pcff_acyclic_cho_v1"):
    result = PCFFAutomaticTypingResult(
        pack(typing_data(chemical_graph(system), source, profile)), source
    )
    result.validate_integrity(system)
    return result


def bridge_data(automatic, source, resolution_policy=None):
    from .expanded import TYPING_SCHEMA
    from .expanded import charge_data as expanded_charges

    if automatic["schema"] in (
        TYPING_SCHEMA,
        "island_pcff_source_typing_v2",
        "island_pcff_source_typing_v3",
        "island_pcff_source_typing_v4",
    ):
        return expanded_charges(automatic, source, resolution_policy=resolution_policy)
    require(resolution_policy is None, "Fallback policy requires expanded typing")
    require(
        automatic["coverage"]["complete"],
        "Incomplete automatic typing cannot assign charges",
    )
    canonical = graph_system(automatic["graph"])
    provenance = f"ISLAND automatic typing {automatic['profile']['name']}; record {identity(automatic)}"
    explicit = assign_pcff_types(
        canonical, source, automatic["assignments"], provenance=provenance
    )
    native = assign_pcff_charges(canonical, explicit)
    return {
        "schema": BRIDGE_SCHEMA,
        "automatic_typing": automatic,
        "automatic_typing_identity": identity(automatic),
        "bridge": "canonical_chemical_graph_explicit_types_v1",
        "native_charge_record": native.payload,
        "native_charge_identity": native.identity,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFAutomaticChargeResult:
    """Automatic evidence plus the unchanged native resolver on a chemical snapshot."""

    json_text: str
    source: PCFFSource

    @boundary
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        auto = PCFFAutomaticTypingResult(pack(p["automatic_typing"]), self.source)
        auto.validate_integrity(system)
        expected = bridge_data(auto.payload, self.source, p.get("resolution_policy"))
        require(pack(p) == pack(expected), "Contradictory automatic charge bridge")

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @property
    def complete(self):
        return self.payload["native_charge_record"]["complete"]

    @property
    def charges(self):
        p = self.payload
        from .expanded import CHARGE_SCHEMA

        if p["schema"] in (CHARGE_SCHEMA, "island_pcff_source_charges_v2"):
            require(
                p["native_charge_record"]["complete"],
                "Incomplete native charges; inspect diagnostics",
            )
            return p["native_charge_record"]["partial_charges"]
        return PCFFChargeResult(pack(p["native_charge_record"]), self.source).charges


@boundary
def assign_automatic_pcff_charges(system, typing, *, resolution_policy=None):
    require(
        type(typing) is PCFFAutomaticTypingResult, "Expected automatic typing record"
    )
    typing.validate_integrity(system)
    result = PCFFAutomaticChargeResult(
        pack(bridge_data(typing.payload, typing.source, resolution_policy)),
        typing.source,
    )
    result.validate_integrity(system)
    return result


@boundary
def save_pcff_automatic_record(result, path):
    require(
        type(result) in (PCFFAutomaticTypingResult, PCFFAutomaticChargeResult),
        "Expected automatic PCFF record",
    )
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_automatic_record(path, source, *, system=None):
    raw = Path(path).read_text()
    p = unpack(raw)
    from .expanded import CHARGE_SCHEMA, TYPING_SCHEMA

    cls = {
        "island_pcff_source_typing_v2": PCFFAutomaticTypingResult,
        "island_pcff_source_typing_v3": PCFFAutomaticTypingResult,
        "island_pcff_source_typing_v4": PCFFAutomaticTypingResult,
        "island_pcff_source_charges_v2": PCFFAutomaticChargeResult,
        TYPING_SCHEMA: PCFFAutomaticTypingResult,
        CHARGE_SCHEMA: PCFFAutomaticChargeResult,
        AUTO_SCHEMA: PCFFAutomaticTypingResult,
        BRIDGE_SCHEMA: PCFFAutomaticChargeResult,
    }.get(p["schema"])
    require(cls is not None, "Unsupported automatic PCFF schema")
    result = cls(raw, source)
    result.validate_integrity(system)
    return result
