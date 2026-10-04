"""Explicit graph typing and bounded source-native bond-increment assignment."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from math import fsum, isfinite
from pathlib import Path

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.workflows.storage import publish

from .source import FLAGS, PIN, PCFFSource, boundary, digest, records, require, select

TYPING_SCHEMA = "island_pcff_explicit_typing_v1"
CHARGE_SCHEMA = "island_pcff_bond_increment_charges_v1"


def identity(value):
    return digest(pack(value).encode())


@boundary
def graph_data(system):
    require(
        system.representation == "atomistic" and system.box is None,
        "PCFF charge domain requires finite nonperiodic atomistic graphs",
    )
    system.topology.validate_bond_graph()
    topology = system.topology
    require(bool(topology.sites), "Empty graph")
    neighbors = {sid: [] for sid in topology.sites}
    bonds = []
    for key, bond in sorted(topology.bonds.items()):
        require(
            key == bond.key
            and type(bond.order) in (int, float)
            and bond.order == 1
            and bond.aromatic is False,
            "Initial PCFF charge domain supports saturated single bonds only",
        )
        a, b = sorted(key)
        neighbors[a].append(b)
        neighbors[b].append(a)
        bonds.append([a, b])
    sites = []
    elements = {"C": (6, 4), "H": (1, 1), "O": (8, 2)}
    for sid, atom in sorted(topology.sites.items()):
        require(
            type(sid) is int and isinstance(atom, AtomSite) and atom.id == sid,
            "Exact integer AtomSite IDs required",
        )
        require(
            atom.element in elements, "Initial PCFF charge domain supports C/H/O only"
        )
        number, valence = elements[atom.element]
        require(
            type(atom.atomic_number) is int and atom.atomic_number == number,
            "Element/atomic-number mismatch",
        )
        require(
            type(atom.formal_charge) is int and atom.formal_charge == 0,
            "Only individually neutral sites supported",
        )
        require(
            type(atom.mass) in (int, float) and isfinite(atom.mass) and atom.mass > 0,
            "Invalid mass",
        )
        require(
            not atom.metadata.get("isotope")
            and not atom.metadata.get("radical_electrons")
            and not atom.metadata.get("aromatic"),
            "Isotope/radical/aromatic scope unsupported",
        )
        require(
            len(neighbors[sid]) == valence,
            f"Site {sid}: incomplete explicit hydrogens or unsupported valence",
        )
        if atom.element == "H":
            require(
                topology.sites[neighbors[sid][0]].element in ("C", "O"),
                "Only carbon/oxygen-bound explicit hydrogen supported",
            )
        sites.append(
            {
                "id": sid,
                "element": atom.element,
                "atomic_number": number,
                "formal_charge": 0,
                "mass": atom.mass,
                "metadata": deepcopy(atom.metadata),
            }
        )
    # Metadata is bound, but coordinates, names, insertion order and interaction caches are not.
    return {
        "sites": sites,
        "bonds": bonds,
        "topology_metadata": deepcopy(topology.metadata)
        if hasattr(topology, "metadata")
        else {},
        "system_metadata": deepcopy(system.metadata),
    }


def system_from_graph(graph):
    require(
        set(graph) == {"sites", "bonds", "topology_metadata", "system_metadata"},
        "Malformed graph",
    )
    topology = Topology()
    for site in graph["sites"]:
        require(
            set(site)
            == {"id", "element", "atomic_number", "formal_charge", "mass", "metadata"},
            "Malformed site",
        )
        topology.add_site(AtomSite(name=str(site["id"]), **deepcopy(site)))
    for a, b in graph["bonds"]:
        require((a, b) not in topology.bonds, "Duplicate bond")
        topology.add_bond(a, b, order=1)
    if graph["topology_metadata"]:
        topology.metadata = deepcopy(graph["topology_metadata"])
    system = MolecularSystem(
        topology=topology,
        coordinates=Coordinates(),
        metadata=deepcopy(graph["system_metadata"]),
    )
    require(
        pack(graph_data(system)) == pack(graph), "Noncanonical/invalid graph inventory"
    )
    return system


def type_data(graph, supplied, provenance, source):
    source.require_assignment()
    system_from_graph(graph)
    require(
        type(supplied) is dict
        and all(type(k) is int and type(v) is str for k, v in supplied.items()),
        "Typing requires integer stable IDs and exact source labels",
    )
    require(
        set(supplied) == {s["id"] for s in graph["sites"]},
        "Typing site coverage mismatch",
    )
    require(
        type(provenance) is str and bool(provenance.strip()),
        "Explicit typing provenance required",
    )
    inventory = source.inventory
    rows = records(inventory, "atom_types")
    selected = {}
    adjacency = {s["id"]: [] for s in graph["sites"]}
    atoms = {s["id"]: s for s in graph["sites"]}
    for a, b in graph["bonds"]:
        adjacency[a].append(b)
        adjacency[b].append(a)
    # Narrow documented type domain, selected by chemical environments, not molecule name.
    supported = {"c", "c1", "c2", "c3", "h", "hc", "h*", "ho", "o", "oc", "oh"}
    for sid, label in sorted(supplied.items()):
        chosen = select([r for r in rows if r["data"]["type"] == label])
        require(chosen is not None, f"Unknown PCFF type {label!r} at site {sid}")
        require(
            label in supported,
            f"Type {label!r} parsed but outside audited charge domain",
        )
        data = chosen["record"]["data"]
        require(
            data["element"] == atoms[sid]["element"]
            and data["connections"] == len(adjacency[sid]),
            f"Type element/connectivity incompatible at {sid}",
        )
        hs = sum(atoms[n]["element"] == "H" for n in adjacency[sid])
        if label in ("c1", "c2", "c3"):
            require(hs == int(label[1]), f"Type hydrogen environment mismatch at {sid}")
        if label in ("h", "hc", "h*", "ho"):
            parent = atoms[adjacency[sid][0]]["element"]
            require(
                parent == ("C" if label in ("h", "hc") else "O"),
                f"Hydrogen parent environment mismatch at {sid}",
            )
        if label in ("oc", "oh"):
            require(
                hs == (1 if label == "oh" else 0),
                f"Oxygen environment mismatch at {sid}",
            )
        selected[sid] = chosen
    return {
        "schema": TYPING_SCHEMA,
        "source": source.identity,
        "graph": graph,
        "graph_identity": identity(graph),
        "types": dict(sorted(supplied.items())),
        "provenance": provenance,
        "selected_atom_records": selected,
        **FLAGS,
    }


def bond_selection(labels, rows):
    """Handle both orientations; inconsistent reversed rows are ambiguous, not precedence."""
    matches = []
    for oriented in (labels, labels[::-1]) if labels[0] != labels[1] else (labels,):
        chosen = select([r for r in rows if r["data"]["types"] == oriented])
        if chosen:
            values = chosen["record"]["data"]["increments"]
            reverse = oriented != labels
            matches.append(
                {
                    "selection": chosen,
                    "reversed": reverse,
                    "increments": values[::-1] if reverse else values,
                }
            )
    if not matches:
        return None
    require(
        all(m["increments"] == matches[0]["increments"] for m in matches),
        "Conflicting forward/reverse bond increment rows",
    )
    require(
        labels[0] != labels[1]
        or matches[0]["increments"][0] == matches[0]["increments"][1],
        "Identical endpoint types have ambiguous unequal increments",
    )
    return {"matches": matches, "increments": matches[0]["increments"]}


def charge_data(typing, source):
    inventory = source.inventory
    incs = records(inventory, "bond_increments")
    equivs = records(inventory, "equivalence")
    autos = records(inventory, "auto_equivalence")
    buckets = {sid: [] for sid in typing["types"]}
    contributions, diagnostics = [], []
    for a, b in typing["graph"]["bonds"]:
        labels = [typing["types"][a], typing["types"][b]]
        try:
            match = bond_selection(labels, incs)
            path, resolved, equivalences = "direct", labels, []
            if match is None:
                equivalences = [
                    select([r for r in equivs if r["data"]["type"] == t])
                    for t in labels
                ]
                if all(equivalences):
                    resolved = [
                        s["record"]["data"]["families"]["bond"] for s in equivalences
                    ]
                    match = bond_selection(resolved, incs)
                    path = "equivalence.bond"
            if match is None:
                auto = [
                    select([r for r in autos if r["data"]["type"] == t]) for t in labels
                ]
                diagnostics.append(
                    {
                        "bond": [a, b],
                        "types": labels,
                        "reason": "missing_increment",
                        "message": "No direct/bond-equivalence increment; auto fallback unsupported",
                        "auto_equivalence_evidence": auto,
                    }
                )
                continue
            values = match["increments"]
            buckets[a].append(values[0])
            buckets[b].append(values[1])
            contributions.append(
                {
                    "sites": [a, b],
                    "supplied_types": labels,
                    "resolved_types": resolved,
                    "path": path,
                    "equivalence_records": equivalences,
                    "source_matches": match["matches"],
                    "increments": values,
                }
            )
        except PCFFError as error:
            diagnostics.append(
                {
                    "bond": [a, b],
                    "types": labels,
                    "reason": "ambiguous_increment",
                    "message": str(error),
                }
            )
    charges = {sid: fsum(values) for sid, values in buckets.items()}
    neighbors = {sid: set() for sid in buckets}
    for a, b in typing["graph"]["bonds"]:
        neighbors[a].add(b)
        neighbors[b].add(a)
    remaining = set(buckets)
    components = []
    while remaining:
        todo, found = [min(remaining)], set()
        while todo:
            sid = todo.pop()
            if sid not in found:
                found.add(sid)
                todo.extend(neighbors[sid] - found)
        remaining -= found
        total = fsum(charges[s] for s in sorted(found))
        components.append({"sites": sorted(found), "formal_charge": 0, "charge": total})
        if abs(total) > PIN["tolerance_e"]:
            diagnostics.append(
                {
                    "sites": sorted(found),
                    "reason": "component_charge_mismatch",
                    "charge": total,
                }
            )
    total = fsum(charges.values())
    if abs(total) > PIN["tolerance_e"]:
        diagnostics.append({"reason": "molecular_charge_mismatch", "charge": total})
    return {
        "schema": CHARGE_SCHEMA,
        "typing": typing,
        "typing_identity": identity(typing),
        "source": source.identity,
        "policy": PIN["charge_policy"],
        "tolerance_e": PIN["tolerance_e"],
        "base_charge": 0.0,
        "contributions": contributions,
        "partial_charges": charges,
        "components": components,
        "total_charge": total,
        "formal_charge": 0,
        "unit": "elementary_charge",
        "complete": not diagnostics,
        "diagnostics": diagnostics,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFTypingResult:
    """Owned explicit labels, not evidence of automatic or scientifically correct typing."""

    json_text: str
    source: PCFFSource

    @boundary
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        expected = type_data(p["graph"], p["types"], p["provenance"], self.source)
        require(pack(p) == pack(expected), "Contradictory PCFF typing record")
        if system is not None:
            require(
                pack(graph_data(system)) == pack(p["graph"]),
                "PCFF bound graph mismatch",
            )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@dataclass(frozen=True)
class PCFFChargeResult:
    """Complete or diagnostic charge evidence. Never a parameterized system."""

    json_text: str
    source: PCFFSource

    @boundary
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        typing = PCFFTypingResult(pack(p["typing"]), self.source)
        typing.validate_integrity(system)
        expected = charge_data(typing.payload, self.source)
        require(pack(p) == pack(expected), "Contradictory PCFF charge record")

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @property
    def complete(self):
        return self.payload["complete"]

    @property
    def charges(self):
        p = self.payload
        require(
            p["complete"],
            "Incomplete PCFF assignment; inspect diagnostics, no usable charges",
        )
        return p["partial_charges"]


@boundary
def assign_pcff_types(system, source, types, *, provenance):
    require(
        isinstance(types, Mapping),
        "Explicit types must be a mapping, not a pair sequence",
    )
    result = PCFFTypingResult(
        pack(type_data(graph_data(system), dict(types), provenance, source)), source
    )
    result.validate_integrity(system)
    return result


@boundary
def assign_pcff_charges(system, typing):
    typing.validate_integrity(system)
    result = PCFFChargeResult(
        pack(charge_data(typing.payload, typing.source)), typing.source
    )
    result.validate_integrity(system)
    return result


@boundary
def save_pcff_record(result, path):
    require(
        type(result) in (PCFFTypingResult, PCFFChargeResult), "Expected PCFF record"
    )
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_record(path, source, *, system=None):
    raw = Path(path).read_text()
    p = unpack(raw)
    cls = {TYPING_SCHEMA: PCFFTypingResult, CHARGE_SCHEMA: PCFFChargeResult}.get(
        p["schema"]
    )
    require(cls is not None, "Unsupported PCFF record schema")
    result = cls(raw, source)
    result.validate_integrity(system)
    return result
