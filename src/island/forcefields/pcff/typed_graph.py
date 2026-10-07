"""External PCFF inputs bound to an exact graph, independently of perception.

Validation certifies the declared structural constraints, never the caller's
chemical authority. No coordinate, molecule name or repeat-unit name is read.
"""

from dataclasses import dataclass
from math import fsum, isfinite
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .automatic import _components, chemical_graph, graph_system
from .charges import identity
from .expanded import MASSES, NUMBERS, environments
from .source import FLAGS, PIN, PCFFSource, boundary, records, require, select
from .validation_cache import record_identity, validated_record

TYPING_SCHEMA = "island_pcff_external_typed_graph_v1"
CHARGE_SCHEMA = "island_pcff_typed_graph_charges_v1"
VALIDATION_POLICY = "island_pcff_external_structural_constraints_v1"
PROFILE_NAME = "island_pcff_external_source_graph_v1"
# This is a constraint registry, not a label-selection rule or an alias table.
SUPPORTED = {
    "c",
    "c1",
    "c2",
    "c3",
    "h",
    "hc",
    "o",
    "oc",
    "oh",
    "ho",
    "cp",
    "ct",
    "nt",
    "na",
    "hn",
    "hn2",
    "n+",
    "n4",
}


def _evidence(provenance, evidence):
    require(
        type(provenance) is str and bool(provenance.strip()),
        "Explicit provenance required",
    )
    require(
        type(evidence) in (list, tuple)
        and bool(evidence)
        and all(type(e) is str and bool(e.strip()) for e in evidence),
        "Nonempty evidence references required",
    )
    return list(evidence)


def _mapping(values, ids, kind):
    require(
        type(values) is dict and all(type(i) is int for i in values),
        f"Exact integer {kind} mapping required",
    )
    require(set(values) == set(ids), f"Exact site coverage required for {kind}")


def typed_data(graph, source, types, provenance, evidence, validation_policy):
    source.require_assignment()
    require(source.identity["sha256"] == PIN["sha256"], "Typed graph source mismatch")
    require(validation_policy == VALIDATION_POLICY, "Unsupported validation policy")
    evidence = _evidence(provenance, evidence)
    graph_system(graph)
    require(
        graph["representation"] == "atomistic" and not graph["has_box"],
        "Finite atomistic graph required",
    )
    atoms, adj, bonds, env = environments(graph)
    _mapping(types, atoms, "types")
    require(all(type(t) is str for t in types.values()), "String type labels required")
    for b in graph["bonds"]:
        require(
            b["aromatic"] == (b["order"] == 1.5)
            and (
                not b["aromatic"]
                or all(atoms[i]["metadata"]["aromatic"] for i in b["sites"])
            ),
            f"Inconsistent aromatic bond semantics at {b['sites']}",
        )
    entries = {}
    atom_rows = records(source.inventory, "atom_types")
    for i, label in types.items():
        row = select([r for r in atom_rows if r["data"]["type"] == label])
        require(row is not None, f"Unknown source label at site {i}: {label}")
        a, e = atoms[i], env[i]
        data = row["record"]["data"]
        require(
            data["element"] == a["element"], f"Element conflict at site {i}: {label}"
        )
        require(
            data["connections"] == e["degree"],
            f"Source connectivity conflict at site {i}: {label}",
        )
        require(
            label in SUPPORTED,
            f"Structural constraints not established at site {i}: {label}",
        )
        require(
            a["atomic_number"] == NUMBERS.get(a["element"])
            and abs(a["mass"] - MASSES[a["element"]]) <= 0.02
            and not a["metadata"]["isotope"]
            and not a["metadata"]["radical_electrons"],
            f"Element/mass/isotope/radical constraint at site {i}: {label}",
        )
        require(
            a["formal_charge"] == (1 if label in {"n+", "n4"} else 0),
            f"Charge-state conflict at site {i}: {label}",
        )
        ns = [atoms[j]["element"] for j in adj[i]]
        orders = e["bond_orders"]
        ok = False
        if label in {"c", "c1", "c2", "c3"}:
            ok = orders == [1, 1, 1, 1] and not e["aromatic"]
            if label != "c":
                ok = ok and e["hydrogens"] == int(label[1])
        elif label in {"h", "hc", "ho", "hn", "hn2"}:
            parents = {
                "h": {"C", "Si", "H"},
                "hc": {"C"},
                "ho": {"O"},
                "hn": {"N"},
                "hn2": {"N"},
            }[label]
            ok = orders == [1] and set(ns) <= parents and not e["aromatic"]
        elif label in {"o", "oc", "oh"}:
            ok = orders == [1, 1] and not e["aromatic"]
            if label == "oc":
                ok = (
                    ok
                    and sorted(ns) == ["C", "C"]
                    and all(
                        not any(
                            bonds[frozenset((j, k))]["order"] == 2
                            and atoms[k]["element"] == "O"
                            for k in adj[j]
                        )
                        for j in adj[i]
                    )
                )
            elif label == "oh":
                ok = ok and sorted(ns) == ["C", "H"]
        elif label == "cp":
            ok = e["aromatic"] and orders in ([1, 1.5, 1.5], [1.5, 1.5, 1.5])
        elif label == "ct":
            ok = not e["aromatic"] and orders == [1, 3]
        elif label == "nt":
            ok = orders == [3] and ns == ["C"] and not e["aromatic"]
        elif label in {"na", "n+", "n4"}:
            ok = orders == [1] * e["degree"] and set(ns) <= {"C", "H"}
            ok = (
                ok
                and not e["aromatic"]
                and all(
                    not atoms[j]["metadata"]["aromatic"]
                    and not any(bonds[frozenset((j, k))]["order"] != 1 for k in adj[j])
                    for j in adj[i]
                )
            )
        require(ok, f"Chemical/bond-order constraint conflict at site {i}: {label}")
        entries[i] = {
            "type": label,
            "environment": e,
            "source_atom_record": row,
            "structural_consistency": "validated",
        }
    return {
        "schema": TYPING_SCHEMA,
        "source": source.identity,
        "graph": graph,
        "graph_identity": identity(graph),
        "profile": {"name": PROFILE_NAME, "validation_policy": validation_policy},
        "assignments": dict(types),
        "entries": entries,
        "origin": "externally_established",
        "provenance": provenance,
        "evidence_references": evidence,
        "validation_policy": validation_policy,
        "validation": {
            "structural_consistency": "validated",
            "external_chemical_authority": "caller_asserted",
            "automatic_perception": "not_performed",
            "unresolved_chemical_interpretation": [
                "External evidence authority and scientific suitability not certified"
            ],
        },
        "coverage": {"total": len(atoms), "typed": len(types), "complete": True},
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFTypedGraph:
    """Immutable source-bound external types; payload properties return copies."""

    json_text: str
    source: PCFFSource

    @boundary
    @validated_record
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        expected = typed_data(
            p["graph"],
            self.source,
            p["assignments"],
            p["provenance"],
            p["evidence_references"],
            p["validation_policy"],
        )
        require(pack(p) == pack(expected), "Contradictory external typed graph")
        if system is not None:
            require(chemical_graph(system) == p["graph"], "Typed graph mismatch")

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return record_identity(self)

    @property
    def complete(self):
        return self.payload["coverage"]["complete"]

    @property
    def assignments(self):
        return self.payload["assignments"]


@boundary
def bind_pcff_types(
    system,
    source,
    types,
    *,
    provenance,
    evidence_references,
    validation_policy=VALIDATION_POLICY,
):
    """Bind external types without requiring automatic perception or agreement."""
    result = PCFFTypedGraph(
        pack(
            typed_data(
                chemical_graph(system),
                source,
                types,
                provenance,
                evidence_references,
                validation_policy,
            )
        ),
        source,
    )
    result.validate_integrity(system)
    return result


def charge_data(typed, source, origin, policy, supplied=None):
    from .expanded import native_increment_data
    from .fallbacks import validate_policy

    if policy is not None:
        validate_policy(policy)
    if origin == "native_increments":
        require(supplied is None, "Native mode cannot contain provided charges")
        data = native_increment_data(typed, source, resolution_policy=policy)
    else:
        require(origin == "provided", "Unknown charge origin")
        require(type(supplied) is dict, "Provided charge input required")
        require(
            set(supplied)
            == {
                "partial_charges",
                "unit",
                "provenance",
                "evidence_references",
                "component_totals",
                "total_charge",
            },
            "Malformed provided charge input",
        )
        evidence = _evidence(supplied["provenance"], supplied["evidence_references"])
        require(supplied["unit"] == "elementary_charge", "Invalid charge unit")
        atoms, adj, _, _ = environments(typed["graph"])
        values = supplied["partial_charges"]
        _mapping(values, atoms, "charges")
        require(
            all(type(q) in (int, float) and isfinite(q) for q in values.values()),
            "Finite nonboolean numeric charges required",
        )
        totals = supplied["component_totals"]
        components = _components(adj)
        _mapping(totals, [c[0] for c in components], "component totals")
        require(
            all(type(q) in (int, float) and isfinite(q) for q in totals.values())
            and type(supplied["total_charge"]) in (int, float)
            and isfinite(supplied["total_charge"]),
            "Finite charge totals required",
        )
        checked = []
        for c in components:
            formal = sum(atoms[i]["formal_charge"] for i in c)
            total = fsum(values[i] for i in c)
            require(
                totals[c[0]] == formal,
                f"Declared component/formal charge mismatch: {c}",
            )
            require(
                abs(total - totals[c[0]]) <= PIN["tolerance_e"],
                f"Provided component charge mismatch: {c}",
            )
            checked.append({"sites": c, "formal_charge": formal, "charge": total})
        formal = sum(a["formal_charge"] for a in atoms.values())
        require(
            supplied["total_charge"] == formal
            and abs(fsum(values.values()) - formal) <= PIN["tolerance_e"],
            "Provided system charge mismatch",
        )
        data = {
            "partial_charges": dict(values),
            "unit": "elementary_charge",
            "provenance": supplied["provenance"],
            "evidence_references": evidence,
            "components": checked,
            "total_charge": fsum(values.values()),
            "formal_charge": formal,
            "tolerance_e": PIN["tolerance_e"],
            "complete": True,
            "diagnostics": [],
            "native_increment_availability": "not_evaluated",
            "validation": "coverage_finiteness_and_declared_formal_totals",
            "chemical_authority": "caller_asserted",
            **FLAGS,
        }
    return {
        "schema": CHARGE_SCHEMA,
        "source": source.identity,
        "typed_graph": typed,
        "typing_identity": identity(typed),
        "origin": origin,
        "resolution_policy": policy,
        "provided_input": supplied,
        "charge_data": data,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFGraphCharges:
    """Native or provided charges, with separate truthful origins."""

    json_text: str
    source: PCFFSource

    @boundary
    @validated_record
    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        typed = PCFFTypedGraph(pack(p["typed_graph"]), self.source)
        typed.validate_integrity(system)
        expected = charge_data(
            typed.payload,
            self.source,
            p["origin"],
            p["resolution_policy"],
            p["provided_input"],
        )
        require(pack(p) == pack(expected), "Contradictory typed graph charges")

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return record_identity(self)

    @property
    def complete(self):
        return self.payload["charge_data"]["complete"]

    @property
    def charges(self):
        p = self.payload
        require(p["charge_data"]["complete"], "Incomplete typed graph charges")
        return p["charge_data"]["partial_charges"]


def _charges(system, typing, origin, policy, supplied=None):
    require(type(typing) is PCFFTypedGraph, "External PCFFTypedGraph required")
    typing.validate_integrity(system)
    result = PCFFGraphCharges(
        pack(charge_data(typing.payload, typing.source, origin, policy, supplied)),
        typing.source,
    )
    result.validate_integrity(system)
    return result


@boundary
def assign_typed_pcff_charges(system, typing, *, resolution_policy=None):
    """Resolve native increments from validated types, without perception."""
    return _charges(system, typing, "native_increments", resolution_policy)


@boundary
def provide_pcff_charges(
    system,
    typing,
    charges,
    *,
    unit,
    provenance,
    evidence_references,
    component_totals,
    total_charge,
    resolution_policy=None,
):
    """Preserve charges exactly; totals are keyed by each component's lowest ID."""
    return _charges(
        system,
        typing,
        "provided",
        resolution_policy,
        {
            "partial_charges": charges,
            "unit": unit,
            "provenance": provenance,
            "evidence_references": evidence_references,
            "component_totals": component_totals,
            "total_charge": total_charge,
        },
    )


def assignment_inputs(charge):
    """Read compatible records without relabelling their origin or serializing adapters."""
    if charge["schema"] == CHARGE_SCHEMA:
        return charge["typed_graph"], charge["charge_data"], charge["typing_identity"]
    return (
        charge["automatic_typing"],
        charge["native_charge_record"],
        charge["automatic_typing_identity"],
    )


def charge_result(payload, source):
    from .automatic import PCFFAutomaticChargeResult

    cls = (
        PCFFGraphCharges
        if payload["schema"] == CHARGE_SCHEMA
        else PCFFAutomaticChargeResult
    )
    return cls(pack(payload), source)


@boundary
def save_pcff_graph_record(result, path):
    require(
        type(result) in (PCFFTypedGraph, PCFFGraphCharges), "PCFF graph record required"
    )
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_graph_record(path, source, *, system=None):
    raw = Path(path).read_text()
    cls = {TYPING_SCHEMA: PCFFTypedGraph, CHARGE_SCHEMA: PCFFGraphCharges}.get(
        unpack(raw)["schema"]
    )
    require(cls is not None, "Unsupported PCFF graph schema")
    result = cls(raw, source)
    result.validate_integrity(system)
    return result
