"""Force-field-neutral immutable final chemical graphs and transformations.

The graph identity excludes coordinates. Coordinates remain in the system bundle
and are checked separately by evaluators. This module performs no force-field
typing, charge fitting, parameter lookup, or scientific readiness attestation.
"""

import json
from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
from math import fsum, isfinite
from pathlib import Path
from uuid import uuid4

from island.core import AtomSite, Coordinates, MolecularSystem, SimulationBox, Topology
from island.exceptions import ValidationError

SCHEMA = "island_final_chemical_graph_v1"
TRANSFORMATION_SCHEMA = "island_final_graph_transformation_v1"
CHARGE_SCHEMA = "island_force_neutral_charges_v1"
BUNDLE_SCHEMA = "island_final_graph_bundle_v1"


def _require(ok, message):
    if not ok:
        raise ValidationError(message)


def _json_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def pack(value):
    return json.dumps(
        {"payload": value, "sha256": sha256(_json_bytes(value)).hexdigest()},
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def unpack(raw):
    envelope = json.loads(raw)
    _require(
        type(envelope) is dict and set(envelope) == {"payload", "sha256"},
        "Malformed record envelope",
    )
    _require(
        envelope["sha256"] == sha256(_json_bytes(envelope["payload"])).hexdigest(),
        "Record checksum changed",
    )
    return envelope["payload"]


def _digest(value):
    return sha256(pack(value).encode()).hexdigest()


def _checksum(raw):
    return sha256(raw).hexdigest()


def _system_data(system):
    payload = system.to_dict()
    if system.topology.impropers:
        payload["impropers"] = [
            {name: getattr(item, name) for name in item.__dataclass_fields__}
            for item in system.topology.impropers
        ]
    return payload


def _finite(value, label):
    _require(
        type(value) in (int, float) and not isinstance(value, bool) and isfinite(value),
        f"{label} must be finite",
    )


def _valence(element, value, aromatic=False):
    allowed = {
        "H": (1,),
        "C": (2, 3, 4),
        "N": (1, 2, 3, 4, 5),
        "O": (1, 2),
        "F": (1,),
        "P": (3, 5),
        "S": (2, 4, 6),
        "Cl": (1,),
        "Br": (1,),
        "I": (1,),
    }
    return value in allowed.get(element, ()) or (
        aromatic and element in {"C", "N"} and value in {3, 4, 4.5}
    )


def _graph_payload(system, *, transformations=(), molecule_membership=None):
    system.validate()
    _require(
        system.representation == "atomistic", "Final graph requires atomistic system"
    )
    ids = sorted(system.topology.sites)
    _require(
        ids and all(type(i) is int and not isinstance(i, bool) for i in ids),
        "Stable integer atom IDs required",
    )
    sites = []
    for i in ids:
        atom = system.topology.sites[i]
        _require(
            isinstance(atom, AtomSite) and atom.id == i, "AtomSite identity mismatch"
        )
        _require(
            type(atom.element) is str and atom.element, f"Invalid element at site {i}"
        )
        _require(
            type(atom.formal_charge) is int
            and not isinstance(atom.formal_charge, bool),
            f"Invalid formal charge at site {i}",
        )
        _require(type(atom.metadata) is dict, f"Invalid atom metadata at site {i}")
        aromatic = bool(atom.metadata.get("aromatic", False))
        sites.append(
            {
                "id": i,
                "element": atom.element,
                "atomic_number": atom.atomic_number,
                "mass": atom.mass,
                "formal_charge": atom.formal_charge,
                "aromatic": aromatic,
                "metadata": deepcopy(atom.metadata),
            }
        )
    bonds = []
    for key, bond in sorted(system.topology.bonds.items()):
        _require(key == bond.key and set(key) <= set(ids), f"Invalid bond IDs: {key}")
        _require(
            type(bond.order) in (int, float)
            and not isinstance(bond.order, bool)
            and isfinite(bond.order)
            and bond.order > 0,
            f"Invalid bond order: {key}",
        )
        _require(bond.aromatic == (bond.order == 1.5), f"Aromatic bond mismatch: {key}")
        bonds.append(
            {"sites": list(key), "order": bond.order, "aromatic": bond.aromatic}
        )
    adjacency = {i: set() for i in ids}
    for b in bonds:
        a, z = b["sites"]
        adjacency[a].add(z)
        adjacency[z].add(a)
    valences = {i: fsum(b["order"] for b in bonds if i in b["sites"]) for i in ids}
    for site, value in valences.items():
        _require(
            _valence(
                system.topology.sites[site].element,
                value,
                bool(system.topology.sites[site].metadata.get("aromatic")),
            ),
            f"Invalid valence at site {site}: {value}",
        )
    components = [sorted(c) for c in system.topology.connected_components()]
    components.sort(key=lambda c: c[0])
    component_of = {i: n for n, c in enumerate(components) for i in c}
    if molecule_membership is None:
        molecule_membership = {
            i: str(system.topology.sites[i].metadata.get("chain_id", component_of[i]))
            for i in ids
        }
    _require(
        type(molecule_membership) is dict and set(molecule_membership) == set(ids),
        "Exact molecule membership coverage required",
    )
    _require(
        all(type(v) is str and v for v in molecule_membership.values()),
        "Molecule membership labels must be nonempty text",
    )
    membership = {i: molecule_membership[i] for i in ids}
    box = None
    if system.box is not None:
        _require(type(system.box) is SimulationBox, "Unsupported simulation box")
        box = {
            "lengths": list(system.box.lengths),
            "periodic": list(system.box.periodic),
            "boundary": ["p" if x else "f" for x in system.box.periodic],
        }
        for value in box["lengths"]:
            _finite(value, "Box length")
    return {
        "schema": SCHEMA,
        "sites": sites,
        "bonds": bonds,
        "connectivity": {str(i): sorted(adjacency[i]) for i in ids},
        "components": components,
        "component_membership": {str(i): component_of[i] for i in ids},
        "molecule_membership": {str(i): membership[i] for i in ids},
        "repeat_provenance": deepcopy(system.metadata.get("polymer", {})),
        "crosslink_provenance": deepcopy(system.metadata.get("crosslinks", [])),
        "periodic": box,
        "transformations": deepcopy(list(transformations)),
    }


def _base_payload(payload):
    p = deepcopy(payload)
    p.pop("graph_identity", None)
    p.pop("history_identity", None)
    return p


def _chemical_payload(payload):
    p = _base_payload(payload)
    p["transformations"] = []
    return p


def _validate_payload(payload):
    _require(
        type(payload) is dict and payload.get("schema") == SCHEMA,
        "Unsupported final graph schema",
    )
    expected = {
        "schema",
        "sites",
        "bonds",
        "connectivity",
        "components",
        "component_membership",
        "molecule_membership",
        "repeat_provenance",
        "crosslink_provenance",
        "periodic",
        "transformations",
        "graph_identity",
        "history_identity",
    }
    _require(set(payload) == expected, "Malformed final graph fields")
    ids = [s.get("id") for s in payload["sites"]]
    _require(
        ids == sorted(ids)
        and len(ids) == len(set(ids))
        and all(type(i) is int for i in ids),
        "Final graph site IDs must be sorted unique integers",
    )
    idset = set(ids)
    bonds = payload["bonds"]
    seen = set()
    for b in bonds:
        _require(set(b) == {"sites", "order", "aromatic"}, "Malformed final graph bond")
        a, z = b["sites"]
        key = tuple(sorted((a, z)))
        _require(
            a != z and key not in seen and set(key) <= idset,
            "Invalid or duplicate final graph bond",
        )
        seen.add(key)
        _finite(b["order"], "Bond order")
        _require(
            b["order"] > 0 and b["aromatic"] == (b["order"] == 1.5),
            "Invalid aromatic bond",
        )
    _require(
        set(payload["connectivity"]) == {str(i) for i in ids},
        "Connectivity ID coverage mismatch",
    )
    for i in ids:
        neighbors = payload["connectivity"][str(i)]
        _require(
            neighbors == sorted(neighbors) and set(neighbors) <= idset,
            "Invalid connectivity",
        )
    _require(
        type(payload["periodic"]) in (dict, type(None)), "Malformed periodic metadata"
    )
    if payload["periodic"] is not None:
        q = payload["periodic"]
        _require(
            set(q) == {"lengths", "periodic", "boundary"}, "Malformed periodic box"
        )
        _require(
            len(q["lengths"]) == len(q["periodic"]) == len(q["boundary"]) == 3,
            "Periodic box dimensions required",
        )
        for value in q["lengths"]:
            _finite(value, "Box length")
            _require(value > 0, "Positive box lengths required")
        _require(
            all(type(v) is bool for v in q["periodic"])
            and q["boundary"] == ["p" if x else "f" for x in q["periodic"]],
            "Boundary metadata mismatch",
        )
    _require(
        payload["graph_identity"] == _digest(_chemical_payload(payload)),
        "Final graph identity mismatch",
    )
    _require(
        payload["history_identity"] == _digest(payload["transformations"]),
        "Final graph history identity mismatch",
    )
    for item in payload["transformations"]:
        transform = GraphTransformation(pack(item))
        transform.validate_integrity()
        _require(
            item["output_graph_identity"] == payload["graph_identity"]
            or item["output_graph_identity"]
            in {t["output_graph_identity"] for t in payload["transformations"]},
            "Transformation output is not bound to graph history",
        )


def graph_identity(payload):
    _validate_payload(payload)
    return payload["graph_identity"]


@dataclass(frozen=True)
class FinalChemicalGraph:
    """Immutable chemical/topological graph; coordinates are intentionally excluded."""

    json_text: str

    def validate_integrity(self, system=None):
        p = unpack(self.json_text)
        _validate_payload(p)
        if system is not None:
            actual = _graph_payload(
                system,
                transformations=p["transformations"],
                molecule_membership={
                    int(i): v for i, v in p["molecule_membership"].items()
                },
            )
            _require(
                _digest(_chemical_payload(actual)) == p["graph_identity"],
                "Final graph chemical/topology/box mismatch",
            )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return unpack(self.json_text)["graph_identity"]

    @property
    def graph_identity(self):
        return self.identity

    @property
    def periodic(self):
        return self.payload["periodic"] is not None


def final_graph(system, *, transformations=(), molecule_membership=None):
    p = _graph_payload(
        system, transformations=transformations, molecule_membership=molecule_membership
    )
    p["graph_identity"] = _digest(_chemical_payload(p))
    p["history_identity"] = _digest(p["transformations"])
    return FinalChemicalGraph(pack(p))


def prepare_final_graph(system, graph=None):
    """Validate and return the neutral final graph for a concrete system.

    This stage owns graph/history/box identity. It deliberately does not select
    a force field or claim that a periodic/crosslinked model is executable.
    """
    result = final_graph(system) if graph is None else graph
    result.validate_integrity(system)
    return result


@dataclass(frozen=True)
class GraphTransformation:
    json_text: str

    def validate_integrity(self, input_graph=None, output_graph=None):
        p = unpack(self.json_text)
        _require(
            set(p)
            == {
                "schema",
                "operation",
                "input_graph_identity",
                "output_graph_identity",
                "parameters",
                "seed",
                "provenance",
                "evidence",
                "identity",
            },
            "Malformed graph transformation",
        )
        _require(
            p["schema"] == TRANSFORMATION_SCHEMA
            and p["identity"] == _digest({k: p[k] for k in p if k != "identity"}),
            "Graph transformation identity mismatch",
        )
        if input_graph is not None:
            _require(
                p["input_graph_identity"] == input_graph.identity,
                "Transformation input mismatch",
            )
        if output_graph is not None:
            _require(
                p["output_graph_identity"] == output_graph.identity,
                "Transformation output mismatch",
            )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        self.validate_integrity()
        return unpack(self.json_text)["identity"]


def transformation(
    input_graph,
    output_graph,
    operation,
    *,
    parameters=None,
    seed=None,
    provenance,
    evidence=None,
):
    input_graph.validate_integrity()
    output_graph.validate_integrity()
    _require(type(operation) is str and operation, "Transformation operation required")
    _require(
        type(provenance) is str and provenance.strip(),
        "Transformation provenance required",
    )
    evidence = [] if evidence is None else list(evidence)
    _require(
        all(type(x) is str and x.strip() for x in evidence),
        "Invalid transformation evidence",
    )
    p = {
        "schema": TRANSFORMATION_SCHEMA,
        "operation": operation,
        "input_graph_identity": input_graph.identity,
        "output_graph_identity": output_graph.identity,
        "parameters": deepcopy(parameters or {}),
        "seed": seed,
        "provenance": provenance,
        "evidence": evidence,
    }
    p["identity"] = _digest(p)
    return GraphTransformation(pack(p))


def _system_from_payload(payload):
    _require(
        type(payload) is dict
        and set(payload)
        >= {"representation", "box", "sites", "bonds", "coordinates", "metadata"},
        "Malformed system bundle",
    )
    topology = Topology()
    for row in payload["sites"]:
        _require(
            row.get("kind") == "AtomSite", "Final graph bundle requires atom sites"
        )
        topology.add_site(AtomSite(**{k: v for k, v in row.items() if k != "kind"}))
    for row in payload["bonds"]:
        topology.add_bond(**row)
    box = payload["box"]
    simbox = (
        None
        if box is None
        else SimulationBox(*box["lengths"], periodic=tuple(box["periodic"]))
    )
    coordinates = {int(k): value for k, value in payload["coordinates"].items()}
    system = MolecularSystem(
        topology,
        Coordinates(coordinates),
        box=simbox,
        metadata=deepcopy(payload["metadata"]),
        representation=payload["representation"],
    )
    system.validate()
    return system


def _apply(system, operation, parameters, *, seed, provenance, evidence):
    prior = deepcopy(system.metadata.get("final_graph_transformations", []))
    before = final_graph(system, transformations=prior)
    result = system.copy()
    if operation == "crosslink_bonds":
        bonds = [tuple(x) for x in parameters["bonds"]]
        _require(
            len(bonds) == len({tuple(sorted(x)) for x in bonds}),
            "Duplicate crosslink bonds",
        )
        for a, z in bonds:
            _require(
                tuple(sorted((a, z))) not in result.topology.bonds,
                "Crosslink duplicates existing bond",
            )
            result.topology.add_bond(a, z, order=1)
        result.topology.rebuild_derived_interactions()
        result.metadata["crosslinks"] = deepcopy(
            result.metadata.get("crosslinks", [])
        ) + [{"sites": list(x), "seed": seed} for x in bonds]
    elif operation == "periodic_box":
        lengths = tuple(parameters["lengths"])
        periodic = tuple(parameters.get("periodic", (True, True, True)))
        result.box = SimulationBox(*lengths, periodic=periodic)
    elif operation == "end_group_completion":
        for a, z, order in parameters.get("bonds", []):
            result.topology.add_bond(a, z, order=order)
        result.topology.rebuild_derived_interactions()
    else:
        _require(False, f"Unsupported graph transformation: {operation}")
    after = final_graph(result)
    tr = transformation(
        before,
        after,
        operation,
        parameters=parameters,
        seed=seed,
        provenance=provenance,
        evidence=evidence,
    )
    # The final graph history binds the transformation itself, so build again with
    # the transformation payload and re-sign the graph.
    history = [tr.payload]
    graph = final_graph(result, transformations=prior + history)
    tr = transformation(
        before,
        graph,
        operation,
        parameters=parameters,
        seed=seed,
        provenance=provenance,
        evidence=evidence,
    )
    graph = final_graph(result, transformations=prior + [tr.payload])
    result.metadata["final_graph_transformations"] = deepcopy(
        graph.payload["transformations"]
    )
    return result, graph, tr


def build_psmiles_graph(psmiles, *, dp, random_seed=2026, **kwargs):
    from island.builders import build_linear_polymer

    system = build_linear_polymer(psmiles, dp=dp, random_seed=random_seed, **kwargs)
    graph = final_graph(system)
    p = {
        "operation": "psmiles_chain_generation",
        "psmiles": psmiles,
        "dp": dp,
        "random_seed": random_seed,
        "builder_kwargs": deepcopy(kwargs),
    }
    tr = transformation(
        final_graph(system),
        graph,
        "psmiles_chain_generation",
        parameters=p,
        seed=random_seed,
        provenance="ISLAND PSMILES linear polymer builder",
        evidence=["builder API and final graph validation"],
    )
    graph = final_graph(system, transformations=[tr.payload])
    tr = transformation(
        final_graph(system),
        graph,
        "psmiles_chain_generation",
        parameters=p,
        seed=random_seed,
        provenance="ISLAND PSMILES linear polymer builder",
        evidence=["builder API and final graph validation"],
    )
    graph = final_graph(system, transformations=[tr.payload])
    system.metadata["final_graph_transformations"] = deepcopy(
        graph.payload["transformations"]
    )
    return system, graph, tr


def complete_end_groups(
    system,
    *,
    bonds=(),
    seed=2026,
    provenance="Explicit end-group completion",
    evidence=None,
):
    return _apply(
        system,
        "end_group_completion",
        {"bonds": [list(x) for x in bonds]},
        seed=seed,
        provenance=provenance,
        evidence=evidence,
    )


def create_crosslinks(
    system,
    bonds,
    *,
    seed=2026,
    provenance="Explicit controlled crosslink fixture",
    evidence=None,
):
    return _apply(
        system,
        "crosslink_bonds",
        {"bonds": [list(x) for x in bonds]},
        seed=seed,
        provenance=provenance,
        evidence=evidence,
    )


def construct_periodic_box(
    system,
    lengths,
    *,
    periodic=(True, True, True),
    seed=2026,
    provenance="Deterministic orthorhombic box construction",
    evidence=None,
):
    return _apply(
        system,
        "periodic_box",
        {"lengths": list(lengths), "periodic": list(periodic)},
        seed=seed,
        provenance=provenance,
        evidence=evidence,
    )


@dataclass(frozen=True)
class FinalGraphChargeRecord:
    json_text: str

    def validate_integrity(self, graph, *, system=None):
        graph.validate_integrity(system)
        p = unpack(self.json_text)
        _require(
            p["schema"] == CHARGE_SCHEMA and p["graph_identity"] == graph.identity,
            "Charge graph identity mismatch",
        )
        ids = {s["id"] for s in graph.payload["sites"]}
        _require(
            set(p["charges"]) == {str(i) for i in ids}
            and all(type(i) is str and i.lstrip("-").isdigit() for i in p["charges"]),
            "Exact charge site coverage required",
        )
        _require(p["unit"] == "elementary_charge", "Unsupported charge unit")
        for q in p["charges"].values():
            _finite(q, "Charge")
        _require(
            type(p["source"]) is dict and p["source"].get("sha256"),
            "Charge source identity required",
        )
        _require(type(p["method"]) is str and p["method"], "Charge method required")
        _require(type(p["origin"]) is str and p["origin"], "Charge origin required")
        _require(
            type(p["provenance"]) is str and p["provenance"].strip(),
            "Charge provenance required",
        )
        _require(
            p["identity"] == _digest({k: p[k] for k in p if k != "identity"}),
            "Charge identity mismatch",
        )

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    @property
    def charges(self):
        return {int(i): value for i, value in self.payload["charges"].items()}


def force_neutral_charges(
    graph,
    charges,
    *,
    force_field,
    source,
    method,
    origin,
    provenance,
    evidence=None,
    unit="elementary_charge",
):
    graph.validate_integrity()
    _require(type(charges) is dict, "Charge mapping required")
    ids = {s["id"] for s in graph.payload["sites"]}
    _require(set(charges) == ids, "Exact charge coverage required")
    _require(type(source) is dict and source.get("sha256"), "Source identity required")
    for q in charges.values():
        _finite(q, "Charge")
    p = {
        "schema": CHARGE_SCHEMA,
        "graph_identity": graph.identity,
        "force_field": force_field,
        "source": deepcopy(source),
        "method": method,
        "origin": origin,
        "charges": {str(i): value for i, value in charges.items()},
        "unit": unit,
        "provenance": provenance,
        "evidence": list(evidence or []),
        "formal_charge_total": sum(s["formal_charge"] for s in graph.payload["sites"]),
        "total_charge": fsum(charges.values()),
    }
    p["identity"] = _digest(p)
    result = FinalGraphChargeRecord(pack(p))
    result.validate_integrity(graph)
    return result


def save_final_graph_bundle(system, graph, path, *, transformations=(), charges=None):
    graph.validate_integrity(system)
    target = Path(path)
    _require(
        not target.exists() and not target.is_symlink(),
        "Bundle destination already exists",
    )
    _require(target.parent.is_dir(), "Bundle parent must exist")
    raw = {
        "system.json": _json_bytes(_system_data(system)),
        "final-graph.json": graph.json_text.encode(),
    }
    if transformations:
        raw["transformations.json"] = _json_bytes(
            [
                t.payload if isinstance(t, GraphTransformation) else t
                for t in transformations
            ]
        )
    if charges is not None:
        charges.validate_integrity(graph, system=system)
        raw["charges.json"] = charges.json_text.encode()
    manifest = {
        "schema": BUNDLE_SCHEMA,
        "graph_identity": graph.identity,
        "history_identity": graph.payload["history_identity"],
        "files": {n: {"sha256": _checksum(v)} for n, v in raw.items()},
    }
    staging = target.parent / ("." + target.name + "." + uuid4().hex)
    staging.mkdir()
    try:
        for n, v in raw.items():
            (staging / n).write_bytes(v)
        (staging / "manifest.json").write_bytes(_json_bytes(manifest))
        target.mkdir()
        for child in staging.iterdir():
            child.replace(target / child.name)
        staging.rmdir()
    except Exception:
        if staging.exists():
            for child in staging.iterdir():
                child.unlink()
            staging.rmdir()
        raise
    return target


def load_final_graph_bundle(path):
    root = Path(path)
    manifest = json.loads((root / "manifest.json").read_text())
    _require(manifest.get("schema") == BUNDLE_SCHEMA, "Unsupported final graph bundle")
    raw = {}
    for name, entry in manifest["files"].items():
        data = (root / name).read_bytes()
        _require(
            _checksum(data) == entry["sha256"],
            f"Bundle checksum mismatch: {name}",
        )
        raw[name] = data
    system = _system_from_payload(json.loads((root / "system.json").read_text()))
    graph = FinalChemicalGraph(raw["final-graph.json"])
    graph.validate_integrity(system)
    _require(
        graph.identity == manifest["graph_identity"]
        and graph.payload["history_identity"] == manifest["history_identity"],
        "Bundle graph/history identity mismatch",
    )
    charges = None
    if "charges.json" in raw:
        charges = FinalGraphChargeRecord(raw["charges.json"])
        charges.validate_integrity(graph, system=system)
    return system, graph, charges


# Stable descriptive aliases for the public contract.
FinalGraph = FinalChemicalGraph
FinalGraphTransformation = GraphTransformation
create_final_graph = final_graph
build_final_graph_from_psmiles = build_psmiles_graph
add_crosslinks = create_crosslinks
build_periodic_box = construct_periodic_box
