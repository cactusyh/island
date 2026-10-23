"""Exact-graph OPLS-compatible terms with explicitly external RESP/Psi4 charges.

Native OPLS charge assignment is untouched. This narrow reviewed profile accepts
only the aromatic benchmark's authoritative C3 native binary64 RESP vector.
It resolves every bonded/LJ term from the original pinned XML, never the charges.
"""

import json
from dataclasses import dataclass
from hashlib import sha256
from math import fsum
from xml.etree import ElementTree as ET

from island.forcefields.ambertools.resp_policy import FiniteRESPCharges
from island.forcefields.oplsaa.adapter import type_atoms
from island.forcefields.oplsaa.models import OPLSTypingResult
from island.forcefields.oplsaa.parameters import inventory, selected_record
from island.forcefields.oplsaa.source import FoyerOPLSSource, require
from island.graph import (
    FinalChemicalGraph,
    unified_assignment_diagnostics,
    unified_graph_charges,
    unified_parameter_assignment,
    unified_typed_graph,
)
from island.graph.final import _system_data, _system_from_payload, pack, unpack
from island.periodic import (
    BACKEND_SCHEMA,
    FLAGS,
    PeriodicForceFieldConfig,
    _signed,
    periodic_chemical_view,
)

SCHEMA = "island.opls-compatible-resp.v1"
GRAPH = "f7730225baafb0a2164cdcf3db2f8813d83d3b8d145923053b0fc4ffd193f052"
CHARGES = "0fcbcbcab2d8cb73a89f943ae4bf5b0daade08a480eb24d232e982d8d7315bf3"


def numerical(system, source, typing, charges):
    """Resolve source terms without invoking or repairing native charge rows."""
    types, _ = source.entries
    assigned = typing.payload["data"]["assignments"]
    paths, pairs = inventory(system)
    root = ET.fromstring(source.xml)
    data = {
        "family": "OPLS-AA",
        "sites": [],
        "bonds": [],
        "angles": [],
        "torsions": [],
        "impropers": [],
        "pairs": [],
    }
    witness = {}
    for kind, paths_for_kind in paths.items():
        witness[kind] = []
        for ids in paths_for_kind:
            names = tuple(assigned[i] for i in ids)
            index, attrs, candidates = selected_record(root, types, names, kind)
            witness[kind].append(
                {
                    "sites": ids,
                    "index": index,
                    "attributes": attrs,
                    "candidates": candidates,
                }
            )
            row = {"sites": ids}
            if kind == "bonds":
                row.update(
                    length=float(attrs["length"]) * 10, k=float(attrs["k"]) / 100
                )
            elif kind == "angles":
                row.update(theta=float(attrs["angle"]), k=float(attrs["k"]))
            else:
                row.update(rb=[float(attrs["c" + str(i)]) for i in range(6)])
            data[{"bonds": "bonds", "angles": "angles", "rb": "torsions"}[kind]].append(
                row
            )
    nb = {r.attrib["type"]: r.attrib for r in root.find("NonbondedForce")}
    for i in sorted(assigned):
        t = assigned[i]
        a = nb[t]
        data["sites"].append(
            {
                "id": i,
                "mass": system.topology.sites[i].mass,
                "charge": charges[i],
                "sigma": float(a["sigma"]) * 10,
                "epsilon": float(a["epsilon"]),
                "atom_type": t,
            }
        )
    data["pairs"] = [
        {"sites": r["sites"], "lj": r["lj_scale"], "coulomb": r["coulomb_scale"]}
        for r in pairs
    ]
    return data, witness


def definition(system, graph, record):
    graph.validate_integrity(system)
    require(
        set(record)
        == {
            "schema",
            "charge_origin",
            "native_opls_charges",
            "label",
            "source_xml",
            "typing_json",
            "resp_charges",
        },
        "Unexpected provenance fields",
    )
    require(
        record["schema"] == SCHEMA and graph.identity == GRAPH,
        "Unsupported exact-graph RESP profile",
    )
    require(
        record["charge_origin"] == "RESP/Psi4"
        and record["native_opls_charges"] is False,
        "Misrepresented charge origin",
    )
    require(
        record["label"] == "OPLS-compatible parameters with RESP/Psi4 charges",
        "Misrepresented model label",
    )
    source = FoyerOPLSSource(record["source_xml"].encode())
    source.validate_integrity()
    chemical = periodic_chemical_view(system)
    typing = OPLSTypingResult(record["typing_json"])
    typing.validate_integrity(chemical, source)
    # Exact graph-specific review also fixes the independent Foyer assignments.
    require(
        typing.payload["data"]["assignments"]
        == {
            s["id"]: "opls_145" if s["element"] == "C" else "opls_146"
            for s in graph.payload["sites"]
        },
        "Changed reviewed aromatic typing",
    )
    provided = FiniteRESPCharges(json.dumps(record["resp_charges"], allow_nan=False))
    q = provided.validate(graph, expected_reference_identity=CHARGES)
    source_id = {
        **source.identity,
        "charge_policy": "RESP/Psi4",
        "resp_charge_identity": CHARGES,
    }
    source_id["sha256"] = sha256(pack(source_id).encode()).hexdigest()
    data, witness = numerical(chemical, source, typing, q)
    config = PeriodicForceFieldConfig(
        family="OPLS-AA",
        source=source_id,
        box_lengths=system.box.lengths,
        nonbonded_cutoff=5.0,
        lj_mixing_rule="geometric_12_6",
        one_four_scaling=(0.5, 0.5),
        electrostatics_method="ewald",
        pme_tolerance=1e-5,
        provenance=record["label"],
        evidence=(graph.identity, CHARGES),
    )
    evidence = [CHARGES, typing.identity]
    typed = unified_typed_graph(
        graph,
        force_field="OPLS-AA",
        atom_types=typing.payload["data"]["assignments"],
        source=source_id,
        typing_profile=SCHEMA,
        typing_method="pinned_foyer",
        automatic_perception="performed",
        provenance=record["label"],
        evidence=evidence,
    )
    charge = unified_graph_charges(
        graph,
        typed,
        force_field="OPLS-AA",
        charges=q,
        charge_method="RESP",
        charge_origin="RESP/Psi4",
        source=source_id,
        component_totals={
            str(c[0]): fsum(q[i] for i in c) for c in graph.payload["components"]
        },
        total_charge=fsum(q.values()),
        provenance=record["label"],
        evidence=evidence,
    )
    diagnostics = unified_assignment_diagnostics(
        typing_complete=True,
        charges_complete=True,
        native_charge_available=False,
        parameter_rows_complete=True,
        executable_potential_complete=True,
        periodic_backend_support=True,
        reasons=[],
    )
    assignment = unified_parameter_assignment(
        graph,
        typed,
        charge,
        force_field="OPLS-AA",
        source=source_id,
        native_payload={"record": record, "source_parameter_rows": witness},
        diagnostics=diagnostics,
    )
    return unpack(
        _signed(
            dict(
                schema=BACKEND_SCHEMA,
                system=_system_data(system),
                final_graph=graph.payload,
                config=config.payload,
                typed_graph=typed.payload,
                graph_charges=charge.payload,
                parameter_assignment=assignment.payload,
                numerical_data=data,
                native_preparation_identity=assignment.identity,
                final_graph_identity=graph.identity,
                typed_graph_identity=typed.identity,
                charge_identity=charge.identity,
                assignment_identity=assignment.identity,
                config_identity=config.identity,
                **FLAGS,
            )
        )
    )


def validate_backend(payload):
    record = payload["parameter_assignment"]["native_payload"]["record"]
    expected = definition(
        _system_from_payload(payload["system"]),
        FinalChemicalGraph(pack(payload["final_graph"])),
        record,
    )
    require(
        expected == payload,
        "OPLS/RESP backend contradicts source, graph, charges or provenance",
    )


@dataclass(frozen=True)
class OPLSRESPBackend:
    json_text: str

    @property
    def payload(self):
        return json.loads(self.json_text)

    def validate_integrity(self):
        validate_backend(self.payload)

    @property
    def system(self):
        return _system_from_payload(self.payload["system"])


def prepare(system, graph, source, provided):
    """Create a separately named model; never reports native OPLS/AM1-BCC success."""
    typing = type_atoms(periodic_chemical_view(system), source)
    record = {
        "schema": SCHEMA,
        "charge_origin": "RESP/Psi4",
        "native_opls_charges": False,
        "label": "OPLS-compatible parameters with RESP/Psi4 charges",
        "source_xml": source.xml.decode(),
        "typing_json": typing.json_text,
        "resp_charges": provided.payload,
    }
    return OPLSRESPBackend(
        json.dumps(definition(system, graph, record), allow_nan=False)
    )
