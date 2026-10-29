"""Coordinate-independent final graph and exact unresolved dependency inventory."""

import json
from copy import deepcopy
from dataclasses import asdict
from hashlib import sha256


def canonical(value):
    # Normalize integer-keyed source maps before sorting; JSON reloads their
    # keys as strings, so signing the pre-normalized map is not round-trip safe.
    normalized = json.loads(json.dumps(value, allow_nan=False))
    return json.dumps(
        normalized, sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def digest(value):
    return sha256(canonical(value).encode()).hexdigest()


def final_graph(system):
    from island.chemistry.rdkit_graph import system_to_rdkit_graph
    from island.forcefields.opls_external.graph import chemical_metadata

    system.validate()
    system.topology.validate_bond_graph()
    system_to_rdkit_graph(system)
    sites = [asdict(a) for a in system.topology.sites.values()]
    for a in sites:
        a["metadata"] = chemical_metadata(a["metadata"])
    if not sites or any(type(a["id"]) is not int for a in sites):
        raise ValueError("Nonempty graph with stable integer atom IDs required")
    return {
        "representation": system.representation,
        "atom_order": list(system.topology.sites),
        "sites": sites,
        "bonds": [asdict(b) for _, b in sorted(system.topology.bonds.items())],
        "impropers": [asdict(i) for i in system.topology.impropers],
        "components": sorted(sorted(c) for c in system.topology.connected_components()),
        "molecules": [
            [i, a.metadata.get("chain_id")] for i, a in system.topology.sites.items()
        ],
        "system_metadata": chemical_metadata(deepcopy(system.metadata)),
        "box": asdict(system.box) if system.box is not None else None,
    }


def requirements(system, family):
    from island.forcefields.oplsaa.parameters import inventory
    from island.forcefields.pcff.automatic import chemical_graph
    from island.forcefields.pcff.class2 import interaction_requests

    rows = []

    def add(kind, sites, dependencies=(), applicability=None):
        rows.append(
            {
                "id": key(kind, sites),
                "family": kind,
                "sites": list(sites),
                "dependencies": list(dependencies),
                "applicability": applicability,
            }
        )

    for i in sorted(system.topology.sites):
        add("atom_type", [i])
        add("charge", [i], [key("atom_type", [i])])
        add("mass", [i], [key("atom_type", [i])])
    if family == "pcff":
        requests, _ = interaction_requests(chemical_graph(system), expanded=True)
        for kind, sites, deps, reason in requests:
            add(kind, sites, [key(k, s) for k, s in deps], reason)
    else:
        paths, pairs = inventory(system)
        for kind, paths_for_kind in paths.items():
            for sites in paths_for_kind:
                add({"bonds": "bond", "angles": "angle", "rb": "proper"}[kind], sites)
        for i in sorted(system.topology.sites):
            add("nonbonded", [i])
        for p in pairs:
            add("special_pair", p["sites"])
        for item in system.topology.impropers:
            add("improper", [item.site1, item.site2, item.site3, item.site4])
    for kind in ("mixing_rule", "scaling_14", "exclusions", "electrostatics"):
        add(kind, [])
    return sorted(rows, key=lambda r: r["id"])


def key(kind, sites):
    return kind + ":" + ",".join(map(str, sites))
