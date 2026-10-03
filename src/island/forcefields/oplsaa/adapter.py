"""Foyer's pinned SMARTS engine, without full force-field application."""

import importlib.util
from importlib.metadata import version
from pathlib import Path

from island.charge_references.records import pack
from island.exceptions import OPLSDependencyError, OPLSTypingError
from island.forcefields.oplsaa.models import (
    OPLSTypingResult,
    chemical_graph,
    typing_data,
)
from island.forcefields.oplsaa.source import PIN, boundary, digest, require


@boundary
def check_foyer_installation():
    spec = importlib.util.find_spec("foyer")
    if spec is None or spec.origin is None:
        raise OPLSDependencyError("Install the pinned optional Foyer environment")
    root = Path(spec.origin).parent
    for relative, sha in PIN["implementation_files"].items():
        path = root / relative
        if not path.is_file() or digest(path.read_bytes()) != sha:
            raise OPLSDependencyError(
                f"Pinned Foyer implementation mismatch: {relative}"
            )
    try:
        from foyer.atomtyper import AtomTypingRulesProvider, find_atomtypes
        from foyer.topology_graph import TopologyGraph

        env = {
            n: version(n)
            for n in (
                "foyer",
                "gmso",
                "networkx",
                "lark",
                "lxml",
                "numpy",
                "parmed",
                "openmm",
                "ele",
            )
        }
    except Exception as error:
        raise OPLSDependencyError(f"Foyer dependency import failed: {error}") from error
    return AtomTypingRulesProvider, find_atomtypes, TopologyGraph, env


@boundary
def type_atoms(system, source):
    graph = chemical_graph(system)
    types, _ = source.entries  # Check source before importing optional dependencies.
    Provider, find, Graph, environment = check_foyer_installation()
    provider = Provider(
        {name: row["def"] for name, row in types.items() if row.get("def")},
        {
            name: {x.strip() for x in row.get("overrides", "").split(",") if x.strip()}
            for name, row in types.items()
        },
        set(),
    )
    external = Graph()
    index_to_site = {i: atom["id"] for i, atom in enumerate(graph["sites"])}
    site_to_index = {sid: i for i, sid in index_to_site.items()}
    for i, atom in enumerate(graph["sites"]):
        external.add_atom(
            i,
            atom["element"],
            atomic_number=atom["atomic_number"],
            symbol=atom["element"],
        )
    for bond in graph["bonds"]:
        a, b = bond["sites"]
        external.add_bond(site_to_index[a], site_to_index[b], bond_order=bond["order"])
    try:
        # No residue-name cache, Forcefield.apply, OpenMM Context or coordinates.
        # Monotonic rule labels admit at most (#sites * #rules) additions.
        matches = find(external, provider, max_iter=len(types) * len(index_to_site) + 1)
    except Exception as error:
        raise OPLSTypingError(
            f"Foyer untyped/ambiguous/unsupported graph: {error}; "
            f"index_to_site_id={index_to_site}",
            index_to_site_id=index_to_site,
            upstream_error=error,
        ) from error
    require(set(matches) == set(index_to_site), "Upstream typing coverage changed")
    rows = {
        index_to_site[i]: {
            "atomtype": row["atomtype"],
            "whitelist": sorted(row["whitelist"]),
            "blacklist": sorted(row["blacklist"]),
        }
        for i, row in matches.items()
    }
    payload = {
        "schema": "island_foyer_typing_v1",
        "source": source.identity,
        "graph": graph,
        "index_to_site_id": index_to_site,
        "matches": rows,
        "environment": environment,
        "data": {},
    }
    payload["data"] = typing_data(payload, system, source)
    return OPLSTypingResult(pack(payload))
