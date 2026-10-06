"""Data-only operational inspection, including source lookups after charge failure.

An inspection is not a parameter assignment, model or prepared force field.
Existing native validators retain exclusive authority over those contracts.
"""

from copy import deepcopy

from .automatic import PCFFAutomaticTypingResult, assign_automatic_pcff_charges
from .catalog import inspect_pcff_full_source
from .class2 import assign_pcff_parameters, canonical, interaction_requests
from .fallbacks import LOWER, resolve, validate_policy
from .model import define_pcff_model, special_pair_policy
from .source import FLAGS, boundary, require


@boundary
def inspect_pcff_operational_support(
    system, typing, *, resolution_policy, special_pairs
):
    """Return owned diagnostic data with all attempted source paths and dependencies.

    Unlike parameterization, structural source inspection can continue after a
    failed native charge check. No partial charge vector is passed to assignment
    or model construction. Missing/ambiguous dependencies remain explicit.
    """
    require(
        type(typing) is PCFFAutomaticTypingResult,
        "Expected native automatic/checked typing",
    )
    typing.validate_integrity(system)
    require(
        typing.complete, "Complete compatible typing required for structural inspection"
    )
    validate_policy(resolution_policy)
    require(
        special_pairs
        == special_pair_policy(
            lj=special_pairs["lj"], coulomb=special_pairs["coulomb"]
        ),
        "Contradictory special-pair policy",
    )
    data = typing.payload
    require(
        data["profile"]["name"].startswith("island_pcff_source_graph_v"),
        "Operational source inspection requires expanded source typing",
    )
    source = typing.source
    catalog = inspect_pcff_full_source(source)
    indexed = {r["id"]: r for r in catalog["records"]}
    inventory = source.inventory
    requests, inventories = interaction_requests(data["graph"], expanded=True)
    queries = []
    base = {}
    cache = {}
    for family, sites, deps, reason in requests:
        labels = tuple(data["assignments"][i] for i in sites)
        if reason:
            query = {"status": "not_applicable", "reason": reason}
        else:
            key = family, labels
            if key not in cache:
                cache[key] = resolve(
                    family,
                    labels,
                    indexed,
                    inventory,
                    policy=resolution_policy,
                    trace=True,
                )
            query = deepcopy(cache[key])
        dependencies = []
        for f, ids in deps:
            b = base[f, canonical(ids)]
            dependencies.append(
                {
                    "request_id": b["request_id"],
                    "sites": list(ids),
                    "status": b["resolution"]["status"],
                    "equilibrium_value": b["resolution"].get(
                        "normalized_values", [None]
                    )[0],
                    "source_rows": [
                        r["record_id"] for r in b["resolution"].get("selected", [])
                    ],
                }
            )
        entry = {
            "request_id": family + ":" + ",".join(map(str, sites)),
            "family": family,
            "sites": list(sites),
            "supplied_types": list(labels),
            "resolution": query,
            "equilibrium_dependencies": dependencies,
            "dependencies_complete": all(
                d["status"] == "assigned" for d in dependencies
            ),
            "automatic_route": "source-declared lower-order family " + LOWER[family]
            if family in LOWER
            else "no automatic cross-term or other family inferred",
        }
        queries.append(entry)
        if family in ("quartic_bond", "quartic_angle"):
            base[family, canonical(sites)] = entry
    charge = assign_automatic_pcff_charges(
        system, typing, resolution_policy=resolution_policy
    )
    result = {
        "schema": "island_pcff_operational_inspection_v1",
        "diagnostic_only": True,
        "source": source.identity,
        "typing_identity": typing.identity,
        "graph_identity": data["graph_identity"],
        "profile_identity": data["profile_identity"],
        "resolution_policy": resolution_policy,
        "special_pairs": deepcopy(special_pairs),
        "typing_entries": data["entries"],
        "native_charge_record": charge.payload,
        "charge_identity": charge.identity,
        "charges_complete": charge.complete,
        "inventories": inventories,
        "interaction_queries": queries,
        "assignment_identity": None,
        "model_identity": None,
        "model_complete": False,
        "model_diagnostics": [
            {
                "reason": "native_charge_incomplete; parameterization and model construction not attempted"
            }
        ],
        "sections": [
            {
                "family": s["name"],
                "namespace": s["namespace"],
                "numerical_records": s["numerical_records"],
            }
            for s in catalog["sections"]
        ],
        **FLAGS,
    }
    # Exact raw numerical rows are indexed once, rather than copied per site.
    candidates = {
        r["record_id"] for q in queries for r in q["resolution"].get("candidates", [])
    }
    result["source_rows"] = {i: deepcopy(indexed[i]) for i in sorted(candidates)}
    if charge.complete:
        assignment = assign_pcff_parameters(
            system, typing, charge, resolution_policy=resolution_policy
        )
        model = define_pcff_model(assignment, special_pairs=special_pairs)
        result.update(
            assignment_identity=assignment.identity,
            model_identity=model.identity,
            model_complete=model.payload["model_definition_complete"],
            model_diagnostics=model.payload["diagnostics"],
        )
    return deepcopy(result)
