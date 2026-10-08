"""Deterministic numerical adapters over validated family-owned assignments."""

from math import fsum, radians

from island.graph import (
    unified_assignment_diagnostics,
    unified_graph_charges,
    unified_parameter_assignment,
    unified_typed_graph,
)
from island.periodic import _require


def _paths(system):
    adj = {i: set() for i in system.topology.sites}
    for bond in system.topology.bonds.values():
        adj[bond.site1].add(bond.site2)
        adj[bond.site2].add(bond.site1)
    pairs = {}
    for start in sorted(adj):
        seen, front = {start}, {start}
        for distance in range(1, 4):
            front = set().union(*(adj[i] for i in front)) - seen
            for end in sorted(front):
                if start < end:
                    pairs[(start, end)] = distance
            seen.update(front)
    return pairs


def binding(system, graph, prepared, config):
    family = {"pcff": "PCFF", "oplsaa": "OPLS-AA", "gaff": "GAFF", "gaff2": "GAFF2"}[
        prepared._family
    ]
    _require(family == config["family"], "Periodic native force-field family mismatch")
    native = prepared.native_result
    data = {
        "family": family,
        "sites": [],
        "bonds": [],
        "angles": [],
        "torsions": [],
        "impropers": [],
        "pairs": [],
    }
    if family == "PCFF":
        model = native.payload
        source = model["source"]
        cr = native.assignment.payload["charge_record"]
        if "typed_graph" in cr:
            typing = cr["typed_graph"]
            perception = typing["validation"]["automatic_perception"]
            origin = cr["origin"]
            charge_method = "provided" if origin == "provided" else "native_increments"
        else:
            typing = cr["automatic_typing"]
            perception = "performed"
            origin = charge_method = "native_increments"
        types = typing["assignments"]
        cr.get("native_charge_record", cr).get(
            "partial_charges", cr.get("charge_data", {}).get("partial_charges", {})
        )
        _require(
            model["model_definition_complete"], "Missing PCFF Class-II interactions"
        )
        # Cartesian Wilson expression is not periodic invariant. Reject before construction.
        _require(
            not any(t["family"] == "wilson_out_of_plane" for t in model["terms"]),
            "Unsupported periodic PCFF interaction family: wilson_out_of_plane",
        )
        data["pcff_terms"] = model["terms"]
        for row in model["nonbonded"]:
            i = row["site"]
            data["sites"].append(
                {
                    "id": i,
                    "mass": system.topology.sites[i].mass,
                    "charge": row["charge"],
                    "sigma": row["rmin"],
                    "epsilon": row["epsilon"],
                    "atom_type": types[i],
                }
            )
        special = model["special_pairs"]
        _require(
            special["lj"][:2] == [0, 0] and special["coulomb"][:2] == [0, 0],
            "Periodic PCFF requires explicit excluded 1-2/1-3 policy",
        )
        scales = [special["lj"][2], special["coulomb"][2]]
        data["pairs"] = [
            {"sites": p["sites"], "lj": p["lj_weight"], "coulomb": p["coulomb_weight"]}
            for p in model["special_pair_inventory"]
        ]
        provenance = (
            "Retained native PCFF typing and charge record; see native assignment"
        )
        evidence = [native.identity, native.assignment.identity]
        native_payload = model
    elif family == "OPLS-AA":
        p = native.payload
        r = p["resolved"]
        source = {
            **prepared.source.identity,
            "sha256": prepared.source.identity["xml_sha256"],
        }
        types = p["charges"]["typing"]["data"]["assignments"]
        perception = "performed"
        origin, charge_method = "native_nonbonded", "native_nonbonded_atom_charge_v1"
        scales = [r["policy"]["lj_14"], r["policy"]["coulomb_14"]]
        _require(not r["impropers"], "Unsupported periodic OPLS improper family")
        for i, row in sorted(r["sites"].items()):
            data["sites"].append(
                {
                    "id": i,
                    "mass": row["mass_dalton"],
                    "charge": row["charge_e"],
                    "sigma": row["sigma_angstrom"],
                    "epsilon": row["epsilon_kj_mol"],
                    "atom_type": types[i],
                }
            )
        for row in r["bonds"]:
            p0 = row["converted"]
            data["bonds"].append(
                {
                    "sites": row["sites"],
                    "length": p0["length_angstrom"],
                    "k": p0["k_kj_mol_angstrom2"],
                }
            )
        for row in r["angles"]:
            p0 = row["converted"]
            data["angles"].append(
                {
                    "sites": row["sites"],
                    "theta": radians(p0["angle_degrees"]),
                    "k": p0["k_kj_mol_radian2"],
                }
            )
        for row in r["rb"]:
            data["torsions"].append(
                {"sites": row["sites"], "rb": row["converted"]["coefficients_kj_mol"]}
            )
        data["pairs"] = [
            {"sites": p["sites"], "lj": p["lj_scale"], "coulomb": p["coulomb_scale"]}
            for p in r["pairs"]
        ]
        provenance = (
            "Retained pinned Foyer typing, native charges and resolved XML rows"
        )
        evidence = [native.identity]
        native_payload = p
    else:
        r = native.imported_result
        source = {"family": family, "sha256": r.source_sha256, "source": r.source}
        types = dict(r.atom_types)
        perception = "performed"
        charge_method = native.record["charge_method"]
        origin = "provided" if charge_method == "provided" else "am1bcc"
        policy = r.nonbonded_policy
        _require(
            policy.mixing_rule == "lorentz_berthelot",
            "Unsupported periodic Amber mixing rule",
        )
        scales = [policy.lj_scale_14, policy.coulomb_scale_14]
        for i in sorted(types):
            p0 = r.site_assignments[i].parameter
            data["sites"].append(
                {
                    "id": i,
                    "mass": system.topology.sites[i].mass,
                    "charge": r.charge_result.assignments[i].charge,
                    "sigma": 10 * p0.sigma,
                    "epsilon": p0.epsilon,
                    "atom_type": types[i],
                }
            )
        for sites, selection in sorted(r.bond_assignments.items()):
            p0 = selection.parameter
            data["bonds"].append(
                {
                    "sites": list(sites),
                    "length": 10 * p0.equilibrium_length,
                    "k": p0.force_constant / 100,
                }
            )
        for sites, selection in sorted(r.angle_assignments.items()):
            p0 = selection.parameter
            data["angles"].append(
                {
                    "sites": list(sites),
                    "theta": radians(p0.equilibrium_angle),
                    "k": p0.force_constant,
                }
            )
        for key, assignments in (
            ("torsions", r.proper_torsion_assignments),
            ("impropers", r.improper_assignments),
        ):
            for sites, selection in sorted(assignments.items()):
                data[key].append(
                    {
                        "sites": list(sites),
                        "fourier": [
                            [t.force_constant, t.periodicity, radians(t.phase)]
                            for t in selection.parameter.terms
                        ],
                    }
                )
        scaled = set(r.source_14_pairs)
        data["pairs"] = [
            {
                "sites": list(p),
                "lj": scales[0] if p in scaled else 0,
                "coulomb": scales[1] if p in scaled else 0,
            }
            for p in sorted(r.source_exclusions)
        ]
        provenance = str(r.provenance)
        evidence = [native.record_signature, r.result_signature]
        native_payload = {
            "record": dict(native.record),
            "imported_result_signature": r.result_signature,
        }
    _require(source == config["source"], "Periodic source identity mismatch")
    _require(
        scales == config["one_four_scaling"],
        "Periodic 1-4 scaling contradicts native assignment",
    )
    _require(
        all(type(s["id"]) is int and s["id"] > 0 for s in data["sites"]),
        "Periodic backend requires positive stable integer atom IDs",
    )
    paths = _paths(system)
    expected = {
        p: (scales if distance == 3 else [0, 0]) for p, distance in paths.items()
    }
    actual = {tuple(sorted(p["sites"])): [p["lj"], p["coulomb"]] for p in data["pairs"]}
    _require(
        len(actual) == len(data["pairs"]) and actual == expected,
        "Native exclusions/1-4 pairs differ from supported shortest-path policy",
    )
    data["sites"].sort(key=lambda p: p["id"])
    charges = {p["id"]: p["charge"] for p in data["sites"]}
    typed = unified_typed_graph(
        graph,
        force_field=family,
        atom_types=types,
        source=source,
        typing_profile="retained_native_periodic_binding_v1",
        typing_method="retained_native_assignment",
        provenance=provenance,
        evidence=evidence,
        automatic_perception=perception,
    )
    charge = unified_graph_charges(
        graph,
        typed,
        force_field=family,
        charges=charges,
        charge_method=charge_method,
        charge_origin=origin,
        source=source,
        component_totals={
            str(c[0]): fsum(charges[i] for i in c) for c in graph.payload["components"]
        },
        total_charge=fsum(charges.values()),
        provenance=provenance,
        evidence=evidence,
    )
    diagnostics = unified_assignment_diagnostics(
        typing_complete=True,
        charges_complete=True,
        native_charge_available=origin in {"native_increments", "native_nonbonded"},
        parameter_rows_complete=True,
        executable_potential_complete=True,
        periodic_backend_support=True,
        reasons=[],
    )
    assignment = unified_parameter_assignment(
        graph,
        typed,
        charge,
        force_field=family,
        source=source,
        native_payload=native_payload,
        diagnostics=diagnostics,
    )
    return typed, charge, assignment, data
