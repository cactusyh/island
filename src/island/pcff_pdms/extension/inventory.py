"""Unique extension terms versus stable graph instances; no parameter lookup."""

from collections import Counter
from itertools import combinations

from island.pcff_pdms import domain_diagnostics
from island.workflows.bundle import system_from

from . import FAMILY, FIXTURES, FLAGS, PROFILE, digest, require, signed

# Equations describe the proposed DERIVED contract, not recovered FRC authority.
# Native source phases/AAT conventions remain unresolved in historical evidence.
FORMS = {
    "quartic_bond": ("sum(n=2..4,Kn*(r-r0)^n)", ["r0", "K2", "K3", "K4"]),
    "quartic_angle": ("sum(n=2..4,Kn*(theta-theta0)^n)", ["theta0", "K2", "K3", "K4"]),
    "torsion_3": (
        "sum(n=1..3,Vn*(1-cos(n*phi-phase_n)))",
        ["V1", "phase1", "V2", "phase2", "V3", "phase3"],
    ),
    "bond-bond": ("Kbb*dr12*dr23", ["Kbb"]),
    "bond-angle": ("(Kleft*dr12+Kright*dr23)*dtheta123", ["Kleft", "Kright"]),
    "angle-angle": ("Kaa*dtheta123*dtheta324", ["Kaa"]),
    "bond-bond_1_3": ("Kbb13*dr12*dr34", ["Kbb13"]),
    "end_bond-torsion_3": (
        "dr12*sum(Ln*cos(n*phi))+dr34*sum(Rn*cos(n*phi))",
        ["L1", "L2", "L3", "R1", "R2", "R3"],
    ),
    "middle_bond-torsion_3": ("dr23*sum(Fn*cos(n*phi))", ["F1", "F2", "F3"]),
    "angle-torsion_3": (
        "dtheta123*sum(Ln*cos(n*phi))+dtheta234*sum(Rn*cos(n*phi))",
        ["L1", "L2", "L3", "R1", "R2", "R3"],
    ),
    "angle-angle-torsion_1": ("Kaat*dtheta123*dtheta234*cos(phi)", ["Kaat"]),
    "partial_charge": ("sum(i<j,k_e*qi*qj/r_ij*explicit_coulomb_weight_ij)", ["q"]),
    "lj_9_6": (
        "epsilon_ij*(2*(rmin_ij/r)^9-3*(rmin_ij/r)^6)*explicit_LJ_weight_ij",
        ["rmin", "epsilon"],
    ),
}
TORSIONAL = tuple(f for f in FORMS if "torsion" in f or f == "bond-bond_1_3")
POLICIES = (
    "atom_typing",
    "charge_origin",
    "mixing_rule",
    "one_two_rules",
    "one_three_rules",
    "one_four_rules",
    "exclusions_and_scaling",
    "long_range_electrostatics",
    "cutoff",
    "tail",
    "torsion_torsion_applicability",
)


def scope_check(payload):
    require(
        payload.get("box") is None,
        "outside_derived_PDMS_scope",
        diagnostics=[{"reason": "nonperiodic_finite_molecules_required"}],
    )
    system = system_from(payload)
    reasons = domain_diagnostics(system)
    si = sum(a.element == "Si" for a in system.topology.sites.values())
    # Exact graph sizes in this phase, not extrapolation along a homologous series.
    if si not in (2, 3, 5, 7):
        reasons.append({"reason": "only_HMDS_and_declared_capped_DP2_DP4_DP6"})
    if payload.get("metadata", {}).get("multiplicity", 1) != 1:
        reasons.append({"reason": "closed_shell_singlet_required"})
    require(not reasons, "outside_derived_PDMS_scope", diagnostics=reasons)
    return system


def chemical_identity(payload):
    system = scope_check(payload)
    return digest(
        {
            "sites": [
                {
                    "id": i,
                    "element": a.element,
                    "formal_charge": a.formal_charge,
                    "radical_electrons": a.metadata.get("radical_electrons", 0),
                }
                for i, a in sorted(system.topology.sites.items())
            ],
            "bonds": [
                {"sites": list(b.key), "order": b.order, "aromatic": b.aromatic}
                for _, b in sorted(system.topology.bonds.items())
            ],
        }
    )


def environment_types(system):
    """New namespace; no FRC type/parameter authority is implied."""
    topo = system.topology
    types = {}
    for i, a in topo.sites.items():
        if a.element == "Si":
            types[i] = (
                "ext_Si_T"
                if sum(topo.sites[n].element == "O" for n in topo.neighbors(i)) == 1
                else "ext_Si_D"
            )
        elif a.element == "O":
            types[i] = "ext_O_bridge"
    for i, a in topo.sites.items():
        if a.element == "C":
            parent = next(n for n in topo.neighbors(i) if topo.sites[n].element == "Si")
            types[i] = "ext_C_" + types[parent].removeprefix("ext_Si_")
    for i, a in topo.sites.items():
        if a.element == "H":
            types[i] = "ext_H_" + types[next(iter(topo.neighbors(i)))].removeprefix(
                "ext_C_"
            )
    return types


def canonical_term(family, types):
    alternatives = [tuple(types)]
    if family == "angle-angle":
        alternatives.append((types[3], types[1], types[2], types[0]))
    elif len(types) > 1:
        alternatives.append(tuple(reversed(types)))
    return min(alternatives)


def units(family, field):
    if field == "r0" or field == "rmin":
        return "angstrom"
    if field == "theta0" or field.startswith("phase"):
        return "radian"
    if field == "q":
        return "e"
    if family == "quartic_bond":
        return "kcal/mol/angstrom^" + field[1:]
    if family == "quartic_angle":
        return "kcal/mol/radian^" + field[1:]
    return {
        "bond-bond": "kcal/mol/angstrom^2",
        "bond-bond_1_3": "kcal/mol/angstrom^2",
        "bond-angle": "kcal/mol/angstrom/radian",
        "angle-angle": "kcal/mol/radian^2",
        "end_bond-torsion_3": "kcal/mol/angstrom",
        "middle_bond-torsion_3": "kcal/mol/angstrom",
        "angle-torsion_3": "kcal/mol/radian",
        "angle-angle-torsion_1": "kcal/mol/radian^2",
    }.get(family, "kcal/mol")


def inventory(graphs):
    require(tuple(graphs) == FIXTURES, "exact_declared_graph_inventory_required")
    terms, instances, graphs_out = {}, [], []
    for name, payload in graphs.items():
        system = scope_check(payload)
        topo = system.topology
        topo.rebuild_derived_interactions()  # authoritative core graph enumeration
        types = environment_types(system)
        graph_id = chemical_identity(payload)
        graphs_out.append(
            {
                "fixture": name,
                "graph_identity": graph_id,
                "atom_order": sorted(types),
                "sites": len(types),
                "types": [[i, types[i]] for i in sorted(types)],
            }
        )

        def add(
            family,
            sites,
            dependencies=(),
            *,
            types=types,
            name=name,
            graph_id=graph_id,
            topo=topo,
        ):
            labels = [types[i] for i in sites]
            ordered = canonical_term(family, labels)
            key = family + ":" + ",".join(ordered)
            equation, fields = FORMS[family]
            aliases = {}
            # Reversal-equivalent ends are ONE variable, not duplicate aliases.
            if ordered == tuple(reversed(ordered)) and family in (
                "bond-angle",
                "end_bond-torsion_3",
                "angle-torsion_3",
            ):
                half = len(fields) // 2
                aliases = dict(zip(fields[half:], fields[:half], strict=True))
            if key not in terms:
                terms[key] = signed(
                    {
                        "id": key,
                        "family": family,
                        "ordered_types": list(ordered),
                        "equation": equation,
                        "force_equation": "F_i=-grad_x_i U; analytic feature derivative required before physical fitting",
                        "fields": [f for f in fields if f not in aliases],
                        "field_aliases": aliases,
                        "units": {f: units(family, f) for f in fields},
                        "origin": "unqualified_extension_target",
                        "value": None,
                        "status": "unresolved",
                        "native_PCFF_row": False,
                        "ordering": "AA swaps outer 1/4 only; proper reversal swaps left/right blocks",
                        "prefactor": "no one-half for Class-II polynomials",
                        "convention_status": "proposed_extension_equation; not scientifically_qualified",
                    }
                )
            instance_id = name + ":" + family + ":" + ",".join(map(str, sites))
            instances.append(
                signed(
                    {
                        "id": instance_id,
                        "parameter_id": key,
                        "fixture": name,
                        "graph_identity": graph_id,
                        "sites": list(sites),
                        "source_order_types": labels,
                        "orientation_reversed": tuple(labels) != ordered,
                        "equilibrium_dependencies": [list(d) for d in dependencies],
                        "equilibrium_parameter_dependencies": [
                            ("quartic_bond" if len(d) == 2 else "quartic_angle")
                            + ":"
                            + ",".join(
                                canonical_term(
                                    "quartic_bond" if len(d) == 2 else "quartic_angle",
                                    [types[i] for i in d],
                                )
                            )
                            + ("/r0" if len(d) == 2 else "/theta0")
                            for d in dependencies
                        ],
                        "dependency_status": "unresolved_reference_geometries"
                        if dependencies
                        else "not_applicable",
                        "local_environment": [
                            {
                                "site": i,
                                "type": types[i],
                                "neighbors": sorted(topo.neighbors(i)),
                            }
                            for i in sites
                        ],
                    }
                )
            )

        for b in sorted(topo.bonds):
            add("quartic_bond", b)
        for a in sorted(topo.angles):
            add("quartic_angle", a)
            add("bond-bond", a, [a[:2], a[1:]])
            add("bond-angle", a, [a[:2], a[1:], a])
        for d in sorted(topo.dihedrals):
            for family in TORSIONAL:
                dep = (
                    [d[:3], d[1:]]
                    if family in ("angle-torsion_3", "angle-angle-torsion_1")
                    else [d[:2], d[2:]]
                    if family in ("end_bond-torsion_3", "bond-bond_1_3")
                    else [d[1:3]]
                    if family == "middle_bond-torsion_3"
                    else []
                )
                add(family, d, dep)
        for center in sorted(types):
            for triple in combinations(sorted(topo.neighbors(center)), 3):
                for shared in triple:
                    outer = sorted(set(triple) - {shared})
                    aa = (outer[0], center, shared, outer[1])
                    add("angle-angle", aa, [aa[:3], (aa[2], aa[1], aa[3])])
        for i in sorted(types):
            add("partial_charge", (i,))
            add("lj_9_6", (i,))
    validate_inventory(list(terms.values()), instances)
    return signed(
        {
            "schema": "island.pdms-extension-inventory.v1",
            "family": FAMILY,
            "profile": PROFILE,
            "graphs": graphs_out,
            "terms": [terms[k] for k in sorted(terms)],
            "instances": sorted(instances, key=lambda x: x["id"]),
            "unique_terms": len(terms),
            "site_instances": len(instances),
            "family_counts": dict(Counter(t["family"] for t in terms.values())),
            "required_policies": [
                {"id": p, "status": "unresolved", "value": None} for p in POLICIES
            ],
            "structural_nonapplicability": [
                "Si-C-C absent: all carbon is methyl bonded to Si",
                "Wilson improper: no three-connected center",
            ],
            "torsion_torsion": "applicability/convention unresolved; no zero inferred",
            "native_PCFF_unresolved_instances": 378,
            **FLAGS,
        }
    )


def validate_inventory(terms, instances):
    keys = [t["id"] for t in terms]
    require(len(keys) == len(set(keys)), "duplicate_parameter_alias")
    require(
        len(instances) == len({i["id"] for i in instances}), "duplicate_site_instance"
    )
    require(
        all(i["parameter_id"] in keys for i in instances), "missing_Class_II_parameter"
    )
    references = {t["id"] + "/" + f for t in terms for f in t["fields"]}
    require(
        all(
            d in references
            for i in instances
            for d in i["equilibrium_parameter_dependencies"]
        ),
        "missing_Class_II_equilibrium_dependency",
    )
    for t in terms:
        require(
            t["origin"] == "unqualified_extension_target" and t["value"] is None,
            "foreign_family_or_analog_substitution",
        )
