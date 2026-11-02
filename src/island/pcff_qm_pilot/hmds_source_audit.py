"""Read-only source closure for the retained HMDS research Hamiltonian.

Missing source coefficients remain absent. This module provides no assignment,
QM, dynamics, parameter-fitting, or source-row editing interface.
"""

from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations

import numpy as np

from island.forcefields.pcff.class2 import FAMILIES, inspect_pcff_class2, resolve
from island.forcefields.pcff.source import PCFFSource, records
from island.pcff_completion import require, signed
from scripts.pcff_j5_reference import equivalents, match, read_source

from .hamiltonian import Jet, bend, cos_phi, harmonic, length
from .hmds_hamiltonian import FORMS, graph_inventory, validate
from .hmds_local_baseline import fixed_families
from .molecular_control import cartesian

FRC = "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c"
FORCE = "F_a = -sum_q (partial E/partial q) (partial q/partial x_a)"
EQUATIONS = {
    **FORMS,
    "nonbond(9-6)": "eps[2(rmin/r)^9-3(rmin/r)^6]",
}
CONVENTIONS = {
    "source_energy": "kcal/mol",
    "internal_energy": "kJ/mol",
    "execution_energy": "kcal/mol",
    "kcal_to_kJ": 4.184,
    "equilibrium_angle_source": "degree",
    "displacement_angle": "radian",
    "length": "angstrom",
    "polynomial_prefactor": 1.0,
    "torsion_source": "sum Vn[1+cos(n*phi-phase_n)]",
    "LAMMPS_phase_mapping": "(source phase + 180) modulo 360 degree",
    "AA_roles": "I-J-K and K-J-L; J center K shared; swap I/L only",
    "phase_geometry_authority": "research convention; native FRC Phi geometry not independently supplied",
}


def source_graph(model, raw):
    """Cross-check independent Decimal ASCII oracle and actual production parser."""
    validate(model)
    require(sha256(raw).hexdigest() == FRC, "J69_source_SHA_mutation")
    oracle = read_source(raw)
    # Historical policy is exact then ordinary equivalence, without wildcard or
    # lower-order substitution. Do not import a more permissive fallback policy.
    exact = defaultdict(list)
    for k, rows in oracle.items():
        exact[k] = [
            r for r in rows if not any(t.startswith("*") for t in r.get("types", []))
        ]
    source = PCFFSource(raw, FRC)
    catalog = inspect_pcff_class2(source)["records"]
    eq = records(source.inventory, "equivalence")
    graph = graph_inventory(model)
    require(graph["all_required_graph_roles_present"], "J69_graph_interaction_mapping")
    nodes, missing, selected = [], [], {}
    for row in model["interaction_inventory"]:
        family, sites = row["family"], row["sites"]
        if family in ("atom_type", "mass", "charge"):
            nodes.append(
                {
                    "id": row["id"],
                    "family": family,
                    "stable_sites": sites,
                    "site_types": [model["typing"][str(i)] for i in sites],
                    "kind": "source_metadata",
                    "assignment": row["assignment"],
                    "dependencies": [
                        {"interaction": d, "role": "source_metadata"}
                        for d in row["dependencies"]
                    ],
                }
            )
            continue
        labels = [model["typing"][str(i)] for i in sites]
        trial = match(exact, family, "cff91", labels)
        searches = [
            {"types": labels, "status": trial["status"], "equivalence_rows": []}
        ]
        if trial["status"] == "missing":
            mapped, evidence = equivalents(
                oracle, labels, "equivalence", [FAMILIES[family][1]] * len(labels)
            )
            if mapped:
                trial = match(exact, family, "cff91", mapped)
                searches.append(
                    {
                        "types": mapped,
                        "status": trial["status"],
                        "equivalence_rows": evidence,
                    }
                )
        parsed = resolve(family, labels, catalog, eq)
        a = row["assignment"]
        require(
            (a is None)
            == (trial["status"] == "missing")
            == (parsed["status"] == "missing"),
            "J69_missing_cross_term_mutation",
            interaction=row["id"],
        )
        if a is not None:
            require(
                trial["status"] == parsed["status"] == "assigned",
                "J69_source_ambiguity",
            )
            require(
                set(a["row_ids"]) == set(trial["selected_ids"]),
                "J69_source_row_mapping",
                interaction=row["id"],
            )
            require(
                np.max(abs(np.array(a["values"]) - trial["values"])) < 1e-11
                and np.array_equal(a["values"], parsed["normalized_values"]),
                "J69_source_values_units_order",
                interaction=row["id"],
            )
            selected[row["id"]] = trial["values"]
        else:
            missing.append(row["id"])
        dependencies = []
        for dep in row["dependencies"]:
            parent = next(r for r in model["interaction_inventory"] if r["id"] == dep)
            require(parent["assignment"] is not None, "J69_missing_equilibrium")
            dependencies.append(
                {
                    "interaction": dep,
                    "source_rows": parent["assignment"]["row_ids"],
                    "equilibrium": parent["assignment"]["values"][0],
                }
            )
        values = None if a is None else list(a["values"])
        if values is not None:
            # Match explicit field dimensions, never divide equilibrium by4.184.
            dims = FAMILIES[family][2]
            execution = [
                v / 4.184 if u.startswith("E") else v
                for u, v in zip(dims, values, strict=True)
            ]
        else:
            execution = None
        lammps = None if execution is None else list(execution)
        if lammps is not None and family == "quartic_angle":
            lammps[0] *= 180 / np.pi
        if lammps is not None and family == "torsion_3":
            for j in (1, 3, 5):
                lammps[j] = (lammps[j] * 180 / np.pi + 180) % 360
        if lammps is not None and family == "nonbond(9-6)":
            lammps = [lammps[1], lammps[0]]  # pair_coeff epsilon rmin
        if lammps is not None and family in ("bond-bond", "bond-angle", "angle-angle"):
            needed = [
                d
                for d in dependencies
                if family != "bond-angle"
                or d["interaction"].startswith("quartic_bond:")
            ]
            lammps += [
                d["equilibrium"]
                * (180 / np.pi if d["interaction"].startswith("quartic_angle:") else 1)
                for d in needed
            ]
        nodes.append(
            {
                "id": row["id"],
                "family": family,
                "stable_sites": sites,
                "site_types": labels,
                "source_rows": [] if a is None else a["row_ids"],
                "source_status": "source_unavailable"
                if a is None
                else "native_row_available",
                "kind": "explicit_research_unknown"
                if a is None
                else "converted_source_term",
                "searches": searches,
                "internal_values": values,
                "execution_values": execution,
                "equation": EQUATIONS[family],
                "force_equation": FORCE,
                "dependencies": dependencies,
                "source_units": FAMILIES[family][2],
                "execution_units": "kcal/mol,angstrom,radian; coefficients carry equation dimensions",
                "multiplicity": "one per enumerated interaction; AA shared-arm slots explicit; n=1,2,3 torsion harmonics; no class-I wildcard divisor",
                "reference_order": "J center K shared for AA; path order and reversal exchange end blocks otherwise",
                "parser_transform": "source kcal coefficients *4.184; degree equilibria/phases *pi/180; symmetric half repeat and oriented end exchange where applicable",
                "LAMMPS_coefficients": lammps,
                "LAMMPS_origin": "generated research representation from audited rows; NOT a successful full msi2lmp export",
                "LAMMPS_cross_references": [
                    {
                        **d,
                        "LAMMPS_reference": d["equilibrium"]
                        * (
                            180 / np.pi
                            if d["interaction"].startswith("quartic_angle:")
                            else 1
                        ),
                    }
                    for d in dependencies
                ],
                "evaluator": "residual.fixed_rows + hmds_hamiltonian._kernels; unknown columns require explicit inputs",
            }
        )
    return signed(
        {
            "schema": "island.j69-source-dependencies.v1",
            "graph_identity": model["graph_identity"],
            "atom_order": model["atom_order"],
            "conventions": CONVENTIONS,
            "nodes": nodes,
            "missing_interactions": missing,
            "family_counts": dict(Counter(n["family"] for n in nodes)),
            "source_complete": not missing,
            "source_present_not_evaluated": [
                {
                    "section": "torsion-torsion_1",
                    "status": "empty section; no numeric zero authority; outside declared research sum",
                },
                {
                    "section": "angle-angle-torsion_1",
                    "status": "literal linear-Phi form not evaluated; no applicable numeric row; research cos(Phi) unknown explicitly distinct",
                },
            ],
            "Wilson": "degree-three converter form; no degree-three centers; degree-four AA retained",
        }
    ), selected


def validate_graph(record, expected):
    require(record == expected, "J69_source_dependency_semantic_mutation")


def family_label(family, sites):
    if family == "quartic_bond" and set(sites) in ({2, 5}, {5, 6}):
        return "quartic_bond:SiO"
    if family == "quartic_angle" and min(tuple(sites), tuple(sites[::-1])) == (2, 5, 6):
        return "quartic_angle:SiOSi"
    return family


def source_sum(model, selected, xyz, jets=False, literal_converter=False):
    """Independent raw-row scalar sum or differentiated scalar potential.

    No production coefficient, equilibrium, family evaluator or pair helper is
    read here. The scalar dihedral uses projections; AD uses existing Jet geometry.
    """
    x = cartesian(xyz) if jets else np.asarray(xyz)
    zero = lambda: Jet(0, np.zeros(81)) if jets else 0.0
    out, interactions, cache = defaultdict(zero), {}, {}

    def coord(kind, s):
        key = kind, tuple(s)
        if key not in cache:
            if jets:
                v = (length if kind == "r" else bend if kind == "t" else cos_phi)(x, s)
            else:
                p = np.array([x[i - 1] for i in s])
                if kind == "r":
                    v = np.linalg.norm(p[0] - p[1])
                elif kind == "t":
                    u, w = p[0] - p[1], p[2] - p[1]
                    v = np.arccos(
                        np.clip(u @ w / (np.linalg.norm(u) * np.linalg.norm(w)), -1, 1)
                    )
                else:
                    axis = p[2] - p[1]
                    axis /= np.linalg.norm(axis)
                    u, w = p[0] - p[1], p[3] - p[2]
                    u, w = u - (u @ axis) * axis, w - (w @ axis) * axis
                    v = (u @ w) / (np.linalg.norm(u) * np.linalg.norm(w))
            cache[key] = v
        return cache[key]

    def delta(f, s):
        key = f + ":" + ",".join(map(str, min(tuple(s), tuple(s[::-1]))))
        require(key in selected, "J69_missing_source_equilibrium")
        return coord("r" if f == "quartic_bond" else "t", s) - selected[key][0]

    for row in model["interaction_inventory"]:
        if row["family"] in ("atom_type", "mass", "charge"):
            continue
        if row["id"] not in selected:
            continue
        f, s, v = row["family"], row["sites"], selected[row["id"]]
        if f == "nonbond(9-6)":
            continue
        if f in ("quartic_bond", "quartic_angle"):
            q = delta(f, s)
            value = sum(v[k] / 4.184 * q ** (k + 1) for k in (1, 2, 3))
        elif f in ("bond-bond", "bond-angle"):
            i, j, k = s
            a, b = delta("quartic_bond", [i, j]), delta("quartic_bond", [j, k])
            value = (
                v[0] / 4.184 * a * b
                if f == "bond-bond"
                else (v[0] * a + v[1] * b) / 4.184 * delta("quartic_angle", s)
            )
        elif f == "angle-angle":
            i, j, k, l = s
            value = (
                v[0]
                / 4.184
                * delta("quartic_angle", [i, j, k])
                * delta("quartic_angle", [k, j, l])
            )
        elif f == "torsion_3":
            require(all(v[j] == 0 for j in (1, 3, 5)), "J69_nonzero_phase_unsupported")
            c = coord("c", s)
            sign = -1 if literal_converter else 1
            value = sum(
                v[2 * (n - 1)]
                / 4.184
                * (
                    1
                    + sign
                    * (
                        harmonic(c, n)
                        if jets
                        else np.cos(n * np.arccos(np.clip(c, -1, 1)))
                    )
                )
                for n in (1, 2, 3)
            )
        else:
            raise ValueError("source-present family not evaluated: " + f)
        out[family_label(f, s)] += value
        interactions[row["id"]] = value
    adj = {i: set() for i in model["atom_order"]}
    for a, b, _ in model["bonds"]:
        adj[a].add(b)
        adj[b].add(a)
    for i, j in combinations(model["atom_order"], 2):
        queue, seen = deque([(i, 0)]), set()
        while queue:
            k, d = queue.popleft()
            if k == j:
                break
            if k not in seen:
                seen.add(k)
                queue.extend((v, d + 1) for v in adj[k] - seen)
        if d < 3:
            continue
        ri, ei = selected[f"nonbond(9-6):{i}"]
        rj, ej = selected[f"nonbond(9-6):{j}"]
        rm = ((ri**6 + rj**6) / 2) ** (1 / 6)
        ep = 2 * np.sqrt(ei * ej) / 4.184 * ri**3 * rj**3 / (ri**6 + rj**6)
        r = coord("r", [i, j])
        lj = ep * (2 * (rm / r) ** 9 - 3 * (rm / r) ** 6)
        coul = 332.06371 * model["charges"][str(i)] * model["charges"][str(j)] / r
        out["LJ96"] += lj
        out["Coulomb"] += coul
        interactions[f"LJ96:{i},{j}"], interactions[f"Coulomb:{i},{j}"] = lj, coul
    return out, interactions


def closure(model, selected, rows, step=1e-5):
    reports, phase = [], []
    for row in rows:
        xyz = np.array(row["coordinates"])
        family, terms = source_sum(model, selected, xyz, jets=True)
        scalar, scalar_terms = source_sum(model, selected, xyz)
        target = fixed_families(model, xyz)["families"]
        require(set(family) == set(target), "J69_missing_evaluated_family")
        numerical = {k: np.zeros((27, 3)) for k in family}
        term_errors = {
            k: {
                "scalar_energy_difference": float(abs(e.v - scalar_terms[k])),
                "finite_difference_component_difference": 0.0,
                "force_max_atom": None,
                "force_max_component": None,
            }
            for k, e in terms.items()
        }
        for i in range(27):
            for j in range(3):
                xp, xm = xyz.copy(), xyz.copy()
                xp[i, j] += step
                xm[i, j] -= step
                plus, pt = source_sum(model, selected, xp)
                minus, mt = source_sum(model, selected, xm)
                for k in family:
                    numerical[k][i, j] = -(plus[k] - minus[k]) / (2 * step)
                for k, e in terms.items():
                    difference = float(
                        abs(e.d[3 * i + j] - (pt[k] - mt[k]) / (2 * step))
                    )
                    if (
                        difference
                        > term_errors[k]["finite_difference_component_difference"]
                    ):
                        term_errors[k].update(
                            finite_difference_component_difference=difference,
                            force_max_atom=i + 1,
                            force_max_component=j,
                        )
        families = {}
        for k, e in family.items():
            force = -e.d.reshape(27, 3)
            diff = force - np.array(target[k]["forces"])
            fd = force - numerical[k]
            local = np.linalg.norm(diff, axis=1)
            families[k] = {
                "energy_difference": abs(e.v - target[k]["energy"]),
                "scalar_energy_difference": abs(e.v - scalar[k]),
                "force_component_difference": float(abs(diff).max()),
                "force_vector_difference": float(local.max()),
                "force_max_atom": int(np.argmax(local)) + 1,
                "finite_difference_component_difference": float(abs(fd).max()),
                "finite_difference_vector_difference": float(
                    np.linalg.norm(fd, axis=1).max()
                ),
                "finite_difference_max_atom": int(np.argmax(np.linalg.norm(fd, axis=1)))
                + 1,
            }
        literal, _ = source_sum(model, selected, xyz, jets=True, literal_converter=True)
        f = literal["torsion_3"] - family["torsion_3"]
        phase.append(
            {
                "id": row["id"],
                "literal_converter_minus_vs_source_plus_energy": f.v,
                "force_component_difference": float(abs(f.d).max()),
                "force_vector_difference": float(
                    np.linalg.norm(f.d.reshape(27, 3), axis=1).max()
                ),
                "scope": "same declared Phi; converter unchanged phases versus source printed equation; historical difference, not a new internal evaluator mismatch",
            }
        )
        largest = max(terms, key=lambda k: np.linalg.norm(terms[k].d))
        reports.append(
            {
                "id": row["id"],
                "families": families,
                "interactions": term_errors,
                "largest_source_force_interaction": largest,
            }
        )
    require(
        max(v["energy_difference"] for r in reports for v in r["families"].values())
        < 1e-9,
        "J69_energy_closure",
    )
    require(
        max(
            v["force_component_difference"]
            for r in reports
            for v in r["families"].values()
        )
        < 1e-9,
        "J69_force_closure",
    )
    require(
        max(
            v["finite_difference_component_difference"]
            for r in reports
            for v in r["families"].values()
        )
        < 2e-5,
        "J69_force_derivatives",
    )
    return signed(
        {
            "schema": "island.j69-source-closure.v1",
            "energy_unit": "kcal/mol",
            "force_unit": "kcal/mol/angstrom",
            "step_A": step,
            "configurations": reports,
            "phase_controls": phase,
            "complete_source_claim": False,
            "source_available_subset_closed": True,
            "engine_executions": 0,
        }
    )


def projection(vector, residual):
    v, r = np.asarray(vector), np.asarray(residual)
    n, rn = np.linalg.norm(v), np.linalg.norm(r)
    return {
        "vector": v.tolist(),
        "norm": float(n),
        "residual_dot": float(v @ r),
        "along_residual": float(v @ r / rn) if rn else 0.0,
        "cosine": float(v @ r / (n * rn)) if n * rn else None,
    }
