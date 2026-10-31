"""Complete evaluable, source-bound HMDS research Class-II Hamiltonian.

Functional accounting is complete. Missing numerical coefficients remain explicit
variables; no default physical values, native-PCFF or parameter-complete claim.
"""

import json
from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed

from .hamiltonian import COULOMB, Jet, bend, cos_phi, harmonic, length
from .hmds_paths import GRAPH, adjacency
from .molecular_control import cartesian
from .residual import fixed_rows

EXPECTED_MODEL = "d9ae6279735fb219c7cdaf447971dde5995980292f2c019b71710213d5b3851c"
IDENTITY = "island_research_class2_hmds_j63_v1"
SOURCE = "01b9535013896a1b53c6764381565b521d0a472c33b41b2b0e60e7e70b7980fb"
ROWS = "a5b6f034ca7ce7c2513dcc3a6331e8224cd1af7e4bba141a040fa48e61f3d92d"
TYPES = "474c02becdb9025e5d992f04387d6b0221c0f992c588162bfe5870aaf0bc6603"
CHARGES = "8a7156b45bf67e2e5c768e767b189d12410dc2e6d46ca6f525b4625248fc343b"
BONDS = (
    (1, 2),
    (2, 3),
    (2, 4),
    (2, 5),
    (5, 6),
    (6, 7),
    (6, 8),
    (6, 9),
    (1, 10),
    (1, 11),
    (1, 12),
    (3, 13),
    (3, 14),
    (3, 15),
    (4, 16),
    (4, 17),
    (4, 18),
    (7, 19),
    (7, 20),
    (7, 21),
    (8, 22),
    (8, 23),
    (8, 24),
    (9, 25),
    (9, 26),
    (9, 27),
)
UNITS = {
    "angle-angle": "kcal/mol/radian^2",
    "bond-bond_1_3": "kcal/mol/angstrom^2",
    "end_bond-torsion_3": "kcal/mol/angstrom",
    "middle_bond-torsion_3": "kcal/mol/angstrom",
    "angle-torsion_3": "kcal/mol/radian",
    "angle-angle-torsion_1": "kcal/mol/radian^2",
}
FORMS = {
    "quartic_bond": "K2 dr^2 + K3 dr^3 + K4 dr^4",
    "quartic_angle": "K2 dt^2 + K3 dt^3 + K4 dt^4",
    "bond-bond": "K dr_IJ dr_JK",
    "bond-angle": "(Kleft dr_IJ + Kright dr_JK) dt_IJK",
    "angle-angle": "K dt_IJK dt_KJL; J center K shared; I/L exchange only",
    "torsion_3": "sum Vn[1+cos(n phi)] from zero source phases; equivalent engine phase pi; signed amplitudes retained",
    "bond-bond_1_3": "K dr_IJ dr_KL",
    "end_bond-torsion_3": "sum_n [Fn_left dr_IJ + Fn_right dr_KL] cos(n phi)",
    "middle_bond-torsion_3": "sum_n Fn dr_JK cos(n phi)",
    "angle-torsion_3": "sum_n [Fn_left dt_IJK + Fn_right dt_JKL] cos(n phi)",
    "angle-angle-torsion_1": "K dt_IJK dt_JKL cos(phi); selected research form, not literal FRC linear-Phi expression",
}


def load_model(repo):
    repo = Path(repo)
    p = repo / "docs/evidence/phase_4j54/audit/hmds/report.json"
    require(pin(p) == SOURCE, "HMDS_source_archive_mutation")
    r = json.loads(p.read_text())
    model = signed(
        {
            "schema": "island.j63-hmds-Hamiltonian.v1",
            "internal_identity": IDENTITY,
            "graph_identity": GRAPH,
            "atom_order": list(range(1, 28)),
            "bonds": [[i, j, 1.0] for i, j in BONDS],
            "typing": r["typing"],
            "charges": r["charge_vector"],
            "interaction_inventory": r["interaction_inventory"],
            "source_hashes": r["source_hashes"],
            "source_record_sha256": SOURCE,
            "source_record_identity": r["identity"],
            "source_typing_identity": r["typing_identity"],
            "charge_identity": r["charge_identity"],
            "charge_origin": r["charge_origin"],
            "improper_inventory": r["improper_structural_inventory"],
            "units": {
                "energy": "kcal/mol",
                "length": "angstrom",
                "angle": "radian",
                "charge": "e",
            },
            "conventions": {
                "torsion": "source plus cos retained, source phases zero; engine minus cos phase pi; signed coefficients preserved",
                "AAT": "LAMMPS class2 cosine research variant, explicitly distinct from printed FRC linear-Phi",
                "BB13": "all declared dihedral paths retained as unknowns; converter cp-only lookup is not physical zero authority",
                "Wilson": "degree-three converter form, no such center; all degree-four AA contributions retained",
                "torsion_torsion": "outside selected class2 sum, not an assigned numerical zero",
                "LJ": "9-6 sixth-power mixing; r_min, not sigma",
                "pairs": "exclude distances1,2; full weights at distance>=3 for LJ and Coulomb; no arbitrary pair weights",
                "electrostatics": "isolated nonperiodic dielectric1, Coulomb332.06371 kcal A/mol/e2, all nonexcluded pairs within100A, no shifting/tail/long-range summation",
            },
            "J60_equation_authority_sha256": pin(
                repo / "docs/evidence/phase_4j60/hamiltonian-definition.json"
            ),
            "origin": "source-bound fixed rows plus explicitly unfilled research-extension variables",
            "original_author_artifact_status": "not_yet_supplied",
            "parameter_completeness": "incomplete",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    validate(model)
    return model


def validate(model):
    require(model["identity"] == EXPECTED_MODEL, "HMDS_external_binding_mutation")
    require(
        model["identity"]
        == digest({k: v for k, v in model.items() if k != "identity"}),
        "HMDS_model_identity_mutation",
    )
    require(
        model["internal_identity"] == IDENTITY
        and model["graph_identity"] == GRAPH
        and model["source_record_sha256"] == SOURCE,
        "unsupported_HMDS_model",
    )
    require(
        digest(model["interaction_inventory"]) == ROWS
        and digest(model["typing"]) == TYPES
        and digest(model["charges"]) == CHARGES,
        "HMDS_source_type_charge_mutation",
    )
    require(
        model["bonds"] == [[i, j, 1.0] for i, j in BONDS]
        and model["units"]
        == {
            "energy": "kcal/mol",
            "length": "angstrom",
            "angle": "radian",
            "charge": "e",
        },
        "HMDS_graph_or_unit_mutation",
    )
    adjacency(model)


def oriented(row, types):
    s = list(row["sites"])
    labels = tuple(types[str(i)] for i in s)
    if row["family"] == "angle-angle":
        reverse = (labels[3], labels[1], labels[2], labels[0])
        if reverse < labels:
            s = [s[3], s[1], s[2], s[0]]
            labels = reverse
    elif labels[::-1] < labels:
        s = s[::-1]
        labels = labels[::-1]
    return s, labels


def _kernels(model, xyz):
    x = cartesian(xyz)
    require(len(x) == 27, "HMDS_atom_mapping")
    rows = model["interaction_inventory"]
    lookup = {r["id"]: r for r in rows}
    types = model["typing"]

    def delta(f, s):
        key = f + ":" + ",".join(map(str, min(tuple(s), tuple(s[::-1]))))
        row = lookup[key]
        require(
            row["assignment"] is not None, "HMDS_equilibrium_dependency_missing", id=key
        )
        return (length(x, s) if f == "quartic_bond" else bend(x, s)) - row[
            "assignment"
        ]["values"][0]

    dr = lambda s: delta("quartic_bond", s)
    dt = lambda s: delta("quartic_angle", s)
    zero = lambda: Jet(0, np.zeros(81))
    fixed = fixed_rows(rows, xyz)
    energy = fixed["energy"]
    force = np.array(fixed["forces"])
    components = dict(fixed["components"])
    extra = defaultdict(zero)
    columns = defaultdict(zero)
    mapping = []
    for row in rows:
        f, s, a = row["family"], row["sites"], row["assignment"]
        if f == "torsion_3":
            v = a["values"]
            require(
                all(abs(v[i]) < 1e-14 for i in (1, 3, 5)),
                "unsupported_HMDS_source_phase",
            )
            c = cos_phi(x, s)
            extra[f] += sum(
                v[2 * (n - 1)] / 4.184 * (1 + harmonic(c, n)) for n in (1, 2, 3)
            )
        elif f == "angle-angle" and a is not None:
            i, j, k, l = s
            extra[f] += a["values"][0] / 4.184 * dt([i, j, k]) * dt([k, j, l])
        elif a is None:
            require(f in UNITS, "unresolved_HMDS_interaction", family=f)
            s, labels = oriented(row, types)
            i, j, k, l = s
            key = f + ":" + ",".join(labels)
            used = []

            def add(name, value, used=used):
                columns[name] += value
                used.append(name)

            if f == "angle-angle":
                add(key, dt([i, j, k]) * dt([k, j, l]))
            elif f == "bond-bond_1_3":
                add(key, dr([i, j]) * dr([k, l]))
            elif f == "angle-angle-torsion_1":
                add(key, dt([i, j, k]) * dt([j, k, l]) * cos_phi(x, s))
            else:
                ds = (
                    [dr([i, j]), dr([k, l])]
                    if f == "end_bond-torsion_3"
                    else [dr([j, k])]
                    if f == "middle_bond-torsion_3"
                    else [dt([i, j, k]), dt([j, k, l])]
                )
                c = cos_phi(x, s)
                for b, d in enumerate(ds):
                    for n in (1, 2, 3):
                        add(key + f":block{b}:n{n}", d * harmonic(c, n))
            mapping.append(
                {
                    "id": row["id"],
                    "directed_query": f
                    + ":"
                    + ",".join(types[str(i)] for i in row["sites"]),
                    "oriented_sites": s,
                    "canonical_types": labels,
                    "columns": used,
                    "units": UNITS[f],
                }
            )
    pair = zero()
    pair_rows = []
    adj = adjacency(model)
    for i, j in combinations(model["atom_order"], 2):
        todo = deque([(i, 0)])
        seen = set()
        while todo:
            node, d = todo.popleft()
            if node == j:
                break
            if node in seen:
                continue
            seen.add(node)
            todo.extend((k, d + 1) for k in adj[node] if k not in seen)
        weight = 0 if d < 3 else 1
        r = length(x, [i, j])
        require(r.v < 100, "isolated_HMDS_cutoff_exceeded")
        ri, ei = lookup[f"nonbond(9-6):{i}"]["assignment"]["values"]
        rj, ej = lookup[f"nonbond(9-6):{j}"]["assignment"]["values"]
        ei /= 4.184
        ej /= 4.184
        rm = ((ri**6 + rj**6) / 2) ** (1 / 6)
        ep = 2 * np.sqrt(ei * ej) * ri**3 * rj**3 / (ri**6 + rj**6)
        lj = weight * ep * (2 * (rm / r) ** 9 - 3 * (rm / r) ** 6)
        coul = (
            weight * COULOMB * model["charges"][str(i)] * model["charges"][str(j)] / r
        )
        pair += lj + coul
        pair_rows.append(
            {
                "sites": [i, j],
                "graph_distance": d,
                "LJ_weight": weight,
                "Coulomb_weight": weight,
                "LJ_energy": lj.v,
                "Coulomb_energy": coul.v,
            }
        )
    extra["LJ96_Coulomb"] = pair
    for f, e in extra.items():
        components[f] = e.v
        energy += e.v
        force -= e.d.reshape(27, 3)
    names = sorted(columns)
    return {
        "fixed_energy": energy,
        "fixed_forces": force,
        "components": components,
        "columns": names,
        "energy_jacobian": np.array([columns[n].v for n in names]),
        "force_jacobian": np.stack(
            [-columns[n].d.reshape(27, 3) for n in names], axis=-1
        ),
        "query_to_columns": mapping,
        "pairs": pair_rows,
        "functional_accounting_complete": True,
        "physical_coefficients_supplied": False,
    }


def parts(model, xyz):
    validate(model)
    return _kernels(model, xyz)


def evaluate(model, xyz, coefficients=None):
    p = parts(model, xyz)
    require(
        coefficients is not None and set(coefficients) == set(p["columns"]),
        "explicit_all_HMDS_unknown_coefficients_required",
        required=p["columns"],
    )
    v = np.array([coefficients[n] for n in p["columns"]], float)
    require(np.isfinite(v).all(), "nonfinite_HMDS_coefficient")
    return {
        "energy": float(p["fixed_energy"] + p["energy_jacobian"] @ v),
        "forces": p["fixed_forces"] + p["force_jacobian"] @ v,
        "Hamiltonian": IDENTITY,
        "parameter_origin": "explicit_research_inputs_not_source_rows",
    }


def accounting(model, xyz):
    p = parts(model, xyz)
    table = []
    for r in model["interaction_inventory"]:
        table.append(
            {
                **r,
                "classification": "verified_fixed_source_contribution"
                if r["assignment"] is not None
                else "explicit_fitted_unknown",
                "functional_form": FORMS.get(r["family"], "source metadata row"),
                "roles": "J center K shared I/L exchange only"
                if r["family"] == "angle-angle"
                else "path reversal exchanges end blocks; no arbitrary permutations",
                "source_equilibrium_dependencies": r["dependencies"],
                "numerical_source_values": None
                if r["assignment"] is None
                else r["assignment"].get("values"),
                "execution_units": UNITS.get(
                    r["family"],
                    "energy kcal/mol; distance angstrom; angle radian; charge e; masses source metadata",
                ),
                "source_to_execution_conversion": "Energy-bearing source coefficients kJ/mol divided by4.184; source equilibrium lengths/angles retained in angstrom/radian; no charge normalization",
            }
        )
    for r in model["improper_inventory"]:
        table.append(
            {
                **r,
                "classification": "structural_inapplicability",
                "family": "wilson_out_of_plane",
                "authority": "degree-three msi2lmp convention explicitly selected; no degree-three center in this graph",
            }
        )
    adj = adjacency(model)
    paths = set()

    def walk(s):
        if len(s) == 5:
            paths.add(min(tuple(s), tuple(s[::-1])))
            return
        for j in adj[s[-1]]:
            if j not in s:
                walk([*s, j])

    for i in adj:
        walk([i])
    for s in sorted(paths):
        table.append(
            {
                "id": "torsion-torsion_1:" + ",".join(map(str, s)),
                "sites": s,
                "family": "torsion-torsion_1",
                "classification": "model_inapplicability",
                "reason": "outside declared LAMMPS class2 research sum; not a numerical zero or native-PCFF omission claim",
            }
        )
    for r in p["pairs"]:
        i, j = r["sites"]
        lookup = {row["id"]: row for row in model["interaction_inventory"]}
        ri, ei = lookup[f"nonbond(9-6):{i}"]["assignment"]["values"]
        rj, ej = lookup[f"nonbond(9-6):{j}"]["assignment"]["values"]
        table.append(
            {
                **r,
                "family": "LJ96_Coulomb",
                "classification": "verified_fixed_source_contribution",
                "coefficient_source": "unaltered source atom LJ rows and bond-increment charges; declared isolated research pair policy",
                "functional_form": "w*epsilon_ij*(2*(rmin_ij/r)^9-3*(rmin_ij/r)^6) + w*332.06371*qi*qj/r; dielectric1; no cutoff shifting",
                "mixing": "rmin_ij^6=(ri^6+rj^6)/2; epsilon_ij=2*sqrt(ei*ej)*ri^3*rj^3/(ri^6+rj^6)",
                "numerical_coefficients": {
                    "rmin_i_A": ri,
                    "rmin_j_A": rj,
                    "epsilon_i_kcal_mol": ei / 4.184,
                    "epsilon_j_kcal_mol": ej / 4.184,
                    "qi_e": model["charges"][str(i)],
                    "qj_e": model["charges"][str(j)],
                    "Coulomb_kcal_A_mol_e2": COULOMB,
                    "dielectric": 1,
                },
                "source_equilibrium_dependencies": [],
                "units": "kcal/mol, angstrom, e",
            }
        )
    return signed(
        {
            "internal_identity": IDENTITY,
            "model_identity": model["identity"],
            "table": table,
            "query_to_columns": p["query_to_columns"],
            "directed_missing_query_count": len(
                {r["directed_query"] for r in p["query_to_columns"]}
            ),
            "independent_coefficient_count": len(p["columns"]),
            "columns": p["columns"],
            "family_coefficient_counts": dict(
                Counter(n.split(":")[0] for n in p["columns"])
            ),
            "structural_symmetry_reduction": "AA swap I/L only with J,K fixed. Dihedral reversal exchanges end blocks. Identical source types share coefficients; source-equation symmetry only, no invented cross-environment equality.",
            "all_graph_interactions_accounted": True,
            "independent_graph_enumeration": graph_inventory(model),
            "unresolved_interpretation_dependencies": [],
            "physical_coefficients_missing": p["columns"],
            "parameter_completeness": "incomplete",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def graph_inventory(model):
    """Re-enumerate path/center roles from the graph, not the source row list."""
    adj = adjacency(model)
    bonds, angles, propers, aa = set(), set(), set(), set()
    canonical = lambda s: min(tuple(s), tuple(s[::-1]))
    for j, neighbors in adj.items():
        for i in neighbors:
            bonds.add(canonical([i, j]))
            for k in neighbors:
                if i == k:
                    continue
                angles.add(canonical([i, j, k]))
                for l in adj[k]:
                    if l not in (i, j):
                        propers.add(canonical([i, j, k, l]))
        for k in neighbors:
            for i, l in combinations(sorted(n for n in neighbors if n != k), 2):
                aa.add((i, j, k, l))
    expected = {
        "quartic_bond": bonds,
        "quartic_angle": angles,
        "bond-angle": angles,
        "bond-bond": angles,
        "angle-angle": aa,
    }
    for family in (
        "torsion_3",
        "bond-bond_1_3",
        "end_bond-torsion_3",
        "middle_bond-torsion_3",
        "angle-torsion_3",
        "angle-angle-torsion_1",
    ):
        expected[family] = propers
    counts = {}
    for family, wanted in expected.items():
        found = [
            tuple(r["sites"])
            for r in model["interaction_inventory"]
            if r["family"] == family
        ]
        found = [
            (min(s[0], s[3]), s[1], s[2], max(s[0], s[3]))
            if family == "angle-angle"
            else canonical(s)
            for s in found
        ]
        require(
            len(found) == len(set(found)) and set(found) == wanted,
            "HMDS_graph_interaction_inventory_incomplete",
            family=family,
        )
        counts[family] = len(wanted)
    return {
        "counts": counts,
        "all_required_graph_roles_present": True,
        "enumeration": "All graph bonds, simple3-/4-site paths and every shared-neighbor angle-angle role; independent of source row inventory and missing-query count.",
    }
