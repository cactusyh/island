"""J59 source-complete fixed subset and explicit unbound physical remainder.

This is an audit evaluator, never a complete force-field or residual fit. All
energies are kcal/mol, lengths angstrom, angles radians. Unknowns never enter
fixed sums. Conditional kernels do not select exclusions or phase conventions.
"""

import json
from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path

import numpy as np

from island.forcefields.pcff.class2 import FAMILIES
from island.pcff_completion import digest, require, signed
from island.pcff_qm.geometry import measure
from island.pcff_qm.records import confirmation_inventory

from .campaign import selected

FIXED = ("quartic_bond", "quartic_angle", "bond-bond", "bond-angle")
FORMS = {
    "quartic_bond": "K2 dr^2 + K3 dr^3 + K4 dr^4",
    "quartic_angle": "K2 dt^2 + K3 dt^3 + K4 dt^4",
    "bond-bond": "K dr_left dr_right",
    "bond-angle": "(Kleft dr_left + Kright dr_right) dt",
    "torsion_3": "sum Vn (1+cos(n phi-phase_n)); Phi geometry mapping remains unbound",
    "angle-angle": "K dt_IJK dt_KJL; center J shared arm K; I/L exchange only",
    "bond-bond_1_3": "K dr_IJ dr_KL; path reversal only",
    "end_bond-torsion_3": "dr_IJ sum Fn_left cos(n phi) + dr_KL sum Fn_right cos(n phi)",
    "middle_bond-torsion_3": "dr_JK sum Fn cos(n phi)",
    "angle-torsion_3": "dt_IJK sum Fn_left cos(n phi) + dt_JKL sum Fn_right cos(n phi)",
    "angle-angle-torsion_1": "K dt_IJK dt_JKL (phi-Phi0); source linear phi; main torsion phase1 candidate dependency retained; no cos replacement",
    "nonbond(9-6)": "eps_ij [2 (rmin_ij/r)^9 -3 (rmin_ij/r)^6]; sixth-power mixing; pair weights not supplied",
}


def bond(x, ids):
    a, b = [i - 1 for i in ids]
    v = x[a] - x[b]
    r = float(np.linalg.norm(v))
    require(r > 1e-10, "degenerate_bond")
    d = np.zeros_like(x)
    d[a] = v / r
    d[b] = -d[a]
    return r, d


def angle(x, ids):
    a, b, c = [i - 1 for i in ids]
    u, v = x[a] - x[b], x[c] - x[b]
    ru, rv = np.linalg.norm(u), np.linalg.norm(v)
    require(min(ru, rv) > 1e-10, "degenerate_angle")
    uh, vh = u / ru, v / rv
    co = float(np.clip(uh @ vh, -1, 1))
    si = np.sqrt(1 - co * co)
    require(si > 1e-7, "singular_angle")
    d = np.zeros_like(x)
    d[a] = (co * uh - vh) / (ru * si)
    d[c] = (co * vh - uh) / (rv * si)
    d[b] = -d[a] - d[c]
    return float(np.arccos(co)), d


def fixed_rows(rows, xyz):
    """Analytical energies/forces for the four fully bound source families only."""
    x = np.asarray(xyz, dtype=float)
    require(
        x.ndim == 2 and x.shape[1] == 3 and np.isfinite(x).all(), "coordinate_mutation"
    )
    lookup = {r["id"]: r for r in rows}
    energy = defaultdict(float)
    gradient = np.zeros_like(x)
    details = []

    def delta(kind, sites, rid):
        value, derivative = (bond if kind == "bond" else angle)(x, sites)
        return value - lookup[rid]["assignment"]["values"][0], derivative

    for row in rows:
        f = row["family"]
        a = row["assignment"]
        s = row["sites"]
        if f not in FIXED:
            continue
        require(a is not None, "fixed_coefficient_missing", id=row["id"])
        v = np.asarray(a["values"], float) / 4.184
        if f in ("quartic_bond", "quartic_angle"):
            value, d = (bond if f == "quartic_bond" else angle)(x, s)
            q = value - a["values"][0]
            e = sum(v[j] * q ** (j + 1) for j in (1, 2, 3))
            g = d * sum((j + 1) * v[j] * q**j for j in (1, 2, 3))
        else:
            i, j, k = s
            r, gr = delta("bond", [i, j], row["dependencies"][0])
            t, gt = delta("bond", [j, k], row["dependencies"][1])
            if f == "bond-bond":
                e = v[0] * r * t
                g = v[0] * (gr * t + gt * r)
            else:
                theta, gtheta = delta("angle", s, row["dependencies"][2])
                e = (v[0] * r + v[1] * t) * theta
                g = (v[0] * gr + v[1] * gt) * theta + (v[0] * r + v[1] * t) * gtheta
        energy[f] += float(e)
        gradient += g
        details.append({"id": row["id"], "energy": float(e)})
    return {
        "components": dict(energy),
        "energy": float(sum(energy.values())),
        "forces": (-gradient).tolist(),
        "rows": details,
        "meaning": "verified fixed subset only; complete physical remainder is explicit and unbound",
    }


def validate_model(model, expected_identity=None):
    require(
        model["identity"]
        == digest({k: v for k, v in model.items() if k != "identity"}),
        "residual_model_mutation",
    )
    if expected_identity:
        require(
            model["identity"] == expected_identity, "residual_external_binding_mutation"
        )
    require(
        model["energy_unit"] == "kcal/mol"
        and model["length_unit"] == "angstrom"
        and model["angle_unit"] == "radian",
        "residual_units_mutation",
    )
    require(
        model["family"] == "pcff" and model["atom_order"] == list(range(1, 10)),
        "residual_mapping_mutation",
    )
    return model


def evaluate(model, xyz, *, expected_identity=None, require_complete=False):
    validate_model(model, expected_identity)
    require(
        not require_complete,
        "unbound_physical_remainder",
        unresolved=[r["id"] for r in model["table"] if r["status"] == "unresolved"],
    )
    require(np.asarray(xyz).shape == (9, 3), "residual_mapping_mutation")
    return fixed_rows(model["interaction_inventory"], xyz)


def lj96(r, ri, rj, ei, ej):
    """Unweighted source pair kernel; never an assigned pair energy."""
    rm = ((ri**6 + rj**6) / 2) ** (1 / 6)
    eps = 2 * np.sqrt(ei * ej) * ri**3 * rj**3 / (ri**6 + rj**6)
    e = eps * (2 * (rm / r) ** 9 - 3 * (rm / r) ** 6)
    dedr = eps * (-18 * (rm / r) ** 9 + 18 * (rm / r) ** 6) / r
    return float(e), float(dedr)


def build(repo, distribution):
    old = Path(repo) / "docs/evidence/phase_4j58"
    scan = json.loads((old / "plan/angle-angle.json").read_text())
    _, assessment, _ = selected(repo, distribution, scan["center_coordinates"])
    p = assessment.payload
    rows = p["interaction_inventory"]
    missing = {r["id"]: r for r in p["missing_terms"]}
    source_equations = p["conventions"]["functional_forms"]["value"]
    table = []
    for row in rows:
        f = row["family"]
        a = row["assignment"]
        fixed = f in FIXED or f in ("atom_type", "mass", "charge")
        table.append(
            {
                **row,
                "status": "fixed" if fixed and a is not None else "unresolved",
                "functional_form": FORMS.get(f, f),
                "coefficient": a.get("values", a.get("value")) if a else None,
                "units": a.get("units", a.get("unit"))
                if a
                else FAMILIES.get(f, (None, None, ["not_applicable"]))[2],
                "roles": row["sites"],
                "permutations": "I/L exchange with J center K shared"
                if f == "angle-angle"
                else "path reversal with end-block exchange; no arbitrary permutation",
                "query": (a or {}).get(
                    "query", missing.get(row["id"], {}).get("query")
                ),
                "source_evidence": source_equations.get(f),
                "provenance": (a or {}).get("row_ids", []),
                "source_equilibrium_dependencies": row["dependencies"],
                "reason": "bound source coefficients and converter equation"
                if fixed
                else "numerical coefficient absent; interpretation confirmation is not a value"
                if a is None
                else "coefficient present; torsion Phi/phase geometry unbound"
                if f == "torsion_3"
                else "atom LJ rows present; pair policy unbound",
            }
        )
    for row in table:
        if row["family"] == "angle-angle-torsion_1":
            main = next(
                r
                for r in rows
                if r["family"] == "torsion_3" and r["sites"] == row["sites"]
            )
            row["source_phase_dependency"] = {
                "source_symbol": "Phi1(0)",
                "candidate_main_torsion_id": main["id"],
                "main_phase1_radian": main["assignment"]["values"][1],
                "source_rows": main["assignment"]["row_ids"],
                "status": "numerical main phase exists; dependency and Phi geometry mapping must be explicit, not invented absent number",
            }
    for r in p["improper_structural_inventory"]:
        table.append(
            {
                "id": r["key"],
                "family": "wilson_out_of_plane",
                "sites": [int(i) for i in r["key"].split(":")[1].split(",")],
                "status": "not_applicable",
                "coefficient": None,
                "units": "kcal/mol/radian^2",
                "functional_form": "K(chi-chi0)^2",
                "reason": r["reason"],
                "authority": "existing structural converter, not invented author omission; degree-four centers separately contain all AA terms",
            }
        )
    for left in (4, 5, 6):
        for right in (7, 8, 9):
            table.append(
                {
                    "id": f"torsion-torsion_1:{left},1,2,3,{right}",
                    "family": "torsion-torsion_1",
                    "sites": [left, 1, 2, 3, right],
                    "status": "unresolved",
                    "coefficient": None,
                    "units": "kcal/mol",
                    "functional_form": "F cos(phi_HSiOSi) cos(phi_SiOSiH)",
                    "reason": "five-site path exists; FRC equation with empty section and no omission authority",
                    "source_line": 5521,
                }
            )
    table.append(
        {
            "id": "improper_domain_applicability",
            "family": "improper_applicability",
            "sites": None,
            "status": "unresolved",
            "coefficient": None,
            "units": "applicability rule",
            "functional_form": "Wilson out-of-plane K(chi-chi0)^2 if source requires it",
            "reason": "eight converter structural not_applicable rows are not source omission authority; primary domain-specific Wilson inclusion/omission rule still needed",
            "provenance": p["conventions"]["improper_applicability"],
        }
    )
    graph = scan["target"]["graph"]
    adj = {i: [] for i in p["atom_order"]}
    for a, b, *_ in graph["bonds"]:
        adj[a].append(b)
        adj[b].append(a)
    pair_rows = []
    for i, j in combinations(p["atom_order"], 2):
        todo = deque([(i, 0)])
        seen = set()
        while todo:
            k, d = todo.popleft()
            if k == j:
                distance = d
                break
            if k in seen:
                continue
            seen.add(k)
            todo.extend((q, d + 1) for q in adj[k] if q not in seen)
        assignment_i = next(
            r["assignment"]["values"] for r in rows if r["id"] == f"nonbond(9-6):{i}"
        )
        assignment_j = next(
            r["assignment"]["values"] for r in rows if r["id"] == f"nonbond(9-6):{j}"
        )
        pairs = {
            "id": f"pair:{i},{j}",
            "family": "LJ9-6+Coulomb",
            "sites": [i, j],
            "status": "unresolved",
            "topological_distance": distance,
            "special_pair_class": f"1-{distance + 1}",
            "coefficient": {
                "ri_ei": assignment_i,
                "rj_ej": assignment_j,
                "qi": p["charge_vector"][str(i)],
                "qj": p["charge_vector"][str(j)],
                "LJ_weight": None,
                "Coulomb_weight": None,
                "Coulomb_prefactor": None,
                "dielectric": None,
            },
            "units": {
                "LJ_rows": "angstrom,kJ/mol",
                "charge": "e",
                "Coulomb_kernel": "e^2/angstrom",
            },
            "functional_form": "wLJ*LJ96(r) + wQ*prefactor*qi*qj/(dielectric*r); cutoff/switch/image policy unbound",
            "exclusion": None,
            "reason": "FRC supplies no source pair exclusions/scalings or Coulomb execution convention; not assigned zero",
            "provenance": [
                f"nonbond(9-6):{i}",
                f"nonbond(9-6):{j}",
                f"charge:{i}",
                f"charge:{j}",
            ],
        }
        pair_rows.append(pairs)
    table += pair_rows
    confirmation = confirmation_inventory()
    confirmed_queries = {r["query"] for r in confirmation["confirmations"]}
    for row in table:
        if row.get("query") in confirmed_queries:
            row["interpretation_status"] = "author_confirmed_audited"
            row["interpretation_authority"] = "user_reported_author_confirmation"
            row["original_author_artifact_status"] = "not_yet_supplied"
            row["numerical_coefficient_status"] = (
                "missing" if row["coefficient"] is None else "supplied"
            )
    result = {
        "schema": "island.j59-residual-audit.v1",
        "family": "pcff",
        "graph_identity": p["graph_identity"],
        "graph": graph,
        "center_coordinates": scan["center_coordinates"],
        "atom_order": p["atom_order"],
        "source_hashes": p["source_hashes"],
        "source_rows": p["source_rows"],
        "assessment_identity": p["identity"],
        "charge_identity": p["charge_identity"],
        "charge_vector": p["charge_vector"],
        "typing": p["typing"],
        "energy_unit": "kcal/mol",
        "length_unit": "angstrom",
        "angle_unit": "radian",
        "interaction_inventory": rows,
        "table": table,
        "status_counts": dict(Counter(r["status"] for r in table)),
        "author_confirmations": confirmation_inventory(),
        "complete_residual": False,
        "unknowns_set_to_zero": False,
        "parameter_completeness": "incomplete",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    return signed(result)


def reconstruct(model, repo, distribution):
    validate_model(model)
    require(model == build(repo, distribution), "source_mapping_or_charge_mutation")
    return model["identity"]


def affected(model, samples):
    """Every table row retains coordinate/dependency ranges for both scans."""
    out = []
    lookup = {r["id"]: r for r in model["interaction_inventory"]}
    for r in model["table"]:
        ids = r.get("sites")
        ranges = {}
        coordinates = []
        if ids and len(ids) > 1:
            if r["family"] == "angle-angle":
                i, j, k, l = ids
                coordinates = [("angle", [i, j, k]), ("angle", [k, j, l])]
            elif r["family"] == "LJ9-6+Coulomb":
                coordinates = [("bond", ids)]
            elif len(ids) == 5:
                coordinates = [("torsion", ids[:4]), ("torsion", ids[1:])]
            elif len(ids) in (2, 3, 4):
                coordinates = [({2: "bond", 3: "angle", 4: "torsion"}[len(ids)], ids)]
        for kind, sites in coordinates:
            values = [
                measure(
                    s["coordinates"],
                    model["atom_order"],
                    {"kind": kind, "sites": sites},
                )
                for s in samples
            ]
            ranges[kind + ":" + ",".join(map(str, sites))] = float(
                np.ptp(np.unwrap(values) if kind == "torsion" else values)
            )
        for dep in r.get("dependencies", []):
            d = lookup[dep]
            if len(d["sites"]) < 2:
                continue
            kind = "bond" if len(d["sites"]) == 2 else "angle"
            ranges[dep] = float(
                np.ptp(
                    [
                        measure(
                            s["coordinates"],
                            model["atom_order"],
                            {"kind": kind, "sites": d["sites"]},
                        )
                        for s in samples
                    ]
                )
            )
        out.append(
            {
                "id": r["id"],
                "status": r["status"],
                "ranges": ranges,
                "affected": any(v > 1e-8 for v in ranges.values()),
                "force_contribution_may_be_nonzero_even_if_energy_constant": True,
            }
        )
    return out


def source_torsion_kernel(phi, values):
    """Literal FRC scalar equation, not an authorized atom-to-Phi mapping."""
    v = np.asarray(values, float)
    return float(
        sum(v[2 * (n - 1)] * (1 + np.cos(n * phi - v[2 * n - 1])) for n in (1, 2, 3))
        / 4.184
    )


def conditional_columns(model, xyz):
    """Source-equation nuisance diagnostic under explicitly conditional Phi.

    Includes all missing coupling families in this graph. Phase/absence/pair
    policy choices remain unbound; this is never evidence of physical closure.
    Pair grouping by graph separation is a *lower-bound diagnostic hypothesis*,
    not authority to tie source pair coefficients or select an exclusion rule.
    """
    x = np.asarray(xyz, float)
    order = model["atom_order"]
    rows = model["interaction_inventory"]
    lookup = {r["id"]: r for r in rows}
    out = defaultdict(float)
    types = model["typing"]

    def coord(kind, sites):
        value = measure(x, order, {"kind": kind, "sites": sites})
        if kind == "torsion":
            reference = measure(
                model["center_coordinates"], order, {"kind": kind, "sites": sites}
            )
            value = reference + np.arctan2(
                np.sin(value - reference), np.cos(value - reference)
            )
        return value

    def equilibrium(family, sites):
        ids = min(tuple(sites), tuple(sites[::-1]))
        r = lookup[family + ":" + ",".join(map(str, ids))]
        return r["assignment"]["values"][0]

    def dr(sites):
        return coord("bond", sites) - equilibrium("quartic_bond", sites)

    def dt(sites):
        return coord("angle", sites) - equilibrium("quartic_angle", sites)

    for row in rows:
        f = row["family"]
        s = row["sites"]
        a = row["assignment"]
        if f == "angle-angle" and a is None:
            i, j, k, l = s
            labels = tuple(types[str(q)] for q in s)
            rev = (labels[3], labels[1], labels[2], labels[0])
            key = f + ":" + ",".join(min(labels, rev))
            out[key] += dt([i, j, k]) * dt([k, j, l])
        elif (
            f
            in (
                "bond-bond_1_3",
                "end_bond-torsion_3",
                "middle_bond-torsion_3",
                "angle-torsion_3",
                "angle-angle-torsion_1",
            )
            and a is None
        ):
            labels = tuple(types[str(q)] for q in s)
            if labels[::-1] < labels:
                s = s[::-1]
                labels = labels[::-1]
            i, j, k, l = s
            phi = coord("torsion", s)
            prefix = f + ":" + ",".join(labels)
            if f == "bond-bond_1_3":
                out[prefix] += dr([i, j]) * dr([k, l])
            elif f == "angle-angle-torsion_1":
                product = dt([i, j, k]) * dt([j, k, l])
                out[prefix + ":Kphi"] += product * phi
                out[prefix + ":Kphase_product"] += product
            else:
                deformations = (
                    [dr([i, j]), dr([k, l])]
                    if f == "end_bond-torsion_3"
                    else [dr([j, k])]
                    if f == "middle_bond-torsion_3"
                    else [dt([i, j, k]), dt([j, k, l])]
                )
                for end, d in enumerate(deformations):
                    for n in (1, 2, 3):
                        out[prefix + f":block{end}:n{n}"] += d * np.cos(n * phi)
        elif f == "torsion_3":
            out["conditional_main_torsion_source_phi_amplitude"] += (
                source_torsion_kernel(coord("torsion", s), a["values"])
            )
    for r in model["table"]:
        if r["family"] == "torsion-torsion_1":
            s = r["sites"]
            out["conditional_torsion_torsion_empty_source"] += np.cos(
                coord("torsion", s[:4])
            ) * np.cos(coord("torsion", s[1:]))
        elif r["family"] == "LJ9-6+Coulomb":
            c = r["coefficient"]
            v = coord("bond", r["sites"])
            ri, ei = c["ri_ei"]
            rj, ej = c["rj_ej"]
            d = r["topological_distance"]
            out[f"conditional_LJ_unit_weight_path{d}"] += lj96(
                v, ri, rj, ei / 4.184, ej / 4.184
            )[0]
            out[f"conditional_Coulomb_qiqj_over_r_path{d}"] += (
                float(c["qi"]) * float(c["qj"]) / v
            )
    return dict(out)
