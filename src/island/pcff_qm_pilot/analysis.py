"""Source-bound pilot identifiability, sensitivity and projected Hessian analysis.

Partial bonded residuals are explicitly incomplete; no omitted term is zero.
"""

import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

from island.pcff_completion import _design, require
from island.pcff_qm.geometry import measure
from island.pcff_qm.reference import CONVERSION

from .campaign import NAMES, selected
from .importer import native


def spectrum(matrix):
    a = np.asarray(matrix, dtype=float)
    norms = np.linalg.norm(a, axis=0)
    active = norms > 1e-12
    scaled = a[:, active] / norms[active] if active.any() else a[:, :0]
    sv = np.linalg.svd(scaled, compute_uv=False)
    rank = int(np.sum(sv > (sv[0] * 1e-10 if len(sv) else 0)))
    return {
        "shape": list(a.shape),
        "column_norms": norms.tolist(),
        "scaled_singular_values": sv.tolist(),
        "rank": rank,
        "active_columns": int(active.sum()),
        "conditioning": float(sv[0] / sv[-1])
        if len(sv) and rank == a.shape[1]
        else None,
        "scaling": "unit 2-norm columns; rank threshold 1e-10 relative; zero-column threshold 1e-12",
    }


def hessian_projection(xyz, hessian):
    """Remove rigid Cartesian translations/rotations, retaining all eigenvalues."""
    xyz, h = np.asarray(xyz, dtype=float), np.asarray(hessian, dtype=float)
    n = len(xyz)
    r = xyz - xyz.mean(axis=0)
    modes = np.column_stack(
        [np.tile(a, (n, 1)).ravel() for a in np.eye(3)]
        + [np.cross(np.tile(a, (n, 1)), r).ravel() for a in np.eye(3)]
    )
    u, sv, _ = np.linalg.svd(modes, full_matrices=True)
    rank = int(np.sum(sv > sv[0] * 1e-10))
    q = u[:, rank:]
    projected = np.linalg.eigvalsh(q.T @ ((h + h.T) / 2) @ q)
    return {
        "method": "Euclidean Cartesian orthogonal projection of rigid translations and rotations; not a mass-weighted frequency calculation",
        "rigid_mode_rank": rank,
        "vibrational_dimension": len(projected),
        "unprojected_eigenvalues_hartree_per_bohr2": np.linalg.eigvalsh(h).tolist(),
        "projected_eigenvalues_hartree_per_bohr2": projected.tolist(),
        "negative_projected_modes_below_minus_1e_minus5": int(
            (projected < -1e-5).sum()
        ),
        "strict_full_positive_definiteness_required": False,
    }


def symmetry_key(query):
    family, labels = query.split(":")
    a, b, c, d = labels.split(",")
    reverse = (d, b, c, a) if family == "angle-angle" else (d, c, b, a)
    return family + ":" + ",".join(min((a, b, c, d), reverse))


def partial_bonded(p, xyz):
    """Verified quartic and adjacent coupling rows only, kcal/mol.

    Torsions and nonbonded interactions are not included and are not zero.
    """
    order = p["atom_order"]
    rows = {r["id"]: r for r in p["interaction_inventory"]}

    def coord(kind, ids):
        return measure(xyz, order, {"kind": kind, "sites": ids})

    result = defaultdict(float)
    for r in rows.values():
        assignment, family, sites = r["assignment"], r["family"], r["sites"]
        if assignment is None:
            continue
        values = assignment.get("values")
        if family in ("quartic_bond", "quartic_angle"):
            kind = "bond" if family == "quartic_bond" else "angle"
            delta = coord(kind, sites) - values[0]
            result[family] += (
                sum(values[j] * delta ** (j + 1) for j in (1, 2, 3)) / 4.184
            )
        elif family in ("bond-bond", "bond-angle"):
            i, j, k = sites
            dr = [
                coord("bond", ids) - rows[key]["assignment"]["values"][0]
                for ids, key in zip(
                    ([i, j], [j, k]), r["dependencies"][:2], strict=True
                )
            ]
            if family == "bond-bond":
                energy = values[0] * dr[0] * dr[1]
            else:
                dt = (
                    coord("angle", sites)
                    - rows[r["dependencies"][2]]["assignment"]["values"][0]
                )
                energy = (values[0] * dr[0] + values[1] * dr[1]) * dt
            result[family] += energy / 4.184
    return dict(result)


def partial_forces(p, xyz):
    x = np.array(xyz, dtype=float)
    out = np.zeros_like(x)
    step = 1e-5
    for i in range(len(x)):
        for j in range(3):
            plus, minus = x.copy(), x.copy()
            plus[i, j] += step
            minus[i, j] -= step
            out[i, j] = -(
                sum(partial_bonded(p, plus).values())
                - sum(partial_bonded(p, minus).values())
            ) / (2 * step)
    return out


def identifiability(repo, distribution, evidence, normalized):
    evidence, normalized = Path(evidence), Path(normalized)
    scan = json.loads((evidence / "plan/angle-angle.json").read_text())
    system, assessment, _ = selected(repo, distribution, scan["center_coordinates"])
    p = assessment.payload
    missing = {r["id"]: r for r in p["missing_terms"]}
    rows = p["interaction_inventory"]
    numerical = [
        r for r in p["missing_terms"] if r["family"] in ("angle-angle", "bond-bond_1_3")
    ]
    columns = sorted({symmetry_key(r["query"]) for r in numerical})
    records, matrices, force_matrices, residuals, affected = [], [], [], [], []
    for name in NAMES:
        scan = json.loads((evidence / "plan" / (name + ".json")).read_text())
        samples = json.loads((normalized / name / "samples.json").read_text())
        matrix = np.zeros((25, len(columns)))
        force_matrix = np.zeros((25, 27, len(columns)))
        for row in numerical:
            idx = columns.index(symmetry_key(row["query"]))
            matrix[:, idx] += np.asarray(_design(system, p, row, samples)).ravel()
            # Force Jacobian from this exact source-equation column, not QM data.
            for atom in range(9):
                for axis in range(3):
                    plus, minus = [], []
                    for sample in samples:
                        xp, xm = (
                            np.array(sample["coordinates"]),
                            np.array(sample["coordinates"]),
                        )
                        xp[atom, axis] += 1e-5
                        xm[atom, axis] -= 1e-5
                        plus.append({"coordinates": xp.tolist()})
                        minus.append({"coordinates": xm.tolist()})
                    force_matrix[:, atom * 3 + axis, idx] += (
                        -(
                            np.asarray(_design(system, p, row, plus)).ravel()
                            - np.asarray(_design(system, p, row, minus)).ravel()
                        )
                        / 2e-5
                    )
        matrices.append(matrix - matrix[12])
        force_matrices.append(force_matrix.reshape(-1, len(columns)))
        fixed = [partial_bonded(p, s["coordinates"]) for s in samples]
        partial = np.array([sum(v.values()) for v in fixed])
        energy = np.array([s["energy"] for s in samples])
        unexplained = (energy - energy[12]) - (partial - partial[12])
        forces = np.array([s["forces"] for s in samples])
        partial_force = np.array([partial_forces(p, s["coordinates"]) for s in samples])
        residual_force = forces - partial_force
        residuals.append(
            {
                "scan": name,
                "meaning": "QM minus ONLY verified quartic/adjacent-coupling contribution; all other terms and model discrepancy remain, NOT a fitted cross-term residual",
                "QM_relative_energy_kcal_mol": (energy - energy[12]).tolist(),
                "verified_partial_energy_components_kcal_mol": fixed,
                "partial_relative_energy_kcal_mol": (partial - partial[12]).tolist(),
                "unexplained_energy_kcal_mol": unexplained.tolist(),
                "unexplained_energy_rms_kcal_mol": float(
                    np.sqrt(np.mean(unexplained**2))
                ),
                "unexplained_force_rms_kcal_mol_angstrom": float(
                    np.sqrt(np.mean(residual_force**2))
                ),
                "unexplained_forces": residual_force.tolist(),
                "fitted_energy_error": None,
                "fitted_force_error": None,
            }
        )
        for row in rows:
            family, sites = row["family"], row["sites"]
            if family in ("atom_type", "mass", "charge", "nonbond(9-6)"):
                continue
            kind = (
                "bond"
                if len(sites) == 2
                else "angle"
                if len(sites) == 3
                else "torsion"
                if family != "angle-angle"
                else "angle_pair"
            )
            if kind == "angle_pair":
                i, j, k, l = sites
                values = [
                    [
                        measure(
                            s["coordinates"],
                            p["atom_order"],
                            {"kind": "angle", "sites": ids},
                        )
                        for ids in ([i, j, k], [k, j, l])
                    ]
                    for s in samples
                ]
            else:
                values = [
                    [
                        measure(
                            s["coordinates"],
                            p["atom_order"],
                            {"kind": kind, "sites": sites},
                        )
                    ]
                    for s in samples
                ]
            dependency_ranges = {}
            for dependency in row["dependencies"]:
                dep = next(r for r in rows if r["id"] == dependency)
                dep_kind = "bond" if len(dep["sites"]) == 2 else "angle"
                dep_values = [
                    measure(
                        sample["coordinates"],
                        p["atom_order"],
                        {"kind": dep_kind, "sites": dep["sites"]},
                    )
                    for sample in samples
                ]
                dependency_ranges[dependency] = float(np.ptp(dep_values))
            if (
                max(
                    float(np.max(np.ptp(values, axis=0))),
                    *dependency_ranges.values(),
                    0.0,
                )
                > 1e-8
            ):
                affected.append(
                    {
                        "scan": name,
                        "id": row["id"],
                        "family": family,
                        "sites": sites,
                        "coordinate_ranges": np.ptp(values, axis=0).tolist(),
                        "dependency_coordinate_ranges": dependency_ranges,
                        "source_assignment": row["assignment"],
                        "query": missing.get(row["id"], {}).get("query"),
                        "numerical_row_status": "source_supplied"
                        if row["assignment"]
                        else "not_supplied",
                        "confirmation_scope": "user-attested interpretation; not a supplied numerical coefficient"
                        if row["assignment"] is None and "torsion" in family
                        else None,
                    }
                )
        for i, j in combinations(range(9), 2):
            distances = [
                float(
                    np.linalg.norm(
                        np.array(s["coordinates"])[i] - np.array(s["coordinates"])[j]
                    )
                )
                for s in samples
            ]
            if np.ptp(distances) > 1e-8:
                affected.append(
                    {
                        "scan": name,
                        "id": f"nonbond_pair:{i + 1},{j + 1}",
                        "family": "LJ9-6_and_Coulomb",
                        "sites": [i + 1, j + 1],
                        "distance_range_angstrom": float(np.ptp(distances)),
                        "status": "source atom rows/charges exist; exclusions, special pairs, electrostatic execution contract must be bound before a residual subtraction; not assumed zero",
                    }
                )
        records.append(
            {
                "scan": name,
                "single_target_design": spectrum(
                    np.asarray(
                        _design(system, p, scan["target"]["target_site"], samples)
                    )
                ),
                "joint_numerical_energy_design": spectrum(matrices[-1]),
                "joint_numerical_energy_force_design": spectrum(
                    np.vstack([matrices[-1], force_matrices[-1]])
                ),
                "held_out_v2": scan["held_out_v2"],
            }
        )
    matrix = np.vstack(matrices)
    # Formal counterexample: unknown model discrepancy per observation can absorb
    # any target coefficient. This is an identifiability demonstration, not a fit.
    augmented = np.column_stack([matrix, np.eye(len(matrix))])
    grouped = defaultdict(list)
    inventory = json.loads(
        (evidence.parent / "phase_4j57/audit/qm-target-inventory.json").read_text()
    )
    for entry in inventory["targets"]:
        grouped[symmetry_key(entry["query"])].append(entry["query"])
    result = {
        "schema": "island.j58-identifiability.v1",
        "model": "QM = fixed verified MM + explicitly fitted unknown MM + model discrepancy",
        "directed_queries": 19,
        "source_role_symmetry_groups": dict(grouped),
        "symmetry_group_count": len(grouped),
        "symmetry_policy": "angle-angle swap I/L with J center and K shared arm; 1-3 bond coupling full path reversal only; no cross-fixture transfer inferred",
        "source_equations": [
            scan["target"]["functional_form"]
            for scan in [
                json.loads((evidence / "plan" / (n + ".json")).read_text())
                for n in NAMES
            ]
        ],
        "source_units": "Kaa kcal/mol/radian^2; Kbb13 kcal/mol/angstrom^2; no half prefactor; equilibrium angle degrees converted once to radians by J54",
        "coefficient_columns_for_exact_graph": columns,
        "joint_numerical_site_counts": dict(Counter(r["query"] for r in numerical)),
        "joint_design_values": matrix.tolist(),
        "joint_energy_design": spectrum(matrix),
        "joint_energy_force_design": spectrum(np.vstack([matrix, *force_matrices])),
        "per_scan": records,
        "discrepancy_counterexample": {
            "matrix_shape": list(augmented.shape),
            "rank": int(np.linalg.matrix_rank(augmented)),
            "nullity": augmented.shape[1] - int(np.linalg.matrix_rank(augmented)),
            "interpretation": "nonzero target changes offset by unknown discrepancy leave QM predictions identical; numerical subdesign rank cannot prove a complete residual is identified",
        },
        "affected_interactions": affected,
        "affected_counts": dict(Counter(r["family"] for r in affected)),
        "source_structural_improper_inventory": p["improper_structural_inventory"],
        "fixed_partial_force_method": "central differences of the verified source polynomials, 1e-5 angstrom; these are diagnostic partial contributions, not a complete backend force",
        "known_partial_families": [
            "quartic_bond",
            "quartic_angle",
            "bond-bond",
            "bond-angle",
        ],
        "other_contributions_not_zero": [
            "source-supplied main torsions require an explicit interpreted evaluator",
            "all unprovided angle-angle and 1-3 bond couplings at symmetry-related sites",
            "attested torsion-coupling interpretations still supply no missing numeric coefficients",
            "LJ9-6 and Coulomb including all pair exclusions/special-pair policies",
            "electronic/model discrepancy; pilot QM center need not equal FRC equilibria",
        ],
        "partial_residuals": residuals,
        "fit_readiness": "qm_data_available",
        "identifiability": "not_identifiable_as_one_missing_term",
        "fit_completed": False,
        "fit_validated": False,
        "created_parameters": [],
        "next_data": [
            "source-bound full numerical torsion-coupling and nonbonded execution contract, with independent evaluator checks",
            "joint, symmetry-reduced model for all affected unknown coefficients; versioned targets before fitting",
            "independent third-arm/angular and central-Si-O deformations, torsion rotations and mixed bond/angle scans to separate nuisance cross terms",
            "held-out independent conformers and HMDS/DP2 methyl environments before transfer assessment",
            "quantified QM method/basis/grid and model discrepancy uncertainty; predeclared fit error criteria",
        ],
        "held_out_review": "historical column split is not independent validation; J58 eight corners/edge midpoints is a new coverage split, still not independent conformers; no fit or transferable validation claimed",
        "parameter_completeness": "incomplete",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    return result


def diagnostics(evidence, work):
    evidence, work = Path(evidence), Path(work)

    def read(key):
        req = json.loads((work / "requests" / (key + ".json")).read_text())
        return native(work / "outputs" / key, req)

    basekeys = {
        "center": "angle-angle/point-0012",
        "angle-corner": "angle-angle/point-0000",
        "bond-corner": "bond-bond/point-0024",
    }
    sensitivity = []
    for method, basis in (("pbe0", "def2-tzvp"), ("b3lyp", "def2-svp")):
        center, _ = read("sensitivity/center-" + method + "-" + basis)
        bc, _ = read(basekeys["center"])
        for name, basekey in basekeys.items():
            r, a = read("sensitivity/" + name + "-" + method + "-" + basis)
            b, ba = read(basekey)
            sensitivity.append(
                {
                    "geometry": name,
                    "method": method,
                    "basis": basis,
                    "energy_hartree": r["energy_hartree"],
                    "relative_energy_kcal_mol": (
                        r["energy_hartree"] - center["energy_hartree"]
                    )
                    * CONVERSION["energy_factor"],
                    "relative_energy_difference_to_pilot_kcal_mol": (
                        (r["energy_hartree"] - center["energy_hartree"])
                        - (b["energy_hartree"] - bc["energy_hartree"])
                    )
                    * CONVERSION["energy_factor"],
                    "gradient_difference_rms_kcal_mol_angstrom": float(
                        np.sqrt(np.mean((a["gradient"] - ba["gradient"]) ** 2))
                    )
                    * CONVERSION["energy_factor"]
                    / CONVERSION["distance_factor"],
                }
            )
    _, reference = read(basekeys["angle-corner"])
    fd = []
    for h in (0.001, 0.0005):
        plus, _ = read(f"finite-difference/h{h:g}-sign+1")
        minus, _ = read(f"finite-difference/h{h:g}-sign-1")
        force = (
            -(plus["energy_hartree"] - minus["energy_hartree"])
            * CONVERSION["energy_factor"]
            / (2 * h)
        )
        analytic = (
            -reference["gradient"][3, 0]
            * CONVERSION["energy_factor"]
            / CONVERSION["distance_factor"]
        )
        fd.append(
            {
                "step_angstrom": h,
                "atom_row_zero_based": 3,
                "axis": "x",
                "positive_energy_hartree": plus["energy_hartree"],
                "negative_energy_hartree": minus["energy_hartree"],
                "central_difference_force_kcal_mol_angstrom": force,
                "analytical_force_kcal_mol_angstrom": float(analytic),
                "absolute_error": abs(float(force - analytic)),
                "sign_reversed_error": abs(float(force + analytic)),
                "predeclared_tolerance": 0.01,
                "passed": bool(abs(force - analytic) <= 0.01),
            }
        )
    hessians = {}
    for name in NAMES:
        raw, arrays = read(name + "/hessian")
        hessians[name] = hessian_projection(raw["geometry_bohr"], arrays["hessian"])
    require(all(row["passed"] for row in fd), "qm_finite_difference_failed")
    return {
        "schema": "island.j58-numerical-checks.v1",
        "sensitivity": sensitivity,
        "finite_difference": fd,
        "hessians": hessians,
        "method_adequacy": "exploratory only; examine deformation-energy and gradient differences, not SCF convergence alone; three geometries cannot establish final QM adequacy",
        "independent_claims": {
            "QM_execution": "passed",
            "QM_data_integrity": "passed",
            "parameter_completeness": "incomplete",
            "fit_completion": "not_run",
            "fit_validation": "not_run",
            "lammps_runtime": "not_run",
            "openmm_build": "not_run",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
        },
    }
