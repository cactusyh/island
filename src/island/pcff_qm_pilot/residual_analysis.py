"""Reconstruct the fixed subset, nuisance sensitivities and quantitative no-fit gate."""

import json
from collections import Counter
from pathlib import Path

import numpy as np

from island.pcff_completion import require, signed

from .analysis import spectrum
from .residual import (
    affected,
    conditional_columns,
    evaluate,
    fixed_rows,
    lj96,
    source_torsion_kernel,
)


def controls():
    """Small independent analytic controls, never physical parameter insertion."""
    rows = [
        {
            "id": "quartic_bond:1,2",
            "family": "quartic_bond",
            "sites": [1, 2],
            "assignment": {"values": [1.0, 2 * 4.184, 3 * 4.184, 4 * 4.184]},
            "dependencies": [],
        }
    ]
    x = [[0, 0, 0], [1.1, 0, 0]]
    r = fixed_rows(rows, x)
    expected = 2 * 0.1**2 + 3 * 0.1**3 + 4 * 0.1**4
    force = 4 * 0.1 + 9 * 0.1**2 + 16 * 0.1**3
    require(
        abs(r["energy"] - expected) < 1e-14 and abs(r["forces"][0][0] - force) < 1e-14,
        "analytic_bond_control",
    )
    theta0 = np.pi / 2 - 0.1
    a = {
        "id": "quartic_angle:1,2,3",
        "family": "quartic_angle",
        "sites": [1, 2, 3],
        "assignment": {"values": [theta0, 2 * 4.184, 3 * 4.184, 4 * 4.184]},
        "dependencies": [],
    }
    ar = fixed_rows([a], [[1, 0, 0], [0, 0, 0], [0, 1, 0]])
    require(
        abs(ar["energy"] - expected) < 1e-14
        and abs(ar["forces"][0][1] - force) < 1e-14,
        "analytic_angle_control",
    )
    ev, der = lj96(4, 4, 4, 0.5, 0.5)
    require(abs(ev + 0.5) < 1e-14 and abs(der) < 1e-14, "analytic_LJ_rmin_control")
    torsion = source_torsion_kernel(0, [4.184, 0, 0, 0, 0, 0])
    require(abs(torsion - 2) < 1e-14, "literal_source_torsion_prefactor_control")
    # Explicit test-only coefficients, including zero base polynomials to
    # isolate both coupling forces. These are never source parameter rows.
    testrows = [
        {
            "id": "b1",
            "family": "quartic_bond",
            "sites": [1, 2],
            "assignment": {"values": [1.0, 0.0, 0.0, 0.0]},
            "dependencies": [],
        },
        {
            "id": "b2",
            "family": "quartic_bond",
            "sites": [2, 3],
            "assignment": {"values": [1.0, 0.0, 0.0, 0.0]},
            "dependencies": [],
        },
        {
            "id": "a",
            "family": "quartic_angle",
            "sites": [1, 2, 3],
            "assignment": {"values": [np.pi / 2, 0.0, 0.0, 0.0]},
            "dependencies": [],
        },
        {
            "id": "bb",
            "family": "bond-bond",
            "sites": [1, 2, 3],
            "assignment": {"values": [7 * 4.184]},
            "dependencies": ["b1", "b2"],
        },
        {
            "id": "ba",
            "family": "bond-angle",
            "sites": [1, 2, 3],
            "assignment": {"values": [11 * 4.184, 13 * 4.184]},
            "dependencies": ["b1", "b2", "a"],
        },
    ]
    theta = np.pi / 2 + 0.1
    xy = np.array(
        [[1.1, 0, 0], [0, 0, 0], [1.2 * np.cos(theta), 1.2 * np.sin(theta), 0]]
    )
    coupling = fixed_rows(testrows, xy)
    expected_first = [-2.5, 3.7 / 1.1, 0]
    expected_last = -(
        2 * np.array([np.cos(theta), np.sin(theta), 0])
        + 3.7 * np.array([-np.sin(theta) / 1.2, np.cos(theta) / 1.2, 0])
    )
    require(
        abs(coupling["components"]["bond-bond"] - 0.14) < 1e-14
        and abs(coupling["components"]["bond-angle"] - 0.37) < 1e-14
        and np.max(np.abs(np.array(coupling["forces"])[0] - expected_first)) < 1e-13
        and np.max(np.abs(np.array(coupling["forces"])[2] - expected_last)) < 1e-13,
        "analytic_cross_coupling_control",
    )
    return {
        "passed": True,
        "synthetic_control_only_not_source_rows": True,
        "bond_bond_energy": coupling["components"]["bond-bond"],
        "bond_angle_energy": coupling["components"]["bond-angle"],
        "coupling_forces": coupling["forces"],
        "quartic_bond_energy": r["energy"],
        "quartic_bond_force": r["forces"][0][0],
        "quartic_angle_energy": ar["energy"],
        "quartic_angle_force": ar["forces"][0][1],
        "LJ_rmin_energy": ev,
        "LJ_rmin_derivative": der,
        "source_plus_cos_torsion_at_zero": torsion,
    }


def checks(model, plan):
    from scripts.independent_phase_4j59 import numerical_force, source_energy

    records = []
    for name, xyz in plan["geometries"].items():
        actual = evaluate(model, xyz)
        independent = source_energy(model["interaction_inventory"], xyz)
        de = max(abs(r["energy"] - independent[r["id"]]) for r in actual["rows"])
        force = np.array(actual["forces"])
        steps = []
        for step in (1e-4, 5e-5):
            error = float(
                np.max(
                    np.abs(
                        force
                        - numerical_force(model["interaction_inventory"], xyz, step)
                    )
                )
            )
            require(
                error < 2e-5,
                "fixed_force_FD_failed",
                geometry=name,
                step=step,
                error=error,
            )
            steps.append({"step_A": step, "maximum_force_error_kcal_mol_A": error})
        torque = np.sum(np.cross(np.asarray(xyz) - np.mean(xyz, axis=0), force), axis=0)
        require(
            de < 1e-10
            and np.max(np.abs(force.sum(axis=0))) < 1e-10
            and np.max(np.abs(torque)) < 1e-10,
            "fixed_independent_or_rigid_invariance_failed",
        )
        records.append(
            {
                "geometry": name,
                "independent_row_energy_max_error": de,
                "FD": steps,
                "net_force": force.sum(axis=0).tolist(),
                "net_torque": torque.tolist(),
            }
        )
    return signed(
        {
            "schema": "island.j59-fixed-checks.v1",
            "controls": controls(),
            "source_geometry_checks": records,
            "FD_tolerance_kcal_mol_A": 2e-5,
            "complete_residual": False,
        }
    )


def analyze(model, repo, normalized):
    """Energy + Cartesian-force design with an explicit nuisance hypothesis."""
    repo, normalized = Path(repo), Path(normalized)
    groups = []
    energy_matrices = []
    force_matrices = []
    impacts = []
    tables = {}
    fixed_results = []
    all_samples = []
    for name in ("angle-angle", "bond-bond"):
        samples = json.loads((normalized / name / "samples.json").read_text())
        all_samples += samples
        tables[name] = affected(model, samples)
        vals = [conditional_columns(model, s["coordinates"]) for s in samples]
        columns = sorted(vals[0])
        groups = columns
        matrix = np.array([[r[k] for k in columns] for r in vals])
        relative = matrix - matrix[12]
        energy_matrices.append(relative)
        force = np.zeros((25, 27, len(columns)))
        for row, s in enumerate(samples):
            x = np.array(s["coordinates"])
            for i in range(9):
                for j in range(3):
                    a = x.copy()
                    b = x.copy()
                    a[i, j] += 1e-5
                    b[i, j] -= 1e-5
                    va, vb = (
                        conditional_columns(model, a),
                        conditional_columns(model, b),
                    )
                    force[row, i * 3 + j, :] = [
                        -(va[k] - vb[k]) / 2e-5 for k in columns
                    ]
        force_matrices.append(force.reshape(-1, len(columns)))
        fixed = [evaluate(model, s["coordinates"]) for s in samples]
        qe = np.array([s["energy"] for s in samples])
        fe = np.array([r["energy"] for r in fixed])
        qf = np.array([s["forces"] for s in samples])
        ff = np.array([r["forces"] for r in fixed])
        unexplained = (qe - qe[12]) - (fe - fe[12])
        fixed_results.append(
            {
                "scan": name,
                "fixed_components": [r["components"] for r in fixed],
                "fixed_energy_kcal_mol": fe.tolist(),
                "fixed_forces_kcal_mol_A": ff.tolist(),
                "QM_minus_verified_fixed_relative_energy_kcal_mol": unexplained.tolist(),
                "unexplained_energy_RMS": float(np.sqrt(np.mean(unexplained**2))),
                "unexplained_force_RMS": float(np.sqrt(np.mean((qf - ff) ** 2))),
                "meaning": "not fitting residual/error: all explicit unknown terms and discrepancy remain",
            }
        )
        for i, k in enumerate(columns):
            impacts.append(
                {
                    "scan": name,
                    "column": k,
                    "relative_column_min": float(relative[:, i].min()),
                    "relative_column_max": float(relative[:, i].max()),
                    "relative_column_RMS": float(np.sqrt(np.mean(relative[:, i] ** 2))),
                    "force_column_RMS": float(np.sqrt(np.mean(force[:, :, i] ** 2))),
                    "unit": "kcal/mol per dimensionless unit weight"
                    if "LJ_unit" in k
                    else "e^2/angstrom per effective Coulomb prefactor"
                    if "Coulomb_" in k
                    else "source equation coefficient coordinate; radian/angstrom powers as table",
                    "not_an_assigned_energy": True,
                }
            )
    e = np.vstack(energy_matrices)
    f = np.vstack(force_matrices)
    target_cols = [
        i
        for i, k in enumerate(groups)
        if k.startswith(("angle-angle:", "bond-bond_1_3:"))
    ]
    _u, sv, vh = np.linalg.svd(
        np.vstack([e, f])
        / np.maximum(np.linalg.norm(np.vstack([e, f]), axis=0), 1e-30),
        full_matrices=False,
    )
    design = np.vstack([e, f])
    scaled = design / np.maximum(np.linalg.norm(design, axis=0), 1e-30)
    gram = scaled.T @ scaled
    correlations = sorted(
        [
            {"a": groups[i], "b": groups[j], "column_cosine": float(gram[i, j])}
            for i in range(len(groups))
            for j in range(i + 1, len(groups))
        ],
        key=lambda r: -abs(r["column_cosine"]),
    )[:12]
    nulls = [
        {
            "singular_value": float(value),
            "largest_scaled_combinations": [
                {"column": groups[int(i)], "weight": float(vh[n, i])}
                for i in np.argsort(np.abs(vh[n]))[-6:][::-1]
            ],
        }
        for n, value in enumerate(sv)
        if value < sv[0] * 1e-10
    ]
    return signed(
        {
            "schema": "island.j59-residual-identifiability.v1",
            "row_metric": "energy kcal/mol and force multiplied by 1 angstrom; equal numerical row weights only for rank diagnostic, not a fitting objective",
            "next_bounded_design": {
                "deformations": [
                    {
                        "kind": "torsion",
                        "offsets_radian": [-0.1, 0.1],
                        "sites": [4, 1, 2, 3],
                    },
                    {
                        "kind": "torsion",
                        "offsets_radian": [-0.1, 0.1],
                        "sites": [1, 2, 3, 7],
                    },
                    {
                        "kind": "angle",
                        "offsets_radian": [-0.04, 0.04],
                        "sites": [1, 2, 3],
                    },
                ],
                "executed_in_J59": False,
                "independent_held_out_conformers_required": 2,
                "initial_screening_not_guaranteed_sufficient": True,
                "new_gradient_geometries": 6,
                "prerequisite": "bind numerical source rows, intramolecular pair policies and torsion geometry; establish reference protocol first",
                "selection": "screen the complete-model Jacobian before execution to reduce the identified end/middle bond-torsion correlation; independent held-outs never tune bounds/weights",
            },
            "model_identity": model["identity"],
            "affected_interactions": tables,
            "family_counts": dict(Counter(r["family"] for r in model["table"])),
            "conditional_columns": groups,
            "conditional_independent_coefficient_count": len(groups),
            "conditional_equations": "source literal linear-Phi AAT plus independent K*Phi0 product; cos Phi convention hypothesis for other torsions; pair weights grouped ONLY as lower-bound nuisance demonstration; empty-source TT represented; none bound or assigned",
            "energy_design": spectrum(e),
            "energy_force_design": spectrum(np.vstack([e, f])),
            "four_target_energy_design": spectrum(e[:, target_cols]),
            "four_target_energy_force_design": spectrum(
                np.vstack([e[:, target_cols], f[:, target_cols]])
            ),
            "source_phase1_linked_alternative": {
                "hypothesis": "interpret source AAT Phi1(0) as the existing main torsion first phase (0 radians for all six sites); separate Kphase product then absent, not a new coefficient",
                "column_count": 30,
                "energy_design": spectrum(
                    e[
                        :,
                        [
                            i
                            for i, k in enumerate(groups)
                            if not k.endswith(":Kphase_product")
                        ],
                    ]
                ),
                "energy_force_design": spectrum(
                    np.vstack([e, f])[
                        :,
                        [
                            i
                            for i, k in enumerate(groups)
                            if not k.endswith(":Kphase_product")
                        ],
                    ]
                ),
                "physical_identifiability": False,
            },
            "null_combinations": nulls,
            "largest_conditional_column_correlations": correlations,
            "coefficient_uncertainty": "not estimated; nuisance conventions unbound and reference adequacy not established",
            "impacts": impacts,
            "fixed_subtraction": fixed_results,
            "unknown_numerical_values": None,
            "physical_identifiability": False,
            "fit_completed": False,
            "fit_validated": False,
            "fitted_candidates": [],
            "decision": "no_fit",
            "quantitative_reasons": [
                "The four-column full-rank target subproblem excludes the other 27 explicit nuisance columns; use the augmented rank/conditioning, not its partial rank.",
                "Conditional pair grouping, torsion geometry, linear-Phi phase mapping, empty torsion-torsion applicability remain unbound; even a full-rank conditional matrix would not close the physical model.",
                "The two J58 5x5 grids share one conformer; neither column split nor corners are independent held-out conformers.",
            ],
            "smallest_next_information": [
                "Numerical author/FRC rows (or independently constrained QM joint targets) for end/middle bond-torsion, angle-torsion, angle-angle-torsion in hsi-sio-osi-sio, including equilibrium/phase and omission authority for H-Si-O-Si-H torsion-torsion.",
                "Explicit intramolecular pair weights/exclusions for path lengths 1..4, Coulomb prefactor/dielectric and finite isolated cutoff policy; these cannot be estimated as the target Kaa/Kbb13.",
                "After the model is bound, at least independent Si-O-Si torsion and center-angle deformations, plus a distinct held-out conformer; use rank-revealing design selection on the complete model before launching scans.",
            ],
            "scope": "9-atom disiloxane Si-H control only",
            "methyl_transfer_required": {
                "HMDS": [
                    "sio with osi,c3,c3,c3; c3 with sio,hc,hc,hc",
                    "Si-C and C-H source-equilibrium/charge dependencies",
                    "c3-sio-c3-hc, c3-sio-osi-sio, osi-sio-c3-hc couplings",
                ],
                "DP2": [
                    "interior sio with osi,osi,c3,c3; distinct terminal/central Si",
                    "osi-sio-osi-sio and sio-osi-sio-osi couplings",
                    "independent methyl rotamers and chain conformers with reference protocol/held-out validation",
                ],
            },
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
