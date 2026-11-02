"""Targeted HMDS intake and protected-subset conic research analysis.

Historical adapters and evidence are immutable. This module extends their
constraint assembly, never their force-field or numerical acceptance limits.
"""

from itertools import combinations

import numpy as np
from scipy import sparse

from island.pcff_completion import digest, require

from .hmds_constrained import canonical, verify
from .hmds_paths import equivalent, leakage

SELECTED = ("joint-032", "joint-002", "joint-024", "joint-034")
PARENT_DECLARATIONS = {
    "schema": "island.j66-receipt.v1",
    "training_feasible": True,
    "candidate_qualified": False,
    "parameter_completeness": "incomplete",
    "explicit_research_trial_vector_complete": True,
    "scientific_validation": "not_established",
    "chemical_transfer": "not_established",
    "cross_engine_equivalence": "not_claimed",
    "automatic_assignment": False,
    "new_QM": 0,
    "validation_QM": 0,
    "MD_steps": 0,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "original_author_artifact_status": "not_yet_supplied",
}
PARENT_RECEIPT_FIELDS = set(PARENT_DECLARATIONS) | {
    "J54_J65_regressions",
    "J65_prerequisite",
    "base_commit",
    "branch",
    "fit_identity",
    "focused_tests",
    "historical_files_preserved",
    "identity",
    "manifest_sha256",
    "numerical_reconstruction",
    "software_counts",
    "uncertainty_identity",
}


def identity(record):
    require(
        record.get("identity")
        == digest({k: v for k, v in record.items() if k != "identity"}),
        "J67_record_identity",
    )


def parent_semantics(receipt, manifest, decision):
    """Scientific declarations remain checked even with valid recomputed seals."""
    for r in (receipt, manifest, decision):
        identity(r)
    require(set(receipt) == PARENT_RECEIPT_FIELDS, "J67_parent_receipt_schema")
    for key, value in PARENT_DECLARATIONS.items():
        require(
            type(receipt.get(key)) is type(value) and receipt[key] == value,
            "J67_parent_semantic_declaration",
            key=key,
        )
    require(
        manifest.get("schema") == "island.numerical-evidence-manifest.v1"
        and set(manifest)
        == {
            "schema",
            "base_commit",
            "historical_files",
            "implementation_files",
            "artifacts",
            "identity",
        }
        and all(
            isinstance(manifest[k], dict) and manifest[k]
            for k in ("historical_files", "implementation_files", "artifacts")
        ),
        "J67_parent_manifest_schema",
    )
    require(
        receipt["base_commit"]
        == manifest["base_commit"]
        == "712e588a8ca9b3c0d3ba03114a51064d3a2a83bb"
        and receipt["J65_prerequisite"] == "1b03f835718599b7f6e53e3331400518112882ec",
        "J67_parent_ancestors",
    )
    expected = {
        "training_feasible": True,
        "constrained_trial_verified": True,
        "nominal_bound_active": True,
        "uncertainty_qualified": False,
        "accepted_candidates": [],
        "new_QM": 0,
        "validation_QM": 0,
        "MD_steps": 0,
        "automatic_assignment": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    require(set(decision) == set(expected) | {"identity"}, "J67_parent_decision_schema")
    for key, value in expected.items():
        require(
            type(decision.get(key)) is type(value) and decision[key] == value,
            "J67_parent_decision_semantics",
            key=key,
        )
    require(
        receipt["numerical_reconstruction"]
        == {
            "MD_steps": 0,
            "atoms": 27,
            "bootstrap_draws": 64,
            "configurations": 22,
            "eligible": 2,
            "native_reconstructed": True,
            "new_QM": 0,
            "offline": True,
            "physical_columns": 58,
            "training_feasible": True,
            "validation_untouched": True,
        },
        "J67_parent_reconstruction_semantics",
    )
    return True


def exact_intake(model, rows, expected, old, reserved):
    require(
        tuple(r["id"] for r in rows) == SELECTED and rows == expected,
        "J67_exact_selected_geometry",
    )
    xyz = [r["coordinates"] for r in rows]
    for a, b in combinations(xyz, 2):
        require(
            not equivalent(model, a, b, 1e-6), "J67_duplicate_or_symmetry_equivalent"
        )
    leakage(model, [r["coordinates"] for r in old] + list(reserved), xyz)
    return True


def bootstrap_draws(observations, seed=4067, count=64):
    blocks = sorted({r["block"] for r in observations})
    require(len(blocks) == 14, "J67_four_new_blocks_required")
    rng = np.random.default_rng(seed)
    return [
        {
            "replicate": k,
            "sampled_blocks": rng.choice(blocks, len(blocks), replace=True).tolist(),
        }
        for k in range(count)
    ]


def protected_system(problem, lower, upper, kind="rho", protected=None, interior=0):
    """Append protected cones after conversion to the one nominal variable scale."""
    system = canonical(problem, lower, upper, kind=kind, interior=interior)
    if protected is None:
        return system
    extra = canonical(protected, lower, upper, kind=kind, bounded=False)
    # The optional second objective epigraph is excluded. Only the full-data
    # weighted LS objective selects a vector; the original subset adds constraints.
    stop = len(extra["b"])
    soc = list(extra["soc"])
    if kind == "least_squares":
        stop -= soc.pop()
    transform = sparse.diags(np.r_[extra["scale"] / system["scale"], 1.0])
    ea = (extra["A"] @ transform)[:stop]
    l, m = system["linear"], extra["linear"]
    system["A"] = sparse.vstack(
        [system["A"][:l], ea[:m], system["A"][l:], ea[m:]], format="csc"
    )
    system["b"] = np.r_[
        system["b"][:l], extra["b"][:m], system["b"][l:], extra["b"][m:stop]
    ]
    system["linear"] += m
    system["soc"] += soc
    system["protected"] = True
    system["labels"] = {
        "full": system["labels"],
        "protected": extra["labels"],
        "assembly": "full linear,protected linear,full SOC,protected SOC; no second objective",
    }
    return system


def verified_subsets(problem, protected, values, lo, hi, rho=1):
    checks = {"full": verify(problem, values, lo, hi, rho)}
    if protected is not None:
        checks["protected22"] = verify(protected, values, lo, hi, rho)
    for name, check in checks.items():
        selected_problem = problem if name == "full" else protected
        check["original_limit_slacks"] = verify(
            selected_problem, values, lo, hi, rho=1
        )["slacks"]
    return checks


def physically_feasible(checks):
    return all(r["original_limits_satisfied"] for r in checks.values())


def solver_qualified(result):
    """Inherited accurate status plus explicit residual/gap checks; no waiver."""
    return (
        result["status"] == "Solved"
        and max(result["primal_residual"], result["dual_residual"]) <= 1e-8
        and result["absolute_gap"] <= 1e-7
    )


def candidate_gates(trial, uncertainty, design):
    checks = trial["checks"]
    return {
        "training": physically_feasible(checks),
        "numerical": solver_qualified(trial["solution"]),
        "bounds": not (
            checks["full"]["active_bound_indices"]
            or trial.get("raw_active_original_bounds", [])
        ),
        "identifiability": design["rank"] == 59 and design["scaled_condition"] <= 10000,
        "uncertainty": uncertainty["qualified"],
    }


def validation_activation(gates, previous_consumption=False):
    require(not previous_consumption, "J67_validation_already_consumed")
    return all(gates.values())


def exact_request(actual, expected):
    """A newly signed producer cannot authorize altered graph/protocol/units."""
    require(actual == expected, "J67_exact_request_semantics")
    require(
        actual["atom_order"] == list(range(1, 28))
        and actual["charge"] == 0
        and actual["multiplicity"] == 1
        and np.asarray(actual["coordinates"]).shape == (27, 3),
        "J67_HMDS_request_scope",
    )
    return True


def preserved_producers(rows, original):
    require(
        rows[: len(original)] == original, "J67_original_observation_producer_mutation"
    )
    return True


def proposal_sensitivities(model, proposal, diagnosis):
    """Rebuild retained weak-direction sensitivities at exact stored coordinates."""
    from .hmds_baseline import extended_parts

    directions = np.array(
        [r["physical_direction"] for r in diagnosis["nominal_weak"]]
        + [r["weak"][-1]["physical_direction"] for r in diagnosis["weak_draws"]]
    )
    require(
        np.array_equal(directions, proposal["weak_directions"]),
        "J67_proposed_weak_direction_mapping",
    )
    results = []
    for row in proposal["selected"]:
        part = extended_parts(model, row["coordinates"])
        e = np.r_[part["energy_jacobian"], 1.0] / (0.1 * np.sqrt(22))
        f = np.c_[part["force_jacobian"].reshape(81, 58), np.zeros(81)] / (
            0.5 * np.sqrt(81 * 22)
        )
        measured = np.linalg.norm(np.vstack([e, f]) @ directions.T, axis=0)
        error = float(np.max(abs(measured - row["target_weak_direction_norms"])))
        require(error < 1e-12, "J67_proposal_sensitivity_reconstruction")
        results.append(
            {
                "id": row["id"],
                "geometry_identity": row["geometry_identity"],
                "maximum_absolute_sensitivity_error": error,
                "weak_direction_norms": measured.tolist(),
            }
        )
    return results


def solver_state_diagnostics(system, result):
    """Recompute raw conic state closure independently of solver status text."""
    x, dual, slack = (np.array(result[k], float) for k in ("x", "dual", "slack"))
    a, b, p, q = (system[k] for k in ("A", "b", "P", "q"))
    require(
        x.shape == q.shape and dual.shape == slack.shape == b.shape,
        "J67_solver_state_dimensions",
    )
    require(
        np.allclose(x[:-1] / system["scale"], result["values"], atol=1e-12, rtol=0),
        "J67_solver_physical_mapping",
    )
    primal = a @ x + slack - b
    stationarity = p @ x + q + a.T @ dual
    cone_margins = {}
    for name, values in [("primal", slack), ("dual", dual)]:
        margins = [float(np.min(values[: system["linear"]]))]
        offset = system["linear"]
        for size in system["soc"]:
            margins.append(
                float(
                    values[offset] - np.linalg.norm(values[offset + 1 : offset + size])
                )
            )
            offset += size
        cone_margins[name] = min(margins)
    po = float(0.5 * x @ p @ x + q @ x)
    do = float(-0.5 * x @ p @ x - b @ dual)
    return {
        "primal_maximum_absolute_closure": float(np.max(abs(primal))),
        "dual_maximum_absolute_stationarity": float(np.max(abs(stationarity))),
        "primal_normalized_L2": float(
            np.linalg.norm(primal)
            / max(1, np.linalg.norm(b), np.linalg.norm(a @ x), np.linalg.norm(slack))
        ),
        "dual_normalized_L2": float(
            np.linalg.norm(stationarity)
            / max(1, np.linalg.norm(q), np.linalg.norm(a.T @ dual))
        ),
        "minimum_cone_margins": cone_margins,
        "recomputed_primal_objective": po,
        "recomputed_dual_objective": do,
        "recomputed_absolute_gap": abs(po - do),
        "qualification_scope": "Independent raw-state diagnostics; original-unit physical checks and frozen solver status/residual/gap gates remain separate.",
    }
