"""Freeze the J66 numerical contract; no QM/MD execution entry point."""

import argparse
import importlib.metadata
import platform
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_precision import read

BASE = "712e588a8ca9b3c0d3ba03114a51064d3a2a83bb"
PREREQUISITE = "1b03f835718599b7f6e53e3331400518112882ec"
BRANCH = "codex/phase-4j66-hmds-constrained-fit-feasibility"
PARENTS = (
    "plan/scientific-plan.json",
    "plan/analysis-source-bindings.json",
    "fitting-dataset.json",
    "actual-fitting-system.npz",
    "trial-fit.json",
    "block-uncertainty.json",
    "receipt.json",
    "native-outputs.tar.xz",
)


def freeze(repo, root):
    require(
        not (root / "plan/scientific-plan.json").exists(), "J66_plan_already_frozen"
    )
    old = repo / "docs/evidence/phase_4j65"
    plan = read(old / "plan/scientific-plan.json")
    result = signed(
        {
            "schema": "island.hmds-constrained-feasibility.v1",
            "base_commit": BASE,
            "J65_prerequisite": PREREQUISITE,
            "J65_in_fetched_main": True,
            "branch": BRANCH,
            "parent_files": {name: pin(old / name) for name in PARENTS},
            "graph_identity": plan["graph_identity"],
            "atom_order": plan["atom_order"],
            "configuration_order": [
                r["id"] for r in read(old / "fitting-dataset.json")["observations"]
            ],
            "protocol": plan["protocol"],
            "parameter_specifications": plan["parameter_specifications"],
            "inherited_fit": plan["fit"],
            "reconstruction_tolerance": {
                "metrics_absolute": 1e-9,
                "physical_arrays_absolute": 1e-10,
                "native_replay": "exact original J65 replay in a relocated directory",
            },
            "problem": "58 unchanged physical variables plus one global energy-only intercept; r=MM-QM. Minimize rho>=0 under original bounds, ||rE||<=.1sqrt(N)rho, |rEi|<=.3rho, ||rF||<=.5sqrt(3NA)rho, ||rFij||<=2rho.",
            "conditional_fit": "If independently verified feasible: minimize original J65 mean(rE^2)/.1^2+mean(rF^2)/.5^2 subject to rho=1 limits and original bounds. Strict convex objective/full rank provides deterministic unique selection; no regularization. Do not tighten/widen physical gates for numerical tolerance.",
            "scaling": "Column norms of actual J65 weighted matrix including offset; z=norm*p. Each residual cone divided by its frozen physical threshold; bound inequalities divided by original bound width. Recompute scaling for each exact bootstrap sample.",
            "solvers": {
                "primary": "clarabel direct, max_iter300,time_limit120s,tol_gap_abs/rel1e-9,tol_feas1e-10",
                "cross_check": "SCS independent backend, max_iters100000,time_limit120s,eps_abs/rel1e-7; same declared physical problem",
                "diagnostic": "If needed: repeat rho with bounds removed, and force-only rho with original finite bounds for a support-function dual lower bound. No unrestricted result is a candidate. Original-bound interior margin1e-6 tested separately when feasible.",
            },
            "uncertainty": "Exact 64 retained J65 sampled_blocks, no resampling/replacement. For each draw solve rho then conditional LS if feasible. All four constraints apply to sampled training rows with multiplicities. Evaluate original22 predictions separately. Eligibility additionally rank59,condition<=10000,no bound activity,solver success,physical verification. Inherit80% and prediction spread limits. Survivor-only intervals explicitly conditional; retain all failures.",
            "diagnosis": "SVD of actual weighted scaled matrices; omitted-block/failure overlaps, bound recurrence, physical weak-direction and predicted force sensitivity. No validation observations.",
            "geometry_only": "If robustness fails, rank a bounded existing J63 graph-controlled candidate pool against measured weak directions; exclude all existing training and reserved validation by graph-aware proper alignment. At most4 proposals; no execution.",
            "resources": {
                "new_QM": 0,
                "MD_steps": 0,
                "workers": 1,
                "BLAS_threads": 1,
                "memory_GiB": 4,
                "per_solve_seconds": 120,
                "campaign_seconds": 21600,
                "scratch": "external temporary reconstruction only; no native engine invoked",
            },
            "failure_handling": "Preserve each solve/status/state and exception; no inaccurate/incomplete status is infeasibility. Independent original-unit witness can establish feasibility. Numerical lower bounds must include checked dual-cone/support residual error. No replacement bootstrap draws.",
            "environment": {
                "python": platform.python_version(),
                "packages": {
                    k: importlib.metadata.version(k)
                    for k in ("numpy", "scipy", "clarabel", "scs")
                },
            },
            "reserved_validation": "all3 remain untouched in this phase regardless of outcome; zero new QM",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "accepted_candidates": [],
            "automatic_assignment": False,
            "author_confirmations": "28 user-attested; original_author_artifact_status not_yet_supplied",
        }
    )
    write(root / "plan/scientific-plan.json", result)
    return result["identity"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j66"))
    a = p.parse_args()
    print(freeze(a.repo.resolve(), a.evidence.resolve()))
