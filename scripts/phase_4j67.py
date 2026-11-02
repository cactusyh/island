"""Freeze four exact J66 proposals, then use the existing counted native capture."""

import argparse
import json
import shutil
import time
from pathlib import Path

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.hmds_augmentation import (
    bootstrap_draws,
    exact_intake,
    parent_semantics,
)
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import read
from scripts.phase_4j64 import run_job

BASE = "fbf093d7a75262a652192a57981cee0872c7c372"
PREREQUISITE = "6e1e3ae7dd20e4d153c6e23dc9107379a475350f"
BRANCH = "codex/phase-4j67-hmds-targeted-information-augmentation"
EXECUTION = (
    "scripts/phase_4j64.py",
    "src/island/pcff_qm_pilot/hmds_engine.py",
    "src/island/pcff_qm_pilot/hmds_precision.py",
    "src/island/pcff_qm_pilot/control_engine.py",
)


def resources(work):
    cg = Path("/sys/fs/cgroup")
    stat = {
        k: int(v)
        for k, v in (s.split() for s in (cg / "memory.stat").read_text().splitlines())
    }
    maximum = int((cg / "memory.max").read_text())
    current = int((cg / "memory.current").read_text())
    result = {
        "memory_max": maximum,
        "memory_current": current,
        "memory_stat": stat,
        "cpu_max": (cg / "cpu.max").read_text().strip(),
        "disk_free": shutil.disk_usage(work.parent).free,
    }
    # Retained TZVP used4GiB configured; allow6GiB RSS and2GiB margin.
    # File cache is reclaimable, anonymous and kernel memory are not assumed free.
    require(
        maximum
        - stat["anon"]
        - sum(stat.get(k, 0) for k in ("kernel_stack", "pagetables", "sock", "slab"))
        > 8 * 1024**3
        and result["disk_free"] > 64 * 1024**3,
        "J67_native_resource_preflight",
    )
    quota, period = result["cpu_max"].split()
    require(quota == "max" or int(quota) / int(period) >= 4, "J67_CPU_preflight")
    return result


def prepare(repo, root, work):
    require(
        not (root / "plan/scientific-plan.json").exists(), "J67_plan_already_frozen"
    )
    parent = repo / "docs/evidence/phase_4j66"
    parent_semantics(
        read(parent / "receipt.json"),
        read(parent / "manifest.json"),
        read(parent / "decision.json"),
    )
    require(
        "training_feasible" in (root / "validation/J66-reconstruction.log").read_text(),
        "J67_J66_replay_required",
    )
    p = read(repo / "docs/evidence/phase_4j65/plan/scientific-plan.json")
    old = read(repo / "docs/evidence/phase_4j65/fitting-dataset.json")["observations"]
    proposed = read(parent / "geometry-only-next-design.json")["selected"]
    model = load_model(repo)
    exact_intake(model, proposed, proposed, old, p["validation_coordinates"])
    target = dict(p["protocol"])
    # Explicitly freeze effective native settings, avoiding stale nested defaults.
    target["settings"] = request(p, proposed[0]["coordinates"])["settings"]
    rows = [{**r, "block": r["id"]} for r in proposed]
    draws = bootstrap_draws(old + rows)
    oldblocks = sorted({r["block"] for r in old})
    require(
        all(any(r["id"] not in d["sampled_blocks"] for d in draws) for r in rows),
        "J67_new_blocks_not_omittable",
    )
    work.mkdir(parents=True, exist_ok=True)
    preflight = resources(work)
    timings = read(repo / "docs/evidence/phase_4j65/compute-accounting.json")[
        "elapsed_seconds"
    ]
    write(
        root / "parent-semantic-audit.json",
        signed(
            {
                "schemas_and_scientific_declarations_checked": True,
                "J66_receipt_sha256": pin(parent / "receipt.json"),
                "J66_native_and64_bootstrap_replayed": True,
                "historical_mutation_tests_claimed": False,
                "graph_identity": model["graph_identity"],
            }
        ),
    )
    plan = signed(
        {
            "schema": "island.hmds-targeted-augmentation.v1",
            "frozen_unix": time.time(),
            "base_commit": BASE,
            "required_J66": PREREQUISITE,
            "J66_in_fetched_main": True,
            "branch": BRANCH,
            "required_ancestors": {
                "J65": "1b03f835718599b7f6e53e3331400518112882ec",
                "J66": PREREQUISITE,
            },
            "parent_receipt_sha256": pin(parent / "receipt.json"),
            "proposal_sha256": pin(parent / "geometry-only-next-design.json"),
            "graph_identity": model["graph_identity"],
            "elements": p["elements"],
            "atom_order": list(range(1, 28)),
            "charge": 0,
            "multiplicity": 1,
            "protocol": target,
            "training_requests": rows,
            "validation_coordinates": p["validation_coordinates"],
            "parameter_specifications": p["parameter_specifications"],
            "inherited_fit": p["fit"],
            "original22_ids": [r["id"] for r in old],
            "configuration_order": [r["id"] for r in old + rows],
            "hypotheses": [
                "Four selected observations may constrain recurring bound/weak directions.",
                "Better nominal matrix conditioning need not improve omission robustness.",
                "New local residuals may expose incompatibility of unchanged model/data/bounds.",
            ],
            "geometry_scope": "Exact center deformations; not independently sampled conformational basins.",
            "intake": "Exact request/protocol/geometry/hash; native convergence,399 orbital functions, all27 gradients; independent decimal energy and ASCII force conversion within1e-8 energy/1e-9 force. Existing J65 protocol precision remains mandatory, unchanged.",
            "limits": {
                "maximum_new_gradients": 8,
                "maximum_failed_attempt_retries": 1,
                "campaign_wall_seconds": 18000,
                "workers": 1,
                "threads": 4,
                "memory_GiB": 4,
                "per_job_wall_seconds": 1800,
                "scratch_GiB": 64,
            },
            "resources": preflight,
            "retained_TZVP_timings_seconds": {
                k: v for k, v in timings.items() if k.startswith("training/")
            },
            "resource_justification": "Retained dense-TZVP training470-568s;1800s timeout,4GiB configured,6GiB RSS envelope+2GiB headroom; one native worker,4threads. File cache separately reported. No QZVP.",
            "execution_sources": {n: pin(repo / n) for n in EXECUTION},
            "objective": "Full26 mean(rE^2)/.1^2 + mean(rF^2)/.5^2; one unbounded energy-only intercept; no regularization/duplicate curvature.",
            "constraints": "Original bounds and four limits simultaneously on protected22 and full26: E_RMSE<=.1,F_RMS<=.5,max|E|<=.3,max atom-vector norm<=2. Rho scales all limits only for feasibility, never acceptance.",
            "conic_scaling": "Actual weighted full59column norms. Protected cones transformed to same variables. Residual cones divided by physical limits; bound rows by widths. Bootstrap recomputes sample scaling.",
            "solver": {
                "backend": "clarabel",
                "version": "0.11.1",
                "settings": {
                    "max_iter": 300,
                    "time_limit": 120.0,
                    "tol_gap_abs": 1e-9,
                    "tol_gap_rel": 1e-9,
                    "tol_feas": 1e-10,
                },
                "qualification": "Inherited Solved required, additionally primal/dual residual<=1e-8 and absolute gap<=1e-7. AlmostSolved unqualified; report physical feasibility independently. No policy relaxation.",
                "crosscheck": "Independent SCS rho if nominal not physically feasible; checked dual support bound if rho>1; no solver failure is infeasibility.",
            },
            "selection": {
                "rule": "Raw solution; if LS physically violates, first strict-feasible convex blend toward verified rho witness in frozen sequence. Does not alter raw solver classification or erase original-bound activity.",
                "alphas": [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
            },
            "uncertainty": {
                "seed": 4067,
                "replicates": 64,
                "old_blocks": oldblocks,
                "new_blocks": [r["id"] for r in rows],
                "rule": "14 blocks sampled with replacement14times. Optimize and constrain ONLY sampled observations with multiplicity; no protected/full post-fit constraints on omitted observations. Rho then conditional LS. Retain every failure; no replacement. Evaluate22/26 after fitting.",
                "rank": 59,
                "condition_max": 10000,
                "bound_activity_fraction": 1e-6,
                "minimum_eligible_fraction": 0.8,
                "prediction_energy_RMS95_max": 0.1,
                "prediction_force_RMS95_max": 0.5,
                "intervals": "Absolute2.5/97.5percentiles, correlations and atom predictive effects; survivor-conditional/unqualified if overall failure.",
            },
            "validation_activation": "All nominal physical/numerical/original-bound/rank/condition and uncertainty gates pass; frozen vector/intercept hash, graph leakage rechecked, one-time activation and remaining budget>=3. No activation on any failure.",
            "failure_handling": "Count startup/interrupted attempts; archive commands/codes/stdout/stderr/native before cleanup. At most one additive retry, never overwrite failed directories; resume exact hash-verified completed requests only. Stop if required stage incomplete.",
            "stopping": "After four observations and numerical assessment stop this allocation; no automatic further QM or model expansion.",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "automatic_assignment": False,
            "original_author_artifact_status": "not_yet_supplied",
            "author_confirmations": "28 user-attested",
            "MD_steps": 0,
        }
    )
    write(root / "plan/scientific-plan.json", plan)
    write(
        root / "plan/expanded-bootstrap-draws.json",
        signed({"seed": 4067, "draws": draws, "replacement_resamples": 0}),
    )
    write(work / "scientific-plan.json", plan)
    jobs = {}
    for row in rows:
        key = "training/" + row["id"]
        r = request(plan, row["coordinates"])
        write(work / "requests" / (key + ".json"), r)
        write(root / "plan/requests" / (key + ".json"), r)
        jobs[key] = {
            "stage": "training",
            "request_identity": digest(r),
            "memory_GiB": 4,
            "wall_seconds": 1800,
        }
    write(work / "jobs.json", {"jobs": jobs})
    write(root / "plan/jobs.json", {"jobs": jobs})
    write(
        root / "plan/frozen-file-pins.json",
        signed(
            {
                "files": {
                    str(f.relative_to(root / "plan")): pin(f)
                    for f in sorted((root / "plan").rglob("*.json"))
                }
            }
        ),
    )
    print(
        json.dumps(
            {"plan": plan["identity"], "preflight": "passed", "requests": list(jobs)}
        )
    )


def execute(repo, root, work, python):
    for key in read(work / "jobs.json")["jobs"]:
        write(
            work / "resources" / ("preflight-" + key.replace("/", "-") + ".json"),
            resources(work),
        )
        code = run_job(repo, root, work, python, key)
        if code:
            return code
    return 0


def conditional_validation(repo, root, work, python):
    """One-time execution, solely after all preregistered gates pass."""
    import numpy as np

    from island.pcff_qm_pilot.hmds_augmentation import identity, validation_activation
    from island.pcff_qm_pilot.hmds_baseline import evaluate
    from island.pcff_qm_pilot.hmds_paths import leakage
    from scripts.analyze_phase_4j63 import sample
    from scripts.analyze_phase_4j65 import independent_value, relative
    from scripts.analyze_phase_4j67 import retain

    activation = read(root / "validation-activation.json")
    identity(activation)
    if not activation["activated"]:
        require(
            not any(
                k.startswith("validation/")
                for k in read(work / "attempts.json")["calls"]
            ),
            "J67_unauthorized_validation_attempt",
        )
        return 0
    require(validation_activation(activation["gates"]), "J67_validation_prerequisite")
    require(
        not (root / "validation-consumption.json").exists(),
        "J67_validation_already_consumed",
    )
    plan = read(root / "plan/scientific-plan.json")
    trial = read(root / "constrained-trial.json")
    values = trial["selected_values"]
    require(
        digest(values) == activation["candidate_values_sha256"]
        and pin(root / "constrained-trial.json") == activation["trial_sha256"]
        and pin(root / "fitting-dataset.json") == activation["dataset_sha256"],
        "J67_candidate_changed_after_activation",
    )
    attempts = read(work / "attempts.json")["calls"]
    require(len(attempts) + 3 <= 8, "J67_validation_attempt_budget")
    observations = read(root / "fitting-dataset.json")["observations"]
    model = load_model(repo)
    leakage_check = leakage(
        model, [r["coordinates"] for r in observations], plan["validation_coordinates"]
    )
    write(
        root / "candidate-freeze.json",
        signed(
            {
                "values": values,
                "activation_identity": activation["identity"],
                "trial_sha256": pin(root / "constrained-trial.json"),
                "plan_sha256": pin(root / "plan/scientific-plan.json"),
                "leakage_check": leakage_check,
                "offset_refitting": False,
            }
        ),
    )
    # Consumption is durable BEFORE even a startup attempt. A failed attempt
    # cannot leave the data labelled untouched.
    write(
        root / "validation-consumption.json",
        signed(
            {
                "consumed": True,
                "activation_identity": activation["identity"],
                "frozen_candidate_sha256": pin(root / "candidate-freeze.json"),
            }
        ),
    )
    jobs = read(work / "jobs.json")
    rows = []
    for number, xyz in enumerate(plan["validation_coordinates"]):
        key = f"validation/reserved-{number}"
        r = request(plan, xyz)
        write(work / "requests" / (key + ".json"), r)
        jobs["jobs"][key] = {
            "stage": "validation",
            "request_identity": digest(r),
            "memory_GiB": 4,
            "wall_seconds": 1800,
        }
        write(work / "jobs.json", jobs)
        write(
            work / "resources" / ("preflight-validation-" + str(number) + ".json"),
            resources(work),
        )
        code = run_job(repo, root, work, python, key)
        if code:
            write(
                root / "validation-result.json",
                signed(
                    {
                        "completed": False,
                        "consumed": True,
                        "failed_key": key,
                        "returncode": code,
                        "candidate_qualified": False,
                    }
                ),
            )
            retain(root, work)
            return code
        native = sample(work, key)
        independent_value(work, key, native)
        mm = evaluate(
            model,
            xyz,
            dict(
                zip(
                    [r["name"] for r in plan["parameter_specifications"]],
                    values[:-1],
                    strict=True,
                )
            ),
        )
        from island.pcff_qm_pilot.hmds_precision import decimal_energy

        energy = relative(
            str(decimal_energy(work / "outputs" / key)[0]),
            observations[0]["native_center_energy_hartree_decimal"],
        )
        de = mm["energy"] + values[-1] - energy
        df = mm["forces"] - np.array(native["forces"])
        rows.append(
            {
                "id": key,
                "energy_error_kcal_mol": de,
                "Cartesian_force_errors": df.tolist(),
                "atom_vector_errors": np.linalg.norm(df, axis=1).tolist(),
            }
        )
    de = np.array([r["energy_error_kcal_mol"] for r in rows])
    df = np.array([r["Cartesian_force_errors"] for r in rows])
    passed = bool(
        np.sqrt(np.mean(de**2)) <= 0.1
        and np.sqrt(np.mean(df**2)) <= 0.5
        and abs(de).max() <= 0.3
        and np.linalg.norm(df, axis=2).max() <= 2
    )
    write(
        root / "validation-result.json",
        signed(
            {
                "completed": True,
                "consumed": True,
                "rows": rows,
                "candidate_qualified": passed,
                "parameters_and_intercept_refitted": False,
            }
        ),
    )
    decision = read(root / "decision.json")
    decision.pop("identity", None)
    decision.update(
        candidate_qualified=passed,
        accepted_coefficients=58 if passed else 0,
        validation_QM=3,
        new_QM=7,
    )
    write(root / "decision.json", signed(decision))
    retain(root, work)
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=["prepare", "execute", "validate"])
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j67"))
    p.add_argument(
        "--work", type=Path, default=Path("../island-validation/phase4j67-qm")
    )
    p.add_argument(
        "--python",
        type=Path,
        default=Path("../island-validation/phase4j36-native-env/bin/python"),
    )
    a = p.parse_args()
    if a.operation == "prepare":
        prepare(a.repo.resolve(), a.evidence.resolve(), a.work.resolve())
    elif a.operation == "execute":
        raise SystemExit(
            execute(
                a.repo.resolve(),
                a.evidence.resolve(),
                a.work.resolve(),
                a.python.resolve(),
            )
        )

    if a.operation == "validate":
        raise SystemExit(
            conditional_validation(
                a.repo.resolve(),
                a.evidence.resolve(),
                a.work.resolve(),
                a.python.resolve(),
            )
        )
