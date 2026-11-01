"""Freeze and execute a protocol-specific HMDS experiment with existing capture."""

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.campaign import verify_completed
from island.pcff_qm_pilot.hmds_baseline import NAMES, extended_parts
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_hamiltonian import UNITS, load_model
from island.pcff_qm_pilot.hmds_paths import equivalent, leakage, path
from island.pcff_qm_pilot.hmds_precision import protocol_contrast, read
from scripts.curvature_phase_4j63 import unknown_null_check
from scripts.phase_4j64 import run_job

BASE = "76db72c0ab4c83f6af711d897b33761efbe12e6a"
PREREQUISITE = "71777122f12179f6c5ce2ca6828d55d92a2fe635"
BRANCH = "codex/phase-4j65-hmds-protocol-specific-joint-fit"


def prepare(repo, root, work):
    require(not (root / "plan/scientific-plan.json").exists(), "frozen_plan_exists")
    old = repo / "docs/evidence/phase_4j64"
    original = repo / "docs/evidence/phase_4j63"
    require(
        read(old / "receipt.json")["reconstruction"]["native_gradients"] == 17,
        "J64_reconstruction_required",
    )
    model = load_model(repo)
    p = read(original / "plan/scientific-plan.json")
    x = np.array(read(original / "stationary-center.json")["center"]["coordinates"])
    target = {**p["protocol"], "grid": [250, 1454], "thresholds": [1e-12, 1e-10]}
    contract = {**p, "protocol": target}
    names = extended_parts(model, x)["columns"]
    specifications = []
    # Broad local research bounds, not source values or a stability certificate.
    limits = {
        "angle-angle": 300.0,
        "bond-bond_1_3": 400.0,
        "angle-angle-torsion_1": 300.0,
        "end_bond-torsion_3": 100.0,
        "middle_bond-torsion_3": 100.0,
        "angle-torsion_3": 100.0,
    }
    for name in names:
        family = name.split(":")[0]
        if name in NAMES:
            types = ["c3", "sio"] if name == NAMES[0] else ["c3", "hc"]
            row = next(
                r
                for r in model["interaction_inventory"]
                if r["family"] == "quartic_bond"
                and sorted(model["typing"][str(i)] for i in r["sites"]) == types
            )
            k2 = row["assignment"]["values"][1] / 4.184
            low, high, unit = -0.5 * k2, k2, "kcal/mol/angstrom^2"
            reason = "Keep total quadratic bond stiffness between0.5 and2 times source K2; retain source r0/K3/K4."
        else:
            low, high, unit = -limits[family], limits[family], UNITS[family]
            reason = "Generous local research envelope: mixed quadratic caps comparable to2*sqrt(source diagonal stiffnesses); torsion-linear amplitudes capped at100 so0.04A/0.06rad changes contribute at most4/6kcal per harmonic/site. Multiplicity and global stability are not certified; bound activity rejects eligibility."
        specifications.append(
            {
                "name": name,
                "units": unit,
                "lower": low,
                "upper": high,
                "justification": reason,
            }
        )
    primary = [{"id": "center", "coordinates": x.tolist(), "block": "center"}]
    refs = read(original / "normalized-reference.json")["samples"]
    for key, value in sorted(refs.items()):
        if key.startswith("paths/"):
            primary.append(
                {
                    "id": key,
                    "coordinates": value["coordinates"],
                    "block": key.split("/")[1].split("-h")[0],
                }
            )
    joints = [
        {**r, "block": r["id"]}
        for r in read(old / "updated-joint-design.json")["selected_geometry_only"]
    ]
    validation = read(original / "plan/reserved-validation.json")["coordinates"]
    # Reuse exact requests only; protocol, graph and bitwise geometry are checked.
    archive = work / "parents/J64"
    unpack(old / "native-outputs.tar.xz", archive)
    reusable = {}
    for row in primary:
        wanted = request(contract, row["coordinates"])
        matches = []
        for file in sorted((archive / "requests").rglob("*.json")):
            r = read(file)
            if r == wanted:
                key = str(file.relative_to(archive / "requests"))[:-5]
                if (archive / "outputs" / key / "result.json").exists():
                    matches.append(key)
        if matches:
            require(len(matches) == 1, "ambiguous_exact_reuse")
            key = matches[0]
            reusable[row["id"]] = {
                "parent": "phase_4j64",
                "key": key,
                "request_identity": digest(wanted),
                "parent_archive_sha256": pin(old / "native-outputs.tar.xz"),
                "producer_result_sha256": pin(
                    archive / "outputs" / key / "result.json"
                ),
            }
    require(
        len(primary) == 17 and len(reusable) == 7, "explicit_reuse_inventory_changed"
    )
    specs = read(old / "plan/proposed-independent-site-check.json")["proposals"]
    calibration = read(old / "baseline-calibration.json")["conditional_candidates"][
        "QZVP_standard"
    ]
    predictions = []
    control_geometries = []
    for proposal in specs:
        spec = proposal["spec"]
        null = unknown_null_check(model, x, spec, 0.003)
        column = NAMES[0] if "SiC" in spec["name"] else NAMES[1]
        oldspec = next(
            s
            for s in p["paths"]
            if s["name"] == ("SiC" if "SiC" in spec["name"] else "CH")
        )
        newpoint = path(model, x, spec, 0.003)["coordinates"]
        oldpoint = path(model, x, oldspec, 0.003)["coordinates"]
        predictions.append(
            {
                "spec": spec,
                "unknown_null": null,
                "fixed_total_curvature": null["verified_fixed_response"][
                    "curvature_from_path_gradient"
                ],
                "predicted_QZVP_total_curvature": null["verified_fixed_response"][
                    "curvature_from_path_gradient"
                ]
                + 2 * calibration[column],
                "fixed_correction": calibration[column],
                "correction_refitted": False,
                "same_site_path_equivalent_at1e_6A": equivalent(
                    model, oldpoint, newpoint, 1e-6
                ),
                "same_site_path_equivalent_at1e_3A": equivalent(
                    model, oldpoint, newpoint, 1e-3
                ),
            }
        )
        control_geometries.extend(
            path(model, x, spec, s * 0.003)["coordinates"] for s in (1, -1)
        )
    # Geometry-only derivative diagnostic; no QM values informed its choice.
    direction = np.random.default_rng(4065).normal(size=(27, 3))
    direction -= direction.mean(axis=0)
    direction /= np.linalg.norm(direction)
    fdcenter = np.array(joints[0]["coordinates"])
    fdgeos = [
        fdcenter + sign * h * direction for h in (0.001, 0.0005) for sign in (1, -1)
    ]
    leak = leakage(
        model,
        [r["coordinates"] for r in primary + joints] + control_geometries + fdgeos,
        validation,
    )
    stat = {
        k: int(v)
        for k, v in (
            s.split()
            for s in Path("/sys/fs/cgroup/memory.stat").read_text().splitlines()
        )
    }
    resources = {
        "memory_max_bytes": int(Path("/sys/fs/cgroup/memory.max").read_text()),
        "memory_current_bytes": int(Path("/sys/fs/cgroup/memory.current").read_text()),
        "memory_stat_bytes": stat,
        "cpu_max": Path("/sys/fs/cgroup/cpu.max").read_text().strip(),
        "disk_available_bytes": shutil.disk_usage(work.parent).free,
        "retained_native_timings": read(old / "compute-accounting.json")[
            "gradient_elapsed_seconds"
        ],
        "QZVP_basis": p["basis_composition"]["def2-qzvp"],
        "QZVP_allocation": "One worker only. J64 QZVP399->1079 orbital functions peaked~9.5GiB RSS with8GiB configured; allow12GiB RSS plus retained~6GiB unrelated anonymous memory within32GiB cgroup. Dense QZVP may cost more than682s center;7200s bound. Never concurrent QZVP/TZVP. Reclaimable cache distinguished from anonymous memory.",
    }
    plan = signed(
        {
            "schema": "island.hmds-method-specific-experiment.v1",
            "base_commit": BASE,
            "J64_prerequisite": PREREQUISITE,
            "J64_in_main": True,
            "branch": BRANCH,
            "J64_receipt_sha256": pin(old / "receipt.json"),
            "J64_reconstructed_before_new_execution": True,
            "graph_identity": model["graph_identity"],
            "atom_order": p["atom_order"],
            "elements": p["elements"],
            "charge": 0,
            "multiplicity": 1,
            "center_coordinates": x.tolist(),
            "protocol": target,
            "paths": p["paths"],
            "basis_composition": p["basis_composition"],
            "scope": "Exact27atom HMDS method-specific research target, no physical/basis/transfer qualification. SiH stopped. Historical gates/author provenance unchanged.",
            "hypotheses": [
                "QZVP cross-site curvatures test frozen shared corrections, not independent basins.",
                "Matched QZVP center grid contrast separates integration from orbital-basis effects.",
                "A complete fixed+58unknown MM model may reproduce one dense-TZVP target; source-baseline incompatibility remains a competing hypothesis.",
            ],
            "limits": {
                "maximum_new_gradients": 32,
                "planned_maximum": 30,
                "campaign_wall_seconds": 36000,
                "maximum_parallel_workers": 1,
                "threads": 4,
                "TZVP_memory_GiB": 4,
                "QZVP_memory_GiB": 8,
                "TZVP_wall_seconds": 1800,
                "QZVP_wall_seconds": 7200,
                "scratch_GiB": 64,
            },
            "resources": resources,
            "primary": primary,
            "reuse": reusable,
            "joints": joints,
            "validation_coordinates": validation,
            "validation_leakage": leak,
            "cross_site_predictions": predictions,
            "parameter_specifications": specifications,
            "fit": {
                "parameters": names,
                "energy_RMSE_scale_kcal_mol": 0.1,
                "force_RMS_scale_kcal_mol_A": 0.5,
                "objective": "mean((E_QM-E_fixed-Ap-offset)^2)/0.1^2 + mean(Cartesian_force_residual^2)/0.5^2; all81 components, one offset on training only",
                "offset": "one unbounded additive training intercept in energy only; frozen unchanged on validation; subtract center electronic energy before kcal conversion for precision",
                "optimizer": "scipy.optimize.lsq_linear TRF bounded linear least squares, column-scaled including intercept, tolerance1e-12,maxiter1000; no regularization",
                "rank_relative_tolerance": 1e-10,
                "maximum_scaled_condition": 10000.0,
                "condition_scope": "Actual weighted59column matrix including offset; 10000 is an experiment identifiability screen, not prior987 geometry-only result or a physical certificate.",
                "energy_RMSE_max": 0.1,
                "force_RMS_max": 0.5,
                "maximum_absolute_energy_error": 0.3,
                "maximum_local_force_error": 2.0,
                "bound_activity_fraction": 1e-6,
                "curvature_constraints": "Diagnostics only; selected-protocol baseline estimates not QZVP. No duplicate gradient-derived observations in objective.",
                "uncertainty": {
                    "replicates": 64,
                    "seed": 4065,
                    "unit": "whole training configuration/path blocks; all signs/steps of each primary path travel together; each joint and center separate blocks",
                    "minimum_identifiable_fraction": 0.8,
                    "prediction_energy_RMS95_max": 0.1,
                    "prediction_force_RMS95_max": 0.5,
                    "interval": "2.5/97.5percentiles only on explicitly eligible bootstrap replicates; retain every rank/conditioning/bound failure, no replacement resampling; failed fraction blocks eligibility. Absolute intervals for all coefficients incl near zero, covariance/correlations and prediction spread; limited small-block/model-discrepancy uncertainty, no physical confidence claim.",
                },
                "candidate_activation": "Full matching primary+joints, numerical gates, rank59/condition<=10000, solver success, no active bound, training thresholds and uncertainty all pass before untouched validation QM. A failed historical complete-basis gate alone never blocks this method-specific trial.",
                "acceptance": "Same frozen energy/force/max gates on untouched validation plus predictive uncertainty; method-specific research only, never automatic assignment.",
            },
            "precision": {
                "inherited": p["precision"],
                "FD_direction": direction.tolist(),
                "FD_center_id": joints[0]["id"],
                "FD_steps_A": [0.001, 0.0005],
                "FD_absolute_floor_kcal_mol_A": 0.01,
                "FD_relative": 0.001,
                "grid_energy_floor_kcal_mol": 0.002,
                "SCF_energy_floor_kcal_mol": 0.00005,
                "force_RMS_floor": 0.005,
                "force_local_floor": 0.02,
                "relative": 0.005,
                "policy": "Mandatory native convergence,mapping,units and energy-gradient checks for fitting observations; two planned common-geometry numerical controls at center/joint040. Preserve original basis failures separately.",
            },
            "activation": {
                "controls": "Five required QZVP controls, prediction-only sites; no refit.",
                "primary": "Ten missing exact dense-TZVP primary observations after controls native success; QZVP basis/grid failures do not gate a distinct target fit.",
                "joints": "Five geometry-screened joints if matching primary mapping/units/convergence and inherited path consistency pass.",
                "precision": "Two300x1454 grid comparisons, one tighterSCF center, four finite-difference points around joint040; within32cap. No other protocol ladder.",
                "validation": "Three reserved configurations only after eligible training candidate; recheck leakage against every observed geometry. Otherwise remain untouched.",
            },
            "failure_policy": "Every process attempt counts, including startup/timeout/interruption. Use existing run_job/capture hashes, retain commands/stdout/stderr/native errors/scratch. Resume only exact verified completed requests. A retry needs additive identity/amendment and remaining32cap; never erase failure.",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "execution_sources": {
                n: pin(repo / n)
                for n in [
                    "scripts/phase_4j65.py",
                    "scripts/phase_4j64.py",
                    "src/island/pcff_qm_pilot/hmds_precision.py",
                    "src/island/pcff_qm_pilot/hmds_engine.py",
                    "src/island/pcff_qm_pilot/control_engine.py",
                ]
            },
        }
    )
    jobs = {}

    def add(key, xyz, stage, **changes):
        r = request(plan, xyz, **changes)
        qz = r["basis"] == "def2-qzvp"
        jobs[key] = {
            "stage": stage,
            "request_identity": digest(r),
            "memory_GiB": 8 if qz else 4,
            "wall_seconds": 7200 if qz else 1800,
        }
        write(work / "requests" / (key + ".json"), r)

    for proposal in specs:
        for sign, label in ((1, "plus"), (-1, "minus")):
            add(
                "controls/" + proposal["spec"]["name"] + "-" + label,
                path(model, x, proposal["spec"], sign * 0.003)["coordinates"],
                "controls",
                basis="def2-qzvp",
                grid=[150, 974],
                thresholds=[1e-10, 1e-8],
            )
    add(
        "controls/QZVP-dense-center",
        x,
        "controls",
        basis="def2-qzvp",
        thresholds=[1e-10, 1e-8],
    )
    for row in primary:
        if row["id"] not in reusable:
            add("training/" + row["id"], row["coordinates"], "primary")
    for row in joints:
        add("training/" + row["id"], row["coordinates"], "joints")
    for key, geo in (("center", x), ("joint040", fdcenter)):
        add("precision/grid-" + key, geo, "precision", grid=[300, 1454])
    add("precision/scf-center", x, "precision", thresholds=[1e-13, 1e-11])
    for h in (0.001, 0.0005):
        for sign, label in ((1, "plus"), (-1, "minus")):
            add(
                f"precision/FD-{h}-{label}",
                fdcenter + sign * h * direction,
                "precision",
            )
    for i, geo in enumerate(validation):
        add(f"validation/config-{i}", geo, "validation")
    require(len(jobs) == 30, "bounded30_requests")
    inventory = signed({"plan_identity": plan["identity"], "jobs": jobs})
    for directory in (root / "plan", work):
        write(directory / "scientific-plan.json", plan)
        write(directory / "jobs.json", inventory)
    write(
        root / "plan/requests.json",
        {k: read(work / "requests" / (k + ".json")) for k in jobs},
    )
    write(
        root / "plan/base-verification.json",
        signed(
            {
                "base": BASE,
                "prerequisite": PREREQUISITE,
                "commands": [
                    {
                        "argv": args,
                        "output": subprocess.check_output(args, text=True).strip(),
                    }
                    for args in (
                        ["git", "rev-parse", "origin/main"],
                        ["git", "branch", "--show-current"],
                        ["git", "merge-base", BASE, PREREQUISITE],
                    )
                ],
            }
        ),
    )
    for row in primary:
        if row["id"] in reusable:
            protocol_contrast(
                request(plan, row["coordinates"]),
                read(archive / "requests" / (reusable[row["id"]]["key"] + ".json")),
                [],
            )
    print(
        json.dumps(
            {
                "plan_sha256": pin(root / "plan/scientific-plan.json"),
                "planned": 30,
                "reuse": len(reusable),
            }
        )
    )


def execute(repo, root, work, python, stage):
    plan = read(work / "scientific-plan.json")
    jobs = read(work / "jobs.json")["jobs"]
    if stage in ("primary", "joints", "precision", "validation"):
        for key, job in jobs.items():
            if job["stage"] == "controls":
                verify_completed(work / "outputs" / key, job["request_identity"])
    if stage in ("joints", "validation"):
        decision = read(
            root
            / (
                "primary-numerical-gates.json"
                if stage == "joints"
                else "trial-fit.json"
            )
        )
        require(
            decision["scientific_plan_identity"] == plan["identity"],
            "decision_plan_binding",
        )
        require(
            decision["joints_authorized"]
            if stage == "joints"
            else decision["validation_authorized"],
            "method_specific_stage_not_authorized",
        )
    if stage == "precision":
        key = "training/" + plan["precision"]["FD_center_id"]
        verify_completed(work / "outputs" / key, jobs[key]["request_identity"])
    for key, job in jobs.items():
        if job["stage"] == stage:
            code = run_job(repo, root, work, python, key)
            if code:
                return code
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=["prepare", "controls", "primary", "joints", "precision", "validation"],
    )
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j65")
    )
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--python", type=Path)
    a = parser.parse_args()
    repo, root, work = a.repo.resolve(), a.evidence.resolve(), a.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if a.stage == "prepare":
        prepare(repo, root, work)
    else:
        raise SystemExit(execute(repo, root, work, a.python.resolve(), a.stage))
