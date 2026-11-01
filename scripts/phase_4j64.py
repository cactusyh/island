"""Prepare and execute the bounded HMDS precision experiment using J63 capture."""

import argparse
import json
import os
import shutil
import signal
import subprocess
import time
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.campaign import verify_completed
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_hamiltonian import load_model, parts
from island.pcff_qm_pilot.hmds_paths import path
from island.pcff_qm_pilot.hmds_precision import check_budget, read

BASE = "3fd458b33fe92e9bd206d82e2b7d9cdbb7cb75e3"
PREREQUISITE = "c10e76f05316c05e66d1d4c0df6581af775572cd"


def prepare(repo, root, work):
    from scripts.replay_phase_4j63 import reconstruct

    require(not (root / "plan/scientific-plan.json").exists(), "frozen_plan_exists")
    previous = repo / "docs/evidence/phase_4j63"
    reconstruction = reconstruct(repo, previous)
    p = read(previous / "plan/scientific-plan.json")
    model = load_model(repo)
    x = np.array(read(previous / "stationary-center.json")["center"]["coordinates"])
    columns = parts(model, x)["columns"]
    require(
        columns == read(previous / "interaction-accounting.json")["columns"]
        and len(columns) == 56,
        "retained56_column_mapping_mutation",
    )
    write(
        root / "J63-reconstruction.json",
        signed(
            {
                **reconstruction,
                "columns": columns,
                "unknown_null_checks": read(
                    previous / "source-baseline-path-curvature.json"
                )["paths"],
                "prerequisite": PREREQUISITE,
            }
        ),
    )
    stat = dict(
        line.split()
        for line in Path("/sys/fs/cgroup/memory.stat").read_text().splitlines()
    )
    resource = {
        "memory_max_bytes": int(Path("/sys/fs/cgroup/memory.max").read_text()),
        "memory_current_bytes": int(Path("/sys/fs/cgroup/memory.current").read_text()),
        "memory_stat_bytes": {k: int(v) for k, v in stat.items()},
        "cpu_max": Path("/sys/fs/cgroup/cpu.max").read_text().strip(),
        "disk_available_bytes": shutil.disk_usage(work.parent).free,
        "historical_timing_seconds": {
            k: v["elapsed_seconds"]
            for k, v in read(
                repo.parent / "island-validation/phase4j63-qm/ledger.json"
            )["calls"].items()
            if k.startswith("controls/")
        },
        "QZVP_full_DF_three_index_tensor_bytes": 1079**2 * 1107 * 8,
        "QZVP_allocation_reason": "32GiB cgroup; about6GiB unrelated anonymous memory, reclaimable file/slab cache. One process only,8GiB Psi4 with12GiB process allowance; integral blocking/disk needed for full9.6GiB DF tensor. No shared-process modification. Native RSS and DF storage retained. Four threads within8CPU quota.553-function TZVPP calls~280-300s; cubic dimension extrapolation~2200s;7200s cap permits controlled overhead, not an unlimited campaign.",
    }
    write(root / "plan/resource-preflight.json", signed(resource))
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    cross = np.array(
        [[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]]
    )
    angle = 0.37
    rotation = (
        np.eye(3) * np.cos(angle)
        + (1 - np.cos(angle)) * np.outer(axis, axis)
        + np.sin(angle) * cross
    )
    plan = signed(
        {
            "schema": "island.hmds-precision-plan.v1",
            "base_commit": BASE,
            "prerequisite": PREREQUISITE,
            "prerequisite_in_main": True,
            "graph_identity": model["graph_identity"],
            "atom_order": p["atom_order"],
            "elements": p["elements"],
            "charge": 0,
            "multiplicity": 1,
            "center_coordinates": x.tolist(),
            "center_identity": read(previous / "stationary-center.json")[
                "geometry_identity"
            ],
            "protocol": p["protocol"],
            "basis_composition": p["basis_composition"],
            "paths": p["paths"],
            "precision": p["precision"],
            "proper_rigid_rotation": rotation.tolist(),
            "hypotheses": [
                "Matched250x1454 grid resolves or bounds retained siloxane energy/gradient mismatch; rounding alone is quantified independently.",
                "Common-center QZVP separates basis-force shifts from geometry/grid/SCF effects; not a stationary or all-path basis claim.",
                "Unknown-null SiC/CH curvatures support only two conditional quadratic corrections if their objective-specific precision passes.",
            ],
            "limits": {
                "maximum_new_gradients": 24,
                "planned_maximum": 21,
                "campaign_wall_seconds": 36000,
                "maximum_parallel_workers": 1,
                "threads": 4,
                "TZVP_memory_GiB": 4,
                "QZVP_memory_GiB": 8,
                "TZVP_wall_seconds": 1800,
                "QZVP_wall_seconds": 7200,
                "scratch_GiB": 64,
                "unallocated": 3,
            },
            "activation": {
                "dense": "Required3 jobs; identical retained SCF-tight150x974 triple except grid250x1454.",
                "qz-center": "Required1 job; identical retained primary/TZVPP center except orbital basis; one worker and separate8GiB/7200s allocation.",
                "stretch": "8 matched tight-SCF150x974/250x1454 stretch +/-0.003A jobs, after dense and QZ center success. Separate SCF and grid contrasts; no equilibrium release.",
                "qz-stretch": "4 standard-SCF150x974 QZVP +/-0.003A jobs after center native success, within resource/campaign cap; objective is curvature/bond-force basis sensitivity, not repair soft torsion.",
                "torsion-control": "Only if dense energy/gradient consistency still fails:2 dense +/-0.04rad points and3 properly whole-molecule-rotated copies of dense center/+/-0.02rad. Isolates step/grid-orientation effects, no grid ladder.",
            },
            "baseline": {
                "identity": "island_research_class2_hmds_j64_quadratic_v1",
                "equation": "E_corr = delta_K2_SiC sum_SiC(r-r0_source)^2 + delta_K2_CH sum_CH(r-r0_source)^2; no1/2, source K3/K4 and equilibrium unchanged",
                "units": "kcal/mol/angstrom^2",
                "curvature_sensitivity": [[2, 0], [0, 2]],
                "selection": "Conditional curvature calibration only; require retained two-step and matched grid/SCF/basis CURVATURE gates. Full Cartesian-force gates reported independently and may fail. No physical56-vector, complete force fit, automatic assignment, or validation claim.",
                "uncertainty": "Report each step/protocol difference as sensitivity, not a statistical confidence interval or complete-basis bound.",
            },
            "rounding": "Decimal native total energies at50-digit subtraction before conversion, independent ASCII gradients; separately bounded text-rounding, binary serialization and offset-before/after conversion errors.",
            "design": "Geometry-only58-column screen using original17 primary geometries and five selected J63 joints; diagnostic scales unchanged. No joint or reserved validation QM.",
            "failure_policy": "Count all subprocess attempts including initialization/timeout failures; native capture also counts started gradients. Preserve requests,stdout,stderr,time-v,native outputs and scratch before cleanup. Nonzero unexpected/timeout statuses remain failures; exact completed/hash-verified jobs only resume. Retry requires additive identity and remaining24-cap.",
            "new_SiH_QM": False,
            "joint_QM": False,
            "validation_QM": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "J63_reconstruction_sha256": pin(root / "J63-reconstruction.json"),
            "execution_sources": {
                n: pin(repo / n)
                for n in [
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

    siloxane = next(s for s in p["paths"] if s["name"] == "siloxane")
    tight = {"thresholds": [1e-12, 1e-10]}
    for name, s in [("center", 0), ("plus", 0.02), ("minus", -0.02)]:
        xyz = path(model, x, siloxane, s)["coordinates"]
        add("dense/" + name, xyz, "dense", grid=[250, 1454], **tight)
        add(
            "rigid/" + name,
            xyz @ rotation.T,
            "torsion-control",
            grid=[250, 1454],
            **tight,
        )
    add("qz/center", x, "qz-center", basis="def2-qzvp")
    for spec in p["paths"]:
        if spec["kind"] != "stretch":
            continue
        for sign, label in ((1, "plus"), (-1, "minus")):
            xyz = path(model, x, spec, sign * 0.003)["coordinates"]
            add(f"tight/{spec['name']}-{label}", xyz, "stretch", **tight)
            add(
                f"dense/{spec['name']}-{label}",
                xyz,
                "stretch",
                grid=[250, 1454],
                **tight,
            )
            add(f"qz/{spec['name']}-{label}", xyz, "qz-stretch", basis="def2-qzvp")
    for sign, label in ((1, "plus"), (-1, "minus")):
        add(
            "dense/coarse-" + label,
            path(model, x, siloxane, sign * 0.04)["coordinates"],
            "torsion-control",
            grid=[250, 1454],
            **tight,
        )
    require(len(jobs) == 21, "bounded_precision_inventory")
    inventory = signed({"plan_sha256": digest(plan), "jobs": jobs})
    for directory in (root / "plan", work):
        write(directory / "scientific-plan.json", plan)
        write(directory / "jobs.json", inventory)
    write(
        root / "plan/request-inventory.json",
        {k: read(work / "requests" / (k + ".json")) for k in jobs},
    )
    print(
        json.dumps(
            {"prepared": len(jobs), "plan_sha256": pin(work / "scientific-plan.json")}
        )
    )


def run_job(repo, root, work, python, key):
    plan = read(work / "scientific-plan.json")
    require(
        pin(work / "scientific-plan.json") == pin(root / "plan/scientific-plan.json"),
        "frozen_precision_plan_mutation",
    )
    for name, sha in plan["execution_sources"].items():
        if pin(repo / name) != sha:
            amendment = read(root / "plan/startup-maintenance.json")
            binding = amendment["sources"][name]
            require(
                binding["old_sha256"] == sha
                and binding["new_sha256"] == pin(repo / name)
                and pin(root / "plan" / binding["original_archive"])
                == binding["archive_sha256"],
                "executed_precision_source_mutation",
                file=name,
            )
    job = read(work / "jobs.json")["jobs"][key]
    r = read(work / "requests" / (key + ".json"))
    require(digest(r) == job["request_identity"], "frozen_precision_request_mutation")
    process = work / "processes" / key
    if process.exists():
        result = read(process / "result.json")
        for name, sha in result["files"].items():
            require(pin(process / name) == sha, "precision_process_record_mutation")
        if result["returncode"] == 0:
            verify_completed(work / "outputs" / key, digest(r))
            return 0
        amendment = read(root / "plan/startup-maintenance.json")
        require(key == amendment["retry_job"], "additive_retry_identity_required")
        process = work / "processes" / amendment["retry_attempt"]
        if process.exists():
            retry = read(process / "result.json")
            require(retry["returncode"] == 0, "additive_retry_identity_required")
            for name, sha in retry["files"].items():
                require(pin(process / name) == sha, "precision_process_record_mutation")
            verify_completed(work / "outputs" / key, digest(r))
            return 0
        require(
            not (work / "outputs" / key).exists(), "additive_retry_identity_required"
        )
    attempt_key = process.relative_to(work / "processes").as_posix()
    ledger = work / "attempts.json"
    a = read(ledger) if ledger.exists() else {"started_unix": time.time(), "calls": {}}
    check_budget(plan, a["calls"], time.time() - a["started_unix"])
    require(attempt_key not in a["calls"], "additive_retry_identity_required")
    a["calls"][attempt_key] = {
        "job_key": key,
        "request_identity": digest(r),
        "status": "reserved_or_interrupted",
        "started_unix": time.time(),
    }
    write(ledger, a)
    process.mkdir(parents=True)
    scratch = work / "scratch" / key
    scratch.mkdir(parents=True, exist_ok=True)
    require(
        sum(p.stat().st_size for p in (work / "scratch").rglob("*") if p.is_file())
        < 64 * 1024**3,
        "precision_scratch_cap",
    )
    command = [
        str(python),
        "-m",
        "island.pcff_qm_pilot.hmds_precision",
        "--work",
        str(work),
        "--key",
        key,
    ]
    env = {
        **os.environ,
        "PYTHONPATH": f"{repo}:{repo / 'src'}",
        "OMP_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "1",
        "PSI_SCRATCH": str(scratch),
    }
    write(
        process / "command.json",
        {
            "argv": command,
            "cwd": str(scratch),
            "allowed_environment": {
                k: env[k]
                for k in (
                    "PYTHONPATH",
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "PSI_SCRATCH",
                )
            },
        },
    )
    start = time.monotonic()
    with (
        (process / "stdout.log").open("w") as out,
        (process / "stderr.log").open("w") as err,
    ):
        try:
            worker = subprocess.Popen(
                command,
                cwd=scratch,
                env=env,
                stdout=out,
                stderr=err,
                start_new_session=True,
            )
            code = worker.wait(timeout=job["wall_seconds"])
        except subprocess.TimeoutExpired:
            # Retain timeout status and terminate the complete worker group.
            os.killpg(worker.pid, signal.SIGTERM)
            try:
                worker.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(worker.pid, signal.SIGKILL)
                worker.wait()
            code = 124
        except OSError as exc:
            code = 127
            err.write(f"{type(exc).__name__}: {exc}\n")
    result = signed(
        {
            "returncode": code,
            "elapsed_seconds": time.monotonic() - start,
            "status": "completed"
            if code == 0
            else "timeout"
            if code == 124
            else "unexpected_failure",
            "files": {p.name: pin(p) for p in process.iterdir() if p.is_file()},
            "scratch_retained": True,
        }
    )
    write(process / "result.json", result)
    a = read(ledger)
    a["calls"][attempt_key].update(
        status=result["status"],
        returncode=code,
        elapsed_seconds=result["elapsed_seconds"],
    )
    write(ledger, a)
    print(
        json.dumps(
            {
                "key": key,
                "returncode": code,
                "elapsed_seconds": result["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return code


def execute(repo, root, work, python, stage):
    jobs = read(work / "jobs.json")["jobs"]
    if stage not in ("dense", "qz-center"):
        for key in ("dense/center", "dense/plus", "dense/minus", "qz/center"):
            verify_completed(work / "outputs" / key, jobs[key]["request_identity"])
    if stage == "torsion-control":
        decision = read(root / "dense-grid-decision.json")
        require(
            decision["energy_gradient_consistency"]["passed"] is False,
            "unjustified_extra_torsion_jobs",
        )
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
        choices=[
            "prepare",
            "dense",
            "qz-center",
            "stretch",
            "qz-stretch",
            "torsion-control",
        ],
    )
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/evidence/phase_4j64")
    )
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--python", type=Path)
    args = parser.parse_args()
    repo, root, work = args.repo.resolve(), args.evidence.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if args.stage == "prepare":
        prepare(repo, root, work)
    else:
        raise SystemExit(execute(repo, root, work, args.python.resolve(), args.stage))
