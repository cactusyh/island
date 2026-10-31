"""Execute/resume only the frozen J62 control and graph-bound HMDS entry."""

import argparse
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.campaign import verify_completed
from island.pcff_qm_pilot.control_engine import request
from island.pcff_qm_pilot.hamiltonian import length
from island.pcff_qm_pilot.molecular_control import cartesian, rigid_projection, torsion
from island.pcff_qm_pilot.molecular_native import normalized


def read(p):
    return json.loads(Path(p).read_text())


def run_process(work, repo, python, key, operation="gradient"):
    out = work / "outputs" / key
    req = work / "requests" / (key + ".json")
    if operation == "gradient" and out.exists():
        return verify_completed(out, digest(read(req)))
    if operation == "optimize" and (work / "hmds-optimization/result.json").exists():
        return read(work / "hmds-optimization/result.json")
    log = work / "processes" / key
    require(not log.exists(), "additive_process_retry_required")
    log.mkdir(parents=True)
    scratch = work / "scratch" / key
    scratch.mkdir(parents=True, exist_ok=True)
    require(
        sum(p.stat().st_size for p in (work / "scratch").rglob("*") if p.is_file())
        < 20 * 1024**3,
        "scratch_cap",
    )
    cmd = [
        str(python),
        "-m",
        "island.pcff_qm_pilot.control_engine",
        operation,
        "--work",
        str(work),
    ]
    if operation == "gradient":
        cmd += ["--key", key]
    write(
        log / "command.json",
        {
            "argv": cmd,
            "cwd": str(scratch),
            "allowed_environment": {
                "PYTHONPATH": str(repo) + ":" + str(repo / "src"),
                "OMP_NUM_THREADS": "4",
                "OPENBLAS_NUM_THREADS": "1",
                "PSI_SCRATCH": str(scratch),
            },
        },
    )
    env = {
        **os.environ,
        "PYTHONPATH": str(repo) + ":" + str(repo / "src"),
        "OMP_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "1",
        "PSI_SCRATCH": str(scratch),
    }
    start = time.monotonic()
    with (
        (log / "stdout.log").open("w") as stdout,
        (log / "stderr.log").open("w") as stderr,
    ):
        try:
            r = subprocess.run(
                cmd,
                cwd=scratch,
                env=env,
                stdout=stdout,
                stderr=stderr,
                timeout=9000 if operation == "optimize" else 1200,
                check=False,
            )
            code = r.returncode
        except subprocess.TimeoutExpired:
            code = 124
    result = {
        "returncode": code,
        "elapsed_seconds": time.monotonic() - start,
        "stdout_sha256": pin(log / "stdout.log"),
        "stderr_sha256": pin(log / "stderr.log"),
        "native_output_directory": (
            "hmds-optimization"
            if operation == "optimize"
            else str(out.relative_to(work))
        ),
        "temporary_cleanup_performed": False,
    }
    write(log / "result.json", signed(result))
    print(json.dumps({"key": key, **result}), flush=True)
    require(code == 0, "native_subprocess_failure_retained", key=key, returncode=code)
    if operation == "gradient":
        normalized(out, read(req))
    return result


def hmds(work, root, repo, python):
    run_process(work, repo, python, "hmds-optimization", "optimize")
    result = read(work / "hmds-optimization/result.json")
    require(result["calls"], "HMDS_no_successful_gradient")
    key = result["calls"][-1]
    sample = normalized(
        work / "outputs" / key, read(work / "requests" / (key + ".json"))
    )
    x = np.array(sample["coordinates"])
    spec = read(work / "scientific-plan.json")["HMDS"]
    native = read(work / "outputs" / key / "result.json")
    g = np.load(work / "outputs" / key / "gradient.npy", allow_pickle=False)
    stationary = (
        result["status"] == "completed"
        and np.max(abs(g)) <= 2.5e-5
        and np.sqrt(np.mean(g * g)) <= 1e-5
    )
    domain = read(root / "HMDS-domain-audit.json")
    lengths = [np.linalg.norm(x[i - 1] - x[j - 1]) for i, j, _ in domain["bonds"]]
    require(min(lengths) > 0.8 and max(lengths) < 2.2, "HMDS_connectivity_lost")
    u, proj = rigid_projection(x)
    directions = []
    for seed in spec["response_seeds"]:
        jets = cartesian(x)
        d = (
            length(jets, seed["sites"])
            if seed["kind"] == "bond"
            else torsion(jets, seed["sites"])
        ).d
        d = proj @ d
        directions.append(d / np.linalg.norm(d))
    write(
        root / "HMDS-center.json",
        signed(
            {
                "center": sample,
                "geometry_identity": digest(x.tolist()),
                "optimization_result": result,
                "stationary_candidate": bool(stationary),
                "full_minimum_verified": False,
                "gradient_max_hartree_bohr": float(np.max(abs(g))),
                "gradient_RMS_hartree_bohr": float(np.sqrt(np.mean(g * g))),
                "basis_functions": native["basis_functions"],
                "basis_shells": native["basis_shells"],
                "directions": np.array(directions).T.tolist(),
                "rigid_overlap_max": float(np.max(abs(u.T @ np.array(directions).T))),
                "interpretation": "Exploratory SVP stationary candidate only if optimizer and gradient criteria pass. Responses are matched to the final successful geometry even if unfinished.",
            }
        ),
    )
    jobs = {"hmds-reference/tzvp-center": request(spec, x, basis="def2-tzvp")}
    for d, vector in enumerate(directions):
        for sign in (-1, 1):
            xyz = x + (sign * 0.003 * vector).reshape(27, 3)
            jobs[f"hmds-response/d{d}-s{sign:+d}"] = request(spec, xyz)
            if d == 0:
                jobs[f"hmds-reference/tzvp-d{d}-s{sign:+d}"] = request(
                    spec, xyz, basis="def2-tzvp"
                )
    write(
        root / "plan/HMDS-response-jobs.json",
        signed(
            {
                "center_identity": digest(x.tolist()),
                "jobs": {k: digest(v) for k, v in jobs.items()},
                "predeclared_rules_sha256": pin(root / "plan/scientific-plan.json"),
            }
        ),
    )
    for k, v in jobs.items():
        p = work / "requests" / (k + ".json")
        if p.exists():
            require(read(p) == v, "HMDS_resume_request_mutation")
        else:
            write(p, v)
        run_process(work, repo, python, k)


def execute(work, root, repo, python):
    require(
        pin(work / "scientific-plan.json") == pin(root / "plan/scientific-plan.json"),
        "frozen_plan_mutation",
    )
    for path, sha in read(root / "plan/implementation-pins.json")["files"].items():
        require(pin(repo / path) == sha, "executed_implementation_mutation", file=path)
    jobs = read(root / "plan/control-jobs.json")["jobs"]

    def control(item):
        k, sha = item
        require(
            digest(read(work / "requests" / (k + ".json"))) == sha,
            "control_request_mutation",
        )
        return run_process(work, repo, python, k)

    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(hmds, work, root, repo, python)] + [
            pool.submit(control, item) for item in jobs.items()
        ]
        errors = []
        for i, f in enumerate(futures):
            try:
                f.result()
            except Exception as exc:  # noqa: BLE001 -- retain every worker failure
                errors.append(
                    {"task": i, "error_type": type(exc).__name__, "error": str(exc)}
                )
    write(
        root / "execution-outcome.json",
        signed(
            {
                "failures": errors,
                "all_tasks_returned": True,
                "ledger_sha256": pin(work / "ledger.json"),
            }
        ),
    )
    require(not errors, "preserved_campaign_failures", errors=errors)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j62"))
    p.add_argument("--python", type=Path, required=True)
    a = p.parse_args()
    execute(
        a.work.resolve(), a.evidence.resolve(), a.repo.resolve(), a.python.resolve()
    )
