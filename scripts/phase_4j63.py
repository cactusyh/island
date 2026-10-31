"""Execute only the frozen J63 HMDS stages; each numerical failure stays nonzero."""

import argparse
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.campaign import verify_completed
from island.pcff_qm_pilot.molecular_native import normalized


def read(p):
    return json.loads(Path(p).read_text())


def run_process(repo, work, python, key, operation="gradient"):
    plan = read(work / "scientific-plan.json")
    directory = work / ("optimization" if operation == "optimize" else "outputs/" + key)
    result = directory / "result.json"
    process = work / "processes" / key
    if result.exists():
        r = read(result)
        if operation == "gradient":
            verify_completed(
                directory, digest(read(work / "requests" / (key + ".json")))
            )
        else:
            require(
                r["identity"]
                == digest({k: v for k, v in r.items() if k != "identity"}),
                "optimizer_status_mutation",
            )
            for name, sha in r["native_files"].items():
                require(pin(directory / name) == sha, "optimizer_native_mutation")
            for k in r["calls"]:
                normalized(
                    work / "outputs" / k, read(work / "requests" / (k + ".json"))
                )
        return read(process / "result.json")["returncode"]
    require(not process.exists(), "additive_retry_identity_required")
    process.mkdir(parents=True)
    scratch = work / "scratch" / key
    scratch.mkdir(parents=True)
    require(
        sum(p.stat().st_size for p in (work / "scratch").rglob("*") if p.is_file())
        < 24 * 1024**3,
        "scratch_cap",
    )
    cmd = [
        str(python),
        "-m",
        "island.pcff_qm_pilot.hmds_engine",
        operation,
        "--work",
        str(work),
    ]
    if operation == "gradient":
        cmd += ["--key", key]
    env = {
        **os.environ,
        "PYTHONPATH": str(repo) + ":" + str(repo / "src"),
        "OMP_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "1",
        "PSI_SCRATCH": str(scratch),
    }
    write(
        process / "command.json",
        {
            "argv": cmd,
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
            "coordinator_sha256": pin(__file__),
        },
    )
    start = time.monotonic()
    timeout = plan["limits"][
        "optimizer_wall_seconds" if operation == "optimize" else "gradient_wall_seconds"
    ]
    timed_out = False
    with (
        (process / "stdout.log").open("w") as stdout,
        (process / "stderr.log").open("w") as stderr,
    ):
        try:
            code = subprocess.run(
                cmd,
                cwd=scratch,
                env=env,
                stdout=stdout,
                stderr=stderr,
                timeout=timeout,
                check=False,
            ).returncode
        except subprocess.TimeoutExpired:
            code = 124
            timed_out = True
    write(
        process / "result.json",
        signed(
            {
                "returncode": code,
                "timed_out": timed_out,
                "elapsed_seconds": time.monotonic() - start,
                "stdout_sha256": pin(process / "stdout.log"),
                "stderr_sha256": pin(process / "stderr.log"),
                "native_directory": str(directory.relative_to(work)),
                "cleanup_performed": False,
            }
        ),
    )
    print(json.dumps({"key": key, "returncode": code}), flush=True)
    return code


def execute(repo, evidence, work, python, stage):
    require(
        pin(work / "scientific-plan.json")
        == pin(evidence / "plan/scientific-plan.json"),
        "frozen_plan_mutation",
    )
    from scripts.validate_phase_4j63_receipt import execution_source

    execution_source(repo, evidence)
    if stage == "optimize":
        return run_process(repo, work, python, "optimization", "optimize")
    from scripts.prepare_phase_4j63 import deformation_jobs

    jobs = deformation_jobs(repo, evidence, work)
    with ThreadPoolExecutor(max_workers=2) as pool:
        codes = list(pool.map(lambda k: run_process(repo, work, python, k), jobs))
    outcome = signed(
        {
            "stage": "deformations",
            "process_codes": dict(zip(jobs, codes, strict=True)),
            "all_completed": not any(codes),
            "ledger_sha256": pin(work / "ledger.json"),
        }
    )
    write(evidence / "deformation-execution.json", outcome)
    return int(any(codes))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["optimize", "deform"])
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--python", type=Path, required=True)
    a = p.parse_args()
    raise SystemExit(
        execute(
            a.repo.resolve(),
            a.evidence.resolve(),
            a.work.resolve(),
            a.python.resolve(),
            a.stage,
        )
    )
