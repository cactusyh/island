"""J60 immutable bounded reference campaign; reuses the native J59 importer."""

import argparse
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from island.pcff_completion import digest, pin, require, write

from .campaign import verify_completed
from .hamiltonian_engine import OPTIONS


def requests(plan, plan_hash):
    jobs = {}
    for name, points in plan["protocol_points"].items():
        for point in points:
            protocol = plan["protocols"][name]
            job = {
                "schema": "island.psi4-job.v1",
                **protocol,
                "scf_type": "df",
                "operation": "gradient",
                "threads": 4,
                "elements": plan["elements"],
                "atom_order": plan["atom_order"],
                "graph_identity": plan["graph_identity"],
                "charge": 0,
                "multiplicity": 1,
                "coordinates": plan["geometries"][point],
                "plan_sha256": plan_hash,
            }
            job["settings"] = {
                **OPTIONS,
                "basis": job["basis"],
                "scf_type": "df",
                "df_basis_scf": "def2-universal-jkfit",
                "dft_radial_points": job["grid"][0],
                "dft_spherical_points": job["grid"][1],
                "e_convergence": job["thresholds"][0],
                "d_convergence": job["thresholds"][1],
            }
            jobs[name + "/" + point] = job
    require(len(jobs) <= plan["budget"]["maximum_gradient_jobs"], "job_budget_exceeded")
    return jobs


def execute(evidence, work, python):
    evidence, work = Path(evidence).resolve(), Path(work).resolve()
    plan_path = evidence / "plan/reference-plan.json"
    plan = json.loads(plan_path.read_text())
    require(
        plan["identity"] == digest({k: v for k, v in plan.items() if k != "identity"}),
        "plan_mutation",
    )
    jobs = requests(plan, pin(plan_path))
    adapter_module = plan.get(
        "adapter_module", "island.pcff_qm_pilot.hamiltonian_engine"
    )
    require(
        adapter_module
        in (
            "island.pcff_qm_pilot.hamiltonian_engine",
            "island.pcff_qm_pilot.hamiltonian_refinement_engine",
        ),
        "unsupported_adapter",
    )
    adapter = Path(__file__).with_name(adapter_module.rsplit(".", 1)[1] + ".py")
    manifest = {
        "plan_sha256": pin(plan_path),
        "adapter_sha256": pin(adapter),
        "jobs": {k: digest(v) for k, v in jobs.items()},
        "budget": plan["budget"],
    }
    mp = evidence / "plan/campaign.json"
    if mp.exists():
        require(json.loads(mp.read_text()) == manifest, "campaign_mutation")
    else:
        write(mp, manifest)
    started = time.monotonic()

    def one(item):
        key, job = item
        request = work / "requests" / (key + ".json")
        output = work / "outputs" / key
        if request.exists():
            require(json.loads(request.read_text()) == job, "request_mutation")
        else:
            write(request, job)
        if output.exists():
            verify_completed(output, digest(job))
            return {"job": key, "status": "verified_resume_skip"}
        require(
            time.monotonic() - started < plan["budget"]["campaign_wall_seconds"],
            "campaign_wall_budget",
        )
        scratch = work / "scratch" / key
        scratch.mkdir(parents=True, exist_ok=True)
        log = work / "logs" / (key + ".log")
        log.parent.mkdir(parents=True, exist_ok=True)
        require(
            not log.exists() and not (work / "attempts" / (key + ".json")).exists(),
            "prior_attempt_requires_additive_retry",
        )
        env = {
            **os.environ,
            "PSI_SCRATCH": str(scratch),
            "OMP_NUM_THREADS": "4",
            "OPENBLAS_NUM_THREADS": "1",
            "PYTHONPATH": str(Path.cwd()) + ":" + str(Path.cwd() / "src"),
        }
        with log.open("w") as stream:
            try:
                p = subprocess.run(
                    [
                        str(Path(python).resolve()),
                        "-m",
                        adapter_module,
                        str(request),
                        str(output),
                    ],
                    env=env,
                    cwd=scratch,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    timeout=plan["budget"]["wall_seconds_per_job"],
                    check=False,
                )
                code = p.returncode
            except subprocess.TimeoutExpired:
                code = 124
        row = {
            "job": key,
            "status": "completed" if code == 0 else "failed_preserved",
            "returncode": code,
        }
        if code == 0:
            verify_completed(output, digest(job))
        print(json.dumps(row), flush=True)
        write(work / "attempts" / (key + ".json"), row)
        return row

    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(one, jobs.items()))
    write(work / "execution-summary.json", rows)
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j60"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--python", type=Path, required=True)
    a = p.parse_args()
    execute(a.evidence, a.work, a.python)
