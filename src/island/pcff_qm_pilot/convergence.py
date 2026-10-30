"""Bounded, plan-bound J59 convergence extension of the executed J58 pilot."""

import argparse
import json
import lzma
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from pathlib import Path

from island.pcff_completion import digest, pin, require, write

from .campaign import verify_completed
from .convergence_fast_engine import OPTIONS, VERSION


def prepare(repo, evidence, work):
    repo, evidence, work = Path(repo), Path(evidence), Path(work)
    plan = json.loads((evidence / "plan/convergence-plan.json").read_text())
    require(
        plan["identity"] == digest({k: v for k, v in plan.items() if k != "identity"}),
        "plan_mutation",
    )
    amendment = json.loads((evidence / "plan/df-amendment.json").read_text())
    jobs = {}
    for name, points in plan["protocol_points"].items():
        protocol = plan["protocols"][name]
        for label in points:
            job = {
                "schema": "island.psi4-job.v1",
                "scf_type": "df",
                "operation": "gradient",
                "method": protocol["method"],
                "basis": protocol["basis"],
                "grid": protocol["grid"],
                "thresholds": protocol.get("thresholds", plan["default_thresholds"]),
                "elements": plan["elements"],
                "coordinates": plan["geometries"][label],
                "atom_order": plan["atom_order"],
                "graph_identity": plan["graph_identity"],
                "charge": 0,
                "multiplicity": 1,
                "plan_sha256": pin(evidence / "plan/convergence-plan.json"),
            }
            job["settings"] = {
                **OPTIONS,
                "scf_type": "df",
                "df_basis_scf": "def2-universal-jkfit",
                "basis": job["basis"],
                "dft_radial_points": job["grid"][0],
                "dft_spherical_points": job["grid"][1],
                "e_convergence": job["thresholds"][0],
                "d_convergence": job["thresholds"][1],
            }
            jobs[name + "/" + label] = job
    for label in amendment["PK_control_points"]:
        control = {**jobs["tz/" + label], "scf_type": "pk"}
        control["settings"] = {
            k: v for k, v in control["settings"].items() if k != "df_basis_scf"
        }
        control["settings"]["scf_type"] = "pk"
        jobs["pk-control/" + label] = control
    optimize = {**jobs["grid/center"], "operation": "optimize"}
    optimize["settings"] = {
        **optimize["settings"],
        "geom_maxiter": 20,
        "opt_coordinates": "cartesian",
    }
    jobs["optimize/qz-dense"] = optimize
    manifest = {
        "schema": "island.j59-campaign.v1",
        "plan_sha256": pin(evidence / "plan/convergence-plan.json"),
        "adapter_version": VERSION,
        "adapter_sha256": pin(
            repo / "src/island/pcff_qm_pilot/convergence_fast_engine.py"
        ),
        "jobs": {k: digest(v) for k, v in jobs.items()},
        "budget": plan["budget"],
        "amendment_sha256": pin(evidence / "plan/df-amendment.json"),
    }
    parent = Path(
        json.loads((evidence / "plan/v2-location.json").read_text())["v2_work"]
    )
    for key in ("grid/angle-angle-0000", "grid/angle-angle-0004"):
        original = json.loads((parent / "requests" / (key + ".json")).read_text())
        verify_completed(parent / "outputs" / key, digest(original))
        jobs[key] = original
    jobs["thread-control/angle-angle-0000"] = dict(jobs["grid/angle-angle-0000"])
    for key, job in jobs.items():
        if key not in ("grid/angle-angle-0000", "grid/angle-angle-0004"):
            job["threads"] = 8
    manifest["jobs"] = {k: digest(v) for k, v in jobs.items()}
    manifest["cpu_amendment_sha256"] = pin(evidence / "plan/cpu-amendment.json")
    import shutil

    for key, job in jobs.items():
        request_path = work / "requests" / (key + ".json")
        if request_path.exists():
            require(
                json.loads(request_path.read_text()) == job, "request_resume_mutation"
            )
        else:
            write(request_path, job)
        if (
            key in ("grid/angle-angle-0000", "grid/angle-angle-0004")
            and not (work / "outputs" / key).exists()
        ):
            shutil.copytree(parent / "outputs" / key, work / "outputs" / key)
    path = evidence / "plan/campaign-fast.json"
    if path.exists():
        require(json.loads(path.read_text()) == manifest, "campaign_resume_mutation")
    else:
        write(path, manifest)
    return manifest


def validate_request(key, job, plan):
    """Reject self-consistent request/manifest mutations outside the saved plan."""
    group, label = key.split("/", 1)
    protocol_name = {
        "pk-control": "tz",
        "thread-control": "grid",
        "optimize": "grid",
    }.get(group, group)
    require(protocol_name in plan["protocols"], "undeclared_protocol")
    protocol = plan["protocols"][protocol_name]
    label = "center" if group == "optimize" else label
    require(
        label in plan["geometries"] and job["coordinates"] == plan["geometries"][label],
        "plan_geometry_mutation",
    )
    require(
        job["graph_identity"] == plan["graph_identity"]
        and job["atom_order"] == plan["atom_order"]
        and job["elements"] == plan["elements"],
        "plan_mapping_mutation",
    )
    require(
        job["method"] == protocol["method"]
        and job["basis"] == protocol["basis"]
        and job["grid"] == protocol["grid"]
        and job["thresholds"] == protocol.get("thresholds", plan["default_thresholds"]),
        "plan_protocol_mutation",
    )
    require(job["charge"] == 0 and job["multiplicity"] == 1, "plan_charge_mutation")
    require(
        job["operation"] == ("optimize" if group == "optimize" else "gradient"),
        "plan_operation_mutation",
    )
    require(
        job["scf_type"] == ("pk" if group == "pk-control" else "df"),
        "plan_algorithm_mutation",
    )
    require(
        job.get("threads", 2)
        == (2 if key in ("grid/angle-angle-0000", "grid/angle-angle-0004") else 8),
        "plan_threads_mutation",
    )


def execution_adapter(evidence, work, python, job):
    """Reproduce the two retained v2 jobs with their exact executed adapter."""
    if job.get("threads", 2) == 2:
        archive = (Path(evidence) / "engine/convergence_df_engine.py.xz").read_bytes()
        source = lzma.decompress(archive)
        pins = json.loads(
            (Path(evidence) / "engine/preflight-adapter-hashes.json").read_text()
        )
        manifest = json.loads((Path(evidence) / "plan/campaign-df.json").read_text())
        require(
            sha256(archive).hexdigest() == pins["convergence_df_engine"]
            and sha256(source).hexdigest() == manifest["adapter_sha256"],
            "historical_execution_adapter_mutation",
        )
        path = Path(work) / "execution-adapters/convergence_df_engine_v2.py"
        if path.exists():
            require(path.read_bytes() == source, "execution_adapter_copy_mutation")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(source)
        return [str(Path(python).resolve()), str(path.resolve())]
    require(job["threads"] == 8, "execution_threads_mutation")
    return [
        str(Path(python).resolve()),
        "-m",
        "island.pcff_qm_pilot.convergence_fast_engine",
    ]


def execute(repo, evidence, work, python):
    repo, evidence, work = (
        Path(repo).resolve(),
        Path(evidence).resolve(),
        Path(work).resolve(),
    )
    manifest = json.loads((evidence / "plan/campaign-fast.json").read_text())
    require(
        pin(repo / "src/island/pcff_qm_pilot/convergence_fast_engine.py")
        == manifest["adapter_sha256"],
        "adapter_mutation",
    )
    require(
        pin(evidence / "plan/convergence-plan.json") == manifest["plan_sha256"],
        "plan_mutation",
    )
    plan = json.loads((evidence / "plan/convergence-plan.json").read_text())
    expected_keys = {
        name + "/" + point
        for name, points in plan["protocol_points"].items()
        for point in points
    }
    expected_keys |= {
        "pk-control/" + p for p in ("center", "angle-angle-0000", "bond-bond-0024")
    }
    expected_keys |= {"thread-control/angle-angle-0000", "optimize/qz-dense"}
    require(
        set(manifest["jobs"]) == expected_keys and len(expected_keys) == 40,
        "campaign_job_budget_mutation",
    )
    require(manifest["budget"] == plan["budget"], "campaign_resource_mutation")
    start = time.monotonic()

    def one(item):
        key, identity = item
        request = work / "requests" / (key + ".json")
        output = work / "outputs" / key
        require(
            digest(json.loads(request.read_text())) == identity, "request_hash_mutation"
        )
        job = json.loads(request.read_text())
        validate_request(key, job, plan)
        if output.exists():
            verify_completed(output, identity)
            return {"job": key, "status": "verified_resume_skip"}
        require(
            time.monotonic() - start < manifest["budget"]["campaign_wall_seconds"],
            "campaign_budget_exhausted",
        )
        scratch = work / "scratch" / key
        scratch.mkdir(parents=True, exist_ok=True)
        log = work / "logs" / (key + ".log")
        log.parent.mkdir(parents=True, exist_ok=True)
        env = {
            **os.environ,
            "PYTHONPATH": str(repo) + ":" + str(repo / "src"),
            "PSI_SCRATCH": str(scratch),
            "OMP_NUM_THREADS": str(job.get("threads", 2)),
            "OPENBLAS_NUM_THREADS": "1",
        }
        with log.open("w") as f:
            try:
                p = subprocess.run(
                    execution_adapter(evidence, work, python, job)
                    + [str(request), str(output)],
                    env=env,
                    cwd=scratch,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    timeout=manifest["budget"]["wall_seconds_per_job"],
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
            verify_completed(output, identity)
        print(json.dumps(row), flush=True)
        return row

    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(one, manifest["jobs"].items()))
    write(work / "execution-summary.json", rows)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=("prepare", "execute"))
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j59"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--python", type=Path)
    a = p.parse_args()
    print(
        json.dumps(
            (
                prepare(a.repo, a.evidence, a.work)
                if a.operation == "prepare"
                else execute(a.repo, a.evidence, a.work, a.python)
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
