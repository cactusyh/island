"""Prepare/resume only the explicitly bounded, hash-bound J58 pilot jobs."""

import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from island.pcff_completion import digest, pin, require, write
from island.pcff_qm.__main__ import load_cases
from island.pcff_qm.geometry import geometries
from island.pcff_qm.targets import targets

from .engine import OPTIONS, VERSION

QUERIES = ("angle-angle:hsi,sio,hsi,hsi", "bond-bond_1_3:sio,osi,sio,hsi")
NAMES = ("angle-angle", "bond-bond")
# Psi4 1.9.1 uses QCElemental's CODATA-2014 molecular length constant.
# This is explicit engine coordinate reconstruction, not a change to J57 units.
ENGINE_BOHR_ANGSTROM = 0.52917721067


def selected(repo, distribution, center=None):
    system, assessment = load_cases(repo, distribution)["disiloxane"]
    if center is not None:
        for i, xyz in zip(assessment.payload["atom_order"], center, strict=True):
            system.coordinates.set(i, xyz)
        from island.pcff_pdms import prepare_polymer_forcefield

        assessment = prepare_polymer_forcefield(
            system, forcefield="pcff", assets=assessment.assets, audit_si_h=True
        )
    chosen = {t.payload["query"]: t for t in targets(system, assessment, "disiloxane")}
    return system, assessment, [chosen[q] for q in QUERIES]


def prepare(repo, distribution, evidence, work):
    """Create a new, additive scan revision; historical J57 records stay intact."""
    repo, evidence, work = Path(repo), Path(evidence), Path(work)
    plan = json.loads((evidence / "plan/pilot-plan.json").read_text())
    opt = json.loads((work / "optimize-restart/result.json").read_text())
    require(
        opt["status"] == "completed" and opt["operation"] == "optimize",
        "optimized_center_required",
    )
    center = [[v * ENGINE_BOHR_ANGSTROM for v in row] for row in opt["geometry_bohr"]]
    system, assessment, chosen = selected(repo, distribution, center)
    scans, jobs = {}, {}
    for index, (name, target) in enumerate(zip(NAMES, chosen, strict=True)):
        t = target.payload
        points = geometries(system, assessment, t["scan"])
        scans[name] = {
            "schema": "island.j58-optimized-scan.v1",
            "parent_J57_target_identity": plan["selected_targets"][index]["identity"],
            "center_native_result_sha256": pin(work / "optimize-restart/result.json"),
            "center_coordinates": center,
            "target": t,
            "points": points,
            "held_out_v2": [
                p["id"]
                for p in points
                if (
                    p["id"]
                    in (
                        "point-0000",
                        "point-0002",
                        "point-0004",
                        "point-0010",
                        "point-0014",
                        "point-0020",
                        "point-0022",
                        "point-0024",
                    )
                )
            ],
            "fit_authorized": False,
        }
        for point in points:
            key = name + "/" + point["id"]
            jobs[key] = job(t, point["coordinates"], "gradient", evidence)
        jobs[name + "/hessian"] = job(t, center, "hessian", evidence)
    sens = [
        ("center", center),
        ("angle-corner", scans["angle-angle"]["points"][0]["coordinates"]),
        ("bond-corner", scans["bond-bond"]["points"][-1]["coordinates"]),
    ]
    for label, xyz in sens:
        for method, basis in (("pbe0", "def2-tzvp"), ("b3lyp", "def2-svp")):
            jobs["sensitivity/" + label + "-" + method + "-" + basis] = job(
                chosen[0].payload, xyz, "gradient", evidence, method=method, basis=basis
            )
    xyz = scans["angle-angle"]["points"][0]["coordinates"]
    for step in (0.001, 0.0005):
        for sign in (-1, 1):
            x = [r[:] for r in xyz]
            x[3][0] += sign * step
            jobs[f"finite-difference/h{step:g}-sign{sign:+d}"] = job(
                chosen[0].payload, x, "energy", evidence
            )
    require(len(jobs) == 62, "pilot_budget_mutation")
    manifest = {
        "schema": "island.j58-campaign.v1",
        "plan_sha256": pin(evidence / "plan/pilot-plan.json"),
        "adapter_version": VERSION,
        "adapter_sha256": pin(Path(__file__).with_name("engine.py")),
        "jobs": {k: digest(v) for k, v in jobs.items()},
        "counts": {
            "scan_gradients": 50,
            "hessians": 2,
            "sensitivity_gradients": 6,
            "finite_difference_energies": 4,
        },
        "maximum_concurrency": 2,
        "wall_seconds_per_job": 900,
        "campaign_wall_seconds": 14400,
        "optimizer_job": "separately retained; max 30 steps",
    }
    if (evidence / "plan/campaign.json").exists():
        require(
            json.loads((evidence / "plan/campaign.json").read_text()) == manifest,
            "campaign_resume_identity_mutation",
        )
    else:
        write(evidence / "plan/campaign.json", manifest)
        for name, scan in scans.items():
            write(evidence / "plan" / (name + ".json"), scan)
    for key, value in jobs.items():
        p = work / "requests" / (key + ".json")
        if p.exists():
            require(json.loads(p.read_text()) == value, "request_resume_mutation")
        else:
            write(p, value)
    return manifest


def job(target, xyz, operation, evidence, *, method="pbe0", basis="def2-svp"):
    return {
        "schema": "island.psi4-job.v1",
        "operation": operation,
        "method": method,
        "basis": basis,
        "settings": {**OPTIONS, "basis": basis},
        "elements": [a["element"] for a in target["graph"]["atoms"]],
        "coordinates": xyz,
        "atom_order": target["atom_order"],
        "graph_identity": target["graph_identity"],
        "charge": 0,
        "multiplicity": 1,
        "plan_sha256": pin(Path(evidence) / "plan/pilot-plan.json"),
    }


def verify_completed(directory, job_identity):
    r = json.loads((directory / "result.json").read_text())
    if r["status"] == "failed" and (directory / "metadata-recovery.json").exists():
        from .recovery import recovered

        r = recovered(directory)
    require(
        r["status"] == "completed" and r["job_identity"] == job_identity,
        "unfinished_or_changed_job",
    )
    require(
        all(pin(directory / name) == sha for name, sha in r["native_files"].items()),
        "native_output_hash_mutation",
    )
    require(
        digest(json.loads((directory / "input.json").read_text())) == job_identity,
        "native_input_mutation",
    )
    return r


def execute(repo, evidence, work, python):
    repo, evidence, work = (
        Path(repo).resolve(),
        Path(evidence).resolve(),
        Path(work).resolve(),
    )
    manifest = json.loads((evidence / "plan/campaign.json").read_text())
    require(
        pin(repo / "src/island/pcff_qm_pilot/engine.py") == manifest["adapter_sha256"],
        "execution_adapter_mutation",
    )
    python = Path(python).resolve()
    start = time.monotonic()

    def one(item):
        key, identity = item
        request, output = work / "requests" / (key + ".json"), work / "outputs" / key
        require(
            digest(json.loads(request.read_text())) == identity, "request_hash_mutation"
        )
        if output.exists():
            # Never overwrite a failed or interrupted job, or silently rerun it.
            verify_completed(output, identity)
            return {"job": key, "status": "verified_resume_skip"}
        require(
            time.monotonic() - start < manifest["campaign_wall_seconds"],
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
            "OMP_NUM_THREADS": "2",
            "OPENBLAS_NUM_THREADS": "1",
        }
        with log.open("w") as handle:
            try:
                p = subprocess.run(
                    [
                        str(python),
                        "-m",
                        "island.pcff_qm_pilot.engine",
                        str(request),
                        str(output),
                    ],
                    cwd=scratch,
                    env=env,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    timeout=manifest["wall_seconds_per_job"],
                    check=False,
                )
                code = p.returncode
            except subprocess.TimeoutExpired:
                code = 124
        status = {
            "job": key,
            "returncode": code,
            "status": "completed" if code == 0 else "failed_preserved",
        }
        if code == 0:
            verify_completed(output, identity)
        print(json.dumps(status), flush=True)
        return status

    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(one, manifest["jobs"].items()))
    write(work / "execution-summary.json", rows)
    require(
        all(r["status"] in ("completed", "verified_resume_skip") for r in rows),
        "campaign_incomplete",
    )
    return rows
