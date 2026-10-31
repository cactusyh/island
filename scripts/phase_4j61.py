"""Bounded J61 preparation/execution/reconstruction using the existing importer."""

import argparse
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm.geometry import _adjust
from island.pcff_qm.reference import CONVERSION
from island.pcff_qm_pilot.campaign import selected, verify_completed
from island.pcff_qm_pilot.convergence_importer import native
from island.pcff_qm_pilot.stationary_engine import job, read
from island.pcff_qm_pilot.stretch_baseline import stretching_projection


def immutable(path, value):
    path = Path(path)
    if path.exists():
        require(read(path) == value, "immutable_J61_record_changed", path=str(path))
    else:
        write(path, value)


def imported(work, key):
    work = Path(work)
    request = read(work / "requests" / (key + ".json"))
    result, gradient = native(work / "outputs" / key, request)
    return {
        "key": key,
        "coordinates": request["coordinates"],
        "graph_identity": request["graph_identity"],
        "atom_order": request["atom_order"],
        "request_identity": digest(request),
        "native_result_sha256": pin(work / "outputs" / key / "result.json"),
        "energy": result["energy_hartree"] * CONVERSION["energy_factor"],
        "forces": (
            -gradient * CONVERSION["energy_factor"] / CONVERSION["distance_factor"]
        ).tolist(),
        "gradient_hartree_bohr": gradient.tolist(),
        "protocol": {
            k: request[k]
            for k in ("method", "basis", "grid", "thresholds", "scf_type", "settings")
        },
    }


def stationary(work, evidence):
    work, evidence = Path(work), Path(evidence)
    result = read(work / "optimization/result.json")
    require(result["status"] == "completed", "stationary_optimization_not_complete")
    rows = [imported(work, key) for key in result["calls"]]
    center = rows[-1]
    g = np.asarray(center["gradient_hartree_bohr"])
    plan = read(evidence / "plan/scientific-plan.json")
    criteria = plan["optimizer"]
    require(
        np.max(abs(g)) <= criteria["stationary_gradient_max_hartree_bohr"]
        and np.sqrt(np.mean(g * g)) <= criteria["stationary_gradient_RMS_hartree_bohr"],
        "stationary_gradient_gate_failed",
    )
    x = np.array(center["coordinates"])
    bonds = plan["projection"]["bonds"]
    lengths = [float(np.linalg.norm(x[i - 1] - x[j - 1])) for i, j in bonds]
    require(
        all(1.3 < v < 2.0 for v in lengths[:2])
        and all(1.15 < v < 1.8 for v in lengths[2:]),
        "stationary_connectivity_changed",
    )
    report = signed(
        {
            "schema": "island.j61-stationary-candidate.v1",
            "geometry_identity": digest(center["coordinates"]),
            "target_identity": digest(
                {
                    "graph": center["graph_identity"],
                    "geometry": center["coordinates"],
                    "protocol": center["protocol"],
                    "role": "J61 stationary stretch target",
                }
            ),
            "center": center,
            "optimization_steps": rows,
            "gradient_calls": len(rows),
            "gradient_max_hartree_bohr": float(np.max(abs(g))),
            "gradient_RMS_hartree_bohr": float(np.sqrt(np.mean(g * g))),
            "bond_lengths_angstrom": lengths,
            "stationary_candidate": True,
            "all_direction_minimum_verified": False,
            "mapping": "Original ordered Si/O/Si/H/H/H/H/H/H identity retained; no remapping to a symmetry-equivalent atom order.",
        }
    )
    immutable(evidence / "stationary-center.json", report)
    return report


def prepare(repo, work, evidence):
    repo, work, evidence = Path(repo), Path(work), Path(evidence)
    plan_path = evidence / "plan/scientific-plan.json"
    require(
        pin(plan_path) == pin(work / "scientific-plan.json"), "execution_plan_changed"
    )
    plan = read(plan_path)
    center = stationary(work, evidence)["center"]
    x = np.asarray(center["coordinates"])
    projection = stretching_projection(x, plan["projection"]["bonds"])
    q = projection["Q"]
    immutable(
        evidence / "stretching-subspace.json",
        signed(
            {
                "schema": "island.j61-Cartesian-stretch-subspace.v1",
                "center_geometry_identity": digest(x.tolist()),
                "Q": q.tolist(),
                "rigid_basis": projection["rigid"].tolist(),
                "rigid_projector": projection["projector"].tolist(),
                "bond_singular_values": projection["bond_singular_values"].tolist(),
                "rank": 8,
                "rigid_rank": 6,
                "units": "unweighted orthonormal Cartesian directions; displacement in angstrom",
                "orthonormal_error": float(np.max(abs(q.T @ q - np.eye(8)))),
                "rigid_overlap_max": float(np.max(abs(projection["rigid"].T @ q))),
            }
        ),
    )
    jobs = {}
    for h in plan["projection"]["steps_angstrom"]:
        for d in range(8):
            for sign in (-1, 1):
                jobs[f"curvature/h{h:g}-d{d}-s{sign:+d}"] = job(
                    plan_path, x + (sign * h * q[:, d]).reshape(9, 3)
                )
    protocols = {
        "basis": {"basis": "cc-pv(q+d)z"},
        "grid": {"grid": [99, 590]},
        "scf": {"thresholds": [1e-12, 1e-10]},
    }
    for name, protocol in protocols.items():
        jobs[f"sensitivity/{name}-center"] = job(plan_path, x, **protocol)
        for d in plan["projection"]["sensitivity_directions"]:
            for sign in (-1, 1):
                jobs[f"sensitivity/{name}-d{d}-s{sign:+d}"] = job(
                    plan_path, x + (sign * 0.002 * q[:, d]).reshape(9, 3), **protocol
                )
    # Energy/gradient checks at two additional step sizes, explicitly counted.
    for h in (0.001, 0.0005):
        for sign in (-1, 1):
            jobs[f"finite-difference/h{h:g}-s{sign:+d}"] = job(
                plan_path, x + (sign * h * q[:, 0]).reshape(9, 3)
            )
    system, _, _ = selected(
        repo,
        repo.parent
        / "island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files",
        x,
    )
    rules = plan["new_geometry_rules"]
    rng = np.random.default_rng(rules["training_seed"])
    coordinates = [
        ("bond", [1, 4], 0.025),
        ("bond", [1, 2], 0.035),
        ("bond", [2, 3], 0.035),
        ("bond", [3, 7], 0.025),
        ("angle", [4, 1, 2], 0.07),
        ("angle", [4, 1, 5], 0.07),
        ("angle", [1, 2, 3], 0.10),
        ("angle", [2, 3, 7], 0.07),
    ]
    selections = []
    for split in ("training", "validation"):
        for k, torsions in enumerate(rules[split + "_torsion_offsets_radian"]):
            xyz = x.copy()
            offsets = []
            for kind, sites, extent in coordinates:
                offset = float(rng.uniform(-extent, extent))
                _adjust(
                    system,
                    xyz,
                    list(range(1, 10)),
                    {"kind": kind, "sites": sites},
                    offset,
                )
                offsets.append(offset)
            for sites, offset in zip(
                ([4, 1, 2, 3], [1, 2, 3, 7]), torsions, strict=True
            ):
                _adjust(
                    system,
                    xyz,
                    list(range(1, 10)),
                    {"kind": "torsion", "sites": sites},
                    offset,
                )
            key = f"fresh/{split}-{k}"
            jobs[key] = job(plan_path, xyz)
            if split == "validation":
                jobs[f"fresh/basis-validation-{k}"] = job(
                    plan_path, xyz, basis="cc-pv(q+d)z"
                )
            selections.append(
                {
                    "id": key,
                    "split": split,
                    "coordinates": xyz.tolist(),
                    "geometry_identity": digest(xyz.tolist()),
                    "local_offsets": offsets,
                    "torsion_offsets": torsions,
                    "conformational_sampling_claim": "none; predeclared deformed configurations",
                }
            )
    immutable(
        evidence / "plan/fresh-configurations.json",
        signed(
            {
                "selections": selections,
                "generation_rules_sha256": pin(plan_path),
                "center_identity": digest(x.tolist()),
                "QM_values_used_for_selection": False,
            }
        ),
    )
    require(len(jobs) == 63, "postoptimization_budget_mutation")
    for key, request in jobs.items():
        immutable(work / "requests" / (key + ".json"), request)
    immutable(
        evidence / "plan/gradient-jobs.json",
        signed(
            {
                "plan_sha256": pin(plan_path),
                "jobs": {k: digest(v) for k, v in jobs.items()},
                "new_gradient_calls": len(jobs),
                "optimization_calls_counted_separately": True,
            }
        ),
    )
    return jobs


def execute(repo, work, evidence, python, group=None):
    repo, work, evidence = (
        Path(repo).resolve(),
        Path(work).resolve(),
        Path(evidence).resolve(),
    )
    manifest = read(evidence / "plan/gradient-jobs.json")
    require(
        pin(repo / "src/island/pcff_qm_pilot/stationary_engine.py")
        == read(evidence / "plan/execution-engine.json")["adapter_sha256"],
        "executed_adapter_changed",
    )
    plan = read(evidence / "plan/scientific-plan.json")

    def one(item):
        key, identity = item
        if group and not key.startswith(group):
            return {"key": key, "status": "not_selected"}
        request = work / "requests" / (key + ".json")
        require(digest(read(request)) == identity, "request_mutation")
        output = work / "outputs" / key
        if output.exists():
            verify_completed(output, identity)
            return {"key": key, "status": "verified_resume_skip"}
        scratch = work / "scratch" / key
        scratch.mkdir(parents=True, exist_ok=True)
        scratch_bytes = sum(
            p.stat().st_size for p in (work / "scratch").rglob("*") if p.is_file()
        )
        require(scratch_bytes < 20 * 1024**3, "scratch_budget_exhausted")
        log = work / "logs" / (key + ".log")
        require(not log.exists(), "prior_attempt_requires_additive_retry")
        log.parent.mkdir(parents=True, exist_ok=True)
        env = {
            **os.environ,
            "PYTHONPATH": str(repo) + ":" + str(repo / "src"),
            "PSI_SCRATCH": str(scratch),
            "OMP_NUM_THREADS": "4",
            "OPENBLAS_NUM_THREADS": "1",
        }
        start = time.monotonic()
        with log.open("w") as stream:
            try:
                result = subprocess.run(
                    [
                        str(Path(python).resolve()),
                        "-m",
                        "island.pcff_qm_pilot.stationary_engine",
                        "gradient",
                        "--work",
                        str(work),
                        "--key",
                        key,
                    ],
                    cwd=scratch,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    timeout=plan["budget"]["wall_seconds_per_gradient"],
                    check=False,
                )
                code = result.returncode
            except subprocess.TimeoutExpired:
                code = 124
        row = {
            "key": key,
            "returncode": code,
            "elapsed_seconds": time.monotonic() - start,
            "status": "completed" if code == 0 else "failed_preserved",
        }
        write(work / "attempts" / (key + ".json"), row)
        print(json.dumps(row), flush=True)
        if code == 0:
            imported(work, key)
        return row

    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(one, manifest["jobs"].items()))
    write(work / ("execution-" + (group or "all").replace("/", "_") + ".json"), rows)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=("prepare", "execute", "import"))
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j61"))
    p.add_argument("--python", type=Path)
    p.add_argument("--group")
    a = p.parse_args()
    if a.operation == "prepare":
        prepare(a.repo, a.work, a.evidence)
    elif a.operation == "execute":
        execute(a.repo, a.work, a.evidence, a.python, a.group)
    else:
        jobs = read(a.evidence / "plan/gradient-jobs.json")["jobs"]
        result = signed(
            {
                "samples": {key: imported(a.work, key) for key in jobs},
                "units": {"energy": "kcal/mol", "forces": "kcal/mol/angstrom"},
                "no_synthetic_QM": True,
            }
        )
        write(a.evidence / "normalized-reference.json", result)


if __name__ == "__main__":
    main()
