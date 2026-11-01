"""Check actual optimizer status and freeze graph-controlled HMDS path requests."""

import argparse
import json
import re
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write
from island.pcff_qm_pilot.hmds_engine import request
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_paths import affected, leakage, path
from island.pcff_qm_pilot.molecular_native import normalized


def read(p):
    return json.loads(Path(p).read_text())


def native_convergence(text):
    rows = re.findall(
        r"^\s+\d+\s+-\d+\.\d+\s+[-+\d.eE]+[^\n]*~?\s*$", text, re.MULTILINE
    )
    table = []
    for row in rows:
        values = re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", row)
        if len(values) == 7:
            table.append(
                dict(
                    zip(
                        (
                            "iteration",
                            "energy_hartree",
                            "delta_energy",
                            "maximum_internal_force",
                            "RMS_internal_force",
                            "maximum_internal_displacement",
                            "RMS_internal_displacement",
                        ),
                        map(float, values),
                        strict=True,
                    )
                )
            )
    require(table, "native_optimizer_convergence_table_missing")
    return {
        "history": table,
        "final": table[-1],
        "native_internal_threshold_checks": {
            k: table[-1][k] <= limit
            for k, limit in (
                ("maximum_internal_force", 1.5e-5),
                ("RMS_internal_force", 1e-5),
                ("maximum_internal_displacement", 6e-5),
                ("RMS_internal_displacement", 4e-5),
            )
        },
        "units_scope": "Optking redundant internal-coordinate dq/fq; bond displacements in bohr, angular displacements in radian. Not a Cartesian displacement test. Native GAU_TIGHT permits its documented flat-potential alternative; final Cartesian gradients are independently gated below.",
    }


def optimizer_outcome(work):
    """A killed process cannot manufacture a native successful terminal record."""
    process = read(work / "processes/optimization/result.json")
    terminal = work / "optimization/result.json"
    if terminal.exists():
        result = read(terminal)
        require(
            result["identity"]
            == digest({k: v for k, v in result.items() if k != "identity"}),
            "optimizer_status_identity",
        )
        return result, process
    require(
        process["returncode"] != 0, "optimizer_success_without_native_terminal_record"
    )
    calls = sorted(
        k for k in read(work / "ledger.json")["calls"] if k.startswith("optimization/")
    )
    return signed(
        {
            "status": "unexpected_failure",
            "calls": calls,
            "native_terminal_record_supplied": False,
            "collection_status": "partial_archive_retained",
            "numerical_exit_code": process["returncode"],
            "reason": "Process exited or timed out without native terminal optimizer record; no numerical success inferred from archive collection.",
        }
    ), process


def stationary(repo, evidence, work):
    result, process = optimizer_outcome(work)
    samples = {
        k: normalized(work / "outputs" / k, read(work / "requests" / (k + ".json")))
        for k in result["calls"]
        if (work / "outputs" / k / "result.json").exists()
        and read(work / "outputs" / k / "result.json")["status"] == "completed"
    }
    require(samples, "no_successful_optimizer_gradient_available")
    last = next(k for k in result["calls"][::-1] if k in samples)
    s = samples[last]
    raw = read(work / "outputs" / last / "result.json")
    g = np.load(work / "outputs" / last / "gradient.npy", allow_pickle=False)
    plan = read(evidence / "plan/scientific-plan.json")
    criteria = plan["optimizer"]
    maximum = float(abs(g).max())
    rms = float(np.sqrt(np.mean(g * g)))
    graph = load_model(repo)
    x = np.array(s["coordinates"])
    distances = []
    for i, j, _ in graph["bonds"]:
        types = sorted([graph["typing"][str(i)], graph["typing"][str(j)]])
        d = np.linalg.norm(x[i - 1] - x[j - 1])
        distances.append(float(d))
        low, high = (
            (0.8, 1.4)
            if "hc" in types
            else (1.3, 2.2)
            if "osi" in types
            else (1.5, 2.3)
        )
        require(low < d < high, "HMDS_stationary_graph_geometry_changed", sites=[i, j])
    pair = np.linalg.norm(x[:, None] - x[None, :], axis=2)
    np.fill_diagonal(pair, np.inf)
    require(np.min(pair) > 0.7, "HMDS_geometry_clash")
    gates = {
        "native_optimizer_converged": result["status"] == "converged",
        "numerical_process_exit_success": process["returncode"] == 0,
        "maximum_gradient": maximum <= criteria["stationary_gradient_max_hartree_bohr"],
        "RMS_gradient": rms <= criteria["stationary_gradient_RMS_hartree_bohr"],
    }
    r = signed(
        {
            "schema": "island.j63-HMDS-stationary-center.v1",
            "center": s,
            "new_target_identity": digest(
                {
                    "graph": s["graph_identity"],
                    "geometry": s["coordinates"],
                    "protocol": s["protocol"],
                    "phase": "J63-TZVP-stationary",
                }
            ),
            "geometry_identity": digest(s["coordinates"]),
            "optimizer_status": result["status"],
            "numerical_process_returncode": process["returncode"],
            "successful_gradients": len(samples),
            "attempted_callbacks": len(result["calls"]),
            "gradient_max_hartree_bohr": maximum,
            "gradient_RMS_hartree_bohr": rms,
            "stationary_gates": gates,
            "stationary_candidate": all(gates.values()),
            "minimum_verified": False,
            "basis_functions": raw["basis_functions"],
            "basis_shells": raw["basis_shells"],
            "native_optimizer_sha256": pin(work / "optimization/result.json")
            if (work / "optimization/result.json").exists()
            else None,
            "optimizer_outcome": result,
            "native_optimizer_log_sha256": pin(work / "optimization/psi4.out"),
            "criterion_scope": "Native GAU_TIGHT displacement/force convergence plus independently checked final ordered gradient; no molecular minimum claim.",
            "native_convergence_table": native_convergence(
                (work / "optimization/psi4.out").read_text()
            )
            if "Convergence Check" in (work / "optimization/psi4.out").read_text()
            else {"status": "not_available", "history": []},
        }
    )
    write(evidence / "stationary-center.json", r)
    return r


def deformation_jobs(repo, evidence, work):
    center = stationary(repo, evidence, work)
    require(center["stationary_candidate"], "HMDS_stationary_stage_prerequisite_failed")
    plan = read(evidence / "plan/scientific-plan.json")
    model = load_model(repo)
    x = np.array(center["center"]["coordinates"])
    jobs = {}
    paths = {}
    for spec in plan["paths"]:
        for h in spec["steps"]:
            for sign in (1, -1):
                key = f"paths/{spec['name']}-h{h:g}-s{sign:+d}"
                p = path(model, x, spec, sign * h)
                changes = affected(model, x, p["coordinates"])
                if spec["kind"] == "rotation":
                    require(
                        changes
                        and all(
                            r["family"] == "torsion_3"
                            and set(r["sites"][1:3]) == set(spec["bond"])
                            and abs(abs(r["change"]) - h) < 1e-9
                            for r in changes
                        ),
                        "rotation_uncontrolled_internal_changes",
                    )
                else:
                    require(
                        len(changes) == 1
                        and changes[0]["family"] == "quartic_bond"
                        and abs(changes[0]["change"] - sign * h) < 1e-9,
                        "stretch_uncontrolled_internal_changes",
                    )
                paths[key] = {
                    k: v.tolist() if isinstance(v, np.ndarray) else v
                    for k, v in p.items()
                }
                paths[key]["affected_internal_coordinates"] = changes
                jobs[key] = request(plan, p["coordinates"])
    for name, spec in plan["controls"].items():
        changes = {k: spec[k] for k in ("grid", "basis", "thresholds") if k in spec}
        jobs[f"controls/{name}-center"] = request(plan, x, **changes)
        for target in spec["paths"]:
            p = next(s for s in plan["paths"] if s["name"] == target)
            h = min(p["steps"])
            for sign in (1, -1):
                jobs[f"controls/{name}-{target}-s{sign:+d}"] = request(
                    plan, path(model, x, p, sign * h)["coordinates"], **changes
                )
    # Reserve new configurations before inspecting any new path/control QM results.
    rng = np.random.default_rng(4063)
    fresh = []
    for i in range(3):
        v = x + rng.normal(0, 0.015, x.shape)
        for name, angle in [
            ("siloxane", [-0.37, 0.29, 0.48][i]),
            ("methyl", [0.23, -0.31, 0.41][i]),
        ]:
            spec = next(s for s in plan["paths"] if s["name"] == name)
            v = path(model, v, spec, angle)["coordinates"]
        fresh.append(v)
    development = [x, *[np.array(r["coordinates"]) for r in jobs.values()]]
    old = read(repo / "docs/evidence/phase_4j62/normalized-reference.json")["samples"]
    development += [
        np.array(v["coordinates"]) for k, v in old.items() if k.startswith("hmds")
    ]
    development += [
        np.array(read(work / "requests" / (k + ".json"))["coordinates"])
        for k in read(work / "optimization/result.json")["calls"]
    ]
    check = leakage(model, development, fresh)
    reservation = signed(
        {
            "seed": 4063,
            "coordinates": [v.tolist() for v in fresh],
            "geometry_identities": [digest(v.tolist()) for v in fresh],
            "graph_identity": model["graph_identity"],
            "atom_order": model["atom_order"],
            "screening": check,
            "QM_outputs_supplied": False,
            "validation_used_for_weights_or_bounds": False,
            "scope": "new deformed configurations, not independently sampled conformational basins",
        }
    )
    write(evidence / "plan/reserved-validation.json", reservation)
    write(
        evidence / "path-geometries.json",
        signed(
            {
                "center_identity": center["identity"],
                "paths": paths,
                "source_graph_identity": model["graph_identity"],
                "all_bond_angle_checks_passed": True,
                "torsions_are_graph_fragment_rotations": True,
            }
        ),
    )
    for key, r in jobs.items():
        file = work / "requests" / (key + ".json")
        if file.exists():
            require(read(file) == r, "HMDS_resume_request_mutation")
        else:
            write(file, r)
    write(
        evidence / "plan/deformation-jobs.json",
        signed(
            {
                "jobs": {k: digest(r) for k, r in jobs.items()},
                "maximum_calls": 35,
                "activation": "native optimizer convergence and independently checked final gradient",
                "center_identity": center["identity"],
            }
        ),
    )
    pins = {
        "src/island/pcff_qm_pilot/hmds_paths.py": pin(
            repo / "src/island/pcff_qm_pilot/hmds_paths.py"
        ),
        "scripts/prepare_phase_4j63.py": pin(__file__),
    }
    f = evidence / "plan/path-implementation-pins.json"
    if f.exists():
        require(read(f)["files"] == pins, "path_implementation_changed_after_execution")
    else:
        write(f, signed({"files": pins, "before_deformation_execution": True}))
    return list(jobs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j63"))
    p.add_argument("--work", type=Path, required=True)
    a = p.parse_args()
    print(stationary(a.repo, a.evidence, a.work)["stationary_gates"])
