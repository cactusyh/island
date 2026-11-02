"""Execute additive, zero-QM J68 span diagnostics with retained solver attempts."""

import argparse
import os
import resource
import subprocess
import sys
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_augmentation import solver_state_diagnostics
from island.pcff_qm_pilot.hmds_constrained import (
    dual_lower_bound,
    physical_problem,
    save_system,
    solve,
)
from island.pcff_qm_pilot.hmds_fit import features
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.hmds_span import (
    authorize_attempt,
    certificate,
    check,
    classification,
    intake,
    projection,
    system,
)
from scripts.analyze_phase_4j66 import full_evaluator_check
from scripts.analyze_phase_4j67 import report_predictions

HYPOTHESES = "40dc0193e635883c7e86ccfc04f61cce3333c04da6b1cb34b37917509b68085c"

ROOT = Path("docs/evidence/phase_4j68")
FROZEN = "147b5785f64fd859cdacff207065c484b442da9f5b77a4efc0c3b26f0f7d0775"


def reconstruct(repo, evidence=None):
    parent = repo / "docs/evidence/phase_4j67"
    root = repo / ROOT if evidence is None else evidence
    plan = read(root / "plan/scientific-plan.json")
    require(pin(root / "plan/scientific-plan.json") == FROZEN, "J68_frozen_plan")
    require(
        all(pin(parent / n) == sha for n, sha in plan["parent_files"].items()),
        "J68_parent_pins",
    )
    rows = read(parent / "fitting-dataset.json")["observations"]
    model = load_model(repo)
    data = features(model, rows)
    with np.load(parent / "physical-arrays.npz", allow_pickle=False) as saved:
        require(
            saved["names"].tolist() == data["names"]
            and all(
                np.array_equal(saved[k], v) for k, v in data.items() if k != "names"
            ),
            "J68_matrix_reconstruction",
        )
    spec = read(parent / "plan/scientific-plan.json")["parameter_specifications"]
    require([s["name"] for s in spec] == data["names"], "J68_column_order")
    with np.load(parent / "physical-arrays.npz", allow_pickle=False) as saved:
        original_arrays = {
            k: (saved[k].tolist() if k == "names" else saved[k]) for k in saved.files
        }
    intake(
        model,
        rows,
        data,
        spec,
        read(parent / "fitting-dataset.json")["observations"],
        original_arrays,
        read(parent / "plan/scientific-plan.json")["parameter_specifications"],
    )
    return (
        plan,
        model,
        rows,
        data,
        np.array([s["lower"] for s in spec]),
        np.array([s["upper"] for s in spec]),
    )


def complete_check(model, rows, data, values, variant=None):
    if variant is None:
        return full_evaluator_check(model, rows, data, values)
    from island.pcff_qm_pilot.hmds_baseline import evaluate
    from island.pcff_qm_pilot.hmds_local_baseline import jet_correction

    energy, forces = [], []
    for row in rows:
        original = evaluate(
            model,
            row["coordinates"],
            dict(zip(data["names"][:58], values[:58], strict=True)),
        )
        e, f = jet_correction(model, row["coordinates"], variant)
        energy.append(original["energy"] + e @ values[58:60] + values[-1])
        forces.append(original["forces"] + f @ values[58:60])
    expected_e = data["fixed_energy"] + data["energy"] @ values[:-1] + values[-1]
    expected_f = data["fixed_force"] + data["force"] @ values[:-1]
    de, df = np.array(energy) - data["QM_energy"], np.array(forces) - data["QM_force"]
    delta_e = float(np.max(abs(np.array(energy) - expected_e)))
    delta_f = float(np.max(abs(np.array(forces) - expected_f)))
    require(max(delta_e, delta_f) < 1e-9, "J68_complete_amended_MM_disagreement")
    return {
        "maximum_energy_array_difference": delta_e,
        "maximum_force_array_difference": delta_f,
        "all60_explicitly_supplied": True,
        "energy_RMSE": float(np.sqrt(np.mean(de**2))),
        "force_RMS": float(np.sqrt(np.mean(df**2))),
        "maximum_absolute_energy_error": float(np.max(abs(de))),
        "maximum_local_force_error": float(np.max(np.linalg.norm(df, axis=2))),
        "original_four_limits_satisfied": bool(
            np.sqrt(np.mean(de**2)) <= 0.1
            and np.sqrt(np.mean(df**2)) <= 0.5
            and np.max(abs(de)) <= 0.3
            and np.max(np.linalg.norm(df, axis=2)) <= 2
        ),
    }


def prediction_report(data, values, ids):
    if len(data["names"]) == 58:
        return report_predictions(data, values, ids)
    groups = {
        "protected22": np.arange(22),
        "new4": np.arange(22, 26),
        "full26": np.arange(26),
    }
    return {
        key: check(
            physical_problem(data, index), None, values, np.array([]), np.array([])
        )["subsets"]["full26"]
        for key, index in groups.items()
    }


def worker(repo, name, backend, bounded, recover=False):
    resource.setrlimit(resource.RLIMIT_AS, (4294967296, 4294967296))
    authorize_attempt(
        read(repo / ROOT / "solver-attempts.json")["attempts"],
        name,
        backend,
        bounded,
        recover,
    )
    plan, model, rows, data, lo, hi = reconstruct(repo)
    variant = next(
        (v for v in ("bridge_angle", "SiO_stretch") if name.startswith(v + "-")), None
    )
    if variant:
        from island.pcff_qm_pilot.hmds_local_baseline import extend

        require(
            pin(repo / ROOT / "plan/baseline-hypotheses.json") == HYPOTHESES,
            "J68_frozen_hypotheses",
        )
        data = extend(model, rows, data, variant)
    p, protected = physical_problem(data), physical_problem(data, np.arange(22))
    s = system(p, lo, hi, protected, bounded)
    root = repo / ROOT / "solves"
    if recover:
        from scripts.replay_phase_4j66 import compare_system

        compare_system(root / (name + "-canonical.npz"), s)
        result = read(root / (name + "-raw.json"))
        result = {k: v for k, v in result.items() if k != "identity"}
    else:
        save_system(root / (name + "-canonical.npz"), s)
        result = solve(s, backend, plan[backend + "_settings"])
        write(root / (name + "-raw.json"), signed(result))
    checks = check(p, protected, result["values"], lo, hi, result["rho"])
    cert = (
        dual_lower_bound(
            s,
            result,
            protected
            if checks["subsets"]["protected22"]["required_rho"]
            > checks["subsets"]["full26"]["required_rho"]
            else p,
            lo,
            hi,
        )
        if bounded
        else certificate(s, result, p, protected, lo, hi)
    )
    if bounded:
        cert["certified"] = True
    out = {
        "solution": result,
        "checks": checks,
        "certificate": cert,
        "classification": classification(checks, cert),
        "raw_state_diagnostics": solver_state_diagnostics(s, result),
        "complete_evaluator_check": complete_check(
            model, rows, data, np.array(result["values"]), variant
        ),
        "predictions": prediction_report(
            data, result["values"], [r["id"] for r in rows]
        ),
        "parameter_names": data["names"] + ["energy_offset"],
        "resources": {
            "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        },
    }
    require(not (root / (name + ".json")).exists(), "J68_no_outcome_overwrite")
    write(root / (name + ".json"), signed(out))
    if recover:
        write(
            root / (name + "-postprocessing-recovery.json"),
            signed(
                {
                    "original_failure_log_sha256": pin(root / (name + ".log")),
                    "raw_state_sha256": pin(root / (name + "-raw.json")),
                    "canonical_sha256": pin(root / (name + "-canonical.npz")),
                    "solver_reexecuted": False,
                    "reason": "Historical J65/J67 prediction helper requires exactly58 parameters. J68 amended model needs60; new generic physical-array report plus independent complete scalar-potential evaluation. Original failed process retained.",
                }
            ),
        )
    print(name, result["status"], result["rho"], cert["lower_bound"], flush=True)


def attempt(repo, name, backend, bounded=False):
    root = repo / ROOT
    ledger_path = root / "solver-attempts.json"
    ledger = read(ledger_path)["attempts"] if ledger_path.exists() else []
    require(
        len(ledger) < 12 and all(r["name"] != name for r in ledger),
        "J68_additive_budget",
    )
    command = [
        sys.executable,
        "scripts/analyze_phase_4j68.py",
        "--worker",
        name,
        "--backend",
        backend,
    ] + (["--bounded"] if bounded else [])
    ledger.append(
        {
            "number": len(ledger) + 1,
            "name": name,
            "backend": backend,
            "bounded": bounded,
            "command": command,
            "status": "started",
        }
    )
    write(ledger_path, signed({"attempts": ledger}))
    with (root / "solves" / (name + ".log")).open("x") as log:
        try:
            r = subprocess.run(
                command,
                env=os.environ.copy(),
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=210,
                check=False,
            )
            ledger[-1].update(
                returncode=r.returncode,
                status="completed" if r.returncode == 0 else "failed",
            )
        except subprocess.TimeoutExpired:
            ledger[-1].update(status="interrupted_timeout", returncode=None)
    write(ledger_path, signed({"attempts": ledger}))
    require(ledger[-1]["status"] == "completed", "J68_solve_failed_retained")
    print((root / "solves" / (name + ".log")).read_text(), flush=True)


def main(repo):
    _plan, model, rows, data, lo, hi = reconstruct(repo)
    p, protected = physical_problem(data), physical_problem(data, np.arange(22))
    projected = projection(p)
    projected["checks"] = check(p, protected, projected["values"], lo, hi)
    projected["predictions"] = report_predictions(
        data, projected["values"], [r["id"] for r in rows]
    )
    projected["complete_evaluator_check"] = full_evaluator_check(
        model, rows, data, projected["values"]
    )
    write(repo / ROOT / "model-span-projection.json", signed(projected))
    for name, backend, bounded in [
        ("bounded-Clarabel", "clarabel", True),
        ("bounded-SCS", "scs", True),
        ("bound-free-Clarabel", "clarabel", False),
        ("bound-free-SCS", "scs", False),
    ]:
        attempt(repo, name, backend, bounded)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recover", action="store_true")
    parser.add_argument("--variants", action="store_true")
    parser.add_argument("--worker")
    parser.add_argument("--backend", default="clarabel")
    parser.add_argument("--bounded", action="store_true")
    a = parser.parse_args()
    if a.worker:
        worker(Path.cwd(), a.worker, a.backend, a.bounded, a.recover)
    elif a.variants:
        for variant in ("bridge_angle", "SiO_stretch"):
            for backend in ("clarabel", "scs"):
                attempt(Path.cwd(), variant + "-" + backend, backend)
    else:
        main(Path.cwd())
