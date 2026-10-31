"""J63 counted HMDS Psi4 execution with distinct numerical/collection outcomes."""

import argparse
import fcntl
import importlib
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write

from .campaign import ENGINE_BOHR_ANGSTROM, verify_completed
from .control_engine import configure, molecule

VERSION = "island.psi4-hmds-j63.v1.1"


def read(p):
    return json.loads(Path(p).read_text())


def request(plan, xyz, **changes):
    s = {**plan["protocol"], **changes}
    settings = dict(s["settings"])
    settings.update(
        basis=s["basis"],
        dft_radial_points=s["grid"][0],
        dft_spherical_points=s["grid"][1],
        e_convergence=s["thresholds"][0],
        d_convergence=s["thresholds"][1],
    )
    return {
        "schema": "island.psi4-job.v1",
        "operation": "gradient",
        "method": "pbe0",
        "basis": s["basis"],
        "grid": s["grid"],
        "thresholds": s["thresholds"],
        "scf_type": "df",
        "threads": 4,
        "settings": settings,
        "elements": plan["elements"],
        "atom_order": plan["atom_order"],
        "graph_identity": plan["graph_identity"],
        "charge": 0,
        "multiplicity": 1,
        "coordinates": np.asarray(xyz).tolist(),
    }


def reserve(work, key, r):
    with (work / "ledger.lock").open("a+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        p = work / "ledger.json"
        a = read(p) if p.exists() else {"started_unix": time.time(), "calls": {}}
        require(key not in a["calls"], "additive_retry_identity_required")
        require(
            len(a["calls"]) < 80 and time.time() - a["started_unix"] < 36000,
            "J63_campaign_budget_exhausted",
        )
        a["calls"][key] = {
            "request_identity": digest(r),
            "started_unix": time.time(),
            "status": "reserved_or_interrupted",
        }
        write(p, a)


def capture(work, key, r, callback):
    import psi4

    reserve(work, key, r)
    out = work / "outputs" / key
    out.mkdir(parents=True)
    write(work / "requests" / (key + ".json"), r)
    write(out / "input.json", r)
    psi4.core.set_output_file(str(out / "psi4.out"), False)
    start = time.monotonic()
    result = {
        "adapter": VERSION,
        "adapter_sha256": pin(__file__),
        "job_identity": digest(r),
        "input_sha256": pin(out / "input.json"),
        "engine": "Psi4",
        "version": psi4.__version__,
        "operation": "gradient",
        "status": "failed",
    }
    try:
        g, wfn = callback()
        np.save(out / "gradient.npy", g.np, allow_pickle=False)
        g.save(str(out / "gradient.dat"), False, False)
        mol = wfn.molecule()
        result.update(
            status="completed",
            energy_hartree=wfn.energy(),
            geometry_bohr=mol.geometry().np.tolist(),
            elements=[mol.symbol(i) for i in range(mol.natom())],
            charge=int(mol.molecular_charge()),
            multiplicity=mol.multiplicity(),
            functional=wfn.functional().name(),
            basis=wfn.basisset().name(),
            units={
                "geometry": "bohr",
                "energy": "hartree",
                "gradient": "hartree/bohr",
                "hessian": "hartree/bohr^2",
            },
            axis_order="atom-major Cartesian xyz",
            scf_converged=True,
            scalar_variables=wfn.scalar_variables(),
            basis_functions=wfn.basisset().nbf(),
            basis_shells=wfn.basisset().nshell(),
        )
        return g, wfn
    except BaseException as exc:
        result.update(error_type=type(exc).__name__, error=str(exc))
        (out / "traceback.txt").write_text(traceback.format_exc())
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - start
        psi4.core.flush_outfile()
        result["native_files"] = {
            p.name: pin(p) for p in sorted(out.iterdir()) if p.is_file()
        }
        write(out / "result.json", result)
        with (work / "ledger.lock").open("a+") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            a = read(work / "ledger.json")
            a["calls"][key].update(
                status=result["status"],
                elapsed_seconds=result["elapsed_seconds"],
                result_sha256=pin(out / "result.json"),
            )
            write(work / "ledger.json", a)


def gradient(work, key):
    import psi4

    r = read(work / "requests" / (key + ".json"))
    out = work / "outputs" / key
    if out.exists():
        return verify_completed(out, digest(r))
    configure(psi4)
    psi4.set_options(r["settings"])
    mol = molecule(psi4, r)
    try:
        capture(
            work, key, r, lambda: psi4.gradient("pbe0", molecule=mol, return_wfn=True)
        )
    finally:
        psi4.core.clean()


def optimization_exit(status):
    return {"converged": 0, "budget_exhausted": 3, "unexpected_failure": 1}[status]


def optimization_budget_exception(exc):
    """Require a specific cap signal; unrelated failures at callback40 stay failures."""
    if getattr(exc, "reason", None) == "J63_optimization_callback_cap":
        return True
    if (
        type(exc).__name__ != "OptimizationConvergenceError"
        or type(exc).__module__ != "psi4.driver.p4util.exceptions"
    ):
        return False
    trace = exc.__traceback__
    while trace is not None:
        frame = trace.tb_frame
        if (
            frame.f_globals.get("__name__") == "psi4.driver.driver"
            and frame.f_code.co_name == "optimize"
        ):
            context = frame.f_locals
            data = context.get("opt_data")
            error = data.get("error") if isinstance(data, dict) else None
            message = error.get("error_message") if isinstance(error, dict) else None
            return message == "Maximum number of steps exceeded: 40." or (
                context.get("n", 0) > 40 and context.get("opt_status") != "FAILED"
            )
        trace = trace.tb_next
    return False


def optimize(work):
    import optking
    import psi4

    plan = read(work / "scientific-plan.json")
    out = work / "optimization"
    require(not out.exists(), "optimization_attempt_exists")
    out.mkdir()
    configure(psi4)
    psi4.core.set_output_file(str(out / "psi4.out"), False)
    r = request(plan, plan["initial_coordinates"])
    settings = {**r["settings"], **plan["optimizer"]["settings"]}
    psi4.set_options(settings)
    mol = molecule(psi4, r)
    driver = importlib.import_module("psi4.driver.driver")
    original = driver.gradient
    calls = []

    def recorded(method, **kwargs):
        require(len(calls) < 40, "J63_optimization_callback_cap")
        rr = request(plan, kwargs["molecule"].geometry().np * ENGINE_BOHR_ANGSTROM)
        key = f"optimization/step-{len(calls):03d}"
        try:
            return capture(work, key, rr, lambda: original(method, **kwargs))
        finally:
            calls.append(key)
            write(out / "calls.json", calls)
            psi4.core.set_output_file(str(out / "psi4.out"), True)

    driver.gradient = recorded
    write(
        out / "instrumentation.json",
        {
            "settings": settings,
            "optking_version": optking.__version__,
            "driver_sha256": pin(driver.__file__),
            "adapter_sha256": pin(__file__),
            "max_callbacks": 40,
            "implicit_QM_Hessian": False,
        },
    )
    start = time.monotonic()
    result = {"status": "unexpected_failure", "calls": calls}
    try:
        e, _, history = psi4.optimize(
            "pbe0", molecule=mol, return_wfn=True, return_history=True
        )
        result.update(status="converged", energy_hartree=e)
        write(
            out / "history.json",
            {k: [np.asarray(v).tolist() for v in vs] for k, vs in history.items()},
        )
    except BaseException as exc:  # noqa: BLE001 -- numerical failure retained and returned nonzero
        exhausted = (
            optimization_budget_exception(exc)
            and len(calls) >= 40
            and all(
                read(work / "outputs" / k / "result.json")["status"] == "completed"
                for k in calls
            )
        )
        result.update(
            status="budget_exhausted" if exhausted else "unexpected_failure",
            error_type=type(exc).__name__,
            error=str(exc),
        )
        (out / "traceback.txt").write_text(traceback.format_exc())
    finally:
        driver.gradient = original
        psi4.core.flush_outfile()
        result.update(
            calls=calls,
            elapsed_seconds=time.monotonic() - start,
            numerical_exit_code=optimization_exit(result["status"]),
        )
        result["native_files"] = {
            p.name: pin(p) for p in sorted(out.iterdir()) if p.is_file()
        }
        write(out / "result.json", signed(result))
        psi4.core.clean()
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=["optimize", "gradient"])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--key")
    a = p.parse_args()
    w = a.work.resolve()
    if a.operation == "optimize":
        sys.exit(optimize(w)["numerical_exit_code"])
    gradient(w, a.key)
