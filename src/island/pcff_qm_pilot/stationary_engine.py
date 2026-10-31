"""J61 counted Psi4 calls; native gradients reuse the verified J59 importer.

Optimization instrumentation is process-local. No installed engine file or
historical adapter is modified. Every attempted gradient consumes a reservation.
"""

import argparse
import fcntl
import json
import os
import signal
import time
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write

from .campaign import ENGINE_BOHR_ANGSTROM, verify_completed
from .hamiltonian_refinement_engine import OPTIONS

VERSION = "island.psi4-j61-counted-stationary.v1"


def read(path):
    return json.loads(Path(path).read_text())


def job(plan_path, coordinates, **protocol):
    plan = read(plan_path)
    p = {**plan["protocol"], **protocol}
    settings = {
        **OPTIONS,
        "basis": p["basis"],
        "scf_type": "df",
        "df_basis_scf": p["df_basis_scf"],
        "dft_radial_points": p["grid"][0],
        "dft_spherical_points": p["grid"][1],
        "e_convergence": p["thresholds"][0],
        "d_convergence": p["thresholds"][1],
    }
    return {
        "schema": "island.psi4-job.v1",
        "method": p["method"],
        "basis": p["basis"],
        "grid": p["grid"],
        "thresholds": p["thresholds"],
        "scf_type": "df",
        "operation": "gradient",
        "threads": 4,
        "elements": plan["elements"],
        "atom_order": plan["atom_order"],
        "graph_identity": plan["graph_identity"],
        "charge": 0,
        "multiplicity": 1,
        "coordinates": np.asarray(coordinates).tolist(),
        "settings": settings,
        "plan_sha256": pin(plan_path),
    }


def reserve(work, key, request, plan):
    """Reserve before entering the engine, including unsuccessful attempts."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    with (work / "ledger.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = work / "ledger.json"
        ledger = (
            read(path) if path.exists() else {"started_unix": time.time(), "calls": {}}
        )
        require(key not in ledger["calls"], "prior_attempt_requires_additive_retry")
        require(len(ledger["calls"]) < 120, "gradient_budget_exhausted")
        require(
            time.time() - ledger["started_unix"]
            < plan["budget"]["campaign_wall_seconds"],
            "campaign_wall_budget_exhausted",
        )
        ledger["calls"][key] = {
            "request_identity": digest(request),
            "started_unix": time.time(),
            "status": "reserved_or_interrupted",
        }
        write(path, ledger)


def finish(work, key, result):
    work = Path(work)
    with (work / "ledger.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = work / "ledger.json"
        ledger = read(path)
        ledger["calls"][key].update(
            status=result["status"],
            elapsed_seconds=result["elapsed_seconds"],
            result_sha256=pin(work / "outputs" / key / "result.json"),
        )
        write(path, ledger)


def gradient(request_path, work, key):
    """One gradient through the historical, unchanged, hash-bound adapter."""
    from .hamiltonian_refinement_engine import run

    work = Path(work)
    request = read(request_path)
    plan = read(work / "scientific-plan.json")
    directory = work / "outputs" / key
    if directory.exists():
        return verify_completed(directory, digest(request))
    reserve(work, key, request, plan)
    try:
        return run(request_path, directory)
    finally:
        if (directory / "result.json").exists():
            finish(work, key, read(directory / "result.json"))


def optimize(work):
    import importlib

    import optking
    import psi4

    work = Path(work).resolve()
    plan_path = work / "scientific-plan.json"
    plan = read(plan_path)
    out = work / "optimization"
    require(not out.exists(), "optimization_attempt_already_exists")
    out.mkdir(parents=True)
    psi4.core.set_output_file(str(out / "psi4.out"), False)
    psi4.set_num_threads(4)
    psi4.set_memory("4 GiB")
    psi4.core.IOManager.shared_object().set_default_path(os.environ["PSI_SCRATCH"])
    initial = job(plan_path, plan["initial_geometry"])
    settings = {
        **initial["settings"],
        "geom_maxiter": 40,
        "opt_coordinates": "redundant",
        "full_hess_every": -1,
        "g_convergence": "gau_tight",
    }
    psi4.set_options(settings)
    body = "\n".join(
        e + " " + " ".join(format(float(v), ".17g") for v in row)
        for e, row in zip(initial["elements"], initial["coordinates"], strict=True)
    )
    molecule = psi4.geometry(
        "0 1\n" + body + "\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n"
    )
    driver = importlib.import_module("psi4.driver.driver")
    original = driver.gradient
    calls = []

    def recorded(method, **kwargs):
        require(len(calls) < 40, "optimization_gradient_budget_exhausted")
        mol = kwargs["molecule"]
        request = job(plan_path, mol.geometry().np * ENGINE_BOHR_ANGSTROM)
        key = f"optimization/step-{len(calls):03d}"
        directory = work / "outputs" / key
        reserve(work, key, request, plan)
        directory.mkdir(parents=True)
        write(work / "requests" / (key + ".json"), request)
        write(directory / "input.json", request)
        psi4.core.set_output_file(str(directory / "psi4.out"), False)
        start = time.monotonic()
        result = {
            "adapter": VERSION,
            "adapter_sha256": pin(__file__),
            "job_identity": digest(request),
            "input_sha256": pin(directory / "input.json"),
            "engine": "Psi4",
            "version": psi4.__version__,
            "operation": "gradient",
            "status": "failed",
        }

        def timeout(_signum, _frame):
            raise TimeoutError("predeclared 1200-second gradient limit")

        previous_alarm = signal.signal(signal.SIGALRM, timeout)
        signal.alarm(1200)
        try:
            g, wfn = original(method, **kwargs)
            np.save(directory / "gradient.npy", g.np, allow_pickle=False)
            g.save(str(directory / "gradient.dat"), False, False)
            result.update(
                status="completed",
                energy_hartree=wfn.energy(),
                geometry_bohr=wfn.molecule().geometry().np.tolist(),
                elements=[mol.symbol(i) for i in range(9)],
                charge=0,
                multiplicity=1,
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
            )
            return g, wfn
        except Exception as exc:
            result.update(error_type=type(exc).__name__, error=str(exc))
            raise
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, previous_alarm)
            result["elapsed_seconds"] = time.monotonic() - start
            psi4.core.flush_outfile()
            result["native_files"] = {
                p.name: pin(p) for p in sorted(directory.iterdir()) if p.is_file()
            }
            write(directory / "result.json", result)
            finish(work, key, result)
            calls.append(key)
            psi4.core.set_output_file(str(out / "psi4.out"), True)
            write(out / "calls.json", calls)

    write(
        out / "instrumentation.json",
        {
            "adapter_sha256": pin(__file__),
            "driver_sha256": pin(driver.__file__),
            "optking_version": optking.__version__,
            "settings": settings,
            "mechanism": "Process-local driver.gradient wrapper; every analytic callback recorded; FULL_HESS_EVERY=-1 prevents hidden QM Hessians.",
        },
    )
    driver.gradient = recorded
    start = time.monotonic()
    result = {"status": "failed"}
    try:
        energy, _, history = psi4.optimize(
            "pbe0", molecule=molecule, return_wfn=True, return_history=True
        )
        result.update(
            status="completed", energy_hartree=energy, last_call=calls[-1], calls=calls
        )
        write(
            out / "history.json",
            {k: [np.asarray(v).tolist() for v in vs] for k, vs in history.items()},
        )
    except Exception as exc:
        result.update(error_type=type(exc).__name__, error=str(exc), calls=calls)
        raise
    finally:
        driver.gradient = original
        result["elapsed_seconds"] = time.monotonic() - start
        psi4.core.flush_outfile()
        write(out / "result.json", signed(result))
        psi4.core.clean()
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=("optimize", "gradient"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--key")
    a = p.parse_args()
    if a.operation == "optimize":
        optimize(a.work)
    else:
        gradient(a.work / "requests" / (a.key + ".json"), a.work, a.key)
