"""Count every J62 Psi4 gradient; retain failures before scratch cleanup."""

import argparse
import fcntl
import importlib
import json
import os
import time
import traceback
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write

from .campaign import ENGINE_BOHR_ANGSTROM, verify_completed
from .hamiltonian_refinement_engine import OPTIONS

VERSION = "island.psi4-j62-molecular-control.v1"


def read(p):
    return json.loads(Path(p).read_text())


def request(spec, xyz, **overrides):
    s = {**spec, **overrides}
    settings = {
        **OPTIONS,
        "basis": s["basis"],
        "scf_type": "df",
        "df_basis_scf": "def2-universal-jkfit",
        "dft_radial_points": s["grid"][0],
        "dft_spherical_points": s["grid"][1],
        "e_convergence": s["thresholds"][0],
        "d_convergence": s["thresholds"][1],
    }
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
        "elements": s["elements"],
        "atom_order": s["atom_order"],
        "graph_identity": s["graph_identity"],
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
            len(a["calls"]) < 48 and time.time() - a["started_unix"] < 28800,
            "campaign_budget_exhausted",
        )
        a["calls"][key] = {
            "request_identity": digest(r),
            "started_unix": time.time(),
            "status": "reserved_or_interrupted",
        }
        write(p, a)


def complete(work, key, r):
    with (work / "ledger.lock").open("a+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        p = work / "ledger.json"
        a = read(p)
        a["calls"][key].update(
            status=r["status"],
            elapsed_seconds=r["elapsed_seconds"],
            result_sha256=pin(work / "outputs" / key / "result.json"),
        )
        write(p, a)


def molecule(psi4, r):
    return psi4.geometry(
        "0 1\n"
        + "\n".join(
            e + " " + " ".join(format(float(v), ".17g") for v in x)
            for e, x in zip(r["elements"], r["coordinates"], strict=True)
        )
        + "\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n"
    )


def configure(psi4):
    psi4.set_num_threads(4)
    psi4.set_memory("4 GiB")
    psi4.core.IOManager.shared_object().set_default_path(os.environ["PSI_SCRATCH"])


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
        complete(work, key, result)


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
            work,
            key,
            r,
            lambda: psi4.gradient(r["method"], molecule=mol, return_wfn=True),
        )
    finally:
        psi4.core.clean()


def optimize(work):
    import optking
    import psi4

    plan = read(work / "scientific-plan.json")
    spec = plan["HMDS"]
    out = work / "hmds-optimization"
    require(not out.exists(), "optimization_attempt_exists")
    out.mkdir()
    configure(psi4)
    psi4.core.set_output_file(str(out / "psi4.out"), False)
    initial = request(spec, spec["initial_coordinates"])
    settings = {
        **initial["settings"],
        "geom_maxiter": 12,
        "g_convergence": "gau_tight",
        "opt_coordinates": "redundant",
        "full_hess_every": -1,
    }
    psi4.set_options(settings)
    mol = molecule(psi4, initial)
    driver = importlib.import_module("psi4.driver.driver")
    original = driver.gradient
    calls = []

    def recorded(method, **kwargs):
        require(len(calls) < 12, "HMDS_optimization_call_cap")
        r = request(spec, kwargs["molecule"].geometry().np * ENGINE_BOHR_ANGSTROM)
        key = f"hmds-opt/step-{len(calls):03d}"
        try:
            return capture(work, key, r, lambda: original(method, **kwargs))
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
            "hidden_QM_Hessians": False,
        },
    )
    start = time.monotonic()
    result = {"status": "unfinished", "calls": calls}
    try:
        e, _, hist = psi4.optimize(
            "pbe0", molecule=mol, return_wfn=True, return_history=True
        )
        result.update(status="completed", energy_hartree=e)
        write(
            out / "history.json",
            {k: [np.asarray(v).tolist() for v in vs] for k, vs in hist.items()},
        )
    except BaseException as exc:  # noqa: BLE001 -- retain interrupted optimizer diagnostics
        result.update(error_type=type(exc).__name__, error=str(exc))
        (out / "traceback.txt").write_text(traceback.format_exc())
    finally:
        driver.gradient = original
        psi4.core.flush_outfile()
        result.update(calls=calls, elapsed_seconds=time.monotonic() - start)
        write(out / "result.json", signed(result))
        psi4.core.clean()
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=["gradient", "optimize"])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--key")
    a = p.parse_args()
    a.work = a.work.resolve()
    optimize(a.work) if a.operation == "optimize" else gradient(a.work, a.key)
