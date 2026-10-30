"""Versioned Psi4 1.9.1 API runner; no other QM adapter is executed.

Each invocation is one immutable job. Scratch is external; outputs survive errors.
"""

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, write

VERSION = "island.psi4-native-pilot.v3"
OPTIONS = {
    "basis": "def2-svp",
    "scf_type": "pk",
    "reference": "rhf",
    "e_convergence": 1e-10,
    "d_convergence": 1e-8,
    "maxiter": 200,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
    "ints_tolerance": 1e-12,
    "fail_on_maxiter": True,
    "geom_maxiter": 30,
    "g_convergence": "gau_tight",
    "findif__fd_project": False,
    "findif__points": 3,
    "findif__disp_size": 0.005,
}


def run(request_path, directory):
    import psi4

    request_path, directory = Path(request_path), Path(directory)
    job = json.loads(request_path.read_text())
    require(job["schema"] == "island.psi4-job.v1", "job_schema_mutation")
    expected_settings = {**OPTIONS, "basis": job["basis"]}
    if job.get("restart_parent_result_sha256"):
        require(job["operation"] == "optimize", "restart_operation_mutation")
        expected_settings.update(geom_maxiter=10, opt_coordinates="cartesian")
    require(job["settings"] == expected_settings, "settings_mutation")
    require(job["method"] in ("pbe0", "b3lyp"), "unsupported_pilot_method")
    require(job["basis"] in ("def2-svp", "def2-tzvp"), "unsupported_pilot_basis")
    xyz = np.asarray(job["coordinates"], dtype=float)
    require(
        xyz.shape == (len(job["elements"]), 3) and np.isfinite(xyz).all(),
        "invalid_geometry",
    )
    require(job["charge"] == 0 and job["multiplicity"] == 1, "pilot_state_mutation")
    require(not directory.exists(), "job_output_already_exists")
    directory.mkdir(parents=True)
    write(directory / "input.json", job)
    psi4.core.set_output_file(str(directory / "psi4.out"), False)
    psi4.set_num_threads(2)
    psi4.set_memory("2 GiB")
    scratch = Path(os.environ["PSI_SCRATCH"])
    require(scratch.is_dir(), "external_scratch_required")
    psi4.core.IOManager.shared_object().set_default_path(str(scratch))
    body = "\n".join(
        e + " " + " ".join(format(v, ".17g") for v in row)
        for e, row in zip(job["elements"], xyz, strict=True)
    )
    molecule = psi4.geometry(
        "0 1\n" + body + "\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n"
    )
    psi4.set_options(job["settings"])
    start = time.monotonic()
    result = {
        "adapter": VERSION,
        "adapter_sha256": pin(__file__),
        "job_identity": digest(job),
        "input_sha256": pin(directory / "input.json"),
        "engine": "Psi4",
        "version": psi4.__version__,
        "status": "failed",
        "operation": job["operation"],
    }
    try:
        op = job["operation"]
        if op == "optimize":
            _, wfn = psi4.optimize(job["method"], molecule=molecule, return_wfn=True)
        elif op == "hessian":
            # Psi4 selects its supported derivative route and prints displacements.
            hessian, wfn = psi4.hessian(
                job["method"], molecule=molecule, return_wfn=True
            )
            np.save(directory / "hessian.npy", hessian.np, allow_pickle=False)
            hessian.save(str(directory / "hessian.dat"), False, False)
        elif op == "gradient":
            _, wfn = psi4.gradient(job["method"], molecule=molecule, return_wfn=True)
        elif op == "energy":
            _, wfn = psi4.energy(job["method"], molecule=molecule, return_wfn=True)
        else:
            raise ValueError("unsupported operation")
        gradient = wfn.gradient()
        if gradient is not None:
            np.save(directory / "gradient.npy", gradient.np, allow_pickle=False)
            gradient.save(str(directory / "gradient.dat"), False, False)
        mol = wfn.molecule()
        result.update(
            {
                "status": "completed",
                "energy_hartree": wfn.energy(),
                "geometry_bohr": mol.geometry().np.tolist(),
                "elements": [mol.symbol(i) for i in range(mol.natom())],
                "charge": int(mol.molecular_charge()),
                "multiplicity": mol.multiplicity(),
                "functional": wfn.functional().name(),
                "basis": wfn.basisset().name(),
                "scalar_variables": wfn.scalar_variables(),
                "units": {
                    "geometry": "bohr",
                    "energy": "hartree",
                    "gradient": "hartree/bohr",
                    "hessian": "hartree/bohr^2",
                },
                "axis_order": "atom-major Cartesian xyz",
                "scf_converged": True,
                "convergence_evidence": "Psi4 returned with FAIL_ON_MAXITER=True; native output retained",
            }
        )
    except Exception as exc:  # noqa: BLE001 -- preserve any native engine failure
        result.update({"error_type": type(exc).__name__, "error": str(exc)})
    finally:
        result["elapsed_seconds"] = time.monotonic() - start
        psi4.core.flush_outfile()
        result["native_files"] = {
            p.name: pin(p) for p in sorted(directory.iterdir()) if p.is_file()
        }
        write(directory / "result.json", result)
        psi4.core.clean()
    require(result["status"] == "completed", "qm_calculation_failed", result=result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.request, args.output)
