"""Additive recovery of completed Psi4 FINDIF metadata, without rerunning QM.

The original adapter failure/result and all numerical bytes remain unchanged.
Only the exactly identified Wavefunction.functional accessor failure is eligible.
"""

import ast
import json
import re
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed, write

VERSION = "island.psi4-findif-metadata-recovery.v1"
ERROR = "'psi4.core.Wavefunction' object has no attribute 'functional'"


def reconstruction(directory):
    directory = Path(directory)
    old = json.loads((directory / "result.json").read_text())
    require(
        old["status"] == "failed"
        and old["operation"] == "hessian"
        and old["error_type"] == "AttributeError"
        and old["error"] == ERROR,
        "metadata_recovery_not_applicable",
    )
    require(
        all(pin(directory / name) == sha for name, sha in old["native_files"].items()),
        "native_output_hash_mutation",
    )
    request = json.loads((directory / "input.json").read_text())
    require(old["job_identity"] == digest(request), "recovery_input_mutation")
    out = (directory / "psi4.out").read_text()
    blocks = re.findall(r"==> Input QCSchema <==\s*\n-+\n(.*?)\n-+", out, re.DOTALL)
    require(
        len(blocks) == 49
        and "49 gradients passed in, including the reference geometry." in out,
        "unfinished_hessian_cannot_be_recovered",
    )
    schemas = [ast.literal_eval(block) for block in blocks]
    require(
        out.count("Energy and wave function converged.") == 49,
        "unconverged_hessian_cannot_be_recovered",
    )
    for schema in schemas:
        require(
            schema["driver"] == "gradient"
            and schema["model"]
            == {"basis": request["basis"], "method": request["method"]},
            "findif_protocol_mutation",
        )
        molecule = schema["molecule"]
        require(
            molecule["symbols"] == request["elements"]
            and molecule["molecular_charge"] == request["charge"]
            and molecule["molecular_multiplicity"] == request["multiplicity"],
            "findif_atom_or_state_mutation",
        )
    first = schemas[0]["molecule"]
    energies = re.findall(r"@RKS Final Energy:\s*([-+\d.]+)", out)
    require(len(energies) == 49, "hessian_energy_coverage_incomplete")
    geometry = np.array(first["geometry"]).reshape(-1, 3)
    require(
        np.allclose(
            geometry * 0.52917721067, request["coordinates"], rtol=0, atol=1e-9
        ),
        "hessian_reference_geometry_mutation",
    )
    metadata = {
        **old,
        "status": "completed",
        "original_adapter_status": "failed_metadata_export",
        "numerical_calculation_status": "completed",
        "energy_hartree": float(energies[0]),
        "geometry_bohr": geometry.tolist(),
        "elements": first["symbols"],
        "charge": int(first["molecular_charge"]),
        "multiplicity": first["molecular_multiplicity"],
        "basis": schemas[0]["model"]["basis"],
        "functional": schemas[0]["model"]["method"],
        "scf_converged": True,
        "axis_order": "atom-major Cartesian xyz",
        "units": {
            "geometry": "bohr",
            "energy": "hartree",
            "gradient": "hartree/bohr",
            "hessian": "hartree/bohr^2",
        },
        "metadata_recovery": VERSION,
        "hessian_derivative_route": "49 actual gradients including reference; 24 translation-free Cartesian SALCs; rotations unprojected; 3-point 0.005 bohr stencil; full 27x27 Cartesian matrix reconstructed by Psi4",
    }
    return signed(
        {
            "schema": "island.psi4-findif-metadata-recovery.v1",
            "implementation": VERSION,
            "original_result_sha256": pin(directory / "result.json"),
            "native_output_sha256": pin(directory / "psi4.out"),
            "original_failure_preserved": True,
            "QM_rerun": False,
            "field_sources": {
                "energy": "first/reference @RKS Final Energy native output line",
                "geometry_order_state_method_basis": "first/reference Input QCSchema native input block; all 49 blocks independently checked",
                "convergence": "49 native convergence markers and completed FINDIF assembly",
                "derivatives": "unchanged native Matrix.save ASCII and lossless arrays",
            },
            "metadata": metadata,
        }
    )


def recover(directory):
    directory = Path(directory)
    result = reconstruction(directory)
    path = directory / "metadata-recovery.json"
    if path.exists():
        require(json.loads(path.read_text()) == result, "metadata_recovery_mutation")
    else:
        write(path, result)
    return result


def recovered(directory):
    directory = Path(directory)
    actual = json.loads((directory / "metadata-recovery.json").read_text())
    require(actual == reconstruction(directory), "metadata_recovery_mutation")
    return actual["metadata"]
