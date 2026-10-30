"""Settings-aware extension using the existing native matrices and unit contract.

Imports only the executed Psi4 gradient adapter. Convergence subsets are not
complete J57 fitting scans and never bypass the fitting intake coverage gate.
"""

import json
import lzma
import re
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, require, signed
from island.pcff_qm.reference import CONVERSION

from .campaign import ENGINE_BOHR_ANGSTROM, verify_completed
from .importer import UNITS, matrix_ascii


def native(directory, expected):
    directory = Path(directory)
    r = verify_completed(directory, digest(expected))
    require(
        r["engine"] == "Psi4" and r["version"] == "1.9.1", "unsupported_native_engine"
    )
    require(
        r["units"] == UNITS and r["axis_order"] == "atom-major Cartesian xyz",
        "native_unit_or_axis_mutation",
    )
    require(
        expected["operation"] == "gradient" and r["operation"] == "gradient",
        "gradient_operation_required",
    )
    require(r["scf_converged"] is True, "unconverged_scf")
    require(
        [e.upper() for e in r["elements"]] == [e.upper() for e in expected["elements"]]
        and expected["atom_order"] == list(range(1, 10)),
        "native_mapping_mutation",
    )
    require(
        r["charge"] == expected["charge"] == 0
        and r["multiplicity"] == expected["multiplicity"] == 1,
        "native_state_mutation",
    )
    require(
        r["functional"].lower() == expected["method"]
        and r["basis"].lower() == expected["basis"],
        "native_method_basis_mutation",
    )
    xyz = np.asarray(r["geometry_bohr"]) * ENGINE_BOHR_ANGSTROM
    require(
        xyz.shape == (9, 3)
        and np.isfinite(xyz).all()
        and np.allclose(xyz, expected["coordinates"], rtol=0, atol=1e-9),
        "native_geometry_mutation",
    )
    out = (directory / "psi4.out").read_text()
    # Analytic gradient calls print the actual SCF model, not an Input QCSchema
    # block (the optimizer does print those). Bind native text independently of
    # the adapter's metadata; do not fabricate absent QCSchema evidence.
    bases = re.findall(r"==> Primary Basis <==\s*Basis Set:\s*(\S+)", out)
    functionals = re.findall(r"Composite Functional:\s*(\S+)", out)
    require(
        bases
        and all(b.lower() == expected["basis"] for b in bases)
        and functionals
        and all(f.lower() == expected["method"] for f in functionals),
        "native_input_protocol_mutation",
    )
    require("Energy and wave function converged." in out, "native_convergence_missing")
    for label, value in zip(
        ("Radial Points", "Spherical Points"), expected["grid"], strict=True
    ):
        matches = re.findall(label + r"\s*=\s*(\d+)", out)
        require(
            matches and all(int(m) == value for m in matches), "native_grid_mutation"
        )
    for label, value in (
        ("Energy threshold", expected["thresholds"][0]),
        ("Density threshold", expected["thresholds"][1]),
        ("Integral threshold", 1e-12),
    ):
        matches = re.findall(label + r"\s*=\s*([\d.Ee+-]+)", out)
        require(
            matches and all(float(m) <= value for m in matches),
            "native_threshold_mutation",
        )
    require(
        f"SCF Algorithm Type is {expected['scf_type'].upper()}." in out
        and f"{expected.get('threads', 2)} Threads" in out,
        "native_algorithm_mutation",
    )
    if expected["scf_type"] == "df":
        require(
            "DEF2-UNIVERSAL-JKFIT" in out.upper(), "native_auxiliary_basis_mutation"
        )
    require(
        f"charge = {expected['charge']}, multiplicity = {expected['multiplicity']}"
        in out,
        "printed_native_state_mutation",
    )
    geometry_rows = re.findall(
        r"^\s+(SI|O|H)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s+[-+\d.]+\s*$",
        out,
        re.MULTILINE,
    )[:9]
    require(
        [v[0] for v in geometry_rows] == [v.upper() for v in expected["elements"]],
        "printed_mapping_mutation",
    )
    printed = np.array([[float(v) for v in row[1:]] for row in geometry_rows])
    unit = re.search(r"Geometry \(in (Angstrom|Bohr)\)", out).group(1)
    if unit == "Bohr":
        printed *= ENGINE_BOHR_ANGSTROM
    require(
        np.allclose(printed, expected["coordinates"], atol=1e-9, rtol=0),
        "printed_geometry_mutation",
    )
    energies = re.findall(r"@(?:DF-)?RKS Final Energy:\s*([-+\d.]+)", out)
    require(
        energies
        and np.isfinite(r["energy_hartree"])
        and abs(float(energies[-1]) - r["energy_hartree"]) < 2e-10,
        "native_energy_mutation",
    )
    gradient = np.load(directory / "gradient.npy", allow_pickle=False)
    require(
        gradient.shape == (9, 3) and np.isfinite(gradient).all(),
        "native_gradient_mutation",
    )
    require(
        np.allclose(
            gradient,
            matrix_ascii(directory / "gradient.dat", (9, 3)),
            rtol=0,
            atol=1e-15,
        ),
        "native_matrix_mutation",
    )
    return r, gradient


def normalize(evidence, work):
    evidence, work = Path(evidence), Path(work)
    manifest = json.loads((evidence / "plan/campaign-fast.json").read_text())
    rows = {}
    plan = json.loads((evidence / "plan/convergence-plan.json").read_text())
    from .convergence import validate_request

    for key, identity in manifest["jobs"].items():
        if key.startswith("optimize/"):
            continue
        expected = json.loads((work / "requests" / (key + ".json")).read_text())
        require(digest(expected) == identity, "campaign_mapping_mutation")
        validate_request(key, expected, plan)
        r, g = native(work / "outputs" / key, expected)
        if expected.get("threads", 2) == 2:
            adapter_bytes = lzma.decompress(
                (evidence / "engine/convergence_df_engine.py.xz").read_bytes()
            )
            adapter_version = "island.psi4-convergence-df.v2"
        else:
            adapter_bytes = (
                Path(__file__).with_name("convergence_fast_engine.py").read_bytes()
            )
            adapter_version = "island.psi4-convergence-df.v3"
        require(
            r["adapter"] == adapter_version
            and r["adapter_sha256"] == sha256(adapter_bytes).hexdigest(),
            "executed_adapter_mutation",
        )
        rows[key] = {
            "engine": {
                "name": r["engine"],
                "version": r["version"],
                "adapter": r["adapter"],
                "adapter_sha256": r["adapter_sha256"],
            },
            "scf_converged": r["scf_converged"],
            "native_array_precision": str(g.dtype),
            "graph_identity": expected["graph_identity"],
            "atom_order": expected["atom_order"],
            "charge": 0,
            "multiplicity": 1,
            "coordinates": expected["coordinates"],
            "threads": expected.get("threads", 2),
            "protocol": {
                k: expected[k]
                for k in (
                    "method",
                    "basis",
                    "grid",
                    "thresholds",
                    "scf_type",
                    "settings",
                )
            },
            "energy": r["energy_hartree"] * CONVERSION["energy_factor"],
            "forces": (
                -g * CONVERSION["energy_factor"] / CONVERSION["distance_factor"]
            ).tolist(),
            "energy_hartree": r["energy_hartree"],
            "gradient_hartree_bohr": g.tolist(),
            "native_job_identity": identity,
            "native_files": r["native_files"],
            "elapsed_seconds": r["elapsed_seconds"],
        }
    return signed(
        {
            "schema": "island.j59-normalized-convergence.v1",
            "conversion": CONVERSION,
            "coordinate_bohr_angstrom": ENGINE_BOHR_ANGSTROM,
            "units": {
                "energy": "kcal/mol",
                "force": "kcal/mol/angstrom",
                "coordinates": "angstrom",
                "raw_energy": "hartree",
                "raw_gradient": "hartree/bohr",
            },
            "coverage": "selected convergence geometries only; not a complete fitting dataset; strict J57 coverage gate remains active",
            "samples": rows,
        }
    )
