"""Additive Psi4 molecular importer using the existing hash and matrix contracts."""

import re
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, require
from island.pcff_qm.reference import CONVERSION

from .campaign import ENGINE_BOHR_ANGSTROM, verify_completed
from .importer import UNITS, matrix_ascii


def native(directory, expected):
    directory = Path(directory)
    r = verify_completed(directory, digest(expected))
    n = len(expected["elements"])
    require(
        expected["atom_order"] == list(range(1, n + 1)) and n in (9, 27),
        "graph_bound_atom_order_required",
    )
    require(
        r["units"] == UNITS and r["axis_order"] == "atom-major Cartesian xyz",
        "native_unit_or_axis_mutation",
    )
    require(
        r["engine"] == "Psi4"
        and r["version"] == "1.9.1"
        and r["scf_converged"] is True,
        "native_engine_or_convergence",
    )
    require(
        r["charge"] == expected["charge"] == 0
        and r["multiplicity"] == expected["multiplicity"] == 1,
        "native_state_mutation",
    )
    require(
        [e.upper() for e in r["elements"]] == [e.upper() for e in expected["elements"]],
        "native_mapping_mutation",
    )
    require(
        r["functional"].lower() == expected["method"]
        and r["basis"].lower() == expected["basis"],
        "native_protocol_mutation",
    )
    x = np.array(r["geometry_bohr"]) * ENGINE_BOHR_ANGSTROM
    require(
        x.shape == (n, 3)
        and np.isfinite(x).all()
        and np.allclose(x, expected["coordinates"], atol=1e-9, rtol=0),
        "native_geometry_mutation",
    )
    out = (directory / "psi4.out").read_text()
    for pattern, wanted in [
        (r"==> Primary Basis <==\s*Basis Set:\s*(\S+)", expected["basis"]),
        (r"Composite Functional:\s*(\S+)", expected["method"]),
    ]:
        values = re.findall(pattern, out)
        require(
            values and all(v.lower() == wanted for v in values),
            "printed_protocol_mutation",
        )
    require(
        "Energy and wave function converged." in out
        and "DEF2-UNIVERSAL-JKFIT" in out.upper(),
        "printed_convergence_or_auxiliary_missing",
    )
    require(
        f"SCF Algorithm Type is {expected['scf_type'].upper()}." in out
        and f"{expected['threads']} Threads" in out,
        "printed_algorithm_mutation",
    )
    require(
        f"charge = {expected['charge']}, multiplicity = {expected['multiplicity']}"
        in out,
        "printed_state_mutation",
    )
    for label, value in [
        ("Radial Points", expected["grid"][0]),
        ("Spherical Points", expected["grid"][1]),
        ("Energy threshold", expected["thresholds"][0]),
        ("Density threshold", expected["thresholds"][1]),
        ("Integral threshold", 1e-12),
    ]:
        matches = re.findall(label + r"\s*=\s*([\d.Ee+-]+)", out)
        require(
            matches and all(float(v) == value for v in matches),
            "printed_setting_mutation",
        )
    rows = re.findall(
        r"^\s+(SI|C|O|H)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s+[-+\d.]+\s*$",
        out,
        re.MULTILINE,
    )[:n]
    require(
        [v[0] for v in rows] == [v.upper() for v in expected["elements"]],
        "printed_mapping_mutation",
    )
    printed = np.array([[float(v) for v in row[1:]] for row in rows])
    unit = re.search(r"Geometry \(in (Angstrom|Bohr)\)", out).group(1)
    if unit == "Bohr":
        printed *= ENGINE_BOHR_ANGSTROM
    require(
        np.allclose(printed, expected["coordinates"], atol=1e-9, rtol=0),
        "printed_geometry_mutation",
    )
    energy = float(re.findall(r"@(?:DF-)?RKS Final Energy:\s*([-+\d.]+)", out)[-1])
    require(
        np.isfinite(r["energy_hartree"]) and abs(energy - r["energy_hartree"]) < 2e-10,
        "native_energy_mutation",
    )
    g = np.load(directory / "gradient.npy", allow_pickle=False)
    require(
        g.shape == (n, 3)
        and np.isfinite(g).all()
        and np.allclose(
            g, matrix_ascii(directory / "gradient.dat", (n, 3)), atol=1e-15, rtol=0
        ),
        "native_gradient_mutation",
    )
    return r, g


def normalized(directory, request):
    r, g = native(directory, request)
    return {
        "coordinates": request["coordinates"],
        "atom_order": request["atom_order"],
        "elements": request["elements"],
        "graph_identity": request["graph_identity"],
        "charge": 0,
        "multiplicity": 1,
        "protocol": {
            k: request[k]
            for k in ("method", "basis", "grid", "thresholds", "scf_type", "settings")
        },
        "energy": r["energy_hartree"] * CONVERSION["energy_factor"],
        "forces": (
            -g * CONVERSION["energy_factor"] / CONVERSION["distance_factor"]
        ).tolist(),
        "native_job_identity": digest(request),
    }
