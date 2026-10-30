"""Strict importer for executed Psi4 1.9.1 native API exports (no pickle).

Native text/Matrix.save evidence and lossless arrays are kept separately from
normalized ISLAND samples. Imported values never authorize a fitted term.
"""

import json
import re
import shutil
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, inventory, pin, require, write
from island.pcff_qm.reference import CONVERSION, intake

from .campaign import ENGINE_BOHR_ANGSTROM, NAMES, selected, verify_completed

VERSION = "island.psi4-native-import.v2"
UNITS = {
    "geometry": "bohr",
    "energy": "hartree",
    "gradient": "hartree/bohr",
    "hessian": "hartree/bohr^2",
}


def matrix_ascii(path, shape):
    """Psi4 Matrix.save full indexed ASCII; native rows/columns are explicit."""
    lines = Path(path).read_text().splitlines()
    require(
        lines[1] == "symmetry 0" and int(lines[2]) == int(np.prod(shape)),
        "native_matrix_header",
    )
    values, seen = np.zeros(shape), set()
    for line in lines[3:]:
        if not line.strip():
            continue
        i, j, v = line.split()
        key = int(i), int(j)
        require(
            key not in seen and 0 <= key[0] < shape[0] and 0 <= key[1] < shape[1],
            "native_matrix_index",
        )
        values[key] = float(v)
        seen.add(key)
    require(
        len(seen) == int(np.prod(shape)) and np.isfinite(values).all(),
        "native_matrix_incomplete",
    )
    return values


def native(directory, expected):
    directory = Path(directory)
    r = verify_completed(directory, digest(expected))
    require(
        r["engine"] == "Psi4" and r["version"] == "1.9.1", "unsupported_native_engine"
    )
    require(r.get("scf_converged") is True, "unconverged_scf")
    require(r["units"] == UNITS, "native_unit_mutation")
    require(
        r["axis_order"] == "atom-major Cartesian xyz", "hessian_axis_order_mutation"
    )
    require(
        [e.upper() for e in r["elements"]] == [e.upper() for e in expected["elements"]],
        "native_atom_order_mutation",
    )
    require(
        r["charge"] == expected["charge"]
        and r["multiplicity"] == expected["multiplicity"],
        "native_state_mutation",
    )
    require(r["basis"].lower() == expected["basis"], "native_basis_mutation")
    require(r["functional"].lower() == expected["method"], "native_method_mutation")
    xyz = np.asarray(r["geometry_bohr"], dtype=float) * ENGINE_BOHR_ANGSTROM
    require(
        np.isfinite(xyz).all()
        and np.allclose(xyz, expected["coordinates"], rtol=0, atol=1e-9),
        "native_geometry_or_order_mutation",
    )
    out = (directory / "psi4.out").read_text()
    require(
        "Energy and wave function converged." in out, "native_scf_convergence_missing"
    )
    geometry_rows = re.findall(
        r"^\s+(SI|O|H)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s+[-+\d.]+\s*$",
        out,
        re.MULTILINE,
    )[: len(expected["elements"])]
    require(
        len(geometry_rows) == len(expected["elements"])
        and [r[0] for r in geometry_rows] == [e.upper() for e in expected["elements"]],
        "printed_atom_order_mutation",
    )
    printed_xyz = np.array([[float(v) for v in row[1:]] for row in geometry_rows])
    first_geometry_unit = re.search(r"Geometry \(in (Angstrom|Bohr)\)", out).group(1)
    if first_geometry_unit == "Bohr":
        printed_xyz *= ENGINE_BOHR_ANGSTROM
    require(
        np.allclose(
            printed_xyz,
            expected["coordinates"],
            atol=1e-9,
            rtol=0,
        ),
        "printed_native_geometry_mutation",
    )
    require(
        f"charge = {expected['charge']}, multiplicity = {expected['multiplicity']}"
        in out,
        "printed_native_state_mutation",
    )
    for label, value in (("Radial Points", 99), ("Spherical Points", 590)):
        matches = re.findall(label + r"\s*=\s*(\d+)", out)
        require(
            matches and all(int(m) == value for m in matches), "native_grid_mutation"
        )
    for label, value in (
        ("Energy threshold", 1e-10),
        ("Density threshold", 1e-8),
        ("Integral threshold", 1e-12),
    ):
        matches = re.findall(label + r"\s*=\s*([\d.Ee+-]+)", out)
        require(
            matches and all(float(m) <= value for m in matches),
            "native_scf_threshold_mutation",
        )
    require(
        "SCF Algorithm Type is PK." in out and "2 Threads" in out,
        "native_execution_settings_mutation",
    )
    energies = re.findall(r"@RKS Final Energy:\s*([-+\d.]+)", out)
    require(
        energies
        and min(
            abs(float(v) - r["energy_hartree"])
            for v in (energies if expected["operation"] == "hessian" else energies[-1:])
        )
        < 2e-10,
        "native_energy_disagreement",
    )
    require(np.isfinite(r["energy_hartree"]), "nonfinite_native_energy")
    arrays = {}
    n = len(expected["atom_order"])
    for key, shape in (("gradient", (n, 3)), ("hessian", (3 * n, 3 * n))):
        path = directory / (key + ".npy")
        if path.exists():
            a = np.load(path, allow_pickle=False)
            require(
                a.shape == shape and np.isfinite(a).all(),
                "native_derivative_shape_or_value",
                channel=key,
            )
            if key == "hessian":
                require(
                    np.allclose(a, a.T, atol=1e-8, rtol=0),
                    "nonsymmetric_native_hessian",
                )
            require(
                np.allclose(
                    a,
                    matrix_ascii(directory / (key + ".dat"), shape),
                    rtol=0,
                    atol=1e-15,
                ),
                "native_matrix_export_disagreement",
                channel=key,
            )
            arrays[key] = a
    require(
        expected["operation"] != "gradient" or "gradient" in arrays,
        "native_gradient_missing",
    )
    require(
        expected["operation"] != "hessian" or "hessian" in arrays,
        "native_hessian_missing",
    )
    return r, arrays


def normalize(repo, distribution, evidence, work, output):
    evidence, work, output = Path(evidence), Path(work), Path(output)
    require(not output.exists(), "normalized_dataset_already_exists")
    manifest = json.loads((evidence / "plan/campaign.json").read_text())
    reports = {}
    for name in NAMES:
        scan = json.loads((evidence / "plan" / (name + ".json")).read_text())
        system, assessment, candidates = selected(
            repo, distribution, scan["center_coordinates"]
        )
        target = next(
            t for t in candidates if t.payload["identity"] == scan["target"]["identity"]
        )
        t = target.payload
        root = output / name
        root.mkdir(parents=True)
        (root / "LICENSE.txt").write_text(
            "ISLAND-generated numerical QM pilot data and inputs: CC0-1.0. Psi4 native logs are retained as computation evidence; upstream software/documentation retains its own license.\n"
        )
        shutil.copyfile(__file__, root / "importer.py")
        rows, proofs = [], []
        for point in scan["points"]:
            key = name + "/" + point["id"]
            request = json.loads((work / "requests" / (key + ".json")).read_text())
            require(digest(request) == manifest["jobs"][key], "campaign_job_mutation")
            src = work / "outputs" / key
            raw, arrays = native(src, request)
            dest = root / "native" / point["id"]
            shutil.copytree(src, dest)
            sample = {
                "id": point["id"],
                "graph_identity": t["graph_identity"],
                "graph": t["graph"],
                "atom_order": t["atom_order"],
                "charge_state": t["charge_state"],
                "spin_state": t["spin_state"],
                "coordinates": point["coordinates"],
                "offsets": point["offsets"],
                "split": point["split"],
                "method_id": "pbe0-svp",
                "energy": raw["energy_hartree"] * CONVERSION["energy_factor"],
                "forces": (
                    -arrays["gradient"]
                    * CONVERSION["energy_factor"]
                    / CONVERSION["distance_factor"]
                ).tolist(),
            }
            raw_values = {
                "units": {
                    "energy": "hartree",
                    "gradient": "hartree/bohr",
                    "hessian": "hartree/bohr^2",
                },
                "energy": raw["energy_hartree"],
                "gradient": arrays["gradient"].tolist(),
            }
            if all(v == 0 for v in point["offsets"]):
                hkey = name + "/hessian"
                hrequest = json.loads(
                    (work / "requests" / (hkey + ".json")).read_text()
                )
                _, harrays = native(work / "outputs" / hkey, hrequest)
                shutil.copytree(work / "outputs" / hkey, root / "native/hessian")
                raw_values["hessian"] = harrays["hessian"].tolist()
                sample["hessian"] = (
                    harrays["hessian"]
                    * CONVERSION["energy_factor"]
                    / CONVERSION["distance_factor"] ** 2
                ).tolist()
            rows.append(sample)
            files = {
                p.relative_to(root).as_posix(): pin(p)
                for p in dest.iterdir()
                if p.is_file()
            }
            if "hessian" in sample:
                files.update(
                    {
                        p.relative_to(root).as_posix(): pin(p)
                        for p in (root / "native/hessian").iterdir()
                        if p.is_file()
                    }
                )
            proofs.append(
                {
                    "sample_id": sample["id"],
                    "output_files": files,
                    "normalized_sample_checksum": digest(sample),
                    "status": "verified_by_supplied_importer",
                    "raw_observables": raw_values,
                }
            )
        write(root / "samples.json", rows)
        proof = {
            "importer_file": "importer.py",
            "importer_sha256": pin(root / "importer.py"),
            "implementation_version": VERSION,
            "normalization": {
                "energy": "hartree to kcal/mol by declared factor",
                "forces": "minus gradient; hartree/bohr to kcal/mol/angstrom",
                "conversion_source": CONVERSION["source"],
                "energy_factor": CONVERSION["energy_factor"],
                "distance_factor": CONVERSION["distance_factor"],
            },
            "records": proofs,
            "engine_coordinate_conversion": {
                "bohr_to_angstrom": ENGINE_BOHR_ANGSTROM,
                "explanation": "engine CODATA-2014 geometry constant; J57 observable normalization retains declared CODATA-2022 factors; relative difference 2.38e-10 is recorded, never silently changed",
            },
        }
        write(root / "verification.json", proof)
        req = {
            "schema": "island.pcff-qm-reference.v2",
            "family": "pcff",
            "target_identity": t["identity"],
            **{
                k: t[k]
                for k in (
                    "graph_identity",
                    "typing_identity",
                    "charge_identity",
                    "atom_order",
                    "charge_vector",
                    "source_hashes",
                    "domain",
                )
            },
            "graph": t["graph"],
            "functional_form_identity": t["functional_form"]["identity"],
            "units": t["units"],
            "files": inventory(root),
            "license_file": "LICENSE.txt",
            "permission_status": "licensed",
            "methods": {
                "pbe0-svp": {
                    "method": "PBE0",
                    "basis_set": "def2-SVP",
                    "dispersion_treatment": "none",
                    "software": "Psi4",
                    "version": "1.9.1",
                    "calculation_settings": t["qm_protocol"],
                }
            },
            "verification_file": "verification.json",
            "samples_file": "samples.json",
            "conversion_source": CONVERSION["source"],
            "energy_factor": CONVERSION["energy_factor"],
            "distance_factor": CONVERSION["distance_factor"],
        }
        result = intake(
            system,
            assessment,
            target,
            req,
            root,
            verification_sha256=pin(root / "verification.json"),
        )
        # The recursive data root is frozen. Requests/results live outside it.
        write(output / (name + "-request.json"), req)
        write(output / (name + "-intake.json"), result.payload)
        reports[name] = {
            "status": "qm_data_available",
            "sample_count": len(rows),
            "intake_identity": result.payload["identity"],
            "target_identity": t["identity"],
        }
    return reports
