"""Standalone native-output/unit/FD audit. Does NOT import the production importer.

Use only the numerical raw archive. No engine installation or scratch is needed.
"""

import argparse
import hashlib
import json
import tarfile
import tempfile
from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np


def check(native_root, normalized_root):
    root, normalized = Path(native_root), Path(normalized_root)
    getcontext().prec = 40
    ef = float(
        Decimal("4.3597447222060e-18") * Decimal("6.02214076e23") / Decimal(4184)
    )
    df = float(Decimal("5.29177210544e-11") * Decimal("1e10"))
    geometry_factor = 0.52917721067
    checks = []
    for name in ("angle-angle", "bond-bond"):
        samples = json.loads((normalized / name / "samples.json").read_text())
        for sample in samples:
            p = root / "outputs" / name / sample["id"]
            r = json.loads((p / "result.json").read_text())
            inp = json.loads((p / "input.json").read_text())
            # Independent native ASCII parsing and complete row/column mapping.
            g = np.zeros((9, 3))
            for line in (p / "gradient.dat").read_text().splitlines()[3:]:
                i, j, v = line.split()
                g[int(i), int(j)] = float(v)
            assert np.allclose(
                g, np.load(p / "gradient.npy", allow_pickle=False), rtol=0, atol=1e-15
            )
            assert np.allclose(-g * ef / df, sample["forces"], rtol=0, atol=2e-12)
            assert abs(r["energy_hartree"] * ef - sample["energy"]) < 1e-9
            assert np.allclose(
                np.array(r["geometry_bohr"]) * geometry_factor,
                sample["coordinates"],
                atol=1e-9,
                rtol=0,
            )
            assert (
                inp["atom_order"] == sample["atom_order"]
                and inp["charge"] == sample["charge_state"]
            )
            assert r["multiplicity"] == 1 and r["scf_converged"]
            if sample.get("hessian") is not None:
                hpath = root / "outputs" / name / "hessian"
                h = np.zeros((27, 27))
                for line in (hpath / "hessian.dat").read_text().splitlines()[3:]:
                    i, j, v = line.split()
                    h[int(i), int(j)] = float(v)
                assert np.allclose(
                    h,
                    np.load(hpath / "hessian.npy", allow_pickle=False),
                    rtol=0,
                    atol=1e-15,
                )
                assert np.allclose(
                    h * ef / df**2, sample["hessian"], rtol=0, atol=3e-12
                )
            checks.append(
                {"scan": name, "sample": sample["id"], "native_to_normalized": "passed"}
            )
    fd = []
    analytic = (
        -np.load(
            root / "outputs/angle-angle/point-0000/gradient.npy", allow_pickle=False
        )[3, 0]
        * ef
        / df
    )
    for h in (0.001, 0.0005):
        values = [
            json.loads(
                (
                    root / f"outputs/finite-difference/h{h:g}-sign{sign:+d}/result.json"
                ).read_text()
            )["energy_hartree"]
            for sign in (1, -1)
        ]
        force = -(values[0] - values[1]) * ef / (2 * h)
        assert abs(force - analytic) <= 0.01
        fd.append(
            {
                "step_angstrom": h,
                "positive_energy_hartree": values[0],
                "negative_energy_hartree": values[1],
                "numerical_force": force,
                "analytical_force": float(analytic),
                "absolute_error": float(abs(force - analytic)),
            }
        )
    return {
        "schema": "island.j58-independent-unit-audit.v1",
        "production_importer_used": False,
        "samples": checks,
        "finite_difference": fd,
        "energy_factor": ef,
        "length_factor": df,
        "engine_geometry_factor": geometry_factor,
        "factor_difference_relative": geometry_factor / df - 1,
        "scope": "independent native indexed ASCII reconstruction of energies, all force vectors and central Hessians; not scientific parameter validation",
    }


def extract(path, dest):
    with tarfile.open(path, "r:xz") as archive:
        assert all(
            i.isfile() and (dest / i.name).resolve().is_relative_to(dest.resolve())
            for i in archive.getmembers()
        )
        archive.extractall(dest, filter="data")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--json", type=Path, required=True)
    a = p.parse_args()
    with tempfile.TemporaryDirectory(prefix="j58-independent-") as temp:
        temp = Path(temp)
        extract(a.evidence / "native-outputs.tar.xz", temp / "native")
        extract(a.evidence / "reference-datasets.tar.xz", temp / "normalized")
        result = check(temp / "native", temp / "normalized")
    result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    a.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
