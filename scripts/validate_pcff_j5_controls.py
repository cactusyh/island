"""J5 isolated source-row controls and converter failures; not whole-model acceptance."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from validate_pcff_fallbacks import compare, compiled, lammps, raw_row
from validate_pcff_operational_support import write
from validate_pcff_singlepoint import reference_inputs

from island.workflows.bundle import system_from
from island.workflows.storage import decode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "lammps", "converter", "declaration", "audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        run(args)
    except Exception as error:  # retain diagnostics, then preserve the original error
        write(
            args.output / "failure.json",
            {"error": f"{type(error).__name__}: {error}", "whole_system_gate": False},
        )
        raise


def run(args):
    declaration = json.loads(args.declaration.read_text())
    raw = args.source.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == declaration["source"]["sha256"]
    term_inputs = []
    for term in declaration["term_controls"]:
        row = raw_row(raw, term["family"], term["types"])
        term_inputs.append(
            dict(
                term,
                source=row,
                coordinates_angstrom=[
                    [0.1, 0.2, 0.3],
                    [row["values"][0] + 0.17, 0.31, -0.12],
                ],
            )
        )
    write(
        args.output / "declaration.json",
        {
            "terms": term_inputs,
            "tolerances": declaration["tolerances"],
            "scope": "Isolated raw source quartic bonds, never completed domain models",
            "source_sha256": hashlib.sha256(raw).hexdigest(),
            "lammps_sha256": hashlib.sha256(args.lammps.read_bytes()).hexdigest(),
            "converter_sha256": hashlib.sha256(args.converter.read_bytes()).hexdigest(),
        },
    )
    report = {"terms": [], "converter_controls": [], "whole_system_gate": False}
    for index, item in enumerate(term_inputs):
        p = item["source"]["values"]
        xyz = np.array(item["coordinates_angstrom"])
        coeff = [p[0], *[v * 4.184 for v in p[1:]]]
        energy, force = compiled("quartic_bond", coeff, xyz)
        ref, rf, meta = lammps(
            args.lammps.resolve(),
            args.output / f"term-{index}",
            xyz,
            "bond",
            "class2",
            " ".join(format(v, ".17g") for v in p),
        )
        comparison = compare(energy, force, ref, rf)
        # Independently differentiated radial polynomial, angstrom input.
        vec = xyz[1] - xyz[0]
        r = np.linalg.norm(vec)
        dr = r - p[0]
        exact_e = 4.184 * sum(p[n - 1] * dr**n for n in (2, 3, 4))
        deriv = 4.184 * sum(n * p[n - 1] * dr ** (n - 1) for n in (2, 3, 4))
        exact_f = np.array([deriv * vec / r, -deriv * vec / r])
        analytic = compare(energy, force, exact_e, exact_f)
        fd = []

        def E(x, p=p):
            d = np.linalg.norm(x[1] - x[0]) - p[0]
            return 4.184 * sum(p[n - 1] * d**n for n in (2, 3, 4))

        for h in declaration["tolerances"]["finite_difference_angstrom"]:
            f = np.zeros_like(xyz)
            for i in range(2):
                for j in range(3):
                    plus, minus = xyz.copy(), xyz.copy()
                    plus[i, j] += h
                    minus[i, j] -= h
                    f[i, j] = -(E(plus) - E(minus)) / (2 * h)
            fd.append(
                {
                    "displacement_angstrom": h,
                    "maximum_force_error": float(np.max(np.abs(f - force))),
                }
            )
        assert (
            fd[-1]["maximum_force_error"]
            < declaration["tolerances"]["force_atol_kj_mol_angstrom"]
        )
        report["terms"].append(
            {
                "source": item,
                "lammps": comparison,
                "analytic": analytic,
                "finite_differences": fd,
                **meta,
            }
        )
    for c in declaration["cases"]:
        directory = args.output / c["name"]
        directory.mkdir()
        observation = decode(
            json.loads((args.audit / c["name"] / "inspection.json").read_text())
        )
        if not observation["charges_complete"]:
            report["converter_controls"].append(
                {
                    "case": c["name"],
                    "executed": False,
                    "reason": "Native charge failure; no substitute charge vector",
                }
            )
            continue
        path = Path(c["retained_input"]["path"])
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == c["retained_input"]["sha256"]
        )
        system = system_from(decode(json.loads(path.read_text())))
        labels = {
            i: entry["type"] for i, entry in observation["typing_entries"].items()
        }
        charges = observation["native_charge_record"]["native_charge_record"][
            "partial_charges"
        ]
        xyz = np.array([system.coordinates.get(i) for i in sorted(labels)])
        reference_inputs(system, labels, charges, xyz, directory)
        cmd = [
            str(args.converter.resolve()),
            "reference",
            "-class",
            "II",
            "-frc",
            str(args.source.resolve()),
            "-p",
            "3",
            "-nocenter",
        ]
        run = subprocess.run(
            cmd, cwd=directory, capture_output=True, text=True, timeout=180, check=False
        )
        (directory / "stdout.txt").write_text(run.stdout)
        (directory / "stderr.txt").write_text(run.stderr)
        report["converter_controls"].append(
            {
                "case": c["name"],
                "executed": True,
                "command": cmd,
                "returncode": run.returncode,
                "output_sha256": hashlib.sha256(run.stdout.encode()).hexdigest(),
                "not_an_auto_equivalence_or_charge_oracle": True,
            }
        )
    report["hashes"] = {
        str(p.relative_to(args.output)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in args.output.rglob("*")
        if p.is_file()
    }
    write(args.output / "report.json", report)
    print("Isolated controls recorded; whole-system acceptance remains unmet")


if __name__ == "__main__":
    main()
