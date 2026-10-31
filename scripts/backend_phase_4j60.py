"""Bounded synthetic LAMMPS class2 static controls, not a physical PCFF fit/MD."""

import argparse
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hamiltonian import bend, cos_phi, geometry, harmonic, length


def control(xyz):
    x = geometry(xyz)
    s = [1, 2, 3, 4]
    c = cos_phi(x, s)
    r = [
        length(x, p) - v
        for p, v in zip(([1, 2], [2, 3], [3, 4]), (1.1, 1.2, 1.3), strict=True)
    ]
    t = [
        bend(x, p) - v for p, v in zip(([1, 2, 3], [2, 3, 4]), (1.7, 1.8), strict=True)
    ]
    e = sum(v * (1 + harmonic(c, n)) for n, v in enumerate((0.3, 0.15, -0.25), 1))
    for d, vs in zip(
        (r[1], r[0], r[2], t[0], t[1]),
        (
            (0.4, -0.2, 0.1),
            (0.15, 0.24, -0.33),
            (-0.14, 0.22, 0.31),
            (0.2, 0.3, -0.4),
            (-0.1, 0.5, 0.2),
        ),
        strict=True,
    ):
        e += d * sum(v * harmonic(c, n) for n, v in enumerate(vs, 1))
    e += 0.31 * t[0] * t[1] * c - 0.29 * r[0] * r[2]
    return e.v, -e.d.reshape(9, 3)[:4]


def run(engine, work):
    engine, work = Path(engine).resolve(), Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    results = []
    for index, phi in enumerate((0.3, -1.1, np.pi - 1e-6, np.pi + 1e-6)):
        xyz = np.array(
            [
                [1.3, 0, 0.2],
                [0, 0, 0],
                [0, 0, 1.5],
                [1.2 * np.cos(phi), 1.2 * np.sin(phi), 1.3],
            ]
            + [[5 + i, 6, 7] for i in range(5)],
            float,
        )
        d = work / str(index)
        d.mkdir(exist_ok=True)
        data = "J60 synthetic kernel control\n\n4 atoms\n1 dihedrals\n\n1 atom types\n1 dihedral types\n\n-20 20 xlo xhi\n-20 20 ylo yhi\n-20 20 zlo zhi\n\nMasses\n\n1 1.0\n\nAtoms # full\n\n"
        data += "\n".join(
            f"{i + 1} 1 1 0 " + " ".join(format(v, ".17g") for v in row)
            for i, row in enumerate(xyz[:4])
        )
        data += "\n\nDihedrals\n\n1 1 1 2 3 4\n"
        (d / "control.data").write_text(data)
        body = """units real
atom_style full
boundary f f f
dihedral_style class2
read_data control.data
pair_style zero 10
pair_coeff * *
dihedral_coeff 1 0.3 180 0.15 180 -0.25 180
dihedral_coeff 1 mbt 0.4 -0.2 0.1 1.2
dihedral_coeff 1 ebt 0.15 0.24 -0.33 -0.14 0.22 0.31 1.1 1.3
dihedral_coeff 1 at 0.2 0.3 -0.4 -0.1 0.5 0.2 THETA1 THETA2
dihedral_coeff 1 aat 0.31 THETA1 THETA2
dihedral_coeff 1 bb13 -0.29 1.1 1.3
thermo 1
run 0
print "J60_E $(edihed:%.17g)"
write_dump all custom forces.dump id fx fy fz modify sort id format float %.17g
""".replace("THETA1", format(1.7 * 180 / np.pi, ".17g")).replace(
            "THETA2", format(1.8 * 180 / np.pi, ".17g")
        )
        (d / "control.in").write_text(body)
        with (d / "stdout.log").open("w") as stream:
            p = subprocess.run(
                [str(engine), "-in", "control.in"],
                cwd=d,
                env={**os.environ, "OMP_NUM_THREADS": "1"},
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
                timeout=60,
            )
        require(
            p.returncode == 0, "static_kernel_control_failed", log=str(d / "stdout.log")
        )
        text = (d / "stdout.log").read_text()
        en = float(re.findall(r"^J60_E ([-+\deE.]+)", text, re.MULTILINE)[-1])
        force = np.array(
            [
                [float(v) for v in row.split()[1:]]
                for row in (d / "forces.dump").read_text().splitlines()[-4:]
            ]
        )
        ee, ff = control(xyz)
        de = abs(ee - en)
        df = float(np.max(abs(ff - force)))
        require(
            de < 1e-11 and df < 1e-9,
            "native_class2_kernel_mismatch",
            energy_error=de,
            force_error=df,
        )
        results.append(
            {
                "phi": phi,
                "engine_energy": en,
                "independent_energy": ee,
                "energy_error": de,
                "force_max_error": df,
                "files": {p.name: pin(p) for p in d.iterdir() if p.is_file()},
            }
        )
    return signed(
        {
            "schema": "island.j60-synthetic-backend-controls.v1",
            "engine_path": str(engine),
            "engine_sha256": pin(engine),
            "controls": results,
            "static_kernel_check": "passed",
            "physical_parameter_values_assigned": False,
            "PCFF_MD_executed": False,
            "physical_backend_readiness": "not_established",
        }
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    write(a.output, run(a.engine, a.work))
