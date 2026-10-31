"""Actual isolated LAMMPS LJ96/Coulomb and special-pair numerical control."""

import argparse
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write


def run(engine, work):
    engine, work = Path(engine).resolve(), Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    x = np.array([[0, 0, 0], [1, 0, 0], [1.3, 1.2, 0], [1.3, 1.4, 1.5]], float)
    data = """J60 explicit special-pair control

4 atoms
3 bonds

2 atom types
1 bond types

-60 60 xlo xhi
-60 60 ylo yhi
-60 60 zlo zhi

Masses

1 1
2 1

Atoms # full

1 1 1 0.13 0 0 0
2 1 1 -0.13 1 0 0
3 1 2 -0.13 1.3 1.2 0
4 1 2 0.13 1.3 1.4 1.5

Bonds

1 1 1 2
2 1 2 3
3 1 3 4
"""
    (work / "control.data").write_text(data)
    (work / "control.in").write_text("""units real
atom_style full
boundary f f f
bond_style zero
read_data control.data
bond_coeff *
pair_style lj/class2/coul/cut 100
dielectric 1
pair_coeff 1 1 0.04 3.7
pair_coeff 2 2 0.03 3.1
special_bonds lj/coul 0 0 1 angle no dihedral no
pair_modify shift no tail no
thermo 1
run 0
print "J60_E $(pe:%.17g)"
write_dump all custom forces.dump id fx fy fz modify sort id format float %.17g
""")
    with (work / "stdout.log").open("w") as stream:
        p = subprocess.run(
            [str(engine), "-in", "control.in"],
            cwd=work,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=60,
        )
    require(p.returncode == 0, "pair_control_failed")
    out = (work / "stdout.log").read_text()
    engine_energy = float(re.findall(r"^J60_E ([-+\d.eE]+)", out, re.MULTILINE)[-1])
    force = np.array(
        [
            [float(v) for v in s.split()[1:]]
            for s in (work / "forces.dump").read_text().splitlines()[-4:]
        ]
    )
    r = float(np.linalg.norm(x[0] - x[3]))
    rm = ((3.7**6 + 3.1**6) / 2) ** (1 / 6)
    eps = 2 * np.sqrt(0.04 * 0.03) * 3.7**3 * 3.1**3 / (3.7**6 + 3.1**6)
    lj = eps * (2 * (rm / r) ** 9 - 3 * (rm / r) ** 6)
    coul = 332.06371 * 0.13**2 / r
    derivative = (
        eps * (-18 * (rm / r) ** 9 + 18 * (rm / r) ** 6) / r
        - 332.06371 * 0.13**2 / r**2
    )
    expected = np.zeros((4, 3))
    expected[0] = -derivative * (x[0] - x[3]) / r
    expected[3] = -expected[0]
    de = abs(engine_energy - lj - coul)
    df = float(np.max(abs(force - expected)))
    require(
        de < 1e-11 and df < 1e-10,
        "pair_units_scaling_mismatch",
        energy_error=de,
        force_error=df,
    )
    return signed(
        {
            "schema": "island.j60-special-pair-native-control.v1",
            "engine_sha256": pin(engine),
            "units": "LAMMPS real",
            "boundary": "f f f",
            "cutoff_A": 100,
            "dielectric": 1,
            "prefactor": 332.06371,
            "special_bonds": [0, 0, 1],
            "included_pairs": [[1, 4]],
            "excluded_pairs": [[1, 2], [2, 3], [3, 4], [1, 3], [2, 4]],
            "mixed_rmin": rm,
            "mixed_epsilon": eps,
            "LJ_energy": lj,
            "Coulomb_energy": coul,
            "native_energy": engine_energy,
            "energy_error": de,
            "force_max_error": df,
            "test_only_values_not_source_rows": True,
            "PCFF_MD_executed": False,
            "files": {p.name: pin(p) for p in work.iterdir() if p.is_file()},
        }
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    write(a.output, run(a.engine, a.work))
