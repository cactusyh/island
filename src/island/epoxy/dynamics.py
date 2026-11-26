"""LAMMPS real-units segments. Retain velocities, stable IDs, box and native inputs."""

import json
import os
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem

from .parameters import assign_network, data_text


def invoke(executable, args, out, log, stdin=None):
    if os.path.dirname(executable):
        executable = os.path.abspath(executable)
    with (out / log).open("w") as stream:
        p = subprocess.run(
            [executable, *args],
            check=False,
            cwd=out,
            stdin=stdin,
            text=True,
            stdout=stream,
            stderr=subprocess.STDOUT,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        )
    if p.returncode:
        raise RuntimeError((out / log).read_text()[-2400:])


def pack(monomers, counts, density, seed, packmol, out):
    mass = sum(
        n * sum(a.GetMass() for a in m.GetAtoms())
        for m, n in zip(monomers, counts, strict=True)
    )
    side = (mass / (0.602214076 * density)) ** (1 / 3)
    script = f"filetype xyz\ntolerance 1.9\nseed {seed}\noutput packed.xyz\npbc {side} {side} {side}\n"
    symbols = []
    for k, (m, n) in enumerate(zip(monomers, counts, strict=True)):
        (out / f"monomer{k}.xyz").write_text(Chem.MolToXYZBlock(m))
        script += f"structure monomer{k}.xyz\n number {n}\n inside box 0 0 0 {side} {side} {side}\nend structure\n"
        symbols.extend([a.GetSymbol() for a in m.GetAtoms()] * n)
    (out / "pack.in").write_text(script)
    with (out / "pack.in").open() as inp:
        invoke(packmol, [], out, "pack.log", inp)
    if "Success!" not in (out / "pack.log").read_text():
        raise RuntimeError("Packing failed to satisfy requested tolerances")
    rows = [line.split() for line in (out / "packed.xyz").read_text().splitlines()[2:]]
    if [r[0] for r in rows] != symbols:
        raise ValueError("PACKMOL atom count/order changed")
    return np.array([r[1:4] for r in rows], float), np.full(3, side)


def last_frame(path):
    frame = None
    with Path(path).open() as h:
        for line in h:
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError("Invalid dump")
            step = int(next(h))
            if next(h).strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("Invalid particle count header")
            n = int(next(h))
            if next(h).strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ValueError("Invalid cell")
            bounds = np.array([next(h).split() for _ in range(3)], float)
            if next(h).split()[2:] != [
                "id",
                "x",
                "y",
                "z",
                "vx",
                "vy",
                "vz",
                "fx",
                "fy",
                "fz",
            ]:
                raise ValueError("Invalid dump fields")
            table = np.array([next(h).split() for _ in range(n)], float)
            if (
                not np.array_equal(table[:, 0], np.arange(1, n + 1))
                or not np.isfinite(table).all()
            ):
                raise ValueError("Invalid particle IDs or nonfinite dynamics")
            frame = (
                step,
                table[:, 1:4] - bounds[:, 0],
                bounds[:, 1] - bounds[:, 0],
                table[:, 4:7],
                table[:, 7:],
            )
    if frame is None:
        raise ValueError("Empty trajectory")
    return frame


def segment(
    out, mol, lib, xyz, box, vel, *, steps, temp, seed, lammps, npt=True, minimize=True
):
    if lib.get("backend") == "radonpy-graph":
        from .gaff import parameterize

        data, params, report = parameterize(mol, lib, xyz, box, vel)
        (out / "parameterization.json").write_text(json.dumps(report, indent=2) + "\n")
    else:
        atoms, terms = assign_network(mol, lib)
        membership = [a.GetIntProp("monomer") + 1 for a in mol.GetAtoms()]
        data = data_text(atoms, terms, lib["counts"], xyz, box, membership, vel)
        params = lib["parameters"]
    (out / "system.data").write_text(data)
    (out / "system.params").write_text(params)
    ensemble = (
        f"npt temp {temp} {temp} 50 iso 1 1 500"
        if npt
        else f"nvt temp {temp} {temp} 50"
    )
    init = (
        f"velocity all create {temp} {seed} mom yes dist gaussian"
        if vel is None
        else ""
    )
    script = f"""units real
atom_style full
boundary p p p
read_data system.data
include system.params
kspace_style pppm 1e-5
neighbor 2 bin
neigh_modify every 1 delay 0 check yes
thermo {max(1, steps // 10)}
thermo_style custom step temp pe ke press density vol
thermo_modify lost error
min_style cg
min_modify dmax 0.08
{"minimize 1e-7 1e-5 1000 4000" if minimize else ""}
reset_timestep 0
timestep 0.5
{init}
fix bath all {ensemble}
dump sample all custom {max(1, steps // 5)} traj.dump id x y z vx vy vz fx fy fz
dump_modify sample sort id format float %.16g
run {steps}
write_dump all custom final.dump id x y z vx vy vz fx fy fz modify sort id format float %.16g
print "$(step:%.0f) $(temp:%.12g) $(pe:%.12g) $(press:%.12g) $(density:%.12g)" file thermo.txt
"""
    (out / "md.in").write_text(script)
    invoke(lammps, ["-in", "md.in"], out, "screen.log")
    end, xyz, box, vel, _ = last_frame(out / "final.dump")
    thermo = np.loadtxt(out / "thermo.txt")
    if end != steps or not np.isfinite(thermo).all():
        raise RuntimeError("Incomplete MD segment")
    return (
        xyz,
        box,
        vel,
        {
            "steps": end,
            "temperature_K": thermo[1],
            "potential_kcal_mol": thermo[2],
            "pressure_atm": thermo[3],
            "density_g_cm3": thermo[4],
        },
    )
