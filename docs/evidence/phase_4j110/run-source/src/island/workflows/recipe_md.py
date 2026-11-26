"""Native LAMMPS input export and checkpointed recipe segments."""

import json
import re
import shutil

import numpy as np

from island.epoxy.dynamics import invoke, last_frame
from island.forcefields.emc_graph import _section
from island.forcefields.emc_md import _check_dump


def radonpy_input(assignment, box):
    """Use RadonPy's native topology writer and retain its coefficient precision."""
    from radonpy.core.utils import Cell
    from radonpy.sim.lammps import MolToLAMMPSdataBlock

    mol = assignment.mol
    mol.cell = Cell(box[0], 0, box[1], 0, box[2], 0)
    data = "\n".join(MolToLAMMPSdataBlock(mol, velocity=False)) + "\n"
    params = [
        "pair_style lj/cut/coul/long 9.5",
        "pair_modify mix arithmetic",
        f"special_bonds lj 0 0 {assignment.lj14_scale:.16g} coul 0 0 {assignment.coulomb14_scale:.16g}",
        "bond_style harmonic",
        "angle_style harmonic",
        "dihedral_style fourier",
        "improper_style cvff",
    ]
    for name, command in [
        ("Pair", "pair"),
        ("Bond", "bond"),
        ("Angle", "angle"),
        ("Dihedral", "dihedral"),
        ("Improper", "improper"),
    ]:
        section = name + " Coeffs"
        values = _section(data, section)
        if values:
            for row in values.splitlines():
                i, rest = row.split(maxsplit=1)
                ids = f"{i} {i}" if command == "pair" else i
                params.append(f"{command}_coeff {ids} {rest}")
            data = data.replace("\n" + section + "\n\n" + values + "\n", "\n", 1)
    return data, "\n".join(params) + "\n"


def run_segment(
    source, out, md, steps, ensemble, seed, lammps, *, restart=False, velocities=False
):
    """Replayable segment with persistent velocities and periodic PPPM setup."""
    for name in ("system.data", "system.params"):
        shutil.copyfile(source / name, out / name)
    if restart:
        shutil.copyfile(source / "final.restart", out / "start.restart")
    params = (out / "system.params").read_text()
    kspace = (
        "kspace_style pppm 1e-5"
        if re.search(r"(?m)^pair_style\s+\S*/long\b", params)
        else ""
    )
    read = (
        "read_restart start.restart"
        if restart
        else "read_data system.data\ninclude system.params"
    )
    initialize = (
        ""
        if restart or velocities
        else f"velocity all create {md.temperature_k} {seed} mom yes rot no dist gaussian"
    )
    minimize = (
        ""
        if restart or velocities
        else "min_style cg\nmin_modify dmax 0.08\nminimize 1e-7 1e-5 2000 10000"
    )
    thermostat = f"{ensemble} temp {md.temperature_k} {md.temperature_k} 100"
    if ensemble == "npt":
        thermostat += f" iso {md.pressure_atm} {md.pressure_atm} 1000"
    columns = "id mol type xu yu zu vx vy vz"
    endpoint = (
        f"write_dump all custom trajectory.lammpstrj {columns} modify append yes sort id"
        if steps % md.sample_every
        else ""
    )
    script = f"""units real
atom_style full
boundary p p p
{read}
{kspace}
neighbor 2 bin
neigh_modify every 1 delay 0 check yes
thermo {md.sample_every}
thermo_style custom step temp pe ke etotal press density vol
thermo_modify lost error flush yes
{minimize}
reset_timestep 0
timestep {md.timestep_fs}
{initialize}
fix bath all {thermostat}
variable island_t equal temp
variable island_pe equal pe
variable island_p equal press
variable island_rho equal density
variable island_v equal vol
fix island_samples all ave/time 1 1 {md.sample_every} v_island_t v_island_pe v_island_p v_island_rho v_island_v file thermo.dat format " %.16g"
dump traj all custom {md.sample_every} trajectory.lammpstrj {columns}
dump_modify traj sort id format float %.16g
run 0
print "$(pe:%.16g)" file start.txt
run {steps} every {md.segment_steps} NULL
undump traj
{endpoint}
write_restart final.restart
write_dump all custom final.dump id x y z vx vy vz fx fy fz modify sort id format float %.16g
print "$(step:%.0f) $(temp:%.16g) $(pe:%.16g) $(press:%.16g) $(density:%.16g)" file end.txt
"""
    (out / "md.in").write_text(script)
    invoke(lammps, ["-nonbuf", "-in", "md.in"], out, "screen.log")
    values = np.loadtxt(out / "end.txt")
    step, xyz, box, _, _ = last_frame(out / "final.dump")
    if (
        step != steps
        or values.shape != (5,)
        or values[0] != steps
        or not np.isfinite(values).all()
        or not (out / "final.restart").stat().st_size
    ):
        raise ValueError("Incomplete recipe dynamics segment")
    frames = _check_dump(
        out / "trajectory.lammpstrj",
        set(range(1, len(xyz) + 1)),
        steps,
        md.sample_every,
    )
    result = {
        "steps": steps,
        "time_fs": steps * md.timestep_fs,
        "ensemble": ensemble,
        "temperature_k": float(values[1]),
        "potential_kcal_mol": float(values[2]),
        "pressure_atm": float(values[3]),
        "density_g_cm3": float(values[4]),
        "box_A": box.tolist(),
        "particles": len(xyz),
        "frames": frames,
    }
    (out / "result.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    return result
