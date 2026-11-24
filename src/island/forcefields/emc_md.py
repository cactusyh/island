"""Run finite-graph EMC assignments with their native LAMMPS potential."""

import json
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from math import isfinite
from pathlib import Path

from island.exceptions import ForceFieldError
from island.forcefields.emc_graph import EMCGraphAssignment


@dataclass(frozen=True)
class EMCMDResult:
    directory: Path
    steps: int
    time_fs: float
    initial_potential_kcal_mol: float
    final_potential_kcal_mol: float
    final_temperature_k: float
    frames: int


def _check_dump(path, particles, steps, stride):
    """Stream and verify every expected frame, including an off-stride endpoint."""
    frames = 0
    last = -1
    with path.open() as f:
        for header in f:
            if header.strip() != "ITEM: TIMESTEP":
                raise ForceFieldError("Invalid LAMMPS trajectory header")
            step = int(f.readline())
            if last == steps or step != (0 if last < 0 else min(last + stride, steps)):
                raise ForceFieldError("Unexpected LAMMPS trajectory timestep")
            if f.readline().strip() != "ITEM: NUMBER OF ATOMS" or int(
                f.readline()
            ) != len(particles):
                raise ForceFieldError("Particle count changed during dynamics")
            if f.readline().strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ForceFieldError("Unexpected trajectory cell")
            for _ in range(3):
                lo, hi = map(float, f.readline().split())
                if not (isfinite(lo) and isfinite(hi) and lo < hi):
                    raise ForceFieldError("Invalid trajectory cell")
            if f.readline().split() != [
                "ITEM:",
                "ATOMS",
                "id",
                "mol",
                "type",
                "xu",
                "yu",
                "zu",
                "vx",
                "vy",
                "vz",
            ]:
                raise ForceFieldError("Unexpected trajectory columns")
            seen = set()
            for _ in particles:
                values = f.readline().split()
                if len(values) != 9 or not all(isfinite(float(x)) for x in values[3:]):
                    raise ForceFieldError("Incomplete or non-finite trajectory")
                seen.add(int(values[0]))
            if seen != particles:
                raise ForceFieldError("Particle identities changed during dynamics")
            last, frames = step, frames + 1
    if last != steps:
        raise ForceFieldError("LAMMPS trajectory ended before the requested step")
    return frames


def run_emc_md(
    assignment: EMCGraphAssignment,
    directory,
    *,
    steps=1000,
    timestep_fs=0.5,
    temperature_k=300.0,
    damping_fs=100.0,
    seed=2026,
    sample_every=100,
    lammps="lmp",
):
    """Minimize and run fixed-volume Nosé–Hoover dynamics in a fresh output directory.

    The directory contains a replayable input, native coefficients, site map,
    unwrapped trajectory, log and binary restart. result.json marks success.
    """
    if not isinstance(assignment, EMCGraphAssignment):
        raise TypeError("Expected the result of assign_emc_graph")
    if any(type(n) is not int or n < 1 for n in (steps, sample_every)):
        raise ValueError("steps and sample_every must be positive integers")
    if type(seed) is not int or not 1 <= seed <= 900_000_000:
        raise ValueError("LAMMPS velocity seed must be in [1, 900000000]")
    if not all(
        not isinstance(v, bool) and isfinite(v) and v > 0
        for v in (temperature_k, timestep_fs, damping_fs)
    ):
        raise ValueError(
            "Temperature, timestep and damping time must be finite and positive"
        )
    exe = shutil.which(os.fspath(lammps))
    if not exe:
        raise FileNotFoundError(f"LAMMPS executable not found: {lammps}")
    exe = str(Path(exe).resolve())
    out = Path(directory).resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "chain.data").write_text(assignment.data)
    (out / "chain.params").write_text(assignment.parameters)
    (out / "assignment.json").write_text(
        json.dumps(
            {
                "field": assignment.field,
                "emc_version": assignment.emc_version,
                "field_sha256": assignment.field_sha256,
                "emc_groups": assignment.emc_groups,
                "input_smiles": assignment.input_smiles,
                "input_smiles_site_ids": assignment.input_smiles_site_ids,
                "site_groups": assignment.site_groups,
                "added_hydrogen_ids": assignment.added_hydrogen_ids,
                "coordinate_source": assignment.coordinate_source,
            },
            indent=2,
        )
        + "\n"
    )
    columns = "id mol type xu yu zu vx vy vz"
    long_range = (
        "kspace_style pppm 1.0e-4"
        if re.search(r"(?m)^pair_style\s+\S*/long\b", assignment.parameters)
        else ""
    )
    final_dump = (
        f"write_dump all custom trajectory.lammpstrj {columns} modify append yes sort id"
        if steps % sample_every
        else ""
    )
    script = f"""units real
atom_style full
boundary p p p
read_data chain.data
include chain.params
{long_range}
thermo {sample_every}
thermo_style custom step temp pe ke etotal press
thermo_modify lost error flush yes
min_style cg
minimize 1.0e-6 1.0e-6 2000 10000
reset_timestep 0
timestep {timestep_fs}
velocity all create {temperature_k} {seed} mom yes rot no dist gaussian
fix ensemble all nvt temp {temperature_k} {temperature_k} {damping_fs}
dump traj all custom {sample_every} trajectory.lammpstrj {columns}
dump_modify traj sort id
run 0
print "ISLAND_START $(pe:%.16g)" file endpoints.txt
run {steps}
print "ISLAND_END $(step:%.0f) $(pe:%.16g) $(temp:%.16g)" append endpoints.txt
undump traj
{final_dump}
write_restart final.restart
"""
    (out / "md.in").write_text(script)
    with (out / "screen.log").open("w") as screen:
        proc = subprocess.run(
            [exe, "-in", "md.in"],
            cwd=out,
            stdout=screen,
            stderr=subprocess.STDOUT,
            env={
                **os.environ,
                "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "1"),
            },
            check=False,
        )
    if proc.returncode:
        with (out / "screen.log").open("rb") as screen:
            screen.seek(max(0, screen.seek(0, 2) - 2500))
            detail = screen.read().decode(errors="replace")
        raise ForceFieldError(f"LAMMPS exited {proc.returncode} in {out}:\n{detail}")
    log = (out / "endpoints.txt").read_text()
    start = re.search(r"(?m)^ISLAND_START (\S+)$", log)
    end = re.search(r"(?m)^ISLAND_END (\d+) (\S+) (\S+)$", log)
    if (
        not start
        or not end
        or int(end[1]) != steps
        or not (out / "final.restart").stat().st_size
    ):
        raise ForceFieldError("LAMMPS did not complete the requested dynamics")
    energies = tuple(map(float, (start[1], end[2], end[3])))
    if not all(isfinite(x) for x in energies):
        raise ForceFieldError("Non-finite LAMMPS thermodynamics")
    frames = _check_dump(
        out / "trajectory.lammpstrj", set(assignment.site_groups), steps, sample_every
    )
    result = EMCMDResult(out, steps, steps * timestep_fs, *energies, frames)
    (out / "result.json").write_text(
        json.dumps(
            {k: v for k, v in asdict(result).items() if k != "directory"}, indent=2
        )
        + "\n"
    )
    return result
