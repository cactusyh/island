"""Matched-state NPT density/shrinkage, cooling and small-strain tensile protocols."""

from pathlib import Path

import numpy as np

from island.forcefields.emc_graph import _section

from .dynamics import invoke, last_frame


def _snapshot(stage, xyz, box, vel):
    data = (stage / "system.data").read_text()
    for section, values, start in [("Atoms", xyz, 4), ("Velocities", vel, 1)]:
        old = _section(data, section)
        rows = (
            [r.split() for r in old.splitlines()]
            if old
            else [[str(i + 1)] for i in range(len(xyz))]
        )
        for row in rows:
            row[start : start + 3] = [
                format(x, ".16g") for x in values[int(row[0]) - 1]
            ]
        new = "\n".join(" ".join(row) for row in rows)
        if old:
            data = data.replace(old + "\n", new + "\n", 1)
        else:
            data += "\n" + section + "\n\n" + new + "\n"
    import re

    for ax, L in zip("xyz", box, strict=True):
        data = re.sub(rf"(?m)^.* {ax}lo {ax}hi$", f"0 {L:.16g} {ax}lo {ax}hi", data)
    return data


def conditions(stage, out, lammps, temperature=300, steps=10000, *, segment_steps=None):
    stage, out = Path(stage), Path(out)
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be finite and positive")
    if (stage / "state.npz").is_file():
        with np.load(stage / "state.npz") as s:
            data = _snapshot(stage, s["xyz"], s["box"], s["vel"])
    else:
        data, _, _ = _from_condition(stage)
    text = f"fix bath all npt temp {temperature} {temperature} 50 iso 1 1 500\n"
    values = _measure(
        out,
        data,
        (stage / "system.params").read_text(),
        text,
        steps,
        lammps,
        segment_steps=segment_steps,
    )
    rho = values[len(values) // 2 :, 2]
    volume = values[len(values) // 2 :, 3]
    return {
        "temperature_K": temperature,
        "pressure_atm": 1,
        "steps": steps,
        "mean_density": float(rho.mean()),
        "mean_volume": float(volume.mean()),
        "density_temporal_sd": float(rho.std()),
    }


def _measure(out, data, params, protocol, steps, lammps, *, segment_steps=None):
    if type(steps) is not int or steps < 200 or steps % 100:
        raise ValueError(
            "Use a positive multiple of 100 steps with at least two sampling blocks"
        )
    run_command = f"run {steps}"
    if segment_steps is not None:
        if (
            type(segment_steps) is not int
            or segment_steps < 100
            or segment_steps % 100
            or steps % segment_steps
        ):
            raise ValueError(
                "Segment length must divide the run and align with the 100-step samples"
            )
        # LAMMPS reinitializes PPPM for the current box between NPT segments.
        run_command += f" every {segment_steps} NULL"
    out.mkdir(parents=True, exist_ok=False)
    (out / "system.data").write_text(data)
    (out / "system.params").write_text(params)
    script = f"""units real
atom_style full
read_data system.data
include system.params
kspace_style pppm 1e-5
timestep 0.5
thermo 1000
thermo_style custom step temp pe density vol press
{protocol}
variable tt equal temp
variable rho equal density
variable vv equal vol
variable xx equal lx
variable yy equal ly
variable zz equal lz
variable px equal pxx
variable py equal pyy
variable pz equal pzz
fix output all ave/time 10 10 100 v_tt v_rho v_vv v_xx v_yy v_zz v_px v_py v_pz file samples.dat
{run_command}
write_dump all custom final.dump id x y z vx vy vz fx fy fz modify sort id format float %.16g
"""
    (out / "property.in").write_text(script)
    invoke(lammps, ["-nonbuf", "-in", "property.in"], out, "screen.log")
    values = np.loadtxt(out / "samples.dat", ndmin=2)
    if values[-1, 0] != steps or not np.isfinite(values).all():
        raise ValueError("Incomplete property samples")
    return values


def _from_condition(stage):
    _, xyz, box, vel, _ = last_frame(stage / "final.dump")
    return _snapshot(stage, xyz, box, vel), (stage / "system.params").read_text(), box


def cooling(hot, out, lammps, rate, high=550, low=250):
    """Estimate a specific-volume breakpoint. Reject nonphysical or unsupported fits."""
    hot, out = Path(hot), Path(out)
    if (
        not 0 < low < high
        or not np.isfinite([low, high, rate]).all()
        or not 0 < rate <= (high - low) * 2000 / 2400
    ):
        raise ValueError(
            "Cooling requires 24 sampling blocks and a positive finite rate and range"
        )
    data, params, _ = _from_condition(hot)
    steps = round(2000 * (high - low) / rate / 100) * 100
    values = _measure(
        out,
        data,
        params,
        f"fix bath all npt temp {high} {low} 50 iso 1 1 500",
        steps,
        lammps,
    )
    bins = np.array_split(values, 24)
    x = np.array([b[:, 1].mean() for b in bins])
    y = np.array([(1 / b[:, 2]).mean() for b in bins])
    residual = lambda pred: float(np.sum((y - pred) ** 2))
    base = residual(np.polyval(np.polyfit(x, y, 1), x))
    fits = []
    for t in np.linspace(low + 0.2 * (high - low), high - 0.2 * (high - low), 81):
        design = np.stack([np.ones_like(x), x - t, np.maximum(x - t, 0)], axis=1)
        c = np.linalg.lstsq(design, y, rcond=None)[0]
        if c[1] > 0 and c[2] > 0:
            fits.append((residual(design @ c), t))
    best = min(fits) if fits else None
    # Penalize both added parameters (slope and fitted breakpoint).
    bic_gain = (
        float(len(x) * np.log(base / max(best[0], 1e-30)) - 2 * np.log(len(x)))
        if best and base > 0
        else None
    )
    accepted = bool(
        best
        and bic_gain is not None
        and bic_gain > 6
        and low + 0.2 * (high - low) < best[1] < high - 0.2 * (high - low)
    )
    return {
        "rate_K_ps": rate,
        "steps": steps,
        "tg_K": float(best[1]) if accepted else None,
        "bic_gain": bic_gain,
        "status": "resolved" if accepted else "no_supported_breakpoint",
        "bins_K": x.tolist(),
        "bins_specific_volume": y.tolist(),
    }


def tension(cold, out, lammps, rate, axis):
    """Pull at fixed engineering strain rate with the two lateral pressures controlled."""
    cold, out = Path(cold), Path(out)
    if axis not in ("x", "y", "z") or not np.isfinite(rate) or not 0 < rate <= 0.04:
        raise ValueError(
            "Use x/y/z and a positive strain rate with enough sampling blocks"
        )
    data, params, box = _from_condition(cold)
    i = "xyz".index(axis)
    transverse = [a for a in "xyz" if a != axis]
    steps = round(2000 * 0.04 / rate / 100) * 100
    protocol = (
        "fix bath all npt temp 300 300 50 "
        + " ".join(f"{a} 1 1 500" for a in transverse)
        + f" couple {''.join(transverse)}\nfix stretch all deform 1 {axis} erate {rate / 1000} units box remap x"
    )
    # LAMMPS only accepts xz, xy, yz coupling spellings.
    values = _measure(out, data, params, protocol, steps, lammps)
    strain = values[:, 4 + i] / box[i] - 1
    stress = -values[:, 7 + i] * 0.000101325
    interval = (strain >= 0.001) & (strain <= 0.03)
    blocks = np.array_split(np.stack([strain[interval], stress[interval]], axis=1), 10)
    x, y = np.array([b.mean(axis=0) for b in blocks]).T
    coef = np.polyfit(x, y, 1)
    r2 = float(1 - np.sum((y - np.polyval(coef, x)) ** 2) / np.sum((y - y.mean()) ** 2))
    return {
        "rate_per_ps": rate,
        "axis": axis,
        "steps": steps,
        "slope_GPa": float(coef[0]),
        "r_squared": r2,
        "modulus_GPa": float(coef[0]) if r2 >= 0.8 and coef[0] > 0 else None,
    }
