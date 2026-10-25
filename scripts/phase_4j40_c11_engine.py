"""Fresh engine single points. Study settings never modify a production backend."""

import math
import re
import subprocess
from pathlib import Path

import numpy as np

from scripts.phase_4j40_c11_numeric import (
    ATM_TO_BAR,
    COULOMB_OPENMM,
    KCAL_TO_KJ,
    analytic_ewald_terms,
    parse_dump,
    require,
    sha,
    write,
)


def openmm_point(backend, xyz, lengths, setting, frame_identity, directory):
    import openmm as mm
    from openmm import unit

    from island._periodic_openmm import build

    directory.mkdir(parents=True, exist_ok=False)
    cfg = dict(
        backend["config"],
        box_lengths=list(np.asarray(lengths) * 10),
        electrostatics_method=setting["method"],
        nonbonded_cutoff=setting["cutoff_nm"] * 10,
        pme_tolerance=setting["tolerance"],
    )
    model = build(backend["numerical_data"], cfg, mm)
    coul = next(f for f in model.getForces() if isinstance(f, mm.NonbondedForce))
    coul.setReciprocalSpaceForceGroup(29)
    if "alpha_nm" in setting:
        coul.setPMEParameters(setting["alpha_nm"], *setting["mesh"])
    integrator = mm.VerletIntegrator(0.0005)
    platform = mm.Platform.getPlatformByName(setting.get("platform", "Reference"))
    context = mm.Context(model, integrator, platform)
    try:
        context.setPositions(xyz)
        state = context.getState(getEnergy=True, getForces=True)
        components = {
            k: 0.0 for k in ("bond", "angle", "proper", "improper", "lj", "coulomb")
        }
        raw_groups = {}
        for group in sorted({f.getForceGroup() for f in model.getForces()} | {29}):
            e = (
                context.getState(getEnergy=True, groups=1 << group)
                .getPotentialEnergy()
                .value_in_unit(unit.kilojoule_per_mole)
            )
            names = {
                f.getName() for f in model.getForces() if f.getForceGroup() == group
            }
            name = (
                "coulomb_reciprocal_including_self_background"
                if group == 29
                else "coulomb_real_plus_exception"
                if group == 30
                else "lj"
                if group == 31
                else next(iter(names))
            )
            raw_groups[name] = float(e)
            key = "coulomb" if group in (29, 30) else name
            components[key] += float(e)
        if setting["method"] == "pme":
            alpha, nx, ny, nz = coul.getPMEParametersInContext(context)
            mesh = [nx, ny, nz]
        else:
            alpha = (
                math.sqrt(-math.log(2 * setting["tolerance"])) / setting["cutoff_nm"]
            )
            mesh = None  # Direct Ewald has reciprocal vectors, not a PME mesh.
        row = {
            "engine": "OpenMM",
            "version": mm.version.version,
            "platform": platform.getName(),
            "precision": "Reference double"
            if platform.getName() == "Reference"
            else "platform-native",
            "setting": setting,
            "frame_identity": frame_identity,
            "energy_kj_mol": float(
                state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            ),
            "components_kj_mol": components,
            "raw_groups_kj_mol": raw_groups,
            "forces_kj_mol_nm": state.getForces(asNumpy=True)
            .value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            .tolist(),
            "alpha_nm": float(alpha),
            "mesh": mesh,
            "pme_interpolation_order": 5 if setting["method"] == "pme" else None,
            "constraints": model.getNumConstraints(),
            "analytic_terms": analytic_ewald_terms(
                [s["charge"] for s in backend["numerical_data"]["sites"]],
                float(alpha),
                float(np.prod(lengths)),
                COULOMB_OPENMM,
            ),
        }
        (directory / "system.xml").write_text(mm.XmlSerializer.serialize(model))
        row["system_xml_sha256"] = sha(directory / "system.xml")
        write(directory / "point.json", row)
        return row
    finally:
        del context, integrator


def lammps_point(
    template, xyz, lengths, setting, frame_identity, directory, executable
):
    directory.mkdir(parents=True, exist_ok=False)
    lines = template.splitlines()
    in_atoms = False
    for i, line in enumerate(lines):
        if re.search(r"\b[xyz]lo\s+[xyz]hi\s*$", line):
            axis = "xyz".index(line.split()[-2][0])
            lines[i] = (
                f"0 {lengths[axis] * 10:.17g} {line.split()[-2]} {line.split()[-1]}"
            )
        if line.startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms and line and line[0].isalpha():
            in_atoms = False
        if in_atoms and line.strip():
            values = line.split()
            sid = int(values[0])
            values[4:7] = [format(v * 10, ".17g") for v in xyz[sid - 1]]
            lines[i] = " ".join(values)
    (directory / "system.data").write_text("\n".join(lines) + "\n")
    extra = []
    if "alpha_nm" in setting:
        extra.append(f"kspace_modify gewald {setting['alpha_nm'] / 10:.17g}")
    if "mesh" in setting:
        extra.append("kspace_modify mesh " + " ".join(map(str, setting["mesh"])))
    if "order" in setting:
        extra.append(f"kspace_modify order {setting['order']}")
    if "kmax" in setting:
        extra.append(
            f"kspace_modify kmax/ewald {setting['kmax']} {setting['kmax']} {setting['kmax']}"
        )
    script = f"""units real
atom_style full
boundary p p p
bond_style harmonic
angle_style harmonic
dihedral_style fourier
improper_style cvff
pair_style lj/cut/coul/long {setting["cutoff_nm"] * 10:.17g}
pair_modify mix arithmetic shift no tail no table 0
special_bonds lj 0 0 0.5 coul 0 0 0.83333333333333337 angle no dihedral no
read_data system.data
kspace_style {setting["method"]} {setting["tolerance"]:.17g}
{chr(10).join(extra)}
neighbor 2 bin
neigh_modify every 1 delay 0 check yes
compute cvir all pressure NULL virial
thermo 1
thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong c_cvir c_cvir[1] c_cvir[2] c_cvir[3] c_cvir[4] c_cvir[5] c_cvir[6]
thermo_modify format float %.17g
dump forces all custom 1 forces.dump id type q x y z fx fy fz
dump_modify forces sort id format float %.17g
run 0
"""
    (directory / "in.study").write_text(script)
    argv = [str(Path(executable).resolve()), "-in", "in.study"]
    proc = subprocess.run(
        argv, cwd=directory, text=True, capture_output=True, check=False
    )
    (directory / "execution.log").write_text(proc.stdout + proc.stderr)
    write(
        directory / "command.json",
        {"argv": argv, "cwd": str(directory.resolve()), "exit_code": proc.returncode},
    )
    require(proc.returncode == 0, "LAMMPS single point failed: " + str(directory))
    lines = proc.stdout.splitlines()
    h = next(i for i, line in enumerate(lines) if "PotEng" in line)
    values = list(map(float, lines[h + 1].split()))
    components = dict(
        zip(
            (
                "total",
                "bond",
                "angle",
                "proper",
                "improper",
                "lj",
                "real",
                "reciprocal",
            ),
            [v * KCAL_TO_KJ for v in values[1:9]],
            strict=True,
        )
    )
    alpha = (
        float(re.search(r"G vector \(1/distance\)\s*=\s*([\deE+.-]+)", proc.stdout)[1])
        * 10
    )
    dump = parse_dump(directory / "forces.dump")
    delta = np.asarray(dump["positions_nm"]) - xyz
    image_shifts = np.rint(delta / lengths)
    require(
        np.allclose(delta - image_shifts * lengths, 0, rtol=0, atol=3e-15),
        "LAMMPS changed atom coordinates beyond periodic image wrapping",
    )
    require(
        np.allclose(dump["box_lengths_nm"], lengths, rtol=0, atol=3e-15),
        "LAMMPS changed box",
    )
    row = {
        "engine": "LAMMPS",
        "setting": setting,
        "frame_identity": frame_identity,
        "energy_kj_mol": components["total"],
        "components_kj_mol": {
            **{k: components[k] for k in ("bond", "angle", "proper", "improper", "lj")},
            "coulomb": components["real"] + components["reciprocal"],
        },
        "raw_groups_kj_mol": components,
        "forces_kj_mol_nm": dump["forces_kj_mol_nm"],
        "atom_ids": dump["atom_ids"],
        "charges_e": dump["charges_e"],
        "periodic_image_shifts": image_shifts.astype(int).tolist(),
        "coordinate_policy": "same exported coordinates; dump wrapping checked modulo actual orthorhombic box",
        "configurational_pressure_bar": values[9] * ATM_TO_BAR,
        "configurational_pressure_tensor_bar": [v * ATM_TO_BAR for v in values[10:16]],
        "kinetic_pressure_bar": 0.0,
        "temperature": "zero velocities for single-point control; not a thermodynamic observation",
        "alpha_nm_from_8_digit_log": alpha,
        "analytic_terms": analytic_ewald_terms(
            dump["charges_e"], alpha, float(np.prod(lengths)), 138.935456264
        ),
        "analytic_alpha_precision_note": "native initialization log rounds alpha; these self/background terms are descriptive, not authoritative component subtraction",
        "raw_files": {p.name: sha(p) for p in directory.iterdir() if p.is_file()},
    }
    write(directory / "point.json", row)
    return row


def finite_strain(
    backend,
    xyz,
    lengths,
    setting,
    frame_identity,
    directory,
    engine,
    template,
    executable,
):
    rows = []
    for h in (1e-4, 1e-5, 1e-6):
        energies = []
        for sign in (-1, 1):
            scale = (1 + sign * h) ** (1 / 3)
            path = directory / f"h-{h}-{sign}"
            if engine == "openmm":
                point = openmm_point(
                    backend, xyz * scale, lengths * scale, setting, frame_identity, path
                )
            else:
                point = lammps_point(
                    template,
                    xyz * scale,
                    lengths * scale,
                    setting,
                    frame_identity,
                    path,
                    executable,
                )
            energies.append(point["components_kj_mol"])
        from scripts.phase_4j40_c11_numeric import KJ_MOL_NM3_TO_BAR

        pressure = {
            k: -(energies[1][k] - energies[0][k])
            / (2 * h * np.prod(lengths))
            * KJ_MOL_NM3_TO_BAR
            for k in energies[0]
        }
        rows.append(
            {
                "relative_volume_step": h,
                "configurational_pressure_components_bar": pressure,
                "total_bar": math.fsum(pressure.values()),
            }
        )
    return rows
