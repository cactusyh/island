"""J40 same-graph KSPACE Ewald single points using retained J22 comparison bounds.

This is independent-engine evidence, separate from trajectory convergence.
The unmodified production cutoff, precision, mixing and 1-4 settings are exported.
"""

import subprocess

import numpy as np

from island._periodic_openmm import build
from island.periodic import export_lammps_data, export_lammps_input


def evaluate(backend, mm, unit, directory, executable):
    model = build(backend.payload["numerical_data"], backend.payload["config"], mm)
    xyz = np.array(
        [
            backend.system.coordinates.get(i) / 10
            for i in sorted(backend.system.topology.sites)
        ]
    )
    integrator = mm.VerletIntegrator(0.001)
    ctx = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(xyz)
    state = ctx.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    force = (
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        / 10
    )
    components = {}
    for group in sorted({f.getForceGroup() for f in model.getForces()}):
        names = {f.getName() for f in model.getForces() if f.getForceGroup() == group}
        value = (
            ctx.getState(getEnergy=True, groups=1 << group)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )
        name = "lj" if group == 31 else "coulomb" if group == 30 else next(iter(names))
        components[name] = value
    h = 1e-6
    errors = []
    # Two separated atoms, all Cartesian directions, including periodic electrostatics.
    for atom in (0, len(xyz) - 1):
        for axis in range(3):
            energies = []
            for sign in (-1, 1):
                x = xyz.copy()
                x[atom, axis] += sign * h
                ctx.setPositions(x)
                energies.append(
                    ctx.getState(getEnergy=True)
                    .getPotentialEnergy()
                    .value_in_unit(unit.kilojoule_per_mole)
                )
            fd = -(energies[1] - energies[0]) / (2 * h) / 10
            errors.append(abs(fd - force[atom, axis]))
    del ctx, integrator
    directory.mkdir(parents=True, exist_ok=False)
    data = export_lammps_data(backend)
    script = export_lammps_input(backend)
    (directory / "system.data").write_text(data)
    (directory / "in.periodic").write_text(script)
    proc = subprocess.run(
        [str(executable), "-in", "in.periodic"],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    (directory / "execution.log").write_text(proc.stdout + proc.stderr)
    result = {
        "backend_identity": backend.payload["identity"],
        "final_graph_identity": backend.payload["final_graph_identity"],
        "source": backend.payload["config"]["source"],
        "configuration": backend.payload["config"],
        "openmm_energy_kj_mol": energy,
        "openmm_components_kj_mol": components,
        "finite_difference_max_error_kj_mol_angstrom": float(max(errors)),
        "lammps_exit_code": proc.returncode,
    }
    if proc.returncode:
        result["passed"] = False
        return result
    lines = proc.stdout.splitlines()
    header = next(n for n, l in enumerate(lines) if "PotEng" in l)
    values = list(map(float, lines[header + 1].split()))
    lammps = dict(
        zip(
            (
                "energy",
                "bond",
                "angle",
                "proper",
                "improper",
                "lj",
                "coulomb_real",
                "coulomb_reciprocal",
            ),
            [v * 4.184 for v in values[1:]],
            strict=True,
        )
    )
    lines = (directory / "forces.dump").read_text().splitlines()
    start = next(n for n, l in enumerate(lines) if l.startswith("ITEM: ATOMS")) + 1
    rows = np.array([[float(v) for v in row.split()] for row in lines[start:]])
    lf = rows[:, -3:] * 4.184
    expected = {
        k: 0.0 for k in ("bond", "angle", "proper", "improper", "lj", "coulomb")
    }
    for key, value in components.items():
        if key in {"quartic_bond"}:
            name = "bond"
        elif key in {"quartic_angle", "bond-bond", "bond-angle"}:
            name = "angle"
        elif key == "angle-angle":
            name = "improper"
        elif key in {"lj", "coulomb", "bond", "angle", "proper", "improper"}:
            name = key
        else:
            name = "proper"
        expected[name] += value
    reference = {k: lammps[k] for k in ("bond", "angle", "proper", "improper", "lj")}
    reference["coulomb"] = lammps["coulomb_real"] + lammps["coulomb_reciprocal"]
    component_errors = {k: abs(reference[k] - expected[k]) for k in expected}
    checks = {
        "energy": abs(energy - lammps["energy"]) <= 1e-5,
        "components": max(component_errors.values()) <= 1e-5,
        "forces": np.allclose(force, lf, atol=1e-5, rtol=2e-10),
        "finite_difference": max(errors) <= 1e-5,
        "charges": np.array_equal(
            rows[:, 2],
            [
                next(
                    s["charge"]
                    for s in backend.payload["numerical_data"]["sites"]
                    if s["id"] == i
                )
                for i in sorted(backend.system.topology.sites)
            ],
        ),
        "atom_order": np.array_equal(rows[:, 0], sorted(backend.system.topology.sites)),
        "box": all(
            np.allclose(
                list(map(float, lines[5 + i].split())),
                [0, backend.payload["config"]["box_lengths"][i]],
                atol=0,
                rtol=0,
            )
            for i in range(3)
        ),
    }
    result.update(
        lammps_components_kj_mol=lammps,
        component_errors_kj_mol=component_errors,
        energy_error_kj_mol=abs(energy - lammps["energy"]),
        force_max_error_kj_mol_angstrom=float(np.max(np.abs(force - lf))),
        checks={k: bool(v) for k, v in checks.items()},
        passed=bool(all(checks.values())),
    )
    return result
