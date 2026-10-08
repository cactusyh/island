"""Four-family J22 periodic single points, components, forces and finite differences.

Requires retained real-source bundles, OpenMM, and an explicitly chosen LAMMPS
binary with KSPACE. No unavailable execution is converted into a pass.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from periodic_j22_fixtures import retained_fixture

from island.periodic import (
    build_openmm_periodic_system,
    export_lammps_data,
    export_lammps_input,
)


def evaluate(backend, mm, unit, directory, executable):
    model = build_openmm_periodic_system(backend)
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
        "backend_identity": backend.identity,
        "final_graph_identity": backend.final_graph.identity,
        "source": backend.config.payload["source"],
        "configuration": backend.config.payload,
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
                backend.graph_charges.charges[i]
                for i in sorted(backend.system.topology.sites)
            ],
        ),
        "atom_order": np.array_equal(rows[:, 0], sorted(backend.system.topology.sites)),
        "box": all(
            np.allclose(
                list(map(float, lines[5 + i].split())),
                [0, backend.config.payload["box_lengths"][i]],
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--method", choices=["pme", "ewald"], default="pme")
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    import openmm as mm
    from openmm import unit

    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island_periodic_numerical_verification_v1",
        "lammps_sha256": hashlib.sha256(args.lammps.read_bytes()).hexdigest(),
        "openmm_version": mm.__version__,
        "energy_atol_kj_mol": 1e-5,
        "force_atol_kj_mol_angstrom": 1e-5,
        "force_rtol": 2e-10,
        "cases": {},
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    for family in ("pcff", "oplsaa", "gaff", "gaff2"):
        backend, _ = retained_fixture(
            family, method=args.method, tolerance=args.tolerance
        )
        result = evaluate(
            backend, mm, unit, args.output / family, args.lammps.resolve()
        )
        report["cases"][family] = result
        (args.output / "receipt.json").write_text(
            json.dumps(report, sort_keys=True, indent=2) + "\n"
        )
        print(
            family,
            json.dumps(result.get("checks", {"lammps_failed": True})),
            flush=True,
        )
    return 0 if all(r["passed"] for r in report["cases"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
