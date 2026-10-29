"""Offline CHARMM-LJ/periodic-Coulomb software controls; never a PDMS profile.

The eight artificial particles and their coefficients are mathematical test
inputs, not a physical model or missing PDMS parameters. This command takes no
external parameter data and cannot execute a blocked force-field preparation.
"""

import argparse
import itertools
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.pdms_literature import pinned_bytes

# Fixed before comparison; tighter reciprocal controls establish numerical error.
BOUNDS = {
    "author_accuracy_energy_kj_mol": 0.03,
    "author_accuracy_force_kj_mol_angstrom": 0.01,
    "tight_energy_kj_mol": 0.0003,
    "tight_force_kj_mol_angstrom": 0.0003,
    "finite_difference_force_kj_mol_angstrom": 0.0001,
}
GROUPS = ("bond", "angle", "proper", "lj", "coulomb")
L = 34.0
XYZ = np.array(
    [
        [1, 1, 1],
        [4, 2, 1],
        [7, 1, 2],
        [10.6, 2, 2],
        [17, 14, 12],
        [19, 15, 11],
        [22, 14, 12],
        [25, 16, 13],
    ],
    dtype=float,
)
Q = np.array([0.3, -0.4, 0.2, -0.1] * 2)
SIGMA = np.array([2.0, 2.3, 2.6, 2.9] * 2)
EPS = np.array([0.01, 0.02, 0.03, 0.04] * 2)
BONDS = [(i, i + 1) for start in (0, 4) for i in range(start, start + 3)]
ANGLES = [(i, i + 1, i + 2) for start in (0, 4) for i in range(start, start + 2)]
TORSIONS = [(0, 1, 2, 3), (4, 5, 6, 7)]
# LAMMPS real-unit electrostatic constant; retain difference from OpenMM in report.
COULOMB = 332.06371 * 4.184


def specials():
    return {
        (i, j): (0.5 if j - i == 3 else 0.0)
        for start in (0, 4)
        for i in range(start, start + 4)
        for j in range(i + 1, start + 4)
    }


def switch(r):
    if r <= 9:
        return 1.0
    if r >= 11:
        return 0.0
    return (121 - r * r) ** 2 * (121 + 2 * r * r - 243) / 40**3


def independent_energy(xyz, *, alpha=0.22, reciprocal_shell=12):
    """Direct Ewald sum + explicit equations, independent of either backend."""
    from scipy.special import erfc

    values = {k: 0.0 for k in GROUPS}
    for a, b in BONDS:
        values["bond"] += 4.184 * 10 * (np.linalg.norm(xyz[a] - xyz[b]) - 3) ** 2
    for a, b, c in ANGLES:
        u, v = xyz[a] - xyz[b], xyz[c] - xyz[b]
        theta = math.acos(
            float(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))
        )
        values["angle"] += 4.184 * 2 * (theta - math.radians(110)) ** 2
    for a, b, c, d in TORSIONS:
        n1 = np.cross(xyz[b] - xyz[a], xyz[c] - xyz[b])
        n2 = np.cross(xyz[c] - xyz[b], xyz[d] - xyz[c])
        phi = math.acos(
            float(
                np.clip(np.dot(n1, n2) / np.linalg.norm(n1) / np.linalg.norm(n2), -1, 1)
            )
        )
        values["proper"] += (
            4.184
            * 0.5
            * sum(
                k * (1 + (-1) ** n * math.cos((n + 1) * phi))
                for n, k in enumerate([0.1, 0.2, 0.3, 0.4])
            )
        )
    sp = specials()
    for i, j in itertools.combinations(range(8), 2):
        delta = xyz[i] - xyz[j]
        delta -= L * np.rint(delta / L)
        r = np.linalg.norm(delta)
        sigma = math.sqrt(SIGMA[i] * SIGMA[j])
        values["lj"] += (
            sp.get((i, j), 1.0)
            * 4
            * 4.184
            * math.sqrt(EPS[i] * EPS[j])
            * ((sigma / r) ** 12 - (sigma / r) ** 6)
            * switch(r)
        )
    n = np.array(list(itertools.product(range(-2, 3), repeat=3)))
    real = 0.0
    for i in range(8):
        for j in range(8):
            r = np.linalg.norm(xyz[i] - xyz[j] + L * n, axis=1)
            r = r[r > 1e-12]
            real += 0.5 * Q[i] * Q[j] * np.sum(erfc(alpha * r) / r)
    h = np.array(
        [
            v
            for v in itertools.product(
                range(-reciprocal_shell, reciprocal_shell + 1), repeat=3
            )
            if v != (0, 0, 0)
        ]
    )
    k = 2 * math.pi / L * h
    k2 = np.sum(k * k, axis=1)
    structure = np.exp(1j * (k @ xyz.T)) @ Q
    reciprocal = (
        2
        * math.pi
        / L**3
        * np.sum(np.exp(-k2 / (4 * alpha**2)) / k2 * np.abs(structure) ** 2)
    )
    self_energy = -alpha / math.sqrt(math.pi) * np.dot(Q, Q)
    background = -math.pi * np.sum(Q) ** 2 / (2 * alpha**2 * L**3)
    correction = 0.0
    for (i, j), scale in sp.items():
        delta = xyz[i] - xyz[j]
        delta -= L * np.rint(delta / L)
        correction -= (1 - scale) * Q[i] * Q[j] / np.linalg.norm(delta)
    values["coulomb"] = COULOMB * (
        real + reciprocal + self_energy + background + correction
    )
    return values


def openmm_context(accuracy):
    import openmm as mm

    model = mm.System()
    for _ in range(8):
        model.addParticle(12)
    model.setDefaultPeriodicBoxVectors(
        *[mm.Vec3(*(L * 0.1 if j == i else 0 for j in range(3))) for i in range(3)]
    )
    bond = mm.HarmonicBondForce()
    for a, b in BONDS:
        bond.addBond(a, b, 0.3, 10 * 836.8)
    angle = mm.HarmonicAngleForce()
    for a, b, c in ANGLES:
        angle.addAngle(a, b, c, math.radians(110), 2 * 8.368)
    proper = mm.CustomTorsionForce(
        "0.5*(0.4184*(1+cos(theta))+0.8368*(1-cos(2*theta))+1.2552*(1+cos(3*theta))+1.6736*(1-cos(4*theta)))"
    )
    for t in TORSIONS:
        proper.addTorsion(*t, [])
    # CHARMM potential switch is a polynomial in r^2, NOT OpenMM's built-in switch.
    switching = "select(step(r-0.9),(1.21-r*r)^2*(1.21+2*r*r-3*0.81)/(1.21-0.81)^3,1)"
    lj = mm.CustomNonbondedForce(
        "4*sqrt(e1*e2)*((sqrt(s1*s2)/r)^12-(sqrt(s1*s2)/r)^6)*" + switching
    )
    lj.addPerParticleParameter("s")
    lj.addPerParticleParameter("e")
    lj.setNonbondedMethod(lj.CutoffPeriodic)
    lj.setCutoffDistance(1.1)
    lj.setUseLongRangeCorrection(False)
    coul = mm.NonbondedForce()
    coul.setNonbondedMethod(coul.PME)
    coul.setCutoffDistance(1.1)
    coul.setEwaldErrorTolerance(accuracy)
    coul.setUseDispersionCorrection(False)
    coul.setExceptionsUsePeriodicBoundaryConditions(True)
    for i in range(8):
        lj.addParticle([SIGMA[i] * 0.1, EPS[i] * 4.184])
        coul.addParticle(Q[i], 1.0, 0)
    lj14 = mm.CustomBondForce("step(1.1-r)*4*e*((s/r)^12-(s/r)^6)*" + switching)
    lj14.addPerBondParameter("s")
    lj14.addPerBondParameter("e")
    lj14.setUsesPeriodicBoundaryConditions(True)
    for (i, j), scale in specials().items():
        lj.addExclusion(i, j)
        coul.addException(i, j, scale * Q[i] * Q[j], 1.0, 0)
        if scale:
            lj14.addBond(
                i,
                j,
                [
                    math.sqrt(SIGMA[i] * SIGMA[j]) * 0.1,
                    scale * math.sqrt(EPS[i] * EPS[j]) * 4.184,
                ],
            )
    for group, force in enumerate([bond, angle, proper, lj, coul]):
        force.setForceGroup(group)
        model.addForce(force)
    lj14.setForceGroup(3)
    model.addForce(lj14)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    return context, integrator


def mm_evaluate(context, xyz):
    from openmm import unit

    context.setPositions(xyz * 0.1)
    components = {
        name: context.getState(getEnergy=True, groups={g})
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
        for g, name in enumerate(GROUPS)
    }
    state = context.getState(getForces=True)
    forces = (
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        * 0.1
    )
    return components, np.array(forces)


def write_control(directory):
    lines = [
        "ISLAND artificial longrange software control; not PDMS",
        "",
        "8 atoms",
        "6 bonds",
        "4 angles",
        "2 dihedrals",
        "0 impropers",
        "",
        "4 atom types",
        "1 bond types",
        "1 angle types",
        "1 dihedral types",
        "",
        f"0 {L} xlo xhi",
        f"0 {L} ylo yhi",
        f"0 {L} zlo zhi",
        "",
        "Masses",
        "",
    ]
    lines += [f"{i + 1} 12" for i in range(4)]
    lines += ["", "Atoms # full", ""]
    lines += [
        f"{i + 1} {i // 4 + 1} {i % 4 + 1} {Q[i]:.17g} "
        + " ".join(format(x, ".17g") for x in XYZ[i])
        for i in range(8)
    ]
    for name, terms in [("Bonds", BONDS), ("Angles", ANGLES), ("Dihedrals", TORSIONS)]:
        lines += ["", name, ""] + [
            f"{n + 1} 1 " + " ".join(str(i + 1) for i in t) for n, t in enumerate(terms)
        ]
    (directory / "control.data").write_text("\n".join(lines) + "\n")


def lammps_evaluate(executable, directory, accuracy, *, table_zero=True):
    commands = [
        "units real",
        "atom_style full",
        "boundary p p p",
        "bond_style harmonic",
        "angle_style harmonic",
        "dihedral_style opls",
        "pair_style lj/charmm/coul/long 9.0 11.0",
        "pair_modify mix geometric" + (" table 0" if table_zero else ""),
        "special_bonds lj/coul 0.0 0.0 0.5",
        "read_data control.data",
        "bond_coeff 1 10 3",
        "angle_coeff 1 2 110",
        "dihedral_coeff 1 0.1 0.2 0.3 0.4",
    ]
    commands += [f"pair_coeff {i + 1} {i + 1} {EPS[i]} {SIGMA[i]}" for i in range(4)]
    commands += [
        f"kspace_style pppm {accuracy}",
        "thermo_style custom step pe ebond eangle edihed evdwl ecoul elong",
        "thermo_modify format float %.17g",
        "run 0",
        "write_dump all custom forces.dump id fx fy fz modify sort id format float %.17g",
    ]
    (directory / "control.in").write_text("\n".join(commands) + "\n")
    result = subprocess.run(
        [str(executable), "-in", "control.in", "-log", "log.lammps", "-nocite"],
        cwd=directory,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    (directory / "stdout.txt").write_text(result.stdout + result.stderr)
    require(result.returncode == 0, "lammps_control_failed:" + result.stdout[-1000:])
    rows = result.stdout.splitlines()
    index = next(
        n
        for n, line in enumerate(rows)
        if line.split()
        == [
            "Step",
            "PotEng",
            "E_bond",
            "E_angle",
            "E_dihed",
            "E_vdwl",
            "E_coul",
            "E_long",
        ]
    )
    values = list(map(float, rows[index + 1].split()))
    components = dict(
        zip(
            GROUPS,
            [values[2], values[3], values[4], values[5], values[6] + values[7]],
            strict=True,
        )
    )
    components = {k: v * 4.184 for k, v in components.items()}
    dump = (directory / "forces.dump").read_text().splitlines()
    start = next(n for n, l in enumerate(dump) if l.startswith("ITEM: ATOMS")) + 1
    force = (
        np.array(
            [[float(x) for x in line.split()[1:]] for line in dump[start : start + 8]]
        )
        * 4.184
    )
    return components, force


def qualify(executable, expected_sha256, output):
    executable = Path(executable).resolve()
    pinned_bytes(executable, expected_sha256)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    help_result = subprocess.run(
        [str(executable), "-h"], capture_output=True, text=True, timeout=30, check=False
    )
    require(help_result.returncode == 0, "lammps_capability_query_failed")
    for style in ["lj/charmm/coul/long", "pppm", "opls", "harmonic"]:
        require(style in help_result.stdout.split(), "backend_style_missing:" + style)
    (output / "lammps-help.txt").write_text(help_result.stdout)
    truth = independent_energy(XYZ)
    second = independent_energy(XYZ, alpha=0.26, reciprocal_shell=14)
    direct_convergence = max(abs(truth[k] - second[k]) for k in GROUPS)
    require(direct_convergence < 1e-8, "direct_Ewald_convergence_failed")
    cases = []
    for name, accuracy in [
        ("author_accuracy", 1e-4),
        ("intermediate", 1e-7),
        ("tight", 1e-10),
    ]:
        directory = output / name
        directory.mkdir()
        write_control(directory)
        la, lf = lammps_evaluate(
            executable, directory, accuracy, table_zero=name != "author_accuracy"
        )
        context, integrator = openmm_context(
            min(accuracy, 1e-4) if name == "author_accuracy" else 1e-7
        )
        mm, mf = mm_evaluate(context, XYZ)
        bound_name = "author_accuracy" if name == "author_accuracy" else "tight"
        e_bound = BOUNDS[bound_name + "_energy_kj_mol"]
        f_bound = BOUNDS[bound_name + "_force_kj_mol_angstrom"]
        errors = {k: abs(mm[k] - la[k]) for k in GROUPS}
        direct_errors = {
            k: {"LAMMPS": abs(la[k] - truth[k]), "OpenMM": abs(mm[k] - truth[k])}
            for k in GROUPS
        }
        force_error = float(np.max(np.abs(mf - lf)))
        passed = bool(
            max(errors.values()) < e_bound
            and force_error < f_bound
            and max(v for row in direct_errors.values() for v in row.values()) < e_bound
        )
        if name == "tight":
            require(passed, "tight_backend_comparison_failed")
        fd_errors = {}
        if name == "tight":
            for step in (1e-4, 5e-5):
                fd = np.zeros((8, 3))
                for i in range(8):
                    for j in range(3):
                        plus, minus = XYZ.copy(), XYZ.copy()
                        plus[i, j] += step
                        minus[i, j] -= step
                        fd[i, j] = -(
                            sum(independent_energy(plus).values())
                            - sum(independent_energy(minus).values())
                        ) / (2 * step)
                fd_errors[str(step)] = {
                    "OpenMM": float(np.max(np.abs(fd - mf))),
                    "LAMMPS": float(np.max(np.abs(fd - lf))),
                }
                require(
                    max(fd_errors[str(step)].values())
                    < BOUNDS["finite_difference_force_kj_mol_angstrom"],
                    "independent_force_difference_failed",
                )
        cases.append(
            {
                "name": name,
                "LAMMPS_electrostatic_accuracy": accuracy,
                "OpenMM_electrostatic_accuracy": accuracy
                if name == "author_accuracy"
                else 1e-7,
                "within_predeclared_bounds": passed,
                "status": "verified" if passed else "scientifically_unqualified",
                "LAMMPS_components": la,
                "OpenMM_components": mm,
                "independent_direct_Ewald_components": truth,
                "component_absolute_errors": errors,
                "errors_against_independent_equations": direct_errors,
                "max_force_difference": force_error,
                "finite_difference_errors": fd_errors,
            }
        )
        del context, integrator
    report = {
        "schema": "island_charmm_longrange_backend_qualification_control_v1",
        "scope": "eight artificial particles; software equations only; NOT a PDMS execution",
        "LAMMPS_executable_sha256": expected_sha256,
        "LAMMPS_executable": str(executable),
        "OpenMM_version": __import__("openmm").__version__,
        "CHARMM_LJ_switch": "polynomial_in_r_squared; 9 to 11 angstrom; no dispersion tail",
        "electrostatics": "LAMMPS PPPM / OpenMM PME / independent direct Ewald; 3D periodic tinfoil",
        "special_bonds": {"lj": [0, 0, 0.5], "coulomb": [0, 0, 0.5]},
        "mixing": "geometric_sigma_and_epsilon",
        "LAMMPS_pair_table": "author_accuracy uses default table; intermediate/tight disable table for independent equation control",
        "LAMMPS_coulomb_constant_kj_angstrom": COULOMB,
        "OpenMM_coulomb_constant_kj_angstrom": 1389.35457644382,
        "bounds_predeclared": BOUNDS,
        "independent_Ewald_alpha_and_shell_convergence": direct_convergence,
        "controls": cases,
        "software_equations_verified": True,
        "author_default_accuracy_control_verified": cases[0][
            "within_predeclared_bounds"
        ],
        "physical_PDMS_executable": False,
        "PDMS_MD_run": False,
        "execution_profile": None,
        "author_model_qualified": False,
        "periodic_scientific_closure": False,
        "limitations": LIMITS,
    }
    report["identity"] = digest(report)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--executable-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = qualify(args.executable, args.executable_sha256, args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
