"""Predeclared isolated Class II comparisons against an actual LAMMPS executable.

All fixtures are synthetic numerical tests, not chemical parameter references.
No production model conversion or interaction-generation helper builds the oracle.
"""

import argparse
import json
import math
import subprocess
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np

from island.forcefields.pcff.terms import (
    Class2Term,
    finite_difference_forces,
    pair_energy,
)
from island.workflows.storage import json_bytes, publish

XYZ = [[0.1, 1.3, 0.2], [0.0, 0.0, 0.0], [1.5, 0.2, 0.1], [2.0, 1.1, 1.4]]
# Raw kcal coefficients and degree equilibria are hand-authored, independent of FRC.
FIXTURES = [
    ("quartic_bond", [1.1, 2.3, -0.7, 0.4], []),
    ("quartic_angle", [107, 2.3, -0.7, 0.4], []),
    ("bond-bond", [1.7], [1.1, 1.4]),
    ("bond-angle", [1.7, -0.8], [1.1, 1.4, 107]),
    ("torsion_3", [1.2, 17, -0.7, -31, 0.4, 54], []),
    ("middle_bond-torsion_3", [1.1, -0.5, 0.3], [1.4]),
    ("end_bond-torsion_3", [1.1, -0.5, 0.3, -0.2, 0.7, -0.4], [1.1, 1.6]),
    ("angle-torsion_3", [1.1, -0.5, 0.3, -0.2, 0.7, -0.4], [107, 119]),
    ("angle-angle-torsion_1", [-1.3], [107, 119]),
    ("bond-bond_1_3", [1.8], [1.1, 1.6]),
    ("angle-angle", [1.1], [101, 107]),
]
STEPS = [1e-3, 3e-4, 1e-4]
TOLERANCES = {
    "energy_atol_kj_mol": 1e-10,
    "energy_rtol": 1e-12,
    "force_atol_kj_mol_angstrom": 3e-5,
    "force_rtol": 1e-7,
}
PAIR_WEIGHTS = {"lj": [0.25, 0.5, 0.75], "coulomb": [0.2, 0.4, 0.6]}


def fmt(values):
    return " ".join(format(v, ".17g") for v in values)


def prepare(family, p, q, reverse=False):
    """Independent handwritten LAMMPS coefficient commands plus kernel fixture."""
    if family == "quartic_bond":
        n = 2
        kind = "bond"
        commands = ["bond_coeff 1 " + fmt(p)]
        coeff = [p[0]] + [v * 4.184 for v in p[1:]]
        eq = []
    elif family in ("quartic_angle", "bond-bond", "bond-angle"):
        n = 3
        kind = "angle"
        commands = [
            "angle_coeff 1 107 0 0 0",
            "angle_coeff 1 bb 0 1.1 1.4",
            "angle_coeff 1 ba 0 0 1.1 1.4",
        ]
        if family == "quartic_angle":
            commands[0] = "angle_coeff 1 " + fmt(p)
            coeff = [math.radians(p[0])] + [v * 4.184 for v in p[1:]]
            eq = []
        elif family == "bond-bond":
            commands[1] = "angle_coeff 1 bb " + fmt(p + q)
            coeff = [v * 4.184 for v in p]
            eq = q[:]
        else:
            commands[2] = "angle_coeff 1 ba " + fmt(p + q[:2])
            coeff = [v * 4.184 for v in p]
            eq = q[:2] + [math.radians(q[2])]
    elif family == "angle-angle":
        n = 4
        kind = "improper"
        commands = [
            "improper_coeff 1 0 0",
            f"improper_coeff 1 aa {p[0]} 0 0 {q[0]} 113 {q[1]}",
        ]
        coeff = [p[0] * 4.184]
        eq = list(map(math.radians, q))
    else:
        n = 4
        kind = "dihedral"
        commands = [
            "dihedral_coeff 1 0 0 0 0 0 0",
            "dihedral_coeff 1 mbt 0 0 0 1.4",
            "dihedral_coeff 1 ebt 0 0 0 0 0 0 1.1 1.6",
            "dihedral_coeff 1 at 0 0 0 0 0 0 107 119",
            "dihedral_coeff 1 aat 0 107 119",
            "dihedral_coeff 1 bb13 0 1.1 1.6",
        ]
        if family == "torsion_3":
            commands[0] = "dihedral_coeff 1 " + fmt(p)
            coeff = [
                v * 4.184 if i % 2 == 0 else math.radians(v) for i, v in enumerate(p)
            ]
            eq = []
        else:
            token, index = {
                "middle_bond-torsion_3": ("mbt", 1),
                "end_bond-torsion_3": ("ebt", 2),
                "angle-torsion_3": ("at", 3),
                "angle-angle-torsion_1": ("aat", 4),
                "bond-bond_1_3": ("bb13", 5),
            }[family]
            commands[index] = "dihedral_coeff 1 " + token + " " + fmt(p + q)
            coeff = [v * 4.184 for v in p]
            eq = list(map(math.radians, q)) if token in ("at", "aat") else q[:]
    indices = list(range(n))
    if reverse:
        indices.reverse()
    term = Class2Term(family, tuple(indices), tuple(coeff), tuple(eq))
    return XYZ[:n], [term], kind, [indices], commands


def write_data(path, xyz, kind, interactions, pair=False):
    counts = {
        k: len(interactions) if k == kind else 0
        for k in ("bond", "angle", "dihedral", "improper")
    }
    text = "ISLAND synthetic isolated-term fixture\n\n" + str(len(xyz)) + " atoms\n"
    text += "".join(f"{n} {k}s\n" for k, n in counts.items())
    text += f"{3 if pair else 1} atom types\n" + "".join(
        f"{int(n > 0)} {k} types\n" for k, n in counts.items()
    )
    text += "\n-50 50 xlo xhi\n-50 50 ylo yhi\n-50 50 zlo zhi\n\nMasses\n\n1 1\n"
    if pair:
        text += "2 1\n3 1\n"
    text += "\nAtoms # full\n\n"
    for i, pos in enumerate(xyz):
        typ = 1 if not pair or i == 0 else 2 if i == len(xyz) - 1 else 3
        charge = (0.17 if i == 0 else -0.17 if i == len(xyz) - 1 else 0) if pair else 0
        text += f"{i + 1} 1 {typ} {charge} " + fmt(pos) + "\n"
    if interactions:
        text += "\n" + kind.capitalize() + "s\n\n"
        for i, indices in enumerate(interactions, 1):
            text += f"{i} 1 " + " ".join(str(j + 1) for j in indices) + "\n"
    path.write_text(text)


def run_lammps(
    exe, directory, xyz, kind, interactions, commands, *, pair=False, weights=None
):
    write_data(directory / "data.lmp", xyz, kind, interactions, pair)
    pre = f"{kind}_style class2\n" if interactions else ""
    if pair:
        pre += "pair_style lj/class2/coul/cut 30.0\n"
        pre += (
            "special_bonds lj "
            + fmt(weights["lj"])
            + " coul "
            + fmt(weights["coulomb"])
            + "\n"
        )
    else:
        pre += "pair_style zero 30.0\nspecial_bonds lj/coul 1 1 1\n"
        commands = ["pair_coeff * *"] + commands
    script = (
        "units real\natom_style full\nboundary f f f\n"
        + pre
        + "read_data data.lmp\n"
        + "\n".join(commands)
    )
    script += "\npair_modify shift no tail no\nthermo 1\nthermo_style custom step pe ebond eangle edihed eimp evdwl ecoul\nthermo_modify format float %.17g\nrun 0\n"
    script += 'print "$(pe:%.17g) $(evdwl:%.17g) $(ecoul:%.17g)" file energy.txt\nwrite_dump all custom forces.dump id fx fy fz modify sort id format float %.17g\n'
    (directory / "in.lmp").write_text(script)
    done = subprocess.run(
        [str(exe), "-in", "in.lmp", "-log", "log.lammps"],
        cwd=directory,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    (directory / "stdout.txt").write_text(done.stdout)
    (directory / "stderr.txt").write_text(done.stderr)
    if done.returncode:
        raise RuntimeError(
            f"LAMMPS exit {done.returncode}: {done.stdout[-1200:]} {done.stderr}"
        )
    energy = (
        np.array(list(map(float, (directory / "energy.txt").read_text().split())))
        * 4.184
    )
    lines = (directory / "forces.dump").read_text().splitlines()
    start = (
        next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS")) + 1
    )
    forces = (
        np.array([[float(x) for x in line.split()[1:]] for line in lines[start:]])
        * 4.184
    )
    return energy, forces


def check(energy, forces, expected, force_series):
    errors = [float(np.max(np.abs(f - forces))) for f in force_series]
    assert np.isclose(
        energy,
        expected,
        atol=TOLERANCES["energy_atol_kj_mol"],
        rtol=TOLERANCES["energy_rtol"],
    ), (energy, expected)
    assert np.allclose(
        force_series[-1],
        forces,
        atol=TOLERANCES["force_atol_kj_mol_angstrom"],
        rtol=TOLERANCES["force_rtol"],
    ), errors
    assert errors[-1] < errors[0] * 0.2 or errors[-1] < 1e-8, errors
    return {
        "energy_error_kj_mol": float(abs(energy - expected)),
        "force_max_errors_kj_mol_angstrom": errors,
    }


def term_case(exe, directory, family, p, q, reverse=False):
    xyz, terms, kind, interactions, commands = prepare(family, p, q, reverse)
    energy, forces = run_lammps(exe, directory, xyz, kind, interactions, commands)
    expected = sum(t.energy({i: xyz[i] for i in t.sites}) for t in terms)
    series = []
    for h in STEPS:
        f = np.zeros((len(xyz), 3))
        for t in terms:
            for i, v in finite_difference_forces(
                t, {i: xyz[i] for i in t.sites}, displacement=h
            ).items():
                f[i] += v
        series.append(f)
    return check(energy[0], forces, expected, series)


def tetrahedral_case(exe, directory):
    xyz = [
        [1.3, 0.2, 0.1],
        [0, 0, 0],
        [-0.2, 1.1, 0.3],
        [-0.4, -0.5, 1.4],
        [-1, -0.7, -0.8],
    ]
    quads = [(i, 1, k, l) for i, k, l in combinations((0, 2, 3, 4), 3)]
    terms = []
    for i, j, k, l in quads:
        for sites, value, degrees in [
            ((i, j, k, l), 1.1, (101, 107)),
            ((k, j, i, l), -0.4, (101, 113)),
            ((i, j, l, k), 0.7, (113, 107)),
        ]:
            terms.append(
                Class2Term(
                    "angle-angle",
                    sites,
                    (value * 4.184,),
                    tuple(map(math.radians, degrees)),
                )
            )
    commands = ["improper_coeff 1 0 0", "improper_coeff 1 aa 1.1 -0.4 0.7 101 113 107"]
    energy, forces = run_lammps(exe, directory, xyz, "improper", quads, commands)
    expected = sum(t.energy({i: xyz[i] for i in t.sites}) for t in terms)
    series = []
    for h in STEPS:
        f = np.zeros((5, 3))
        for t in terms:
            for i, v in finite_difference_forces(
                t, {i: xyz[i] for i in t.sites}, displacement=h
            ).items():
                f[i] += v
        series.append(f)
    return {
        **check(energy[0], forces, expected, series),
        "angle_angle_couplings": len(terms),
    }


def pair_case(exe, directory, path_length, zero_lj=False, excluded=False):
    n = path_length + 1 if path_length else 2
    xyz = np.array(
        [[0, 0, 0], [3.7, 0.2, 0.1], [5.2, 1.1, 0.4], [6.1, 2.0, 1.3]][:n], dtype=float
    )
    weights = {"lj": [0, 0, 0], "coulomb": [0, 0, 0]} if excluded else PAIR_WEIGHTS
    edges = [(i, i + 1) for i in range(n - 1)] if path_length else []
    commands = [
        f"pair_coeff 1 1 {0 if zero_lj else 0.15} 3.2",
        "pair_coeff 2 2 .35 4.1",
        "pair_coeff 3 3 0 2.0",
    ]
    if edges:
        commands.append("bond_coeff 1 1.1 0 0 0")
    energy, forces = run_lammps(
        exe, directory, xyz, "bond", edges, commands, pair=True, weights=weights
    )

    def kernel(x):
        return pair_energy(
            float(np.linalg.norm(x[0] - x[-1])),
            rmin_i=3.2,
            epsilon_i=0 if zero_lj else 0.15 * 4.184,
            rmin_j=4.1,
            epsilon_j=0.35 * 4.184,
            charge_i=0.17,
            charge_j=-0.17,
            lj_weight=weights["lj"][path_length - 1] if path_length else 1,
            coulomb_weight=weights["coulomb"][path_length - 1] if path_length else 1,
        )

    expected = kernel(xyz)
    series = []
    for h in STEPS:
        f = np.zeros_like(xyz)
        for i in range(n):
            for axis in range(3):
                plus = xyz.copy()
                minus = xyz.copy()
                plus[i, axis] += h
                minus[i, axis] -= h
                f[i, axis] = -(kernel(plus)["total"] - kernel(minus)["total"]) / (2 * h)
        series.append(f)
    assert np.allclose(
        energy[1:], [expected["lj_9_6"], expected["coulomb"]], atol=1e-10, rtol=1e-12
    )
    return check(energy[0], forces, expected["total"], series)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    exe = Path(args.lammps).resolve()
    declaration = {
        "schema": "island_pcff_term_acceptance_v1",
        "fixtures": FIXTURES,
        "xyz": XYZ,
        "fd_displacements_angstrom": STEPS,
        "tolerances": TOLERANCES,
        "special_pairs": PAIR_WEIGHTS,
        "special_pairs_origin": "explicit test model choice, not recovered PCFF fact",
        "extra_cases": [
            "all proper terms reversed",
            "tetrahedral 12 angle-angle couplings",
            "pair unbonded/1-2/1-3/1-4",
            "zero LJ charged site",
            "excluded pair",
        ],
        "timeout_seconds_per_case": 60,
        "retries": 0,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(out / "declaration.json", json_bytes(declaration))
    outcomes = []
    try:
        info = subprocess.run(
            [str(exe), "-help"], capture_output=True, text=True, check=True, timeout=30
        )
        publish(
            out / "executable.json",
            json_bytes(
                {
                    "path": str(exe),
                    "sha256": sha256(exe.read_bytes()).hexdigest(),
                    "help": info.stdout,
                    "revision_claim": "see separately verified build provenance; help alone does not prove source revision",
                }
            ),
        )
        cases = []
        for family, p, q in FIXTURES:
            cases.append(
                (family, lambda d, f=family, p=p, q=q: term_case(exe, d, f, p, q))
            )
            if family in (
                "torsion_3",
                "middle_bond-torsion_3",
                "end_bond-torsion_3",
                "angle-torsion_3",
                "angle-angle-torsion_1",
                "bond-bond_1_3",
            ):
                cases.append(
                    (
                        family + "_reversed",
                        lambda d, f=family, p=p, q=q: term_case(exe, d, f, p, q, True),
                    )
                )
        cases.append(("tetrahedral", lambda d: tetrahedral_case(exe, d)))
        for path in range(4):
            cases.append(
                (f"pair_path_{path}", lambda d, path=path: pair_case(exe, d, path))
            )
        cases.extend(
            [
                ("zero_lj_charged", lambda d: pair_case(exe, d, 0, True)),
                ("excluded_pair", lambda d: pair_case(exe, d, 1, excluded=True)),
            ]
        )
        for name, fn in cases:
            directory = out / name
            directory.mkdir()
            try:
                result = {"name": name, "status": "passed", **fn(directory)}
            except Exception as error:  # noqa: BLE001 -- retain every declared failure
                result = {
                    "name": name,
                    "status": "failed",
                    "error": f"{type(error).__name__}: {error}",
                }
            result["files"] = {
                p.name: sha256(p.read_bytes()).hexdigest()
                for p in sorted(directory.iterdir())
                if p.is_file()
            }
            outcomes.append(result)
            print(json.dumps(result), flush=True)
            publish(
                out / "report.json", json_bytes({"outcomes": outcomes}), replace=True
            )
    except (OSError, subprocess.SubprocessError) as error:
        outcomes.append(
            {"name": "executable", "status": "unavailable", "error": str(error)}
        )
        publish(out / "report.json", json_bytes({"outcomes": outcomes}), replace=True)
    return int(not outcomes or any(o["status"] != "passed" for o in outcomes))


if __name__ == "__main__":
    raise SystemExit(main())
