"""Synthetic complete nine-atom Hamiltonian export/run-0 control, never a profile."""

import argparse
import itertools
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hamiltonian import canonical, evaluate, load_model, parts


def run(repo, engine, work):
    import json

    repo, engine, work = Path(repo), Path(engine).resolve(), Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    model = load_model(repo)
    xyz = json.loads(
        (repo / "docs/evidence/phase_4j60/plan/reference-plan.json").read_text()
    )["geometries"]["center"]
    rows = model["interaction_inventory"]
    lookup = {r["id"]: r for r in rows}
    columns = parts(model, xyz)["columns"]
    values = {name: (i + 1) * 0.31 - 2.7 for i, name in enumerate(columns)}
    require(all(v != 0 for v in values.values()), "zero_synthetic_control")
    expected = evaluate(model, xyz, values)

    def row(family, sites):
        key = family + ":" + ",".join(map(str, min(tuple(sites), tuple(sites[::-1]))))
        return lookup[key]

    def eq(family, sites):
        return row(family, sites)["assignment"]["values"][0]

    def coeff(family, sites, block=None, n=None):
        _, types = canonical({"family": family, "sites": sites}, model["typing"])
        key = family + ":" + ",".join(types)
        if block is not None:
            key += f":block{block}:n{n}"
        return values[key]

    def line(command, *numbers):
        return command + " " + " ".join(format(float(n), ".17g") for n in numbers)

    bonds = [r for r in rows if r["family"] == "quartic_bond"]
    angles = [r for r in rows if r["family"] == "quartic_angle"]
    dihedrals = [r for r in rows if r["family"] == "torsion_3"]
    neighbors = {i: [] for i in range(1, 10)}
    for r in bonds:
        i, j = r["sites"]
        neighbors[i].append(j)
        neighbors[j].append(i)
    impropers = [
        [i, j, k, l]
        for j, ns in neighbors.items()
        for i, k, l in itertools.combinations(sorted(ns), 3)
    ]
    data = [
        "J60 synthetic full-graph static control",
        "",
        "9 atoms",
        f"{len(bonds)} bonds",
        f"{len(angles)} angles",
        f"{len(dihedrals)} dihedrals",
        f"{len(impropers)} impropers",
        "",
        "3 atom types",
        f"{len(bonds)} bond types",
        f"{len(angles)} angle types",
        f"{len(dihedrals)} dihedral types",
        f"{len(impropers)} improper types",
        "",
        "-60 60 xlo xhi",
        "-60 60 ylo yhi",
        "-60 60 zlo zhi",
        "",
        "Masses",
        "",
    ]
    labels = ["sio", "osi", "hsi"]
    for t, label in enumerate(labels, 1):
        atom = next(i for i in range(1, 10) if model["typing"][str(i)] == label)
        data.append(line(str(t), lookup[f"mass:{atom}"]["assignment"]["value"]))
    data += ["", "Atoms # full", ""]
    for i, x in enumerate(xyz, 1):
        data.append(
            line(
                f"{i} 1 {labels.index(model['typing'][str(i)]) + 1}",
                model["charge_vector"][str(i)],
                *x,
            )
        )
    for section, entries in (("Bonds", bonds), ("Angles", angles)):
        data += ["", section, ""]
        data += [
            f"{i} {i} " + " ".join(map(str, r["sites"]))
            for i, r in enumerate(entries, 1)
        ]
    ordered = [canonical(r, model["typing"])[0] for r in dihedrals]
    for section, entries in (("Dihedrals", ordered), ("Impropers", impropers)):
        data += ["", section, ""]
        data += [f"{i} {i} " + " ".join(map(str, s)) for i, s in enumerate(entries, 1)]
    (work / "control.data").write_text("\n".join(data) + "\n")
    commands = [
        "units real",
        "atom_style full",
        "boundary f f f",
        "bond_style class2",
        "angle_style class2",
        "dihedral_style class2",
        "improper_style class2",
        "read_data control.data",
        "pair_style lj/class2/coul/cut 100",
        "dielectric 1",
        "pair_modify mix sixthpower shift no tail no",
        "special_bonds lj/coul 0 0 1 angle no dihedral no",
    ]
    for t, label in enumerate(labels, 1):
        atom = next(i for i in range(1, 10) if model["typing"][str(i)] == label)
        radius, epsilon = lookup[f"nonbond(9-6):{atom}"]["assignment"]["values"]
        commands.append(line(f"pair_coeff {t} {t}", epsilon / 4.184, radius))
    for i, r in enumerate(bonds, 1):
        equilibrium, *ks = r["assignment"]["values"]
        commands.append(line(f"bond_coeff {i}", equilibrium, *(k / 4.184 for k in ks)))
    for i, r in enumerate(angles, 1):
        s = r["sites"]
        theta, *ks = r["assignment"]["values"]
        r1, r2 = eq("quartic_bond", s[:2]), eq("quartic_bond", s[1:])
        bb = row("bond-bond", s)["assignment"]["values"][0]
        ba = row("bond-angle", s)["assignment"]["values"]
        commands += [
            line(f"angle_coeff {i}", np.degrees(theta), *(k / 4.184 for k in ks)),
            line(f"angle_coeff {i} bb", bb / 4.184, r1, r2),
            line(f"angle_coeff {i} ba", *(k / 4.184 for k in ba), r1, r2),
        ]
    for i, (r, s) in enumerate(zip(dihedrals, ordered, strict=True), 1):
        ks = r["assignment"]["values"]
        commands.append(
            line(
                f"dihedral_coeff {i}",
                *[
                    v
                    for n in range(3)
                    for v in (ks[2 * n] / 4.184, np.degrees(ks[2 * n + 1]) + 180)
                ],
            )
        )
        rs = [eq("quartic_bond", s[j : j + 2]) for j in range(3)]
        ts = [np.degrees(eq("quartic_angle", s[j : j + 3])) for j in range(2)]
        for flag, family, blocks, equilibria in (
            ("mbt", "middle_bond-torsion_3", 1, [rs[1]]),
            ("ebt", "end_bond-torsion_3", 2, [rs[0], rs[2]]),
            ("at", "angle-torsion_3", 2, ts),
        ):
            commands.append(
                line(
                    f"dihedral_coeff {i} {flag}",
                    *[
                        coeff(family, s, block, n)
                        for block in range(blocks)
                        for n in (1, 2, 3)
                    ],
                    *equilibria,
                )
            )
        commands += [
            line(f"dihedral_coeff {i} aat", coeff("angle-angle-torsion_1", s), *ts),
            line(f"dihedral_coeff {i} bb13", coeff("bond-bond_1_3", s), rs[0], rs[2]),
        ]
    for number, (i, j, k, l) in enumerate(impropers, 1):
        # Explicit AA-only carrier for structurally absent Wilson interaction;
        # K=0 is an encoding of declared model inapplicability, not a missing row.
        commands.append(f"improper_coeff {number} 0 0")
        commands.append(
            line(
                f"improper_coeff {number} aa",
                coeff("angle-angle", [i, j, k, l]),
                coeff("angle-angle", [k, j, i, l]),
                coeff("angle-angle", [i, j, l, k]),
                *(
                    np.degrees(eq("quartic_angle", s))
                    for s in ([i, j, k], [i, j, l], [k, j, l])
                ),
            )
        )
    commands += [
        "thermo 1",
        "run 0",
        'print "J60_E $(pe:%.17g)"',
        "write_dump all custom forces.dump id fx fy fz modify sort id format float %.17g",
    ]
    (work / "control.in").write_text("\n".join(commands) + "\n")
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
    require(p.returncode == 0, "full_graph_control_failed")
    energy = float(
        re.findall(
            r"^J60_E ([-+\d.eE]+)", (work / "stdout.log").read_text(), re.MULTILINE
        )[-1]
    )
    forces = np.array(
        [
            [float(v) for v in s.split()[1:]]
            for s in (work / "forces.dump").read_text().splitlines()[-9:]
        ]
    )
    de, df = (
        abs(energy - expected["energy"]),
        float(np.max(abs(forces - expected["forces"]))),
    )
    require(
        de < 1e-9 and df < 1e-8, "full_graph_energy_force_mismatch", energy=de, force=df
    )
    return signed(
        {
            "schema": "island.j60-full-graph-synthetic-control.v1",
            "hamiltonian": expected["hamiltonian"],
            "graph_identity": model["graph_identity"],
            "synthetic_coefficients": values,
            "source_physical_coefficients_fitted": False,
            "PCFF_MD_executed": False,
            "static_only": True,
            "atom_count": 9,
            "bond_count": len(bonds),
            "angle_count": len(angles),
            "dihedral_count": len(dihedrals),
            "AA_only_carrier_count": len(impropers),
            "Wilson_carrier_encoding": "Explicit 0 0 encoding of structurally inapplicable Wilson; all24 AA products active and nonzero synthetic coefficients. This does not assign an absent physical Wilson row.",
            "engine_sha256": pin(engine),
            "expected_energy": expected["energy"],
            "native_energy": energy,
            "energy_error": de,
            "force_max_error": df,
            "files": {p.name: pin(p) for p in work.iterdir() if p.is_file()},
        }
    )


def reconstruct(repo, evidence):
    import json

    repo, evidence = Path(repo), Path(evidence)
    result = json.loads((evidence / "backend-full-graph-control.json").read_text())
    root = evidence / "backend-static/full-graph"
    require(
        all(pin(root / p) == h for p, h in result["files"].items()),
        "full_graph_native_mutation",
    )
    xyz = json.loads((evidence / "plan/reference-plan.json").read_text())["geometries"][
        "center"
    ]
    expected = evaluate(load_model(repo), xyz, result["synthetic_coefficients"])
    energy = float(
        re.findall(
            r"^J60_E ([-+\d.eE]+)", (root / "stdout.log").read_text(), re.MULTILINE
        )[-1]
    )
    forces = np.array(
        [
            [float(v) for v in s.split()[1:]]
            for s in (root / "forces.dump").read_text().splitlines()[-9:]
        ]
    )
    require(
        abs(energy - expected["energy"]) == result["energy_error"]
        and float(np.max(abs(forces - expected["forces"]))) == result["force_max_error"]
        and result["energy_error"] < 1e-9
        and result["force_max_error"] < 1e-8,
        "full_graph_reconstruction_mismatch",
    )
    return {
        "energy_error": result["energy_error"],
        "force_max_error": result["force_max_error"],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    write(a.output, run(a.repo, a.engine, a.work))
