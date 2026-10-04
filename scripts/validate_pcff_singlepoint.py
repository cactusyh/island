"""Declared finite PCFF acceptance against independent pinned msi2lmp + LAMMPS.

Validation-only CAR/MDF writer; not a supported general LAMMPS exporter.
No production model interactions are used to construct the reference.
"""

import argparse
import json
import subprocess
import time
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
from validate_pcff_automatic import expected_labels, reference_charges
from validate_pcff_class2 import CASES, oracle_check

from island.charge_references.records import unpack
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    load_pcff_source,
    save_pcff_model,
    save_pcff_parameters,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.source import PIN
from island.workflows.storage import json_bytes, publish

TOL = {"energy_atol_kj_mol": 1e-5, "force_atol_kj_mol_angstrom": 1e-5, "rtol": 2e-10}
GROUPS = {
    "bond": ["quartic_bond"],
    "angle": ["quartic_angle", "bond-bond", "bond-angle"],
    "dihedral": [
        "torsion_3",
        "middle_bond-torsion_3",
        "end_bond-torsion_3",
        "angle-torsion_3",
        "angle-angle-torsion_1",
        "bond-bond_1_3",
    ],
    "improper": ["angle-angle"],
    "lj": ["lj_9_6"],
    "coulomb": ["coulomb"],
}


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def run(command, directory, label):
    start = time.perf_counter()
    proc = subprocess.run(
        list(map(str, command)),
        cwd=directory,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )
    (directory / (label + ".log")).write_text(proc.stdout + proc.stderr)
    (directory / (label + ".command.json")).write_text(
        json.dumps(list(map(str, command)))
    )
    if proc.returncode:
        raise ValueError(f"{label} exit {proc.returncode}; see retained log")
    return time.perf_counter() - start


def build(case):
    if "smiles" in case:
        from island.chemistry import from_smiles

        return from_smiles(case["smiles"], random_seed=2026)
    from island.builders import build_linear_polymer

    return build_linear_polymer(
        case["psmiles"],
        dp=case["dp"],
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )


def reference_inputs(system, labels, charges, xyz, directory):
    ids = sorted(labels)
    index = {sid: i + 1 for i, sid in enumerate(ids)}
    names = {sid: system.topology.sites[sid].element + str(index[sid]) for sid in ids}
    neighbors = {sid: [] for sid in ids}
    for b in system.topology.bonds.values():
        neighbors[b.site1].append(b.site2)
        neighbors[b.site2].append(b.site1)
    car = [
        "!BIOSYM archive 3",
        "PBC=OFF",
        "ISLAND independent H5 reference",
        "!DATE fixed 2026-10-04",
    ]
    mdf = ["!BIOSYM molecular_data 4", "#topology"]
    mdf += [
        f"@column {i} {name}"
        for i, name in enumerate(
            [
                "element",
                "atom_type",
                "charge_group",
                "isotope",
                "formal_charge",
                "charge",
                "switching_atom",
                "oop_flag",
                "chirality_flag",
                "occupancy",
                "xray_temp_factor",
                "connections",
            ],
            1,
        )
    ]
    mdf += ["@molecule chain", ""]
    for i, sid in enumerate(ids):
        el = system.topology.sites[sid].element
        car.append(
            f"{names[sid]} "
            + " ".join(format(v, ".17g") for v in xyz[i])
            + f" XXXX 1 {labels[sid]} {el} {charges[sid]:.17g}"
        )
        mdf.append(
            f"XXXX_1:{names[sid]} {el} {labels[sid]} ? 0 0 {charges[sid]:.17g} 0 0 8 1.0 0.0 "
            + " ".join(names[n] for n in sorted(neighbors[sid]))
        )
    car += ["end", "end"]
    mdf += ["", "!", "#end"]
    (directory / "reference.car").write_text("\n".join(car) + "\n")
    (directory / "reference.mdf").write_text("\n".join(mdf) + "\n")
    return ids, neighbors


def data_sections(path):
    sections, current = {}, None
    for line in path.read_text().splitlines():
        s = line.split("#", 1)[0].strip()
        if not s:
            continue
        if s[0].isalpha():
            current = s
            sections[current] = []
        elif current:
            sections[current].append(s.split())
    return sections


def check_inventory(sections, ids, neighbors, charges):
    def canon(row):
        return min(tuple(row), tuple(reversed(row)))

    expected_bonds = {canon((a, b)) for a in ids for b in neighbors[a]}
    expected_angles = {
        canon((a, b, c)) for b in ids for a, c in combinations(neighbors[b], 2)
    }
    expected_propers = {
        canon((a, b, c, d))
        for b, c in expected_bonds
        for a in neighbors[b]
        if a != c
        for d in neighbors[c]
        if d != b and d != a
    }
    expected_impropers = {
        (b, tuple(sorted(arms))) for b in ids for arms in combinations(neighbors[b], 3)
    }
    actual = {}
    for section, expected in (
        ("Bonds", expected_bonds),
        ("Angles", expected_angles),
        ("Dihedrals", expected_propers),
        ("Impropers", expected_impropers),
    ):
        mapped = [[ids[int(v) - 1] for v in row[2:]] for row in sections[section]]
        found = {
            (r[1], tuple(sorted((r[0], r[2], r[3]))))
            if section == "Impropers"
            else canon(r)
            for r in mapped
        }
        if found != expected or len(mapped) != len(expected):
            raise ValueError(f"Independent {section} inventory mismatch")
        actual[section] = len(found)
    for row in sections["Atoms"]:
        if abs(float(row[3]) - charges[ids[int(row[0]) - 1]]) > 1e-12:
            raise ValueError("Converter charge mismatch")
    return actual


def check_model_graph(model, labels, charges, neighbors):
    """Independent explicit walks (not the H4 BFS/pair-generation helper)."""
    distances = {}
    for start in neighbors:
        paths = [(start,)]
        for length in (1, 2, 3):
            paths = [p + (n,) for p in paths for n in neighbors[p[-1]] if n not in p]
            for path in paths:
                pair = tuple(sorted((path[0], path[-1])))
                distances[pair] = min(length, distances.get(pair, 4))
    pairs = []
    for pair, length in sorted(distances.items()):
        pairs.append(
            {
                "sites": list(pair),
                "path_length": length,
                "lj_weight": model["special_pairs"]["lj"][length - 1],
                "coulomb_weight": model["special_pairs"]["coulomb"][length - 1],
            }
        )
    if pairs != model["special_pair_inventory"]:
        raise ValueError("Independent special-pair inventory mismatch")
    actual = {r["site"]: r["charge"] for r in model["nonbonded"]}
    if set(actual) != set(labels) or any(
        abs(actual[s] - charges[s]) > 1e-12 for s in labels
    ):
        raise ValueError("Independent source-native charge mismatch")
    zeros = [r for r in model["terms"] if r["family"] == "bond-bond_1_3"]
    if any(
        r["coefficients"] != [0]
        or r["origin"] != "policy_derived_zero"
        or r["source_rows"]
        or any(labels[s] == "cp" for s in r["sites"])
        for r in zeros
    ):
        raise ValueError("Unexpected policy-zero applicability/provenance")
    return {
        "pairs_checked": len(pairs),
        "sites_checked": len(actual),
        "policy_zeros_checked": len(zeros),
        "charge_atol_e": 1e-12,
    }


def aa_reference_coefficients(sections):
    """Independent dependencies from converter graph/angle rows, not H3 terms.

    Pinned GetParameters.c writes ABC,CBD,ABD while improper_class2.cpp
    consumes ABC,ABD,CBD. Preserve raw file; explicitly override only these
    equilibrium roles in the validation input. Coefficients remain unchanged.
    """
    angles = {
        min(tuple(r[2:]), tuple(reversed(r[2:]))): r[1] for r in sections["Angles"]
    }
    angle_coeff = {r[0]: float(r[1]) for r in sections["Angle Coeffs"]}
    aa = {r[0]: list(map(float, r[1:])) for r in sections["AngleAngle Coeffs"]}
    resolved = {}
    for row in sections["Impropers"]:
        a, b, c, d = row[2:]
        keys = [(a, b, c), (a, b, d), (c, b, d)]
        theta = [angle_coeff[angles[min(k, k[::-1])]] for k in keys]
        old = aa[row[1]]
        if old[3:] != [theta[0], theta[2], theta[1]]:
            raise ValueError("Unexpected converter AA dependency behavior")
        new = old[:3] + theta
        if row[1] in resolved and resolved[row[1]] != new:
            raise ValueError("Conflicting AA type dependencies")
        resolved[row[1]] = new
    return [
        "improper_coeff " + k + " aa " + " ".join(format(v, ".17g") for v in values)
        for k, values in resolved.items()
    ]


def lammps_input(xyz, lj, coul, aa_commands=()):
    span = float(np.max(np.ptp(xyz, axis=0)))
    cutoff = 2 * span + 20
    bound = float(np.max(np.abs(xyz))) + cutoff + 5
    lines = [
        "units real",
        "atom_style full",
        "boundary f f f",
        "bond_style class2",
        "angle_style class2",
        "dihedral_style class2",
        "improper_style class2",
        f"pair_style lj/class2/coul/cut {cutoff:.17g}",
        "special_bonds lj "
        + " ".join(map(str, lj))
        + " coul "
        + " ".join(map(str, coul)),
        "read_data reference.data",
        f"change_box all x final {-bound} {bound} y final {-bound} {bound} z final {-bound} {bound} units box",
        "pair_modify mix sixthpower shift no tail no",
    ]
    lines.extend(aa_commands)
    # Restore declared coordinates at full precision after converter formatting/recentering.
    lines += [
        f"set atom {i + 1} x {r[0]:.17g} y {r[1]:.17g} z {r[2]:.17g}"
        for i, r in enumerate(xyz)
    ]
    lines += [
        "neighbor 2.0 nsq",
        "thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul",
        "run 0",
        'print "ENERGIES $(pe:%.17g) $(ebond:%.17g) $(eangle:%.17g) $(edihed:%.17g) $(eimp:%.17g) $(evdwl:%.17g) $(ecoul:%.17g)"',
        "write_dump all custom forces.dump id x y z fx fy fz modify sort id format float %.17g",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--converter", required=True)
    parser.add_argument("--lammps", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    source_path = Path(args.source).resolve()
    converter = Path(args.converter).resolve()
    lammps = Path(args.lammps).resolve()
    declaration = {
        "schema": "island_pcff_whole_system_acceptance_v1",
        "source": PIN,
        "source_sha256": digest(source_path),
        "converter_sha256": digest(converter),
        "lammps_sha256": digest(lammps),
        "cases": CASES,
        "seed": 2026,
        "polymer_coordinates": "local_template; template/assembly seed=2026",
        "perturbation": "0.015*sin(arange(3*N).reshape(N,3)+0.3) angstrom; small cases only",
        "tolerances": TOL,
        "policies": {
            "all_cases": {"lj": [0, 0, 1], "coulomb": [0, 0, 1]},
            "ethanol_additional": {"lj": [0, 0.2, 0.5], "coulomb": [0, 0.3, 0.7]},
        },
        "reference_AA_policy": "explicit ABC,ABD,CBD equilibrium dependencies from independent converter angle inventory; raw converter ABC,CBD,ABD retained and compared separately",
        "timeouts_seconds": 180,
        "retries": 0,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(out / "declaration.json", json_bytes(declaration))
    source = load_pcff_source(source_path)
    results = []
    for case in CASES:
        folder = out / case["name"]
        folder.mkdir()
        row = {"case": case["name"], "passed": False}
        try:
            start = time.perf_counter()
            system = build(case)
            ids = sorted(system.topology.sites)
            xyz = np.array([system.coordinates.get(s) for s in ids])
            labels = expected_labels(system)
            charges = reference_charges(source.raw, system, labels)
            typing = type_pcff_atoms(system, source)
            if typing.assignments != labels:
                raise ValueError("Independent atom type mismatch")
            native = assign_automatic_pcff_charges(system, typing)
            assignment = assign_pcff_parameters(system, typing, native)
            save_pcff_parameters(assignment, folder / "parameters.json")
            row["source_selection_check"] = oracle_check(
                unpack(assignment.json_text), source.raw
            )
            row["construction_assignment_seconds"] = time.perf_counter() - start
            reference_inputs(system, labels, charges, xyz, folder)
            row["converter_seconds"] = run(
                [
                    converter,
                    "reference",
                    "-class",
                    "II",
                    "-frc",
                    source_path,
                    "-p",
                    "3",
                    "-nocenter",
                ],
                folder,
                "converter",
            )
            sections = data_sections(folder / "reference.data")
            row["independent_inventory"] = check_inventory(
                sections,
                ids,
                {
                    s: [
                        b.site2 if b.site1 == s else b.site1
                        for b in system.topology.bonds.values()
                        if s in (b.site1, b.site2)
                    ]
                    for s in ids
                },
                charges,
            )
            aa_commands = aa_reference_coefficients(sections)
            publish(folder / "aa_reference_commands.json", json_bytes(aa_commands))
            row["comparisons"] = []
            policies = [([0, 0, 1], [0, 0, 1])]
            if case["name"] == "ethanol":
                policies.append(([0, 0.2, 0.5], [0, 0.3, 0.7]))
            for pi, (lj, coul) in enumerate(policies):
                spec = define_pcff_model(
                    assignment, special_pairs=special_pair_policy(lj=lj, coulomb=coul)
                )
                row.setdefault("independent_model_checks", []).append(
                    check_model_graph(
                        unpack(spec.json_text),
                        labels,
                        charges,
                        {
                            sid: [
                                b.site2 if b.site1 == sid else b.site1
                                for b in system.topology.bonds.values()
                                if sid in (b.site1, b.site2)
                            ]
                            for sid in ids
                        },
                    )
                )
                save_pcff_model(spec, folder / f"model-{pi}.json")
                start = time.perf_counter()
                evaluator = PCFFSinglePointEvaluator(system, spec)
                compile_seconds = time.perf_counter() - start
                frames = (
                    [xyz]
                    if case["name"] == "PE_DP50"
                    else [
                        xyz,
                        xyz
                        + 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3),
                    ]
                )
                for fi, frame in enumerate(frames):
                    tag = f"p{pi}-f{fi}"
                    work = folder / tag
                    work.mkdir()
                    (work / "reference.data").write_bytes(
                        (folder / "reference.data").read_bytes()
                    )
                    publish(
                        work / "coordinates.json",
                        json_bytes({str(s): v.tolist() for s, v in zip(ids, frame)}),
                    )
                    (work / "input.lmp").write_text(
                        lammps_input(frame, lj, coul, aa_commands)
                    )
                    duration = run([lammps, "-in", "input.lmp"], work, "lammps")
                    text = (work / "lammps.log").read_text()
                    line = [x for x in text.splitlines() if x.startswith("ENERGIES ")][
                        -1
                    ]
                    energies = np.array([float(v) for v in line.split()[1:]]) * 4.184
                    dump = (
                        (work / "forces.dump")
                        .read_text()
                        .split("ITEM: ATOMS id x y z fx fy fz\n")[1]
                    )
                    forces = (
                        np.array(
                            [
                                [float(v) for v in ln.split()[4:7]]
                                for ln in dump.strip().splitlines()
                            ]
                        )
                        * 4.184
                    )
                    start = time.perf_counter()
                    result = evaluator.evaluate(dict(zip(ids, frame)))
                    elapsed = time.perf_counter() - start
                    calculated = np.array(
                        [result.potential_energy]
                        + [
                            sum(result.energy_components[k] for k in names)
                            for names in GROUPS.values()
                        ]
                    )
                    actual = np.array([result.forces[s] for s in ids])
                    comparison = {
                        "frame": tag,
                        "energy_max_abs_kj_mol": float(max(abs(calculated - energies))),
                        "force_max_abs_kj_mol_angstrom": float(
                            np.max(abs(actual - forces))
                        ),
                        "island_energies": calculated.tolist(),
                        "lammps_energies": energies.tolist(),
                        "compile_seconds": compile_seconds,
                        "evaluation_seconds": elapsed,
                        "lammps_seconds": duration,
                        "coordinate_fingerprint": result.coordinate_fingerprint,
                        "parameter_fingerprint": result.parameter_fingerprint,
                        "passed": bool(
                            np.allclose(calculated, energies, atol=1e-5, rtol=2e-10)
                            and np.allclose(actual, forces, atol=1e-5, rtol=2e-10)
                        ),
                    }
                    row["comparisons"].append(comparison)
            row["passed"] = all(r["passed"] for r in row["comparisons"])
        except Exception as exc:  # noqa: BLE001 -- preserve failed acceptance outcomes
            row["error"] = f"{type(exc).__name__}: {exc}"
        row["artifacts"] = {
            str(p.relative_to(folder)): digest(p)
            for p in sorted(folder.rglob("*"))
            if p.is_file()
        }
        results.append(row)
        publish(folder / "outcome.json", json_bytes(row))
        print(case["name"], row["passed"], row.get("error", ""), flush=True)
    report = {
        "declaration_sha256": digest(out / "declaration.json"),
        "results": results,
        "passed": all(r["passed"] for r in results),
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(out / "report.json", json_bytes(report))
    return int(not report["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
