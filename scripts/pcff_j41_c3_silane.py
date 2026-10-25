"""Finite silane raw-row energy oracle, independent of native PCFF assembly.

This bounded nonperiodic control has no active nonbonded or torsional pairs.
It does not qualify alternate-source chemical typing or intermolecular forces.
"""

import subprocess
from itertools import combinations
from math import acos, pi

import numpy as np

from scripts.pcff_j5_reference import match, read_source, resolve_raw
from scripts.pcff_j41_reference import require
from scripts.phase_4j41_c2 import write


def raw_coefficients(variant):
    rows = read_source(variant.raw)
    si, h, ns = (
        ("si", "hsi", "cff91")
        if variant.filename == "pcff.frc"
        else ("si4", "h1", "compass")
    )
    values, evidence = {}, {}
    for family, labels in (
        ("quartic_bond", [si, h]),
        ("quartic_angle", [h, si, h]),
        ("bond-bond", [h, si, h]),
        ("bond-angle", [h, si, h]),
        ("angle-angle", [h, si, h, h]),
    ):
        q = (
            resolve_raw(rows, family, labels, guarded_msi=True)
            if ns == "cff91"
            else match(rows, family, ns, labels, guarded_msi=True)
        )
        require(q["status"] == "assigned", "Incomplete silane source: " + family)
        indexed = {r["id"]: r for r in rows[family, ns]}
        originals = [indexed[i]["values"] for i in q["selected_ids"]]
        require(all(v == originals[0] for v in originals), "Conflicting source values")
        values[family] = list(map(float, originals[0]))
        evidence[family] = q["selected_ids"]
    if len(values["bond-angle"]) == 1:
        values["bond-angle"] *= 2
    return values, evidence


def raw_energy(xyz, center, coefficients):
    arms = [i for i in range(len(xyz)) if i != center]
    r0, *ks = coefficients["quartic_bond"]
    theta0, *ka = coefficients["quartic_angle"]
    theta0 *= pi / 180
    dr = {i: np.linalg.norm(xyz[i] - xyz[center]) - r0 for i in arms}
    delta = {}
    for i, j in combinations(arms, 2):
        a, b = xyz[i] - xyz[center], xyz[j] - xyz[center]
        delta[i, j] = (
            acos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))
            - theta0
        )
    bond = sum(sum(k * d**p for p, k in enumerate(ks, 2)) for d in dr.values())
    angle = sum(sum(k * d**p for p, k in enumerate(ka, 2)) for d in delta.values())
    bb = coefficients["bond-bond"][0] * sum(dr[i] * dr[j] for i, j in delta)
    left, right = coefficients["bond-angle"]
    ba = sum((left * dr[i] + right * dr[j]) * d for (i, j), d in delta.items())
    aa = 0.0
    for triple in combinations(arms, 3):
        angles = [delta[tuple(sorted(pair))] for pair in combinations(triple, 2)]
        aa += coefficients["angle-angle"][0] * sum(
            a * b for a, b in combinations(angles, 2)
        )
    # LAMMPS real units, source kcal/mol; output kJ/mol.
    return (
        np.array(
            [bond + angle + bb + ba + aa, bond, angle + bb + ba, 0.0, aa, 0.0, 0.0]
        )
        * 4.184
    )


def silane_numerics(system, variant, folder, lammps, *, execute):
    from scripts.pcff_j5_reference import charge_raw
    from scripts.pcff_j41_c3_reference import compass_trial
    from scripts.validate_pcff_singlepoint import (
        check_inventory,
        data_sections,
        lammps_input,
    )

    ids = sorted(system.topology.sites)
    require(
        len(ids) == 5
        and sorted(a.element for a in system.topology.sites.values())
        == ["H", "H", "H", "H", "Si"],
        "Finite silane control only",
    )
    center = next(
        i for i, sid in enumerate(ids) if system.topology.sites[sid].element == "Si"
    )
    xyz = np.array([system.coordinates.get(i) for i in ids])
    coeff, evidence = raw_coefficients(variant)
    sections = data_sections(folder / "reference.data")
    # Silane has no graph dihedrals; the converter omits that empty data section.
    sections.setdefault("Dihedrals", [])
    si, h = ("si", "hsi") if variant.filename == "pcff.frc" else ("si4", "h1")
    labels = {sid: si if i == center else h for i, sid in enumerate(ids)}
    if variant.filename == "pcff.frc":
        q = charge_raw(read_source(variant.raw), system, labels)
        require(q["complete"], "Incomplete independent silane charges")
        charges = {i: float(value) for i, value in q["partial_charges"].items()}
    else:
        q = compass_trial(system, variant)["charges"]
        require(q["complete"], "Incomplete independent silane charges")
        charges = {i: float(value) for i, value in q["vector_e"].items()}
    neighbors = {sid: [] for sid in ids}
    for b in system.topology.bonds.values():
        neighbors[b.site1].append(b.site2)
        neighbors[b.site2].append(b.site1)
    counts = check_inventory(sections, ids, neighbors, charges)
    atom_decls = {r["types"][0]: r for r in variant.payload["atom_types"]}
    type_order = list(dict.fromkeys(labels.values()))
    require(
        [int(row[0]) for row in sections["Atoms"]] == list(range(1, 6)),
        "Converter atom ordering differs",
    )
    require(
        [type_order[int(row[2]) - 1] for row in sections["Atoms"]]
        == list(labels.values()),
        "Converter atom label mapping differs",
    )
    source_rows = read_source(variant.raw)
    ns = "cff91" if variant.filename == "pcff.frc" else "compass"
    masses = {int(row[0]): float(row[1]) for row in sections["Masses"]}
    pairs = {int(row[0]): list(map(float, row[1:])) for row in sections["Pair Coeffs"]}
    for i, label in enumerate(type_order, 1):
        require(
            abs(masses[i] - float(atom_decls[label]["mass"])) <= 1e-12,
            "Converter mass differs from source",
        )
        pair = match(source_rows, "nonbond(9-6)", ns, [label], guarded_msi=True)
        require(pair["status"] == "assigned", "Missing source nonbonded row")
        selected = next(
            r
            for r in source_rows["nonbond(9-6)", ns]
            if r["id"] == pair["selected_ids"][0]
        )
        np.testing.assert_allclose(
            pairs[i], list(map(float, selected["values"][::-1])), atol=1e-12, rtol=1e-13
        )
    r0 = coeff["quartic_bond"][0]
    theta0 = coeff["quartic_angle"][0]
    require(
        len(sections["Atoms"]) == 5
        and len(sections["Bonds"]) == 4
        and len(sections["Angles"]) == 6
        and len(sections["Impropers"]) == 4,
        "Converter interaction count differs",
    )
    expected = {
        "Bond Coeffs": coeff["quartic_bond"],
        "Angle Coeffs": coeff["quartic_angle"],
        "BondBond Coeffs": coeff["bond-bond"] + [r0, r0],
        "BondAngle Coeffs": coeff["bond-angle"] + [r0, r0],
        "AngleAngle Coeffs": coeff["angle-angle"] * 3 + [theta0] * 3,
    }
    for name, values in expected.items():
        require(len(sections[name]) == 1, "Unexpected silane coefficient type count")
        np.testing.assert_allclose(
            list(map(float, sections[name][0][1:])), values, atol=1e-12, rtol=1e-13
        )
    if execute:
        (folder / "lammps.in").write_text(lammps_input(xyz, [0, 0, 1], [0, 0, 1]))
        cmd = [str(lammps), "-in", "lammps.in", "-log", "none"]
        p = subprocess.run(
            cmd, cwd=folder, capture_output=True, text=True, timeout=180, check=False
        )
        (folder / "lammps.log").write_text(p.stdout + p.stderr)
        write(folder / "lammps-command.json", {"argv": cmd, "exit_code": p.returncode})
        require(p.returncode == 0, "Finite LAMMPS control failed")
    text = (folder / "lammps.log").read_text()
    energies = (
        np.array(
            list(
                map(
                    float,
                    next(
                        line
                        for line in text.splitlines()
                        if line.startswith("ENERGIES ")
                    ).split()[1:],
                )
            )
        )
        * 4.184
    )
    force_text = (
        (folder / "forces.dump").read_text().split("ITEM: ATOMS id x y z fx fy fz\n")[1]
    )
    lines = [line.split() for line in force_text.strip().splitlines()]
    require(
        [int(row[0]) for row in lines] == list(range(1, 6)), "Changed force atom order"
    )
    np.testing.assert_allclose(
        np.array([list(map(float, row[1:4])) for row in lines]),
        xyz,
        atol=1e-13,
        rtol=1e-13,
    )
    forces = np.array([list(map(float, row[4:7])) for row in lines]) * 4.184
    raw = raw_energy(xyz, center, coeff)
    np.testing.assert_allclose(raw, energies, atol=1e-5, rtol=2e-10)
    sweep = []
    for h in (0.001, 0.0001, 1e-5):
        numerical = np.zeros_like(xyz)
        for i in range(5):
            for axis in range(3):
                lo, hi = xyz.copy(), xyz.copy()
                lo[i, axis] -= h
                hi[i, axis] += h
                numerical[i, axis] = -(
                    raw_energy(hi, center, coeff)[0] - raw_energy(lo, center, coeff)[0]
                ) / (2 * h)
        sweep.append(
            {
                "h_angstrom": h,
                "max_force_error_kj_mol_angstrom": float(
                    np.max(abs(numerical - forces))
                ),
            }
        )
    require(
        sweep[-1]["max_force_error_kj_mol_angstrom"] < 1e-5,
        "Independent silane finite differences failed",
    )
    return {
        "interaction_counts": counts,
        "atom_order": ids,
        "source_labels": labels,
        "charge_vector_e": charges,
        "source_masses_dalton": {
            label: atom_decls[label]["mass"] for label in type_order
        },
        "source_rows": evidence,
        "source_coefficients": coeff,
        "component_order": [
            "total",
            "bond",
            "angle+bb+ba",
            "dihedral",
            "improper_aa",
            "lj",
            "coulomb",
        ],
        "raw_energies_kj_mol": raw.tolist(),
        "lammps_energies_kj_mol": energies.tolist(),
        "lammps_forces_kj_mol_angstrom": forces.tolist(),
        "finite_difference_sweep": sweep,
        "atol": 1e-5,
        "rtol": 2e-10,
        "scope": "finite isolated silane diagnostic only; no alternate-source execution profile, active nonbonded pairs, periodic or scientific qualification",
    }
