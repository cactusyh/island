"""Deterministic LAMMPS coefficients, retaining source-defined interaction roles."""

from math import degrees

from island.periodic import _require


def _number(value):
    return str(value) if type(value) is int else format(value, ".17g")


def _line(values):
    return " ".join(_number(x) for x in values)


def _pcff(data):
    terms = data["pcff_terms"]
    allowed = {
        "quartic_bond",
        "quartic_angle",
        "bond-bond",
        "bond-angle",
        "torsion_3",
        "middle_bond-torsion_3",
        "end_bond-torsion_3",
        "angle-torsion_3",
        "angle-angle-torsion_1",
        "bond-bond_1_3",
        "angle-angle",
    }
    _require(
        all(t["family"] in allowed for t in terms),
        "LAMMPS PCFF export does not support this bonded interaction family",
    )
    lookup = {(t["family"], tuple(t["sites"])): t for t in terms}
    _require(len(lookup) == len(terms), "Duplicate PCFF numerical terms")
    used = set()

    def get(family, sites):
        key = family, tuple(sites)
        _require(key in lookup, f"Missing PCFF cross term: {family} {sites}")
        used.add(key)
        return lookup[key]

    sections = {}
    inventory = {}

    def add(section, values):
        sections.setdefault(section, []).append(values)

    for kind, family in [
        ("Bonds", "quartic_bond"),
        ("Angles", "quartic_angle"),
        ("Dihedrals", "torsion_3"),
    ]:
        rows = sorted(
            [t for t in terms if t["family"] == family], key=lambda t: t["sites"]
        )
        inventory[kind] = [r["sites"] for r in rows]
        for row in rows:
            sites = row["sites"]
            c = get(family, sites)["coefficients"]
            if kind == "Bonds":
                add("Bond Coeffs", [c[0]] + [v / 4.184 for v in c[1:]])
            elif kind == "Angles":
                add("Angle Coeffs", [degrees(c[0])] + [v / 4.184 for v in c[1:]])
                t = get("bond-bond", sites)
                add(
                    "BondBond Coeffs",
                    [v / 4.184 for v in t["coefficients"]] + t["equilibria"],
                )
                t = get("bond-angle", sites)
                add(
                    "BondAngle Coeffs",
                    [v / 4.184 for v in t["coefficients"]] + t["equilibria"][:2],
                )
            else:
                add(
                    "Dihedral Coeffs",
                    [
                        c[0] / 4.184,
                        degrees(c[1]),
                        c[2] / 4.184,
                        degrees(c[3]),
                        c[4] / 4.184,
                        degrees(c[5]),
                    ],
                )
                for name, f, angular in [
                    ("MiddleBondTorsion Coeffs", "middle_bond-torsion_3", False),
                    ("EndBondTorsion Coeffs", "end_bond-torsion_3", False),
                    ("AngleTorsion Coeffs", "angle-torsion_3", True),
                    ("AngleAngleTorsion Coeffs", "angle-angle-torsion_1", True),
                    ("BondBond13 Coeffs", "bond-bond_1_3", False),
                ]:
                    t = get(f, sites)
                    add(
                        name,
                        [v / 4.184 for v in t["coefficients"]]
                        + [degrees(v) if angular else v for v in t["equilibria"]],
                    )
    # Class-II AA order: K1 ABC*CBD, K2 ABC*ABD, K3 ABD*CBD.
    # The zero out-of-plane stiffness is structural: these tetrahedral AA
    # interactions contain no Wilson term in the validated native model.
    groups = {}
    for t in terms:
        if t["family"] == "angle-angle":
            a, b, c, d = t["sites"]
            groups.setdefault((b, tuple(sorted((a, c, d)))), []).append(t)
    inventory["Impropers"] = []
    for (b, arms), rows in sorted(groups.items()):
        a, c, d = arms
        sites = [a, b, c, d]

        def angle_key(x, y):
            return tuple(sorted((x, y)))

        required = [
            (angle_key(a, c), angle_key(c, d)),
            (angle_key(a, c), angle_key(a, d)),
            (angle_key(a, d), angle_key(c, d)),
        ]
        vals = {}
        eq = {}
        for row in rows:
            x, _, y, z = row["sites"]
            angles = (angle_key(x, y), angle_key(y, z))
            key = tuple(sorted(angles))
            _require(key not in vals, "Duplicate angle-angle role")
            vals[key] = row["coefficients"][0] / 4.184
            for ang, value in zip(angles, row["equilibria"], strict=True):
                _require(
                    ang not in eq or eq[ang] == value,
                    "Conflicting AA equilibrium dependency",
                )
                eq[ang] = value
            used.add(("angle-angle", tuple(row["sites"])))
        _require(
            len(vals) == 3 and set(vals) == {tuple(sorted(x)) for x in required},
            "Missing PCFF angle-angle cross term",
        )
        inventory["Impropers"].append(sites)
        add("Improper Coeffs", [0.0, 0.0])
        add(
            "AngleAngle Coeffs",
            [vals[tuple(sorted(x))] for x in required]
            + [
                degrees(eq[x])
                for x in (angle_key(a, c), angle_key(a, d), angle_key(c, d))
            ],
        )
    _require(used == set(lookup), "Unrepresented PCFF interaction terms")
    return inventory, sections


def render(backend):
    p = backend.payload
    data = p["numerical_data"]
    cfg = p["config"]
    system = backend.system
    pcff = data["family"] == "PCFF"
    sites = {s["id"]: s for s in data["sites"]}
    ids = sorted(sites)
    # Different masses or numerical LJ values get distinct type numbers even
    # when the native label is shared. Stable atom IDs are never renumbered.
    keys = {
        i: (
            sites[i]["atom_type"],
            sites[i]["mass"],
            sites[i]["sigma"],
            sites[i]["epsilon"],
        )
        for i in ids
    }
    unique = sorted(set(keys.values()))
    types = {key: i + 1 for i, key in enumerate(unique)}
    if pcff:
        inventory, sections = _pcff(data)
        styles = {
            "bond": "class2",
            "angle": "class2",
            "dihedral": "class2",
            "improper": "class2",
        }
    else:
        inventory = {
            name: [r["sites"] for r in data[key]]
            for name, key in [
                ("Bonds", "bonds"),
                ("Angles", "angles"),
                ("Dihedrals", "torsions"),
                ("Impropers", "impropers"),
            ]
        }
        sections = {
            "Bond Coeffs": [[r["k"] / 8.368, r["length"]] for r in data["bonds"]],
            "Angle Coeffs": [
                [r["k"] / 8.368, degrees(r["theta"])] for r in data["angles"]
            ],
        }
        styles = {
            "bond": "harmonic",
            "angle": "harmonic",
            "dihedral": "multi/harmonic" if data["family"] == "OPLS-AA" else "fourier",
            "improper": "cvff",
        }
        if data["family"] == "OPLS-AA":
            _require(
                all(r["rb"][5] == 0 for r in data["torsions"]),
                "LAMMPS multi/harmonic cannot represent nonzero RB c5",
            )
            sections["Dihedral Coeffs"] = [
                [v * (-1) ** n / 4.184 for n, v in enumerate(r["rb"][:5])]
                for r in data["torsions"]
            ]
        else:
            sections["Dihedral Coeffs"] = [
                [len(r["fourier"])]
                + [
                    v
                    for k, n, phase in r["fourier"]
                    for v in (k / 4.184, n, degrees(phase))
                ]
                for r in data["torsions"]
            ]
            out = []
            for row in data["impropers"]:
                _require(
                    len(row["fourier"]) == 1,
                    "LAMMPS Amber improper requires one representable periodic term",
                )
                k, n, phase = row["fourier"][0]
                _require(
                    abs(degrees(phase)) < 1e-10 or abs(degrees(phase) - 180) < 1e-10,
                    "Unsupported Amber improper phase",
                )
                out.append([k / 4.184, 1 if abs(phase) < 1e-10 else -1, n])
            sections["Improper Coeffs"] = out
    # Prevent a cutoff-dependent disagreement: scaled pairs are full exceptions
    # in OpenMM, but special pairs use the LJ cutoff in LAMMPS.
    lengths = cfg["box_lengths"]
    for pair in data["pairs"]:
        if pair["lj"] or pair["coulomb"]:
            a, b = pair["sites"]
            delta = system.coordinates.get(a) - system.coordinates.get(b)
            r2 = sum(
                (float(x) - round(float(x) / length) * length) ** 2
                for x, length in zip(delta, lengths, strict=True)
            )
            _require(
                r2 < cfg["nonbonded_cutoff"] ** 2,
                "Scaled special pair exceeds periodic cutoff",
            )
    labels = sorted(set(p["final_graph"]["molecule_membership"].values()))
    molecules = {label: n + 1 for n, label in enumerate(labels)}
    lines = ["ISLAND periodic force-field backend v1", "", f"{len(ids)} atoms"]
    for name, rows in inventory.items():
        lines.append(f"{len(rows)} {name.lower()}")
    lines += ["", f"{len(unique)} atom types"]
    for name, rows in inventory.items():
        lines.append(f"{len(rows)} {name[:-1].lower()} types")
    lines += [""] + [
        f"0 {_number(length)} {axis}lo {axis}hi"
        for axis, length in zip("xyz", lengths, strict=True)
    ]

    def section(name, rows):
        if rows:
            lines.extend(["", name, ""] + [_line(row) for row in rows])

    section("Masses", [[n, key[1]] for key, n in types.items()])
    section("Pair Coeffs", [[n, key[3] / 4.184, key[2]] for key, n in types.items()])
    for name, rows in sections.items():
        section(name, [[n] + row for n, row in enumerate(rows, 1)])
    section(
        "Atoms # full",
        [
            [
                i,
                molecules[p["final_graph"]["molecule_membership"][str(i)]],
                types[keys[i]],
                sites[i]["charge"],
                *map(float, system.coordinates.get(i)),
            ]
            for i in ids
        ],
    )
    for name, rows in inventory.items():
        section(name, [[n, n] + list(row) for n, row in enumerate(rows, 1)])
    mix = (
        "sixthpower"
        if pcff
        else ("geometric" if data["family"] == "OPLS-AA" else "arithmetic")
    )
    lj, coul = cfg["one_four_scaling"]
    script = ["units real", "atom_style full", "boundary p p p"]
    script += [f"{key}_style {style}" for key, style in styles.items()]
    script += [
        f"pair_style {'lj/class2/coul/long' if pcff else 'lj/cut/coul/long'} {_number(cfg['nonbonded_cutoff'])}",
        f"pair_modify mix {mix} shift no tail no table 0",
        f"special_bonds lj 0 0 {_number(lj)} coul 0 0 {_number(coul)} angle no dihedral no",
        "read_data system.data",
        f"kspace_style {'pppm' if cfg['electrostatics_method'] == 'pme' else 'ewald'} {_number(cfg['pme_tolerance'])}",
        "neighbor 2 bin",
        "neigh_modify every 1 delay 0 check yes",
        "thermo 1",
        "thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong",
        "thermo_modify format float %.17g",
        "dump forces all custom 1 forces.dump id type q x y z fx fy fz",
        "dump_modify forces sort id format float %.17g",
        "run 0",
        "",
    ]
    return "\n".join(lines) + "\n", "\n".join(script)
