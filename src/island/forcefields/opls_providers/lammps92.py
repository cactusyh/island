"""Strict LAMMPS text readers. Never execute a command or convert a coefficient."""

import re
from collections import deque

from island.forcefields.assignment_core.schema import identity, require

from .schema import number, sealed

WIDTHS = {
    "Masses": 2,
    "Pair Coeffs": 3,
    "Bond Coeffs": 3,
    "Angle Coeffs": 3,
    "Dihedral Coeffs": 5,
    "Improper Coeffs": 4,
    "Atoms": 7,
    "Bonds": 4,
    "Angles": 5,
    "Dihedrals": 6,
    "Impropers": 6,
}
COUNT_KEYS = {
    "Masses": "atom types",
    "Pair Coeffs": "atom types",
    "Atoms": "atoms",
    "Bond Coeffs": "bond types",
    "Bonds": "bonds",
    "Angle Coeffs": "angle types",
    "Angles": "angles",
    "Dihedral Coeffs": "dihedral types",
    "Dihedrals": "dihedrals",
    "Improper Coeffs": "improper types",
    "Impropers": "impropers",
}
STYLES = {
    "bond": {
        "style": "harmonic",
        "equation": "K*(r-r0)^2",
        "order": ["K", "r0"],
        "units": ["kcal/mol/Angstrom^2", "Angstrom"],
    },
    "angle": {
        "style": "harmonic",
        "equation": "K*(theta-theta0)^2",
        "order": ["K", "theta0"],
        "units": ["kcal/mol/radian^2", "degree"],
        "displacement": "radians",
    },
    "dihedral": {
        "style": "opls",
        "equation": "(K1*(1+cos(phi))+K2*(1-cos(2*phi))+K3*(1+cos(3*phi))+K4*(1-cos(4*phi)))/2",
        "order": ["K1", "K2", "K3", "K4"],
        "units": ["kcal/mol"] * 4,
        "periodicities": [1, 2, 3, 4],
        "atom_order": "source_order_preserved",
    },
    "improper": {
        "style": "cvff",
        "equation": "K*(1+d*cos(n*phi))",
        "order": ["K", "d", "n"],
        "units": ["kcal/mol", "sign", "integer"],
        "atom_order": "source_order_preserved_requires_declared_center",
    },
    "nonbonded": {
        "style": "lj/cut/coul/cut",
        "equation": "4*epsilon*((sigma/r)^12-(sigma/r)^6)",
        "order": ["epsilon", "sigma"],
        "units": ["kcal/mol", "Angstrom"],
        "sigma": "zero_crossing",
        "minimum": "2^(1/6)*sigma",
        "epsilon": "positive_well_depth",
    },
}
EXPANSION_VERSION = "island.j92.shortest-bond-distance-1-2-3.v1"


def integer(token, reason="integer_token_required", *, zero=False):
    require(
        isinstance(token, str) and re.fullmatch(r"[0-9]+", token) is not None, reason
    )
    v = int(token)
    require(v >= (0 if zero else 1), reason)
    return v


def source(file, line, raw, section):
    return {
        "path": file["path"],
        "sha256": file["sha256"],
        "line": line,
        "section": section,
        "row_identity": identity(raw),
    }


def parse_data(raw, file):
    lines = raw.decode("utf-8").splitlines()
    require(lines, "empty_LAMMPS_data")
    match = re.fullmatch(
        r"LAMMPS data file generated with LigParGen ([\w.+-]+) - Jorgensen Lab",
        lines[0],
    )
    require(match is not None, "provider_version_header_required")
    counts, rows, boxes, current = {}, {}, {}, None
    for line_no, original in enumerate(lines[1:], 2):
        content, _, annotation = original.partition("#")
        content = content.strip()
        if not content:
            continue
        if content in WIDTHS:
            require(content not in rows, "duplicate_section:" + content)
            require(
                not annotation.strip()
                or annotation.strip()
                == {
                    "Atoms": "full",
                    "Pair Coeffs": "lj/cut/coul/cut",
                    "Bond Coeffs": "harmonic",
                    "Angle Coeffs": "harmonic",
                    "Dihedral Coeffs": "opls",
                    "Improper Coeffs": "cvff",
                }.get(content),
                "section_style_annotation_mismatch",
            )
            current = content
            rows[current] = []
            continue
        fields = content.split()
        if current is None:
            declaration = re.fullmatch(
                r"([0-9]+) ("
                + "|".join(sorted(set(COUNT_KEYS.values()), key=len, reverse=True))
                + r")",
                content,
            )
            if declaration:
                name = declaration[2]
                require(name not in counts, "duplicate_header_count")
                counts[name] = integer(declaration[1], zero=True)
            elif len(fields) == 4 and fields[2:] in (
                ["xlo", "xhi"],
                ["ylo", "yhi"],
                ["zlo", "zhi"],
            ):
                require(fields[2] not in boxes, "duplicate_box_bounds")
                require(number(fields[0]) < number(fields[1]), "invalid_box_bounds")
                boxes[fields[2]] = fields[:2]
            else:
                raise ValueError("unsupported_LAMMPS_header_or_section")
            continue
        require(len(fields) == WIDTHS[current], "invalid_row_width:" + current)
        rid = integer(fields[0])
        require(
            rid not in {r["id"] for r in rows[current]}, "duplicate_row_id:" + current
        )
        for token in fields[1:]:
            number(token)
        if current in ("Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"):
            integer(fields[1])
            integer(fields[2])
        if current in ("Bonds", "Angles", "Dihedrals", "Impropers"):
            for token in fields[2:]:
                integer(token)
            require(
                len(set(fields[2:])) == len(fields[2:]), "repeated_interaction_site"
            )
        if current == "Masses":
            require(number(fields[1]) > 0, "invalid_mass")
        if current == "Pair Coeffs":
            require(
                number(fields[1]) >= 0 and number(fields[2]) > 0,
                "invalid_pair_coefficients",
            )
        if current in ("Bond Coeffs", "Angle Coeffs"):
            require(
                number(fields[1]) >= 0 and number(fields[2]) > 0,
                "invalid_primary_coefficients",
            )
        if current == "Improper Coeffs":
            require(
                number(fields[1]) >= 0 and fields[2] in ("-1", "1"),
                "invalid_cvff_improper",
            )
            require(
                integer(fields[3], zero=True) in (0, 1, 2, 3, 4, 6),
                "unsupported_cvff_periodicity",
            )
        rows[current].append(
            {
                "id": rid,
                "tokens": fields,
                "source": source(file, line_no, original, current),
            }
        )
    require(set(boxes) == {"xlo", "ylo", "zlo"}, "orthogonal_box_required")
    require(
        all(
            k in counts
            for k in set(COUNT_KEYS.values()) - {"impropers", "improper types"}
        ),
        "complete_header_counts_required",
    )
    for section, key in COUNT_KEYS.items():
        require(
            len(rows.get(section, [])) == counts.get(key, 0),
            "section_count_mismatch:" + section,
        )
    require(counts["atoms"] > 0 and rows.get("Atoms"), "atoms_required")
    # Coefficient and atom type references must be valid; never supply defaults.
    tables = {name: {r["id"]: r for r in rr} for name, rr in rows.items()}
    for section, header in COUNT_KEYS.items():
        if section.endswith("Coeffs") or section == "Masses":
            require(
                set(tables.get(section, {}))
                == set(range(1, counts.get(header, 0) + 1)),
                "coefficient_type_inventory_mismatch:" + section,
            )
    atoms = tables["Atoms"]
    for atom in atoms.values():
        typ = int(atom["tokens"][2])
        require(
            typ in tables.get("Masses", {}) and typ in tables.get("Pair Coeffs", {}),
            "unresolved_atom_type",
        )
    for section, coeff in (
        ("Bonds", "Bond Coeffs"),
        ("Angles", "Angle Coeffs"),
        ("Dihedrals", "Dihedral Coeffs"),
        ("Impropers", "Improper Coeffs"),
    ):
        for row in rows.get(section, []):
            require(
                int(row["tokens"][1]) in tables.get(coeff, {}),
                "missing_coefficient_reference:" + section,
            )
            require(
                all(int(t) in atoms for t in row["tokens"][2:]),
                "unknown_interaction_atom",
            )
    return {
        "provider": "LigParGen",
        "provider_version": match[1],
        "header_source": source(file, 1, lines[0], "header"),
        "counts": counts,
        "rows": rows,
    }


def parse_settings(raw, file, data_name):
    """Small declarative grammar. Recognized execution requests are recorded inertly."""
    commands, inert = {}, []
    for line_no, original in enumerate(raw.decode("utf-8").splitlines(), 1):
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        tokens = line.split()
        name = tokens[0]
        loc = source(file, line_no, original, "settings")
        # Recognize this supplier's run-0/energy-print tail, but NEVER evaluate it.
        if (
            tokens == ["run", "0"]
            or tokens == ["variable", "e", "equal", "pe"]
            or (
                name == "print"
                and re.fullmatch(r'print\s+"LIG: \$\{e\} kcal/mol"', line)
            )
        ):
            inert.append({"command": line, "source": loc, "executed": False})
            continue
        require(name not in commands, "duplicate_or_overridden_setting:" + name)
        require(
            not any(c in line for c in ("$", "&", ";", "\\")),
            "dynamic_settings_forbidden",
        )
        require(
            name
            in (
                "units",
                "atom_style",
                "boundary",
                "pair_style",
                "pair_modify",
                "bond_style",
                "angle_style",
                "dihedral_style",
                "improper_style",
                "special_bonds",
                "read_data",
                "neighbor",
            ),
            "unsupported_settings_command:" + name,
        )
        commands[name] = {"tokens": tokens[1:], "source": loc, "raw_command": line}
    expected = {
        "units": ["real"],
        "atom_style": ["full"],
        "boundary": ["m", "m", "m"],
        "pair_style": ["lj/cut/coul/cut", "40.0", "40.0"],
        "pair_modify": ["mix", "geometric", "tail", "no"],
        "bond_style": ["harmonic"],
        "angle_style": ["harmonic"],
        "dihedral_style": ["opls"],
        "improper_style": ["cvff"],
        "read_data": [data_name],
    }
    blockers = []
    for name, value in expected.items():
        if name not in commands:
            blockers.append("required_setting_missing:" + name)
        elif commands[name]["tokens"] != value:
            # Accept equivalent numeric spellings for the explicitly supplied cutoffs.
            v = commands[name]["tokens"]
            if not (
                name == "pair_style"
                and len(v) == 3
                and v[0] == value[0]
                and all(number(x) == 40 for x in v[1:])
            ):
                blockers.append("setting_mismatch:" + name)
    special = commands.get("special_bonds")
    lj, coul = None, None
    if special:
        tokens = special["tokens"]
        if len(tokens) == 4 and tokens[0] == "lj/coul":
            lj = coul = tokens[1:]
        elif len(tokens) == 8 and tokens[0] == "lj" and tokens[4] == "coul":
            lj, coul = tokens[1:4], tokens[5:8]
        else:
            blockers.append("unsupported_special_bonds_expansion")
        if lj is not None and (
            list(map(number, lj)) != [0, 0, number("0.5")]
            or list(map(number, coul)) != [0, 0, number("0.5")]
        ):
            blockers.append("setting_mismatch:special_bonds")
    else:
        blockers.append("required_setting_missing:special_bonds")
    if "neighbor" in commands:
        v = commands["neighbor"]["tokens"]
        require(
            len(v) == 2 and v[1] == "bin" and number(v[0]) >= 0,
            "unsupported_neighbor_setting",
        )
    if "read_data" in commands:
        for name in (
            "units",
            "atom_style",
            "boundary",
            "pair_style",
            "bond_style",
            "angle_style",
            "dihedral_style",
            "improper_style",
        ):
            if (
                name in commands
                and commands[name]["source"]["line"]
                > commands["read_data"]["source"]["line"]
            ):
                blockers.append("setting_after_read_data:" + name)
    return {
        "commands": commands,
        "inert_commands": inert,
        "blockers": sorted(blockers),
        "conventions": {
            "units": "real",
            "energy": "kcal/mol",
            "length": "Angstrom",
            "charge": "e",
            "forms": STYLES,
            "mixing": "geometric",
            "scale14_lj": lj[2] if lj else None,
            "scale14_coulomb": coul[2] if coul else None,
            "exclusion_distances": [1, 2],
            "lj_cutoff_Angstrom": "40",
            "coulomb_cutoff_Angstrom": "40",
            "long_range": "none_truncated_coulomb_explicit_coul_cut_style",
            "tail": "no",
            "validated": True,
        }
        if not blockers
        else {"validated": False},
    }


def expand_pairs(atom_ids, bonds, graph_identity, settings):
    require(not settings["blockers"], "settings_not_validated_no_pair_expansion")
    adjacency = {i: set() for i in atom_ids}
    edges = set()
    for a, b in bonds:
        require(a in adjacency and b in adjacency and a != b, "invalid_expansion_edge")
        key = tuple(sorted((a, b)))
        require(key not in edges, "duplicate_expansion_edge")
        edges.add(key)
        adjacency[a].add(b)
        adjacency[b].add(a)
    require(adjacency, "empty_expansion_graph")
    distances = {1: set(), 2: set(), 3: set()}
    for start in sorted(adjacency):
        visited = {start: 0}
        queue = deque([start])
        while queue:
            a = queue.popleft()
            if visited[a] == 3:
                continue
            for b in sorted(adjacency[a]):
                if b in visited:
                    continue
                visited[b] = visited[a] + 1
                queue.append(b)
                if start < b:
                    distances[visited[b]].add((start, b))
    # BFS over all atoms establishes the full finite graph, not a fragment.
    seen = set()
    queue = [min(adjacency)]
    while queue:
        a = queue.pop()
        if a not in seen:
            seen.add(a)
            queue.extend(adjacency[a] - seen)
    require(seen == set(adjacency), "disconnected_graph_expansion_blocked")
    command = settings["commands"]["special_bonds"]
    body = {
        "algorithm_version": EXPANSION_VERSION,
        "graph_identity": graph_identity,
        "indexed_bond_graph_hash": identity(
            {"atoms": sorted(atom_ids), "bonds": sorted(edges)}
        ),
        "source_command": command,
        "shortest_distance_policy": "each pair occurs once at its shortest bond distance",
        "pair14": [
            {
                "sites": list(p),
                "distance": 3,
                "lj_scale": settings["conventions"]["scale14_lj"],
                "coulomb_scale": settings["conventions"]["scale14_coulomb"],
            }
            for p in sorted(distances[3])
        ],
        "exclusions": [
            {"sites": list(p), "distance": d, "lj_scale": "0", "coulomb_scale": "0"}
            for d in (1, 2)
            for p in sorted(distances[d])
        ],
        "record_origin": "deterministic_expansion_of_explicit_source_policy",
        "numeric_parameter_generation": False,
    }
    return {**body, "expansion_hash": identity(body)}


def row_metadata(parsed):
    return [
        sealed(
            {
                "section": section,
                "id": row["id"],
                "source": row["source"],
                "claim_hash": identity(row),
            }
        )
        for section, rows in parsed["rows"].items()
        for row in rows
    ]
