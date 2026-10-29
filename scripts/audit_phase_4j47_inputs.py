"""Independent standard-library source/graph/charge/coefficient reconstruction.

No production parser, ISLAND import, parameter default or numerical engine.
Explicit ZIP bytes are independently pinned; no literature assets are copied.
"""

import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from zipfile import ZipFile

TYPES = {"H": 98, "O": 122, "Si": 868, "C": 871}
DFT = {
    "O": ["-1.261", "-1.255"],
    "Si": ["2.121", "1.874"],
    "C": ["-1.102", "-1.083"],
    "H": ["0.224", "0.2225"],
}


def require(condition, stage):
    if not condition:
        raise ValueError(stage)


def digest(value):
    return sha256(
        (
            json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
            + "\n"
        ).encode()
    ).hexdigest()


def read(path):
    e = json.loads(Path(path).read_text())
    require(
        set(e) == {"payload", "sha256"} and digest(e["payload"]) == e["sha256"],
        "independent_envelope_hash",
    )

    def decode(x):
        if x is None or type(x) in (str, int, float, bool):
            return x
        require(type(x) is dict and len(x) == 1, "independent_container")
        tag, values = next(iter(x.items()))
        if tag == "map":
            require(
                len({row[0] for row in values}) == len(values),
                "independent_duplicate_key",
            )
            return {k: decode(v) for k, v in values}
        require(tag in ("list", "tuple"), "independent_tag")
        return [decode(v) for v in values]

    return decode(e["payload"])


def sections(text):
    names = {
        "Masses",
        "Atoms",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        "Velocities",
        "Pair Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
    }
    result, active = {}, None
    for n, raw in enumerate(text.splitlines(), 1):
        t = raw.split("#")[0].strip()
        if t in names:
            active = t
            result[active] = {}
        elif active and t:
            cols = t.split()
            result[active][int(cols[0])] = (n, raw, cols[1:])
    return result


def audit(directory, archive, archive_sha256, *, dps=(0, 2, 4, 6)):
    raw = Path(archive).read_bytes()
    require(sha256(raw).hexdigest() == archive_sha256, "independent_archive_hash")
    with ZipFile(archive) as z:

        def text(name):
            matches = [n for n in z.namelist() if n.endswith("/" + name)]
            require(len(matches) == 1, "independent_archive_membership")
            return z.read(matches[0]).decode()

        source = sections(text("param/OPLS/opls_system.data"))
        settings = text("param/OPLS/opls_system.in.settings")
        raw_files = {
            name: sha256(text(name).encode()).hexdigest()
            for name in (
                "README.md",
                "param/OPLS/opls_system.in",
                "param/OPLS/opls_system.in.init",
                "param/OPLS/opls_system.in.settings",
                "param/OPLS/opls_system.data",
            )
        }
    coefficients = {}
    for n, line in enumerate(settings.splitlines(), 1):
        cols = line.split("#")[0].split()
        if not cols:
            continue
        command, tid = cols[0], int(cols[1])
        key = command, tid
        require(key not in coefficients, "independent_duplicate_coefficient")
        if command == "pair_coeff":
            require(cols[1] == cols[2], "independent_pair_override")
        coefficients[key] = (n, line, cols[3:] if command == "pair_coeff" else cols[2:])
    source_charge = {tid: set() for tid in TYPES.values()}
    for _, _, row in source["Atoms"].values():
        source_charge[int(row[1])].add(row[2])
    results = {}
    for dp in dps:
        root = Path(directory) / f"dp{dp}"
        s = read(root / "system.json")
        p = read(root / "source-preparation.json")
        require(p["archive_sha256"] == archive_sha256, "independent_stale_archive")
        require(
            p["source"]["source_hashes"] == raw_files, "independent_companion_hashes"
        )
        require(
            p["identity"] == digest({k: v for k, v in p.items() if k != "identity"}),
            "independent_assignment_hash",
        )
        atoms = {a["id"]: a for a in s["sites"]}
        require(
            len(atoms) == (27 if dp == 0 else 17 + 10 * dp), "independent_fixture_size"
        )
        order = [a["id"] for a in s["sites"]]
        adjacency = {i: set() for i in atoms}
        for b in s["bonds"]:
            require(b["order"] == 1 and not b["aromatic"], "independent_bond_state")
            a, c = b["site1"], b["site2"]
            adjacency[a].add(c)
            adjacency[c].add(a)
        element = {i: a["element"] for i, a in atoms.items()}
        ends = {
            i
            for i in atoms
            if element[i] == "Si"
            and Counter(element[j] for j in adjacency[i]) == {"C": 3, "O": 1}
        }
        require(len(ends) == 2, "independent_termini")
        for i, a in atoms.items():
            require(
                a["formal_charge"] == 0
                and not a["metadata"].get("radical_electrons", 0),
                "independent_chemical_state",
            )
            counts = Counter(element[j] for j in adjacency[i])
            require(
                (element[i] == "H" and counts == {"C": 1})
                or (element[i] == "C" and counts == {"Si": 1, "H": 3})
                or (element[i] == "O" and counts == {"Si": 2})
                or (element[i] == "Si" and (counts == {"O": 2, "C": 2} or i in ends)),
                "independent_atom_environment",
            )
        anchors = p["e_oxygen_site_ids"]
        eo = anchors[0] if anchors else next(i for i in atoms if element[i] == "O")
        require(
            element[eo] == "O" and bool(adjacency[eo] & ends), "independent_E_oxygen"
        )
        e = set(ends) | {eo}
        for i in ends:
            for c in adjacency[i]:
                if element[c] == "C":
                    e.add(c)
                    e.update(adjacency[c] - {i})
        require(len(e) == 27, "independent_E_fragment")
        require([a["site_id"] for a in p["mapping"]] == order, "independent_atom_order")
        charges = []
        for a in p["mapping"]:
            i = a["site_id"]
            el = element[i]
            tid = TYPES[el]
            q = format(float(DFT[el][int(i in e)]), ".16e")
            require(
                a["element"] == el
                and a["fragment"] == ("E" if i in e else "M")
                and a["molecule_id"] == 1
                and a["author_type_id"] == tid
                and a["charge_e"] == q
                and q in source_charge[tid],
                "independent_type_charge_component",
            )
            require(
                Fraction(a["mass_dalton"]) == Fraction(source["Masses"][tid][2][0]),
                "independent_mass",
            )
            charges.append(q)
        require(
            str(sum(map(Fraction, charges), Fraction())) == p["charge_policy"]["sum_e"]
            and p["charge_policy"]["normalization"] is False
            and p["charge_policy"]["redistribution"] is False,
            "independent_charge_origin",
        )
        rows = {r["evidence_checksum"]: r for r in p["source"]["rows"]}
        for r in rows.values():
            command = {
                "bond": "bond_coeff",
                "angle": "angle_coeff",
                "proper_torsion": "dihedral_coeff",
                "nonbonded": "pair_coeff",
            }[r["family"]]
            n, line, values = coefficients[command, r["author_type_id"]]
            require(
                r["units"]
                == "LAMMPS real: kcal/mol, angstrom, degree equilibrium, radian displacement"
                and r["conversion_version"]
                == "lammps_real_to_nm_kJ_equation_preserving_v1",
                "independent_unit_policy",
            )
            require(
                (n, line)
                == (r["settings_evidence"]["line"], r["settings_evidence"]["raw_text"])
                and values == r["values"],
                "independent_raw_coefficient",
            )
            require(
                r["evidence_checksum"]
                == digest({k: v for k, v in r.items() if k != "evidence_checksum"}),
                "independent_row_checksum",
            )
        paths = {
            "bond": {tuple(sorted((b["site1"], b["site2"]))) for b in s["bonds"]},
            "angle": {
                min((a, b, c), (c, b, a))
                for b in atoms
                for a, c in combinations(sorted(adjacency[b]), 2)
            },
            "proper_torsion": {
                min((a, b, c, d), (d, c, b, a))
                for b in atoms
                for c in adjacency[b]
                for a in adjacency[b] - {c}
                for d in adjacency[c] - {b, a}
            },
            "nonbonded": {(i,) for i in atoms},
        }
        for kind, expected in paths.items():
            selected = p["interactions"][kind]
            require(
                len(selected) == len(expected)
                and {tuple(r["sites"]) for r in selected} == expected,
                "independent_interaction_inventory",
            )
            for r in selected:
                row = rows[r["row_checksum"]]
                labels = tuple(element[i] for i in r["sites"])
                require(
                    row["family"] == kind
                    and tuple(row["labels"]) == min(labels, labels[::-1]),
                    "independent_interaction_labels",
                )
        exported = sections((root / "lammps/system.data").read_text())
        for n, a in enumerate(p["mapping"], 1):
            row = exported["Atoms"][n][2]
            i = a["site_id"]
            require(
                int(row[0]) == a["molecule_id"]
                and int(row[1]) == {"H": 1, "O": 2, "Si": 3, "C": 4}[a["element"]]
                and row[2] == a["charge_e"]
                and list(map(float, row[3:])) == s["coordinates"][i],
                "independent_export_atom_mapping",
            )
        results[str(dp)] = {
            "atoms": len(atoms),
            "source_rows": len(rows),
            "raw_source_hashes": raw_files,
            "ordered_charge_checksum": digest(charges),
            "independent_interaction_counts": {k: len(v) for k, v in paths.items()},
            "graph_type_charge_mass_and_source_rows_reconstructed": True,
        }
    return {
        "schema": "island.j47-independent-input-audit.v1",
        "archive_sha256": archive_sha256,
        "raw_settings_coefficient_rows": len(coefficients),
        "cases": results,
        "numerical_MD_executed": False,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    report = audit(a.directory, a.archive, a.archive_sha256)
    a.output.write_text(
        json.dumps({**report, "identity": digest(report)}, indent=2, sort_keys=True)
        + "\n"
    )


if __name__ == "__main__":
    main()
