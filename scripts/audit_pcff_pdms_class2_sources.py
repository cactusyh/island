"""Independent ASCII/ZIP/PDF evidence for J54; never imports a force-field parser.

COMPASS archive rows are kept separate from FRC and Sun literature records.
A literal author-exported zero is evidence of that input, not proof of a missing
silicon coupling's physical absence or of the generating converter's provenance.
"""

import argparse
import json
import zipfile
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from scripts.audit_pdms_source_completion import strict_resolution
from scripts.pcff_j5_reference import ARITY, charge_raw, raw_inventory, read_source

ARCHIVE_PIN = "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4"
COMPASS_ROOT = "PDMS_benchmark-main/param/COMPASS/"
DOI = "10.1021/acs.jpcb.4c08471"


def coefficient_units(command, subterm):
    """Label author real/class2 cells; no values or missing coefficients generated."""
    e, l, a = "kcal/mol", "angstrom", "radian"
    return {
        ("pair_coeff", "main"): [e, l],
        ("bond_coeff", "main"): [l] + [f"{e}/{l}^{n}" for n in (2, 3, 4)],
        ("angle_coeff", "main"): ["degree"] + [f"{e}/{a}^{n}" for n in (2, 3, 4)],
        ("angle_coeff", "bb"): [f"{e}/{l}^2", l, l],
        ("angle_coeff", "ba"): [f"{e}/({l}*{a})"] * 2 + [l, l],
        ("dihedral_coeff", "main"): [e, "degree"] * 3,
        ("dihedral_coeff", "mbt"): [f"{e}/{l}"] * 3 + [l],
        ("dihedral_coeff", "ebt"): [f"{e}/{l}"] * 6 + [l, l],
        ("dihedral_coeff", "at"): [f"{e}/{a}"] * 6 + ["degree"] * 2,
        ("dihedral_coeff", "aat"): [f"{e}/{a}^2", "degree", "degree"],
        ("dihedral_coeff", "bb13"): [f"{e}/{l}^2", l, l],
        ("improper_coeff", "main"): [f"{e}/{a}^2", "degree"],
        ("improper_coeff", "aa"): [f"{e}/{a}^2"] * 3 + ["degree"] * 3,
    }[command, subterm]


def signed(value):
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return {**value, "evidence_checksum": sha256(raw).hexdigest()}


def archive_audit(archive):
    raw = Path(archive).read_bytes()
    if sha256(raw).hexdigest() != ARCHIVE_PIN:
        raise ValueError("author_archive_hash_mutation")
    with zipfile.ZipFile(archive) as z:
        names = [n for n in z.namelist() if not n.endswith("/")]
        if len(names) != len(set(names)) or any(
            Path(n).is_absolute() or ".." in Path(n).parts for n in names
        ):
            raise ValueError("unsafe_or_duplicate_archive_members")
        inventory = [
            {"file": n, "size": len(z.read(n)), "sha256": sha256(z.read(n)).hexdigest()}
            for n in sorted(names)
        ]
        source = {
            name: z.read(COMPASS_ROOT + name)
            for name in (
                "system.in",
                "system.data",
                "system.in.init",
                "system.in.settings",
                "system.in.charges",
            )
        }
    parts = defaultdict(list)
    section = None
    for n, line in enumerate(source["system.data"].decode().splitlines(), 1):
        w = line.split("#", 1)[0].split()
        if not w:
            continue
        if w[0][0].isalpha():
            section = " ".join(w)
        elif section:
            parts[section].append((n, line, w))
    atoms = {
        int(w[0]): {"type": int(w[2]), "charge": w[3], "line": n, "raw": raw}
        for n, raw, w in parts["Atoms"]
    }
    adj = defaultdict(list)
    for _, _, w in parts["Bonds"]:
        i, j = map(int, w[2:4])
        adj[i].append(j)
        adj[j].append(i)
    active = {
        k: {int(w[1]) for _, _, w in parts[k]}
        for k in ("Bonds", "Angles", "Dihedrals", "Impropers")
    }
    active_types = {a["type"] for a in atoms.values()}
    improper_centers = Counter(
        atoms[int(w[3])]["type"] for _, _, w in parts["Impropers"]
    )
    coefficient_rows = []
    command_sections = {
        "bond_coeff": "Bonds",
        "angle_coeff": "Angles",
        "dihedral_coeff": "Dihedrals",
        "improper_coeff": "Impropers",
    }
    subterms = {"bb", "ba", "mbt", "ebt", "at", "aat", "bb13", "aa"}
    patterns = {}
    for command, section in command_sections.items():
        found = defaultdict(set)
        for _, _, w in parts[section]:
            found[int(w[1])].add(tuple(atoms[int(i)]["type"] for i in w[2:]))
        patterns[command] = found
    for n, line in enumerate(source["system.in.settings"].decode().splitlines(), 1):
        words = line.split("#", 1)[0].split()
        if not words:
            continue
        command = words[0]
        used = (
            command == "pair_coeff"
            and int(words[1]) in active_types
            and int(words[2]) in active_types
        ) or (
            command in command_sections
            and int(words[1]) in active[command_sections[command]]
        )
        if not used:
            continue
        start = 3 if command == "pair_coeff" or words[2] in subterms else 2
        values = words[start:]
        subterm = words[2] if words[2] in subterms else "main"
        units = coefficient_units(command, subterm)
        if len(values) != len(units):
            raise ValueError("coefficient_shape_or_functional_form_changed")
        coefficient_rows.append(
            signed(
                {
                    "source_file": COMPASS_ROOT + "system.in.settings",
                    "source_sha256": sha256(source["system.in.settings"]).hexdigest(),
                    "archive_sha256": ARCHIVE_PIN,
                    "doi": DOI,
                    "citation": "Xiang et al., J. Phys. Chem. B, DOI "
                    + DOI
                    + "; separate author COMPASS inputs",
                    "location": {"line": n},
                    "raw_text": line,
                    "command": command,
                    "type": words[1],
                    "subterm": subterm,
                    "labels": [
                        list(p) for p in sorted(patterns[command][int(words[1])])
                    ]
                    if command != "pair_coeff"
                    else [[int(words[1]), int(words[2])]],
                    "original_values": values,
                    "converted_values": values,
                    "original_units": units,
                    "converted_units": units,
                    "unit_authority": "explicit author units real/styles class2; command dimensions are converter/backend interpretation, not new paper facts",
                    "functional_form": "author " + command + "/class2/" + subterm,
                    "conversion_formula": "identity; author engine input retained without physical reassignment",
                    "implementation_version": "independent-j54-ascii-v1",
                    "explicit_zero_cells": [
                        i for i, v in enumerate(values) if float(v) == 0
                    ],
                    "authority": "author input example only; not a new ISLAND parameter assignment",
                }
            )
        )
    settings = [
        signed(
            {
                "source_file": COMPASS_ROOT + "system.in.init",
                "source_sha256": sha256(source["system.in.init"]).hexdigest(),
                "archive_sha256": ARCHIVE_PIN,
                "doi": DOI,
                "line": n,
                "raw_text": line,
            }
        )
        for n, line in enumerate(source["system.in.init"].decode().splitlines(), 1)
        if line.strip() and not line.startswith("#")
    ]
    types = []
    for t in sorted(active_types):
        motifs = Counter(
            tuple(sorted(atoms[j]["type"] for j in adj[i]))
            for i, a in atoms.items()
            if a["type"] == t
        )
        n, raw, _ = next(r for r in parts["Masses"] if int(r[2][0]) == t)
        types.append(
            {
                "type": t,
                "mass_source_line": n,
                "raw_mass": raw,
                "observed_neighbor_types": [
                    {"neighbors": list(k), "count": v}
                    for k, v in sorted(motifs.items())
                ],
                "source_charges": sorted(
                    {a["charge"] for a in atoms.values() if a["type"] == t}
                ),
                "example": next(
                    {"site": i, **a} for i, a in atoms.items() if a["type"] == t
                ),
            }
        )
    # Type 44 has explicit Si mass and four-connected methyl/bridging-O motifs.
    # The source topology has no improper/angle-angle centered on these Si sites.
    silicon = [i for i, a in atoms.items() if a["type"] == 44]
    absent = []
    for i in silicon:
        if improper_centers[44]:
            raise ValueError("archive_silicon_improper_scope_changed")
        for arms in combinations(sorted(adj[i]), 3):
            absent.append({"center": i, "arms": list(arms)})
    missing_patterns = Counter(
        tuple(sorted(atoms[j]["type"] for j in r["arms"])) for r in absent
    )
    # A graph-specific bulk charge list is not a documented charge method for
    # newly constructed HMDS/DP2/DP4/DP6 graphs.
    charge_commands = source["system.in.charges"].decode().splitlines()
    includes = [
        l
        for l in source["system.in"].decode().splitlines()
        if l.strip().startswith("include")
    ]
    return signed(
        {
            "archive_sha256": ARCHIVE_PIN,
            "inventory": inventory,
            "source_files": {n: sha256(v).hexdigest() for n, v in source.items()},
            "family": "separate_author_COMPASS_class2_example; not OPLS or pinned PCFF FRC",
            "license": "no license/permission file in this archive; article CC-BY does not establish archive redistribution",
            "source_version": "archive content hash; exact generating converter source/version not supplied",
            "settings": settings,
            "active_coefficient_rows": coefficient_rows,
            "atom_types_and_charges": types,
            "atom_count": len(atoms),
            "improper_center_types": dict(improper_centers),
            "silicon_neighbor_triples_without_source_improper_or_AA": len(absent),
            "missing_silicon_AA_patterns": [
                {
                    "center_type": 44,
                    "arm_types": list(k),
                    "source_triples": v,
                    "reason": "author topology omits Si-centered AA; no source absence rule or transferable row",
                }
                for k, v in sorted(missing_patterns.items())
            ],
            "missing_source_sites_checksum": sha256(
                json.dumps(absent, sort_keys=True).encode()
            ).hexdigest(),
            "charge_file_commands": len(charge_commands),
            "main_includes": includes,
            "missing_requirements": [
                "source-authorized transferable Si-centered angle-angle rows or explicit applicability rule",
                "exact generating converter rules/provenance for literal exported zero couplings",
                "source-bound atom typing and terminal charge transfer rule for new bounded final graphs",
                "source authorization/license scope for this parameter archive",
            ],
            "parameter_completeness": "incomplete",
            "registered_profiles": [],
            "cross_engine_equivalence": "not_claimed",
            "MD_executed": False,
        }
    )


def independent_frc(system, report, source):
    """Reconstruct every selected row, graph path and charge without production parsing."""
    raw = Path(source).read_bytes()
    if sha256(raw).hexdigest() != report["source_hashes"]["pcff.frc"]:
        raise ValueError("source_hash_mutation")
    lines = raw.decode().splitlines()
    parsed = read_source(raw)
    labels = {int(i): t for i, t in report["typing"].items()}
    adj = defaultdict(list)
    for b in system.topology.bonds.values():
        adj[b.site1].append(b.site2)
        adj[b.site2].append(b.site1)
    for i, atom in system.topology.sites.items():
        expected = {"Si": "sio", "O": "osi", "C": "c3"}.get(atom.element)
        if atom.element == "H":
            expected = (
                "hsi" if system.topology.sites[adj[i][0]].element == "Si" else "hc"
            )
        assert labels[i] == expected
    for row in report["source_rows"].values():
        line = lines[row["location"]["line"] - 1]
        assert line == row["raw_text"]
        w = line.split("!", 1)[0].split()
        if row["family"] in ARITY:
            n = ARITY[row["family"]]
            assert (
                w[2 : 2 + n] == row["labels"] and w[2 + n :] == row["original_values"]
            )
        elif row["family"] == "atom_typing_mass":
            assert w[3:6] == row["original_values"]
        elif row["family"] == "equivalence":
            assert w[3:] == row["original_values"]
    paths = {
        f + ":" + ",".join(map(str, sites))
        for f, sites in raw_inventory(system)
        if f != "wilson_out_of_plane"
    }
    inventory = {r["id"]: r for r in report["interaction_inventory"]}
    atom_paths = {
        f"{kind}:{i}" for kind in ("atom_type", "mass", "charge") for i in labels
    }
    assert paths | atom_paths == set(inventory)
    for rid, req in inventory.items():
        if rid in atom_paths:
            continue
        q = strict_resolution(parsed, req["family"], [labels[i] for i in req["sites"]])
        assigned = req["assignment"]
        assert (q["status"] == "assigned") == (assigned is not None), rid
        if assigned:
            assert q["row_ids"] == assigned["row_ids"]
            assert q["values"] == assigned["values"]
    charges = charge_raw(parsed, system, labels)
    assert charges["complete"] and {
        i: float(v) for i, v in charges["partial_charges"].items()
    } == {int(i): v for i, v in report["charge_vector"].items()}
    return {
        "status": "passed",
        "raw_source_sha256": sha256(raw).hexdigest(),
        "source_rows_reconstructed": len(report["source_rows"]),
        "graph_requests_reconstructed": len(paths),
        "charge_decimal_e": charges["partial_charges"],
        "component_totals": charges["components"],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.write_text(
        json.dumps(archive_audit(a.archive), indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
