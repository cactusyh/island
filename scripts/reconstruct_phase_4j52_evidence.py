"""Independent source, graph, charge, FD and native-pressure reconstruction.

No production force-field parser is imported. This uses raw source bytes and
the historical independent envelope/section readers, not the assignment API.
"""

import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from zipfile import ZipFile

from scripts.audit_phase_4j47_inputs import DFT, TYPES, read, sections
from scripts.audit_phase_4j50_rows import audit as numerical_audit
from scripts.audit_phase_4j51_sources import reconstruct_partitions
from scripts.audit_phase_4j52_sources import check, corpus, digest, signed


def reconstruct_mapping(directory, sources):
    p = json.loads((directory / "canonical-assignment.json").read_text())
    s = read(directory / "system.json")
    g = p["final_graph"]
    check(p["graph_identity"] == digest(g), "Final graph identity")
    check(g["atom_order"] == [a["id"] for a in s["sites"]], "Final atom order")
    check(
        g["bonds"]
        == sorted(
            [[b["site1"], b["site2"], b["order"], b["aromatic"]] for b in s["bonds"]]
        ),
        "Final bonds",
    )
    adjacent = {i: set() for i in g["atom_order"]}
    for a, b, order, aromatic in g["bonds"]:
        check(order == 1 and not aromatic, "Unsupported bonding")
        adjacent[a].add(b)
        adjacent[b].add(a)
    atoms = {a["id"]: a for a in s["sites"]}
    element = {i: a["element"] for i, a in atoms.items()}
    ends = {
        i
        for i in atoms
        if element[i] == "Si"
        and Counter(element[j] for j in adjacent[i]) == {"C": 3, "O": 1}
    }
    check(len(ends) == 2, "Terminal graph chemistry")
    typed = p["typed_graph"]
    oxygen = typed["e_oxygen_site_ids"]
    if oxygen:
        check(
            len(oxygen) == 1
            and element[oxygen[0]] == "O"
            and adjacent[oxygen[0]] & ends,
            "Terminal oxygen anchor",
        )
        end = ends | set(oxygen)
    else:
        check(sum(e == "O" for e in element.values()) == 1, "HMDS oxygen")
        end = ends | {i for i in atoms if element[i] == "O"}
    for si in ends:
        for carbon in adjacent[si]:
            if element[carbon] == "C":
                end.add(carbon)
                end.update(h for h in adjacent[carbon] if element[h] == "H")
    data = sections(sources["param/OPLS/opls_system.data"].decode())
    rows = p["source_evidence"]["rows"]
    check(len(rows) == 17, "Parameter rows incomplete")
    for row in rows:
        for key in ("data_evidence", "settings_evidence"):
            cite = row[key]
            raw = sources[cite["source_file"]]
            check(
                sha256(raw).hexdigest() == cite["sha256"]
                and raw.decode().splitlines()[cite["line"] - 1] == cite["raw_text"],
                "Raw parameter citation",
            )
        words = row["settings_evidence"]["raw_text"].split()
        raw_values = words[3:] if row["family"] == "nonbonded" else words[2:]
        check(
            [Fraction(v) for v in raw_values] == [Fraction(v) for v in row["values"]],
            "Raw numeric parameter values",
        )
    charges = dict(p["graph_charges"]["vector"])
    serialized = {t: set() for t in TYPES.values()}
    for _, _, row in data["Atoms"].values():
        serialized[int(row[1])].add(row[2])
    for row in typed["atom_types"]:
        i = row["site_id"]
        expected = format(float(DFT[element[i]][int(i in end)]), ".16e")
        check(
            row["author_type_id"] == TYPES[element[i]]
            and row["fragment"] == ("E" if i in end else "M")
            and charges[i] == expected
            and expected in serialized[TYPES[element[i]]],
            "Independent type/E-M/serialized charge mapping",
        )
        check(
            Fraction(str(atoms[i]["mass"]))
            == Fraction(data["Masses"][TYPES[element[i]]][2][0]),
            "Exact source mass",
        )
    check(
        p["graph_charges"]["component_totals"]
        == [
            str(sum((Fraction(charges[i]) for i in c), Fraction()))
            for c in g["components"]
        ],
        "Exact rational component charge totals",
    )
    paths = {
        "bond": {min((a, b), (b, a)) for a in adjacent for b in adjacent[a]},
        "angle": {
            min((a, b, c), (c, b, a))
            for b in adjacent
            for a, c in combinations(sorted(adjacent[b]), 2)
        },
        "proper_torsion": {
            min((a, b, c, d), (d, c, b, a))
            for b in adjacent
            for c in adjacent[b]
            for a in adjacent[b] - {c}
            for d in adjacent[c] - {b, a}
        },
    }
    by_hash = {r["evidence_checksum"]: r for r in rows}
    for family, expected in paths.items():
        assigned = p["assigned_terms"][family]
        check(
            {tuple(r["sites"]) for r in assigned} == expected
            and len(assigned) == len(expected),
            "Complete independent interaction enumeration",
        )
        for term in assigned:
            row = by_hash[term["row_checksum"]]
            labels = tuple(element[i] for i in term["sites"])
            check(
                row["family"] == family
                and tuple(row["labels"]) == min(labels, labels[::-1]),
                "Parameter-to-interaction labels",
            )
    exported = sections((directory / "lammps/system.data").read_text())
    for n, i in enumerate(g["atom_order"], 1):
        row = exported["Atoms"][n][2]
        check(
            row[2] == charges[i] and list(map(float, row[3:])) == s["coordinates"][i],
            "Exact exported charge/coordinate ordering",
        )
    return {
        "source_rows": len(rows),
        "atoms": len(atoms),
        "graph_type_mass_charge_and_parameter_rows_reconstructed": True,
        "interaction_counts": {k: len(v) for k, v in paths.items()},
    }


def reconstruct(args):
    searched = corpus(args)
    root = args.directory
    check(
        searched == json.loads((root / "supplied-corpus.json").read_text()),
        "Independent complete-corpus replay differs",
    )
    numerical = numerical_audit(
        root,
        args.archive,
        searched["source_archive_sha256"],
        args.backend,
        args.references / "jp4c08471.pdf",
    )
    with ZipFile(args.archive) as z:
        sources = {
            n.split("/", 1)[1]: z.read(n) for n in z.namelist() if not n.endswith("/")
        }
    mapping, partitions, lines = {}, {}, 0
    final_table = []
    for result in searched["final_resolutions"]:
        raw = []
        for asset in searched["assets"]:
            raw.extend(
                {
                    k: asset[k]
                    for k in (
                        "source_file",
                        "source_sha256",
                        "doi",
                        "citation",
                        "scope",
                    )
                }
                | cite
                for cite in asset["matches"][result["setting"]]
            )
        final_table.append(
            signed(
                {
                    **{k: v for k, v in result.items() if k != "identity"},
                    "corpus_result_identity": result["identity"],
                    "raw_evidence": raw,
                }
            )
        )
    for dp in (0, 2, 4, 6):
        for policy in ("source_replay_v1", "exact_erfc_diagnostic_v1"):
            directory = root / f"dp{dp}" / policy
            key = f"dp{dp}/{policy}"
            mapping[key] = reconstruct_mapping(directory, sources)
            partitions[key] = reconstruct_partitions(directory)
            final = json.loads((directory / "final-resolution.json").read_text())
            check(
                final["source_audit_identity"] == searched["identity"],
                "Final corpus identity",
            )
            for row in final["final_resolutions"]:
                check(
                    row["resolution"] == "unresolved_permanently"
                    and row["value"] is None,
                    "Invented final author convention",
                )
                for cite in row["backend"]["evidence"]:
                    raw = (args.backend / cite["source_file"]).read_bytes()
                    check(
                        sha256(raw).hexdigest() == cite["source_sha256"]
                        and raw.decode().splitlines()[cite["line"] - 1]
                        == cite["raw_text"],
                        "Independent backend definition citation",
                    )
                    lines += 1
    build = json.loads(
        (args.repo / "docs/evidence/phase_4j49/kernel-build.json").read_text()
    )
    for name, change in build["changes"].items():
        original = (args.backend / "src" / name).read_bytes()
        modified = (args.exact_source / "src" / name).read_bytes()
        check(
            sha256(original).hexdigest() == change["original_sha256"]
            and sha256(modified).hexdigest() == change["modified_sha256"]
            and original.replace(
                change["old_expression"].encode(), change["replacement"].encode()
            )
            == modified,
            "Independent exact-erfc edit reconstruction",
        )
    prior_raw = {}
    prior_errors = {}
    for phase in range(47, 52):
        for path in sorted(
            (args.repo / f"docs/evidence/phase_4j{phase}/audit").rglob("*")
        ):
            if not path.is_file():
                continue
            raw = path.read_bytes()
            name = str(path.relative_to(args.repo))
            prior_raw[name] = sha256(raw).hexdigest()
            if path.suffix == ".json":
                json.loads(raw)
            if path.suffix in (".screen", ".log"):
                errors = [
                    {"line": n, "raw_text": line}
                    for n, line in enumerate(raw.decode().splitlines(), 1)
                    if line.startswith("ERROR:")
                ]
                if errors:
                    prior_errors[name] = errors
    table = signed(
        {
            "schema": "island.j52-final-nine-table.v1",
            "source_audit_identity": searched["identity"],
            "rows": final_table,
        }
    )
    (root / "final-nine-conventions.json").write_text(
        json.dumps(table, indent=2, sort_keys=True) + "\n"
    )
    return signed(
        {
            "schema": "island.j52-independent-reconstruction.v1",
            "outside_production_parser": True,
            "source_audit_identity": searched["identity"],
            "assets_independently_reconstructed": len(searched["assets"]),
            "final_convention_raw_evidence_rows": sum(
                len(r["raw_evidence"]) for r in final_table
            ),
            "numerical": numerical,
            "independent_graph_mapping": mapping,
            "independent_special_pairs_and_pressure": partitions,
            "backend_citation_lines_reconstructed": lines,
            "exact_kernel_edits_reconstructed": len(build["changes"]),
            "J47_J51_raw_files_rechecked": prior_raw,
            "historical_engine_errors_preserved": prior_errors,
            "MD_executed": False,
            "NPT_executed": False,
        }
    )


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in (
        "directory",
        "archive",
        "references",
        "backend",
        "exact-source",
        "output",
    ):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    args = p.parse_args(argv)
    value = reconstruct(args)
    args.output.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "identity": value["identity"],
                "assets_reconstructed": value["assets_independently_reconstructed"],
                "raw_FD_samples": value["numerical"]["raw_FD_samples_reconstructed"],
                "prior_raw_files_read": len(value["J47_J51_raw_files_rechecked"]),
                "backend_citations": value["backend_citation_lines_reconstructed"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
