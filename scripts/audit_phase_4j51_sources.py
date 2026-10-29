"""Independent offline ZIP/PDF/engine-line reconstruction, outside production parsers."""

import argparse
import json
import math
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile

from scripts.audit_phase_4j50_rows import audit as numerical_audit


def check(value, message):
    if not value:
        raise ValueError(message)


def pin(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def reconstruct_partitions(directory):
    """Raw pressure tensor and graph-distance 1-4 equations, no force-field parser."""
    c = json.loads((directory / "components.json").read_text())
    d = json.loads((directory / "detailed-components.json").read_text())
    lines = (directory / "lammps/comparison/engine.screen").read_text().splitlines()
    i = next(
        i
        for i, t in enumerate(lines)
        if t.strip().startswith("Step ") and "E_bond" in t
    )
    pressure = list(map(float, lines[i + 1].split()))[9:]
    virial = [v * 40**3 * 4.184 / 68568.415 for v in pressure]
    check(
        pressure == d["pressure_proxy"]
        and virial
        == d["LAMMPS_virial"]
        == c["LAMMPS"]["virial_kj_mol_xx_yy_zz_xy_xz_yz"],
        "Raw pressure/virial conversion changed",
    )
    assignment = json.loads((directory / "canonical-assignment.json").read_text())
    graph = assignment["final_graph"]
    adjacent = {a: set() for a in graph["atom_order"]}
    for a, b, *_ in graph["bonds"]:
        adjacent[a].add(b)
        adjacent[b].add(a)
    required = set()
    for a in adjacent:
        distances, frontier = {a: 0}, [a]
        for atom in frontier:
            if distances[atom] == 3:
                continue
            for b in sorted(adjacent[atom]):
                if b not in distances:
                    distances[b] = distances[atom] + 1
                    frontier.append(b)
        required.update((a, b) for b, n in distances.items() if a < b and n == 3)
    pairs = d["special_pairs_independently_reconstructed"]["pairs"]
    check(
        {tuple(v) for v in pairs} == required and len(pairs) == len(required),
        "Special-pair topology changed",
    )
    dump = (directory / "lammps/comparison/state.dump").read_text().splitlines()
    i = next(i for i, t in enumerate(dump) if t.startswith("ITEM: ATOMS"))
    rows = [list(map(float, t.split())) for t in dump[i + 1 :]]
    atom = {site: row for site, row in zip(graph["atom_order"], rows)}
    charges = dict(assignment["graph_charges"]["vector"])
    check(
        all(float(charges[a]) == atom[a][3] for a in atom),
        "Dump charge mapping changed",
    )
    nb = {}
    for r in assignment["source_evidence"]["rows"]:
        if r["family"] == "nonbonded":
            nb[r["author_type_id"]] = [
                float(v) for v in r["settings_evidence"]["raw_text"].split()[3:5]
            ]
    types = {
        r["site_id"]: r["author_type_id"]
        for r in assignment["typed_graph"]["atom_types"]
    }
    coulomb, lj = 0.0, 0.0
    for a, b in pairs:
        delta = [atom[a][k] - atom[b][k] for k in (4, 5, 6)]
        delta = [v - 40 * round(v / 40) for v in delta]
        r = math.sqrt((delta[0] ** 2 + delta[1] ** 2) + delta[2] ** 2)
        coulomb += 0.5 * 1389.35457644382 * atom[a][3] * atom[b][3] / r
        ea, sa = nb[types[a]]
        eb, sb = nb[types[b]]
        if r < 11:
            switch = (
                1 if r <= 9 else (121 - r * r) ** 2 * (121 + 2 * r * r - 243) / 40**3
            )
            lj += (
                0.5
                * 4
                * 4.184
                * math.sqrt(ea * eb)
                * ((math.sqrt(sa * sb) / r) ** 12 - (math.sqrt(sa * sb) / r) ** 6)
                * switch
            )
    s = d["special_pairs_independently_reconstructed"]
    # Arithmetic reconstruction roundoff only; numerical acceptance bounds are unchanged.
    check(
        abs(coulomb - s["coulomb14_exact_equation_kj_mol"]) <= 1e-10
        and abs(lj - s["lj14_kj_mol"]) <= 1e-10,
        "Independent 1-4 equations differ",
    )
    return {
        "special_pair_count": len(pairs),
        "Coulomb14_kj_mol": coulomb,
        "LJ14_kj_mol": lj,
        "native_pressure_proxy_reconstructed": True,
        "pressure_scientifically_validated": False,
    }


def review(reference_directory):
    from pypdf import PdfReader, __version__

    # Re-read the relevant original sources, not a new blind inventory.
    spec = (
        (
            "jp4c08471.pdf",
            "10.1021/acs.jpcb.4c08471",
            "selected author paper; Qeq results do not validate DFT/NPA",
        ),
        (
            "jp4c08471_si_001.pdf",
            "10.1021/acs.jpcb.4c08471",
            "OPLS Table S4 lj/cut instruction differs from bound archive; DFT scope unproven",
        ),
        (
            "ma00107a006.pdf",
            "10.1021/ma00107a006",
            "Sun polysilanes CFF93; Si-Si backbone, not OPLS siloxane execution authority",
        ),
        (
            "COMPASS, An ab Initio Force-Field Optimized forCondensed-Phase Applications.pdf",
            "10.1021/jp980939v",
            "Sun COMPASS Class-II functional form; no import into external OPLS",
        ),
        (
            "J Comput Chem - July 1994 - Sun - Force field for computation of conformational energies  structures  and vibrational.pdf",
            "10.1002/jcc.540150708",
            "Sun aromatic polyesters CFF93; no OPLS DFT branch convention authority",
        ),
        (
            "jp047434r.pdf",
            "10.1021/jp047434r",
            "Smith/Borodin/Smith independent PDMS model; no cross-family convention borrowing",
        ),
    )
    result = []
    for name, doi, role in spec:
        path = reference_directory / name
        pages = [page.extract_text() or "" for page in PdfReader(path).pages]
        citations = [
            {"page": 1, "line": 1, "raw_text": "\n".join(pages[0].splitlines()[:8])}
        ]
        if name == "jp4c08471_si_001.pdf":
            raw = "pair_style lj/cut/coul/long 11 11"
            check(raw in pages[5], "SI Table S4 changed")
            citations.append({"page": 6, "table": "S4", "raw_text": raw})
        if name == "ma00107a006.pdf":
            check("CFF93" in pages[2], "Sun functional-form context lost")
            line = next(i for i, t in enumerate(pages[2].splitlines()) if "CFF93" in t)
            citations.append(
                {"page": 3, "line": line + 1, "raw_text": pages[2].splitlines()[line]}
            )
        result.append(
            {
                "source_file": name,
                "source_sha256": pin(path),
                "doi": doi,
                "pages_read": len(pages),
                "role": role,
                "citations": citations,
                "parameter_rows_imported": 0,
                "pypdf_version": __version__,
                "extracted_text_sha256": sha256("\n".join(pages).encode()).hexdigest(),
            }
        )
    return result


def audit(args):
    from pypdf import PdfReader

    root = args.directory
    numerical = numerical_audit(
        root,
        args.archive,
        args.archive_sha256,
        args.lammps_source,
        args.reference_directory / "jp4c08471.pdf",
    )
    reviews = review(args.reference_directory)
    build_path = args.repo / "docs/evidence/phase_4j49/kernel-build.json"
    build = json.loads(build_path.read_text())
    modifications = []
    for file, change in build["changes"].items():
        original = (args.lammps_source / "src" / file).read_bytes()
        modified = (args.exact_kernel_source / "src" / file).read_bytes()
        check(
            sha256(original).hexdigest() == change["original_sha256"]
            and sha256(modified).hexdigest() == change["modified_sha256"]
            and original.count(change["old_expression"].encode())
            == change["occurrences"]
            and original.replace(
                change["old_expression"].encode(), change["replacement"].encode()
            )
            == modified,
            "Exact-kernel change not independently reconstructed",
        )
        modifications.append({"source_file": file, **change})
    pdfs = {
        r["source_file"]: PdfReader(args.reference_directory / r["source_file"])
        for r in reviews
    }
    with ZipFile(args.archive) as z:
        raw = {
            n.split("/", 1)[-1]: z.read(n) for n in z.namelist() if not n.endswith("/")
        }
    backend_lines, source_lines, absent, converted = 0, 0, 0, 0
    partitions = {}
    for dp in (0, 2, 4, 6):
        for policy in ("source_replay_v1", "exact_erfc_diagnostic_v1"):
            p = json.loads((root / f"dp{dp}" / policy / "resolution.json").read_text())
            partitions[f"dp{dp}/{policy}"] = reconstruct_partitions(
                root / f"dp{dp}" / policy
            )
            check(len(p["ten_convention_resolutions"]) == 10, "Lost ten-item table")
            for row in p["rows"]:
                for claim in row["candidates"]:
                    for c in claim["evidence"]:
                        file = c["source_file"]
                        if file in raw:
                            b = raw[file]
                        elif file in pdfs:
                            check(
                                pin(args.reference_directory / file)
                                == c["source_sha256"]
                                and c["raw_text"]
                                in pdfs[file].pages[c["page"] - 1].extract_text(),
                                "PDF source citation changed",
                            )
                            source_lines += 1
                            continue
                        else:
                            b = (args.lammps_source / file).read_bytes()
                        check(
                            sha256(b).hexdigest() == c["source_sha256"]
                            and b.decode().splitlines()[c["line"] - 1] == c["raw_text"],
                            "Convention source line changed: " + file,
                        )
                        source_lines += 1
                if row["resolution"] == "unresolved":
                    check(
                        row["value"] is None and row["missing_evidence"],
                        "Invented convention",
                    )
                    absent += 1
                if row["setting"] == "lj_tail_correction":
                    check(
                        row["resolution"] == "converter_resolved"
                        and row["value"] == "no"
                        and not row["source_authoritative"],
                        "Promoted tail authority",
                    )
                    text = (args.lammps_source / "doc/src/pair_charmm.rst").read_text()
                    check(
                        "None of the *lj/charmm*" in text and "tail option" in text,
                        "Missing explicit functional-form restriction",
                    )
                    converted += 1
                for c in row["backend"]["evidence"]:
                    if c.get("authority") != "backend_only":
                        continue
                    b = (args.lammps_source / c["source_file"]).read_bytes()
                    check(
                        sha256(b).hexdigest() == c["source_sha256"]
                        and b.decode().splitlines()[c["line"] - 1] == c["raw_text"],
                        "Backend citation changed",
                    )
                    backend_lines += 1
            scope = p["source_scope_review"]["evidence"]
            check(
                scope["raw_text"]
                in pdfs[scope["source_file"]].pages[scope["page"] - 1].extract_text(),
                "Missing competing SI instruction",
            )
    prior = []
    for phase in (47, 48, 49, 50):
        for path in sorted(
            (args.repo / f"docs/evidence/phase_4j{phase}").rglob("*.json")
        ):
            if path.name in (
                "cases.json",
                "result.json",
                "summary.json",
                "qualification.json",
                "historical-control.json",
            ) or any(
                t in path.name
                for t in (
                    "component",
                    "comparison",
                    "diagnosis",
                    "reconciliation",
                    "convention-table",
                )
            ):
                v = json.loads(path.read_text())
                prior.append(
                    {
                        "file": str(path.relative_to(args.repo)),
                        "sha256": pin(path),
                        "top_level_fields": sorted(v),
                    }
                )
    return {
        "schema": "island.j51-independent-reconstruction.v1",
        "outside_production_parser": True,
        "numerical": numerical,
        "source_reviews": reviews,
        "prior_comparisons_read": prior,
        "backend_citation_lines": backend_lines,
        "convention_citation_lines": source_lines,
        "unresolved_rows_preserved": absent,
        "no_tail_rule_reconstructed": converted,
        "physical_parameters_or_charges_imported": False,
        "exact_kernel_independently_reconstructed": modifications,
        "kernel_build_receipt_sha256": pin(build_path),
        "independent_special_pairs_and_pressure": partitions,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--reference-directory", type=Path, required=True)
    p.add_argument("--lammps-source", type=Path, required=True)
    p.add_argument("--exact-kernel-source", type=Path, required=True)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(argv)
    value = audit(args)
    args.output.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {k: v for k, v in value.items() if isinstance(v, (str, int, bool))},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
