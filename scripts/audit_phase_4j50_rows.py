"""Independent reconstruction of source citations and raw numerical samples.

Does not import the production parser, convention catalogue or FD driver.
pypdf is used only for explicit local PDF citation reconstruction.
"""

import argparse
import json
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile


def check(condition, reason):
    if not condition:
        raise ValueError(reason)


def audit(root, archive, expected_sha256, lammps_source, paper):
    from pypdf import PdfReader

    raw = Path(archive).read_bytes()
    check(sha256(raw).hexdigest() == expected_sha256, "Archive mutation")
    with ZipFile(archive) as z:
        sources = {
            n.split("/", 1)[-1]: z.read(n) for n in z.namelist() if not n.endswith("/")
        }
    inventory = {n: sha256(b).hexdigest() for n, b in sources.items()}
    paper_hash = sha256(Path(paper).read_bytes()).hexdigest()
    pdf = PdfReader(paper)
    root, lammps_source = Path(root), Path(lammps_source)
    rows, pdf_rows, samples, cases = 0, 0, 0, {}
    for dp in (0, 2, 4, 6):
        case = {}
        for policy in ("source_replay_v1", "exact_erfc_diagnostic_v1"):
            directory = root / f"dp{dp}" / policy
            p = json.loads((directory / "conventions.json").read_text())
            source = p["source_authoritative_settings"]
            check(
                source["archive_file_inventory"] == inventory,
                "Source inventory mutation",
            )
            for r in source["settings"].values():
                for cite in r["evidence"]:
                    if "source_file" in cite:
                        if cite.get("source_medium") == "PDF":
                            check(
                                paper_hash == cite["source_sha256"]
                                and cite["raw_text"]
                                in pdf.pages[cite["page"] - 1].extract_text(),
                                "PDF raw source citation changed",
                            )
                            pdf_rows += 1
                            continue
                        data = sources[cite["source_file"]]
                        check(
                            sha256(data).hexdigest() == cite["source_sha256"],
                            "Source citation hash",
                        )
                        check(
                            data.decode().splitlines()[cite["line"] - 1]
                            == cite["raw_text"],
                            "Source raw line",
                        )
                        rows += 1
                if r["classification"] == "unresolved":
                    check(
                        r["value"] is None and r["missing_evidence"],
                        "Invented source authority",
                    )
                    check(
                        all(c["scope"].startswith("not_linked") for c in r["evidence"]),
                        "Cross-scope import",
                    )
            for file, e in p["converter_settings"]["engine_definitions"].items():
                data = (lammps_source / file).read_bytes()
                check(
                    sha256(data).hexdigest() == e["sha256"],
                    "Engine definition mutation: " + file,
                )
                check(
                    all(
                        data.decode().splitlines()[c["line"] - 1] == c["raw_text"]
                        for c in e["citations"]
                    ),
                    "Engine raw line",
                )
            c = json.loads((directory / "components.json").read_text())
            fd = c["finite_difference"]
            check(
                fd["bound"] == 1e-4 and c["energy_bound"] == c["force_bound"] == 3e-4,
                "Changed numerical limit",
            )
            values = {}
            log = (
                (directory / "lammps/comparison/engine.screen").read_text().splitlines()
            )
            for line in log:
                if line.startswith("J49FD ") and "$(" not in line:
                    words = line.split()
                    e = [float(t) * 4.184 for t in words[2:]]
                    values[int(words[1])] = {
                        "bonded": (e[0] + e[1]) + e[2],
                        "lj": e[3],
                        "coulomb": e[4] + e[5],
                        "total": e[6],
                    }
            check(len(values) == 2 * len(fd["raw_records"]), "Missing raw FD energies")
            maxima = dict.fromkeys(("bonded", "lj", "coulomb", "total"), 0.0)
            for index, record in enumerate(fd["raw_records"]):
                for component, v in record["components"].items():
                    plus, minus = (
                        values[2 * index][component],
                        values[2 * index + 1][component],
                    )
                    force = -(plus - minus) / (2 * record["h_angstrom"])
                    error = abs(force - v["analytical_force"])
                    check(
                        plus == v["positive_energy_kj_mol"]
                        and minus == v["negative_energy_kj_mol"]
                        and force == v["central_difference_force"]
                        and error == v["absolute_error"],
                        "Raw FD equation or energy changed",
                    )
                    maxima[component] = max(maxima[component], error)
                samples += 1
            check(
                maxima == fd["max_error_by_component"]
                and fd["passed"] == all(v <= 1e-4 for v in maxima.values()),
                "FD classification mutation",
            )
            dump = (directory / "lammps/comparison/state.dump").read_text().splitlines()
            start = next(
                i for i, line in enumerate(dump) if line.startswith("ITEM: ATOMS")
            )
            force = [
                [float(t) * 4.184 for t in line.split()[7:]]
                for line in dump[start + 1 :]
            ]
            check(
                force == c["LAMMPS"]["forces_kj_mol_angstrom"],
                "Native force dump mutation",
            )
            case[policy] = {
                "raw_FD_samples": len(fd["raw_records"]),
                "FD_passed": fd["passed"],
                "native_force_dump_verified": True,
            }
        cases[str(dp)] = case
    return {
        "schema": "island.j50-independent-reconstruction.v1",
        "outside_production_parser": True,
        "source_archive_sha256": expected_sha256,
        "source_files": len(inventory),
        "source_citations_reconstructed": rows,
        "PDF_citations_reconstructed": pdf_rows,
        "paper_sha256": paper_hash,
        "raw_FD_samples_reconstructed": samples,
        "cases": cases,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--lammps-source", type=Path, required=True)
    p.add_argument("--paper", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(argv)
    value = audit(a.directory, a.archive, a.archive_sha256, a.lammps_source, a.paper)
    a.output.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    print(json.dumps(value, sort_keys=True))


if __name__ == "__main__":
    main()
