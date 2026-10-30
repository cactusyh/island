"""Offline exhaustive pinned-corpus search; independent of production FF parsers.

Absence is scoped to these exact supplied bytes, not a claim about future evidence.
The PDF reader is an explicit audit dependency, never a runtime discovery service.
"""

import argparse
import json
import re
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile

NINE = (
    "dielectric",
    "coulomb_prefactor",
    "periodic_image_policy",
    "precision",
    "neighbor_list_policy",
    "pressure_convention",
    "virial_convention",
    "pair_table_bits",
    "erfc_kernel",
)
QUERIES = {
    "dielectric": r"dielectric|permittivity",
    "coulomb_prefactor": r"Coulomb|qqr[de2]|332\.06|1389\.35|epsilon_?0|epsilon_0",
    "periodic_image_policy": r"periodic|boundary|slab|conducting|tin.foil|surface term",
    "precision": r"precision|double|single|LAMMPS|version|compiler|build",
    "neighbor_list_policy": r"neigh\w*|atom_modify.*sort",
    "pressure_convention": r"pressure|thermo_style|compute.*press|virial|pxy|pxz|pyz",
    "virial_convention": r"virial|compute.*press|pressure tensor|pxy|pxz|pyz",
    "pair_table_bits": r"pair_modify|table.*bits|interpolat|ncoultablebits",
    "erfc_kernel": r"erfc|EWALD_P|EWALD_F|A1\s*=|approx.*error function",
}
MISSING = {
    "dielectric": "DFT/NPA-scoped dielectric declaration or executable/restart proving its value",
    "coulomb_prefactor": "author LAMMPS revision/build and its real-unit electrostatic constant",
    "periodic_image_policy": "DFT/NPA-scoped surface/slab/image-sum policy; x/y/z PBC alone is insufficient",
    "precision": "author executable hash, build precision, accelerator and compiler manifest",
    "neighbor_list_policy": "complete DFT/NPA driver linking the generic neighbor/skin/rebuild/sort settings",
    "pressure_convention": "DFT/NPA pressure compute, kinetic/constraint terms and thermodynamic estimator definition",
    "virial_convention": "DFT/NPA virial tensor, image/constraint/kinetic contributions and conversion definition",
    "pair_table_bits": "author pair_modify table setting or exact executable proving interpolation semantics",
    "erfc_kernel": "author engine revision/binary identifying erfc energy and force implementation",
}
OBSERVABLES = {
    "dielectric": [
        "Coulomb real",
        "Coulomb reciprocal",
        "Coulomb14",
        "forces",
        "virial",
        "pressure",
    ],
    "coulomb_prefactor": [
        "Coulomb real",
        "Coulomb reciprocal",
        "Coulomb14",
        "forces",
        "virial",
        "pressure",
    ],
    "periodic_image_policy": ["Coulomb reciprocal", "forces", "virial", "pressure"],
    "precision": [
        "energies",
        "forces",
        "finite differences",
        "checkpoint/continuation bytes",
    ],
    "neighbor_list_policy": [
        "pair energies",
        "forces",
        "trajectory/continuation bytes",
    ],
    "pressure_convention": ["pressure", "NPT", "pressure proxy interpretation"],
    "virial_convention": ["virial", "pressure", "NPT"],
    "pair_table_bits": ["Coulomb real", "forces", "finite differences"],
    "erfc_kernel": ["Coulomb real", "forces", "finite differences"],
}
REASONS = {
    "dielectric": "The only explicit archive dielectric belongs to Dreiding, not the selected OPLS include closure; a backend default would choose an unproven electrostatic model.",
    "coulomb_prefactor": "The source selects real units but does not pin the author engine constant; independently known present-engine constants do not identify the historical executable.",
    "periodic_image_policy": "Three-dimensional PBC is explicit, but it does not specify a surface dielectric, slab correction or electrostatic image convention.",
    "precision": "Present diagnostic double precision and output formatting cannot establish the author's binary/build or bitwise numerical contract.",
    "neighbor_list_policy": "Generic examples specify check yes; the selected OPLS entry point does not include them or specify skin and sorting. Selecting a neighbor policy changes numerical execution authority.",
    "pressure_convention": "Generic press/pxy output is not a complete DFT/NPA pressure estimator; the diagnostic configurational proxy cannot stand in for thermodynamic pressure.",
    "virial_convention": "A current backend tensor implementation is auditable but does not prove the author's image/kinetic/constraint convention; no scientific comparison bound exists.",
    "pair_table_bits": "Source inputs omit table bits. Present-engine interpolation defaults cannot be promoted to the author contract; J49 measured a numerical effect.",
    "erfc_kernel": "The author engine is unpinned. Exact-erfc diagnostic source edits are explicit experiments, not proof of the original implementation or permission to substitute it.",
}
CONVERTERS = (
    "src/island/forcefields/oplsaa/pdms_author.py",
    "src/island/forcefields/pdms_opls.py",
    "src/island/pdms.py",
    "src/island/pdms_lammps.py",
    "src/island/forcefields/pdms_lammps/__init__.py",
    "src/island/forcefields/oplsaa/pdms_conventions.py",
    "src/island/opls_execution.py",
    "src/island/opls_conventions/__init__.py",
    "src/island/opls_conventions/evidence.py",
    "scripts/build_phase_4j49_kernel.py",
)
BACKEND = (
    "doc/src/units.rst",
    "doc/src/dielectric.rst",
    "doc/src/pair_charmm.rst",
    "doc/src/pair_modify.rst",
    "doc/src/special_bonds.rst",
    "doc/src/dihedral_opls.rst",
    "doc/src/kspace_style.rst",
    "doc/src/kspace_modify.rst",
    "doc/src/boundary.rst",
    "doc/src/neighbor.rst",
    "doc/src/neigh_modify.rst",
    "doc/src/compute_pressure.rst",
    "src/update.cpp",
    "src/force.cpp",
    "src/pair.cpp",
    "src/kspace.cpp",
    "src/KSPACE/pppm.cpp",
    "src/KSPACE/pair_lj_charmm_coul_long.cpp",
    "src/ewald_const.h",
)
ARCHIVE_PIN = "8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4"
PDFS = {
    "jp4c08471.pdf": "04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc",
    "jp4c08471_si_001.pdf": "6e86aa64085466e286b5c59878c2288f5f16da1b7a5298d37ba1833a88cd602f",
    "ja11225.pdf": "451cf754daab3fcf6e357a3cdfb64e5d8ab8eb800ce8afdad53897c0bbedcc5a",
}


def canonical(v):
    return json.dumps(v, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(v):
    return sha256(canonical(v).encode()).hexdigest()


def signed(v):
    return {**v, "identity": digest(v)}


def check(v, message):
    if not v:
        raise ValueError(message)


def corpus(args):
    from pypdf import PdfReader, __version__

    assets = []
    patterns = {key: re.compile(query, re.IGNORECASE) for key, query in QUERIES.items()}

    def add(name, raw, kind, scope, pages=None):
        pages = pages if pages is not None else [raw.decode()]
        matches = {key: [] for key in NINE}
        for page, text in enumerate(pages, 1):
            for line, textline in enumerate(text.splitlines(), 1):
                for key, pattern in patterns.items():
                    if pattern.search(textline):
                        matches[key].append(
                            {
                                "page": page if kind == "PDF" else None,
                                "line": line,
                                "raw_text": textline,
                            }
                        )
        assets.append(
            {
                "source_file": name,
                "source_sha256": sha256(raw).hexdigest(),
                "size_bytes": len(raw),
                "kind": kind,
                "scope": scope,
                "citation": "Xiang et al., JPCB 129 (2025), 1864-1873"
                if kind in ("PDF", "archive_member") and name != "ja11225.pdf"
                else "Jorgensen, Maxwell, Tirado-Rives, JACS 118 (1996), 11225-11236"
                if name == "ja11225.pdf"
                else "Pinned local implementation; not author authority",
                "doi": "10.1021/acs.jpcb.4c08471"
                if kind in ("PDF", "archive_member") and name != "ja11225.pdf"
                else "10.1021/ja9621760"
                if name == "ja11225.pdf"
                else None,
                "pages_searched": len(pages),
                "lines_searched": sum(len(p.splitlines()) for p in pages),
                "extracted_text_sha256": sha256("\n".join(pages).encode()).hexdigest(),
                "matches": matches,
            }
        )

    archive = args.archive.read_bytes()
    check(sha256(archive).hexdigest() == ARCHIVE_PIN, "Author archive replacement")
    with ZipFile(args.archive) as z:
        names = [n for n in z.namelist() if not n.endswith("/")]
        check(len(names) == len(set(names)), "Duplicate ZIP member")
        for name in sorted(names):
            relative = name.split("/", 1)[1]
            scope = (
                "selected_DFT_NPA_include_closure"
                if relative.startswith("param/OPLS/")
                else "author_generic_not_linked_to_DFT_NPA"
                if relative.startswith("input/")
                else "other_family_or_readme_not_OPLS_authority"
            )
            add(relative, z.read(name), "archive_member", scope)
    for name, expected in PDFS.items():
        raw = (args.references / name).read_bytes()
        check(
            expected is None or sha256(raw).hexdigest() == expected, "PDF replacement"
        )
        pages = [
            page.extract_text() or ""
            for page in PdfReader(args.references / name).pages
        ]
        add(
            name,
            raw,
            "PDF",
            "selected_paper_or_SI"
            if name != "ja11225.pdf"
            else "original_OPLS_supporting_tables_no_PDMS_execution_authority",
            pages,
        )
    for name in CONVERTERS:
        add(
            name,
            (args.repo / name).read_bytes(),
            "converter",
            "current_pinned_converter_not_original_author_version",
        )
    for name in BACKEND:
        add(
            name,
            (args.backend / name).read_bytes(),
            "backend",
            "explicit_current_backend_definition_not_author_settings",
        )
    assets.sort(key=lambda r: (r["kind"], r["source_file"]))
    corpus_id = digest(assets)
    rows = []
    for key in NINE:
        searched = [
            {
                "source_file": a["source_file"],
                "source_sha256": a["source_sha256"],
                "scope": a["scope"],
                "match_count": len(a["matches"][key]),
            }
            for a in assets
        ]
        rows.append(
            signed(
                {
                    "setting": key,
                    "resolution": "unresolved_permanently",
                    "value": None,
                    "source_authoritative": False,
                    "source_corpus_identity": corpus_id,
                    "query": QUERIES[key],
                    "all_searched_evidence": searched,
                    "exact_missing_artifact": MISSING[key],
                    "affected_execution_modes": [
                        "research_opt_in",
                        "production_candidate",
                    ],
                    "affected_observables": OBSERVABLES[key],
                    "default_scientifically_unjustified": REASONS[key],
                    "raw_evidence_reference": "assets[].matches[setting]; all matches retained including empty results and rejected scopes",
                    "interpreted_value": None,
                    "units": "no author value; see separate backend interpretation",
                    "conversion_formula": "not applicable: no source value to convert",
                    "implementation_version": "island.j52-pinned-corpus-exhaustion.v1",
                    "independent_reconstruction": {
                        "method": "raw ZIP members, complete PDF pages, converter and backend bytes outside production parser",
                        "result": "all searched; no DFT/NPA-scoped declaration proves the missing artifact",
                    },
                    "evidence_checksum": digest(
                        {"corpus": corpus_id, "search": searched, "query": QUERIES[key]}
                    ),
                }
            )
        )
    # Review all prior comparisons by reading their actual bytes; no new claims.
    comparisons = {}
    for phase in range(47, 52):
        root = args.repo / f"docs/evidence/phase_4j{phase}"
        for p in sorted(root.rglob("*.json")):
            if any(
                t in p.name
                for t in (
                    "component",
                    "comparison",
                    "summary",
                    "result",
                    "case",
                    "kernel",
                    "receipt",
                )
            ):
                raw = p.read_bytes()
                json.loads(raw)
                comparisons[str(p.relative_to(args.repo))] = sha256(raw).hexdigest()
    return signed(
        {
            "schema": "island.j52-exhaustive-corpus.v1",
            "nine_conventions": list(NINE),
            "source_archive_sha256": ARCHIVE_PIN,
            "assets": assets,
            "corpus_identity": corpus_id,
            "final_resolutions": rows,
            "prior_comparison_records_read": comparisons,
            "pypdf_version": __version__,
            "no_network_discovery": True,
            "closure_scope": "Permanently unsupported by this exact pinned supplied corpus; new evidence requires a new audited version. Not a claim that future author evidence is impossible.",
            "source_scope_conflict": "SI p6 Table S4 uses lj/cut/coul/long 11 11; archive OPLS init line 8 uses lj/charmm/coul/long 9.0 11.0. No DFT/NPA branch reconciliation is supplied.",
            "original_author_converter_version": "not supplied; the pinned ISLAND converter is not proof of the author's Moltemplate version",
        }
    )


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--references", type=Path, required=True)
    p.add_argument("--backend", type=Path, required=True)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--verify", type=Path)
    a = p.parse_args(argv)
    value = corpus(a)
    if a.verify:
        check(
            value == json.loads(a.verify.read_text()),
            "Independent corpus replay differs",
        )
    a.output.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "identity": value["identity"],
                "assets": len(value["assets"]),
                "unresolved_permanently": list(NINE),
            }
        )
    )


if __name__ == "__main__":
    main()
