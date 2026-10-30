"""Reproducible, non-executable external PCFF/COMPASS PDMS Class-II audit."""

import argparse
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields.pcff.mixed_sip import DISTRIBUTION_PINS
from island.forcefields.selection.inventory import canonical, digest
from island.pcff_pdms import (
    CANDIDATE,
    PCFFPDMSAssets,
    prepare_polymer_forcefield,
    reconstruct,
    save,
)
from scripts.audit_pcff_pdms_class2_sources import archive_audit, independent_frc

BASE = "b296a06488abca891ae46511aba3a55b075beebf"
BRANCH = "codex/phase-4j54-external-pcff-pdms-class2"
REPEAT = "[*:1]O[Si](C)(C)[*:2]"


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def signed(value):
    return {**value, "identity": digest(value)}


def fixtures():
    yield "disiloxane", from_smiles("[SiH3]O[SiH3]", random_seed=2026)
    yield "hmds", from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)
    for dp in (2, 4, 6):
        yield (
            f"dp{dp}",
            build_linear_polymer(
                REPEAT,
                dp=dp,
                head_end_group="[*][Si](C)(C)C",
                tail_end_group="[*]C",
                add_hydrogens=True,
                random_seed=2026,
            ),
        )


def repeat_audit():
    from rdkit import Chem

    mol = Chem.MolFromSmiles(REPEAT)
    graph = {
        "psmiles": REPEAT,
        "atoms": [
            {
                "index": a.GetIdx(),
                "atomic_number": a.GetAtomicNum(),
                "attachment_map": a.GetAtomMapNum(),
            }
            for a in mol.GetAtoms()
        ],
        "bonds": [
            [b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()]
            for b in mol.GetBonds()
        ],
    }
    return signed(
        {
            "family": "pcff",
            "graph_descriptor": graph,
            "graph_identity": digest(graph),
            "candidate_profile": CANDIDATE,
            "selected_profile": None,
            "parameter_completeness": "incomplete",
            "typing_status": "incomplete",
            "charge_status": "incomplete",
            "missing_terms": [
                {
                    "family": "final_graph",
                    "reason": "open_repeat_attachments; explicit_DP_and_terminal_groups_required",
                    "attachment_maps": [1, 2],
                }
            ],
            "lammps_export_status": "blocked",
            "lammps_runtime_status": "blocked",
            "openmm_build_status": "blocked",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def audit(
    distribution, archive, output, *, fresh=True, references=None, pdf_python=None
):
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    assets = PCFFPDMSAssets(distribution)
    assets.load()
    cases, all_rows, grouped = [], {}, defaultdict(list)
    matrix = []
    relocated = root.parent / (root.name + "-relocated")
    if fresh:
        relocated.mkdir(parents=True, exist_ok=False)
        shutil.copytree(
            distribution, relocated / "sources", ignore=shutil.ignore_patterns("*.pdf")
        )
    for name, system in fixtures():
        result = prepare_polymer_forcefield(
            system, forcefield="pcff", assets=assets, audit_si_h=name == "disiloxane"
        )
        p = result.payload
        for family, coverage in sorted(p["coverage"].items()):
            matrix.append(
                {
                    "fixture": name,
                    "family": family,
                    "source_rows": coverage,
                    "dependencies": sorted(
                        {
                            d.split(":")[0]
                            for row in p["interaction_inventory"]
                            if row["family"] == family
                            for d in row["dependencies"]
                        }
                    ),
                }
            )
        case = root / name
        pin = save(system, result, case)
        independent = independent_frc(system, p, Path(distribution) / "pcff.frc")
        write(root / (name + "-independent.json"), independent)
        for r in p["missing_terms"]:
            grouped[r["family"], r["query"]].append(
                {
                    "fixture": name,
                    "id": r["id"],
                    "sites": r["sites"],
                    "dependencies": r["dependencies"],
                    "resolved_types": r["source_lookup"]["resolved_types"],
                }
            )
        all_rows.update(p["source_rows"])
        processes = []
        if fresh:
            shutil.copytree(case, relocated / name)
            for seed in (1, 99397):
                command = [
                    sys.executable,
                    "-m",
                    "scripts.phase_4j54",
                    "reconstruct",
                    "--directory",
                    str((relocated / name).resolve()),
                    "--manifest-sha256",
                    pin,
                    "--distribution",
                    str((relocated / "sources").resolve()),
                ]
                run = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=False,
                    env={**os.environ, "PYTHONHASHSEED": str(seed)},
                )
                if (
                    run.returncode
                    or json.loads(run.stdout)["identity"] != p["identity"]
                ):
                    raise ValueError("relocated reconstruction failed: " + run.stderr)
                processes.append(
                    {
                        "process_hash_seed": seed,
                        "returncode": run.returncode,
                        "record_identity": json.loads(run.stdout)["identity"],
                    }
                )
        case_summary = {
            k: p[k]
            for k in (
                "parameter_completeness",
                "typing_status",
                "charge_status",
                "lammps_export_status",
                "lammps_runtime_status",
                "openmm_build_status",
                "cross_engine_equivalence",
                "scientific_validation",
                "production_readiness",
                "coverage",
                "missing_conventions",
            )
        }
        case_summary.update(
            fixture=name,
            atoms=system.number_of_sites,
            report_identity=p["identity"],
            manifest_sha256=pin,
            missing_requests=len(p["missing_terms"]),
            graph_identity=p["graph_identity"],
            fresh_processes=processes,
        )
        cases.append(case_summary)
    write(root / "repeat.json", repeat_audit())
    author = archive_audit(archive)
    write(root / "author-compass.json", author)
    write(
        root / "source-to-row.json",
        signed(
            {
                "frc_rows": all_rows,
                "author_compass_rows": author["active_coefficient_rows"],
                "overlays_authorized": False,
            }
        ),
    )
    write(
        root / "missing-rows.json",
        signed(
            {
                "queries": [
                    {"family": f, "query": q, "requests": r}
                    for (f, q), r in sorted(grouped.items())
                ],
                "source_specific": "single pinned FRC only; not a statement that author COMPASS omits these numerical coefficients",
            }
        ),
    )
    if references is not None:
        from scripts.pcff_pdms_source_context import review

        write(
            root / "literature-context.json",
            review(Path.cwd(), references, pdf_python=pdf_python),
        )
    write(
        root / "class2-dependency-matrix.json",
        signed(
            {
                "rows": matrix,
                "additional_required_source_policies": cases[0]["missing_conventions"],
                "C-C": "not used by the bounded methyl-only graphs; no extrapolated C-C assignment",
                "improper": "structural Wilson applicability preserved; source omission authority unresolved",
                "torsion-torsion": "source-declared empty section; applicability unresolved",
            }
        ),
    )
    summary = signed(
        {
            "schema": "island.j54-audit.v1",
            "candidate_profile": CANDIDATE,
            "registered_profiles": [],
            "source_hashes": {
                **DISTRIBUTION_PINS,
                "PDMS_benchmark-main.zip": author["archive_sha256"],
            },
            "fixtures": cases,
            "open_repeat": "repeat.json",
            "missing_frc_requests": sum(c["missing_requests"] for c in cases),
            "missing_frc_source_queries": len(grouped),
            "independent_rows": len(all_rows),
            "author_COMPASS": {
                "parameter_completeness": "incomplete",
                "missing_requirements": author["missing_requirements"],
                "source_settings_not_transferable_to_FRC": True,
            },
            "source_selection": "FRC and author COMPASS assessed separately; Sun papers retained as historical independent source catalogs",
            "parameter_completeness_independent_of_cross_engine": True,
            "backend_execution": [],
            "MD_executed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "summary.json", summary)
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest="command", required=True)
    a = commands.add_parser("audit")
    a.add_argument("--distribution", type=Path, required=True)
    a.add_argument("--archive", type=Path, required=True)
    a.add_argument("--output", type=Path, required=True)
    a.add_argument("--references", type=Path, required=True)
    a.add_argument("--pdf-python", type=Path, required=True)
    a.add_argument("--require-parameter-complete", action="store_true")
    r = commands.add_parser("reconstruct")
    r.add_argument("--directory", type=Path, required=True)
    r.add_argument("--distribution", type=Path, required=True)
    r.add_argument("--manifest-sha256", required=True)
    args = p.parse_args(argv)
    if args.command == "reconstruct":
        _, result = reconstruct(
            args.directory,
            args.manifest_sha256,
            assets=PCFFPDMSAssets(args.distribution),
        )
        print(
            canonical(
                {
                    "identity": result.payload["identity"],
                    "parameter_completeness": result.payload["parameter_completeness"],
                }
            )
        )
        return 0
    result = audit(
        args.distribution,
        args.archive,
        args.output,
        references=args.references,
        pdf_python=args.pdf_python,
    )
    print(
        canonical(
            {
                k: result[k]
                for k in (
                    "identity",
                    "registered_profiles",
                    "missing_frc_requests",
                    "missing_frc_source_queries",
                )
            }
        )
    )
    return int(args.require_parameter_complete)


if __name__ == "__main__":
    raise SystemExit(main())
