"""Assign the explicit local author DFT/NPA PDMS source and run research NVT."""

import argparse
import json
from pathlib import Path

import numpy as np

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.core import SimulationBox
from island.forcefields.pdms_opls import export_lammps, save_preparation
from island.forcefields.research import (
    ForceFieldRequest,
    PDMSOPLSError,
    PDMSOPLSOptions,
    prepare_forcefield,
)
from island.pdms import resume_pdms_md, run_pdms_md


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    graph = parser.add_mutually_exclusive_group()
    graph.add_argument("--dp", type=int, default=None)
    graph.add_argument("--smiles")
    parser.add_argument("--e-oxygen-site-id", type=int)
    parser.add_argument("--box-angstrom", type=float, default=40)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--temperature", type=float, default=300)
    parser.add_argument("--timestep-fs", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--report-interval", type=int, default=100)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output must be a new directory")
    if args.resume:
        if not args.manifest_sha256:
            parser.error("Resume requires an independent --manifest-sha256")
        from island.charge_references.records import unpack

        try:
            previous = unpack((args.resume / "preparation.json").read_text())
            options = PDMSOPLSOptions(
                args.archive, args.archive_sha256, tuple(previous["e_oxygen_site_ids"])
            )
            report = resume_pdms_md(
                args.resume,
                args.output,
                options=options,
                manifest_sha256=args.manifest_sha256,
                steps=args.steps,
                report_interval=args.report_interval,
            )
        except PDMSOPLSError as error:
            print(
                json.dumps({"status": error.status, "diagnostics": error.diagnostics})
            )
            return 1
        except (OSError, ValueError, KeyError, TypeError) as error:
            print(
                json.dumps(
                    {
                        "status": "asset_unavailable",
                        "diagnostics": [
                            {
                                "stage": "resume_record_unavailable_or_invalid",
                                "reason": str(error),
                            }
                        ],
                    }
                )
            )
            return 1
        print(json.dumps(report))
        return 0
    if args.smiles:
        system = from_smiles(args.smiles, random_seed=args.seed)
        anchors = () if args.e_oxygen_site_id is None else (args.e_oxygen_site_id,)
    else:
        system = build_linear_polymer(
            "[*:1]O[Si](C)(C)[*:2]",
            dp=args.dp if args.dp is not None else 4,
            head_end_group="[*][Si](C)(C)C",
            tail_end_group="[*]C",
            add_hydrogens=True,
            random_seed=args.seed,
        )
        anchors = (
            system.metadata["polymer"]["head_site_id"]
            if args.e_oxygen_site_id is None
            else args.e_oxygen_site_id,
        )
    system.box = SimulationBox(args.box_angstrom, args.box_angstrom, args.box_angstrom)
    xyz = np.array([system.coordinates.get(i) for i in system.topology.sites])
    # Explicit geometric preflight: avoid self-image overlap in this single-chain demo.
    if max(np.ptp(xyz, axis=0)) + 22 >= args.box_angstrom:
        parser.error(
            "Box must exceed the molecular span plus two 11 A cutoffs; increase --box-angstrom"
        )
    system.coordinates.translate(args.box_angstrom / 2 - np.mean(xyz, axis=0))
    options = PDMSOPLSOptions(args.archive, args.archive_sha256, anchors)
    prepared = prepare_forcefield(
        system, ForceFieldRequest("external_pdms_opls", options)
    )
    try:
        prepared.require_executable(system)
        if args.prepare_only:
            args.output.mkdir(parents=True)
            pin = save_preparation(prepared, args.output / "preparation.json")
            files = export_lammps(system, prepared, args.output / "lammps")
            report = {
                "status": "prepared",
                "preparation_sha256": pin,
                "profile_identity": prepared.payload["profile_identity"],
                "lammps_files": files,
                "MD_run": False,
            }
            (args.output / "result.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
        else:
            report = run_pdms_md(
                system,
                prepared,
                args.output,
                steps=args.steps,
                temperature_kelvin=args.temperature,
                timestep_fs=args.timestep_fs,
                seed=args.seed,
                report_interval=args.report_interval,
            )
    except PDMSOPLSError as error:
        report = {
            "status": error.status,
            "diagnostics": error.diagnostics,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "blocked.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, ensure_ascii=False))
        return 1
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
