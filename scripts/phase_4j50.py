"""Offline four-layer PDMS reconciliation. Exit 1 until source MD gates close."""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from island.forcefields.opls_registry import (
    OPLSAssets,
    OPLSResolutionError,
    use_opls_assets,
)
from island.forcefields.opls_registry.registry import PAPER_SHA256
from island.opls_conventions import LAYERS, MODES, POLICIES, signed
from island.opls_execution import (
    LAMMPSOptions,
    generate_lammps,
    prepare_polymer_forcefield,
    reconstruct,
    run_md,
    run_static,
)
from island.pdms_lammps import pin, seal
from scripts.phase_4j47 import Reference, write
from scripts.phase_4j49 import analyze_mode
from scripts.qualify_pdms_public import molecule

BASE = "b207a3741f0b6c2e1a6c632a283d82fefdd5e0d9"


def fixture(dp):
    system, _ = molecule(dp)
    xyz = np.array([system.coordinates.get(i) for i in system.topology.sites])
    system.coordinates.translate(20 - xyz.mean(axis=0))
    return system


def relocated(root, result, assets, expected_identity):
    destination = root.parent / (root.name + "-moved")
    shutil.copytree(root, destination)
    archive = destination.parent / "relocated-source.zip"
    if not archive.exists():
        shutil.copyfile(assets.archive_path, archive)
    records = []
    for seed in (1, 99397):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "reconstruct",
            "--directory",
            str(destination.resolve()),
            "--manifest-sha256",
            result["manifest_sha256"],
            "--archive",
            str(archive.resolve()),
            "--archive-sha256",
            assets.archive_sha256,
        ]
        proc = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
        )
        if (
            proc.returncode != 0
            or json.loads(proc.stdout)["identity"] != expected_identity
        ):
            raise ValueError("Fresh process relocation failed: " + proc.stderr)
        records.append(
            {
                "hash_seed": seed,
                "identity": expected_identity,
                "returncode": proc.returncode,
            }
        )
    shutil.rmtree(destination)
    archive.unlink()
    return {
        "verified": True,
        "fresh_processes": records,
        "source_archive_relocated": True,
        "independent_manifest_sha256": result["manifest_sha256"],
    }


def blocked_modes(system, assets, directory):
    records = []
    for mode in MODES:
        with use_opls_assets(assets):
            prepared = prepare_polymer_forcefield(
                system, forcefield="opls", execution_mode=mode
            )
        if mode == "research_opt_in":
            write(directory / "research-conventions.json", prepared.payload)
        for opt_in in (False, True):
            target = directory / f"forbidden-{mode}-{opt_in}"
            try:
                run_md(system, prepared, target, research_opt_in=opt_in)
            except OPLSResolutionError as error:
                assert not target.exists()
                records.append(
                    {
                        "execution_mode": mode,
                        "caller_opt_in": opt_in,
                        "status": error.status,
                        "diagnostics": error.diagnostics,
                        "engine_started": False,
                        "checkpoint_created": False,
                        "trajectory_created": False,
                        "MD_executed": False,
                    }
                )
            else:
                raise ValueError("Unqualified MD was authorized")
    return records


def audit(args):
    if pin(args.paper) != PAPER_SHA256:
        raise ValueError("Pinned local paper hash changed")
    assets = OPLSAssets(args.archive, args.archive_sha256)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    engines = {
        POLICIES[0]: LAMMPSOptions(args.lammps, args.lammps_sha256),
        POLICIES[1]: LAMMPSOptions(args.exact_kernel, args.exact_kernel_sha256),
    }
    cases = {}
    for dp in (0, 2, 4, 6):
        system = fixture(dp)
        case_dir = root / f"dp{dp}"
        case_dir.mkdir()
        comparisons = {}
        for policy in POLICIES:
            with use_opls_assets(assets):
                prepared = prepare_polymer_forcefield(
                    system, forcefield="opls", backend_policy=policy
                )
            run_dir = case_dir / policy
            result = run_static(system, prepared, run_dir, engine=engines[policy])
            # Independent repeat export must have exactly identical bytes.
            repeat = case_dir / (policy + "-repeat")
            generated = generate_lammps(system, prepared, repeat)
            assert generated == result["deterministic_inputs"]
            shutil.rmtree(repeat)
            source = prepared.assignment.validate_integrity(system).validate_integrity(
                system
            )
            xyz = np.array(
                [
                    system.coordinates.get(r["site_id"])
                    for r in source.payload["mapping"]
                ]
            )
            probes = [
                (
                    next(
                        n
                        for n, a in enumerate(source.payload["mapping"])
                        if a["element"] == element
                    ),
                    axis,
                )
                for element in ("H", "O", "Si", "C")
                for axis in range(3)
            ]
            prefix = (
                (run_dir / "lammps/system.in")
                .read_text()
                .replace("read_data system.data", "read_data ../system.data")
                .replace("include coefficients.in", "include ../coefficients.in")
            )
            reference = Reference(source.payload, xyz)
            try:
                comparison = analyze_mode(
                    run_dir / "lammps/comparison",
                    prefix,
                    xyz,
                    probes,
                    engines[policy],
                    reference,
                    source,
                )
            finally:
                reference.close()
            write(run_dir / "components.json", comparison)
            # Added raw FD records become part of the checksummed workflow.
            result["manifest_sha256"] = seal(run_dir)
            relocation = relocated(run_dir, result, assets, prepared.identity)
            comparisons[policy] = {
                "convention_identity": prepared.identity,
                "layer_identities": {
                    k: prepared.payload[k]["identity"] for k in LAYERS
                },
                "canonical_assignment_identity": prepared.assignment.identity,
                "energy_error_kj_mol": comparison["energy_error_kj_mol"],
                "force_error_kj_mol_angstrom": comparison[
                    "force_error_kj_mol_angstrom"
                ],
                "finite_difference_passed": comparison["finite_difference"]["passed"],
                "deterministic_export_verified": True,
                "relocation": relocation,
                "MD_executed": False,
            }
            if dp == 0:
                write(
                    root / (policy + "-convention-table.json"),
                    signed(
                        {
                            "rows": [
                                {
                                    "setting": k,
                                    "source": prepared.payload[LAYERS[0]]["settings"][
                                        k
                                    ],
                                    "backend": prepared.payload[LAYERS[2]]["settings"][
                                        k
                                    ],
                                    "converter_record_identity": prepared.payload[
                                        LAYERS[1]
                                    ]["identity"],
                                    "scientific_record_identity": prepared.payload[
                                        LAYERS[3]
                                    ]["identity"],
                                }
                                for k in prepared.payload["required_settings"]
                            ]
                        }
                    ),
                )
        cases[str(dp)] = {
            "policies": comparisons,
            "blocked_modes": blocked_modes(system, assets, case_dir),
        }
        author, exact = (comparisons[k] for k in POLICIES)
        reconciliation = {
            "source_setting": "PPPM 1e-4; source does not specify pair table/kernel/precision",
            "current_backend_setting": "explicit exact-erfc kernel; Ewald 1e-7; table off; matched alpha; kmax 20",
            "physical_parameters_or_charges_changed": False,
            "difference_classification": "representational numerical solver/kernel; source scope and pressure remain unresolved",
            "affected_components": [
                "Coulomb real/reciprocal",
                "total energy",
                "total/Coulomb forces",
                "FD",
            ],
            "source_replay_errors": author["energy_error_kj_mol"],
            "backend_diagnostic_errors": exact["energy_error_kj_mol"],
            "source_MD_authorized": False,
        }
        write(case_dir / "reconciliation.json", signed(reconciliation))
        print(
            f"DP{dp}: diagnostic comparisons recorded; research MD blocked", flush=True
        )
    policies = {
        "diagnostic": {
            "permits": ["static", "component_comparison"],
            "MD": False,
            "requirements": [
                "all selected-backend settings resolved and cited",
                "four-layer integrity",
            ],
        },
        "research_opt_in": {
            "permits_if_qualified": ["short_minimization", "short_NVT"],
            "available": False,
            "requires_explicit_opt_in": True,
            "requirements": [
                "all source settings resolved",
                "all backend settings resolved",
                "independent numerical qualification",
                "author-convention reconciliation",
            ],
        },
        "production_candidate": {
            "available": False,
            "requirements": [
                "Ewald",
                "pressure",
                "cutoff",
                "sampling",
                "cross-engine gates",
            ],
        },
    }
    write(root / "execution-policy.json", signed(policies))
    exact_pass = all(
        all(
            v <= 3e-4
            for v in c["policies"][POLICIES[1]]["energy_error_kj_mol"].values()
        )
        and all(
            v <= 3e-4
            for v in c["policies"][POLICIES[1]]["force_error_kj_mol_angstrom"].values()
        )
        and c["policies"][POLICIES[1]]["finite_difference_passed"]
        for c in cases.values()
    )
    summary = signed(
        {
            "schema": "island.j50-reconciliation-summary.v1",
            "base_commit": BASE,
            "cases": cases,
            "source_assets": {"archive": args.archive_sha256, "paper": PAPER_SHA256},
            "diagnostic_acceptance_gate": 0 if exact_pass else 1,
            "research_opt_in_acceptance_gate": 1,
            "production_candidate_acceptance_gate": 1,
            "source_authority_promoted": False,
            "new_physical_profiles": [],
            "MD_executed": False,
            "NPT_executed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "summary.json", summary)
    return 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    p = subs.add_parser("audit")
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--paper", type=Path, required=True)
    p.add_argument("--lammps", type=Path, required=True)
    p.add_argument("--lammps-sha256", required=True)
    p.add_argument("--exact-kernel", type=Path, required=True)
    p.add_argument("--exact-kernel-sha256", required=True)
    p.add_argument("--output", type=Path, required=True)
    p = subs.add_parser("reconstruct")
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args(argv)
    if args.command == "audit":
        return audit(args)
    _, preparation = reconstruct(
        args.directory,
        manifest_sha256=args.manifest_sha256,
        assets=OPLSAssets(args.archive, args.archive_sha256),
    )
    print(
        json.dumps(
            {
                "identity": preparation.identity,
                "four_layer_identities": {
                    k: v.identity for k, v in preparation.records().items()
                },
                "MD_authorized": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
