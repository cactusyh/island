"""Final nine-convention audit and static replay; exit 1 while MD is blocked."""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from island.forcefields.opls_registry import OPLSAssets, OPLSResolutionError
from island.opls_conventions import MODES, POLICIES, signed
from island.opls_final_execution import (
    LAMMPSOptions,
    generate_lammps,
    prepare_polymer_forcefield,
    reconstruct,
    run_md,
    run_static,
)
from island.opls_final_resolution import ConventionSources, corpus_audit
from island.pdms_lammps import seal
from scripts.phase_4j47 import Reference, write
from scripts.phase_4j49 import analyze_mode
from scripts.phase_4j50 import fixture
from scripts.phase_4j51 import detail

BASE = "78a516cdf9f373b905ad8fd729bc82453ac8a886"
BRANCH = "codex/phase-4j52-opls-pdms-final-convention-resolution"


def relocated(directory, manifest, assets, sources, identity):
    moved = directory.parent / (directory.name + "-moved")
    shutil.copytree(directory, moved)
    files = [moved.parent / n for n in ("relocated.zip", "paper.pdf", "si.pdf")]
    processes = []
    try:
        for target, source in zip(
            files, (assets.archive_path, sources.paper, sources.supplement)
        ):
            shutil.copyfile(source, target)
        for seed in (1, 99397):
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "reconstruct",
                "--directory",
                str(moved),
                "--manifest-sha256",
                manifest,
                "--archive",
                str(files[0]),
                "--archive-sha256",
                assets.archive_sha256,
                "--paper",
                str(files[1]),
                "--supplement",
                str(files[2]),
            ]
            proc = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
                env={**os.environ, "PYTHONHASHSEED": str(seed)},
            )
            if proc.returncode or json.loads(proc.stdout)["identity"] != identity:
                raise ValueError("Fresh J52 relocation failed: " + proc.stderr)
            processes.append({"hash_seed": seed, "returncode": 0, "identity": identity})
    finally:
        shutil.rmtree(moved)
        for p in files:
            p.unlink(missing_ok=True)
    return {
        "verified": True,
        "fresh_processes": processes,
        "ZIP_and_both_PDFs_relocated": True,
        "manifest_sha256": manifest,
    }


def blocked(system, original, directory):
    outcomes = []
    for mode in MODES:
        prepared = prepare_polymer_forcefield(
            system,
            forcefield="opls",
            assets=original.assignment.assets,
            sources=original.sources,
            execution_mode=mode,
        )
        for action in ("minimize", "nvt", "npt"):
            for opt_in in (False, True):
                target = directory / f"blocked-{mode}-{action}-{opt_in}"
                try:
                    run_md(
                        system,
                        prepared,
                        target,
                        engine=object(),
                        action=action,
                        research_opt_in=opt_in,
                    )
                except OPLSResolutionError as error:
                    if target.exists():
                        raise ValueError("Blocked J52 path created output") from error
                    outcomes.append(
                        {
                            "mode": mode,
                            "action": action,
                            "explicit_opt_in": opt_in,
                            "status": error.status,
                            "diagnostics": error.diagnostics,
                            "engine_started": False,
                            "outputs_created": False,
                        }
                    )
                else:
                    raise ValueError("Permanent source blocker reached execution")
    return outcomes


def audit(args):
    assets = OPLSAssets(args.archive, args.archive_sha256)
    sources = ConventionSources(args.paper, args.supplement)
    sources.validate()
    corpus = corpus_audit()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    write(root / "supplied-corpus.json", corpus)
    write(
        root / "permanent-blockers.json",
        signed(
            {
                "closure_scope": corpus["closure_scope"],
                "rows": corpus["final_resolutions"],
                "source_scope_conflict": corpus["source_scope_conflict"],
                "MD_authorized": False,
            }
        ),
    )
    engines = {
        POLICIES[0]: LAMMPSOptions(args.lammps, args.lammps_sha256),
        POLICIES[1]: LAMMPSOptions(args.exact_kernel, args.exact_kernel_sha256),
    }
    cases = {}
    for dp in (0, 2, 4, 6):
        system = fixture(dp)
        case_dir = root / f"dp{dp}"
        case_dir.mkdir()
        policies = {}
        for policy in POLICIES:
            prepared = prepare_polymer_forcefield(
                system,
                forcefield="opls",
                assets=assets,
                sources=sources,
                backend_policy=policy,
            )
            directory = case_dir / policy
            result = run_static(system, prepared, directory, engine=engines[policy])
            repeat = case_dir / (policy + "-repeat")
            if (
                generate_lammps(system, prepared, repeat)
                != result["deterministic_inputs"]
            ):
                raise ValueError("Nondeterministic export")
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
                        if a["element"] == e
                    ),
                    axis,
                )
                for e in ("H", "O", "Si", "C")
                for axis in range(3)
            ]
            prefix = (
                (directory / "lammps/system.in")
                .read_text()
                .replace("read_data system.data", "read_data ../system.data")
                .replace("include coefficients.in", "include ../coefficients.in")
            )
            ref = Reference(source.payload, xyz)
            try:
                comparison = analyze_mode(
                    directory / "lammps/comparison",
                    prefix,
                    xyz,
                    probes,
                    engines[policy],
                    ref,
                    source,
                )
                write(
                    directory / "detailed-components.json",
                    detail(ref, source, xyz, comparison),
                )
            finally:
                ref.close()
            write(directory / "components.json", comparison)
            manifest = seal(directory)
            policies[policy] = {
                "final_resolution_identity": prepared.identity,
                "assignment_identity": prepared.assignment.identity,
                "deterministic_export": True,
                "energy_error_kj_mol": comparison["energy_error_kj_mol"],
                "force_error_kj_mol_angstrom": comparison[
                    "force_error_kj_mol_angstrom"
                ],
                "FD_passed": comparison["finite_difference"]["passed"],
                "relocation": relocated(
                    directory, manifest, assets, sources, prepared.identity
                ),
            }
            if dp == 0:
                write(
                    root / (policy + "-nine-conventions.json"),
                    signed(
                        {
                            "rows": prepared.payload["final_resolutions"],
                            "source_scope_review": prepared.payload[
                                "source_scope_review"
                            ],
                        }
                    ),
                )
        cases[str(dp)] = {
            "policies": policies,
            "blocked_execution": blocked(system, prepared, case_dir),
        }
        print(
            f"DP{dp}: actual statics/FD, relocated fresh processes; research MD blocked",
            flush=True,
        )
    summary = signed(
        {
            "schema": "island.j52-summary.v1",
            "base_commit": BASE,
            "cases": cases,
            "source_assets": sources.validate(),
            "source_archive_sha256": assets.archive_sha256,
            "source_audit_identity": corpus["identity"],
            "diagnostic_acceptance_gate": int(
                not all(
                    all(
                        v <= 3e-4
                        for v in c["policies"][POLICIES[1]][
                            "energy_error_kj_mol"
                        ].values()
                    )
                    and all(
                        v <= 3e-4
                        for v in c["policies"][POLICIES[1]][
                            "force_error_kj_mol_angstrom"
                        ].values()
                    )
                    and c["policies"][POLICIES[1]]["FD_passed"]
                    for c in cases.values()
                )
            ),
            "research_opt_in_acceptance_gate": 1,
            "production_candidate_acceptance_gate": 1,
            "MD_executed": False,
            "NPT_executed": False,
            "source_authority_promoted": False,
            "checkpoint_RNG_continuation": "not run: source prerequisites remain blocked",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "summary.json", summary)
    return 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("audit", "reconstruct"):
        p = sub.add_parser(name)
        p.add_argument("--archive", type=Path, required=True)
        p.add_argument("--archive-sha256", required=True)
        p.add_argument("--paper", type=Path, required=True)
        p.add_argument("--supplement", type=Path, required=True)
        if name == "audit":
            for flag in (
                "lammps",
                "lammps-sha256",
                "exact-kernel",
                "exact-kernel-sha256",
            ):
                p.add_argument("--" + flag, required=True)
            p.add_argument("--output", type=Path, required=True)
        else:
            p.add_argument("--directory", type=Path, required=True)
            p.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args(argv)
    if args.command == "audit":
        return audit(args)
    _, prepared = reconstruct(
        args.directory,
        manifest_sha256=args.manifest_sha256,
        assets=OPLSAssets(args.archive, args.archive_sha256),
        sources=ConventionSources(args.paper, args.supplement),
    )
    print(
        json.dumps(
            {"identity": prepared.identity, "research_MD_available": False},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
