"""Offline, explicit ten-convention closure audit. Exit 1 while source gaps remain."""

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
from island.opls_resolution import ConventionSources, prepare_polymer_forcefield
from island.opls_resolved_execution import (
    LAMMPSOptions,
    generate_lammps,
    reconstruct,
    run_md,
    run_static,
)
from island.pdms_lammps import seal
from scripts.phase_4j47 import Reference, special_energies, write
from scripts.phase_4j49 import analyze_mode
from scripts.phase_4j50 import fixture
from scripts.qualify_pdms_public import LMP_COULOMB, MM_COULOMB

BASE = "72dbb07a0d61a15436a04720357d76f63e41d927"
BRANCH = "codex/phase-4j51-opls-pdms-source-convention-closure"


def relocated(directory, manifest, assets, sources, identity):
    moved = directory.parent / (directory.name + "-moved")
    shutil.copytree(directory, moved)
    archive, paper, si = (moved / n for n in ("relocated.zip", "paper.pdf", "si.pdf"))
    # Relocated sources live outside the workflow membership, never in its manifest.
    for p, original in (
        (archive, assets.archive_path),
        (paper, sources.paper),
        (si, sources.supplement),
    ):
        shutil.copyfile(original, moved.parent / p.name)
    archive, paper, si = (moved.parent / p.name for p in (archive, paper, si))
    processes = []
    try:
        for seed in (1, 99397):
            proc = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "reconstruct",
                    "--directory",
                    str(moved),
                    "--manifest-sha256",
                    manifest,
                    "--archive",
                    str(archive),
                    "--archive-sha256",
                    assets.archive_sha256,
                    "--paper",
                    str(paper),
                    "--supplement",
                    str(si),
                ],
                capture_output=True,
                text=True,
                check=False,
                env={**os.environ, "PYTHONHASHSEED": str(seed)},
            )
            if proc.returncode or json.loads(proc.stdout)["identity"] != identity:
                raise ValueError("Fresh relocation failed: " + proc.stderr)
            processes.append({"hash_seed": seed, "returncode": 0, "identity": identity})
    finally:
        shutil.rmtree(moved)
        for p in (archive, paper, si):
            p.unlink()
    return {
        "verified": True,
        "fresh_processes": processes,
        "ZIP_and_both_PDFs_relocated": True,
        "manifest_sha256": manifest,
    }


def blocked(system, original, directory):
    outcomes = []
    for mode in MODES:
        p = prepare_polymer_forcefield(
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
                        p,
                        target,
                        engine=object(),
                        action=action,
                        research_opt_in=opt_in,
                    )
                except OPLSResolutionError as error:
                    assert not target.exists()
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
                    raise ValueError("Blocked convention reached execution")
    return outcomes


def detail(ref, source, xyz, comparison):
    mm = ref.components()
    lmp = comparison["LAMMPS"]["energy_kj_mol"]
    special = special_energies(source.payload, xyz)
    affine = ref.virial()
    virial = comparison["LAMMPS"]["virial_kj_mol_xx_yy_zz_xy_xz_yz"]
    return {
        "OpenMM": mm,
        "LAMMPS": lmp,
        "special_pairs_independently_reconstructed": special,
        "LAMMPS_coulomb14_physical_equation_kj_mol": special[
            "coulomb14_exact_equation_kj_mol"
        ]
        * LMP_COULOMB
        / MM_COULOMB,
        "special_pair_partition": "independent physical partition, not a native LAMMPS observable; no output corrections",
        "OpenMM_affine_virial": affine,
        "LAMMPS_virial": virial,
        "virial_absolute_error": [
            abs(a["values_at_1e-4_and_5e-5"][-1] - b) for a, b in zip(affine, virial)
        ],
        "virial_acceptance_bound": None,
        "pressure_validated": False,
        "pressure_proxy": comparison["LAMMPS"].get(
            "configurational_pressure_proxy_atm"
        ),
        "component_classification": {
            "bond_angle_torsion_LJ": "precision issue",
            "improper": "not applicable: no declared-domain improper interaction",
            "coulomb_real_reciprocal": "OpenMM/LAMMPS backend mismatch: PPPM vs Ewald, tables and erfc kernel; see selected policy",
            "coulomb_constant": "precision issue; engine-owned constants retained",
            "special_pairs": "same source mapping; explicit equation reconstructed, no native split-force observable",
            "virial_pressure_proxy": "unresolved scientific discrepancy; no inherited acceptance bound",
            "SI_vs_archive_pair_style": "source convention scope mismatch; no DFT/NPA branch mapping supplied",
        },
        "identical_physical_inputs": [
            "graph",
            "coordinates",
            "atom order",
            "box",
            "masses",
            "charges",
            "exclusions",
            "1-4 scaling",
            "LJ mixing",
            "cutoffs",
            "CHARMM LJ switching",
            "unit conversion",
            "double precision",
        ],
        "numerical_contract": {
            "OpenMM": "Reference Ewald 1e-7, exact erfc, engine-owned Coulomb constant",
            "LAMMPS": comparison["numerical_input_addendum"],
            "source_policy_vs_OpenMM_identical_solver": False,
            "diagnostic_control_not_source_replacement": True,
        },
    }


def audit(args):
    assets = OPLSAssets(args.archive, args.archive_sha256)
    sources = ConventionSources(args.paper, args.supplement)
    sources.validate()
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
        policies = {}
        for policy in POLICIES:
            p = prepare_polymer_forcefield(
                system,
                forcefield="opls",
                assets=assets,
                sources=sources,
                backend_policy=policy,
            )
            directory = case_dir / policy
            result = run_static(system, p, directory, engine=engines[policy])
            repeat = case_dir / (policy + "-repeat")
            assert generate_lammps(system, p, repeat) == result["deterministic_inputs"]
            shutil.rmtree(repeat)
            source = p.assignment.validate_integrity(system).validate_integrity(system)
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
                (directory / "lammps/system.in")
                .read_text()
                .replace("read_data system.data", "read_data ../system.data")
                .replace("include coefficients.in", "include ../coefficients.in")
            )
            ref = Reference(source.payload, xyz)
            try:
                c = analyze_mode(
                    directory / "lammps/comparison",
                    prefix,
                    xyz,
                    probes,
                    engines[policy],
                    ref,
                    source,
                )
                write(
                    directory / "detailed-components.json", detail(ref, source, xyz, c)
                )
            finally:
                ref.close()
            write(directory / "components.json", c)
            manifest = seal(directory)
            policies[policy] = {
                "resolution_identity": p.identity,
                "assignment_identity": p.assignment.identity,
                "deterministic_export": True,
                "energy_error_kj_mol": c["energy_error_kj_mol"],
                "force_error_kj_mol_angstrom": c["force_error_kj_mol_angstrom"],
                "FD_passed": c["finite_difference"]["passed"],
                "relocation": relocated(
                    directory, manifest, assets, sources, p.identity
                ),
            }
            if dp == 0:
                write(
                    root / (policy + "-ten-conventions.json"),
                    signed(
                        {
                            "rows": p.payload["ten_convention_resolutions"],
                            "all_conventions": p.payload["rows"],
                            "source_scope_review": p.payload["source_scope_review"],
                        }
                    ),
                )
        cases[str(dp)] = {
            "policies": policies,
            "blocked_execution": blocked(system, p, case_dir),
        }
        print(
            f"DP{dp}: eight-way source/backend conventions retained; research MD blocked",
            flush=True,
        )
    summary = signed(
        {
            "schema": "island.j51-summary.v1",
            "base_commit": BASE,
            "cases": cases,
            "source_assets": sources.validate(),
            "source_archive_sha256": assets.archive_sha256,
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
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "summary.json", summary)
    return 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    for command in ("audit", "reconstruct"):
        p = subs.add_parser(command)
        p.add_argument("--archive", type=Path, required=True)
        p.add_argument("--archive-sha256", required=True)
        p.add_argument("--paper", type=Path, required=True)
        p.add_argument("--supplement", type=Path, required=True)
        if command == "audit":
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
    a = parser.parse_args(argv)
    if a.command == "audit":
        return audit(a)
    _, p = reconstruct(
        a.directory,
        manifest_sha256=a.manifest_sha256,
        assets=OPLSAssets(a.archive, a.archive_sha256),
        sources=ConventionSources(a.paper, a.supplement),
    )
    print(
        json.dumps(
            {"identity": p.identity, "research_MD_available": False}, sort_keys=True
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
