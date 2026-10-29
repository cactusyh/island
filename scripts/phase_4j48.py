"""Reproduce automatic OPLS assignment and raw, fail-closed PDMS statics.

Exit 1 is the numerical acceptance result, not a parser exception. Asset setup
is explicit and no physical parameter, charge or acceptance threshold changes.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from island.forcefields.external_pdms import LIMITS
from island.forcefields.opls_registry import (
    OPLSAssets,
    OPLSResolutionError,
    canonical_registry,
    prepare_polymer_forcefield,
    use_opls_assets,
)
from island.forcefields.opls_registry.registry import PAPER_SHA256
from island.forcefields.pdms_opls import PROFILE, PDMSOPLSError
from island.forcefields.selection.inventory import digest
from island.opls_lammps import LAMMPSOptions, run_md, run_static
from island.pdms_lammps import pin
from scripts.phase_4j47 import compare, lammps_finite_difference, write
from scripts.qualify_pdms_public import molecule

BASE = "fdf117aadf613912786fdfbebe653ad58695bea3"


def audit(args):
    from numpy import asarray

    if pin(args.paper) != PAPER_SHA256:
        raise OPLSResolutionError(
            "source_incompatible", [{"stage": "paper_hash_mismatch"}]
        )
    args.output.mkdir(parents=True, exist_ok=False)
    assets = OPLSAssets(args.archive, args.archive_sha256)
    engine = LAMMPSOptions(args.lammps, args.lammps_sha256)
    cases = {}
    write(args.output / "canonical-registry.json", canonical_registry())
    for dp in (0, 2, 4, 6):
        system, _ = molecule(dp)
        xyz = asarray([system.coordinates.get(i) for i in system.topology.sites])
        system.coordinates.translate(20 - xyz.mean(axis=0))
        with use_opls_assets(assets):
            prepared = prepare_polymer_forcefield(system, forcefield="opls")
        backend = prepared.validate_integrity(system)
        directory = args.output / f"dp{dp}"
        observed = run_static(system, prepared, directory, engine=engine)
        if observed["status"] != "static_executed_unqualified":
            raise OPLSResolutionError(observed["status"], observed["diagnostics"])
        case = compare(system, backend, observed["observed"])
        xyz = asarray(observed["observed"]["coordinates_angstrom"])
        case["LAMMPS_finite_difference"] = lammps_finite_difference(
            directory / "finite-difference-control",
            backend.validate_integrity(system),
            engine,
            xyz,
            observed["observed"]["forces_kj_mol_angstrom"],
        )
        shutil.move(
            str(directory / "finite-difference-control"),
            str(args.output / f"finite-difference-dp{dp}"),
        )
        moved = args.output / f"relocated-dp{dp}"
        shutil.copytree(directory, moved)
        relocated_archive = args.output / "relocated-source.zip"
        shutil.copyfile(args.archive, relocated_archive)
        code = """import sys
from island.forcefields.opls_registry import OPLSAssets
from island.opls_lammps import reconstruct,LAMMPSOptions
s,p,r=reconstruct(sys.argv[1],manifest_sha256=sys.argv[2],assets=OPLSAssets(sys.argv[3],sys.argv[4]),engine=LAMMPSOptions(sys.argv[5],sys.argv[6]))
print(p.identity)
"""
        identities = []
        for seed in ("1", "99397"):
            proc = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    code,
                    str(moved.resolve()),
                    observed["manifest_sha256"],
                    str(relocated_archive.resolve()),
                    args.archive_sha256,
                    str(args.lammps.resolve()),
                    args.lammps_sha256,
                ],
                capture_output=True,
                text=True,
                check=True,
                env={**os.environ, "PYTHONHASHSEED": seed},
            )
            identities.append(proc.stdout.strip())
        assert identities == [prepared.identity] * 2
        shutil.rmtree(moved)
        relocated_archive.unlink()
        # Prove the numerical gate rejects before engine invocation or MD artifacts.
        try:
            run_md(system, prepared, args.output / f"forbidden-MD-dp{dp}")
        except OPLSResolutionError as error:
            assert error.status == "scientifically_unqualified"
            md_diagnostic = error.diagnostics
        else:
            raise AssertionError("Blocked static profile reached MD")
        assert not (args.output / f"forbidden-MD-dp{dp}").exists()
        case.update(
            {
                "public_forcefield": "opls",
                "user_profile_argument": None,
                "user_charge_method_argument": None,
                "selected_profile": prepared.payload["selected_profile"],
                "charge_origin": prepared.payload["charge_origin"],
                "assignment_identity": prepared.identity,
                "graph_identity": prepared.payload["graph_identity"],
                "atoms": len(backend.validate_integrity(system).payload["mapping"]),
                "workflow_manifest_sha256": observed["manifest_sha256"],
                "source_rows_reconstructed": 17,
                "relocation_verified": True,
                "fresh_process_reconstruction_verified": True,
                "process_hash_seeds": [1, 99397],
                "MD_block_diagnostics": md_diagnostic,
                "missing_conventions": prepared.payload["missing_conventions"],
                "minimization": "blocked_by_failed_static_gates",
                "NVT": "blocked_by_failed_static_gates",
                "split_uninterrupted_checkpoint_RNG_equality": "blocked_by_failed_static_gates",
            }
        )
        cases[str(dp)] = case
        write(args.output / "cases.json", cases)
        print(
            f"DP{dp}: automatic OPLS; static delta={case['component_comparison']['total']['absolute_error']:.10g}; MD blocked",
            flush=True,
        )
    summary = {
        "schema": "island.j48-static-acceptance.v1",
        "base_commit": BASE,
        "registry_identity": canonical_registry()["identity"],
        "profile_identity": PROFILE,
        "cases": cases,
        "source_assets": {"archive": args.archive_sha256, "paper": PAPER_SHA256},
        "software_selection_acceptance_gate": 0,
        "lammps_end_to_end_acceptance_gate": 1,
        "scientific_acceptance_gate": 1,
        "new_physical_parameters": [],
        "MD_executed": False,
        "NPT_executed": False,
        "limitations": LIMITS,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "discrepancy_diagnosis": {
            "bond_angle_torsion_LJ": "precision; all within inherited bounds",
            "electrostatics": "long-range solver; author PPPM 1e-4/default pair table versus Ewald 1e-7",
            "historical_exact_control": "erfc polynomial precision and Coulomb unit constant; raw failures remain",
            "matching_reciprocal_contract": "absent; no authorized replacement numerical convention supplied",
            "virial_pressure": "unresolved; no inherited virial bound; pressure unavailable",
            "forbidden_adjustments": [],
        },
    }
    write(args.output / "summary.json", {**summary, "identity": digest(summary)})
    return 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    parser.add_argument("--lammps-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        return audit(args)
    except (OPLSResolutionError, PDMSOPLSError, OSError) as error:
        print(
            json.dumps(
                {
                    "status": getattr(error, "status", "asset_unavailable"),
                    "diagnostics": getattr(
                        error,
                        "diagnostics",
                        [{"stage": "explicit_asset_unavailable", "reason": str(error)}],
                    ),
                    "lammps_end_to_end_acceptance_gate": 1,
                    "MD_executed": False,
                },
                sort_keys=True,
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
