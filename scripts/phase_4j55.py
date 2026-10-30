"""Reproduce independent J55 policy reports without numerical/MD promotion."""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from island.forcefield_policy import (
    OPLSAssets,
    check_lammps,
    check_openmm,
    export_lammps,
    prepare_polymer_forcefield,
)
from island.forcefield_policy.evidence import reconstruct, save
from island.forcefields.opls_registry.registry import ARCHIVE_SHA256
from island.forcefields.pdms_opls import PROFILE
from island.opls_conventions import signed
from island.opls_resolution import ConventionSources
from scripts.phase_4j47 import write
from scripts.phase_4j50 import fixture
from scripts.validate_phase_4j47_receipt import pin

BASE = "6715a9d742b9bf85de37acabb512fbacfa74b9fb"
BRANCH = "codex/phase-4j55-parameter-completeness-policy"


def audit(args):
    from island.chemistry import from_smiles
    from island.pdms_lammps import LAMMPSOptions

    root, repo = args.output, Path(__file__).resolve().parents[1]
    assets = OPLSAssets(args.archive, ARCHIVE_SHA256, args.native, pin(args.native))
    sources = ConventionSources(args.paper, args.supplement)
    cases = []
    for name, system in [("native-ethane", from_smiles("CC", random_seed=2026))] + [
        (f"dp{dp}", fixture(dp)) for dp in (0, 2, 4, 6)
    ]:
        prepared = prepare_polymer_forcefield(
            system, forcefield="opls", assets=assets, sources=sources
        )
        out = root / name
        if name.startswith("dp"):
            prepared = export_lammps(system, prepared, out / "lammps")
            if args.lammps:
                prepared = check_lammps(
                    system,
                    prepared,
                    out / "static",
                    engine=LAMMPSOptions(args.lammps, args.lammps_sha256),
                )
        prepared = check_openmm(system, prepared)
        sha = save(system, prepared, out / "bundle")
        write(out / "report.json", prepared.payload)
        relocation = []
        with tempfile.TemporaryDirectory(prefix="island-j55-relocate-") as tmp:
            moved = Path(tmp)
            shutil.copytree(out / "bundle", moved / "bundle")
            for key, value in [
                ("archive", args.archive),
                ("native", args.native),
                ("paper", args.paper),
                ("supplement", args.supplement),
            ]:
                shutil.copyfile(value, moved / key)
            for seed in ("1", "99397"):
                command = [
                    sys.executable,
                    "-m",
                    "scripts.phase_4j55",
                    "reconstruct",
                    "--bundle",
                    str(moved / "bundle"),
                    "--manifest-sha256",
                    sha,
                    "--archive",
                    str(moved / "archive"),
                    "--native",
                    str(moved / "native"),
                    "--paper",
                    str(moved / "paper"),
                    "--supplement",
                    str(moved / "supplement"),
                    "--output",
                    str(moved / "reconstructed.json"),
                ]
                run = subprocess.run(
                    command,
                    env={**os.environ, "PYTHONHASHSEED": seed},
                    check=False,
                    capture_output=True,
                    text=True,
                )
                if run.returncode:
                    raise RuntimeError(run.stderr)
                if (
                    json.loads((moved / "reconstructed.json").read_text())
                    != prepared.payload
                ):
                    raise AssertionError("Changed relocated report")
                relocation.append(
                    {
                        "seed": seed,
                        "returncode": run.returncode,
                        "report_identity": prepared.identity,
                        "source_files_relocated": True,
                    }
                )
        cases.append(
            {
                "case": name,
                "report_identity": prepared.identity,
                "bundle_manifest_sha256": sha,
                "relocation": relocation,
                "parameter_completeness": prepared.payload["parameter_completeness"],
                "missing_conventions": prepared.payload["missing_conventions"],
                "lammps_export_status": prepared.payload["lammps_export_status"],
                "lammps_runtime_status": prepared.payload["lammps_runtime_status"],
                "openmm_build_status": prepared.payload["openmm_build_status"],
                "openmm_runtime_status": prepared.payload["openmm_runtime_status"],
            }
        )
        print(name, cases[-1], flush=True)
    prior = {
        str(p.relative_to(repo)): pin(p)
        for phase in range(47, 54)
        for p in (repo / f"docs/evidence/phase_4j{phase}").rglob("*.json")
        if any(
            k in p.name
            for k in ("component", "comparison", "convention", "resolution", "blocker")
        )
    }
    write(
        root / "prior-evidence.json",
        signed({"files": prior, "historical_claims_rewritten": False}),
    )
    write(
        root / "summary.json",
        signed(
            {
                "schema": "island.j55-audit.v1",
                "base_commit": BASE,
                "cases": cases,
                "source_hashes": {
                    Path(f).name: pin(f)
                    for f in (args.archive, args.native, args.paper, args.supplement)
                },
                "profile_identity": PROFILE,
                "new_profiles": [],
                "MD_executed": False,
                "cross_engine_diagnostics_requested": False,
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
    )


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["audit", "reconstruct"])
    for name in ("archive", "native", "paper", "supplement"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--lammps", type=Path)
    p.add_argument("--lammps-sha256")
    p.add_argument("--bundle", type=Path)
    p.add_argument("--manifest-sha256")
    args = p.parse_args(argv)
    if args.action == "audit":
        audit(args)
    else:
        assets = OPLSAssets(args.archive, ARCHIVE_SHA256, args.native, pin(args.native))
        _, prep = reconstruct(
            args.bundle,
            manifest_sha256=args.manifest_sha256,
            assets=assets,
            sources=ConventionSources(args.paper, args.supplement),
        )
        write(args.output, prep.payload)


if __name__ == "__main__":
    main()
