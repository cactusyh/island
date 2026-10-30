"""J53: detect genuinely new assets without repeating the exhausted J52 audit."""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from island.forcefields.opls_registry import OPLSAssets
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import canonical
from island.opls_conventions import signed
from island.opls_final_resolution import ConventionSources, use_convention_sources
from island.opls_source_execution import reconstruct, run_md, save_preparation
from island.opls_source_intake import (
    blocker_reopen_report,
    inspect_source,
    inventory,
    prepare_polymer_forcefield,
)
from island.opls_source_intake.schema import source_intake_schema
from island.polymer import use_opls_assets
from scripts.phase_4j47 import write
from scripts.phase_4j50 import fixture
from scripts.validate_phase_4j47_receipt import pin

BASE = "c5a06496a9c27560ef627baf21bb95aeaef40080"
BRANCH = "codex/phase-4j53-authoritative-opls-pdms-source-intake"


def asset_delta(repo, sources):
    """Hashes only: previously recorded assets are not read as new authority."""
    repo, sources = Path(repo), Path(sources)
    names = (
        subprocess.check_output(
            ["git", "ls-tree", "-r", "--name-only", BASE, "docs"],
            cwd=repo,
        )
        .decode()
        .splitlines()
    )
    known = {}
    for name in names:
        if Path(name).suffix not in (".json", ".md"):
            continue
        text = (repo / name).read_text()
        for sha in set(re.findall(r"\b[0-9a-f]{64}\b", text)):
            known.setdefault(sha, []).append(name)
    rows = inventory(sources)
    return signed(
        {
            "schema": "island.j53-local-source-delta.v1",
            "base_commit": BASE,
            "method": "file hashes compared to committed evidence; no repeated convention audit",
            "files": [
                {**r, "previous_evidence": sorted(known.get(r["sha256"], []))}
                for r in rows
            ],
            "new_content_candidates": [r for r in rows if r["sha256"] not in known],
            "explicit_new_source_package": None,
            "no_network_discovery": True,
        }
    )


def audit(args):
    root = args.output
    root.mkdir(parents=True, exist_ok=False)
    delta = asset_delta(args.repo, args.sources)
    write(root / "source-delta.json", delta)
    write(root / "source-intake-schema.json", source_intake_schema())
    result = inspect_source()
    write(root / "source-intake.json", result)
    write(root / "blocker-reopen.json", blocker_reopen_report(result))
    previous = {}
    for phase in range(49, 53):
        for path in sorted(
            (args.repo / f"docs/evidence/phase_4j{phase}").rglob("*.json")
        ):
            if any(
                s in path.name
                for s in (
                    "convention",
                    "component",
                    "resolution",
                    "blocker",
                    "comparison",
                )
            ):
                json.loads(path.read_text())
                previous[str(path.relative_to(args.repo))] = pin(path)
    write(
        root / "prior-comparison-records.json",
        signed(
            {
                "schema": "island.j53-prior-records.v1",
                "files": previous,
                "new_comparisons_run": False,
                "source_model_changed": False,
            }
        ),
    )
    shutil.copyfile(args.archive, root.parent / "relocated-author.zip")
    shutil.copyfile(args.paper, root.parent / "relocated-paper.pdf")
    shutil.copyfile(args.supplement, root.parent / "relocated-supplement.pdf")
    records = []
    assets = OPLSAssets(args.archive, pin(args.archive))
    sources = ConventionSources(args.paper, args.supplement)
    for dp in (0, 2, 4, 6):
        system = fixture(dp)
        with use_opls_assets(assets), use_convention_sources(sources):
            prepared = prepare_polymer_forcefield(system, forcefield="opls")
        directory = root / f"dp{dp}"
        manifest_pin = save_preparation(system, prepared, directory)
        blocked = []
        for action in ("minimize", "nvt", "npt"):
            from island.forcefields.opls_registry import OPLSResolutionError

            output = root / f"forbidden-{dp}-{action}"
            try:
                run_md(system, prepared, output, action=action, research_opt_in=True)
                raise RuntimeError("Source-blocked action was executed")
            except OPLSResolutionError as error:
                blocked.append(
                    {
                        "action": action,
                        "status": error.status,
                        "stage": error.diagnostics[0]["stage"],
                        "outputs_created": output.exists(),
                        "engine_started": False,
                    }
                )
        moved = root.parent / f"relocated-dp{dp}"
        shutil.copytree(directory, moved)
        fresh = []
        for seed in (1, 99397):
            command = [
                sys.executable,
                "-m",
                "scripts.phase_4j53",
                "reconstruct",
                "--directory",
                str(moved),
                "--manifest-sha256",
                manifest_pin,
                "--archive",
                str(root.parent / "relocated-author.zip"),
                "--paper",
                str(root.parent / "relocated-paper.pdf"),
                "--supplement",
                str(root.parent / "relocated-supplement.pdf"),
            ]
            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
                env={**os.environ, "PYTHONHASHSEED": str(seed)},
            )
            if (
                process.returncode
                or json.loads(process.stdout)["identity"] != prepared.identity
            ):
                raise ValueError(
                    "Fresh process intake relocation failed: " + process.stderr
                )
            fresh.append(
                {
                    "hash_seed": seed,
                    "command": command,
                    "identity": prepared.identity,
                    "returncode": process.returncode,
                    "stdout": process.stdout,
                    "stderr": process.stderr,
                }
            )
        records.append(
            {
                "dp": dp,
                "preparation_identity": prepared.identity,
                "bundle_manifest_sha256": manifest_pin,
                "selected_profile": prepared.payload["selected_profile"],
                "source_status": result["status"],
                "blocked_actions": blocked,
                "relocated_fresh_processes": fresh,
            }
        )
    summary = signed(
        {
            "schema": "island.j53-audit.v1",
            "base_commit": BASE,
            "family": "opls",
            "selected_profile": PROFILE,
            "source_status": result["status"],
            "new_source_packages": [],
            "new_profiles": [],
            "new_content_candidates": delta["new_content_candidates"],
            "cases": records,
            "MD_executed": False,
            "NPT_executed": False,
            "new_static_comparisons_executed": False,
            "checkpoint_RNG_continuation": "not run: no new authoritative source",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "summary.json", summary)
    print(canonical(summary))
    return 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    subs = p.add_subparsers(dest="command", required=True)
    a = subs.add_parser("audit")
    a.add_argument("--repo", type=Path, default=Path.cwd())
    a.add_argument("--sources", type=Path, required=True)
    a.add_argument("--output", type=Path, required=True)
    r = subs.add_parser("reconstruct")
    r.add_argument("--directory", type=Path, required=True)
    r.add_argument("--manifest-sha256", required=True)
    for sub in (a, r):
        for arg in ("archive", "paper", "supplement"):
            sub.add_argument("--" + arg, type=Path, required=True)
    args = p.parse_args(argv)
    if args.command == "audit":
        return audit(args)
    with (
        use_opls_assets(OPLSAssets(args.archive, pin(args.archive))),
        use_convention_sources(ConventionSources(args.paper, args.supplement)),
    ):
        _, prepared = reconstruct(args.directory, manifest_sha256=args.manifest_sha256)
    print(canonical({"identity": prepared.identity, "MD_authorized": False}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
