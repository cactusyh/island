"""J93 deterministic PE5 metadata gate and replay receipt."""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

from island.forcefields.assignment_core.schema import is_hash
from island.forcefields.opls_providers.intake93 import SOURCE_DIR, intake
from scripts.phase_4j73_source_intake import (
    exact,
    pin,
    read,
    require,
    signed,
    verify_signed,
    write,
)
from scripts.phase_4j76_source_resolution import scope

BASE = "7c4dbc8070f5d6c76c7e080787a94ab56832ee03"
BRANCH = "codex/phase-4j93-ligpargen-pe5-metadata-gate"
ROOT = Path("docs/evidence/phase_4j93")
REPORTS = ("result.json", "missing-terms-and-blockers.json")
SCOPE = {
    "J93_QM_executions": 0,
    "J93_MD_steps": 0,
    "LAMMPS_commands_executed": 0,
    "energy_force_evaluations": 0,
    "minimization_executions": 0,
    "physical_validations": 0,
    "accepted_coefficients": 0,
    "runtime_registration": False,
    "automatic_assignment_enabled": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "chemical_family_transfer": "not_established",
}


def input_identity(source_dir=SOURCE_DIR, *, data_sha256=None, settings_sha256=None):
    return {
        "source_directory_declared": str(source_dir),
        "data_sha256": data_sha256,
        "settings_sha256": settings_sha256,
    }


def reconstruct(repo=None, **external):
    del repo
    scope()
    result = intake(**external).report
    snapshot = signed(
        {
            "schema": "island.j93.external-input-pins.v1",
            "provided": input_identity(**external),
            "raw_source_records_persisted": False,
        }
    )
    reports = {
        "result.json": signed({**result, "execution_scope": SCOPE}),
        "missing-terms-and-blockers.json": signed(
            {
                "decision": result["decision"],
                "missing_terms": result["missing_terms"],
                "convention_blockers": result["convention_blockers"],
                "blocking_reasons": result["blockers"],
            }
        ),
    }
    return snapshot, reports


def receipt(root, reports):
    return signed(
        {
            "schema": "island.j93.intake-receipt.v1",
            "base": BASE,
            "branch": BRANCH,
            "snapshot_sha256": pin(root / "captured-input.json"),
            "reports": {name: pin(root / name) for name in REPORTS},
            "decision": reports["result.json"]["decision"],
            "scope": SCOPE,
        }
    )


def historical(repo):
    entries = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    count = 0
    for row in entries.splitlines():
        meta, name = row.split("\t", 1)
        path = repo / name
        raw = (
            str(path.readlink()).encode()
            if meta.startswith("120000")
            else path.read_bytes()
        )
        blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        require(blob == meta.split()[2], "historical_byte_mutation:" + name)
        count += 1
    return count


def emit(repo, root, **external):
    scope()
    require(not root.exists(), "do_not_overwrite_intake_attempt")
    snapshot, reports = reconstruct(repo, **external)
    root.mkdir(parents=True)
    write(root / "captured-input.json", snapshot)
    for name, value in reports.items():
        write(root / name, value)
    write(root / "receipt.json", receipt(root, reports))
    write(
        root / "historical-preservation.json",
        signed(
            {
                "base": BASE,
                "files_preserved": historical(Path(repo)),
                "J92_receipt_and_blocker_report": "byte_identical",
            }
        ),
    )
    return reports["result.json"]


def replay(repo, root, snapshot_sha=None, **external):
    scope()
    snapshot = read(root / "captured-input.json")
    require(is_hash(snapshot_sha), "snapshot_hash_required")
    require(pin(root / "captured-input.json") == snapshot_sha, "snapshot_hash_mutation")
    supplied = input_identity(**external)
    for name in ("data_sha256", "settings_sha256"):
        exact(
            snapshot["provided"][name],
            supplied[name],
            "external_source_pin_mismatch:" + name,
        )
    fresh, reports = reconstruct(repo, **external)
    fresh_body = {key: value for key, value in fresh.items() if key != "identity"}
    fresh_body["provided"]["source_directory_declared"] = snapshot["provided"][
        "source_directory_declared"
    ]
    fresh = signed(fresh_body)
    exact(snapshot, fresh, "external_sources_unavailable_or_mutated")
    for name, value in reports.items():
        exact(read(root / name), value, "report_semantic_mutation:" + name)
    stored = read(root / "receipt.json")
    verify_signed(stored)
    exact(stored, receipt(root, reports), "receipt_mutation")
    require(
        read(root / "historical-preservation.json")["files_preserved"]
        == historical(Path(repo)),
        "historical_preservation_mutation",
    )
    return reports["result.json"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("emit", "replay", "validate"))
    parser.add_argument("--directory", type=Path, default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=os.environ.get("ISLAND_OPLS_SOURCE_DIR", str(SOURCE_DIR)),
    )
    parser.add_argument(
        "--data-sha256", default=os.environ.get("ISLAND_OPLS_DATA_SHA256")
    )
    parser.add_argument(
        "--settings-sha256", default=os.environ.get("ISLAND_OPLS_SETTINGS_SHA256")
    )
    parser.add_argument("--snapshot-sha256")
    args = parser.parse_args()
    external = {
        "source_dir": args.source_dir,
        "data_sha256": args.data_sha256,
        "settings_sha256": args.settings_sha256,
    }
    try:
        scope()
        if args.command == "emit":
            value = emit(args.repo, args.directory, **external)
        else:
            value = replay(
                args.repo,
                args.directory,
                args.snapshot_sha256,
                **external,
            )
        print(json.dumps(value, sort_keys=True))
        return 0
    except (ValueError, TypeError, KeyError, OSError, RuntimeError) as error:
        print(
            json.dumps(
                {
                    "status": "blocked",
                    "decision": "parameter_bundle_incomplete",
                    "reason": str(error),
                    "scope": SCOPE,
                },
                sort_keys=True,
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
