"""J92 fenced real-source intake; no default scientific graph or provider."""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

from island.forcefields.assignment_core.schema import is_hash
from island.forcefields.opls_providers.intake92 import READY, SOURCE_DIR, intake

FLAGS = READY


def input_identity(source_dir=SOURCE_DIR, *, manifest_sha256=None):
    return {
        "source_directory_declared": str(source_dir),
        "manifest_sha256": manifest_sha256,
    }


from scripts.phase_4j73_source_intake import exact, pin, read, require, signed, write
from scripts.phase_4j76_source_resolution import scope

BASE = "43945345da32aa86c4217538126af67b2937c24e"
BRANCH = "codex/phase-4j92-lammps-opls-bundle-intake"
ROOT = Path("docs/evidence/phase_4j92")
IMPLEMENTATION = (
    "src/island/forcefields/opls_providers/intake92.py",
    "src/island/forcefields/opls_providers/lammps92.py",
    "scripts/phase_4j92_intake.py",
    "scripts/validate_phase_4j92_receipt.py",
    "tests/test_phase_4j92.py",
    "docs/phase_4j92.md",
)
SCOPE = {
    "J92_QM_executions": 0,
    "J92_MD_steps": 0,
    "converter_executions": 0,
    "physical_fitting_executions": 0,
    "parameter_generation_executions": 0,
    "energy_force_evaluations": 0,
    "minimization_executions": 0,
    "OpenMM_Context_executions": 0,
    "production_executions": 0,
    "reserved_validation_executions": 0,
    **FLAGS,
    "J72_historical_violation_and_C1_quarantine": "preserved_unchanged",
    "synthetic_controls_are_scientific_data": False,
}
REPORTS = ("result.json", "per-graph-report.json", "missing-terms-and-blockers.json")


def reconstruct(repo=None, **external):
    del repo
    scope()
    result = intake(**external).report
    snapshot = {
        "schema": "island.j92.external-input-pins.v1",
        "provided": input_identity(**external),
        "raw_source_records_persisted": False,
    }
    reports = {
        "result.json": signed({**result, "execution_scope": SCOPE}),
        "per-graph-report.json": signed(
            {
                "requested_graph": result["requested_graph"],
                "qualified_graph": result["qualified_graph"],
                "settings": result["requested_settings"],
                "coverage": result["coverage"],
                "interaction_instances": result["interaction_instances"],
                "scope": "one_exact_finite_capped_oligomer_only",
                "chemical_family_transfer": "not_established",
            }
        ),
        "missing-terms-and-blockers.json": signed(
            {
                "missing_terms": result["missing_terms"],
                "convention_blockers": result["conventions"].get("blockers", []),
                "blocking_reasons": result["blockers"],
            }
        ),
    }
    return snapshot, reports


def intake_receipt(root, reports):
    return signed(
        {
            "schema": "island.j92.intake-receipt.v1",
            "base": BASE,
            "snapshot_sha256": pin(root / "captured-input.json"),
            "reports": {n: pin(root / n) for n in REPORTS},
            "decision": reports["result.json"]["decision"],
            "scope": SCOPE,
        }
    )


def emit(repo, root, **external):
    scope()
    require(not root.exists(), "do_not_overwrite_intake_attempt")
    inputs, reports = reconstruct(repo, **external)
    root.mkdir(parents=True)
    write(root / "captured-input.json", inputs)
    for name, report in reports.items():
        write(root / name, report)
    write(root / "qualification-receipt.json", intake_receipt(root, reports))
    return {
        "result": reports["result.json"],
        "snapshot_sha256": pin(root / "captured-input.json"),
    }


def replay(root, snapshot_sha=None, *, repo=None, **external):
    scope()
    inputs = read(root / "captured-input.json")
    if inputs["provided"] != input_identity():
        require(is_hash(snapshot_sha), "external_snapshot_pin_required")
    if snapshot_sha is not None:
        require(
            pin(root / "captured-input.json") == snapshot_sha, "snapshot_hash_mutation"
        )
    exact(
        inputs["provided"]["manifest_sha256"],
        input_identity(**external)["manifest_sha256"],
        "explicit_external_originals_and_pins_required_for_replay",
    )
    fresh, reports = reconstruct(repo, **external)
    fresh["provided"]["source_directory_declared"] = inputs["provided"][
        "source_directory_declared"
    ]
    exact(inputs, fresh, "external_sources_unavailable_or_mutated")
    for name, report in reports.items():
        exact(read(root / name), report, "resigned_semantic_mutation:" + name)
    exact(
        read(root / "qualification-receipt.json"),
        intake_receipt(root, reports),
        "intake_receipt_mutation",
    )
    return reports["result.json"]


def historical(repo):
    entries = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    hashes = {}
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
        hashes[name] = hashlib.sha256(raw).hexdigest()
    return hashes


def external_arguments(parser):
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=os.environ.get("ISLAND_OPLS_SOURCE_DIR", str(SOURCE_DIR)),
    )
    parser.add_argument(
        "--manifest-sha256", default=os.environ.get("ISLAND_OPLS_MANIFEST_SHA256")
    )
    parser.add_argument("--snapshot-sha256")


def arguments(a):
    return {"source_dir": a.source_dir, "manifest_sha256": a.manifest_sha256}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("emit", "ingest", "replay"))
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    external_arguments(p)
    a = p.parse_args()
    try:
        scope()
        external = arguments(a)
        if a.command == "replay":
            value = replay(a.directory, a.snapshot_sha256, repo=a.repo, **external)
        else:
            value = emit(a.repo, a.directory, **external)
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
