"""Strict J71 receipt validation, including re-signed scientific declarations."""

import argparse
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_augmentation import identity
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.literature_intake import semantic_declarations
from scripts.extract_phase_4j71 import BASE, BRANCH, DECLARATIONS, REQUIRED, ROOT
from scripts.replay_phase_4j71 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts

FILES = (
    "src/island/pcff_qm_pilot/literature_intake.py",
    "scripts/extract_phase_4j71.py",
    "scripts/read_phase_4j71_sources.py",
    "scripts/replay_phase_4j71.py",
    "scripts/validate_phase_4j71_receipt.py",
    "tests/test_phase_4j71.py",
    "docs/phase_4j71.md",
)


def semantics(receipt, manifest):
    semantic_declarations(receipt)
    identity(manifest)
    require(
        set(receipt)
        == set(DECLARATIONS)
        | {
            "identity",
            "schema",
            "base_commit",
            "required_J70",
            "branch",
            "manifest_sha256",
            "software_counts",
            "terminal_unique_tests",
            "offline_reconstruction",
            "historical_files_preserved",
        },
        "J71_receipt_schema",
    )
    require(
        receipt["schema"] == "island.j71-receipt.v1"
        and receipt["base_commit"] == BASE
        and receipt["required_J70"] == REQUIRED
        and receipt["branch"] == BRANCH,
        "J71_ancestry",
    )
    require(
        set(manifest)
        == {
            "schema",
            "identity",
            "base_commit",
            "implementation_files",
            "historical_files",
            "artifacts",
            "scientific_scope",
        }
        and manifest["schema"] == "island.j71-evidence-manifest.v1"
        and manifest["base_commit"] == BASE,
        "J71_manifest_schema",
    )
    require(
        manifest["scientific_scope"] == DECLARATIONS
        and all(
            type(manifest["scientific_scope"].get(k)) is type(v)
            for k, v in DECLARATIONS.items()
        ),
        "J71_manifest_scientific_semantics",
    )
    r = receipt["offline_reconstruction"]
    require(
        r
        == {
            "offline": True,
            "extracted_candidate_rows": 467,
            "nonblank_numeric_cells": 1470,
            "J70_instances": 378,
            "resolved_instances": 0,
            "decision": "source_conflict",
            "parameter_completeness": "incomplete",
            "converter_executed": False,
            "static_complete_closure": False,
            "new_QM": 0,
            "QM_fits": 0,
            "MD_steps": 0,
            "accepted_coefficients": 0,
            "validation_untouched": True,
        },
        "J71_nested_semantics",
    )
    for key in (
        "new_QM",
        "QM_fits",
        "MD_steps",
        "accepted_coefficients",
        "extracted_candidate_rows",
        "nonblank_numeric_cells",
        "J70_instances",
        "resolved_instances",
    ):
        require(type(r[key]) is int, "J71_nested_types")
    for key in (
        "offline",
        "converter_executed",
        "static_complete_closure",
        "validation_untouched",
    ):
        require(type(r[key]) is bool, "J71_nested_types")


def audit(repo, root):
    result = replay(repo, root)
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == h for n, h in old.items()),
        "J71_historical_byte_mutation",
    )
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        len(nodes) == len(set(nodes))
        and Counter(nodes) == Counter(read(root / "collection.json")["nodes"]),
        "J71_unique_terminal_coverage",
    )
    require(
        counts["failure"] == counts["error"] == 0
        and all(c["returncode"] == 0 for c in read(root / "validation/checks.json")),
        "J71_required_checks",
    )
    return {
        "software_counts": counts,
        "terminal_unique_tests": len(nodes),
        "offline_reconstruction": result,
        "historical_files_preserved": len(old),
    }, old


def seal(repo, root):
    result, old = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.j71-evidence-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "historical_files": old,
            "artifacts": artifacts(root),
            "scientific_scope": DECLARATIONS,
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j71-receipt.v1",
            **DECLARATIONS,
            "base_commit": BASE,
            "required_J70": REQUIRED,
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantics(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {"status": "sealed", "receipt_sha256": pin(root / "receipt.json"), **result}


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantics(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["implementation_files"] == {n: pin(repo / n) for n in FILES}
        and manifest["artifacts"] == artifacts(root),
        "J71_manifest_inventory",
    )
    result, old = audit(repo, root)
    require(
        manifest["historical_files"] == old
        and all(receipt[k] == v for k, v in result.items()),
        "J71_receipt_reconstruction",
    )
    return {"status": "passed", "receipt_sha256": pin(root / "receipt.json"), **result}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--seal", action="store_true")
    args = parser.parse_args()
    print(
        (seal if args.seal else validate)(args.repo.resolve(), args.directory.resolve())
    )
