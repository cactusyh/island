"""Strict J70 receipt, semantic declarations, exact inventories and byte preservation."""

import argparse
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_augmentation import identity
from island.pcff_qm_pilot.hmds_precision import read
from scripts.intake_phase_4j70 import BASE, REQUIRED, ROOT, reports
from scripts.replay_phase_4j70 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin
from scripts.validate_phase_4j64_receipt import artifacts

FILES = (
    "src/island/pcff_qm_pilot/source_intake.py",
    "scripts/intake_phase_4j70.py",
    "scripts/replay_phase_4j70.py",
    "scripts/validate_phase_4j70_receipt.py",
    "tests/test_phase_4j70.py",
    "docs/phase_4j70.md",
)
DECLARATIONS = {
    "decision": "parameter_incomplete",
    "source_complete": False,
    "converter_complete": False,
    "parameter_assignment_capability": False,
    "execution_profile_created": False,
    "registered_profiles": [],
    "accepted_coefficients": 0,
    "automatic_assignment": False,
    "new_QM": 0,
    "MD_steps": 0,
    "production_fits": 0,
    "reserved_validation_activated": False,
    "reserved_validation_consumed": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "scientific_validation": "not_established",
    "chemical_transfer": "not_established",
    "native_PCFF_equivalence": "not_established",
    "original_author_artifact_status": "not_yet_supplied",
    "author_confirmations": "28 user-attested",
}


def declarations(record):
    identity(record)
    for k, v in DECLARATIONS.items():
        require(
            type(record.get(k)) is type(v) and record[k] == v,
            "J70_scientific_declaration",
            field=k,
        )


def semantics(receipt, manifest):
    declarations(receipt)
    identity(manifest)
    require(
        set(receipt)
        == set(DECLARATIONS)
        | {
            "schema",
            "identity",
            "base_commit",
            "required_J69",
            "branch",
            "manifest_sha256",
            "software_counts",
            "terminal_unique_tests",
            "numerical_reconstruction",
            "historical_files_preserved",
        },
        "J70_receipt_schema",
    )
    require(
        receipt["schema"] == "island.j70-receipt.v1"
        and receipt["base_commit"] == BASE
        and receipt["required_J69"] == REQUIRED
        and receipt["branch"]
        == "codex/phase-4j70-authoritative-pcff-missing-row-intake",
        "J70_ancestry",
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
        }
        and manifest["schema"] == "island.numerical-evidence-manifest.v1"
        and manifest["base_commit"] == BASE,
        "J70_manifest_schema",
    )
    require(
        receipt["numerical_reconstruction"]
        == {
            "offline": True,
            "instances": 378,
            "observations_retained": 26,
            "source_complete": False,
            "converter_complete": False,
            "parameter_assignment_capability": False,
            "new_QM": 0,
            "MD_steps": 0,
            "production_fits": 0,
            "accepted_coefficients": 0,
            "validation_untouched": True,
        },
        "J70_nested_semantics",
    )
    require(
        all(
            type(v)
            is (
                int
                if k
                in (
                    "instances",
                    "observations_retained",
                    "new_QM",
                    "MD_steps",
                    "production_fits",
                    "accepted_coefficients",
                )
                else bool
            )
            for k, v in receipt["numerical_reconstruction"].items()
        ),
        "J70_nested_semantic_types",
    )


def audit(repo, root):
    declarations(read(root / "decision.json"))
    require(
        read(root / "decision.json") == reports(repo, root)["decision.json"],
        "J70_decision_schema_semantics",
    )
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / n) == h for n, h in old.items()),
        "J70_historical_byte_mutation",
    )
    numerical = replay(repo, root)
    counts, nodes = junit(root, read(root / "validation/junit-files.json"))
    require(
        len(nodes) == len(set(nodes))
        and Counter(nodes) == Counter(read(root / "collection.json")["nodes"]),
        "J70_unique_terminal_coverage",
    )
    require(
        counts["failure"] == counts["error"] == 0
        and all(r["returncode"] == 0 for r in read(root / "validation/checks.json")),
        "J70_required_checks",
    )
    return {
        "software_counts": counts,
        "terminal_unique_tests": len(nodes),
        "numerical_reconstruction": numerical,
        "historical_files_preserved": len(old),
    }, old


def seal(repo, root):
    result, old = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.numerical-evidence-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "historical_files": old,
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            **DECLARATIONS,
            "schema": "island.j70-receipt.v1",
            "base_commit": BASE,
            "required_J69": REQUIRED,
            "branch": "codex/phase-4j70-authoritative-pcff-missing-row-intake",
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantics(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {"receipt_sha256": pin(root / "receipt.json"), **result}


def validate(repo, root):
    receipt = read(root / "receipt.json")
    manifest = read(root / "manifest.json")
    semantics(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["implementation_files"] == {n: pin(repo / n) for n in FILES}
        and manifest["artifacts"] == artifacts(root),
        "J70_manifest_inventory",
    )
    result, old = audit(repo, root)
    require(
        manifest["historical_files"] == old
        and all(receipt[k] == v for k, v in result.items()),
        "J70_receipt_reconstruction",
    )
    return {"status": "passed", "receipt_sha256": pin(root / "receipt.json"), **result}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path, nargs="?", default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    a = p.parse_args()
    print((seal if a.seal else validate)(a.repo.resolve(), a.directory.resolve()))
