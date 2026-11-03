"""Strict offline verification of the blocked J74 source-resolution receipt."""

import argparse
import json
import lzma
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j74_source_intake import (
    BASE,
    BRANCH,
    J73,
    MAIN,
    ROOT,
    SCOPE,
    IntakeError,
    exact,
    inputs,
    pin,
    read,
    replay,
    require,
    scope,
    sha_bytes,
    signed,
    verify_signed,
    write,
)

FILES = (
    "scripts/phase_4j74_source_intake.py",
    "scripts/validate_phase_4j74_receipt.py",
    "tests/test_phase_4j74.py",
    "docs/phase_4j74.md",
)
PARSER_PIP_DIAGNOSTIC = (
    "torch 2.5.1 requires fsspec, which is not installed.\n"
    'torch 2.5.1 has requirement sympy==1.13.1; python_version >= "3.9", but you have sympy 1.14.0.\n'
)


def artifacts(root):
    require(
        not any(
            p.is_file()
            and p.suffix.lower() in {".pdf", ".frc", ".zip", ".dat", ".data"}
            for p in root.rglob("*")
        ),
        "source_asset_redistribution_forbidden_in_J74",
    )
    return {
        str(path.relative_to(root)): pin(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
        and str(path.relative_to(root)) not in {"receipt.json", "manifest.json"}
        and "__pycache__" not in path.parts
    }


def history(repo, root):
    inputs(repo, root)
    paths = read(root / "historical-files.json")
    for name, expected in paths.items():
        path = repo / name
        actual = (
            sha_bytes(str(path.readlink()).encode()) if path.is_symlink() else pin(path)
        )
        require(actual == expected, "historical_byte_mutation:" + name)
    return len(paths)


def verification(root):
    from scripts.phase_4j72_c1_integrity import verify_events

    validation = root / "validation"
    records = read(validation / "checks.json")
    expected = {
        "focused_tests",
        "ruff",
        "compileall",
        "pip_supported",
        "pip_parser",
        "J73_receipt",
        "C1_receipt",
        "offline_replay",
        "relocated_replay",
        "source_originals",
        "diff",
    }
    require(
        len(records) == len(expected)
        and {r["name"] for r in records} == expected
        and all(
            type(r["returncode"]) is int
            and r["returncode"] == (1 if r["name"] == "pip_parser" else 0)
            for r in records
        ),
        "required_fenced_verification",
    )
    require(all((root / r["log"]).is_file() for r in records), "missing_check_log")
    parser_log = root / next(r["log"] for r in records if r["name"] == "pip_parser")
    raw = (
        lzma.decompress(parser_log.read_bytes())
        if parser_log.suffix == ".xz"
        else parser_log.read_bytes()
    )
    require(
        raw == PARSER_PIP_DIAGNOSTIC.encode(),
        "parser_pip_failure_must_remain_disclosed",
    )
    with lzma.open(validation / "focused.xml.xz", "rb") as stream:
        cases = list(ET.parse(stream).iter("testcase"))
    actual = [
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
        for c in cases
    ]
    collection = read(validation / "collection.json")["nodes"]
    require(
        len(actual) == len(set(actual)) and set(actual) == set(collection),
        "unique_terminal_coverage",
    )
    counts = dict.fromkeys(("passed", "skipped", "failure", "error"), 0)
    for case in cases:
        category = next(
            (k for k in ("failure", "error", "skipped") if case.find(k) is not None),
            "passed",
        )
        counts[category] += 1
    require(counts["failure"] == counts["error"] == 0, "terminal_software_failure")
    directories = [root / name for name in read(validation / "fenced-runs.json")]
    events = verify_events(directories)
    return {
        "terminal_unique_tests": len(actual),
        "software_counts": counts,
        "fence_events": events,
        "scope_skips": "Live science/converter/evaluator controls explicitly skipped, never reported passed",
        "initial_readonly_inspection_before_J74_preflight": True,
        "initial_sequence_compliant": False,
        "capture_and_terminal_verification_fenced": True,
        "parser_environment_pip_check": {
            "returncode": 1,
            "status": "existing_torch_dependency_conflicts_retained; no_dependency_changes",
            "diagnostic": PARSER_PIP_DIAGNOSTIC,
        },
    }


def semantic(receipt, manifest):
    verify_signed(receipt)
    verify_signed(manifest)
    require(
        set(receipt)
        == {
            "identity",
            "schema",
            "status",
            "branch",
            "base_commit",
            "previous_J73_commit",
            "fetched_origin_main",
            "scope",
            "manifest_sha256",
            "source_resolution",
            "historical_files_preserved",
            "verification",
        },
        "receipt_schema",
    )
    require(
        set(manifest)
        == {
            "identity",
            "schema",
            "base_commit",
            "previous_J73_commit",
            "fetched_origin_main",
            "scope",
            "implementation_files",
            "artifacts",
            "frozen_inputs_sha256",
        },
        "manifest_schema",
    )
    for document, schema in (
        (receipt, "island.j74-receipt.v1"),
        (manifest, "island.j74-manifest.v1"),
    ):
        require(
            document["schema"] == schema
            and document["base_commit"] == BASE
            and document["previous_J73_commit"] == J73
            and document["fetched_origin_main"] == MAIN,
            "receipt_base_or_schema",
        )
        exact(document["scope"], SCOPE, "resigned_scientific_scope")
    require(
        receipt["branch"] == BRANCH
        and receipt["status"] == "blocked"
        and receipt["source_resolution"]["source_resolution_status"] == "blocked",
        "blocked_source_status_required",
    )
    require(
        receipt["source_resolution"]["resolved_rows"] == 0,
        "fabricated_exact_resolution",
    )
    exact(receipt["source_resolution"]["scope"], SCOPE, "nested_scope")
    exact(
        receipt["source_resolution"]["counts"],
        {"analog_only": 294, "convention_blocked": 60, "source_missing": 24},
        "coverage_semantics",
    )
    require(
        receipt["verification"]["initial_readonly_inspection_before_J74_preflight"]
        is True
        and receipt["verification"]["initial_sequence_compliant"] is False,
        "inspection_deviation_must_not_be_erased",
    )
    exact(
        receipt["verification"]["parser_environment_pip_check"],
        {
            "returncode": 1,
            "status": "existing_torch_dependency_conflicts_retained; no_dependency_changes",
            "diagnostic": PARSER_PIP_DIAGNOSTIC,
        },
        "parser_environment_failure_must_not_be_erased",
    )
    require(
        set(manifest["implementation_files"]) == set(FILES),
        "explicit_implementation_inventory",
    )


def audit(repo, root):
    return {
        "source_resolution": replay(repo, root),
        "historical_files_preserved": history(repo, root),
        "verification": verification(root),
    }


def seal(repo, root):
    result = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.j74-manifest.v1",
            "base_commit": BASE,
            "previous_J73_commit": J73,
            "fetched_origin_main": MAIN,
            "scope": SCOPE,
            "implementation_files": {name: pin(repo / name) for name in FILES},
            "artifacts": artifacts(root),
            "frozen_inputs_sha256": pin(root / "frozen-inputs.json"),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j74-receipt.v1",
            "status": "blocked",
            "branch": BRANCH,
            "base_commit": BASE,
            "previous_J73_commit": J73,
            "fetched_origin_main": MAIN,
            "scope": SCOPE,
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantic(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {
        "validation_status": "sealed",
        "scientific_source_status": "blocked",
        **result,
    }


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantic(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["frozen_inputs_sha256"] == pin(root / "frozen-inputs.json")
        and manifest["implementation_files"]
        == {name: pin(repo / name) for name in FILES}
        and manifest["artifacts"] == artifacts(root),
        "manifest_inventory_mutation",
    )
    result = audit(repo, root)
    exact({k: receipt[k] for k in result}, result, "resigned_receipt_reconstruction")
    return {
        "validation_status": "passed",
        "scientific_source_status": "blocked",
        "receipt_sha256": pin(root / "receipt.json"),
        **result,
    }


def main():
    try:
        scope()
    except RuntimeError as error:
        print(json.dumps({"status": "blocked", "reason": str(error), "scope": SCOPE}))
        return 2
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--seal", action="store_true")
    modes.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    function = seal if args.seal else replay if args.replay else validate
    try:
        print(
            json.dumps(
                function(args.repo.resolve(), args.directory.resolve()), sort_keys=True
            )
        )
    except (IntakeError, OSError, ValueError, TypeError, KeyError) as error:
        print(json.dumps({"status": "rejected", "reason": str(error)}))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
