"""Offline, fail-closed J75 receipt reconstruction; never opens source assets."""

import argparse
import json
import lzma
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j75_source_path import (
    BASE,
    BRANCH,
    FREEZE_SHA,
    MAIN,
    ROOT,
    SCOPE,
    IntakeError,
    exact,
    inputs,
    old,
    pin,
    read,
    replay,
    require,
    scope,
    signed,
    verify_signed,
    write,
)

FILES = (
    "scripts/phase_4j75_source_path.py",
    "scripts/validate_phase_4j75_receipt.py",
    "tests/test_phase_4j75.py",
    "docs/phase_4j75.md",
)


def artifacts(root):
    result = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file() or "__pycache__" in p.parts:
            continue
        require(
            p.suffix.lower() not in {".pdf", ".frc", ".zip", ".dat", ".data"},
            "source_asset_redistribution_forbidden",
        )
        n = p.relative_to(root).as_posix()
        if n not in {"receipt.json", "manifest.json"}:
            result[n] = pin(p)
    return result


def history(repo, root):
    inputs(repo, root)
    entries = read(root / "historical-files.json")
    for name, sha in entries.items():
        p = repo / name
        actual = old.sha_bytes(str(p.readlink()).encode()) if p.is_symlink() else pin(p)
        require(actual == sha, "historical_byte_mutation:" + name)
    return len(entries)


def coverage(root):
    with lzma.open(root / "validation/focused.xml.xz", "rb") as stream:
        cases = list(ET.parse(stream).iter("testcase"))
    nodes = [
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
        for c in cases
    ]
    collection = read(root / "validation/collection.json")["nodes"]
    require(
        len(nodes) == len(set(nodes)) and set(nodes) == set(collection),
        "unique_terminal_coverage",
    )
    counts = dict.fromkeys(("passed", "skipped", "failure", "error"), 0)
    for case in cases:
        key = next(
            (k for k in ("failure", "error", "skipped") if case.find(k) is not None),
            "passed",
        )
        counts[key] += 1
    require(counts["failure"] == counts["error"] == 0, "terminal_software_failure")
    return {
        "terminal_unique_tests": len(nodes),
        "counts": counts,
        "scientific_checks_passed": 0,
        "skipped_controls": "Live science, evaluator, converter and validation controls are scope skips; broader historical scientific suites were not executed",
    }


def verification(root):
    from scripts.phase_4j72_c1_integrity import verify_events
    from scripts.validate_phase_4j74_receipt import PARSER_PIP_DIAGNOSTIC

    checks = read(root / "validation/checks.json")
    names = {
        "focused_tests",
        "ruff",
        "compileall",
        "pip_supported",
        "pip_parser",
        "J74_receipt",
        "J73_receipt",
        "C1_receipt",
        "offline_replay",
        "relocated_replay",
        "source_originals",
        "diff",
    }
    require(
        len(checks) == len(names) and {c["name"] for c in checks} == names,
        "required_verification_inventory",
    )
    for c in checks:
        require(
            type(c["returncode"]) is int
            and c["returncode"] == (1 if c["name"] == "pip_parser" else 0),
            "verification_returncode",
        )
        require((root / c["log"]).is_file(), "missing_check_log")
        if c["name"] == "pip_parser":
            require(
                lzma.decompress((root / c["log"]).read_bytes())
                == PARSER_PIP_DIAGNOSTIC.encode(),
                "existing_parser_dependency_failure_must_be_disclosed",
            )
    runs = [root / n for n in read(root / "validation/fenced-runs.json")]
    require(runs and len(runs) == len(set(runs)), "fenced_runs_required")
    events = verify_events(runs)
    for run in runs:
        require((run / "result.json").is_file(), "incomplete_fenced_run")
    return {
        **coverage(root),
        "fence_events": events,
        "parser_pip_check": "existing torch/fsspec/sympy dependency failures; no environment changes",
        "initial_sequence_compliant": True,
        "source_parsing_only": True,
    }


def semantic(receipt, manifest):
    verify_signed(receipt)
    verify_signed(manifest)
    for doc, schema in (
        (receipt, "island.j75-receipt.v1"),
        (manifest, "island.j75-manifest.v1"),
    ):
        require(
            doc["schema"] == schema
            and doc["base_commit"] == BASE
            and doc["fetched_origin_main"] == MAIN,
            "receipt_base_identity",
        )
        exact(doc["scope"], SCOPE, "resigned_scientific_scope")
    require(
        set(receipt)
        == {
            "identity",
            "schema",
            "base_commit",
            "fetched_origin_main",
            "branch",
            "decision",
            "scope",
            "manifest_sha256",
            "source_resolution",
            "historical_files_preserved",
            "verification",
        }
        and set(manifest)
        == {
            "identity",
            "schema",
            "base_commit",
            "fetched_origin_main",
            "scope",
            "frozen_inputs_sha256",
            "implementation_files",
            "artifacts",
        },
        "strict_receipt_manifest_schema",
    )
    require(
        receipt["branch"] == BRANCH
        and receipt["decision"] in {"parameter_incomplete", "source_path_unavailable"},
        "receipt_decision",
    )
    exact(
        receipt["source_resolution"]["counts"],
        {"analog_only": 294, "convention_blocked": 60, "source_missing": 24},
        "coverage_semantics",
    )
    require(
        receipt["source_resolution"]["resolved_native_rows"] == 0,
        "fabricated_native_resolution",
    )
    require(
        manifest["frozen_inputs_sha256"] == FREEZE_SHA
        and set(manifest["implementation_files"]) == set(FILES),
        "explicit_implementation_and_freeze",
    )
    require(
        receipt["verification"]["initial_sequence_compliant"] is True,
        "source_read_before_fence",
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
            "schema": "island.j75-manifest.v1",
            "base_commit": BASE,
            "fetched_origin_main": MAIN,
            "scope": SCOPE,
            "frozen_inputs_sha256": FREEZE_SHA,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j75-receipt.v1",
            "base_commit": BASE,
            "fetched_origin_main": MAIN,
            "branch": BRANCH,
            "decision": result["source_resolution"]["decision"],
            "scope": SCOPE,
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantic(receipt, manifest)
    require(not (root / "receipt.json").exists(), "do_not_overwrite_receipt")
    write(root / "receipt.json", receipt)
    return {"validation_status": "sealed", **result}


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantic(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json"),
        "manifest_hash_mutation",
    )
    exact(
        manifest["implementation_files"],
        {n: pin(repo / n) for n in FILES},
        "implementation_mutation",
    )
    exact(manifest["artifacts"], artifacts(root), "artifact_inventory_mutation")
    result = audit(repo, root)
    exact({k: receipt[k] for k in result}, result, "resigned_receipt_reconstruction")
    require(
        receipt["decision"] == result["source_resolution"]["decision"],
        "resigned_decision",
    )
    return {
        "validation_status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        **result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--seal", action="store_true")
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    try:
        scope()
        result = (
            seal(args.repo, args.directory)
            if args.seal
            else replay(args.repo, args.directory)
            if args.replay
            else validate(args.repo, args.directory)
        )
        print(json.dumps(result, sort_keys=True))
        return 0
    except (
        IntakeError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        RuntimeError,
    ) as error:
        print(json.dumps({"status": "blocked", "reason": str(error)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
