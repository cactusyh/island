"""Strict J76 offline replay, implementation-commit proof and receipt checks."""

import argparse
import json
import lzma
import re
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j76_source_resolution import (
    BASE,
    BRANCH,
    FREEZE_SHA,
    MAIN,
    ROOT,
    SCOPE,
    STATES,
    IntakeError,
    exact,
    inputs,
    j75,
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
    "scripts/phase_4j76_source_resolution.py",
    "scripts/validate_phase_4j76_receipt.py",
    "tests/test_phase_4j76.py",
    "docs/phase_4j76.md",
)
EXPECTED_COUNTS = {
    **dict.fromkeys(STATES, 0),
    "analog_only": 294,
    "source_missing": 24,
    "unresolved": 60,
}


def artifacts(root):
    result = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file() or "__pycache__" in p.parts:
            continue
        require(
            p.suffix.lower() not in {".pdf", ".frc", ".zip", ".dat", ".data"},
            "source_asset_redistribution_forbidden",
        )
        n = str(p.relative_to(root))
        if n not in {"manifest.json", "receipt.json"}:
            result[n] = pin(p)
    return result


def history(repo, root):
    inputs(repo, root)
    files = read(root / "historical-files.json")
    for n, sha in files.items():
        p = repo / n
        actual = (
            j75.old.sha_bytes(str(p.readlink()).encode()) if p.is_symlink() else pin(p)
        )
        require(actual == sha, "historical_byte_mutation:" + n)
    return len(files)


def tested_commit(repo, root):
    record = read(root / "validation/tested-implementation.json")
    sha = record["commit"]
    require(bool(re.fullmatch(r"[0-9a-f]{40}", sha)), "tested_commit_required")
    exact(
        record["implementation_files"],
        {n: pin(repo / n) for n in FILES},
        "tested_implementation_changed",
    )
    # All commands are read-only and independently mediated by C1. A fabricated
    # hash or a real commit with different source bytes cannot pass this proof.
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE, sha],
        cwd=repo,
        check=True,
        capture_output=True,
    )
    for n in FILES:
        raw = subprocess.check_output(["git", "show", sha + ":" + n], cwd=repo)
        require(
            j75.old.sha_bytes(raw) == record["implementation_files"][n],
            "tested_commit_blob_mismatch:" + n,
        )
    return sha


def coverage(root):
    with lzma.open(root / "validation/terminal.xml.xz", "rb") as stream:
        cases = list(ET.parse(stream).iter("testcase"))
    nodes = [
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
        for c in cases
    ]
    require(
        len(nodes) == len(set(nodes))
        and set(nodes) == set(read(root / "validation/collection.json")["nodes"]),
        "unique_terminal_coverage",
    )
    counts = {
        g: dict.fromkeys(("passed", "skipped", "failure", "error"), 0)
        for g in ("J76_focused", "J73_J75_regressions")
    }
    for node, case in zip(nodes, cases, strict=True):
        category = next(
            (k for k in ("failure", "error", "skipped") if case.find(k) is not None),
            "passed",
        )
        group = (
            "J76_focused"
            if node.startswith("tests/test_phase_4j76.py::")
            else "J73_J75_regressions"
        )
        counts[group][category] += 1
    require(
        all(c["failure"] == c["error"] == 0 for c in counts.values()),
        "terminal_software_failure",
    )
    return {
        "unique_terminal_tests": len(nodes),
        "groups": counts,
        "scientific_checks_passed": 0,
        "scope_skips": "Live QM/MD/converter/evaluator/fitting/reserved controls not launched; no scientific pass claim",
    }


def verification(repo, root):
    from scripts.phase_4j72_c1_integrity import verify_events

    checks = read(root / "validation/checks.json")
    required = {
        "terminal_tests",
        "ruff",
        "compileall",
        "pip_supported",
        "J75_receipt",
        "J74_receipt",
        "J73_receipt",
        "C1_receipt",
        "offline_replay",
        "relocated_replay",
        "source_originals",
        "diff",
    }
    require(
        len(checks) == len(required) and {c["name"] for c in checks} == required,
        "required_check_inventory",
    )
    for c in checks:
        require(
            type(c["returncode"]) is int
            and c["returncode"] == 0
            and (root / c["log"]).is_file(),
            "required_terminal_check_failed",
        )
    runs = [root / n for n in read(root / "validation/fenced-runs.json")]
    require(runs and len(runs) == len(set(runs)), "fenced_execution_records_required")
    for run in runs:
        require((run / "result.json").is_file(), "incomplete_fenced_attempt")
    return {
        **coverage(root),
        "tested_implementation_commit": tested_commit(repo, root),
        "fence_events": verify_events(runs),
        "source_fence_preceded_access": True,
        "PDF_parser_used": False,
        "historical_parser_pip_conflicts": "J75 retained; J76 uses only supported standard-library intake environment",
    }


def semantic(receipt, manifest):
    verify_signed(receipt)
    verify_signed(manifest)
    for doc, schema in (
        (receipt, "island.j76-receipt.v1"),
        (manifest, "island.j76-manifest.v1"),
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
            "row_resolution",
            "historical_files_preserved",
            "verification",
            "tested_implementation_commit",
        },
        "receipt_schema",
    )
    require(
        set(manifest)
        == {
            "identity",
            "schema",
            "base_commit",
            "fetched_origin_main",
            "scope",
            "implementation_files",
            "artifacts",
            "frozen_inputs_sha256",
        },
        "manifest_schema",
    )
    require(
        receipt["branch"] == BRANCH and receipt["decision"] == "parameter_incomplete",
        "decision_or_branch_mutation",
    )
    exact(
        receipt["row_resolution"]["status_counts"],
        EXPECTED_COUNTS,
        "resigned_native_row_promotion",
    )
    require(
        receipt["tested_implementation_commit"]
        == receipt["verification"]["tested_implementation_commit"],
        "tested_commit_inconsistency",
    )
    require(
        manifest["frozen_inputs_sha256"] == FREEZE_SHA
        and set(manifest["implementation_files"]) == set(FILES),
        "explicit_implementation_inventory",
    )


def audit(repo, root):
    return {
        "row_resolution": replay(repo, root),
        "historical_files_preserved": history(repo, root),
        "verification": verification(repo, root),
    }


def seal(repo, root):
    require(
        not (root / "receipt.json").exists() and not (root / "manifest.json").exists(),
        "do_not_overwrite_sealed_attempt",
    )
    result = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.j76-manifest.v1",
            "base_commit": BASE,
            "fetched_origin_main": MAIN,
            "scope": SCOPE,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "artifacts": artifacts(root),
            "frozen_inputs_sha256": FREEZE_SHA,
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j76-receipt.v1",
            "base_commit": BASE,
            "fetched_origin_main": MAIN,
            "branch": BRANCH,
            "decision": "parameter_incomplete",
            "scope": SCOPE,
            "manifest_sha256": pin(root / "manifest.json"),
            "tested_implementation_commit": result["verification"][
                "tested_implementation_commit"
            ],
            **result,
        }
    )
    semantic(receipt, manifest)
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
        subprocess.CalledProcessError,
    ) as error:
        print(json.dumps({"status": "blocked", "reason": str(error)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
