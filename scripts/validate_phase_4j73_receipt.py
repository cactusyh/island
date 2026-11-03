"""Offline, fenced J73 request reconstruction and strict semantic validation."""

import argparse
import hashlib
import json
import lzma
import os
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j73_source_intake import (
    BASE,
    BRANCH,
    ROOT,
    SCOPE,
    IntakeError,
    build,
    exact,
    parents,
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
    "scripts/phase_4j73_scope.py",
    "scripts/phase_4j73_source_intake.py",
    "scripts/validate_phase_4j73_receipt.py",
    "tests/test_phase_4j73.py",
    "docs/phase_4j73.md",
)


def artifacts(root):
    return {
        str(p.relative_to(root)): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file()
        and p.name not in {"receipt.json", "manifest.json"}
        and "__pycache__" not in p.parts
    }


def history(repo, root):
    parents(repo, root)
    inventory = read(root / "historical-files.json")
    frozen = read(root / "frozen-inputs.json")
    require(
        pin(root / "historical-files.json") == frozen["historical_files_sha256"],
        "historical_inventory_mutation",
    )
    for name, sha in inventory.items():
        path = repo / name
        # Git symlinks store their link text, not target bytes.
        actual = (
            hashlib.sha256(os.readlink(path).encode()).hexdigest()
            if path.is_symlink()
            else pin(path)
        )
        require(actual == sha, "historical_byte_mutation:" + name)
    return len(inventory)


def verification(root):
    from scripts.phase_4j72_c1_integrity import verify_events

    validation = root / "validation"
    records = read(validation / "checks.json")
    expected = {
        "focused_tests",
        "C1_receipt",
        "ruff",
        "compileall",
        "pip",
        "diff",
        "offline_replay",
        "relocated_replay",
    }
    require(
        len(records) == len(expected)
        and {r["name"] for r in records} == expected
        and all(type(r["returncode"]) is int and r["returncode"] == 0 for r in records),
        "required_fenced_checks",
    )
    for record in records:
        require((root / record["log"]).is_file(), "missing_verification_log")
    with lzma.open(validation / "focused.xml.xz", "rb") as stream:
        cases = list(ET.parse(stream).iter("testcase"))
    collection = read(validation / "collection.json")["nodes"]
    actual = [
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
        for c in cases
    ]
    require(
        len(actual) == len(set(actual)) and set(actual) == set(collection),
        "unique_terminal_test_coverage",
    )
    counts = dict.fromkeys(("passed", "skipped", "failure", "error"), 0)
    for case in cases:
        counts[
            next(
                (
                    k
                    for k in ("error", "failure", "skipped")
                    if case.find(k) is not None
                ),
                "passed",
            )
        ] += 1
    require(counts["failure"] == counts["error"] == 0, "terminal_test_failure")
    runs = [root / n for n in read(validation / "fenced-runs.json")]
    events = verify_events(runs)
    # J73 verification permits only Python, read-only Git and Ruff execs, and
    # no executed scientific imports. The fixtures use marker-only fake code.
    for run in runs:
        data = lzma.decompress((run / "events.jsonl.xz").read_bytes()).decode()
        for line in data.splitlines():
            event = json.loads(line)
            if event["action"] == "exec_allowed":
                args = event["command"]
                require(
                    not any(a in ("-S", "-I", "-E") for a in args[1:2]),
                    "unfenced_Python_bypass",
                )
    return {
        "terminal_unique_tests": len(actual),
        "software_counts": counts,
        "fence_verification": events,
        "scientific_controls": "explicit scope skips; never successful scientific checks",
    }


def semantic(receipt, manifest):
    verify_signed(receipt)
    verify_signed(manifest)
    require(
        set(receipt)
        == {
            "schema",
            "identity",
            "base_commit",
            "branch",
            "scope",
            "manifest_sha256",
            "historical_files_preserved",
            "verification",
            "offline_reconstruction",
        },
        "receipt_schema",
    )
    require(
        set(manifest)
        == {
            "schema",
            "identity",
            "base_commit",
            "scope",
            "implementation_files",
            "artifacts",
            "frozen_inputs_sha256",
        },
        "manifest_schema",
    )
    require(
        receipt["schema"] == "island.j73-receipt.v1"
        and manifest["schema"] == "island.j73-manifest.v1"
        and receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["branch"] == BRANCH,
        "receipt_base_schema",
    )
    exact(receipt["scope"], SCOPE, "receipt_scope_semantics")
    exact(manifest["scope"], SCOPE, "manifest_scope_semantics")
    exact(
        receipt["offline_reconstruction"],
        {
            "offline": True,
            "instances": 378,
            "unresolved_counts": {
                "analog_unqualified": 294,
                "convention_unresolved": 60,
                "source_missing": 24,
            },
            "scope": SCOPE,
        },
        "nested_receipt_scope",
    )


def audit(repo, root):
    return {
        "offline_reconstruction": replay(repo, root),
        "historical_files_preserved": history(repo, root),
        "verification": verification(root),
    }


def seal(repo, root):
    result = audit(repo, root)
    manifest = signed(
        {
            "schema": "island.j73-manifest.v1",
            "base_commit": BASE,
            "scope": SCOPE,
            "implementation_files": {n: pin(repo / n) for n in FILES},
            "artifacts": artifacts(root),
            "frozen_inputs_sha256": pin(root / "frozen-inputs.json"),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j73-receipt.v1",
            "base_commit": BASE,
            "branch": BRANCH,
            "scope": SCOPE,
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantic(receipt, manifest)
    write(root / "receipt.json", receipt)
    return {"status": "sealed", **result}


def validate(repo, root):
    receipt, manifest = read(root / "receipt.json"), read(root / "manifest.json")
    semantic(receipt, manifest)
    require(
        receipt["manifest_sha256"] == pin(root / "manifest.json")
        and manifest["frozen_inputs_sha256"] == pin(root / "frozen-inputs.json")
        and manifest["implementation_files"] == {n: pin(repo / n) for n in FILES}
        and manifest["artifacts"] == artifacts(root),
        "manifest_inventory_mutation",
    )
    result = audit(repo, root)
    exact({k: receipt[k] for k in result}, result, "receipt_reconstruction")
    return {"status": "passed", "receipt_sha256": pin(root / "receipt.json"), **result}


def main():
    try:
        scope()  # A missing C1 installation stops before any intake/check.
    except RuntimeError as error:
        print(json.dumps({"status": "blocked", "reason": str(error)}))
        return 2
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=ROOT)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--build", action="store_true")
    mode.add_argument("--replay", action="store_true")
    mode.add_argument("--seal", action="store_true")
    args = parser.parse_args()
    function = (
        build
        if args.build
        else replay
        if args.replay
        else seal
        if args.seal
        else validate
    )
    try:
        print(
            json.dumps(
                function(args.repo.resolve(), args.directory.resolve()), sort_keys=True
            )
        )
    except (IntakeError, OSError, KeyError, TypeError, ValueError) as error:
        print(json.dumps({"status": "rejected", "reason": str(error)}))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
