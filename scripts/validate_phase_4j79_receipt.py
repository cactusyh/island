"""Strict offline J79 semantic replay, historical byte audit and receipt seal."""

import argparse
import hashlib
import json
import lzma
import re
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j72_c1_integrity import verify_events
from scripts.phase_4j79_literature import (
    BASE,
    BRANCH,
    FREEZE_SHA,
    IMPLEMENTATION,
    ROOT,
    SCOPE,
    digest,
    exact,
    pin,
    read,
    replay,
    require,
    scope,
    signed,
    verify_signed,
    write,
)


def history(repo, root):
    replay(root)
    historical = read(root / "historical-files.json")
    for n, sha in historical.items():
        p = repo / n
        require(
            (digest(str(p.readlink())) if p.is_symlink() else pin(p)) == sha,
            "historical_byte_mutation:" + n,
        )
    for n, sha in read(root / "frozen-inputs.json")["original_inputs"].items():
        require(pin(repo / n) == sha, "original_input_mutation:" + n)
    return len(historical)


def artifacts(root):
    return {
        p.relative_to(root).as_posix(): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file()
        and "__pycache__" not in p.parts
        and p.name not in ("receipt.json", "manifest.json")
    }


def verification(repo, root):
    tested = read(root / "validation/tested-implementation.json")
    require(
        bool(re.fullmatch("[0-9a-f]{40}", tested["commit"])), "tested_commit_required"
    )
    exact(
        tested["files"],
        {n: pin(repo / n) for n in IMPLEMENTATION},
        "tested_implementation_changed",
    )
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE, tested["commit"]],
        cwd=repo,
        check=True,
        capture_output=True,
    )
    for n, sha in tested["files"].items():
        blob = subprocess.check_output(
            ["git", "show", tested["commit"] + ":" + n], cwd=repo
        )
        require(hashlib.sha256(blob).hexdigest() == sha, "tested_commit_blob_mutation")
    with lzma.open(root / "validation/terminal.xml.xz", "rb") as f:
        cases = list(ET.parse(f).iter("testcase"))
    nodes = [
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
        for c in cases
    ]
    require(
        len(nodes) == len(set(nodes))
        and set(nodes) == set(read(root / "validation/collection.json")["nodes"]),
        "unique_terminal_test_coverage",
    )
    groups = {
        n: dict.fromkeys(("passed", "skipped", "failure", "error"), 0)
        for n in ("J79_focused", "regressions")
    }
    for n, c in zip(nodes, cases, strict=True):
        group = (
            "J79_focused"
            if n.startswith("tests/test_phase_4j79.py::")
            else "regressions"
        )
        outcome = next(
            (k for k in ("failure", "error", "skipped") if c.find(k) is not None),
            "passed",
        )
        groups[group][outcome] += 1
    require(
        all(g["failure"] == g["error"] == 0 for g in groups.values()),
        "terminal_test_failure",
    )
    required = {
        "terminal_tests",
        "ruff",
        "compileall",
        "pip_supported",
        "J78_receipt",
        "J76_receipt",
        "C1_receipt",
        "offline_replay",
        "relocated_replay",
        "diff",
    }
    checks = read(root / "validation/checks.json")
    require(
        len(checks) == len(required) and {c["name"] for c in checks} == required,
        "explicit_check_inventory",
    )
    for c in checks:
        require(
            type(c["returncode"]) is int
            and c["returncode"] == 0
            and (root / c["log"]).is_file(),
            "terminal_check_failed",
        )
    runs = read(root / "validation/fenced-runs.json")
    require(runs and len(runs) == len(set(runs)), "fenced_attempt_inventory")
    for n in runs:
        result = read(root / n / "result.json")
        # Failed software probes are retained, never relabeled as science.
        require(type(result["returncode"]) is int, "missing_attempt_outcome")
    return {
        "tested_implementation_commit": tested["commit"],
        "groups": groups,
        "unique_terminal_tests": len(nodes),
        "fence_events": verify_events([root / n for n in runs]),
        "scientific_checks_passed": 0,
        "scientific_controls": "not launched; real QM/MD/converter/fitting/evaluator controls are scope skips",
        "all_software_attempts": [
            {"run": n, "returncode": read(root / n / "result.json")["returncode"]}
            for n in runs
        ],
    }


def semantic(receipt, manifest):
    for obj, schema in (
        (receipt, "island.j79-receipt.v1"),
        (manifest, "island.j79-manifest.v1"),
    ):
        verify_signed(obj)
        require(
            obj["schema"] == schema and obj["base"] == BASE and obj["branch"] == BRANCH,
            "receipt_identity_mutation",
        )
        exact(obj["scope"], SCOPE, "resigned_readiness_or_scope_mutation")
    require(
        receipt["decision"] == "literature_reconciliation_partial"
        and manifest["freeze_sha256"] == FREEZE_SHA,
        "fabricated_qualification_or_source_freeze",
    )


def audit(repo, root):
    return {
        "result": replay(root),
        "historical_files_preserved": history(repo, root),
        "verification": verification(repo, root),
    }


def seal(repo, root):
    require(
        not (root / "receipt.json").exists() and not (root / "manifest.json").exists(),
        "do_not_overwrite_sealed_attempt",
    )
    result = audit(repo, root)
    m = signed(
        {
            "schema": "island.j79-manifest.v1",
            "base": BASE,
            "branch": BRANCH,
            "scope": SCOPE,
            "freeze_sha256": FREEZE_SHA,
            "implementation": {n: pin(repo / n) for n in IMPLEMENTATION},
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", m)
    r = signed(
        {
            "schema": "island.j79-receipt.v1",
            "base": BASE,
            "branch": BRANCH,
            "decision": "literature_reconciliation_partial",
            "scope": SCOPE,
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantic(r, m)
    write(root / "receipt.json", r)
    return {"status": "sealed", **result}


def validate(repo, root):
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    semantic(r, m)
    require(
        r["manifest_sha256"] == pin(root / "manifest.json"), "manifest_hash_mutation"
    )
    exact(
        m["implementation"],
        {n: pin(repo / n) for n in IMPLEMENTATION},
        "implementation_mutation",
    )
    exact(m["artifacts"], artifacts(root), "artifact_inventory_mutation")
    result = audit(repo, root)
    exact({k: r[k] for k in result}, result, "resigned_receipt_semantic_mutation")
    return {"status": "passed", "receipt_sha256": pin(root / "receipt.json"), **result}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    p.add_argument("--replay", action="store_true")
    a = p.parse_args()
    try:
        scope()
        value = (
            seal(a.repo, a.directory)
            if a.seal
            else replay(a.directory)
            if a.replay
            else validate(a.repo, a.directory)
        )
        print(json.dumps(value, sort_keys=True))
        return 0
    except (
        ValueError,
        OSError,
        KeyError,
        TypeError,
        RuntimeError,
        subprocess.CalledProcessError,
    ) as e:
        print(json.dumps({"status": "blocked", "reason": str(e)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
