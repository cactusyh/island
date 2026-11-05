"""Offline strict semantic validation of the separate J78 extension proposal."""

import argparse
import json
import lzma
import re
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j72_c1_integrity import verify_events
from scripts.phase_4j78_extension import (
    BASE,
    BRANCH,
    FREEZE_SHA,
    IMPLEMENTATION,
    ROOT,
    SCOPE,
    digest,
    exact,
    inputs,
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
    inputs(root)
    old = read(root / "historical-files.json")
    for name, sha in old.items():
        p = repo / name
        require(
            (digest(str(p.readlink())) if p.is_symlink() else pin(p)) == sha,
            "historical_native_or_evidence_byte_mutation:" + name,
        )
    frozen = read(root / "frozen-inputs.json")
    for n, sha in frozen["original_inputs"].items():
        require(pin(repo / n) == sha, "historical_input_mutation:" + n)
    return len(old)


def artifacts(root):
    return {
        p.relative_to(root).as_posix(): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file()
        and "__pycache__" not in p.parts
        and p.name not in ("manifest.json", "receipt.json")
    }


def verification(repo, root):
    r = read(root / "validation/tested-implementation.json")
    require(bool(re.fullmatch("[0-9a-f]{40}", r["commit"])), "tested_commit_required")
    exact(
        r["files"],
        {n: pin(repo / n) for n in IMPLEMENTATION},
        "tested_implementation_changed",
    )
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE, r["commit"]],
        cwd=repo,
        check=True,
        capture_output=True,
    )
    for n, sha in r["files"].items():
        import hashlib

        blob = subprocess.check_output(["git", "show", r["commit"] + ":" + n], cwd=repo)
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
        "unique_terminal_coverage",
    )
    groups = {
        n: dict.fromkeys(("passed", "skipped", "failure", "error"), 0)
        for n in ("J78_focused", "regressions")
    }
    for n, c in zip(nodes, cases, strict=True):
        group = (
            "J78_focused"
            if n.startswith("tests/test_phase_4j78.py::")
            else "regressions"
        )
        outcome = next(
            (k for k in ("failure", "error", "skipped") if c.find(k) is not None),
            "passed",
        )
        groups[group][outcome] += 1
    require(
        all(g["failure"] == g["error"] == 0 for g in groups.values()),
        "software_test_failure",
    )
    checks = read(root / "validation/checks.json")
    required = {
        "terminal_tests",
        "ruff",
        "compileall",
        "pip_supported",
        "J76_receipt",
        "J75_receipt",
        "C1_receipt",
        "offline_replay",
        "relocated_replay",
        "diff",
    }
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
    runs = [root / n for n in read(root / "validation/fenced-runs.json")]
    require(
        runs
        and len(runs) == len(set(runs))
        and all((p / "result.json").is_file() for p in runs),
        "fenced_attempt_inventory",
    )
    return {
        "tested_implementation_commit": r["commit"],
        "unique_terminal_tests": len(nodes),
        "groups": groups,
        "fence_events": verify_events(runs),
        "scientific_checks_passed": 0,
        "synthetic_controls": "bounded matrix fitting only; no physical coefficients or observations",
        "scope_skips": "real_QM_MD_converter_evaluator_fitting_reserved_validation_not_launched",
    }


def audit(repo, root):
    return {
        "result": replay(root),
        "historical_files_preserved": history(repo, root),
        "verification": verification(repo, root),
    }


def semantic(receipt, manifest):
    verify_signed(receipt)
    verify_signed(manifest)
    require(
        receipt["schema"] == "island.j78-receipt.v1"
        and manifest["schema"] == "island.j78-manifest.v1",
        "receipt_schema",
    )
    for d in (receipt, manifest):
        require(
            d["base"] == BASE and d["branch"] == BRANCH, "receipt_base_branch_mutation"
        )
        exact(d["scope"], SCOPE, "resigned_scientific_readiness_mutation")
    require(
        receipt["decision"] == "derived_profile_blocked"
        and manifest["freeze_sha256"] == FREEZE_SHA,
        "native_completion_or_qualification_promotion",
    )


def seal(repo, root):
    require(
        not (root / "receipt.json").exists() and not (root / "manifest.json").exists(),
        "do_not_overwrite_sealed_attempt",
    )
    result = audit(repo, root)
    m = signed(
        {
            "schema": "island.j78-manifest.v1",
            "base": BASE,
            "branch": BRANCH,
            "scope": SCOPE,
            "freeze_sha256": FREEZE_SHA,
            "implementation": {n: pin(repo / n) for n in IMPLEMENTATION},
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", m)
    receipt = signed(
        {
            "schema": "island.j78-receipt.v1",
            "base": BASE,
            "branch": BRANCH,
            "decision": "derived_profile_blocked",
            "scope": SCOPE,
            "manifest_sha256": pin(root / "manifest.json"),
            **result,
        }
    )
    semantic(receipt, m)
    write(root / "receipt.json", receipt)
    return {"validation_status": "sealed", **result}


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
    exact(
        {k: r[k] for k in result}, result, "resigned_qualification_or_receipt_mutation"
    )
    return {
        "validation_status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        **result,
    }


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
        RuntimeError,
        OSError,
        KeyError,
        TypeError,
        subprocess.CalledProcessError,
    ) as e:
        print(json.dumps({"status": "blocked", "reason": str(e)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
