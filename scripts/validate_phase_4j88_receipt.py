"""Strict offline J88 reconstruction, scope checks and historical preservation."""

import argparse
import hashlib
import json
import lzma
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j72_c1_integrity import verify_events
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
from scripts.phase_4j88_opls import (
    BASE,
    BRANCH,
    IMPLEMENTATION,
    ROOT,
    SCOPE,
    historical,
    replay,
)


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
    exact(
        tested["files"],
        {n: pin(repo / n) for n in IMPLEMENTATION},
        "tested_implementation_mutation",
    )
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE, tested["commit"]],
        cwd=repo,
        capture_output=True,
        check=True,
    )
    for n, sha in tested["files"].items():
        raw = subprocess.check_output(
            ["git", "show", tested["commit"] + ":" + n], cwd=repo
        )
        require(hashlib.sha256(raw).hexdigest() == sha, "tested_commit_blob_mutation")
    with lzma.open(root / "validation/terminal.xml.xz", "rb") as stream:
        cases = list(ET.parse(stream).iter("testcase"))
    nodes = [
        c.attrib["classname"].replace(".", "/") + ".py::" + c.attrib["name"]
        for c in cases
    ]
    expected = read(root / "validation/collection.json")["nodes"]
    require(
        len(nodes) == len(set(nodes)) and set(nodes) == set(expected),
        "terminal_collection_mutation",
    )
    files = {
        "tests/test_phase_4j" + n + ".py"
        for n in (
            "73",
            "74",
            "75",
            "76",
            "78",
            "79",
            "80",
            "81",
            "82",
            "83",
            "84",
            "85",
            "86",
            "87",
            "88",
        )
    }
    files |= {
        "tests/" + n + ".py"
        for n in (
            "test_psmiles",
            "test_polymer_builder",
            "test_polymer_termination",
            "test_atom_typing",
            "test_atom_typing_rules",
            "test_atom_typing_compatibility",
            "test_parameter_inventory",
            "test_parameter_assignment_graph_only",
        )
    }
    require(
        {n.split("::")[0] for n in nodes} == files,
        "explicit_regression_inventory_required",
    )
    groups = {
        g: {k: 0 for k in ("passed", "skipped", "failure", "error")}
        for g in ("J88_focused", "regressions")
    }
    for c, node in zip(cases, nodes, strict=True):
        outcome = next(
            (k for k in ("failure", "error", "skipped") if c.find(k) is not None),
            "passed",
        )
        group = (
            "J88_focused"
            if node.startswith("tests/test_phase_4j88.py::")
            else "regressions"
        )
        groups[group][outcome] += 1
    require(
        all(g["failure"] == g["error"] == 0 for g in groups.values()),
        "terminal_test_failure",
    )
    checks = read(root / "validation/checks.json")
    required = {
        "terminal_tests",
        "ruff",
        "compileall",
        "pip_supported",
        "J81_receipt",
        "J82_receipt",
        "J83_receipt",
        "J84_receipt",
        "J85_receipt",
        "J86_receipt",
        "J87_receipt",
        "J81_replay",
        "J82_replay",
        "J83_replay",
        "J84_replay",
        "J85_replay",
        "J86_replay",
        "J87_replay",
        "J80_receipt",
        "J79_receipt",
        "J78_receipt",
        "J76_receipt",
        "C1_receipt",
        "offline_replay",
        "relocated_replay",
        "diff",
    }
    require(
        len(checks) == len(required) and {c["name"] for c in checks} == required,
        "required_check_inventory",
    )
    require(
        all(c["returncode"] == 0 and (root / c["log"]).is_file() for c in checks),
        "verification_check_failed",
    )
    runs = read(root / "validation/fenced-runs.json")
    require(runs and len(runs) == len(set(runs)), "fenced_attempt_inventory")
    return {
        "tested_implementation_commit": tested["commit"],
        "groups": groups,
        "unique_terminal_tests": len(nodes),
        "fence_events": verify_events([root / n for n in runs]),
        "software_attempts": [
            {"run": n, "returncode": read(root / n / "result.json")["returncode"]}
            for n in runs
        ],
        "scientific_checks_passed": 0,
    }


def seal(repo, root):
    scope()
    require(
        not (root / "receipt.json").exists() and not (root / "manifest.json").exists(),
        "do_not_overwrite_seal",
    )
    result = replay(root, repo=repo)
    hashes = historical(repo)
    checks = verification(repo, root)
    write(root / "historical-files.json", hashes)
    manifest = signed(
        {
            "schema": "island.j88-manifest.v1",
            "base": BASE,
            "branch": BRANCH,
            "scope": SCOPE,
            "implementation": {n: pin(repo / n) for n in IMPLEMENTATION},
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j88-receipt.v1",
            "base": BASE,
            "branch": BRANCH,
            "scope": SCOPE,
            "decision": result["decision"],
            "result": result,
            "verification": checks,
            "snapshot_sha256": pin(root / "captured-input.json"),
            "manifest_sha256": pin(root / "manifest.json"),
            "historical_files_preserved": len(hashes),
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "status": "sealed",
        "decision": result["decision"],
        "historical_files_preserved": len(hashes),
    }


def validate(repo, root):
    scope()
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    for obj, name in ((r, "receipt"), (m, "manifest")):
        verify_signed(obj)
        require(
            obj["schema"] == "island.j88-" + name + ".v1"
            and obj["base"] == BASE
            and obj["branch"] == BRANCH,
            "receipt_identity_mutation",
        )
        exact(obj["scope"], SCOPE, "resigned_scope_or_readiness_mutation")
    require(
        pin(root / "manifest.json") == r["manifest_sha256"], "manifest_hash_mutation"
    )
    exact(
        m["implementation"],
        {n: pin(repo / n) for n in IMPLEMENTATION},
        "implementation_mutation",
    )
    exact(m["artifacts"], artifacts(root), "artifact_hash_mutation")
    result = replay(root, repo=repo)
    exact(r["result"], result, "resigned_receipt_semantic_mutation")
    require(
        r["decision"] == result["decision"]
        and r["snapshot_sha256"] == pin(root / "captured-input.json"),
        "decision_or_snapshot_mutation",
    )
    exact(
        r["verification"],
        verification(repo, root),
        "resigned_validation_metadata_mutation",
    )
    hashes = historical(repo)
    exact(read(root / "historical-files.json"), hashes, "historical_inventory_mutation")
    require(r["historical_files_preserved"] == len(hashes), "historical_count_mutation")
    return {
        "status": "passed",
        "decision": r["decision"],
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files_preserved": len(hashes),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    a = p.parse_args()
    try:
        print(
            json.dumps(
                (seal if a.seal else validate)(a.repo, a.directory),
                sort_keys=True,
            )
        )
        return 0
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as e:
        print(
            json.dumps(
                {
                    "status": "blocked",
                    "decision": "parameter_bundle_incomplete",
                    "reason": str(e),
                }
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
