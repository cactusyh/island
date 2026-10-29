"""Validate J46 report/manifest, full-suite accounting and historical bytes."""

import argparse
import json
import subprocess
import tarfile
from collections import Counter
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree

from island.forcefields.external_pdms import LIMITS
from island.forcefields.selection import capability_registry
from island.forcefields.selection.inventory import digest

BASE = "fdf8e8058173faa9c013703aa4910b3065820e7c"
HISTORICAL = "8257a2ef25042494d0e4541b8da117cd8e97072a"
BRANCH = "codex/phase-4j46-unified-forcefield-selection"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pin(path):
    path = Path(path)
    # Git stores a symlink's target text, not the contents it points to.
    raw = str(path.readlink()).encode() if path.is_symlink() else path.read_bytes()
    return sha256(raw).hexdigest()


def read(path):
    data = json.loads(Path(path).read_text())
    require(
        data["identity"] == digest({k: v for k, v in data.items() if k != "identity"}),
        "Changed identity: " + str(path),
    )
    return data


def history(repo):
    """Hash Git object bytes independently of checkout filters and file mtimes."""
    raw = subprocess.check_output(["git", "ls-tree", "-r", "-z", HISTORICAL], cwd=repo)
    entries = []
    for item in raw.split(b"\0"):
        if item:
            fields, name = item.split(b"\t", 1)
            _mode, kind, oid = fields.split()
            require(kind == b"blob", "Unsupported historical object")
            entries.append((name.decode(), oid.decode()))
    objects = list(dict.fromkeys(oid for _, oid in entries))
    stream = subprocess.check_output(
        ["git", "cat-file", "--batch"],
        input=("\n".join(objects) + "\n").encode(),
        cwd=repo,
    )
    cursor, hashes = 0, {}
    for oid in objects:
        end = stream.index(b"\n", cursor)
        header = stream[cursor:end].split()
        require(header[:2] == [oid.encode(), b"blob"], "Changed Git object")
        size = int(header[2])
        start = end + 1
        hashes[oid] = sha256(stream[start : start + size]).hexdigest()
        cursor = start + size + 1
    return {name: hashes[oid] for name, oid in entries}


def validate(directory, *, repo=None, expected_sha256=None):
    directory = Path(directory).resolve()
    repo = Path.cwd() if repo is None else Path(repo)
    if expected_sha256 is not None:
        require(
            pin(directory / "receipt.json") == expected_sha256,
            "Independent receipt hash mismatch",
        )
    receipt = read(directory / "receipt.json")
    manifest = read(directory / "manifest.json")
    require(
        receipt["schema"] == "island.j46-receipt.v1"
        and manifest["schema"] == "island.j46-manifest.v1",
        "Wrong schema",
    )
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["branch"] == BRANCH,
        "Wrong base/branch",
    )
    require(receipt["historical_snapshot"] == HISTORICAL, "Wrong historical snapshot")
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["limitations"] == LIMITS,
        "Promoted readiness",
    )
    require(
        receipt["production_md_executed"] is False
        and receipt["new_physical_parameters"] == [],
        "Unexpected physical parameters/MD",
    )
    require(
        pin(directory / "manifest.json") == receipt["manifest_sha256"],
        "Changed manifest",
    )
    for name, expected in manifest["implementation_files"].items():
        require(pin(repo / name) == expected, "Changed implementation: " + name)
    expected_history = history(repo)
    require(
        manifest["historical_files"] == expected_history, "Changed historical inventory"
    )
    for name, expected in expected_history.items():
        require(pin(repo / name) == expected, "Changed historical bytes: " + name)
    for name, expected in manifest["artifacts"].items():
        path = (directory / name).resolve()
        require(
            path.is_relative_to(directory) and pin(path) == expected,
            "Changed artifact: " + name,
        )
    summary = read(directory / "summary.json")
    require(
        summary["identity"] == receipt["summary_identity"],
        "Changed source audit summary",
    )
    require(
        summary["production_validated"] is False
        and summary["simulation_readiness"] == "not_established",
        "Promoted audit readiness",
    )
    require(
        receipt["registry"] == capability_registry(), "Changed registry declaration"
    )
    with tarfile.open(directory / "audit.tar.xz", "r:xz") as archive:
        require(
            archive.getnames() == ["audit.json"], "Changed audit archive membership"
        )
        audit = json.load(archive.extractfile("audit.json"))
    require(
        audit["identity"]
        == digest({k: v for k, v in audit.items() if k != "identity"})
        == summary["audit_identity"]
        == receipt["audit_identity"],
        "Changed full audit",
    )
    require(
        audit["registry"] == receipt["registry"]
        and audit["source_assets"]
        == summary["source_assets"]
        == manifest["source_assets"]
        and audit["numerical_md_executed"] is False,
        "Changed source/backend scope",
    )
    projected = {}
    for name, case in audit["cases"].items():
        require(
            case["identity"]
            == digest({k: v for k, v in case.items() if k != "identity"}),
            "Changed assignment report",
        )
        require(
            case["final_graph_identity"] == digest(case["final_graph"]),
            "Changed final graph",
        )
        require(
            case["assigned_terms"]
            == sorted(
                r["id"] for r in case["requirements"] if r["status"] == "assigned"
            )
            and case["blocked_terms"]
            == sorted(
                r["id"]
                for r in case["requirements"]
                if r["status"] not in ("assigned", "not_applicable")
            ),
            "Changed interaction inventory",
        )
        projected[name] = {
            k: case[k]
            for k in (
                "status",
                "profile_identity",
                "source_kind",
                "stages",
                "diagnostics",
                "blocked_terms",
                "identity",
            )
        }
    require(
        projected == summary["cases"], "Summary differs from full assignment evidence"
    )
    for case in summary["cases"].values():
        s = case["stages"]
        require(
            s["scientifically_validated"] is False
            and s["numerically_verified"] is False,
            "Promoted scientific/numerical status",
        )
        require(
            not s["executable"]
            or case["status"] == "prepared"
            and all(
                s[k]
                for k in ("typing_complete", "charge_complete", "parameter_complete")
            )
            and not case["blocked_terms"],
            "Incomplete execution case",
        )
    counts = Counter()
    observed = []
    for name in receipt["junit_files"]:
        tree = ElementTree.parse(directory / name)
        for case in tree.iter("testcase"):
            classname = case.attrib["classname"]
            require(classname.startswith("tests."), "Unexpected test class")
            # The suite contains module-level test functions only.
            observed.append(classname.replace(".", "/") + ".py::" + case.attrib["name"])
            counts["tests"] += 1
            for flag in ("failure", "error", "skipped"):
                counts[flag] += int(case.find(flag) is not None)
    collection = json.loads((directory / "collection.json").read_text())["nodes"]
    require(
        Counter(observed) == Counter(collection)
        and all(v == 1 for v in Counter(observed).values()),
        "Incomplete or repeated full-suite coverage",
    )
    require(
        dict(counts) == receipt["full_suite_counts"]
        and counts["failure"] == counts["error"] == 0,
        "Changed full-suite outcome",
    )
    require(
        all(v == 0 for v in receipt["check_returncodes"].values()),
        "Failed required check",
    )
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(expected_history),
        "full_suite_counts": dict(counts),
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--sha256")
    a = p.parse_args(argv)
    print(
        json.dumps(
            validate(a.directory, repo=a.repo, expected_sha256=a.sha256), sort_keys=True
        )
    )


if __name__ == "__main__":
    main()
