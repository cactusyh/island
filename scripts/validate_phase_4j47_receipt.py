"""Validate J47 evidence and exact base bytes; a blocked gate must stay nonzero."""

import argparse
import json
import subprocess
from collections import Counter
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.external_pdms import LIMITS, digest
from island.forcefields.pdms_opls import PROFILE
from scripts.validate_phase_4j46_receipt import pin

BASE = "9e7d9765d0248e14320b4214db8831771f06e84e"
BRANCH = "codex/phase-4j47-pdms-lammps-end-to-end"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def history(repo, commit=BASE):
    """Hash the complete base Git tree, including symlink target bytes."""
    entries = []
    for item in subprocess.check_output(
        ["git", "ls-tree", "-r", "-z", commit], cwd=repo
    ).split(b"\0"):
        if item:
            fields, name = item.split(b"\t", 1)
            _, kind, oid = fields.split()
            require(kind == b"blob", "Unsupported historical object")
            entries.append((name.decode(), oid.decode()))
    objects = list(dict.fromkeys(oid for _, oid in entries))
    raw = subprocess.check_output(
        ["git", "cat-file", "--batch"],
        cwd=repo,
        input=("\n".join(objects) + "\n").encode(),
    )
    hashes, cursor = {}, 0
    for oid in objects:
        end = raw.index(b"\n", cursor)
        header = raw[cursor:end].split()
        require(header[:2] == [oid.encode(), b"blob"], "Historical Git object changed")
        size = int(header[2])
        start = end + 1
        hashes[oid] = sha256(raw[start : start + size]).hexdigest()
        cursor = start + size + 1
    return {name: hashes[oid] for name, oid in entries}


def read(path):
    value = json.loads(Path(path).read_text())
    require(
        value["identity"]
        == digest({k: v for k, v in value.items() if k != "identity"}),
        "Changed identity: " + str(path),
    )
    return value


def validate(root, *, repo, expected_sha256=None):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    if expected_sha256 is not None:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "Independent receipt hash mismatch",
        )
    r = read(root / "receipt.json")
    m = read(root / "manifest.json")
    require(
        r["schema"] == "island.j47-receipt.v1"
        and m["schema"] == "island.j47-manifest.v1",
        "Wrong schema",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong base/branch",
    )
    require(
        r["profile_identity"] == PROFILE and m["profile_identities"] == [PROFILE],
        "Changed profile/family",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established"
        and r["limitations"] == LIMITS,
        "Promoted readiness",
    )
    require(
        r["end_to_end_acceptance_gate"] == 1
        and r["scientific_acceptance_gate"] == 1
        and r["native_pcff_acceptance_gate"] == 1,
        "Promoted incomplete workflow",
    )
    require(
        r["MD_executed"] is False
        and r["NPT_executed"] is False
        and r["new_physical_parameters"] == []
        and r["new_profiles"] == [],
        "Changed physical authority",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Changed manifest")
    historical = history(repo)
    require(m["historical_files"] == historical, "Historical inventory changed")
    for name, expected in historical.items():
        require(pin(repo / name) == expected, "Historical bytes changed: " + name)
    for collection, base in ((m["implementation_files"], repo), (m["artifacts"], root)):
        for name, expected in collection.items():
            path = base / name
            require(
                path.resolve().is_relative_to(base) and pin(path) == expected,
                "Changed file: " + name,
            )
    independent = read(root / "independent-inputs.json")
    require(
        set(independent["cases"]) == {"0", "2", "4", "6"}
        and all(
            c["graph_type_charge_mass_and_source_rows_reconstructed"] is True
            and c["source_rows"] == 17
            for c in independent["cases"].values()
        ),
        "Missing independent source/type/charge reconstruction",
    )
    s = read(root / "audit/summary.json")
    require(
        s["identity"] == r["summary_identity"]
        and s["source_assets"] == m["source_assets"],
        "Changed source audit",
    )
    require(set(s["cases"]) == {"0", "2", "4", "6"}, "Missing acceptance case")
    require(
        s["end_to_end_acceptance_gate"] == 1 and s["exact_static_parity"] is False,
        "Promoted parity",
    )
    for name, case in s["cases"].items():
        require(
            case["acceptance_passed"] is False
            and case["MD_executed"] is False
            and case["NPT_executed"] is False,
            "Promoted incomplete case",
        )
        require(
            case["relocation_verified"] is True
            and case["fresh_process_reconstruction_verified"] is True,
            "Missing reconstruction",
        )
        require(
            case["identical_physical_inputs"] is True
            and case["identical_numerical_long_range_settings"] is False,
            "Changed numerical contract",
        )
        require(case["raw_rows_reconstructed"] == 17, "Missing source reconstruction")
        for kind in ("bond", "angle", "proper_torsion", "lj"):
            c = case["component_comparison"][kind]
            require(
                c["bound"] == 3e-4 and c["passed"] is True,
                "Changed bonded/LJ comparison",
            )
        require(
            case["component_comparison"]["total"]["passed"] is False
            and case["component_comparison"]["forces"]["passed"] is False,
            "Lost raw failure",
        )
        require(
            case["component_comparison"]["virial"]["passed"] is False, "Promoted virial"
        )
        require(
            case["OpenMM_finite_difference"]["comparisons"] == 24
            and case["LAMMPS_finite_difference"]["comparisons"] == 24,
            "Missing FD measurements",
        )
    counts = Counter()
    nodes = []
    for name in r["junit_files"]:
        for case in ET.parse(root / name).iter("testcase"):
            klass = case.attrib["classname"]
            require(klass.startswith("tests."), "Wrong test class")
            nodes.append(klass.replace(".", "/") + ".py::" + case.attrib["name"])
            counts["tests"] += 1
            for flag in ("failure", "error", "skipped"):
                counts[flag] += int(case.find(flag) is not None)
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected)
        and all(n == 1 for n in Counter(nodes).values()),
        "Missing/duplicated suite coverage",
    )
    require(
        dict(counts) == r["full_suite_counts"]
        and counts["failure"] == counts["error"] == 0,
        "Changed suite results",
    )
    require(
        all(code == 0 for code in r["check_returncodes"].values()),
        "Failed required validation",
    )
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(historical),
        "full_suite_counts": dict(counts),
        "end_to_end_acceptance_gate": 1,
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
