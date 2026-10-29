"""Validate canonical registry, raw blocked evidence, full suite and historical bytes."""

import argparse
import json
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.external_pdms import LIMITS
from island.forcefields.opls_registry import canonical_registry
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import digest
from scripts.validate_phase_4j47_receipt import history, pin, require

BASE = "fdf117aadf613912786fdfbebe653ad58695bea3"
BRANCH = "codex/phase-4j48-canonical-opls-registry"


def read(path):
    p = json.loads(Path(path).read_text())
    require(
        p["identity"] == digest({k: v for k, v in p.items() if k != "identity"}),
        "Changed identity: " + str(path),
    )
    return p


def validate(root, *, repo, expected_sha256=None):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    if expected_sha256 is not None:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "Independent receipt hash mismatch",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["schema"] == "island.j48-receipt.v1"
        and m["schema"] == "island.j48-manifest.v1",
        "Wrong schema",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong base/branch",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established"
        and r["limitations"] == LIMITS,
        "Readiness promoted",
    )
    require(
        r["software_selection_acceptance_gate"] == 0
        and r["lammps_end_to_end_acceptance_gate"]
        == r["scientific_acceptance_gate"]
        == 1,
        "Altered acceptance",
    )
    require(
        r["MD_executed"] is False
        and r["NPT_executed"] is False
        and not r["new_physical_parameters"]
        and not r["new_profiles"],
        "Altered physical authority",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Changed manifest")
    require(
        read(root / "canonical-registry.json") == canonical_registry(),
        "Changed canonical registry",
    )
    require(
        m["profile_identities"]
        == [e["profile_identity"] for e in canonical_registry()["entries"]]
        and m["registry_identity"] == canonical_registry()["identity"],
        "Changed profile identities",
    )
    historical = history(repo, BASE)
    require(m["historical_files"] == historical, "Historical inventory changed")
    for name, expected in historical.items():
        require(pin(repo / name) == expected, "Historical bytes changed: " + name)
    for collection, base in ((m["implementation_files"], repo), (m["artifacts"], root)):
        for name, expected in collection.items():
            path = base / name
            require(
                not Path(name).is_absolute()
                and ".." not in Path(name).parts
                and path.resolve().is_relative_to(base)
                and pin(path) == expected,
                "Changed file: " + name,
            )
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p.name not in ("receipt.json", "manifest.json")
    }
    # Workflow submanifests are ordinary artifacts and must not escape membership.
    expected_files = set(m["artifacts"])
    actual_all = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "receipt.json", root / "manifest.json")
    }
    require(
        actual_all == expected_files and actual <= actual_all,
        "Artifact membership changed",
    )
    independent = read(root / "independent-inputs.json")
    require(
        independent[
            "canonical_source_graph_type_charge_and_parameter_mapping_reconstructed"
        ]
        is True
        and independent["cases"] == [0, 2, 4, 6],
        "Missing independent reconstruction",
    )
    s = read(root / "audit/summary.json")
    require(
        s["identity"] == r["summary_identity"]
        and s["software_selection_acceptance_gate"] == 0
        and s["lammps_end_to_end_acceptance_gate"] == 1
        and s["source_assets"] == m["source_assets"],
        "Changed matrix acceptance",
    )
    require(set(s["cases"]) == {"0", "2", "4", "6"}, "Missing fixture")
    for dp, c in s["cases"].items():
        require(
            c["selected_profile"] == PROFILE
            and c["public_forcefield"] == "opls"
            and c["user_profile_argument"] is c["user_charge_method_argument"] is None,
            "Automatic OPLS selection missing",
        )
        require(
            c["acceptance_passed"] is False
            and c["MD_executed"] is False
            and c["NPT_executed"] is False,
            "Promoted blocked MD",
        )
        require(
            c["relocation_verified"] is True
            and c["fresh_process_reconstruction_verified"] is True
            and c["process_hash_seeds"] == [1, 99397],
            "Missing process relocation",
        )
        require(
            c["identical_physical_inputs"] is True
            and c["identical_numerical_long_range_settings"] is False,
            "Altered solver contract",
        )
        for kind in ("bond", "angle", "proper_torsion", "lj"):
            item = c["component_comparison"][kind]
            require(
                item["passed"] is True and item["bound"] == 3e-4,
                "Altered component bound",
            )
        for kind in ("total", "coulomb", "forces", "virial"):
            require(
                c["component_comparison"][kind]["passed"] is False, "Lost raw failure"
            )
        for engine in ("OpenMM", "LAMMPS"):
            require(c[engine + "_finite_difference"]["comparisons"] == 24, "Missing FD")
        require(
            c["LAMMPS_finite_difference"]["passed"] is False
            and c["LAMMPS_finite_difference"]["bound"] == 1e-4,
            "Promoted FD",
        )
        assignment = read(root / f"audit/dp{dp}/canonical-assignment.json")
        require(
            assignment["identity"] == c["assignment_identity"]
            and assignment["selected_profile"] == PROFILE
            and assignment["source_kind"] == "external_literature"
            and assignment["charge_method"] == "author_dft_npa"
            and not assignment["graph_charges"]["normalization"]
            and not assignment["graph_charges"]["redistribution"]
            and assignment["missing_conventions"],
            "Changed automatic assignment",
        )
    counts, nodes = Counter(), []
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
        and all(v == 1 for v in Counter(nodes).values()),
        "Incomplete or repeated suite coverage",
    )
    require(
        dict(counts) == r["full_suite_counts"]
        and counts["failure"] == counts["error"] == 0,
        "Changed suite outcomes",
    )
    require(
        all(v == 0 for v in r["check_returncodes"].values()), "Failed required check"
    )
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(historical),
        "full_suite_counts": dict(counts),
        "software_selection_acceptance_gate": 0,
        "lammps_end_to_end_acceptance_gate": 1,
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
