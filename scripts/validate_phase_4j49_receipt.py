"""Validate raw J49 controls without promoting failed author-source MD gates."""

import argparse
import json
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.external_pdms import LIMITS
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import digest
from scripts.validate_phase_4j47_receipt import history, pin, require

BASE = "9e2e5df9cfebf119c0c857c841a401fd98fcb42a"
BRANCH = "codex/phase-4j49-opls-pdms-electrostatic-closure"


def read(path):
    r = json.loads(Path(path).read_text())
    require(
        r["identity"] == digest({k: v for k, v in r.items() if k != "identity"}),
        "Changed identity: " + str(path),
    )
    return r


def validate(root, *, repo, expected_sha256=None):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    if expected_sha256 is not None:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "Independent receipt hash mismatch",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["schema"] == "island.j49-receipt.v1"
        and m["schema"] == "island.j49-manifest.v1",
        "Wrong schema",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong exact base/branch",
    )
    require(
        r["profile_identity"] == PROFILE and m["profile_identities"] == [PROFILE],
        "Profile changed",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established"
        and r["limitations"] == LIMITS,
        "Readiness promoted",
    )
    require(
        r["end_to_end_acceptance_gate"] == r["scientific_acceptance_gate"] == 1
        and r["bounded_diagnostic_acceptance_gate"] == 0,
        "Changed source/diagnostic authority",
    )
    require(
        not r["author_profile_MD_executed"]
        and r["diagnostic_MD_executed"]
        and not r["NPT_executed"],
        "Changed MD scope",
    )
    require(
        r["new_profiles"] == r["new_physical_parameters"] == [],
        "Physical model changed",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Changed manifest")
    old = history(repo, BASE)
    require(old == m["historical_files"], "Historical inventory changed")
    for name, sha in old.items():
        require(pin(repo / name) == sha, "Historical bytes changed: " + name)
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p.name not in ("receipt.json", "manifest.json")
    }
    # Nested workflow manifests are artifacts too.
    actual.update(
        str(p.relative_to(root))
        for p in root.rglob("manifest.json")
        if p != root / "manifest.json"
    )
    require(actual == set(m["artifacts"]), "Evidence membership changed")
    for mapping, base in ((m["artifacts"], root), (m["implementation_files"], repo)):
        for name, sha in mapping.items():
            p = base / name
            require(
                not Path(name).is_absolute()
                and ".." not in Path(name).parts
                and p.resolve().is_relative_to(base)
                and not p.is_symlink()
                and pin(p) == sha,
                "Changed file: " + name,
            )
    s = read(root / "audit/summary.json")
    require(
        s["identity"] == r["summary_identity"]
        and s["end_to_end_acceptance_gate"] == 1
        and s["bounded_diagnostic_acceptance_gate"] == 0,
        "Changed acceptance",
    )
    require(
        s["source_assets"] == m["source_assets"]
        and set(s["cases"]) == {"0", "2", "4", "6"},
        "Changed source or cases",
    )
    require(
        s["source_authority_unchanged"]
        and not s["corrections_applied_to_author_outputs"],
        "Silent correction",
    )
    for dp, c in s["cases"].items():
        require(
            c["selected_profile"] == PROFILE
            and c["public_forcefield"] == "opls"
            and c["user_profile_argument"] is c["user_charge_method_argument"] is None,
            "Not canonical automatic selection",
        )
        require(
            c["relocation"]["verified"]
            and c["relocation"]["process_hash_seeds"] == [1, 99397],
            "Missing fresh-process relocation",
        )
        require(
            c["OpenMM_FD_passed"]
            and not c["acceptance_passed"]
            and not c["MD_executed"],
            "Changed original source failure",
        )
        binding = read(root / f"audit/input-binding-dp{dp}.json")
        require(
            binding["identity"] == c["input_binding_identity"]
            and binding["selected_profile"] == PROFILE
            and binding["charge_origin"] == c["charge_origin"],
            "Changed bound inputs",
        )
        modes = json.loads((root / f"audit/components-dp{dp}.json").read_text())
        require(
            set(modes)
            == {
                "author",
                "no-table",
                "converged-pppm",
                "matched-ewald",
                "exact-erfc-ewald",
            },
            "Missing numerical isolation",
        )
        for name, v in modes.items():
            fd = v["finite_difference"]
            require(
                fd["steps_angstrom"] == [2e-4, 1e-4, 5e-5]
                and fd["bound"] == 1e-4
                and fd["comparisons_per_component"] == 36,
                "Changed FD coverage/bounds",
            )
            require(
                v["energy_bound"] == v["force_bound"] == 3e-4
                and v["diagnostic_only"]
                and not v["MD_authority"],
                "Control promoted",
            )
            maxima = {
                k: max(
                    row["components"][k]["absolute_error"] for row in fd["raw_records"]
                )
                for k in ("bonded", "lj", "coulomb", "total")
            }
            require(
                maxima == fd["max_error_by_component"]
                and fd["passed"] == all(v <= 1e-4 for v in maxima.values()),
                "Changed raw FD gate",
            )
        require(
            not modes["author"]["finite_difference"]["passed"]
            and modes["author"]["energy_error_kj_mol"]["coulomb"] > 3e-4,
            "Author failure lost",
        )
        exact = modes["exact-erfc-ewald"]
        require(
            exact["finite_difference"]["passed"]
            and all(v <= 3e-4 for v in exact["energy_error_kj_mol"].values())
            and all(v <= 3e-4 for v in exact["force_error_kj_mol_angstrom"].values()),
            "Failed diagnostic numerical gate",
        )
        w = s["diagnostic_workflows"][dp]
        directory = root / f"audit/workflow-dp{dp}"
        require(
            w["NVT"]["steps"] == 40
            and w["NVT"]["split_step"] == 20
            and w["NVT"]["dt_fs"] == 0.25
            and w["minimization"]["executed"]
            and w["relocated_separate_process_continuation"],
            "Missing bounded workflow",
        )
        require(
            w["checkpoint_bytes_equal"]
            and w["complete_RNG_state_equal"]
            and w["RNG_reconstructed_in_fresh_process"]
            and w["future_random_stream_probe_equal"]
            and w["full_state_dump_equal"]
            and w["bounded_diagnostic_workflow_passed"],
            "Failed exact continuation",
        )
        require(
            (directory / "whole.restart").read_bytes()
            == (directory / "resumed.restart").read_bytes()
            and pin(directory / "whole.restart")
            == w["whole_checkpoint_sha256"]
            == w["resumed_checkpoint_sha256"],
            "Native checkpoint changed",
        )
        require(
            (directory / "rng-initial.json").read_bytes()
            == (directory / "rng-continuation.json").read_bytes()
            and (directory / "whole.dump").read_bytes()
            == (directory / "resumed.dump").read_bytes(),
            "RNG or complete physical state changed",
        )
        require(
            not w["canonical_author_MD_authority"] and not w["NPT_executed"],
            "Promoted diagnostic MD",
        )
    inv = json.loads((root / "audit/kspace-inventory.json").read_text())
    for name in ("author", "exact_erfc_diagnostic"):
        require(
            "KSPACE" in inv[name]["help_text"]
            and inv[name]["units"] == "real"
            and "double" in inv[name]["help_text"],
            "Missing KSPACE/build precision",
        )
    require(
        inv["author"]["executable_sha256"]
        != inv["exact_erfc_diagnostic"]["executable_sha256"],
        "Hidden kernel replacement",
    )
    independent = read(root / "independent-inputs.json")
    require(
        independent["outside_production_parser"]
        and independent["native_workflows_reconstructed"] == ["0", "2", "4", "6"],
        "Missing independent reconstruction",
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
        "Suite coverage missing/repeated",
    )
    require(
        dict(counts) == r["full_suite_counts"]
        and counts["failure"] == counts["error"] == 0,
        "Test failure",
    )
    require(all(v == 0 for v in r["check_returncodes"].values()), "Failed verification")
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(old),
        "full_suite_counts": dict(counts),
        "bounded_diagnostic_acceptance_gate": 0,
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
