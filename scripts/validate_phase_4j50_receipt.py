"""Validate four-layer J50 evidence without granting source or MD authority."""

import argparse
import json
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.external_pdms import LIMITS
from island.forcefields.opls_registry.registry import ARCHIVE_SHA256, PAPER_SHA256
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import digest
from island.opls_conventions import CLASSIFICATIONS, LAYERS, POLICIES, REQUIRED_SETTINGS
from scripts.phase_4j50 import BASE
from scripts.validate_phase_4j47_receipt import history, pin, require

BRANCH = "codex/phase-4j50-opls-pdms-lammps-convention-reconciliation"


def read(path):
    value = json.loads(Path(path).read_text())
    require(
        value["identity"]
        == digest({k: v for k, v in value.items() if k != "identity"}),
        "Changed identity: " + str(path),
    )
    return value


def numerical(directory):
    c = json.loads((directory / "components.json").read_text())
    require(
        c["diagnostic_only"] and not c["MD_authority"], "Promoted numerical comparison"
    )
    require(c["energy_bound"] == c["force_bound"] == 3e-4, "Changed limits")
    values = {}
    for line in (
        (directory / "lammps/comparison/engine.screen").read_text().splitlines()
    ):
        if line.startswith("J49FD ") and "$(" not in line:
            w = line.split()
            e = [float(v) * 4.184 for v in w[2:]]
            values[int(w[1])] = {
                "bonded": (e[0] + e[1]) + e[2],
                "lj": e[3],
                "coulomb": e[4] + e[5],
                "total": e[6],
            }
    fd = c["finite_difference"]
    require(
        fd["bound"] == 1e-4
        and fd["steps_angstrom"] == [2e-4, 1e-4, 5e-5]
        and fd["comparisons_per_component"] == 36
        and len(values) == 72,
        "Changed FD scope",
    )
    maxima = dict.fromkeys(("bonded", "lj", "coulomb", "total"), 0.0)
    for i, r in enumerate(fd["raw_records"]):
        for k, v in r["components"].items():
            positive, negative = values[2 * i][k], values[2 * i + 1][k]
            force = -(positive - negative) / (2 * r["h_angstrom"])
            error = abs(force - v["analytical_force"])
            require(
                positive == v["positive_energy_kj_mol"]
                and negative == v["negative_energy_kj_mol"]
                and force == v["central_difference_force"]
                and error == v["absolute_error"],
                "Changed raw FD data",
            )
            maxima[k] = max(maxima[k], error)
    require(
        maxima == fd["max_error_by_component"]
        and fd["passed"] == all(v <= 1e-4 for v in maxima.values()),
        "Changed FD gate",
    )
    return c


def validate(root, *, repo, expected_sha256=None):
    root, repo = Path(root).resolve(), Path(repo).resolve()
    if expected_sha256 is not None:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "Independent receipt hash mismatch",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["schema"] == "island.j50-receipt.v1"
        and m["schema"] == "island.j50-manifest.v1",
        "Wrong schema",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong base/branch",
    )
    require(
        r["profile_identity"] == PROFILE and m["profile_identities"] == [PROFILE],
        "Wrong profile",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established"
        and r["limitations"] == LIMITS,
        "Changed scientific readiness",
    )
    require(
        not r["MD_executed"]
        and not r["NPT_executed"]
        and not r["new_physical_profiles"]
        and not r["charges_or_parameters_changed"],
        "Promoted execution or changed model",
    )
    require(
        r["research_opt_in_acceptance_gate"]
        == r["production_candidate_acceptance_gate"]
        == 1
        and r["diagnostic_acceptance_gate"] == 0,
        "Promoted source acceptance",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Changed manifest")
    old = history(repo, BASE)
    require(m["historical_files"] == old, "Changed historical inventory")
    for name, expected in old.items():
        require(pin(repo / name) == expected, "Changed historical bytes: " + name)
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "receipt.json", root / "manifest.json")
    }
    require(actual == set(m["artifacts"]), "Changed evidence membership")
    for collection, base in ((m["artifacts"], root), (m["implementation_files"], repo)):
        for name, expected in collection.items():
            path = base / name
            require(
                not Path(name).is_absolute()
                and ".." not in Path(name).parts
                and not path.is_symlink()
                and path.resolve().is_relative_to(base)
                and pin(path) == expected,
                "Changed evidence/implementation: " + name,
            )
    require(
        m["source_assets"] == {"archive": ARCHIVE_SHA256, "paper": PAPER_SHA256},
        "Changed assets",
    )
    summary = read(root / "audit/summary.json")
    require(
        summary["identity"] == r["summary_identity"]
        and set(summary["cases"]) == {"0", "2", "4", "6"},
        "Changed cases",
    )
    require(
        summary["diagnostic_acceptance_gate"] == 0
        and summary["research_opt_in_acceptance_gate"]
        == summary["production_candidate_acceptance_gate"]
        == 1
        and not summary["source_authority_promoted"]
        and not summary["MD_executed"],
        "Changed authority",
    )
    for dp, case in summary["cases"].items():
        original_assignment = None
        for policy in POLICIES:
            directory = root / f"audit/dp{dp}" / policy
            p = read(directory / "conventions.json")
            result = case["policies"][policy]
            require(
                p["identity"] == result["convention_identity"]
                and p["selected_profile"] == PROFILE
                and p["family"] == "opls"
                and p["execution_mode"] == "diagnostic",
                "Not canonical selection",
            )
            for layer in LAYERS:
                require(
                    p[layer]["identity"]
                    == digest({k: v for k, v in p[layer].items() if k != "identity"})
                    == result["layer_identities"][layer],
                    "Changed layer identity",
                )
            if original_assignment is None:
                original_assignment = p["assignment_identity"]
            require(
                p["assignment_identity"] == original_assignment,
                "Numerical policy changed physical assignment",
            )
            source, backend = p[LAYERS[0]], p[LAYERS[2]]
            require(
                not p[LAYERS[3]]["current_graph_numerically_verified_by_preparation"]
                and p[LAYERS[3]]["prior_verified_cases_only"]
                == ["HMDS", "capped_DP2", "capped_DP4", "capped_DP6"],
                "Preparation promoted numerical qualification",
            )
            require(
                set(source["settings"])
                == set(backend["settings"])
                == set(REQUIRED_SETTINGS),
                "Missing setting",
            )
            for name, setting in source["settings"].items():
                require(
                    setting["classification"] in CLASSIFICATIONS,
                    "Unknown authority class",
                )
                if setting["classification"] == "unresolved":
                    require(
                        setting["value"] is None
                        and setting["missing_evidence"]
                        and backend["settings"][name]["classification"]
                        == "selected_by_backend_policy"
                        and not backend["settings"][name]["source_authority"],
                        "Promoted backend setting",
                    )
            require(
                p["unresolved_settings"] and not p["research_MD_available"],
                "Hidden source gaps",
            )
            require(
                result["deterministic_export_verified"]
                and result["relocation"]["verified"]
                and [x["hash_seed"] for x in result["relocation"]["fresh_processes"]]
                == [1, 99397],
                "Missing relocation",
            )
            c = numerical(directory)
            require(
                result["energy_error_kj_mol"] == c["energy_error_kj_mol"]
                and result["force_error_kj_mol_angstrom"]
                == c["force_error_kj_mol_angstrom"]
                and result["finite_difference_passed"]
                == c["finite_difference"]["passed"],
                "Changed component report",
            )
            if policy == POLICIES[0]:
                require(
                    c["energy_error_kj_mol"]["coulomb"] > 3e-4
                    and not c["finite_difference"]["passed"],
                    "Lost author-setting failure",
                )
            else:
                require(
                    all(v <= 3e-4 for v in c["energy_error_kj_mol"].values())
                    and all(
                        v <= 3e-4 for v in c["force_error_kj_mol_angstrom"].values()
                    )
                    and c["finite_difference"]["passed"],
                    "Failed diagnostic control",
                )
        require(
            len(case["blocked_modes"]) == 6
            and all(
                not x["MD_executed"]
                and not x["engine_started"]
                and not x["checkpoint_created"]
                and not x["trajectory_created"]
                for x in case["blocked_modes"]
            ),
            "MD boundary bypassed",
        )
    independent = json.loads((root / "independent-rows.json").read_text())
    require(
        independent["outside_production_parser"]
        and independent["raw_FD_samples_reconstructed"] == 288
        and independent["source_citations_reconstructed"] > 0,
        "Missing independent reconstruction",
    )
    require(
        independent["PDF_citations_reconstructed"] == 8
        and independent["paper_sha256"] == PAPER_SHA256,
        "Missing independent paper boundary evidence",
    )
    counts, nodes = Counter(), []
    for file in r["junit_files"]:
        for case in ET.parse(root / file).iter("testcase"):
            name = case.attrib["classname"]
            require(name.startswith("tests."), "Wrong test class")
            nodes.append(name.replace(".", "/") + ".py::" + case.attrib["name"])
            counts["tests"] += 1
            for flag in ("failure", "error", "skipped"):
                counts[flag] += int(case.find(flag) is not None)
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected)
        and all(v == 1 for v in Counter(nodes).values()),
        "Missing/repeated test coverage",
    )
    require(
        dict(counts) == r["full_suite_counts"]
        and counts["failure"] == counts["error"] == 0,
        "Test failures",
    )
    require(
        all(v == 0 for v in r["check_returncodes"].values()), "Failed required check"
    )
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(old),
        "full_suite_counts": dict(counts),
        "diagnostic_acceptance_gate": 0,
        "research_opt_in_acceptance_gate": 1,
        "production_candidate_acceptance_gate": 1,
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
