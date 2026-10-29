"""Validate J51 bytes and numerical evidence without authorizing unresolved MD."""

import argparse
import json
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.opls_registry.registry import ARCHIVE_SHA256, PAPER_SHA256
from island.forcefields.pdms_opls import PROFILE
from island.opls_conventions import (
    AUTHOR_ENGINE,
    EXACT_ENGINE,
    LAYERS,
    POLICIES,
    REQUIRED_SETTINGS,
)
from island.opls_resolution import (
    REMAINING_SOURCE_GAPS,
    SUPPLEMENT_SHA256,
    TEN_CONVENTIONS,
)
from scripts.phase_4j51 import BASE, BRANCH
from scripts.validate_phase_4j47_receipt import history, pin, require
from scripts.validate_phase_4j50_receipt import numerical, read


def validate(directory, *, repo, expected_sha256=None):
    root, repo = Path(directory).resolve(), Path(repo).resolve()
    if expected_sha256 is not None:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "Independent receipt hash mismatch",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["schema"] == "island.j51-receipt.v1"
        and m["schema"] == "island.j51-manifest.v1",
        "Wrong schema",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong base or branch",
    )
    require(
        r["profile_identity"] == PROFILE and m["profile_identities"] == [PROFILE],
        "Changed physical profile",
    )
    require(
        m["engine_executable_sha256"]
        == {
            "source_replay_v1": AUTHOR_ENGINE,
            "exact_erfc_diagnostic_v1": EXACT_ENGINE,
        },
        "Changed pinned executables",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established",
        "Changed readiness",
    )
    require(
        not r["MD_executed"]
        and not r["NPT_executed"]
        and not r["charges_parameters_or_tolerances_changed"],
        "Unauthorized execution or model change",
    )
    require(
        r["research_opt_in_acceptance_gate"]
        == r["production_candidate_acceptance_gate"]
        == 1
        and r["diagnostic_acceptance_gate"] == 0,
        "Promoted source authority",
    )
    require(
        r["ten_conventions"] == list(TEN_CONVENTIONS)
        and r["converter_resolved_from_ten"] == ["lj_tail_correction"]
        and r["remaining_source_gaps"] == list(REMAINING_SOURCE_GAPS),
        "Changed closure outcome",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Changed manifest")
    old = history(repo, BASE)
    require(m["historical_files"] == old, "Changed history membership")
    for name, expected in old.items():
        require(pin(repo / name) == expected, "Changed historical bytes: " + name)
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "manifest.json", root / "receipt.json")
    }
    require(actual == set(m["artifacts"]), "Changed evidence membership")
    for pins, base in ((m["artifacts"], root), (m["implementation_files"], repo)):
        for name, expected in pins.items():
            path = base / name
            require(
                not Path(name).is_absolute()
                and ".." not in Path(name).parts
                and not path.is_symlink()
                and path.resolve().is_relative_to(base)
                and pin(path) == expected,
                "Changed evidence/code: " + name,
            )
    summary = read(root / "audit/summary.json")
    require(
        summary["identity"] == r["summary_identity"]
        and summary["diagnostic_acceptance_gate"] == 0
        and summary["research_opt_in_acceptance_gate"]
        == summary["production_candidate_acceptance_gate"]
        == 1
        and not summary["MD_executed"]
        and not summary["source_authority_promoted"],
        "Changed acceptance summary",
    )
    require(
        summary["source_assets"]
        == {"jp4c08471.pdf": PAPER_SHA256, "jp4c08471_si_001.pdf": SUPPLEMENT_SHA256}
        and summary["source_archive_sha256"] == ARCHIVE_SHA256,
        "Changed primary sources",
    )
    require(set(summary["cases"]) == {"0", "2", "4", "6"}, "Changed bounded fixtures")
    for dp, case in summary["cases"].items():
        assignment = None
        for policy in POLICIES:
            directory = root / f"audit/dp{dp}" / policy
            p = read(directory / "resolution.json")
            prior = read(directory / "conventions.json")
            result = case["policies"][policy]
            require(
                p["family"] == "opls"
                and p["selected_profile"] == PROFILE
                and p["identity"] == result["resolution_identity"]
                and p["four_layer_identity"] == prior["identity"],
                "Lost canonical resolution",
            )
            require(
                all(p[layer] == prior[layer] for layer in LAYERS),
                "Changed historical authority layers",
            )
            require(
                p["ten_conventions"] == list(TEN_CONVENTIONS)
                and [v["setting"] for v in p["ten_convention_resolutions"]]
                == list(TEN_CONVENTIONS)
                and {v["setting"] for v in p["rows"]} == set(REQUIRED_SETTINGS)
                and {v["setting"] for v in p["unresolved_conventions"]}
                == set(REMAINING_SOURCE_GAPS),
                "Lost explicit gap table",
            )
            for row in p["rows"]:
                if row["setting"] in REMAINING_SOURCE_GAPS:
                    require(
                        row["resolution"] == "unresolved"
                        and row["value"] is None
                        and row["missing_evidence"]
                        and not row["source_authoritative"],
                        "Invented source proof",
                    )
                if row["setting"] == "lj_tail_correction":
                    require(
                        row["resolution"] == "converter_resolved"
                        and row["value"] == "no"
                        and not row["source_authoritative"]
                        and any(
                            "None of the" in c["raw_text"]
                            for c in row["candidates"][0]["evidence"]
                        ),
                        "Tail default substituted for proof",
                    )
                require(
                    row["backend"]["resolution"] == "backend_resolved"
                    and not row["backend"]["source_authoritative"],
                    "Promoted backend authority",
                )
            require(
                not p["research_MD_available"]
                and not p["conflicts"]
                and p["source_scope_review"]["status"] == "unresolved",
                "Invented execution permission",
            )
            if assignment is None:
                assignment = prior["assignment_identity"]
            require(
                prior["assignment_identity"]
                == assignment
                == result["assignment_identity"],
                "Changed physical assignment",
            )
            require(
                result["deterministic_export"]
                and result["relocation"]["verified"]
                and [v["hash_seed"] for v in result["relocation"]["fresh_processes"]]
                == [1, 99397],
                "Missing fresh relocation",
            )
            c = numerical(directory)
            require(
                c["energy_error_kj_mol"] == result["energy_error_kj_mol"]
                and c["force_error_kj_mol_angstrom"]
                == result["force_error_kj_mol_angstrom"]
                and c["finite_difference"]["passed"] == result["FD_passed"],
                "Changed component results",
            )
            if policy == POLICIES[0]:
                require(
                    c["energy_error_kj_mol"]["coulomb"] > 3e-4
                    and not result["FD_passed"],
                    "Lost source-command discrepancy",
                )
            else:
                require(
                    all(v <= 3e-4 for v in c["energy_error_kj_mol"].values())
                    and all(
                        v <= 3e-4 for v in c["force_error_kj_mol_angstrom"].values()
                    )
                    and result["FD_passed"],
                    "Failed unchanged diagnostic bounds",
                )
            d = json.loads((directory / "detailed-components.json").read_text())
            require(
                d["virial_acceptance_bound"] is None
                and not d["pressure_validated"]
                and d["OpenMM"]["improper"]["status"] == "not_applicable"
                and d["special_pairs_independently_reconstructed"]["pairs"],
                "Missing components or promoted pressure",
            )
        require(
            len(case["blocked_execution"]) == 18
            and all(
                not v["engine_started"] and not v["outputs_created"]
                for v in case["blocked_execution"]
            ),
            "Blocked path reached MD",
        )
    independent = json.loads((root / "independent-reconstruction.json").read_text())
    require(
        independent["outside_production_parser"]
        and independent["numerical"]["raw_FD_samples_reconstructed"] == 288
        and independent["no_tail_rule_reconstructed"] == 8
        and independent["unresolved_rows_preserved"] == 72
        and len(independent["source_reviews"]) == 6
        and len(independent["exact_kernel_independently_reconstructed"]) == 2
        and independent["kernel_build_receipt_sha256"]
        == pin(repo / "docs/evidence/phase_4j49/kernel-build.json")
        and independent["prior_comparisons_read"],
        "Incomplete independent reconstruction",
    )
    require(
        len(independent["independent_special_pairs_and_pressure"]) == 8
        and all(
            v["native_pressure_proxy_reconstructed"]
            and not v["pressure_scientifically_validated"]
            for v in independent["independent_special_pairs_and_pressure"].values()
        ),
        "Missing independent special pairs or promoted pressure",
    )
    controls = json.loads((root / "validation/validator-controls.json").read_text())
    require(
        len(controls["controls"]) == 7
        and all(v["rejected"] for v in controls["controls"])
        and controls["relocated_receipt"]["verified"]
        and controls["relocated_receipt"]["separate_process"],
        "Missing validator rejection or relocation controls",
    )
    expected_sources = {
        v["source_file"]: v["source_sha256"] for v in independent["source_reviews"]
    }
    expected_sources["PDMS_benchmark-main.zip"] = ARCHIVE_SHA256
    require(m["source_assets"] == expected_sources, "Changed source inventory")
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
        "Missing or duplicate tests",
    )
    require(
        dict(counts) == r["full_suite_counts"]
        and counts["failure"] == counts["error"] == 0,
        "Failed suite",
    )
    require(all(v == 0 for v in r["check_returncodes"].values()), "Failed verification")
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(old),
        "full_suite_counts": dict(counts),
        "ten_conventions": 10,
        "converter_closed": 1,
        "source_gaps_remaining": 9,
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
