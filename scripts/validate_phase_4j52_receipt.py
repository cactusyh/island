"""Validate J52 scientific limits, independent raw numbers and every historical byte."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.forcefields.pdms_opls import PROFILE
from island.opls_conventions import (
    AUTHOR_ENGINE,
    EXACT_ENGINE,
    LAYERS,
    POLICIES,
    signed,
)
from island.opls_final_resolution import NINE_CONVENTIONS, corpus_audit
from scripts.phase_4j52 import BASE, BRANCH
from scripts.seal_phase_4j51_receipt import junit
from scripts.seal_phase_4j52_receipt import IMPLEMENTATION
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
        r["schema"] == "island.j52-receipt.v1"
        and m["schema"] == "island.j52-manifest.v1",
        "Wrong schema",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong exact base or branch",
    )
    require(
        r["profile_identity"] == PROFILE
        and m["profile_identities"] == [PROFILE]
        and r["public_family"] == "opls",
        "Wrong family/profile",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established"
        and not r["MD_executed"]
        and not r["NPT_executed"]
        and not r["charges_parameters_or_tolerances_changed"],
        "Unauthorized readiness or MD claim",
    )
    require(
        r["diagnostic_acceptance_gate"] == 0
        and r["research_opt_in_acceptance_gate"]
        == r["production_candidate_acceptance_gate"]
        == 1,
        "Promoted source authority",
    )
    require(
        r["nine_conventions"]
        == r["permanently_unresolved_in_pinned_corpus"]
        == list(NINE_CONVENTIONS),
        "Lost exact nine-item inventory",
    )
    require(
        r["checkpoint_RNG_continuation"]
        == "not run: source prerequisites remain blocked",
        "Unauthorized continuation claim",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Changed manifest")
    old = history(repo, BASE)
    require(m["historical_files"] == old, "Changed historical membership")
    require(all(pin(repo / n) == v for n, v in old.items()), "Changed historical bytes")
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "receipt.json", root / "manifest.json")
    }
    require(actual == set(m["artifacts"]), "Changed evidence membership")
    require(
        set(m["implementation_files"]) == set(IMPLEMENTATION),
        "Changed implementation membership",
    )
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
    audit = corpus_audit()
    require(
        read(root / "audit/supplied-corpus.json") == audit
        and r["source_audit_identity"]
        == m["source_audit_identity"]
        == audit["identity"],
        "Changed exhaustive corpus",
    )
    require(
        m["source_assets"]
        == {a["source_file"]: a["source_sha256"] for a in audit["assets"]}
        and m["source_archive_sha256"] == audit["source_archive_sha256"],
        "Changed source pins",
    )
    require(
        m["engine_executable_sha256"]
        == {POLICIES[0]: AUTHOR_ENGINE, POLICIES[1]: EXACT_ENGINE},
        "Changed engines",
    )
    summary = read(root / "audit/summary.json")
    require(
        summary["identity"] == r["summary_identity"]
        and summary["diagnostic_acceptance_gate"] == 0
        and summary["research_opt_in_acceptance_gate"]
        == summary["production_candidate_acceptance_gate"]
        == 1
        and not summary["MD_executed"]
        and not summary["NPT_executed"]
        and not summary["source_authority_promoted"],
        "Changed acceptance summary",
    )
    require(set(summary["cases"]) == {"0", "2", "4", "6"}, "Changed fixture matrix")
    table = read(root / "audit/final-nine-conventions.json")
    require(
        [v["setting"] for v in table["rows"]] == list(NINE_CONVENTIONS)
        and table["source_audit_identity"] == audit["identity"],
        "Lost final nine table",
    )
    for row, expected in zip(table["rows"], audit["final_resolutions"]):
        raw = [
            {
                **{
                    k: a[k]
                    for k in (
                        "source_file",
                        "source_sha256",
                        "doi",
                        "citation",
                        "scope",
                    )
                },
                **cite,
            }
            for a in audit["assets"]
            for cite in a["matches"][row["setting"]]
        ]
        require(
            row
            == signed(
                {
                    **{k: v for k, v in expected.items() if k != "identity"},
                    "corpus_result_identity": expected["identity"],
                    "raw_evidence": raw,
                }
            ),
            "Changed final raw evidence",
        )
    blockers = read(root / "audit/permanent-blockers.json")
    require(
        blockers["rows"] == audit["final_resolutions"]
        and not blockers["MD_authorized"],
        "Changed permanent blockers",
    )
    for dp, case in summary["cases"].items():
        physical = None
        for policy in POLICIES:
            directory = root / f"audit/dp{dp}" / policy
            final = read(directory / "final-resolution.json")
            historical = read(directory / "resolution.json")
            prior = read(directory / "conventions.json")
            result = case["policies"][policy]
            require(
                final["family"] == "opls"
                and final["selected_profile"] == PROFILE
                and final["identity"] == result["final_resolution_identity"]
                and final["J51_resolution_identity"] == historical["identity"]
                and final["source_audit_identity"] == audit["identity"],
                "Changed canonical final identity",
            )
            require(
                all(
                    final[layer] == historical[layer] == prior[layer]
                    for layer in LAYERS
                ),
                "Changed four authority layers",
            )
            require(
                final["nine_conventions"] == list(NINE_CONVENTIONS)
                and [v["setting"] for v in final["final_resolutions"]]
                == list(NINE_CONVENTIONS)
                and final["source_assets"] == historical["source_assets"]
                and final["charge_origin"] == historical["charge_origin"]
                and final["domain"] == historical["domain"]
                and final["unresolved_conventions"]
                == final["permanent_blockers"]
                == final["final_resolutions"]
                and not final["research_MD_available"]
                and not final["default_MD_enabled"]
                and not final["production_candidate_available"],
                "Promoted final authority",
            )
            for row, source in zip(
                final["final_resolutions"], audit["final_resolutions"]
            ):
                old_row = next(
                    v for v in historical["rows"] if v["setting"] == source["setting"]
                )
                expected = signed(
                    {
                        **{k: v for k, v in source.items() if k != "identity"},
                        "historical_resolution_identity": old_row["identity"],
                        "backend": old_row["backend"],
                        "units": old_row["units"],
                        "conversion_implementation": old_row[
                            "conversion_implementation"
                        ],
                        "backend_conversion_formula": old_row["conversion_formula"],
                    }
                )
                require(row == expected, "Changed final convention or backend")
            # Exact source-bound physical and numerical inputs must match J51.
            for name in (
                "resolution.json",
                "conventions.json",
                "canonical-assignment.json",
                "system.json",
                "lammps/system.data",
                "lammps/system.in",
                "lammps/coefficients.in",
            ):
                require(
                    pin(directory / name)
                    == pin(
                        repo / f"docs/evidence/phase_4j51/audit/dp{dp}" / policy / name
                    ),
                    "Changed J51 fixture contract: " + name,
                )
            require(
                physical is None or prior["assignment_identity"] == physical,
                "Changed assignment by numerical policy",
            )
            physical = prior["assignment_identity"]
            require(
                result["assignment_identity"] == physical
                and result["deterministic_export"]
                and result["relocation"]["verified"]
                and [v["hash_seed"] for v in result["relocation"]["fresh_processes"]]
                == [1, 99397]
                and all(
                    v["returncode"] == 0 and v["identity"] == final["identity"]
                    for v in result["relocation"]["fresh_processes"]
                ),
                "Missing exact relocation reconstruction",
            )
            c = numerical(directory)
            require(
                c["energy_error_kj_mol"] == result["energy_error_kj_mol"]
                and c["force_error_kj_mol_angstrom"]
                == result["force_error_kj_mol_angstrom"]
                and c["finite_difference"]["passed"] == result["FD_passed"],
                "Changed raw component result",
            )
            if policy == POLICIES[0]:
                require(
                    c["energy_error_kj_mol"]["coulomb"] > 3e-4
                    and not result["FD_passed"],
                    "Lost source replay discrepancy",
                )
            else:
                require(
                    all(v <= 3e-4 for v in c["energy_error_kj_mol"].values())
                    and all(
                        v <= 3e-4 for v in c["force_error_kj_mol_angstrom"].values()
                    )
                    and result["FD_passed"],
                    "Failed unchanged diagnostic limits",
                )
            detail = json.loads((directory / "detailed-components.json").read_text())
            require(
                detail["virial_acceptance_bound"] is None
                and not detail["pressure_validated"]
                and detail["special_pairs_independently_reconstructed"]["pairs"],
                "Promoted pressure or missing partitions",
            )
        require(
            len(case["blocked_execution"]) == 18
            and all(
                not v["engine_started"] and not v["outputs_created"]
                for v in case["blocked_execution"]
            ),
            "Blocked MD executed",
        )
    independent = read(root / "independent-reconstruction.json")
    require(
        independent["outside_production_parser"]
        and not independent["MD_executed"]
        and not independent["NPT_executed"]
        and independent["assets_independently_reconstructed"] == 64
        and independent["source_audit_identity"] == audit["identity"]
        and independent["numerical"]["raw_FD_samples_reconstructed"] == 288
        and independent["exact_kernel_edits_reconstructed"] == 2
        and len(independent["independent_graph_mapping"])
        == len(independent["independent_special_pairs_and_pressure"])
        == 8
        and all(
            v["source_rows"] == 17
            and v["graph_type_mass_charge_and_parameter_rows_reconstructed"]
            for v in independent["independent_graph_mapping"].values()
        ),
        "Incomplete independent reconstruction",
    )
    counts, nodes = junit(root, r["junit_files"])
    expected = json.loads((root / "collection.json").read_text())["nodes"]
    base_nodes = json.loads(
        (repo / "docs/evidence/phase_4j51/collection.json").read_text()
    )["nodes"]
    focused, focused_nodes = junit(root, ["validation/focused-final.xml"])
    regressions, _ = junit(root, ["validation/j49-j51-regression.xml"])
    require(
        Counter(expected) == Counter(base_nodes + focused_nodes)
        and all(n.startswith("tests/test_phase_4j52.py::") for n in focused_nodes),
        "Lost historical full-suite nodes",
    )
    require(
        r["focused_tests_passed"] == focused["tests"] - focused["skipped"]
        and r["J49_J51_regressions_passed"]
        == regressions["tests"] - regressions["skipped"],
        "Changed focused or regression counts",
    )
    require(
        Counter(nodes) == Counter(expected)
        and all(n == 1 for n in Counter(nodes).values())
        and counts == r["full_suite_counts"],
        "Incomplete test accounting",
    )
    require(
        {
            "ruff",
            "compileall",
            "diff",
            "pip-phase4e2",
            "pip-phase4g1",
            "pip-phase4j32",
            "pip-phase4j36",
            "pip-md_rg",
            "J49_receipt",
            "J50_receipt",
            "J51_receipt",
            "historical_hashes",
            "independent_reconstruction",
            "fresh_process_relocation",
        }.issubset(r["check_returncodes"])
        and all(v == 0 for v in r["check_returncodes"].values()),
        "Failed required verification",
    )
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(old),
        "full_suite_counts": counts,
        "unresolved_permanently": list(NINE_CONVENTIONS),
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
