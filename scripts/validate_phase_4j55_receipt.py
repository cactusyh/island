"""Validate J55 independent claims, raw runtime evidence, history and test coverage."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.charge_references.records import unpack
from island.forcefield_policy import policy_identity
from island.forcefield_policy.claims import assess
from island.forcefields.opls_registry.registry import NATIVE
from island.forcefields.pdms_opls import PROFILE
from scripts.phase_4j55 import BASE, BRANCH
from scripts.seal_phase_4j51_receipt import junit
from scripts.seal_phase_4j55_receipt import IMPLEMENTATION, SUITES
from scripts.validate_phase_4j47_receipt import history, pin, require
from scripts.validate_phase_4j50_receipt import read


def validate(directory, *, repo, expected_sha256=None):
    root, repo = Path(directory).resolve(), Path(repo).resolve()
    if expected_sha256:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "Receipt independent pin changed",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["schema"] == "island.j55-receipt.v1"
        and m["schema"] == "island.j55-manifest.v1",
        "Wrong schema",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong exact base/branch",
    )
    require(
        r["profile_identities"] == m["profile_identities"] == [NATIVE, PROFILE]
        and r["new_profiles"] == m["new_profiles"] == [],
        "Changed profile family or authority",
    )
    require(
        r["parameter_completeness_independent_of_equivalence"] is True
        and r["ordinary_forcefield"] == "opls"
        and r["ordinary_profile_argument"] is None
        and r["ordinary_charge_method_argument"] is None,
        "Lost family-only selection",
    )
    require(
        r["cross_engine_equivalence"] == "not_claimed"
        and r["scientific_validation"]
        == r["production_readiness"]
        == r["simulation_readiness"]
        == "not_established"
        and all(
            r[k] is False
            for k in (
                "production_validated",
                "MD_executed",
                "NPT_executed",
                "charges_parameters_or_tolerances_changed",
                "historical_claims_rewritten",
            )
        ),
        "Promoted readiness or execution",
    )
    require(
        r["bounded_PDMS_parameter_completeness"] == "incomplete"
        and r["bounded_PDMS_missing_model_conventions"]
        == ["dielectric", "periodic_image_policy"]
        and r["native_finite_control_parameter_completeness"] == "complete",
        "Changed parameter claim",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Manifest changed")
    old = history(repo, BASE)
    require(
        m["historical_files"] == old and r["historical_files_preserved"] == len(old),
        "Changed history membership",
    )
    require(
        all(pin(repo / n) == value for n, value in old.items()),
        "Historical bytes changed",
    )
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "receipt.json", root / "manifest.json")
    }
    require(set(m["artifacts"]) == actual, "Evidence membership changed")
    require(
        set(m["implementation_files"]) == set(IMPLEMENTATION),
        "Implementation membership changed",
    )
    for pins, base in ((m["artifacts"], root), (m["implementation_files"], repo)):
        for name, expected in pins.items():
            p = base / name
            require(
                not Path(name).is_absolute()
                and ".." not in Path(name).parts
                and not p.is_symlink()
                and p.resolve().is_relative_to(base)
                and pin(p) == expected,
                "Changed evidence/code: " + name,
            )
    summary = read(root / "audit/summary.json")
    require(
        summary["identity"] == r["summary_identity"]
        and summary["base_commit"] == BASE
        and summary["source_hashes"] == m["source_hashes"]
        and summary["profile_identity"] == PROFILE
        and summary["new_profiles"] == []
        and summary["cross_engine_diagnostics_requested"] is False
        and summary["MD_executed"] is False
        and summary["production_validated"] is False
        and summary["simulation_readiness"] == "not_established",
        "Changed audit",
    )
    require(
        [c["case"] for c in summary["cases"]]
        == ["native-ethane", "dp0", "dp2", "dp4", "dp6"],
        "Lost fixture",
    )
    for c in summary["cases"]:
        case = root / "audit" / c["case"]
        p = read(case / "report.json")
        require(
            p["selected_profile"] == p["source_assignment"]["selected_profile"]
            and p["source_identity"] == p["source_assignment"]["source_identity"]
            and p["charge_record"] == p["source_assignment"]["graph_charges"]
            and p["assignment_identity"] == p["source_assignment"]["identity"],
            "Changed source/profile/charge binding",
        )
        require(
            p["identity"] == c["report_identity"]
            and p["policy_identity"] == policy_identity(),
            "Changed report/policy",
        )
        require(
            p["cross_engine_equivalence"] == "not_claimed"
            and p["cross_engine_diagnostics"] == []
            and p["scientific_validation"]
            == p["production_readiness"]
            == "not_established"
            and p["production_validated"] is False
            and p["MD_authorized"] is False,
            "Promoted fixture claim",
        )
        inv = p["interaction_inventory"]
        from island.forcefields.selection.inventory import requirements
        from island.pdms import restore_pdms_system
        from island.workflows.bundle import system_from

        physical = unpack((case / "bundle/system.json").read_text())
        system = (
            restore_pdms_system(physical)
            if physical["box"] is not None
            else system_from(physical)
        )
        exact = assess(
            requirements(system, "opls"),
            [row["assignment"] for row in inv if row["assignment"] is not None],
            p["conventions"],
            typing_complete=True,
            charge_complete=True,
            source_valid=True,
        )
        require(
            all(p[k] == v for k, v in exact.items()), "Changed completeness assessment"
        )
        require(
            p["typing_status"] == p["charge_status"] == "complete"
            and not p["missing_terms"],
            "Lost typing/charge/row coverage",
        )
        require(
            p["parameter_completeness"] == c["parameter_completeness"]
            and p["missing_conventions"] == c["missing_conventions"],
            "Changed fixture parameter status",
        )
        for name in (
            "openmm_build_status",
            "openmm_runtime_status",
            "lammps_export_status",
            "lammps_runtime_status",
        ):
            require(p[name] == c[name], "Changed independent backend status")
        require(
            p["openmm_build_status"] == p["openmm_runtime_status"] == "passed"
            and p["openmm_static"]["energy"]
            == p["openmm_static"]["result"]["potential_energy"],
            "Lost OpenMM runtime result",
        )
        if c["case"].startswith("dp"):
            dp = c["case"][2:]
            previous = read(
                repo
                / f"docs/evidence/phase_4j52/audit/dp{dp}/source_replay_v1/canonical-assignment.json"
            )
            require(
                p["source_assignment"] == previous,
                "Changed historical graph/source/types/charges",
            )
            require(
                p["parameter_completeness"] == "incomplete"
                and p["missing_conventions"] == ["dielectric", "periodic_image_policy"]
                and len(p["source_convention_gaps"]) == 9,
                "Erased physical or historical gaps",
            )
            require(
                p["lammps_export_status"] == p["lammps_runtime_status"] == "passed",
                "Lost LAMMPS runtime",
            )
            for n, value in p["lammps_export"]["files"].items():
                require(
                    pin(case / "lammps" / n) == value, "Changed deterministic export"
                )
            raw = json.loads((case / "static/result.json").read_text())
            require(
                raw["observed"] == p["lammps_static"]["observed"]
                and raw["MD_executed"] is False
                and raw["engine"]["executable_sha256"]
                == "6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191",
                "Changed pinned runtime evidence",
            )
        else:
            require(
                p["selected_profile"] == NATIVE
                and p["parameter_completeness"] == "complete",
                "Lost complete control",
            )
        require(
            pin(case / "bundle/manifest.json") == c["bundle_manifest_sha256"],
            "Changed workflow bundle",
        )
        bundle = read(case / "bundle/manifest.json")
        require(
            set(bundle["files"]) == {"system.json", "report.json"},
            "Changed workflow membership",
        )
        require(
            all(pin(case / "bundle" / n) == v for n, v in bundle["files"].items()),
            "Changed workflow file",
        )
        require(read(case / "bundle/report.json") == p, "Changed bundled report")
        from island.forcefields.selection.inventory import digest

        require(
            digest(unpack((case / "bundle/system.json").read_text()))
            == p["execution_input_identity"],
            "Changed serialized system or execution input",
        )
        require(
            [v["seed"] for v in c["relocation"]] == ["1", "99397"]
            and all(
                v["returncode"] == 0
                and v["source_files_relocated"] is True
                and v["report_identity"] == p["identity"]
                for v in c["relocation"]
            ),
            "Lost fresh-process relocation",
        )
    prior = read(root / "audit/prior-evidence.json")
    expected_prior = {
        str(p.relative_to(repo)): pin(p)
        for phase in range(47, 54)
        for p in (repo / f"docs/evidence/phase_4j{phase}").rglob("*.json")
        if any(
            k in p.name
            for k in ("component", "comparison", "convention", "resolution", "blocker")
        )
    }
    require(
        prior["files"] == expected_prior
        and prior["historical_claims_rewritten"] is False,
        "Changed historical diagnostics",
    )
    counts, nodes = junit(root, SUITES)
    expected_nodes = json.loads((root / "collection.json").read_text())["nodes"]
    require(
        Counter(nodes) == Counter(expected_nodes)
        and len(nodes) == len(set(nodes))
        and counts == r["full_suite_counts"]
        and r["junit_files"] == SUITES,
        "Incomplete test coverage",
    )
    focused, _ = junit(root, ["validation/focused.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    require(
        r["focused_tests_passed"] == focused["tests"] - focused["skipped"]
        and r["J48_J53_regression_counts"] == regressions
        and r["J54"] == "absent_from_exact_base",
        "Changed regression accounting",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(
        all(c["returncode"] == 0 for c in checks)
        and r["check_returncodes"] == {c["name"]: c["returncode"] for c in checks},
        "Failed check",
    )
    required_checks = (
        {"ruff", "compileall", "diff", "historical_hashes"}
        | {f"J{n}_receipt" for n in range(48, 54)}
        | {
            "pip-phase4e2",
            "pip-phase4g1",
            "pip-phase4j32",
            "pip-phase4j36",
            "pip-md_rg",
        }
    )
    require(required_checks <= set(r["check_returncodes"]), "Missing required check")
    return {
        "valid": True,
        "historical_files": len(old),
        "full_suite_counts": counts,
        "relocated_fresh_processes": 10,
        "production_validated": False,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--sha256")
    args = p.parse_args(argv)
    print(
        json.dumps(
            validate(args.directory, repo=args.repo, expected_sha256=args.sha256),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
