"""Validate exact J53 source absence, nine blockers, artifacts and historical bytes."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import digest
from island.opls_conventions import signed
from island.opls_final_resolution import NINE_CONVENTIONS, corpus_audit
from island.opls_source_intake import blocker_reopen_report, inspect_source
from island.opls_source_intake.schema import source_intake_schema
from scripts.phase_4j53 import BASE, BRANCH
from scripts.seal_phase_4j51_receipt import junit
from scripts.seal_phase_4j53_receipt import IMPLEMENTATION
from scripts.validate_phase_4j47_receipt import history, pin, require
from scripts.validate_phase_4j50_receipt import read


def validate(directory, *, repo, expected_sha256=None):
    root, repo = Path(directory).resolve(), Path(repo).resolve()
    if expected_sha256:
        require(
            pin(root / "receipt.json") == expected_sha256,
            "Changed independent receipt pin",
        )
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["schema"] == "island.j53-receipt.v1"
        and m["schema"] == "island.j53-manifest.v1",
        "Wrong schemas",
    )
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Wrong exact base/branch",
    )
    require(
        r["public_family"] == "opls"
        and r["profile_identity"] == PROFILE
        and m["profile_identities"] == [PROFILE],
        "Wrong public/internal family",
    )
    require(
        r["source_status"] == "not_supplied"
        and r["new_profiles"] == m["new_profiles"] == [],
        "Invented source or profile",
    )
    require(
        r["permanent_blockers"] == list(NINE_CONVENTIONS), "Lost permanent blockers"
    )
    require(
        r["new_source_executable_acceptance_gate"]
        == r["research_opt_in_acceptance_gate"]
        == r["production_candidate_acceptance_gate"]
        == 1,
        "Promoted authority",
    )
    require(
        all(
            r[k] is False
            for k in (
                "MD_executed",
                "NPT_executed",
                "charges_parameters_or_tolerances_changed",
                "new_static_comparisons_executed",
                "production_validated",
            )
        )
        and r["simulation_readiness"] == "not_established"
        and r["checkpoint_RNG_continuation"] == "not run: no new authoritative source",
        "Promoted readiness or execution",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "Manifest changed")
    old = history(repo, BASE)
    require(m["historical_files"] == old, "Historical membership changed")
    require(all(pin(repo / n) == v for n, v in old.items()), "Historical bytes changed")
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "receipt.json", root / "manifest.json")
    }
    require(set(m["artifacts"]) == actual, "Evidence membership changed")
    require(
        set(m["implementation_files"]) == set(IMPLEMENTATION), "Changed code membership"
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
                "Changed artifact/code: " + name,
            )
    result = read(root / "audit/source-intake.json")
    require(result == inspect_source(), "Changed no-source intake")
    require(
        read(root / "audit/blocker-reopen.json") == blocker_reopen_report(result),
        "Changed exact blocker reopen evidence",
    )
    require(
        json.loads((root / "audit/source-intake-schema.json").read_text())
        == source_intake_schema(),
        "Changed intake contract",
    )
    audit = corpus_audit()
    require(
        r["J52_source_audit_identity"]
        == m["J52_source_audit_identity"]
        == audit["identity"],
        "Changed exhausted corpus identity",
    )
    delta = read(root / "audit/source-delta.json")
    require(
        delta["base_commit"] == BASE
        and not delta["new_content_candidates"]
        and delta["explicit_new_source_package"] is None
        and delta["no_network_discovery"],
        "Invented new source evidence",
    )
    require(
        len({f["path"] for f in delta["files"]}) == len(delta["files"])
        and m["source_hashes"] == {f["path"]: f["sha256"] for f in delta["files"]},
        "Changed source inventory",
    )
    for file in delta["files"]:
        require(file["previous_evidence"], "New asset requires explicit intake")
        for name in file["previous_evidence"]:
            require(
                name in old and file["sha256"] in (repo / name).read_text(),
                "Invented previous source pin",
            )
    previous = read(root / "audit/prior-comparison-records.json")
    expected_previous = {
        str(p.relative_to(repo)): pin(p)
        for phase in range(49, 53)
        for p in (repo / f"docs/evidence/phase_4j{phase}").rglob("*.json")
        if any(
            s in p.name
            for s in ("convention", "component", "resolution", "blocker", "comparison")
        )
    }
    require(
        previous["files"] == expected_previous
        and not previous["new_comparisons_run"]
        and not previous["source_model_changed"],
        "Changed historical comparison references",
    )
    summary = read(root / "audit/summary.json")
    require(
        summary["identity"] == r["summary_identity"]
        and summary["source_status"] == "not_supplied"
        and summary["new_source_packages"]
        == summary["new_profiles"]
        == summary["new_content_candidates"]
        == []
        and all(
            summary[k] is False
            for k in (
                "MD_executed",
                "NPT_executed",
                "new_static_comparisons_executed",
                "production_validated",
            )
        )
        and summary["simulation_readiness"] == "not_established"
        and summary["checkpoint_RNG_continuation"] == r["checkpoint_RNG_continuation"],
        "Changed audit outcome",
    )
    require(
        [v["dp"] for v in summary["cases"]] == [0, 2, 4, 6], "Changed fixture matrix"
    )
    for case in summary["cases"]:
        dp = case["dp"]
        path = root / f"audit/dp{dp}"
        value = read(path / "preparation.json")
        prior = read(
            repo
            / f"docs/evidence/phase_4j52/audit/dp{dp}/source_replay_v1/final-resolution.json"
        )
        assignment = read(
            repo
            / f"docs/evidence/phase_4j52/audit/dp{dp}/source_replay_v1/canonical-assignment.json"
        )
        expected = signed(
            {
                "schema": "island.opls-source-intake-preparation.v1",
                **{k: v for k, v in prior.items() if k not in ("schema", "identity")},
                "J52_resolution_identity": prior["identity"],
                "assignment_binding": {
                    "canonical_assignment_identity": assignment["identity"],
                    "graph_identity": assignment["graph_identity"],
                    "typed_graph_identity": digest(assignment["typed_graph"]),
                    "charge_identity": digest(assignment["graph_charges"]),
                    "source_identity": assignment["source_identity"],
                },
                "source_intake": result,
                "blocker_reopen": blocker_reopen_report(result),
                "new_source_selected": False,
                "selection_reason": "no_new_independently_qualified_authoritative_package",
            }
        )
        require(
            value == expected and value["identity"] == case["preparation_identity"],
            "Changed source/graph/charge/convention preparation",
        )
        require(
            pin(path / "system.json")
            == pin(
                repo
                / f"docs/evidence/phase_4j52/audit/dp{dp}/source_replay_v1/system.json"
            ),
            "Changed physical fixture",
        )
        bundle = read(path / "manifest.json")
        require(
            pin(path / "manifest.json") == case["bundle_manifest_sha256"]
            and bundle["files"]
            == {n: pin(path / n) for n in ("system.json", "preparation.json")}
            and bundle["preparation_identity"] == value["identity"]
            and not bundle["MD_authorized"],
            "Changed preparation bundle",
        )
        require(
            [v["action"] for v in case["blocked_actions"]] == ["minimize", "nvt", "npt"]
            and all(
                not v["outputs_created"]
                and not v["engine_started"]
                and v["stage"] == "new_authoritative_source_required"
                for v in case["blocked_actions"]
            ),
            "Blocked action executed",
        )
        require(
            [v["hash_seed"] for v in case["relocated_fresh_processes"]] == [1, 99397]
            and all(
                v["returncode"] == 0
                and not v["stderr"]
                and json.loads(v["stdout"])
                == {"identity": value["identity"], "MD_authorized": False}
                for v in case["relocated_fresh_processes"]
            ),
            "Missing fresh process reconstruction",
        )
    counts, nodes = junit(root, r["junit_files"])
    collection = json.loads((root / "collection.json").read_text())["nodes"]
    base = json.loads((repo / "docs/evidence/phase_4j52/collection.json").read_text())[
        "nodes"
    ]
    focused, new_nodes = junit(root, ["validation/focused.xml"])
    regression, regression_nodes = junit(root, ["validation/j49-j52-regression.xml"])
    require(
        Counter(nodes) == Counter(collection) == Counter(base + new_nodes)
        and all(v == 1 for v in Counter(nodes).values())
        and all(n.startswith("tests/test_phase_4j53.py::") for n in new_nodes)
        and counts == r["full_suite_counts"],
        "Incomplete ordinary test accounting",
    )
    require(
        Counter(regression_nodes)
        == Counter(
            [
                n
                for n in base
                if any(
                    n.startswith(f"tests/test_phase_4j{phase}.py::")
                    for phase in range(49, 53)
                )
            ]
        ),
        "Incomplete J49-J52 regression accounting",
    )
    require(
        r["focused_tests_passed"] == focused["tests"] - focused["skipped"]
        and r["J49_J52_regressions_passed"]
        == regression["tests"] - regression["skipped"],
        "Changed focused/regression counts",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(
        r["check_returncodes"] == {c["name"]: c["returncode"] for c in checks}
        and all(v == 0 for v in r["check_returncodes"].values())
        and {
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
            "J52_receipt",
            "historical_hashes",
        }
        <= set(r["check_returncodes"]),
        "Missing required checks",
    )
    return {
        "receipt_valid": True,
        "manifest_valid": True,
        "historical_files": len(old),
        "source_status": "not_supplied",
        "permanent_blockers": list(NINE_CONVENTIONS),
        "full_suite_counts": counts,
        "production_validated": False,
        "simulation_readiness": "not_established",
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
