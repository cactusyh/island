"""Validate J54 source-bound incompleteness, independent statuses and exact history."""

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from island.forcefields.selection.inventory import canonical, digest
from island.pcff_pdms import CANDIDATE, PCFFPDMSAssets, reconstruct
from scripts.audit_pcff_pdms_class2_sources import archive_audit, independent_frc
from scripts.phase_4j54 import BASE, BRANCH, repeat_audit
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history, pin, require

IMPLEMENTATION = (
    "src/island/pcff_pdms/__init__.py",
    "src/island/pcff_pdms/__main__.py",
    "scripts/phase_4j54.py",
    "scripts/audit_pcff_pdms_class2_sources.py",
    "scripts/pcff_pdms_source_context.py",
    "scripts/seal_phase_4j54_receipt.py",
    "scripts/validate_phase_4j54_receipt.py",
    "scripts/check_phase_4j54_receipt.py",
    "scripts/select_phase_4j54_validation.py",
    "tests/test_phase_4j54.py",
    "docs/phase_4j54.md",
)
SUITES = [f"validation/full-{i}.xml" for i in range(8) if i != 5] + [
    "validation/full-5-passed.xml",
    "validation/full-5-rerun.xml",
    "validation/focused.xml",
]


def read(path):
    p = json.loads(Path(path).read_text())
    require(
        p["identity"] == digest({k: v for k, v in p.items() if k != "identity"}),
        "Changed signed record: " + str(path),
    )
    return p


def validate_audit(directory, *, distribution, archive):
    root = Path(directory)
    summary = read(root / "summary.json")
    require(
        summary["candidate_profile"] == CANDIDATE
        and summary["registered_profiles"] == []
        and not summary["MD_executed"]
        and not summary["production_validated"]
        and summary["simulation_readiness"] == "not_established"
        and summary["backend_execution"] == []
        and summary["parameter_completeness_independent_of_cross_engine"],
        "Changed audit authority",
    )
    require(
        [c["fixture"] for c in summary["fixtures"]]
        == ["disiloxane", "hmds", "dp2", "dp4", "dp6"],
        "Changed fixture membership",
    )
    require(
        json.loads((root / "repeat.json").read_text()) == repeat_audit(),
        "Changed unfinished repeat",
    )
    all_rows, gaps, matrix = {}, defaultdict(list), []
    assets = PCFFPDMSAssets(distribution)
    for c in summary["fixtures"]:
        name = c["fixture"]
        s, result = reconstruct(root / name, c["manifest_sha256"], assets=assets)
        p = result.validate_integrity(s)
        require(
            p["identity"] == c["report_identity"]
            and p["graph_identity"] == c["graph_identity"]
            and len(p["missing_terms"]) == c["missing_requests"]
            and s.number_of_sites == c["atoms"],
            "Fixture identity/count changed",
        )
        for field in (
            "parameter_completeness",
            "typing_status",
            "charge_status",
            "lammps_export_status",
            "lammps_runtime_status",
            "openmm_build_status",
            "cross_engine_equivalence",
            "scientific_validation",
            "production_readiness",
            "coverage",
            "missing_conventions",
        ):
            require(p[field] == c[field], "Changed independent status: " + field)
        independent = independent_frc(s, p, Path(distribution) / "pcff.frc")
        require(
            canonical(independent)
            == canonical(json.loads((root / (name + "-independent.json")).read_text())),
            "Independent reconstruction changed",
        )
        require(
            [r["process_hash_seed"] for r in c["fresh_processes"]] == [1, 99397]
            and all(
                r["returncode"] == 0 and r["record_identity"] == p["identity"]
                for r in c["fresh_processes"]
            ),
            "Missing fresh-process evidence",
        )
        all_rows.update(p["source_rows"])
        for row in p["missing_terms"]:
            gaps[row["family"], row["query"]].append(
                {
                    "fixture": name,
                    "id": row["id"],
                    "sites": row["sites"],
                    "dependencies": row["dependencies"],
                    "resolved_types": row["source_lookup"]["resolved_types"],
                }
            )
        for family, coverage in sorted(p["coverage"].items()):
            matrix.append(
                {
                    "fixture": name,
                    "family": family,
                    "source_rows": coverage,
                    "dependencies": sorted(
                        {
                            d.split(":")[0]
                            for row in p["interaction_inventory"]
                            if row["family"] == family
                            for d in row["dependencies"]
                        }
                    ),
                }
            )
    require(
        read(root / "missing-rows.json")["queries"]
        == [
            {"family": f, "query": q, "requests": r}
            for (f, q), r in sorted(gaps.items())
        ],
        "Changed missing rows/dependencies",
    )
    require(
        read(root / "class2-dependency-matrix.json")["rows"] == matrix,
        "Changed Class-II dependency matrix",
    )
    author = archive_audit(archive)
    require(
        canonical(author)
        == canonical(json.loads((root / "author-compass.json").read_text())),
        "Changed separate COMPASS evidence",
    )
    mapping = read(root / "source-to-row.json")
    require(
        mapping["frc_rows"] == all_rows
        and mapping["author_compass_rows"] == author["active_coefficient_rows"]
        and mapping["overlays_authorized"] is False,
        "Changed source mapping or authorized borrowing",
    )
    require(
        summary["missing_frc_requests"]
        == sum(c["missing_requests"] for c in summary["fixtures"])
        and summary["missing_frc_source_queries"] == len(gaps)
        and summary["independent_rows"] == len(all_rows),
        "Changed gap totals",
    )
    return summary


def validate(directory, *, repo, distribution, archive, references):
    root, repo = Path(directory), Path(repo)
    r, m = read(root / "receipt.json"), read(root / "manifest.json")
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "Changed exact base/branch",
    )
    require(
        r["schema"] == "island.j54-receipt.v1"
        and m["schema"] == "island.j54-manifest.v1",
        "Wrong schema",
    )
    require(r["manifest_sha256"] == pin(root / "manifest.json"), "Changed manifest")
    old = history(repo, BASE)
    require(
        old == m["historical_files"] and len(old) == r["historical_files_preserved"],
        "Changed historical membership",
    )
    require(all(pin(repo / n) == v for n, v in old.items()), "Historical bytes changed")
    require(
        set(m["implementation_files"]) == set(IMPLEMENTATION),
        "Changed implementation membership",
    )
    actual = {
        str(p.relative_to(root)): pin(p)
        for p in root.rglob("*")
        if p.is_file() and p not in (root / "receipt.json", root / "manifest.json")
    }
    require(m["artifacts"] == actual, "Changed evidence files")
    require(
        all(pin(repo / n) == v for n, v in m["implementation_files"].items()),
        "Implementation changed",
    )
    summary = validate_audit(root / "audit", distribution=distribution, archive=archive)
    require(
        summary["identity"] == r["summary_identity"]
        and summary["source_hashes"] == m["source_hashes"],
        "Changed summary/source identity",
    )
    context = read(root / "audit/literature-context.json")
    from scripts.pcff_pdms_source_context import CATALOG, CATALOG_PIN, SCOPE

    require(
        context["catalog"] == {"file": CATALOG, "sha256": CATALOG_PIN}
        and pin(repo / CATALOG) == CATALOG_PIN,
        "Changed prior literature catalog",
    )
    catalog = json.loads((repo / CATALOG).read_text())
    require(
        context["source_row_mapping"]
        == [
            {
                "id": row["id"],
                "source": row["source"],
                "family": row["family"],
                "labels": row["labels"],
                "location": row["source_evidence"],
                "evidence_checksum": row["evidence_checksum"],
            }
            for row in catalog["rows"]
        ],
        "Changed literature row mapping",
    )
    require(
        {s["file"]: s["scope"] for s in context["sources"]} == SCOPE,
        "Changed source comparison scope",
    )
    require(
        context["independent_original_PDF_reconstruction"][
            "parameter_rows_reconstructed"
        ]
        == 105
        and context["independent_original_PDF_reconstruction"][
            "partial_rows_not_zero_filled"
        ]
        == 8,
        "Changed paper reconstruction",
    )
    require(
        all(
            pin(Path(references) / s["file"]) == s["sha256"] for s in context["sources"]
        ),
        "Literature source replaced",
    )
    require(
        m["literature_sources"] == context["sources"],
        "Changed literature source hashes",
    )
    require(
        r["candidate_profile"] == CANDIDATE
        and r["registered_profiles"] == []
        and r["parameter_completeness"] == "incomplete"
        and r["cross_engine_equivalence"] == "not_claimed"
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
                "parameters_charges_or_tolerances_changed",
                "historical_claims_rewritten",
            )
        ),
        "Promoted readiness/profile authority",
    )
    from scripts.select_phase_4j54_validation import selection

    retry = selection(root)
    require(r["ordinary_suite_retry"] == retry, "Changed historical failure disclosure")
    counts, nodes = junit(root, SUITES)
    require(
        Counter(nodes)
        == Counter(json.loads((root / "collection.json").read_text())["nodes"])
        and len(nodes) == len(set(nodes)),
        "Incomplete or duplicate full-suite membership",
    )
    require(counts == r["full_suite_counts"], "Changed full-suite counts")
    focused, _ = junit(root, ["validation/focused.xml"])
    regressions, _ = junit(root, ["validation/regressions.xml"])
    require(
        r["focused_counts"] == focused and r["J55_regression_counts"] == regressions,
        "Changed focused/regression counts",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    expected = {"ruff", "compileall", "diff", "historical_hashes", "J55_receipt"} | {
        "pip-" + n for n in ("phase4e2", "phase4g1", "phase4j32", "phase4j36", "md_rg")
    }
    require(
        expected <= {c["name"] for c in checks}
        and all(c["returncode"] == 0 for c in checks)
        and r["check_returncodes"] == {c["name"]: c["returncode"] for c in checks},
        "Required check failed/absent",
    )
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(old),
        "full_suite_counts": counts,
        "registered_profiles": [],
        "parameter_completeness": "incomplete",
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--distribution", type=Path, required=True)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--references", type=Path, required=True)
    a = p.parse_args(argv)
    print(
        json.dumps(
            validate(
                a.directory,
                repo=a.repo,
                distribution=a.distribution,
                archive=a.archive,
                references=a.references,
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
