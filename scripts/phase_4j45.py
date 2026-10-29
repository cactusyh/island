"""Reproduce one source-bound J45 gap-closure target; acceptance remains nonzero."""

import argparse
import json
from pathlib import Path

from island.chemistry import from_smiles
from island.forcefields.external_pdms import LIMITS, require
from island.forcefields.literature_audit import fixture_requirements
from island.forcefields.literature_closure import (
    DomainClosure,
    DomainClosureBlockedError,
    DomainSelection,
    close_domain,
)
from island.workflows.bundle import system_data, system_from
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j41_c5 import CASES, fixture

BASE = "fdf8e8058173faa9c013703aa4910b3065820e7c"
J44_SNAPSHOT = "95531e5509b262265583e43c078afc34daa5c54a"
BRANCH = "codex/phase-4j45-complete-external-domain"


def fixtures():
    for case in CASES:
        yield case[0], fixture(case)
    yield "mixed_si_o_c", from_smiles("CO[Si](C)(C)C", random_seed=2026)


def summarize(name, result):
    p = result.payload
    try:
        result.require_executable()
    except DomainClosureBlockedError as error:
        boundary = {"exception": type(error).__name__, "status": error.status}
        require(
            error.result == result and error.diagnostics == p["missing_rows"],
            "lost_diagnostic",
        )
    else:
        raise ValueError("Unqualified literature domain reached numerical execution")
    return {
        "name": name,
        "graph_identity": p["requirements"]["graph_identity"],
        "status": p["status"],
        "missing_assignments": len(p["missing_rows"]),
        "literal_candidate_matches": sum(
            bool(v) for v in p["literal_row_candidates"].values()
        ),
        "unsupported_sites": p["unsupported_environments"],
        "stages": p["stages"],
        "boundary": boundary,
    }


def report(selection, cases, errors):
    return {
        "schema": "island.j45-study.v1",
        "base_commit": BASE,
        "j44_snapshot": J44_SNAPSHOT,
        "branch": BRANCH,
        "domain": selection.payload["selected_domain"],
        "selection_identity": selection.payload["identity"],
        "inventory_reused_from_j44": True,
        "known_source_files": len(selection.catalog["files"]),
        "source_diagnostics": errors,
        "cases": cases,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "native_pcff_acceptance_gate": 1,
        "limitations": LIMITS,
    }


def execute(output, references, selection_path, catalog_path):
    selection = DomainSelection(
        Path(selection_path).read_text(), Path(catalog_path).read_text()
    )
    errors = selection.verify_sources(references)
    root = Path(output)
    root.mkdir(parents=True, exist_ok=False)
    (root / "selection.json").write_text(selection.json_text)
    (root / "j44-source-audit.json").write_text(selection.catalog_text)
    write(
        root / "source-hashes.json",
        {r["filename"]: r["sha256"] for r in selection.catalog["files"]},
    )
    cases = []
    for name, system in fixtures():
        folder = root / name
        folder.mkdir()
        result = close_domain(system, selection, source_errors=errors)
        write(folder / "system.json", system_data(system))
        (folder / "gap-report.json").write_text(result.json_text)
        cases.append(summarize(name, result))
    result = report(selection, cases, errors)
    write(root / "report.json", result)
    return result


def reconstruct(root, references=None):
    root = Path(root)
    selection = DomainSelection(
        (root / "selection.json").read_text(),
        (root / "j44-source-audit.json").read_text(),
    )
    stored_report = read(root / "report.json")
    errors = stored_report["source_diagnostics"]
    if references is not None:
        require(selection.verify_sources(references) == errors, "source_state_changed")
    cases = []
    for name, canonical in fixtures():
        folder = root / name
        system = system_from(read(folder / "system.json"))
        require(
            fixture_requirements(system, "opls_aa")
            == fixture_requirements(canonical, "opls_aa"),
            "fixture_graph_or_atom_membership_changed",
        )
        stored = DomainClosure((folder / "gap-report.json").read_text())
        stored.validate_integrity(system, selection, errors)
        cases.append(summarize(name, stored))
    require(
        stored_report == report(selection, cases, errors),
        "study_status_or_stages_changed",
    )
    require(
        read(root / "source-hashes.json")
        == {r["filename"]: r["sha256"] for r in selection.catalog["files"]},
        "source_manifest_changed",
    )
    return {
        "reconstructed_cases": len(cases),
        "known_sources": len(selection.catalog["files"]),
        "execution_profiles": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--catalog", type=Path)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--reconstruct", type=Path)
    args = parser.parse_args()
    if args.reconstruct:
        print(json.dumps(reconstruct(args.reconstruct, args.references)))
        return 0
    if not all((args.references, args.selection, args.catalog)):
        parser.error(
            "--output requires explicit --references, --selection and --catalog"
        )
    print(
        json.dumps(execute(args.output, args.references, args.selection, args.catalog))
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
