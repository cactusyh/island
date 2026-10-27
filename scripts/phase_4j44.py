"""Inventory only the explicit local reference directory; acceptance exits 1."""

import argparse
import json
from pathlib import Path

from island.chemistry import from_smiles
from island.forcefields.external_pdms import LIMITS
from island.forcefields.literature_audit import (
    LiteratureAudit,
    LiteratureBlockedError,
    LiteratureBundleAudit,
    LiteratureSourceReview,
    audit_parameter_bundle,
    inventory_literature,
)
from island.workflows.bundle import system_data, system_from
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j41_c5 import CASES, fixture

BASE = "fdf8e8058173faa9c013703aa4910b3065820e7c"
J43_SNAPSHOT = "be8d94ec6e526d8cf8251ceb7c1baef86d609ee5"
BRANCH = "codex/phase-4j44-local-literature-parameter-audit"
CONTROLS = (
    ("mixed_silicon_methylsilane", "C[SiH3]", "class2"),
    ("opls_aliphatic_ethane", "CC", "opls_aa"),
    ("opls_aromatic_benzene", "c1ccccc1", "opls_aa"),
    ("pcff_ethanol", "CCO", "class2"),
    ("opls_dimethyl_phosphate_anion", "COP(=O)([O-])OC", "opls_aa"),
)


def fixtures():
    for spec in CASES:
        yield spec[0], fixture(spec), "class2"
    for name, smiles, form in CONTROLS:
        yield name, from_smiles(smiles, random_seed=2026), form


def boundary(result):
    try:
        result.require_executable()
    except LiteratureBlockedError as error:
        return {"exception": type(error).__name__, "status": error.status}
    raise ValueError("Incomplete literature reached numerical execution")


def execute(output, references, catalog=None):
    root = Path(output)
    root.mkdir(parents=True, exist_ok=False)
    audit = inventory_literature(references)
    (root / "inventory.json").write_text(audit.json_text)
    review = None
    if catalog is not None:
        review = LiteratureSourceReview(Path(catalog).read_text())
        review.validate_inventory(audit)
        (root / "source-audit.json").write_text(review.json_text)
    report = {
        "schema": "island.j44-study.v1",
        "base_commit": BASE,
        "j43_snapshot": J43_SNAPSHOT,
        "branch": BRANCH,
        "inventory_identity": audit.identity,
        "source_files": len(audit.payload["files"]),
        "stages": audit.payload["stages"],
        "cases": [],
        "phosphorus": "graph_only_dimethyl_phosphate_anion; no state transfer or execution",
        "source_review_identity": review.identity if review else None,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    for name, system, form in fixtures():
        folder = root / name
        folder.mkdir()
        result = audit_parameter_bundle(system, root=references, audit=audit, form=form)
        (folder / "bundle-audit.json").write_text(result.json_text)
        write(folder / "system.json", system_data(system))
        p = result.payload
        report["cases"].append(
            {
                "name": name,
                "functional_form": form,
                "atoms": system.number_of_sites,
                "graph_identity": p["requirements"]["graph_identity"],
                "status": p["status"],
                "missing_rows": len(p["missing_rows"]),
                "unresolved_dependencies": len(p["unresolved_dependencies"]),
                "boundary": boundary(result),
            }
        )
    write(
        root / "source-hashes.json",
        {f["filename"]: f["sha256"] for f in audit.payload["files"]},
    )
    write(
        root / "extracted-parameter-rows.json",
        {
            "rows": review.payload["extracted_row_candidates"] if review else [],
            "accepted_rows": [],
            "reason": "Candidate decimal evidence is not a complete independently qualified parameter profile",
        },
    )
    write(root / "report.json", report)
    return report


def reconstruct(root, references=None):
    root = Path(root)
    report = read(root / "report.json")
    audit = LiteratureAudit((root / "inventory.json").read_text())
    if references is not None:
        audit.validate_integrity(references)
    review = (
        LiteratureSourceReview((root / "source-audit.json").read_text())
        if (root / "source-audit.json").exists()
        else None
    )
    if review is not None:
        review.validate_inventory(audit)
    expected = []
    for name, _, form in fixtures():
        folder = root / name
        system = system_from(read(folder / "system.json"))
        result = audit_parameter_bundle(system, root=references, audit=audit, form=form)
        stored = LiteratureBundleAudit((folder / "bundle-audit.json").read_text())
        if stored != result:
            raise ValueError("Changed source, graph, requirements or blocked outcome")
        p = result.payload
        expected.append(
            {
                "name": name,
                "functional_form": form,
                "atoms": system.number_of_sites,
                "graph_identity": p["requirements"]["graph_identity"],
                "status": p["status"],
                "missing_rows": len(p["missing_rows"]),
                "unresolved_dependencies": len(p["unresolved_dependencies"]),
                "boundary": boundary(result),
            }
        )
    original = {
        "schema": "island.j44-study.v1",
        "base_commit": BASE,
        "j43_snapshot": J43_SNAPSHOT,
        "branch": BRANCH,
        "inventory_identity": audit.identity,
        "source_files": len(audit.payload["files"]),
        "stages": audit.payload["stages"],
        "cases": expected,
        "phosphorus": "graph_only_dimethyl_phosphate_anion; no state transfer or execution",
        "source_review_identity": review.identity if review else None,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    if report != original:
        raise ValueError("Changed study stages or fixture identities")
    if read(root / "source-hashes.json") != {
        f["filename"]: f["sha256"] for f in audit.payload["files"]
    }:
        raise ValueError("Changed source hash manifest")
    if read(root / "extracted-parameter-rows.json") != {
        "rows": review.payload["extracted_row_candidates"] if review else [],
        "accepted_rows": [],
        "reason": "Candidate decimal evidence is not a complete independently qualified parameter profile",
    }:
        raise ValueError("Forged parameter extraction")
    return {
        "reconstructed_cases": len(expected),
        "source_files": report["source_files"],
        "acceptance_gate": 1,
        "execution_profiles": [],
        "limitations": LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reconstruct", type=Path)
    args = parser.parse_args()
    if args.reconstruct:
        if args.output:
            parser.error("--output and --reconstruct are mutually exclusive")
        print(json.dumps(reconstruct(args.reconstruct, args.references)))
        return 0
    if not args.output:
        parser.error("--output is required")
    print(json.dumps(execute(args.output, args.references, args.catalog)))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
