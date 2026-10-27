"""Execute explicit literature ingestion. Exit 1 retains blocked acceptance."""

import argparse
import json
from pathlib import Path

from island.charge_references.records import pack
from island.forcefields import ExternalPDMSBlockedError, create_evaluator
from island.forcefields.external_pdms import LIMITS
from island.forcefields.external_pdms_adapter import ExternalPDMSPreparation
from island.forcefields.pdms_literature import LiteratureIngestion, ingest_literature
from island.workflows.bundle import system_data, system_from
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j41_c5 import CASES, fixture

BASE = "fdf8e8058173faa9c013703aa4910b3065820e7c"
BRANCH = "codex/phase-4j43-literature-pdms-parameter-ingestion"


def boundary(system, preparation):
    try:
        create_evaluator(system, preparation)
    except ExternalPDMSBlockedError as error:
        return {
            "exception": type(error).__name__,
            "status": error.status,
            "diagnostics": error.diagnostics,
        }
    raise ValueError("Unqualified literature reached a numerical evaluator")


def execute(output, **inputs):
    root = Path(output)
    root.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island.j43-study.v1",
        "base_commit": BASE,
        "branch": BRANCH,
        "source_supplied": inputs.get("source") is not None,
        "bundle_supplied": inputs.get("bundle") is not None,
        "literature_inventory": [],
        "extracted_parameter_rows": {},
        "cases": [],
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "native_pcff_acceptance_gate": 1,
        "limitations": LIMITS,
    }
    for spec in CASES:
        system = fixture(spec)
        result = ingest_literature(system, **inputs)
        folder = root / spec[0]
        folder.mkdir()
        write(folder / "system.json", system_data(system))
        (folder / "ingestion.json").write_text(result.json_text)
        p = result.payload
        prep = ExternalPDMSPreparation(pack(p["preparation"]))
        outcome = boundary(system, prep)
        write(folder / "evaluator-boundary.json", outcome)
        record = p["source_record"]
        if record is not None and not report["literature_inventory"]:
            report["literature_inventory"].append(record)
            report["extracted_parameter_rows"] = record["extracted_rows"]
        report["cases"].append(
            {
                "name": spec[0],
                "graph_identity": prep.payload["requirements"]["graph_identity"],
                "atoms": system.number_of_sites,
                "status": p["status"],
                "missing_interactions": len(p["missing_interactions"]),
                "unresolved_dependencies": len(p["unresolved_dependencies"]),
                "boundary": outcome,
            }
        )
    write(root / "literature-inventory.json", report["literature_inventory"])
    write(
        root / "extracted-parameter-rows.json",
        {
            "rows": report["extracted_parameter_rows"],
            "reason": "No user-supplied literature assets"
            if not report["source_supplied"]
            else "Only independently reconstructed supplied rows; not numerical qualification",
        },
    )
    write(root / "report.json", report)
    return report


def reconstruct(root):
    """Reconstruct absent-asset studies; supplied assets always need explicit pins."""
    root = Path(root)
    report = read(root / "report.json")
    if report["source_supplied"] or report["bundle_supplied"]:
        raise ValueError(
            "Reexecution requires explicit source/bundle paths and independent pins"
        )
    if [r["name"] for r in report["cases"]] != [c[0] for c in CASES]:
        raise ValueError("Changed fixture inventory")
    for row in report["cases"]:
        folder = root / row["name"]
        system = system_from(read(folder / "system.json"))
        stored = LiteratureIngestion((folder / "ingestion.json").read_text())
        fresh = ingest_literature(system)
        if fresh != stored:
            raise ValueError("Changed graph, inventory or ingestion outcome")
        p = stored.payload
        prep = ExternalPDMSPreparation(pack(p["preparation"]))
        prep.validate_integrity(system)
        expected = {
            "name": row["name"],
            "atoms": system.number_of_sites,
            "graph_identity": prep.payload["requirements"]["graph_identity"],
            "status": p["status"],
            "missing_interactions": len(p["missing_interactions"]),
            "unresolved_dependencies": len(p["unresolved_dependencies"]),
            "boundary": boundary(system, prep),
        }
        if row != expected or row["boundary"] != read(
            folder / "evaluator-boundary.json"
        ):
            raise ValueError("Changed boundary result")
    if (
        report["literature_inventory"] != []
        or report["extracted_parameter_rows"] != {}
        or read(root / "literature-inventory.json") != []
        or read(root / "extracted-parameter-rows.json")
        != {"rows": {}, "reason": "No user-supplied literature assets"}
        or report["execution_profiles"] != []
        or report["backend_executions"] != []
        or report["acceptance_gate"] != 1
        or report["native_pcff_acceptance_gate"] != 1
        or report["limitations"] != LIMITS
        or report["base_commit"] != BASE
        or report["branch"] != BRANCH
        or report["schema"] != "island.j43-study.v1"
    ):
        raise ValueError("Promoted or changed study")
    return {
        "reconstructed_cases": len(CASES),
        "acceptance_gate": 1,
        "execution_profiles": [],
        "backend_executions": [],
        "limitations": LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--source-sha256")
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--bundle-sha256")
    parser.add_argument("--reconstruct", type=Path)
    args = parser.parse_args()
    if args.reconstruct:
        if any(
            (
                args.output,
                args.source,
                args.source_sha256,
                args.bundle,
                args.bundle_sha256,
            )
        ):
            parser.error("--reconstruct takes only its explicit study directory")
        print(json.dumps(reconstruct(args.reconstruct)))
        return 0
    if not args.output:
        parser.error("--output is required")
    print(
        json.dumps(
            execute(
                args.output,
                source=args.source,
                source_sha256=args.source_sha256,
                bundle=args.bundle,
                bundle_sha256=args.bundle_sha256,
            )
        )
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
