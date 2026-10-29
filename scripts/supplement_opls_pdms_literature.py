"""Export/query audited OPLS-compatible PDMS data without preparing a force field."""

import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path

from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_literature import (
    OPLSPDMSLiteratureError,
    build_opls_pdms_literature,
    load_opls_pdms_literature,
)
from island.forcefields.pdms_literature import INPUT_ERRORS
from scripts.phase_4j45 import fixtures, summarize


def audit(supplement, references):
    p = supplement.payload
    return {
        "data_identity": p["data_identity"],
        "identity": supplement.identity,
        "row_counts": dict(Counter(row["kind"] for row in p["rows"])),
        "source_diagnostics": supplement.source_diagnostics(references),
        "cases": [
            summarize(name, supplement.assess(system, references=references))
            for name, system in fixtures()
        ],
        "status": p["status"],
        "profile_identity": None,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "limitations": p["limitations"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--reconstruct", type=Path)
    parser.add_argument("--sha256", help="Independent hash of the explicit export")
    parser.add_argument("--require-executable", action="store_true")
    args = parser.parse_args()
    if args.reconstruct and not args.sha256:
        parser.error("--reconstruct requires an independent --sha256")
    try:
        selection = DomainSelection(
            args.selection.read_text(), args.catalog.read_text()
        )
        if args.reconstruct:
            supplement = load_opls_pdms_literature(
                args.reconstruct,
                args.sha256,
                selection=selection,
                references=args.references,
            )
        else:
            supplement = build_opls_pdms_literature(selection)
            supplement.verify_sources(args.references)
        report = audit(supplement, args.references)
        raw = supplement.json_text.encode("utf-8")
        report["supplement_file_sha256"] = sha256(raw).hexdigest()
        if args.output:
            args.output.mkdir(parents=True, exist_ok=False)
            (args.output / "parameters.json").write_bytes(raw)
            (args.output / "report.json").write_text(
                json.dumps(report, ensure_ascii=False, indent=2) + "\n"
            )
        print(json.dumps(report, ensure_ascii=False, sort_keys=True))
        # Export success is distinct from the unchanged numerical acceptance gate.
        return 1 if args.require_executable else 0
    except OPLSPDMSLiteratureError as error:
        print(
            json.dumps(
                {
                    "status": error.status,
                    "diagnostics": error.diagnostics,
                    "execution_profiles": [],
                    "acceptance_gate": 1,
                }
            )
        )
        return 1
    except INPUT_ERRORS as error:
        print(
            json.dumps(
                {
                    "status": "asset_unavailable"
                    if isinstance(error, FileNotFoundError)
                    else "source_incompatible",
                    "diagnostics": [
                        {"stage": "explicit_evidence_invalid", "reason": str(error)}
                    ],
                    "execution_profiles": [],
                    "acceptance_gate": 1,
                }
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
