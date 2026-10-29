"""Explicit offline new-graph mapping; blocked acceptance stays nonzero."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.opls_external.pdms_mapping import (
    load_mapping,
    map_author_pdms,
    source_rights,
)
from island.forcefields.oplsaa.pdms_author import load_author_opls_audit
from island.forcefields.oplsaa.pdms_literature import (
    OPLSPDMSLiteratureError,
    build_opls_pdms_literature,
)
from island.forcefields.pdms_literature import INPUT_ERRORS
from island.workflows.bundle import system_data, system_from
from scripts.assess_opls_pdms_conventions import fixtures


def read_author(args):
    supplement = build_opls_pdms_literature(
        DomainSelection(args.selection.read_text(), args.catalog.read_text())
    )
    return load_author_opls_audit(
        args.author_record,
        args.author_record_sha256,
        sources=args.sources,
        supplement=supplement,
        references=args.references,
        engine_docs=args.engine_docs,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "sources",
        "references",
        "engine-docs",
        "selection",
        "catalog",
        "author-record",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--author-record-sha256", required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--reconstruct", type=Path)
    parser.add_argument("--sha256")
    parser.add_argument("--system", type=Path)
    parser.add_argument("--e-oxygen-site-id", type=int)
    parser.add_argument("--require-executable", action="store_true")
    args = parser.parse_args()
    author = read_author(args)
    if args.reconstruct:
        if not (args.sha256 and args.system):
            parser.error("reconstruction requires --sha256 and --system")
        system = system_from(unpack(args.system.read_text()))
        result = load_mapping(
            args.reconstruct,
            args.sha256,
            system,
            author,
            references=args.references,
            e_oxygen_site_id=args.e_oxygen_site_id,
        )
        if args.require_executable:
            try:
                result.require_executable(system)
            except OPLSPDMSLiteratureError as error:
                print(
                    json.dumps(
                        {"status": error.status, "diagnostics": error.diagnostics}
                    )
                )
                return 1
        print(
            json.dumps(
                {
                    "identity": result.payload["identity"],
                    "status": result.payload["status"],
                    "reconstructed": True,
                }
            )
        )
        return 0
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "source_rights": source_rights(args.references),
        "cases": {},
        "PDMS_executable": False,
        "PDMS_MD_run": False,
        "execution_profiles": [],
        "acceptance_gate": 1,
    }
    for name, system, oxygen in fixtures():
        try:
            result = map_author_pdms(
                system, author, references=args.references, e_oxygen_site_id=oxygen
            )
            raw = result.json_text.encode()
            (args.output / (name + ".json")).write_bytes(raw)
            (args.output / (name + "-system.json")).write_text(
                pack(system_data(system))
            )
            p = result.payload
            report["cases"][name] = {
                "status": p["status"],
                "mapping_identity": p["identity"],
                "file_sha256": sha256(raw).hexdigest(),
                "atoms": len(p["mapping"]),
                "E_oxygen_site_id": p["E_oxygen_site_id"],
                "charge_policy": p["charge_policy"],
                "missing_dependencies": p["missing_dependencies"],
            }
        except OPLSPDMSLiteratureError as error:
            report["cases"][name] = {
                "status": error.status,
                "diagnostics": error.diagnostics,
            }
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 1 if args.require_executable else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except OPLSPDMSLiteratureError as error:
        print(json.dumps({"status": error.status, "diagnostics": error.diagnostics}))
        raise SystemExit(1) from None
    except INPUT_ERRORS as error:
        print(
            json.dumps(
                {
                    "status": "typing_blocked",
                    "diagnostics": [
                        {"stage": "mapping_input_invalid", "reason": str(error)}
                    ],
                }
            )
        )
        raise SystemExit(1) from None
