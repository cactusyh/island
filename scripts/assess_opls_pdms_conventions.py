"""Reproduce explicit Figure 7 E/M partitions; numerical acceptance stays nonzero."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_conventions import (
    assess_opls_pdms_conventions,
    load_opls_pdms_conventions,
    source_evidence,
)
from island.forcefields.oplsaa.pdms_literature import (
    OPLSPDMSLiteratureError,
    build_opls_pdms_literature,
)
from island.forcefields.pdms_literature import INPUT_ERRORS


def fixtures():
    yield "hmds", from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026), None
    for dp in (2, 4):
        system = build_linear_polymer(
            "[*:1]O[Si](C)(C)[*:2]",
            dp=dp,
            random_seed=2026,
            add_hydrogens=True,
            head_end_group="[*][Si](C)(C)C",
            tail_end_group="[*]C",
        )
        # Explicit input from the validated builder, recorded with the exact graph.
        # The library never chooses a terminal oxygen by atom number or adjacency.
        yield f"tms_capped_dp{dp}", system, system.metadata["polymer"]["head_site_id"]
        yield f"tms_capped_dp{dp}_unoriented", system, None
    yield "unsupported_sih3", from_smiles("[SiH3]O[SiH3]", random_seed=2026), None


def summary(result):
    p = result.payload
    return {
        k: p[k]
        for k in (
            "data_identity",
            "identity",
            "status",
            "graph_identity",
            "atom_order",
            "charge_partition",
            "diagnostics",
            "missing_dependencies",
            "stages",
            "profile_identity",
            "execution_profiles",
            "backend_executions",
            "acceptance_gate",
            "limitations",
        )
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--reconstruct", type=Path)
    parser.add_argument("--sha256")
    parser.add_argument("--require-executable", action="store_true")
    args = parser.parse_args()
    if args.reconstruct and not args.sha256:
        parser.error("--reconstruct requires an independent --sha256")
    try:
        supplement = build_opls_pdms_literature(
            DomainSelection(
                args.selection.read_text(),
                args.catalog.read_text(),
            )
        )
        supplement.verify_sources(args.references)
        if args.reconstruct:
            report = summary(
                load_opls_pdms_conventions(
                    args.reconstruct,
                    args.sha256,
                    supplement=supplement,
                    references=args.references,
                )
            )
        else:
            # Complete all validation before publishing any output.
            cases = [
                (
                    name,
                    assess_opls_pdms_conventions(
                        system,
                        supplement,
                        references=args.references,
                        e_oxygen_site_id=oxygen,
                    ),
                )
                for name, system, oxygen in fixtures()
            ]
            args.output.mkdir(parents=True, exist_ok=False)
            report = {
                "cases": {},
                "acceptance_gate": 1,
                "execution_profiles": [],
                "backend_executions": [],
            }
            for name, result in cases:
                raw = result.json_text.encode()
                (args.output / f"{name}.json").write_bytes(raw)
                report["cases"][name] = {
                    **summary(result),
                    "file_sha256": sha256(raw).hexdigest(),
                }
            for name, document in (
                ("source-evidence.json", source_evidence(supplement)),
                ("report.json", report),
            ):
                (args.output / name).write_text(
                    json.dumps(document, ensure_ascii=False, indent=2) + "\n"
                )
        print(json.dumps(report, ensure_ascii=False, sort_keys=True))
        return 1 if args.require_executable else 0
    except INPUT_ERRORS as error:
        status = (
            error.status
            if isinstance(error, OPLSPDMSLiteratureError)
            else "asset_unavailable"
            if isinstance(error, FileNotFoundError)
            else "source_incompatible"
        )
        print(
            json.dumps(
                {
                    "status": status,
                    "diagnostics": getattr(
                        error,
                        "diagnostics",
                        [
                            {
                                "stage": "explicit_conventions_input_invalid",
                                "reason": str(error),
                            },
                        ],
                    ),
                    "acceptance_gate": 1,
                    "execution_profiles": [],
                }
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
