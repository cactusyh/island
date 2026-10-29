#!/usr/bin/env python
"""Prepare an explicitly pinned external OPLS source from SMILES/PSMILES."""

import argparse
import json
from pathlib import Path

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields.extensible import (
    ExternalOPLSError,
    ExternalOPLSOptions,
    ForceFieldRequest,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.opls_external.storage import save_preparation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--smiles")
    inputs.add_argument("--psmiles")
    parser.add_argument("--dp", type=int)
    parser.add_argument("--head-end-group")
    parser.add_argument("--tail-end-group")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--source-sha256")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--evaluate", action="store_true")
    parser.add_argument("--minimize", action="store_true")
    args = parser.parse_args()
    if bool(args.source) != bool(args.source_sha256):
        parser.error(
            "Supply both an explicit source manifest and its independent SHA-256"
        )
    if args.psmiles and not (args.dp and args.head_end_group and args.tail_end_group):
        parser.error("PSMILES requires DP and explicit head/tail end groups")
    if args.output.exists():
        parser.error("Output must be a new directory")
    if args.psmiles:
        system = build_linear_polymer(
            args.psmiles,
            dp=args.dp,
            head_end_group=args.head_end_group,
            tail_end_group=args.tail_end_group,
            add_hydrogens=True,
            random_seed=2026,
        )
    else:
        system = from_smiles(args.smiles, random_seed=2026)
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "external_opls", ExternalOPLSOptions(args.source, args.source_sha256)
        ),
    )
    args.output.mkdir(parents=True)
    saved_pin = save_preparation(prepared, args.output / "preparation.json")
    p = prepared.payload
    report = {
        "status": p["status"],
        "diagnostics": p["diagnostics"],
        "profile_identity": p["profile_identity"],
        "preparation_identity": p["identity"],
        "preparation_file_sha256": saved_pin,
        "native_opls_coverage": False,
        "limits": p["limitations"],
        "MD_run": False,
    }
    try:
        prepared.require_executable(system)
        if args.evaluate or args.minimize:
            evaluator = create_evaluator(system, prepared)
            e = evaluator.evaluate()
            report["evaluation"] = {
                "energy_kJ_mol": e.potential_energy,
                "components": dict(e.energy_components),
                "forces_kJ_mol_A": dict(e.forces),
            }
            if args.minimize:
                from island.minimization import minimize_geometry

                minimum = minimize_geometry(system, evaluator)
                report["minimization"] = {
                    "converged": minimum.converged,
                    "reason": minimum.termination_reason,
                    "final_energy_kJ_mol": minimum.final_evaluation.potential_energy
                    if minimum.final_evaluation is not None
                    else None,
                    "final_evaluation_verified": minimum.final_evaluation_verified,
                }
    except ExternalOPLSError as error:
        report.update(status=error.status, diagnostics=error.diagnostics)
    raw = json.dumps(report, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    (args.output / "report.json").write_text(raw)
    print(raw)
    return 0 if report["status"] == "prepared" else 1


if __name__ == "__main__":
    raise SystemExit(main())
