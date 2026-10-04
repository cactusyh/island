"""Construct a graph, automatically type it, and report separate charge coverage."""

import argparse
import json
from pathlib import Path

from island.exceptions import IslandError
from island.workflows.storage import json_bytes, publish

from .automatic import (
    assign_automatic_pcff_charges,
    save_pcff_automatic_record,
    type_pcff_atoms,
)
from .source import FLAGS, load_pcff_source


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", required=True, help="Pinned local pcff.frc; no download"
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--psmiles")
    inputs.add_argument("--smiles")
    parser.add_argument("--dp", type=int, default=3)
    parser.add_argument("--output", help="New exclusive output directory; optional")
    args = parser.parse_args(argv)
    try:
        source = load_pcff_source(args.source)
        if args.psmiles:
            from island.builders import build_linear_polymer

            system = build_linear_polymer(args.psmiles, dp=args.dp, generate_3d=False)
        else:
            from island.chemistry import from_smiles

            # Construction uses its established embedding API, not the typer.
            system = from_smiles(args.smiles, add_hydrogens=True, random_seed=2026)
        typing = type_pcff_atoms(system, source)
        charge = (
            assign_automatic_pcff_charges(system, typing) if typing.complete else None
        )
        report = {
            "source": source.identity,
            "typing_identity": typing.identity,
            "typing_coverage": typing.payload["coverage"],
            "typing_diagnostics": typing.payload["diagnostics"],
            "assignments": typing.payload["assignments"],
            "charge_coverage": None,
            "class_ii_parameter_coverage": "not_implemented",
            **FLAGS,
        }
        if charge is not None:
            p = charge.payload["native_charge_record"]
            report["charge_coverage"] = {
                k: p[k]
                for k in (
                    "complete",
                    "total_charge",
                    "components",
                    "diagnostics",
                    "unit",
                )
            }
            report["automatic_charge_identity"] = charge.identity
        if args.output:
            root = Path(args.output)
            root.mkdir()
            save_pcff_automatic_record(typing, root / "typing.json")
            if charge is not None:
                save_pcff_automatic_record(charge, root / "charges.json")
            publish(root / "report.json", json_bytes(report))
        print(json.dumps(report, indent=2, allow_nan=False))
        return 0 if charge is not None and charge.complete else 1
    except (IslandError, ImportError, OSError, ValueError) as error:
        print(json.dumps({"status": "failed", "error": str(error), **FLAGS}))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
