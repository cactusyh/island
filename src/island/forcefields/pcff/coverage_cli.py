"""Inspect the full pinned inventory or assign an expanded graph profile.

The default inspection gate is structural; --require-full-source requires every
coverage stage and currently returns nonzero. Assignment never implies coverage.
"""

import argparse
import json
from pathlib import Path

from island.exceptions import IslandError
from island.workflows.storage import json_bytes, publish, read_json

from .automatic import (
    assign_automatic_pcff_charges,
    save_pcff_automatic_record,
    type_pcff_atoms,
)
from .catalog import inspect_pcff_full_source, resolve_pcff_source_record
from .class2 import assign_pcff_parameters, save_pcff_parameters
from .coverage import pcff_coverage_ledger
from .expanded import PROFILE_NAME, assign_pcff_source_types
from .source import FLAGS, load_pcff_source


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--output")
    p.add_argument("--smiles")
    p.add_argument(
        "--types",
        help="JSON object: integer stable IDs as decimal string keys to exact PCFF labels",
    )
    p.add_argument("--provenance", help="Required for explicit types")
    p.add_argument("--family")
    p.add_argument("--labels", nargs="+")
    p.add_argument("--namespace", default="cff91", choices=["cff91", "cff91_auto"])
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args(argv)
    try:
        source = load_pcff_source(a.source)
        if a.family:
            report = resolve_pcff_source_record(
                source, a.family, a.labels or [], namespace=a.namespace
            )
            passed = report["status"] == "assigned"
        elif a.smiles:
            from island.chemistry import from_smiles

            m = from_smiles(a.smiles, random_seed=2026)
            if a.types:
                raw = read_json(Path(a.types))
                if type(raw) is not dict or any(str(int(k)) != k for k in raw):
                    raise ValueError("Canonical decimal integer ID keys required")
                typing = assign_pcff_source_types(
                    m,
                    source,
                    {int(k): v for k, v in raw.items()},
                    provenance=a.provenance,
                )
            else:
                typing = type_pcff_atoms(m, source, profile=PROFILE_NAME)
            charge = (
                assign_automatic_pcff_charges(m, typing) if typing.complete else None
            )
            assignment = (
                assign_pcff_parameters(m, typing, charge)
                if charge and charge.complete
                else None
            )
            report = {
                "typing": typing.payload,
                "charges": charge.payload if charge else None,
                "parameter_coverage": assignment.payload["coverage"]
                if assignment
                else None,
                "full_source_complete": False,
                **FLAGS,
            }
            passed = (
                assignment is not None
                and assignment.payload["parameter_coverage_complete"]
            )
            if a.output:
                out = Path(a.output)
                out.mkdir(parents=True, exist_ok=False)
                save_pcff_automatic_record(typing, out / "typing.json")
                if charge:
                    save_pcff_automatic_record(charge, out / "charges.json")
                if assignment:
                    save_pcff_parameters(assignment, out / "parameters.json")
                publish(out / "report.json", json_bytes(report))
        else:
            if a.types:
                raise ValueError("--types requires --smiles")
            report = pcff_coverage_ledger(source)
            report["numeric_rows"] = len(inspect_pcff_full_source(source)["records"])
            passed = True
        if a.require_full_source:
            passed = False
        if a.output and not a.smiles:
            publish(Path(a.output), json_bytes(report))
        print(json.dumps(report, indent=2, allow_nan=False))
        return 0 if passed else 1
    except (IslandError, OSError, ValueError, TypeError) as error:
        print(json.dumps({"status": "failed", "error": str(error), **FLAGS}))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
