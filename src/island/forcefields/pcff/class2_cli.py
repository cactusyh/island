"""Resolve a saved automatic charge record without construction or embedding."""

import argparse
import json

from island.charge_references.records import pack
from island.exceptions import PCFFError

from .automatic import (
    PCFFAutomaticTypingResult,
    graph_system,
    load_pcff_automatic_record,
)
from .class2 import assign_pcff_parameters, save_pcff_parameters
from .source import load_pcff_source


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, help="Pinned local pcff.frc")
    parser.add_argument(
        "--charges", required=True, help="Saved automatic native-charge record"
    )
    parser.add_argument("--output", required=True, help="Exclusive result JSON path")
    parser.add_argument(
        "--require-complete",
        action="store_true",
        help="Exit nonzero on missing parameter coverage",
    )
    args = parser.parse_args(argv)
    try:
        source = load_pcff_source(args.source)
        charges = load_pcff_automatic_record(args.charges, source)
        data = charges.payload
        auto = data["automatic_typing"]
        typing = PCFFAutomaticTypingResult(pack(auto), source)
        result = assign_pcff_parameters(graph_system(auto["graph"]), typing, charges)
        save_pcff_parameters(result, args.output)
        report = result.payload
        print(
            json.dumps(
                {
                    k: report[k]
                    for k in (
                        "coverage",
                        "parameter_coverage_complete",
                        "physical_model_complete",
                        "production_validated",
                        "simulation_readiness",
                    )
                },
                indent=2,
            )
        )
        return int(args.require_complete and not report["parameter_coverage_complete"])
    except (PCFFError, KeyError) as error:
        print(json.dumps({"error": str(error)}))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
