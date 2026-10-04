"""PSMILES -> automatic PCFF types/native charges -> explicit Class II coverage.

Run with --source /path/to/the/pinned/pcff.frc. No energy model is created.
"""

import argparse
import json

from island.builders import build_linear_polymer
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    load_pcff_source,
    type_pcff_atoms,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--psmiles", default="[*:1]CC[*:2]")
    parser.add_argument("--dp", type=int, default=3)
    args = parser.parse_args()
    system = build_linear_polymer(args.psmiles, dp=args.dp, generate_3d=False)
    source = load_pcff_source(args.source)
    typing = type_pcff_atoms(system, source)
    charges = assign_automatic_pcff_charges(system, typing)
    result = assign_pcff_parameters(system, typing, charges)
    p = result.payload
    print(
        json.dumps(
            {
                k: p[k]
                for k in (
                    "coverage",
                    "parameter_coverage_complete",
                    "physical_model_complete",
                )
            },
            indent=2,
        )
    )
    print("Native records only; no PCFF evaluator or simulation readiness.")


if __name__ == "__main__":
    main()
