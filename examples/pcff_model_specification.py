"""Load an H3 assignment and define a named model with explicit pair choices.

No typing, embedding, QM or simulation is performed. Pair weights are a caller
model choice, not values recovered from pcff.frc.
"""

import argparse
import json

from island.forcefields.pcff import (
    define_pcff_model,
    load_pcff_parameters,
    load_pcff_source,
    save_pcff_model,
    special_pair_policy,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--assignment", required=True)
    parser.add_argument("--lj", required=True, type=float, nargs=3)
    parser.add_argument("--coulomb", required=True, type=float, nargs=3)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    assignment = load_pcff_parameters(args.assignment, load_pcff_source(args.source))
    result = define_pcff_model(
        assignment, special_pairs=special_pair_policy(lj=args.lj, coulomb=args.coulomb)
    )
    save_pcff_model(result, args.output)
    p = result.payload
    print(
        json.dumps(
            {
                k: p[k]
                for k in (
                    "assignment_identity",
                    "raw_parameter_coverage_complete",
                    "model_definition_complete",
                    "term_origins",
                    "numerical_verification",
                    "diagnostics",
                    "production_validated",
                    "simulation_readiness",
                )
            },
            indent=2,
        )
    )
    return int(not p["model_definition_complete"])


if __name__ == "__main__":
    raise SystemExit(main())
