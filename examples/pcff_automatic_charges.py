"""PSMILES -> existing builder -> ISLAND PCFF rules -> native increments.

python examples/pcff_automatic_charges.py --source /local/pcff.frc --dp 3
"""

import argparse
import json

from island.builders import build_linear_polymer
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    load_pcff_source,
    type_pcff_atoms,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--psmiles", default="[*:1]CCO[*:2]")
    parser.add_argument("--dp", type=int, default=3)
    args = parser.parse_args()
    system = build_linear_polymer(args.psmiles, dp=args.dp, generate_3d=False)
    source = load_pcff_source(args.source)
    typing = type_pcff_atoms(system, source)
    charge = assign_automatic_pcff_charges(system, typing) if typing.complete else None
    print(
        json.dumps(
            {
                "typing": typing.payload["coverage"],
                "unsupported": typing.payload["diagnostics"],
                "charge": None
                if charge is None
                else {
                    k: charge.payload["native_charge_record"][k]
                    for k in ("complete", "total_charge", "components", "diagnostics")
                },
                "bonded_class_ii_parameters": "not_implemented",
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
            indent=2,
        )
    )
    return 0 if charge is not None and charge.complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
