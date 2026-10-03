"""Explicit installed source; no embedding, QM or full parameterization."""

import argparse
import json

from island.builders import build_linear_polymer
from island.forcefields.oplsaa import (
    assign_native_charges,
    load_oplsaa_source,
    type_atoms,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml", required=True)
    parser.add_argument("--psmiles", default="[*:1]CC[*:2]")
    parser.add_argument("--dp", type=int, default=3)
    args = parser.parse_args()
    system = build_linear_polymer(args.psmiles, dp=args.dp, generate_3d=False)
    source = load_oplsaa_source(args.xml)
    typing = type_atoms(system, source)
    charges = assign_native_charges(system, typing, source)
    print(
        json.dumps(
            {
                "sites": system.number_of_sites,
                "coverage": len(charges.payload["charges"]),
                "total_charge_e": charges.payload["data"]["total_e"],
                "source": source.identity,
                "capabilities": charges.payload["data"]["capabilities"],
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
