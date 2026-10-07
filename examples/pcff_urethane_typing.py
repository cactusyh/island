"""PSMILES final graph -> opt-in v6 typing, charges and honest model coverage.

The default urethane has missing source couplings and exits 1. No converter or
CAR/MDF is required. Charge completeness does not authorize evaluation.
"""

import argparse
from collections import Counter

from island.builders import build_linear_polymer
from island.forcefields import ForceFieldRequest, PCFFOptions, prepare_forcefield
from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY, MSI_POLICY
from island.forcefields.pcff.urethane_domains import PROFILE_NAME


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--psmiles", default="[*:1]CCOC(=O)NCC[*:2]")
    parser.add_argument("--dp", type=int, default=2)
    parser.add_argument("--converter-compatibility", action="store_true")
    args = parser.parse_args()
    system = build_linear_polymer(args.psmiles, dp=args.dp, random_seed=2026)
    source = load_pcff_source(args.source)
    policy = COMPATIBILITY_POLICY if args.converter_compatibility else MSI_POLICY
    typing = type_pcff_atoms(system, source, profile=PROFILE_NAME)
    if not typing.complete:
        print({"typing_complete": False, "diagnostics": typing.payload["diagnostics"]})
        return 1
    report = inspect_pcff_operational_support(
        system,
        typing,
        resolution_policy=policy,
        special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    print(
        {
            "typing_profile": PROFILE_NAME,
            "source": source.identity,
            "labels": dict(Counter(typing.assignments.values())),
            "resolution_policy": policy,
            **{
                k: report[k]
                for k in ("charges_complete", "model_complete", "model_diagnostics")
            },
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    if not report["model_complete"]:
        return 1
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                args.source,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile=PROFILE_NAME,
                resolution_policy=policy,
            ),
        ),
    )
    print({"prepared_identity": prepared.identity})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
