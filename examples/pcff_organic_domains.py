"""Explicit v4 typing/charge/model coverage; no claim of complete PCFF coverage."""

import argparse
import json

from island.chemistry import from_smiles
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    load_pcff_source,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import DOMAIN_POLICY
from island.forcefields.pcff.organic_domains import PROFILE_NAME


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", required=True, help="Explicit pinned local pcff.frc"
    )
    parser.add_argument("--smiles", default="[H]Cl")
    parser.add_argument("--evaluate", action="store_true")
    args = parser.parse_args()
    system = from_smiles(args.smiles, random_seed=2026)
    source = load_pcff_source(args.source)
    typing = type_pcff_atoms(system, source, profile=PROFILE_NAME)
    report = {
        "typing_complete": typing.complete,
        "typing": typing.payload["assignments"],
        "diagnostics": typing.payload["diagnostics"],
        "source": source.identity,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    if typing.complete:
        charges = assign_automatic_pcff_charges(
            system, typing, resolution_policy=DOMAIN_POLICY
        )
        report.update(
            charges_complete=charges.complete,
            charge_diagnostics=charges.payload["native_charge_record"]["diagnostics"],
        )
    print(json.dumps(report, indent=2))
    if args.evaluate:
        prepared = prepare_forcefield(
            system,
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    source_path=args.source,
                    typing_profile=PROFILE_NAME,
                    resolution_policy=DOMAIN_POLICY,
                    lj=(0, 0, 1),
                    coulomb=(0, 0, 1),
                ),
            ),
        )
        result = create_evaluator(system, prepared).evaluate_fresh()
        print(
            "Energy (kJ/mol):",
            result.potential_energy,
            "; model:",
            result.model_fingerprint,
        )


if __name__ == "__main__":
    main()
