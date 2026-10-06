"""Explicit opt-in J3 chemistry and native charge/parameter coverage.

Run with a verified external pcff.frc. No implicit preparation fallback.
"""

import argparse

from island.chemistry import from_smiles
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    assign_pcff_source_types,
    define_pcff_model,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.domains import PROFILE_NAME
from island.forcefields.pcff.fallbacks import DOMAIN_POLICY


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--smiles", default="S")
    p.add_argument("--explicit-checked", action="store_true")
    args = p.parse_args()
    system = from_smiles(args.smiles, random_seed=2026)
    source = load_pcff_source(args.source)
    typed = type_pcff_atoms(system, source, profile=PROFILE_NAME)
    print("Typing:", typed.complete, typed.assignments, typed.payload["diagnostics"])
    if not typed.complete:
        return 1
    if args.explicit_checked:
        typed = assign_pcff_source_types(
            system,
            source,
            typed.assignments,
            profile=PROFILE_NAME,
            provenance="Example labels from validated graph typing; not independent chemical evidence",
        )
    charges = assign_automatic_pcff_charges(
        system, typed, resolution_policy=DOMAIN_POLICY
    )
    print(
        "Native charges:",
        charges.complete,
        charges.payload["native_charge_record"]["diagnostics"],
    )
    if not charges.complete:
        return 1
    assignment = assign_pcff_parameters(
        system, typed, charges, resolution_policy=DOMAIN_POLICY
    )
    model = define_pcff_model(
        assignment, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    print(
        "Model:",
        model.identity,
        model.payload["model_definition_complete"],
        model.payload["diagnostics"],
    )
    print("production_validated=False; simulation_readiness=not_established")
    return 0 if model.payload["model_definition_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
