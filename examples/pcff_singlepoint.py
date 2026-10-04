"""Explicit-H chain -> pinned native PCFF assignment -> finite single point.

This is not a simulation-ready snapshot or an accuracy/equilibration claim.
"""

import argparse
import json

from island.builders import build_linear_polymer
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, help="Pinned local pcff.frc")
    parser.add_argument("--psmiles", default="[*:1]CC[*:2]")
    parser.add_argument("--dp", type=int, default=3)
    parser.add_argument("--lj", type=float, nargs=3, required=True)
    parser.add_argument("--coulomb", type=float, nargs=3, required=True)
    args = parser.parse_args()
    system = build_linear_polymer(
        args.psmiles,
        dp=args.dp,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    source = load_pcff_source(args.source)
    typing = type_pcff_atoms(system, source)
    charges = assign_automatic_pcff_charges(system, typing)
    assignment = assign_pcff_parameters(system, typing, charges)
    specification = define_pcff_model(
        assignment, special_pairs=special_pair_policy(lj=args.lj, coulomb=args.coulomb)
    )
    evaluator = PCFFSinglePointEvaluator(system, specification)
    result = evaluator.evaluate_fresh()
    print(
        json.dumps(
            {
                "potential_energy_kj_mol": result.potential_energy,
                "components_kj_mol": dict(result.energy_components),
                "forces_kj_mol_angstrom": dict(result.forces),
                "model_fingerprint": result.model_fingerprint,
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
