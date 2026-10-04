"""Explicit-H chain -> pinned PCFF model -> bounded session-backed minimum.

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
from island.minimization import MinimizationOptions, minimize_geometry


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
    options = MinimizationOptions(
        force_tolerance=0.1, max_iterations=5000, max_evaluations=10000
    )
    with evaluator.open_session() as session:
        result = minimize_geometry(system, session, options)
    print(
        json.dumps(
            {
                "converged": result.converged,
                "verified": result.final_evaluation_verified,
                "termination_reason": result.termination_reason,
                "initial_energy_kj_mol": result.initial_energy,
                "final_energy_kj_mol": result.final_energy,
                "fmax_kj_mol_angstrom": result.final_fmax,
                "evaluations": result.evaluations,
                "iterations": result.iterations,
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
            indent=2,
        )
    )
    if not (result.converged and result.final_evaluation_verified):
        return 1
    optimized = result.to_system(system)
    optimized.validate()
    specification.validate_integrity(optimized)
    print("Validated copied minimum; caller unchanged. This is not equilibration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
