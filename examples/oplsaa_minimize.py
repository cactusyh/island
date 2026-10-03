"""Short-chain local minimum; completion does not establish equilibration."""

import argparse
import json

from island.builders import build_linear_polymer
from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.forcefields.oplsaa import load_oplsaa_source, parameterize_oplsaa
from island.minimization import MinimizationOptions, minimize_geometry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xml", required=True, help="Installed checksum-pinned OPLS XML"
    )
    parser.add_argument("--psmiles", default="[*:1]CC[*:2]")
    parser.add_argument("--dp", type=int, default=3)
    args = parser.parse_args()
    system = build_linear_polymer(
        args.psmiles,
        dp=args.dp,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    source = load_oplsaa_source(args.xml)
    parameters = parameterize_oplsaa(system, source)
    evaluator = OPLSSinglePointEvaluator(system, parameters, source)
    options = MinimizationOptions(
        force_tolerance=0.1, max_iterations=5000, max_evaluations=10000
    )
    with evaluator.open_session() as session:
        result = minimize_geometry(system, session, options)
    print(
        json.dumps(
            {
                "converged": result.converged,
                "termination_reason": result.termination_reason,
                "independently_verified": result.final_evaluation_verified,
                "initial_energy_kj_mol": result.initial_energy,
                "final_energy_kj_mol": result.final_energy,
                "final_fmax_kj_mol_angstrom": result.final_fmax,
                "evaluations": result.evaluations,
                "iterations": result.iterations,
                "parameter_identity": evaluator.parameter_fingerprint,
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
    parameters.validate_integrity(optimized, source)
    print(
        "Validated copied system:", optimized.number_of_sites, "sites; caller unchanged"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
