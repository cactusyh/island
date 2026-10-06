"""Final PSMILES graph -> explicitly opt-in converter BB13 policy -> evaluation.

This is converter compatibility, not strict FRC completeness or scientific validation.
CAR/MDF, msi2lmp and LAMMPS are not runtime dependencies.
"""

import argparse

from island.builders import build_linear_polymer
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--psmiles", default="[*:1]CCO[*:2]")
    p.add_argument("--dp", type=int, default=3)
    a = p.parse_args()
    system = build_linear_polymer(a.psmiles, dp=a.dp, random_seed=2026)
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                a.source,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile="island_pcff_source_graph_v4",
                resolution_policy=COMPATIBILITY_POLICY,
            ),
        ),
    )
    model = prepared.native_result.payload
    with create_evaluator(system, prepared).open_session() as session:
        result = session.evaluate(
            {i: system.coordinates.get(i) for i in system.topology.sites}
        )
    print(
        {
            "native_identity": prepared.native_result.identity,
            "policy": COMPATIBILITY_POLICY,
            "raw_source_coverage_complete": model["raw_parameter_coverage_complete"],
            "term_origins": model["term_origins"],
            "energy_kj_mol": result.potential_energy,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


if __name__ == "__main__":
    main()
