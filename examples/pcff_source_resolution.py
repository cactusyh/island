"""PSMILES final graph -> explicit guarded PCFF source policy -> single point.

No CAR/MDF, converter, or LAMMPS is needed. An unresolved type, charge or
coupling prevents preparation. This example is not full-source acceptance.
"""

import argparse

from island.builders import build_linear_polymer
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.pcff.fallbacks import MSI_POLICY


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, help="Exact pinned external pcff.frc")
    a = p.parse_args()
    system = build_linear_polymer("[*:1]c1ccc([*:2])cc1", dp=2, random_seed=2026)
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                a.source,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile="island_pcff_source_graph_v4",
                resolution_policy=MSI_POLICY,
            ),
        ),
    )
    evaluator = create_evaluator(system, prepared)
    coordinates = {i: system.coordinates.get(i) for i in system.topology.sites}
    with evaluator.open_session() as session:
        result = session.evaluate(coordinates)
        fresh = session.evaluate_fresh(coordinates)
    assert result.evaluation_fingerprint == fresh.evaluation_fingerprint
    print("source:", prepared.native_result.assignment.source.identity["sha256"])
    print("resolver:", MSI_POLICY)
    print("model:", prepared.native_result.identity)
    print("energy (kJ/mol):", result.potential_energy)
    print("force coverage:", len(result.forces))
    print("production_validated=False; simulation_readiness=not_established")


if __name__ == "__main__":
    main()
