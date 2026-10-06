"""Explicit profile -> native PCFF -> single point; source library stays local."""

import argparse
from pathlib import Path

from island.chemistry import from_smiles
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.pcff import PCFFOperationalSelection
from island.forcefields.pcff.operational_profile import DEFINITION, NAME


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    a = p.parse_args()
    system = from_smiles("c1ccccc1", random_seed=2026)
    selected = PCFFOperationalSelection(NAME, DEFINITION["source_sha256"])
    request = ForceFieldRequest(
        "pcff",
        PCFFOptions(
            a.source,
            (0, 0, 1),
            (0, 0, 1),
            typing_profile="island_pcff_source_graph_v1",
            source_profile=selected,
        ),
    )
    prepared = prepare_forcefield(system, request)
    evaluator = create_evaluator(system, prepared)
    result = evaluator.evaluate_fresh()
    print("Selected source/profile:", selected.sha256, selected.name)
    print(
        "Prepared/native identities:",
        prepared.identity,
        prepared.native_result.identity,
    )
    print("Energy (kJ/mol):", result.potential_energy)
    print(
        "Maximum force (kJ/(mol*angstrom)):",
        max(sum(x * x for x in f) ** 0.5 for f in result.forces.values()),
    )
    print(
        "Scope: connected neutral explicit-H C6H6 aromatic cycle; zero Wilson equilibrium"
    )
    print("production_validated=False; simulation_readiness=not_established")


if __name__ == "__main__":
    main()
