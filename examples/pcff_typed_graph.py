"""PSMILES polyethylene with explicitly selected source-family c/h types.

Requires only the pinned local FRC. No CAR/MDF or converter at runtime.
The FRC describes c as generic SP3 carbon and h as H bound to C/Si/H.
This bounded example does not establish scientific suitability.
"""

import argparse

from island.builders import build_linear_polymer
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    assign_typed_pcff_charges,
    bind_pcff_types,
    load_pcff_source,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    args = parser.parse_args()
    system = build_linear_polymer("[*:1]CC[*:2]", dp=3, random_seed=2026)
    source = load_pcff_source(args.source)
    # Explicit choice for this saturated C/H graph, not automatic typing output.
    labels = {
        i: {"C": "c", "H": "h"}[a.element] for i, a in system.topology.sites.items()
    }
    typed = bind_pcff_types(
        system,
        source,
        labels,
        provenance="FRC generic SP3 carbon and carbon-bound hydrogen descriptions",
        evidence_references=["pcff.frc atom_types c", "pcff.frc atom_types h"],
    )
    charges = assign_typed_pcff_charges(
        system, typed, resolution_policy=COMPATIBILITY_POLICY
    )
    prepared = prepare_forcefield(
        system,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                args.source,
                lj=(0, 0, 1),
                coulomb=(0, 0, 1),
                resolution_policy=COMPATIBILITY_POLICY,
                typed_graph=typed,
                graph_charges=charges,
            ),
        ),
    )
    result = create_evaluator(system, prepared).evaluate_fresh()
    print(
        {
            "typing_origin": typed.payload["origin"],
            "charge_origin": charges.payload["origin"],
            "energy_kj_mol": result.potential_energy,
            "prepared_identity": prepared.identity,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


if __name__ == "__main__":
    main()
