"""Automatic/explicit ring typing -> native model -> unified finite evaluator.

J1 implements a partial full-source profile; inspect unresolved interactions.
No simulation-readiness or full FRC coverage is asserted.
"""

import argparse

from island.chemistry import from_smiles
from island.forcefields import adopt_forcefield, create_evaluator
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    assign_pcff_source_types,
    define_pcff_model,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument(
        "--explicit",
        action="store_true",
        help="Use independently declared cyclohexane labels",
    )
    args = parser.parse_args()
    system = from_smiles("C1CCCCC1", random_seed=2026)
    source = load_pcff_source(args.source)
    if args.explicit:
        labels = {
            i: ("c2" if a.element == "C" else "hc")
            for i, a in system.topology.sites.items()
        }
        typing = assign_pcff_source_types(
            system,
            source,
            labels,
            provenance="Manual cyclohexane fixture: six sp3 CH2 sites and twelve C-bound H; source lines 71/102",
        )
    else:
        typing = type_pcff_atoms(system, source, profile="island_pcff_source_graph_v1")
    charges = assign_automatic_pcff_charges(system, typing)
    parameters = assign_pcff_parameters(system, typing, charges)
    model = define_pcff_model(
        parameters, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    prepared = adopt_forcefield(system, "pcff", model)
    with create_evaluator(system, prepared).open_session() as session:
        result = session.evaluate()
    print("Typing origin:", typing.payload["origin"])
    print("Charge:", sum(charges.charges.values()), "e")
    print("Energy:", result.potential_energy, "kJ/mol")
    print("Model:", result.model_fingerprint)
    print("production_validated=False; simulation_readiness=not_established")


if __name__ == "__main__":
    main()
