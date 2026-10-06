"""Explicit source-fallback policy: inspect coverage, evaluate only complete models."""

import argparse

from island.chemistry import from_smiles
from island.forcefields import adopt_forcefield, create_evaluator
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import POLICY


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--smiles", default="ClCl")
    args = p.parse_args()
    system = from_smiles(args.smiles, random_seed=2026)
    source = load_pcff_source(args.source)
    typing = type_pcff_atoms(system, source, profile="island_pcff_source_graph_v2")
    print("Typing:", typing.payload["coverage"])
    if not typing.complete:
        return 1
    charges = assign_automatic_pcff_charges(system, typing, resolution_policy=POLICY)
    print("Native charges complete:", charges.complete)
    if not charges.complete:
        print(charges.payload["native_charge_record"]["diagnostics"])
        return 1
    parameters = assign_pcff_parameters(
        system, typing, charges, resolution_policy=POLICY
    )
    model = define_pcff_model(
        parameters, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    print("Source family coverage:", parameters.payload["coverage"])
    print("Model complete:", model.payload["model_definition_complete"])
    if not model.payload["model_definition_complete"]:
        print(model.payload["diagnostics"])
        return 1
    prepared = adopt_forcefield(system, "pcff", model)
    with create_evaluator(system, prepared).open_session() as session:
        result = session.evaluate()
        print("Energy (kJ/mol):", result.potential_energy)
        print("Components:", dict(result.energy_components))
    print(
        "Experimental source-profile fidelity only; scientific readiness not established."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
