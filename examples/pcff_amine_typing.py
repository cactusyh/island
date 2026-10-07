"""PSMILES final graph -> opt-in ref-1 amine typing -> native charges/model.

No CAR/MDF or converter is used. Strict source v3 remains incomplete for these
chains. The default here explicitly selects existing converter compatibility v4.
"""

import argparse
from collections import Counter

from island.builders import build_linear_polymer
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.amine_domains import PROFILE_NAME
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY, MSI_POLICY


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--psmiles", default="[*:1]CCN[*:2]")
    p.add_argument("--dp", type=int, default=2)
    p.add_argument("--strict-source", action="store_true")
    a = p.parse_args()
    s = build_linear_polymer(a.psmiles, dp=a.dp, random_seed=2026)
    source = load_pcff_source(a.source)
    policy = MSI_POLICY if a.strict_source else COMPATIBILITY_POLICY
    typing = type_pcff_atoms(s, source, profile=PROFILE_NAME)
    if not typing.complete:
        print(typing.payload["diagnostics"])
        return 1
    q = assign_automatic_pcff_charges(s, typing, resolution_policy=policy)
    diagnostic = inspect_pcff_operational_support(
        s,
        typing,
        resolution_policy=policy,
        special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    print(
        {
            "typing": dict(Counter(typing.assignments.values())),
            "charges_complete": q.complete,
            "total_charge_e": sum(q.charges.values()) if q.complete else None,
            "model_complete": diagnostic["model_complete"],
            "typing_profile": PROFILE_NAME,
            "resolution_policy": policy,
            "diagnostics": diagnostic["model_diagnostics"],
        }
    )
    if not diagnostic["model_complete"]:
        return 1
    prepared = prepare_forcefield(
        s,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                a.source,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile=PROFILE_NAME,
                resolution_policy=policy,
            ),
        ),
    )
    with create_evaluator(s, prepared).open_session() as session:
        evaluated = session.evaluate(
            {i: s.coordinates.get(i) for i in s.topology.sites}
        )
    print(
        {
            "energy_kj_mol": evaluated.potential_energy,
            "native_identity": prepared.native_result.identity,
            "raw_source_coverage_complete": prepared.native_result.payload[
                "raw_parameter_coverage_complete"
            ],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
