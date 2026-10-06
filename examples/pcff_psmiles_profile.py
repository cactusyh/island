"""Final PSMILES graph -> exact operational PCFF profile -> native single point.

Only linked six-carbon benzenoid rings (and methyl ends) are authorized here.
No CAR/MDF, external converter, embedding during assignment, or QM is involved.
Construction generates coordinates; charges and interactions use the final graph.
"""

import argparse

from island.builders import build_linear_polymer
from island.forcefields import PCFFOptions, create_evaluator
from island.forcefields.pcff import prepare_pcff_polymer, select_pcff_profile
from island.forcefields.pcff.registry import LINKED_BENZENOID
from island.forcefields.pcff.source import PIN


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--dp", type=int, default=2)
    args = parser.parse_args()
    system = build_linear_polymer("[*:1]c1ccc([*:2])cc1", dp=args.dp, random_seed=2026)
    selection = select_pcff_profile(LINKED_BENZENOID, sha256=PIN["sha256"])
    options = PCFFOptions(
        args.source,
        (0, 0, 1),
        (0, 0, 1),
        typing_profile="island_pcff_source_graph_v1",
        source_profile=selection,
    )
    final = prepare_pcff_polymer(system, options)
    evaluator = create_evaluator(system, final.prepared)
    with evaluator.open_session() as session:
        result = session.evaluate(
            {i: system.coordinates.get(i) for i in system.topology.sites}
        )
    print("profile:", selection.name, "source:", selection.sha256)
    print("sites:", system.number_of_sites, "energy kJ/mol:", result.potential_energy)
    print("components:", dict(result.energy_components))
    print("production_validated=False; simulation_readiness=not_established")


if __name__ == "__main__":
    main()
