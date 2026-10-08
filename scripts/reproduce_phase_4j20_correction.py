"""Run the three reviewed reproductions against either the base or corrected code."""

import argparse
import json

from island.core import AtomSite, Coordinates, MolecularSystem, Topology
from island.crosslinking import ReactiveSiteRule, generate_crosslink_network
from island.exceptions import ValidationError
from island.graph import final_graph


def fragments(count=4):
    topology = Topology()
    for i in range(1, 2 * count + 1):
        topology.add_site(
            AtomSite(i, "C", 12, element="C", metadata={"reactive": True})
        )
    for i in range(1, 2 * count + 1, 2):
        topology.add_bond(i, i + 1, order=2)
    system = MolecularSystem(
        topology, Coordinates({i: (i, 0, 0) for i in topology.sites})
    )
    return system, final_graph(system)


def run(system, graph, element="C", seed=2026):
    return generate_crosslink_network(
        system,
        graph,
        ReactiveSiteRule(
            "test",
            "reactive",
            True,
            (element,),
            1,
            1.0,
            "software fixture",
            ("reproduction",),
        ),
        target_crosslinks=1,
        seed=seed,
        provenance="software fixture",
        evidence=("reproduction",),
    )


def reproduce():
    system, graph = fragments()
    first = run(system, graph)
    second = run(first.system, first.graph, seed=0)
    report = {
        "capacity": {
            "first": first.plan.payload["selected_crosslink_bonds"],
            "second": second.plan.payload["selected_crosslink_bonds"],
        }
    }
    system, _ = fragments(2)
    labels = {i: "original-A" if i < 3 else "original-B" for i in system.topology.sites}
    graph = final_graph(system, molecule_membership=labels)
    result = run(system, graph)
    report["membership"] = {
        "input": graph.payload["molecule_membership"],
        "output": result.graph.payload["molecule_membership"],
    }
    topology = Topology()
    for nitrogen in (1, 5):
        topology.add_site(
            AtomSite(nitrogen, "N", 14, element="N", metadata={"reactive": True})
        )
        for hydrogen in range(nitrogen + 1, nitrogen + 4):
            topology.add_site(AtomSite(hydrogen, "H", 1, element="H"))
            topology.add_bond(nitrogen, hydrogen, order=1)
    system = MolecularSystem(
        topology, Coordinates({i: (i, 0, 0) for i in topology.sites})
    )
    try:
        result = run(system, final_graph(system), "N")
    except ValidationError as error:
        report["neutral_nitrogen"] = {"rejected": str(error)}
    else:
        report["neutral_nitrogen"] = {
            "selected": result.plan.payload["selected_crosslink_bonds"],
            "output": [
                {
                    "id": i,
                    "formal_charge": result.system.topology.sites[i].formal_charge,
                    "degree": result.system.topology.degree(i),
                }
                for i in (1, 5)
            ],
        }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expect", choices=("defects", "corrected"), required=True)
    args = parser.parse_args()
    report = reproduce()
    print(json.dumps(report, indent=2))
    reused = bool(
        {i for pair in report["capacity"]["first"] for i in pair}
        & {i for pair in report["capacity"]["second"] for i in pair}
    )
    reset = report["membership"]["input"] != report["membership"]["output"]
    accepted_nitrogen = "output" in report["neutral_nitrogen"]
    expected = args.expect == "defects"
    return (
        0
        if all(value == expected for value in (reused, reset, accepted_nitrogen))
        else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
