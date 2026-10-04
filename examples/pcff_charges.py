"""Illustrative manual typing, no automatic typing/embedding or energy model.

Run: python examples/pcff_charges.py /local/pinned/pcff.frc
"""

import argparse
import json

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    assign_pcff_charges,
    assign_pcff_types,
    load_pcff_source,
)


def molecule(case="ethanol", *, reverse=False):
    """Explicit hand-specified CHO graphs and illustrative source-comment-based types."""
    definitions = {
        "ethane": (["C", "C"], [(0, 1)], [3, 3], ["c3", "c3"]),
        "ethanol": (["C", "C", "O"], [(0, 1), (1, 2)], [3, 2, 1], ["c3", "c2", "oh"]),
        "dimethyl_ether": (
            ["C", "O", "C"],
            [(0, 1), (1, 2)],
            [3, 0, 3],
            ["c3", "oc", "c3"],
        ),
    }
    elements, bonds, hcounts, types = definitions[case]
    elements, bonds, types = list(elements), list(bonds), list(types)
    for parent, count in enumerate(hcounts):
        for _ in range(count):
            bonds.append((parent, len(elements)))
            types.append("ho" if elements[parent] == "O" else "hc")
            elements.append("H")
    ids = [10 + i * 7 for i in range(len(elements))]
    topology = Topology()
    values = {"C": (6, 12.011), "H": (1, 1.008), "O": (8, 15.999)}
    for i in reversed(range(len(ids))) if reverse else range(len(ids)):
        z, mass = values[elements[i]]
        topology.add_site(
            AtomSite(
                ids[i], f"{elements[i]}{i}", mass, element=elements[i], atomic_number=z
            )
        )
    for a, b in reversed(bonds) if reverse else bonds:
        topology.add_bond(
            ids[b] if reverse else ids[a], ids[a] if reverse else ids[b], order=1
        )
    return MolecularSystem(topology=topology, coordinates=Coordinates()), dict(
        zip(ids, types)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source")
    args = parser.parse_args()
    source = load_pcff_source(args.source)
    system, labels = molecule()
    typing = assign_pcff_types(
        system, source, labels, provenance="Illustrative manual source-comment typing"
    )
    result = assign_pcff_charges(system, typing)
    print(
        json.dumps(
            {
                "source": source.identity,
                "charges_e": result.charges,
                "total_e": result.payload["total_charge"],
                "limitations": "No automatic typing or energy model",
            }
        )
    )
    try:
        assign_pcff_types(
            system,
            source,
            {k: v for k, v in labels.items() if k != min(labels)},
            provenance="Intentionally incomplete example",
        )
    except PCFFError as error:
        print("Expected missing coverage rejection:", error)
    return 0 if result.complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
