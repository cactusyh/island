"""Inspect a retained source-typed system without preparing an incomplete model."""

import argparse
import json
from collections import Counter
from pathlib import Path

from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.workflows.bundle import system_from
from island.workflows.storage import decode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--system", type=Path, required=True)
    args = parser.parse_args()
    system = system_from(decode(json.loads(args.system.read_text())))
    source = load_pcff_source(args.source)
    typing = type_pcff_atoms(system, source, profile="island_pcff_source_graph_v4")
    if not typing.complete:
        print(typing.payload["diagnostics"])
        return 1
    report = inspect_pcff_operational_support(
        system,
        typing,
        resolution_policy="island_pcff_positional_fallbacks_v2",
        special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    print("Native charges complete:", report["charges_complete"])
    print("Operational model complete:", report["model_complete"])
    print(
        "Structural queries:",
        dict(Counter(q["resolution"]["status"] for q in report["interaction_queries"])),
    )
    print(
        "Diagnostic only. No evaluator is constructed and no scientific record is published."
    )
    return 0  # Inspection success does not assert scientific/operational acceptance.


if __name__ == "__main__":
    raise SystemExit(main())
