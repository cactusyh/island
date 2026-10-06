"""Inspect explicitly supplied sources and query one declared hash/profile."""

import argparse
from pathlib import Path

from island.forcefields.pcff import (
    PCFFVariantPolicy,
    PCFFVariantSelection,
    compare_pcff_sources,
    resolve_pcff_variant_record,
    save_pcff_source_comparison,
    save_pcff_variant_resolution,
)
from island.forcefields.pcff.variants_cli import declared_sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--family", default="bond_increments")
    parser.add_argument("--namespace", default="cff91_auto")
    parser.add_argument("--types", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = declared_sources(args.sources)
    comparison = compare_pcff_sources(sources)
    policy = PCFFVariantPolicy(PCFFVariantSelection(args.sha256, args.profile))
    query = resolve_pcff_variant_record(
        sources,
        policy=policy,
        family=args.family,
        namespace=args.namespace,
        types=args.types,
    )
    args.output.mkdir(parents=True, exist_ok=False)
    save_pcff_source_comparison(comparison, args.output / "comparison.json")
    save_pcff_variant_resolution(query, args.output / "query.json")
    for variant in sources:
        p = variant.payload
        print(
            p["selection"],
            "atom labels:",
            len(p["atom_labels"]),
            "increments:",
            p["counts"].get("bond_increments", 0),
        )
    print("Comparison identity:", comparison.identity)
    print(
        "Query:",
        query.payload["classification"],
        "source:",
        query.payload["selected_source"],
    )
    print(
        "Operational model assembly authorized:",
        query.payload["model_assembly_authorized"],
    )
    print("production_validated=False; simulation_readiness=not_established")


if __name__ == "__main__":
    main()
