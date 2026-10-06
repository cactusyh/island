"""Explicit local source inventory/comparison; success is inspection only."""

import argparse
from pathlib import Path

from island.charge_references.records import pack
from island.workflows.storage import publish, read_json

from .variants import (
    PCFFVariantPolicy,
    PCFFVariantSelection,
    compare_pcff_sources,
    load_pcff_source_variant,
    resolve_pcff_variant_record,
    save_pcff_source_comparison,
    save_pcff_variant_resolution,
)


def declared_sources(path):
    data = read_json(Path(path))
    sources = []
    for row in data["sources"]:
        source = load_pcff_source_variant(
            row["path"], expected_sha256=row["sha256"], provenance=row["provenance"]
        )
        if source.selection != {"sha256": row["sha256"], "profile": row["profile"]}:
            raise ValueError("Wrong declared source profile")
        sources.append(source)
    if not sources:
        raise ValueError("Nonempty explicit source declaration required")
    return sources


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--declaration", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--family")
    p.add_argument("--namespace", default="cff91")
    p.add_argument("--types", nargs="+")
    p.add_argument("--sha256")
    p.add_argument("--profile")
    p.add_argument(
        "--fallback",
        nargs=2,
        action="append",
        default=[],
        metavar=("SHA256", "PROFILE"),
    )
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    if a.family and (not a.types or not a.sha256 or not a.profile):
        p.error("Queries require --types, --sha256 and --profile")
    if not a.family and any((a.types, a.sha256, a.profile, a.fallback)):
        p.error("Source selections apply only to an explicit --family query")
    a.output.mkdir(parents=True, exist_ok=False)
    sources = declared_sources(a.declaration)
    result = compare_pcff_sources(sources)
    save_pcff_source_comparison(result, a.output / "comparison.json")
    report = {
        "comparison_identity": result.identity,
        "inspection_complete": True,
        "full_source_complete": False,
        "operational_models_authorized": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    if a.family:
        policy = PCFFVariantPolicy(
            PCFFVariantSelection(a.sha256, a.profile),
            tuple(PCFFVariantSelection(*r) for r in a.fallback),
        )
        query = resolve_pcff_variant_record(
            sources,
            policy=policy,
            family=a.family,
            namespace=a.namespace,
            types=a.types,
        )
        save_pcff_variant_resolution(query, a.output / "query.json")
        report["query_identity"] = query.identity
        report["classification"] = query.payload["classification"]
    publish(a.output / "report.json", pack(report).encode())
    print(pack(report))
    return 1 if a.require_full_source else 0


if __name__ == "__main__":
    raise SystemExit(main())
