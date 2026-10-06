"""Explicit source inventories and separate coverage stages; never authorization."""

import argparse
from pathlib import Path

from island.forcefields.pcff import load_pcff_source, load_pcff_source_variant
from island.forcefields.pcff.registry import list_pcff_profiles
from island.forcefields.pcff.registry_cli import pcff_profile_coverage
from island.workflows import storage


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--variants-declaration", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    native = pcff_profile_coverage(load_pcff_source(a.source))
    declaration = storage.read_json(a.variants_declaration)
    variants = []
    seen = set()
    for supplied in declaration["sources"]:
        if (
            supplied["sha256"] == native["source"]["sha256"]
            or supplied["sha256"] in seen
        ):
            continue
        variant = load_pcff_source_variant(
            supplied["path"],
            expected_sha256=supplied["sha256"],
            provenance=supplied["provenance"],
        )
        inv = variant.payload
        seen.add(supplied["sha256"])
        rows = inv["rows"]
        types = []
        for label in inv["atom_labels"]:
            atom = [
                r
                for r in rows
                if r["content"]["section"] == "atom_types"
                and r["semantic"]["type"] == label
            ]
            eq = [
                r["id"]
                for r in rows
                if r["content"]["section"] == "equivalence"
                and r["semantic"]["type"] == label
            ]
            auto = [
                r["id"]
                for r in rows
                if r["content"]["section"] == "auto_equivalence"
                and r["semantic"]["type"] == label
            ]
            types.append(
                {
                    "type": label,
                    "source_rows": [r["id"] for r in atom],
                    "inventory_presence": True,
                    "ordinary_equivalence_coverage": eq,
                    "automatic_equivalence_coverage": auto,
                    **{
                        stage: "candidate_not_authorized"
                        for stage in (
                            "typing_status",
                            "native_charge_status",
                            "lower_order_parameter_status",
                            "cross_term_status",
                            "executable_interpretation_status",
                            "independent_energy_force_status",
                            "operational_profile_status",
                            "final_polymer_graph_coverage",
                        )
                    },
                    "exact_blocking_reason": "Comparison source only; independently audit typing/charges/forms before any runtime authorization",
                }
            )
        variants.append(
            {
                "source": inv["selection"],
                "provenance": inv["provenance"],
                "counts": inv["counts"],
                "type_rows": types,
                "families": [
                    {
                        "name": section["name"],
                        "namespace": section["namespace"],
                        "record_ids": section["record_ids"],
                        "record_count": section["record_count"],
                        "runtime_status": "candidate_not_authorized",
                    }
                    for section in inv["sections"]
                ],
            }
        )
    report = {
        "schema": "island_j9_all_declared_sources_coverage_v1",
        "native": native,
        "comparison_only_sources": variants,
        "registry": list_pcff_profiles(),
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    storage.publish(a.output, storage.json_bytes(report))
    print(
        {
            "native_labels": len(native["type_rows"]),
            "candidate_labels": [len(v["type_rows"]) for v in variants],
            "full_source_complete": False,
        }
    )
    return 1 if a.require_full_source else 0


if __name__ == "__main__":
    raise SystemExit(main())
