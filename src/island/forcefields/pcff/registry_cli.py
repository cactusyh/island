"""Inspect source/profile coverage or prepare a final PSMILES-built graph."""

import argparse
import json

from .registry import inspect_pcff_profile, list_pcff_profiles, select_pcff_profile
from .source import FLAGS, PIN, boundary, load_pcff_source


@boundary
def pcff_profile_coverage(source):
    from .domain_coverage import pcff_domain_coverage

    ledger = pcff_domain_coverage(source, profile="island_pcff_source_graph_v4")
    for row in ledger["type_rows"]:
        row.update(
            inventory_presence=True,
            ordinary_equivalence_coverage=[
                r["id"]
                for r in row["equivalence_rows"]
                if r["section"] == "equivalence"
            ],
            automatic_equivalence_coverage=[
                r["id"]
                for r in row["equivalence_rows"]
                if r["section"] == "auto_equivalence"
            ],
            typing_status=row["graph_predicate_status"],
            native_charge_status="graph-dependent; no global completion claim",
            lower_order_parameter_status="source candidates only; requires oriented final-graph assignment",
            cross_term_status="interaction-dependent; missing remains missing",
            executable_interpretation_status=row["executable_evaluator_status"],
            independent_energy_force_status=row["independent_numerical_verification"],
            operational_profile_status=[
                p["name"]
                for p in list_pcff_profiles()
                if p["state"] == "operational"
                and row["type"] in p["definition"]["authorized_atom_labels"]
            ],
            final_polymer_graph_coverage="bounded linked-benzenoid profile only"
            if row["type"] in ("cp", "hc", "c3")
            else "not independently established",
            msi2lmp_chemical_perception="not implemented by converter; requires supplied types",
            exact_blocking_reason=row["next_action"],
        )
    ledger.update(
        schema="island_pcff_registry_coverage_v1",
        registry=list_pcff_profiles(),
        msi2lmp_capabilities={
            "supplied_type_conversion": True,
            "automatic_typing": False,
            "bond_increments": False,
            "auto_equivalence_supplementation": False,
            "runtime_dependency": False,
        },
        full_source_complete=False,
        **FLAGS,
    )
    return ledger


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile")
    parser.add_argument("--source")
    parser.add_argument("--psmiles")
    parser.add_argument("--dp", type=int, default=2)
    parser.add_argument("--diagnostic", action="store_true")
    parser.add_argument("--require-full-source", action="store_true")
    args = parser.parse_args(argv)
    if args.psmiles:
        from island.builders import build_linear_polymer
        from island.forcefields import PCFFOptions

        from .polymer import prepare_pcff_polymer

        system = build_linear_polymer(args.psmiles, dp=args.dp)
        selection = select_pcff_profile(args.profile, sha256=PIN["sha256"])
        options = PCFFOptions(
            args.source,
            (0, 0, 1),
            (0, 0, 1),
            typing_profile=selection.profile().payload["typing_profile"],
            source_profile=selection,
        )
        result = prepare_pcff_polymer(
            system, options, mode="diagnostic" if args.diagnostic else "strict"
        )
        report = (
            result
            if args.diagnostic
            else {
                "prepared": result.prepared.metadata,
                "final_graph_evidence_identity": result.identity,
            }
        )
    elif args.source:
        report = pcff_profile_coverage(load_pcff_source(args.source))
    else:
        report = (
            inspect_pcff_profile(args.profile) if args.profile else list_pcff_profiles()
        )
    print(json.dumps(report, sort_keys=True, indent=2))
    return 1 if args.require_full_source else 0


if __name__ == "__main__":
    raise SystemExit(main())
