"""Deterministic, offline final-graph force-field inspection; never runs MD."""

import argparse
import json
from pathlib import Path

from island.forcefields.selection import capability_registry, prepare_polymer_forcefield


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    graph = parser.add_mutually_exclusive_group()
    graph.add_argument("--smiles")
    graph.add_argument("--psmiles")
    parser.add_argument("--dp", type=int)
    parser.add_argument("--head-end-group")
    parser.add_argument("--tail-end-group")
    parser.add_argument("--forcefield", choices=("pcff", "opls", "gaff", "gaff2"))
    parser.add_argument("--profile")
    parser.add_argument("--charge-method")
    parser.add_argument(
        "--resolution-policy",
        default="strict_source",
        choices=(
            "strict_source",
            "derived_from_declared_formula",
            "fit_from_external_reference",
        ),
    )
    parser.add_argument("--source", type=Path)
    parser.add_argument("--sha256")
    parser.add_argument("--amberhome", type=Path)
    parser.add_argument("--lj", type=float, nargs=3)
    parser.add_argument("--coulomb", type=float, nargs=3)
    parser.add_argument("--oxygen-anchor", type=int, action="append", default=[])
    parser.add_argument("--registry", action="store_true")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.registry:
        data, code = {"profiles": capability_registry()}, 0
    else:
        if not (args.smiles or args.psmiles) or not args.forcefield:
            parser.error("a final-graph SMILES or PSMILES and forcefield are required")
        try:
            if args.psmiles:
                from island.builders import build_linear_polymer

                if (
                    args.dp is None
                    or not args.head_end_group
                    or not args.tail_end_group
                ):
                    parser.error(
                        "PSMILES requires dp and explicit head/tail end groups"
                    )
                system = build_linear_polymer(
                    args.psmiles,
                    dp=args.dp,
                    head_end_group=args.head_end_group,
                    tail_end_group=args.tail_end_group,
                    add_hydrogens=True,
                    random_seed=2026,
                )
            else:
                from island.chemistry import from_smiles

                system = from_smiles(args.smiles, random_seed=2026)
            from island.forcefields.selection.registry import (
                PDMS,
                USER_OPLS,
                inspect_profile,
            )

            options = None
            if args.profile and not inspect_profile(args.profile).get(
                "selection_adapter", True
            ):
                pass  # The common request boundary returns the exact policy diagnostic.
            elif args.profile == PDMS:
                from island.forcefields.pdms_opls import PDMSOPLSOptions

                if args.source and args.sha256:
                    options = PDMSOPLSOptions(
                        args.source,
                        args.sha256,
                        e_oxygen_site_ids=tuple(args.oxygen_anchor),
                    )
            elif args.profile == USER_OPLS:
                from island.forcefields.opls_external import ExternalOPLSOptions

                if args.source and args.sha256:
                    options = ExternalOPLSOptions(args.source, args.sha256)
            elif args.forcefield == "pcff" and args.source:
                from island.forcefields.preparation import PCFFOptions

                if args.lj is None or args.coulomb is None:
                    parser.error(
                        "PCFF requires explicit LJ and Coulomb special-pair tuples"
                    )
                options = PCFFOptions(
                    args.source,
                    tuple(args.lj),
                    tuple(args.coulomb),
                    typing_profile=args.profile or "island_pcff_acyclic_cho_v1",
                )
            elif args.forcefield == "opls" and args.source:
                from island.forcefields.preparation import OPLSOptions

                options = OPLSOptions(args.source)
            elif args.forcefield in ("gaff", "gaff2") and args.amberhome:
                from island.forcefields.ambertools import AmberToolsOptions

                if args.charge_method == "am1bcc":
                    options = AmberToolsOptions(
                        force_field=args.forcefield,
                        charge_method="am1bcc",
                        amberhome=args.amberhome,
                    )
            result = prepare_polymer_forcefield(
                system,
                args.forcefield,
                source_profile=args.profile,
                charge_method=args.charge_method,
                resolution_policy=args.resolution_policy,
                options=options,
            )
            data = result.payload
            code = 0 if data["stages"]["executable"] else 1
        except Exception as error:  # noqa: BLE001 -- classified public boundary
            data = {
                "status": "typing_ambiguous",
                "diagnostics": [
                    {
                        "stage": "input_construction",
                        "exception": type(error).__name__,
                        "reason": str(error),
                    }
                ],
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
            code = 1
    args.json.write_text(
        json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
