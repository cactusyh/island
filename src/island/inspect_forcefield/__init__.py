"""Canonical OPLS inspection entry point, with the historical expert CLI retained.

The additive package keeps the J46 module bytes intact for historical manifests.
Non-OPLS and explicit legacy expert invocations delegate to that exact module.
"""

import argparse
import json
import runpy
import sys
from pathlib import Path


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    probe = argparse.ArgumentParser(add_help=False)
    probe.add_argument("--forcefield")
    probe.add_argument("--profile")
    probe.add_argument("--charge-method")
    probe.add_argument("--source")
    dispatch, _ = probe.parse_known_args(argv)
    if dispatch.forcefield != "opls" or any(
        (dispatch.profile, dispatch.charge_method, dispatch.source)
    ):
        return runpy.run_path(
            str(Path(__file__).parent.parent / "inspect_forcefield.py")
        )["main"](argv)
    from island.forcefields.opls_registry import ForceFieldRequest, prepare_forcefield

    parser = argparse.ArgumentParser(
        description="Automatic canonical OPLS inspection; no MD",
        epilog="Expert-only --profile/--charge-method/--source use the historical explicit-review CLI; ordinary use requires none of them.",
    )
    graph = parser.add_mutually_exclusive_group()
    graph.add_argument("--smiles")
    graph.add_argument("--psmiles")
    parser.add_argument("--dp", type=int)
    parser.add_argument("--head-end-group")
    parser.add_argument("--tail-end-group")
    parser.add_argument("--forcefield", choices=["opls"], required=True)
    parser.add_argument("--registry", action="store_true")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.registry:
        from island.forcefields.opls_registry import canonical_registry

        args.json.write_text(
            json.dumps(canonical_registry(), sort_keys=True, indent=2) + "\n"
        )
        return 0
    if not (args.psmiles or args.smiles):
        parser.error("a final-graph SMILES or PSMILES is required")
    try:
        if args.psmiles:
            from island.builders import build_linear_polymer

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
        # This CLI inspects a declared periodic PDMS research model, not density.
        # A box is required for export; inspection uses its explicit audit default.
        from island.core import SimulationBox

        system.box = (
            SimulationBox(40, 40, 40)
            if any(a.element == "Si" for a in system.topology.sites.values())
            else None
        )
        result = prepare_forcefield(system, ForceFieldRequest("opls"))
        data = result.payload
        code = 0 if data["stages"]["parameter_complete"] else 1
    except Exception as error:  # noqa: BLE001 -- classified CLI boundary
        data = {
            "family": "opls",
            "status": "typing_ambiguous",
            "diagnostics": [{"stage": "input_construction", "reason": str(error)}],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        code = 1
    args.json.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
    return code
