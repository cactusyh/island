"""Inspect a separate PCFF Class-II PDMS candidate without granting execution."""

import argparse
import json
from pathlib import Path

from island.forcefield_policy import ClaimRequirements
from island.pcff_pdms import PCFFPDMSAssets, prepare_polymer, prepare_polymer_forcefield


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--smiles")
    g.add_argument("--psmiles")
    p.add_argument("--forcefield", choices=("pcff",), default="pcff")
    p.add_argument(
        "--distribution",
        type=Path,
        help="administrator-installed, pinned original sources",
    )
    p.add_argument("--dp", type=int)
    p.add_argument("--head-end-group")
    p.add_argument("--tail-end-group")
    p.add_argument("--require-parameter-complete", action="store_true")
    p.add_argument("--require-cross-engine-equivalence", action="store_true")
    p.add_argument("--json", type=Path, required=True)
    a = p.parse_args(argv)
    claims = ClaimRequirements(
        parameter_complete=a.require_parameter_complete,
        cross_engine_equivalence=a.require_cross_engine_equivalence,
    )
    try:
        assets = PCFFPDMSAssets(a.distribution) if a.distribution else None
        if a.psmiles:
            if a.dp is None or not a.head_end_group or not a.tail_end_group:
                from island.pcff_pdms import Class2AuditError

                raise Class2AuditError(
                    "typing_ambiguous",
                    "repeat_is_not_a_final_graph; explicit_DP_and_terminal_groups_required",
                )
            _, prepared = prepare_polymer(
                psmiles=a.psmiles,
                forcefield="pcff",
                dp=a.dp,
                head_end_group=a.head_end_group,
                tail_end_group=a.tail_end_group,
                assets=assets,
            )
        else:
            from island.chemistry import from_smiles

            system = from_smiles(a.smiles, random_seed=2026)
            prepared = prepare_polymer_forcefield(
                system, forcefield="pcff", assets=assets
            )
        report = prepared.payload
        report["requested_claim_failures"] = claims.failures(report)
        code = int(bool(report["requested_claim_failures"]))
    except (ValueError, RuntimeError) as error:
        report = {
            "family": "pcff",
            "parameter_completeness": "incomplete",
            "selected_profile": None,
            "status": getattr(error, "status", "source_incompatible"),
            "diagnostics": getattr(error, "diagnostics", [str(error)]),
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        code = 1
    a.json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
