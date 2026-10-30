"""Inspect separate parameter, backend and scientific claims (J55).

Asset installation is independent of the ordinary forcefield="opls" choice.
No diagnostics or runtime is executed unless explicitly requested.
"""

import argparse
import json
from pathlib import Path

from island.forcefield_policy import (
    ClaimRequirements,
    OPLSAssets,
    check_lammps,
    check_openmm,
    export_lammps,
    prepare_polymer_forcefield,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    graph = parser.add_mutually_exclusive_group(required=True)
    graph.add_argument("--smiles")
    graph.add_argument("--psmiles")
    parser.add_argument(
        "--forcefield", choices=("opls", "pcff", "gaff", "gaff2"), required=True
    )
    parser.add_argument("--dp", type=int)
    parser.add_argument("--head-end-group")
    parser.add_argument("--tail-end-group")
    parser.add_argument("--box", type=float, nargs=3)
    parser.add_argument("--assets", type=Path, help="administrator asset configuration")
    parser.add_argument("--assets-sha256")
    parser.add_argument("--paper", type=Path)
    parser.add_argument("--supplement", type=Path)
    parser.add_argument("--export-lammps", type=Path)
    parser.add_argument("--lammps", type=Path)
    parser.add_argument("--lammps-sha256")
    parser.add_argument("--runtime-output", type=Path)
    parser.add_argument("--diagnostic", type=Path)
    parser.add_argument("--diagnostic-sha256")
    for name in (
        "parameter-complete",
        "lammps-runtime",
        "openmm-runtime",
        "cross-engine-equivalence",
        "scientific-validation",
    ):
        parser.add_argument("--require-" + name, action="store_true")
    parser.add_argument("--build-openmm", action="store_true")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args(argv)
    wanted = ClaimRequirements(
        args.require_parameter_complete,
        args.require_lammps_runtime,
        args.require_openmm_runtime,
        args.require_cross_engine_equivalence,
        args.require_scientific_validation,
    )
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
        if args.box:
            from island.core import SimulationBox

            system.box = SimulationBox(*args.box)
        assets = None
        if args.assets:
            from hashlib import sha256

            raw = args.assets.read_bytes()
            if sha256(raw).hexdigest() != args.assets_sha256:
                raise ValueError("independent asset configuration hash required")
            assets = OPLSAssets(**json.loads(raw))
        sources = None
        if args.paper and args.supplement:
            from island.opls_resolution import ConventionSources

            sources = ConventionSources(args.paper, args.supplement)
        prepared = prepare_polymer_forcefield(
            system, args.forcefield, assets=assets, sources=sources
        )
        if args.export_lammps:
            prepared = export_lammps(system, prepared, args.export_lammps)
        if args.require_openmm_runtime:
            prepared = check_openmm(system, prepared)
        elif args.build_openmm:
            from island.forcefield_policy import build_openmm

            prepared, evaluator = build_openmm(system, prepared)
            if hasattr(evaluator, "close"):
                evaluator.close()
        if args.require_lammps_runtime:
            from island.pdms_lammps import LAMMPSOptions

            if not (args.lammps and args.lammps_sha256 and args.runtime_output):
                from island.forcefield_policy import _update

                prepared = _update(
                    prepared,
                    system,
                    lammps_runtime_status="blocked",
                    lammps_runtime_diagnostic="LAMMPS runtime requires pinned executable and output directory",
                )
            else:
                prepared = check_lammps(
                    system,
                    prepared,
                    args.runtime_output,
                    engine=LAMMPSOptions(args.lammps, args.lammps_sha256),
                )
        if args.diagnostic:
            from island.forcefield_policy.evidence import record_static_comparison

            prepared = record_static_comparison(
                system,
                prepared,
                args.diagnostic,
                expected_sha256=args.diagnostic_sha256,
            )
        report = prepared.payload
        report["requested_claim_failures"] = wanted.failures(report)
        report["preparation_identity"] = report.pop("identity")
        from island.forcefields.selection.inventory import digest

        report["identity"] = digest(report)
        code = int(
            bool(report["requested_claim_failures"])
            or report["source_integrity"] != "passed"
        )
    except Exception as error:  # noqa: BLE001 -- classified public inspection boundary
        report = {
            "status": "blocked",
            "diagnostics": str(error),
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
        code = 1
    args.json.write_text(
        json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
