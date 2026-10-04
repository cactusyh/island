"""PSMILES -> local-template chain -> explicitly selected native force field."""

import argparse
import json
from pathlib import Path

from island.builders import build_linear_polymer
from island.forcefields import (
    AmberToolsOptions,
    ForceFieldRequest,
    OPLSOptions,
    PCFFOptions,
    create_evaluator,
    prepare_forcefield,
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--family", choices=("gaff", "gaff2", "oplsaa", "pcff"), required=True
    )
    p.add_argument("--psmiles", default="[*:1]CC[*:2]")
    p.add_argument("--dp", type=int, default=3)
    p.add_argument("--source", type=Path)
    p.add_argument("--charge-method", choices=("provided", "am1bcc"))
    p.add_argument(
        "--charges", type=Path, help="JSON object of exact stable ID -> charge in e"
    )
    p.add_argument("--artifacts", type=Path)
    p.add_argument("--lj", type=float, nargs=3)
    p.add_argument("--coulomb", type=float, nargs=3)
    args = p.parse_args()
    if args.family in ("gaff", "gaff2"):
        if args.charge_method is None or args.source or args.lj or args.coulomb:
            p.error("GAFF requires explicit --charge-method and no source/pair options")
        charges = (
            None
            if args.charges is None
            else {int(k): v for k, v in json.loads(args.charges.read_text()).items()}
        )
        options = AmberToolsOptions(
            args.family,
            args.charge_method,
            charges,
            work_root=args.artifacts,
            retain_success_artifacts=True,
        )
    else:
        if args.charge_method or args.charges or args.artifacts or not args.source:
            p.error(
                "Native OPLS/PCFF charges require --source and reject Amber options"
            )
        if args.family == "pcff":
            if args.lj is None or args.coulomb is None:
                p.error("PCFF requires both --lj and --coulomb")
            options = PCFFOptions(args.source, args.lj, args.coulomb)
        else:
            if args.lj or args.coulomb:
                p.error("OPLS policy comes from the pinned source")
            options = OPLSOptions(args.source)
    system = build_linear_polymer(
        args.psmiles,
        dp=args.dp,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    prepared = prepare_forcefield(system, ForceFieldRequest(args.family, options))
    evaluator = create_evaluator(system, prepared)
    with evaluator.open_session() as session:
        result = session.evaluate()
    print(
        json.dumps(
            {
                "preparation": prepared.metadata,
                "wrapper_identity": prepared.identity,
                "potential_energy_kj_mol": result.potential_energy,
                "forces_kj_mol_angstrom": dict(result.forces),
                "native_model_fingerprint": result.model_fingerprint,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
