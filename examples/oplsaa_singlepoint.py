"""Pinned OPLS-AA source -> local-template chain -> owned snapshot -> single point."""

import argparse
import json

from island.builders import build_linear_polymer
from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.forcefields.oplsaa import load_oplsaa_source, parameterize_oplsaa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xml", required=True, help="Installed checksum-pinned oplsaa.xml"
    )
    parser.add_argument("--psmiles", default="[*:1]CC[*:2]")
    parser.add_argument("--dp", type=int, default=3)
    args = parser.parse_args()
    system = build_linear_polymer(
        args.psmiles,
        dp=args.dp,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    source = load_oplsaa_source(args.xml)
    parameters = parameterize_oplsaa(system, source)
    snapshot = parameters.to_parameterized_system(system, source)
    result = OPLSSinglePointEvaluator.from_parameterized_system(
        snapshot, source
    ).evaluate()
    print(
        json.dumps(
            {
                "source": source.identity,
                "parameters": parameters.identity,
                "coverage": parameters.payload["resolved"]["coverage"],
                "energy_kj_mol": result.potential_energy,
                "components_kj_mol": dict(result.energy_components),
                "forces_kj_mol_angstrom": dict(result.forces),
                "coordinate_fingerprint": result.coordinate_fingerprint,
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
