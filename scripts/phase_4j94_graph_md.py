"""Reproduce J94 graph assignment and short NVT runs on actual block and branch graphs."""

import argparse
import json
from pathlib import Path

from island.builders import build_block_copolymer
from island.forcefields.assignment import assign_polymer_parameters, run_emc_md
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_core.schema import Branch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--lammps", default="lmp")
    parser.add_argument("--steps", type=int, default=20_000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    ether = build_block_copolymer([("*CC*", 3), ("*CCO*", 2), ("*CC*", 1)])
    styrene = build_block_copolymer([("*CC*", 3), ("*CC(c1ccccc1)*", 2), ("*CC*", 1)])
    branch = normalize(PolymerSpec((("A", "*CC*"),), 4, branches=(Branch(4, "*CC"),)))
    branched_block = normalize(
        PolymerSpec(
            (("A", "*CC*"), ("B", "*CC(c1ccccc1)*")),
            6,
            sequence=("A", "A", "A", "B", "B", "A"),
            branches=(Branch(2, "*CC"), Branch(4, "*CCC")),
        )
    )
    long_block = build_block_copolymer(
        [("*CC*", 30), ("*CC(c1ccccc1)*", 10)], generate_3d=False
    )
    cases = [
        ("ether_block", ether, ("pcff", "opls/2012/opls-aa", "opls/2024/opls-ua")),
        ("styrene_block", styrene, ("pcff", "opls/2024/opls-aa")),
        ("ethyl_branch", branch, ("pcff", "opls/2024/opls-aa", "opls/2024/opls-ua")),
        ("branched_styrene_block", branched_block, ("pcff", "opls/2024/opls-aa")),
        ("dp40_styrene_block", long_block, ("pcff", "opls/2024/opls-aa")),
    ]
    summary = []
    for label, graph, fields in cases:
        for field in fields:
            assignment = assign_polymer_parameters(graph, field)
            name = label + "_" + field.replace("/", "_")
            result = run_emc_md(
                assignment,
                args.output / name,
                steps=args.steps,
                sample_every=min(1000, args.steps),
                lammps=args.lammps,
            )
            summary.append(
                {
                    "run": name,
                    "field": assignment.field,
                    "input_atom_count": sum(map(len, assignment.site_groups.values())),
                    "lammps_particle_count": len(assignment.site_groups),
                    "coordinate_source": assignment.coordinate_source,
                    "steps": result.steps,
                    "time_fs": result.time_fs,
                    "frames": result.frames,
                    "final_temperature_K": result.final_temperature_k,
                    "final_potential_kcal_mol": result.final_potential_kcal_mol,
                }
            )
            print(name, result.steps, result.final_temperature_k, flush=True)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
