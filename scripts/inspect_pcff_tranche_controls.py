"""Declared J13 protonation and near-miss controls; no model execution on failure."""

import argparse
from collections import Counter
from pathlib import Path

from pcff_j5_reference import charge_raw, read_source

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields import ForceFieldRequest, PCFFOptions, prepare_forcefield
from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY
from island.workflows import storage


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument(
        "--declaration",
        type=Path,
        default=Path("docs/evidence/phase_4j13_declaration.json"),
    )
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    d = storage.read_json(a.declaration)
    source = load_pcff_source(a.source)
    if source.identity["sha256"] != d["source"]["sha256"]:
        raise ValueError("Source mismatch")
    a.output.mkdir(exist_ok=False)
    cases = []
    for c in d["control_graphs"]:
        system = (
            build_linear_polymer(c["psmiles"], dp=c["dp"], random_seed=2026)
            if "psmiles" in c
            else from_smiles(c["smiles"], random_seed=2026)
        )
        storage.publish(
            a.output / (c["name"] + "-system.json"),
            storage.json_bytes(storage.encode(system.to_dict())),
        )
        typing = type_pcff_atoms(system, source, profile=d["typing_profile"])
        result = {
            "case": c["name"],
            "typing_complete": typing.complete,
            "labels": dict(Counter(typing.payload["assignments"].values())),
            "typing_diagnostics": typing.payload["diagnostics"],
        }
        if typing.complete:
            inspection = inspect_pcff_operational_support(
                system,
                typing,
                resolution_policy=COMPATIBILITY_POLICY,
                special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
            )
            independent = charge_raw(
                read_source(source.raw), system, typing.assignments
            )
            if independent["complete"] != inspection["charges_complete"]:
                raise ValueError("Independent charge status mismatch")
            storage.publish(
                a.output / (c["name"] + "-inspection.json"),
                storage.json_bytes(storage.encode(inspection)),
            )
            result.update(
                charges_complete=inspection["charges_complete"],
                model_complete=inspection["model_complete"],
                independent_charge=independent,
                model_diagnostics=inspection["model_diagnostics"],
            )
        try:
            prepare_forcefield(
                system,
                ForceFieldRequest(
                    "pcff",
                    PCFFOptions(
                        a.source,
                        (0, 0, 1),
                        (0, 0, 1),
                        typing_profile=d["typing_profile"],
                        resolution_policy=COMPATIBILITY_POLICY,
                    ),
                ),
            )
        except Exception as e:  # noqa: BLE001 -- retain failed control evidence
            result["public_rejection"] = type(e).__name__ + ": " + str(e)
        else:
            result["public_preparation_complete"] = True
        cases.append(result)
    storage.publish(
        a.output / "report.json",
        storage.json_bytes(
            storage.encode(
                {
                    "source": source.identity,
                    "cases": cases,
                    "full_source_complete": False,
                    "production_validated": False,
                    "simulation_readiness": "not_established",
                }
            )
        ),
    )
    print(
        [
            {
                k: v
                for k, v in c.items()
                if k
                in (
                    "case",
                    "typing_complete",
                    "charges_complete",
                    "model_complete",
                    "public_rejection",
                )
            }
            for c in cases
        ]
    )


if __name__ == "__main__":
    main()
