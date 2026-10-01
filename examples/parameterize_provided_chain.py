"""Inspect exact chain IDs, then prepare using an explicit JSON charge file.

Charges JSON: {"source": "citation/file/method description", "charges": {"1": -0.1, ...}}
No values are filled, copied across repeats, renormalized or inferred.
"""

import argparse
import os
from pathlib import Path

from island.forcefields import AmberToolsParameterizationEngine
from island.workflows import WorkflowConfig, inspect_chain, storage
from island.workflows.bundle import system_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("inspect", "prepare"))
    parser.add_argument("--psmiles", default="[*]CC[*]")
    parser.add_argument("--dp", type=int, default=3)
    parser.add_argument("--max-atoms", type=int, default=100)
    parser.add_argument("--force-field", choices=("gaff", "gaff2"), default="gaff2")
    parser.add_argument("--charges", type=Path)
    parser.add_argument("--amberhome", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    charges, source = {}, None
    if args.operation == "prepare":
        if args.charges is None:
            parser.error(
                "prepare requires --charges; no automatic charge model is used"
            )
        raw = storage.read_json(args.charges)
        if type(raw) is not dict or set(raw) != {"source", "charges"}:
            parser.error("Charge file requires exactly source and charges")
        source = raw["source"]
        values = raw["charges"]
        if type(values) is not dict or any(str(int(k)) != k for k in values):
            parser.error("Charge keys must be canonical integer stable IDs")
        charges = {int(k): q for k, q in values.items()}
    cfg = WorkflowConfig(
        args.psmiles,
        args.dp,
        str(args.output),
        args.force_field,
        "provided",
        provided_charges=charges,
        max_atoms=args.max_atoms,
        charge_source=source,
        amberhome=str(args.amberhome) if args.amberhome else None,
    )
    system = inspect_chain(cfg)
    if args.operation == "inspect":
        storage.publish(
            args.output,
            storage.json_bytes(
                {
                    "actual_atoms": system.number_of_sites,
                    "seeds": {
                        "template": cfg.template_seed,
                        "assembly": cfg.assembly_seed,
                    },
                    "sites": [
                        {
                            "id": s,
                            "element": a.element,
                            "formal_charge": a.formal_charge,
                        }
                        for s, a in sorted(system.topology.sites.items())
                    ],
                }
            ),
        )
        print(f"Recorded {system.number_of_sites} sites; no charges generated")
        return
    if args.amberhome:
        os.environ["PATH"] = (
            str(args.amberhome / "bin") + os.pathsep + os.environ.get("PATH", "")
        )
    args.output.mkdir(parents=True, exist_ok=False)
    prepared = AmberToolsParameterizationEngine().parameterize(
        system, cfg.amber_options(args.output)
    )
    storage.publish(
        args.output / "input-system.json",
        storage.json_bytes(storage.encode(system_data(system))),
    )
    storage.publish(
        args.output / "preparation.json",
        storage.json_bytes(
            storage.encode(
                {
                    "record": dict(prepared.record),
                    "record_signature": prepared.record_signature,
                    "import_source": prepared.imported_result.source,
                    "import_provenance": dict(prepared.imported_result.provenance),
                }
            )
        ),
    )
    print(prepared.record_signature)
    print(
        "Provided-charge preparation completed; scientific suitability remains unestablished"
    )


if __name__ == "__main__":
    main()
