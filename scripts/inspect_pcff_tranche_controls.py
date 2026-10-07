"""Accept only declared J13 chemical rejections; preserve reports on failure.

This is an acceptance CLI, not an inspection-only success code. No forces,
minimization or dynamics are executed.
"""

import argparse
from collections import Counter
from pathlib import Path

from pcff_j5_reference import charge_raw, read_source
from pcff_tranche_expectations import assess, contract_for

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldError,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY
from island.workflows import storage


def execute_control(c, d, source, source_path, output):
    result = {"case": c["name"]}
    stage = "construction"
    try:
        system = (
            build_linear_polymer(c["psmiles"], dp=c["dp"], random_seed=2026)
            if "psmiles" in c
            else from_smiles(c["smiles"], random_seed=2026)
        )
        storage.publish(
            output / (c["name"] + "-system.json"),
            storage.json_bytes(storage.encode(system.to_dict())),
        )
        stage = "typing"
        typing = type_pcff_atoms(system, source, profile=d["typing_profile"])
        result.update(
            typing_complete=typing.complete,
            labels=dict(Counter(typing.payload["assignments"].values())),
            typing_diagnostics=typing.payload["diagnostics"],
        )
        if typing.complete:
            stage = "charge_model_inspection"
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
                output / (c["name"] + "-inspection.json"),
                storage.json_bytes(storage.encode(inspection)),
            )
            result.update(
                charges_complete=inspection["charges_complete"],
                model_complete=inspection["model_complete"],
                independent_charge=independent,
                model_diagnostics=inspection["model_diagnostics"],
            )
        stage = "public_preparation"
        try:
            prepare_forcefield(
                system,
                ForceFieldRequest(
                    "pcff",
                    PCFFOptions(
                        source_path,
                        (0, 0, 1),
                        (0, 0, 1),
                        typing_profile=d["typing_profile"],
                        resolution_policy=COMPATIBILITY_POLICY,
                    ),
                ),
            )
        except (PCFFError, PreparedForceFieldError) as e:
            result["public_exception"] = {"type": type(e).__name__, "message": str(e)}
            result["public_rejection"] = type(e).__name__ + ": " + str(e)
        else:
            result["public_preparation_complete"] = True
    except Exception as e:  # noqa: BLE001 -- retain execution failures, never accept as chemical rejection
        result["execution_error"] = {
            "stage": stage,
            "type": type(e).__name__,
            "message": str(e),
        }
    return result


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
    a.output.mkdir(exist_ok=False)
    report = {
        "schema": "island_j13_negative_acceptance_v1",
        "mode": "acceptance",
        "acceptance_passed": False,
        "cases": [],
        "failures": [],
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        raw = a.declaration.read_bytes()
        d = storage.read_json(a.declaration)
        source = load_pcff_source(a.source)
        report["source"] = source.identity
        name, contract, checksum = contract_for(d, raw, source.identity)
        report.update(
            control_set=name,
            expectation_sha256=checksum,
            declaration_sha256=storage.checksum(raw),
        )
        for c in d["control_graphs"]:
            result = execute_control(c, d, source, a.source, a.output)
            report["cases"].append(result)
            result["expected"] = contract["expectations"][c["name"]]
            result["acceptance_failures"] = assess(result, result["expected"])
            result["acceptance_passed"] = not result["acceptance_failures"]
        report["acceptance_passed"] = bool(report["cases"]) and all(
            r["acceptance_passed"] for r in report["cases"]
        )
    except Exception as e:  # noqa: BLE001 -- malformed contract, source or report is a failed gate
        report["failures"].append({"type": type(e).__name__, "message": str(e)})
    storage.publish(
        a.output / "report.json", storage.json_bytes(storage.encode(report))
    )
    print(
        {
            "control_set": report.get("control_set"),
            "acceptance_passed": report["acceptance_passed"],
            "cases": [
                {
                    k: v
                    for k, v in c.items()
                    if k
                    in (
                        "case",
                        "acceptance_passed",
                        "acceptance_failures",
                        "execution_error",
                    )
                }
                for c in report["cases"]
            ],
            "failures": report["failures"],
        }
    )
    return 0 if report["acceptance_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
