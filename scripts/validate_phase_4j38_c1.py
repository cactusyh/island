"""Reproduce failed J38 cells, publish diagnostics, return nonzero acceptance."""

import argparse
import json
from pathlib import Path

from scripts.phase_4j38_c1_audit import aromatic, require, resp_inventory
from scripts.phase_4j38_c1_contract import FAILED, matrix


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def run(output):
    from island.exceptions import OPLSAssignmentError
    from island.forcefields.oplsaa import load_oplsaa_source, type_atoms
    from island.periodic import periodic_chemical_view
    from scripts.phase_4j32_assignment import prepare
    from scripts.phase_4j32_real_fixture import XML
    from scripts.phase_4j38_fixture import build

    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    write(
        output / "run-state.json",
        {"status": "started", "acceptance_met": False, "dynamics_executed": False},
    )
    try:
        fixtures = {}
        for key in ("pe", "aromatic"):
            fixture, evidence = build(key)
            fixtures[key] = fixture
            write(output / f"{key}-fixture.json", evidence)
        evidence = json.loads((output / "aromatic-fixture.json").read_text())
        report = aromatic(evidence, XML)
        system = periodic_chemical_view(fixtures["aromatic"].packed.system)
        typing = type_atoms(system, load_oplsaa_source(XML))
        assignments = typing.payload["data"]["assignments"]
        require(
            assignments == {r["atom_id"]: r["source_type"] for r in report["rows"]},
            "Native/external typing disagreement",
        )
        write(output / "aromatic-audit.json", report)
        try:
            prepare(fixtures["aromatic"], "oplsaa", output / "native-preparation")
        except OPLSAssignmentError as exc:
            require(
                str(exc) == "Native component neutrality failed: -0.23 e",
                "Unexpected native rejection",
            )
            rejection = {
                "exception": type(exc).__name__,
                "diagnostic": str(exc),
                "stage": "charge_completeness",
                "classification": "source_incomplete",
                "final_graph_identity": fixtures["aromatic"].packed.graph.identity,
                "dynamics_executed": False,
                "typing_identity": typing.identity,
            }
        else:
            raise ValueError(
                "Unexpected native aromatic success; requires separate source review"
            )
        write(output / "native-rejection.json", rejection)
        inventory = resp_inventory(
            {f"gaff2/{key}": f.packed.graph.identity for key, f in fixtures.items()}
        )
        write(output / "resp-inventory.json", inventory)
        history = json.loads(Path("docs/evidence/phase_4j38/matrix.json").read_text())
        result = matrix(history, report, inventory)
        write(output / "matrix.json", result)
        write(
            output / "run-state.json",
            {
                "status": "diagnostic_complete",
                "acceptance_met": False,
                "dynamics_executed": False,
                "unmet_cases": sorted(FAILED),
            },
        )
        return result
    except Exception as exc:
        write(
            output / "run-state.json",
            {
                "status": "failed",
                "acceptance_met": False,
                "dynamics_executed": False,
                "exception": type(exc).__name__,
                "diagnostic": str(exc),
            },
        )
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.output)
    print(
        "Diagnostic complete; acceptance unmet: OPLS source incomplete and exact-graph RESP assets unavailable"
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
