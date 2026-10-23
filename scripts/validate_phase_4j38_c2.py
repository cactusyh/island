"""Source incompatibility and exact-graph asset inventory; unmet acceptance exits 1."""

import argparse
import json
from pathlib import Path

from scripts.phase_4j38_c2 import (
    SourceIncompatibleError,
    aromatic_decision,
    matrix,
    require_native_source,
)
from scripts.validate_phase_4j38_c1 import run as c1_run
from scripts.validate_phase_4j38_c1 import write


def run(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    try:
        c1 = c1_run(output / "c1-recheck")
        read = lambda name: json.loads((output / "c1-recheck" / name).read_text())
        decision = aromatic_decision(
            read("aromatic-audit.json"), read("aromatic-fixture.json")
        )
        write(output / "native-source-decision.json", decision)
        try:
            require_native_source(decision)
        except SourceIncompatibleError as exc:
            write(
                output / "native-rejection.json",
                {
                    "exception": type(exc).__name__,
                    "diagnostic": str(exc),
                    "classification": "source_incompatible",
                    "dynamics_executed": False,
                },
            )
        inventory = read("resp-inventory.json")
        result = matrix(c1, decision, inventory)
        write(output / "matrix.json", result)
        write(
            output / "run-state.json",
            {
                "status": "diagnostic_complete",
                "acceptance_met": False,
                "dynamics_executed": False,
            },
        )
        return result
    except Exception as exc:
        write(
            output / "run-state.json",
            {
                "status": "failed",
                "acceptance_met": False,
                "diagnostic": str(exc),
                "exception": type(exc).__name__,
                "dynamics_executed": False,
            },
        )
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.output)
    print(
        "Diagnostic complete; acceptance unmet: native aromatic source incompatible; PE/aromatic exact RESP assets unavailable"
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
