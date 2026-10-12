"""Validate a J23 workflow, calculate declared observables, and write a record.

The command intentionally writes a blocked/negative result when a workflow is
not physically sufficient.  It never creates labels or regenerates chemistry.
"""

import argparse
import json
from pathlib import Path

from island.analysis import (
    calculate_observables,
    load_verified_trajectory,
    make_analysis_record,
    observable_definitions,
    write_analysis_record,
)
from island.exceptions import AnalysisError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", type=Path)
    parser.add_argument("--observables", nargs="+", required=True, choices=sorted(observable_definitions()))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-size", type=int, default=5)
    parser.add_argument("--expected-identity")
    args = parser.parse_args()
    try:
        trajectory = load_verified_trajectory(args.workflow, expected_identity=args.expected_identity)
        values = calculate_observables(trajectory, tuple(args.observables), block_size=args.block_size)
        record = make_analysis_record(
            trajectory,
            values,
            provenance={"command": "scripts/analyze_periodic_trajectory.py", "block_size": args.block_size},
            evidence_references=("docs/phase_4j23.md", "docs/evidence/phase_4j23.json"),
        )
        path = write_analysis_record(record, args.output)
        print(json.dumps({"path": str(path), "identity": record.identity, "observables": values}, sort_keys=True))
    except AnalysisError as error:
        parser.exit(2, f"analysis: {error}\n")


if __name__ == "__main__":
    main()
