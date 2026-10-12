"""Generate supplementary metadata tables from analysis records only."""

import argparse
import json
from pathlib import Path

from island.analysis import load_analysis_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", nargs="*", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for path in args.records:
        record = load_analysis_record(path)
        metadata = record.payload["trajectory"]["metadata"]
        rows.append({"record": str(path), "record_identity": record.identity, "frame_count": record.payload["trajectory"]["frame_count"], "final_graph_identity": metadata["final_graph_identity"], "packing_identity": metadata["packing_identity"], "force_field_family": metadata["force_field_family"], "workflow_identity": metadata["workflow_identity"], "checkpoint_identity": metadata["checkpoint_identity"]})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"schema": "island.paper-supplementary-table.v1", "rows": rows, "status": "no_records" if not rows else "generated"}, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "rows": len(rows)}, sort_keys=True))


if __name__ == "__main__":
    main()
