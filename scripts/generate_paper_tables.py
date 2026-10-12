"""Generate machine-readable tables from immutable analysis records.

No scientific values are fabricated: an empty or blocked input produces a
table with an explicit status and no labels.
"""

import argparse
import csv
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
        for name, value in record.payload["observables"].items():
            rows.append({"record": str(path), "record_identity": record.identity, "observable": name, "status": value.get("status"), "value": value.get("value"), "uncertainty": value.get("uncertainty"), "units": value.get("units")})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["record", "record_identity", "observable", "status", "value", "uncertainty", "units"])
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"output": str(args.output), "rows": len(rows), "status": "no_records" if not rows else "generated"}, sort_keys=True))


if __name__ == "__main__":
    main()
