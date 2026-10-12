"""Generate provenance-labelled figures from analysis-record tables.

Matplotlib is optional.  Without it the command emits a JSON manifest that
keeps the figure gate blocked instead of substituting a hand-edited image.
"""

import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = list(csv.DictReader(args.table.open()))
    manifest = {"schema": "island.paper-figure.v1", "table": str(args.table), "rows": len(rows), "status": "blocked_no_records" if not rows else "generated", "scientific_readiness": "not_established"}
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        manifest["status"] = "blocked_missing_matplotlib"
        args.output.with_suffix(".json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        return
    numeric = [float(row["value"]) for row in rows if row.get("status") == "ok" and row.get("value") not in (None, "") and row["value"].replace(".", "", 1).replace("-", "", 1).isdigit()]
    if not numeric:
        manifest["status"] = "blocked_no_numeric_observables"
        args.output.with_suffix(".json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        return
    fig, axis = plt.subplots(figsize=(5, 3.2))
    axis.plot(range(len(numeric)), numeric, marker="o")
    axis.set_xlabel("recorded observation index")
    axis.set_ylabel("value (declared units)")
    axis.set_title("ISLAND analysis records")
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
    manifest["figure_sha256"] = __import__("hashlib").sha256(args.output.read_bytes()).hexdigest()
    args.output.with_suffix(".json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
