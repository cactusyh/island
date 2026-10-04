"""Inspect a local FRC; never downloads sources or assigns unknown semantics."""

import argparse
import json
from collections import Counter

from island.exceptions import PCFFError

from .source import load_pcff_source, select


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path")
    parser.add_argument(
        "--sha256",
        help="Explicit checksum for inspection of an unreviewed local source",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Include original section lines and semantic records",
    )
    args = parser.parse_args(argv)
    try:
        source = load_pcff_source(args.path, expected_sha256=args.sha256)
        inventory = source.inventory
        diagnostics = []
        counts = []
        for section in inventory["sections"]:
            groups = {}
            for record in section["records"]:
                key = (
                    tuple(record["data"]["types"])
                    if "types" in record["data"]
                    else record["data"]["type"]
                )
                groups.setdefault(key, []).append(record)
            for candidates in groups.values():
                try:
                    select(candidates)
                except PCFFError as error:
                    diagnostics.append(str(error))
            counts.append(
                {k: section[k] for k in ("name", "namespace", "line", "interpretation")}
                | {
                    "semantic_records": len(section["records"]),
                    "raw_lines": len(section["lines"]),
                    "duplicate_groups": sum(len(g) > 1 for g in groups.values()),
                }
            )
        result = {
            "source": source.identity,
            "declarations": inventory["declarations"],
            "sections": counts,
            "diagnostics": diagnostics,
            "interpretations": dict(Counter(s["interpretation"] for s in counts)),
            "assignment_profile_available": source.identity["profile"] is not None,
        }
        if args.full:
            result["inventory"] = inventory
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except PCFFError as error:
        print(json.dumps({"status": "failed", "error": str(error)}))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
