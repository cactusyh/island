"""Offline CLI: inventory, inspect and report; never execute a supplied package."""

import argparse
import json
from pathlib import Path

from island.opls_source_intake import (
    SourcePackage,
    blocker_reopen_report,
    inspect_source,
    inventory,
)
from island.opls_source_intake.schema import source_intake_schema


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("schema", "inventory", "inspect"))
    parser.add_argument("--package", type=Path)
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--json", type=Path)
    parser.add_argument("--blocker-reopen-json", type=Path)
    args = parser.parse_args(argv)
    if args.command == "schema":
        value, code = source_intake_schema(), 0
    elif args.command == "inventory":
        if args.package is None:
            parser.error("inventory requires --package")
        value, code = inventory(args.package), 0
    else:
        package = (
            SourcePackage(args.package, args.manifest_sha256) if args.package else None
        )
        value = inspect_source(package)
        code = 0 if value["status"] == "source_complete_not_executed" else 1
        if args.blocker_reopen_json:
            args.blocker_reopen_json.write_text(
                json.dumps(blocker_reopen_report(value), indent=2, sort_keys=True)
                + "\n"
            )
    raw = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(raw)
    else:
        print(raw, end="")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
