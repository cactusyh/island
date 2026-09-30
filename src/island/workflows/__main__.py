"""Thin data-only CLI: python -m island.workflows COMMAND ..."""

import argparse
import json

from island.exceptions import IslandError

from . import (
    WorkflowConfig,
    inspect_chain,
    resume_workflow,
    start_workflow,
    workflow_status,
)
from .storage import read_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("inspect", "start"):
        p = sub.add_parser(name)
        p.add_argument("config")
        if name == "start":
            p.add_argument("--segments", type=int, default=1)
    for name in ("resume", "status"):
        p = sub.add_parser(name)
        p.add_argument("directory")
        if name == "resume":
            p.add_argument("--segments", type=int, default=1)
    args = parser.parse_args()
    try:
        if args.command in ("inspect", "start"):
            config = WorkflowConfig.from_dict(read_json(args.config))
            result = (
                inspect_chain(config).to_dict()
                if args.command == "inspect"
                else start_workflow(config, segments=args.segments)
            )
        else:
            result = (
                workflow_status(args.directory)
                if args.command == "status"
                else resume_workflow(args.directory, segments=args.segments)
            )
        print(json.dumps(result, indent=2, allow_nan=False))
        return 1 if result.get("status") == "stage_failed" else 0
    except IslandError as error:
        parser.exit(2, f"{type(error).__name__}: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
