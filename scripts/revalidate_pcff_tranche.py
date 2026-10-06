"""Read-only bundle/workflow reconstruction with scientific dependencies blocked."""

import argparse
import builtins
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--root", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    original = builtins.__import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if level == 0 and name.split(".")[0] in {
            "openmm",
            "rdkit",
            "scipy",
            "foyer",
            "parmed",
        }:
            raise AssertionError(
                "Scientific import during offline reconstruction: " + name
            )
        return original(name, globals, locals, fromlist, level)

    builtins.__import__ = blocked
    try:
        from island.forcefields import (
            PreparedForceFieldSources,
            load_prepared_forcefield,
        )
        from island.workflows import (
            prepared_workflow_status,
            read_prepared_workflow_frames,
            resume_prepared_workflow,
            storage,
        )

        results = []
        for root in args.root:
            before = {
                str(p): storage.checksum(p.read_bytes())
                for p in root.rglob("*")
                if p.is_file()
            }
            sources = PreparedForceFieldSources(pcff_frc=args.source)
            bundle = load_prepared_forcefield(
                root / "bundle-relocated", sources=sources
            )
            run = root / "workflow-relocated"
            status = prepared_workflow_status(run, sources=sources)
            frames = read_prepared_workflow_frames(run, sources=sources)
            resumed = resume_prepared_workflow(run, sources=sources)
            if status != resumed or status["status"] != "completed":
                raise ValueError("Completed no-op/status mismatch")
            if any(
                storage.checksum(Path(p).read_bytes()) != h for p, h in before.items()
            ):
                raise ValueError("Read-only artifacts changed")
            results.append(
                {
                    "root": str(root),
                    "status": status["status"],
                    "frames": len(frames),
                    "native": bundle.prepared.native_result.identity,
                    "facade": bundle.prepared.identity,
                    "historical_hashes_unchanged": True,
                    "scientific_imports_blocked": True,
                }
            )
        storage.publish(args.output, storage.json_bytes(results))
    finally:
        builtins.__import__ = original


if __name__ == "__main__":
    main()
