"""Read-only bundle/workflow inspection with optional science imports blocked."""

import argparse
import builtins
from pathlib import Path

from island.forcefields import PreparedForceFieldSources, load_prepared_forcefield
from island.workflows import (
    prepared_workflow_status,
    read_prepared_workflow_frames,
    resume_prepared_workflow,
)
from island.workflows.storage import checksum, json_bytes, publish


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--bundle", type=Path, action="append", default=[])
    p.add_argument("--workflow", type=Path, action="append", default=[])
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    sources = PreparedForceFieldSources(pcff_frc=a.source)
    importer = builtins.__import__

    def blocked(name, *args, **kw):
        if name.split(".")[0] in ("openmm", "scipy", "rdkit", "parmed", "foyer"):
            raise AssertionError("Forbidden scientific import: " + name)
        return importer(name, *args, **kw)

    def hashes(path):
        return {
            str(f.relative_to(path)): checksum(f.read_bytes())
            for f in path.rglob("*")
            if f.is_file()
        }

    result = {"bundles": [], "workflows": [], "passed": False}
    builtins.__import__ = blocked
    try:
        for path in a.bundle:
            before = hashes(path)
            loaded = load_prepared_forcefield(path, sources=sources)
            assert hashes(path) == before
            result["bundles"].append(
                {
                    "path": str(path),
                    "identity": loaded.prepared.identity,
                    "native": loaded.prepared.native_result.identity,
                    "hashes": before,
                }
            )
        for path in a.workflow:
            before = hashes(path)
            status = prepared_workflow_status(path, sources=sources)
            frames = read_prepared_workflow_frames(path, sources=sources)
            assert status["status"] == "completed"
            assert resume_prepared_workflow(path, sources=sources) == status
            assert hashes(path) == before
            result["workflows"].append(
                {
                    "path": str(path),
                    "status": status["status"],
                    "frames": len(frames),
                    "identity": status["prepared_identity"],
                    "hashes": before,
                }
            )
        result["passed"] = True
    finally:
        builtins.__import__ = importer
    publish(a.output, json_bytes(result))
    print("Read-only checks passed")


if __name__ == "__main__":
    main()
