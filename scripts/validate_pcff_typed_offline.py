"""Separate-process J17 relocation, offline inspection and semantic tamper checks."""

import argparse
import builtins
import json
import os
import shutil
from pathlib import Path
from unittest.mock import patch


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--workflow", type=Path, required=True)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(a.bundle, a.output / "bundle")
    shutil.copytree(a.workflow, a.output / "workflow")
    shutil.copyfile(a.source, a.output / "current.frc")
    original = builtins.__import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if level == 0 and name.split(".")[0] in {
            "openmm",
            "rdkit",
            "scipy",
            "parmed",
            "foyer",
        }:
            raise AssertionError("Scientific import during reconstruction: " + name)
        return original(name, globals, locals, fromlist, level)

    builtins.__import__ = blocked
    try:
        from island.exceptions import WorkflowError
        from island.forcefields import (
            PreparedForceFieldSources,
            load_prepared_forcefield,
            preparation,
        )
        from island.forcefields.pcff import automatic, class2, typed_graph
        from island.workflows import prepared as workflow
        from island.workflows import (
            prepared_workflow_status,
            read_prepared_workflow_frames,
            resume_prepared_workflow,
            storage,
        )

        sources = PreparedForceFieldSources(pcff_frc=a.output / "current.frc")

        def forbidden(*args, **kwargs):
            raise AssertionError(
                "Preparation/automatic perception invoked during native loading"
            )

        with (
            patch.object(preparation, "prepare_forcefield", forbidden),
            patch.object(automatic, "type_pcff_atoms", forbidden),
            patch.object(automatic, "assign_automatic_pcff_charges", forbidden),
            patch.object(class2, "assign_pcff_parameters", forbidden),
            patch.object(typed_graph, "bind_pcff_types", forbidden),
            patch.object(typed_graph, "assign_typed_pcff_charges", forbidden),
            patch.object(typed_graph, "provide_pcff_charges", forbidden),
        ):
            loaded = load_prepared_forcefield(a.output / "bundle", sources=sources)
            status = prepared_workflow_status(a.output / "workflow", sources=sources)
            frames = read_prepared_workflow_frames(
                a.output / "workflow", sources=sources
            )
            assert status == resume_prepared_workflow(
                a.output / "workflow", sources=sources
            )
            assert status["status"] == "completed" and len(frames) == 5
            assert status["prepared_identity"] == loaded.prepared.identity
        # Rechecksummed workflow identity mutation cannot attach old checkpoints.
        tamper = a.output / "tampered-workflow"
        shutil.copytree(a.output / "workflow", tamper)
        altered = dict(status)
        altered["prepared_identity"] = "0" * 64
        workflow._manifest(tamper, altered)
        rejected = []
        for operation in (
            prepared_workflow_status,
            read_prepared_workflow_frames,
            resume_prepared_workflow,
        ):
            try:
                operation(tamper, sources=sources)
            except WorkflowError as error:
                rejected.append({"operation": operation.__name__, "error": str(error)})
            else:
                raise AssertionError("Changed workflow identity accepted")
        result = {
            "passed": True,
            "pid": os.getpid(),
            "scientific_imports_blocked": True,
            "public_preparation_forbidden": True,
            "original_source_path_not_used": True,
            "source_sha256": storage.checksum((a.output / "current.frc").read_bytes()),
            "prepared": loaded.prepared.identity,
            "model": loaded.prepared.native_result.identity,
            "frames": len(frames),
            "tampering_rejected": rejected,
        }
        storage.publish(a.output / "outcome.json", storage.json_bytes(result))
        print(json.dumps(result, indent=2))
    finally:
        builtins.__import__ = original


if __name__ == "__main__":
    main()
