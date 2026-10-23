"""Fresh-process NPT continuation that retains genuine partial failure histories."""

import argparse
import base64
import traceback
from pathlib import Path

from island.analysis.npt import FLAGS, NPTConfiguration, _signed, _write_json
from island.analysis.npt_continuation import capture, load_continuation, require
from island.periodic_md import _engine


def resume(source, destination, *, expected_restart, expected_parent, expected_config):
    import openmm as mm
    from openmm import unit

    workflow, restart = load_continuation(
        source,
        expected_restart=expected_restart,
        expected_parent=expected_parent,
        expected_config=expected_config,
    )
    require(not destination.exists(), "Resume destination exists")
    manifest, frames = workflow["manifest"], workflow["frames"]
    cfg = NPTConfiguration.from_payload(manifest["config"])
    boundary = restart["split_boundary_npt_step"]
    require(
        boundary is not None and cfg.production_steps == 0, "No resumable NPT boundary"
    )
    require(_engine(mm) == manifest["engine"], "Changed execution engine")
    model = mm.XmlSerializer.deserialize((source / "system.xml").read_text())
    integrator = mm.XmlSerializer.deserialize((source / "integrator.xml").read_text())
    ctx = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    status = "completed"
    failure = None
    attempted_npt_step = None
    try:
        ctx.loadCheckpoint(
            base64.b64decode(
                (source / "npt-initial-checkpoint.txt").read_text(), validate=True
            )
        )
        for _ in range(boundary):
            integrator.step(1)
        require(
            capture(ctx, cfg, len(manifest["atom_ids"]), frames[-1]["step"], mm, unit)
            == frames[-1],
            "Independent barostat prefix/RNG replay differs",
        )
        ctx.loadCheckpoint(
            base64.b64decode(frames[-1]["checkpoint_base64"], validate=True)
        )
        require(
            capture(ctx, cfg, len(manifest["atom_ids"]), frames[-1]["step"], mm, unit)
            == frames[-1],
            "Restored checkpoint/RNG differs",
        )
        for i in range(boundary + 1, cfg.npt_steps + 1):
            attempted_npt_step = i
            integrator.step(1)
            if i % cfg.checkpoint_interval == 0 or i == cfg.npt_steps:
                frames.append(
                    capture(
                        ctx, cfg, len(manifest["atom_ids"]), cfg.nvt_steps + i, mm, unit
                    )
                )
    except Exception as error:  # noqa: BLE001 -- preserve failed real continuation
        status = "failed"
        failure = {
            "exception": type(error).__name__,
            "diagnostic": str(error),
            "attempted_npt_step": attempted_npt_step,
        }
        destination.mkdir(parents=True)
        (destination / "failure.log").write_text(traceback.format_exc())
    finally:
        del ctx, integrator
    destination.mkdir(parents=True, exist_ok=True)
    record = _signed(
        dict(
            schema="island.j40-c2-continuation-attempt.v1",
            parent_restart=restart,
            config=cfg.payload,
            backend_identity=manifest["backend_identity"],
            heating=manifest["heating"],
            frames=frames,
            status=status,
            failure=failure,
            **FLAGS,
        )
    )
    _write_json(destination / "attempt.json", record)
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--resume", type=Path, required=True)
    p.add_argument("--destination", type=Path, required=True)
    p.add_argument("--restart", required=True)
    p.add_argument("--parent", required=True)
    p.add_argument("--config", required=True)
    a = p.parse_args()
    r = resume(
        a.resume,
        a.destination,
        expected_restart=a.restart,
        expected_parent=a.parent,
        expected_config=a.config,
    )
    return int(r["status"] != "completed")


if __name__ == "__main__":
    raise SystemExit(main())
