"""Save a loaded preparation, relocate it, then evaluate in another process.

Start from a validated I2 bundle. Source libraries are resolved explicitly.
No parameterization, typing, coordinate construction or QM runs here.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from island.forcefields import (
    AmberBundleArtifacts,
    PreparedForceFieldSources,
    create_evaluator,
    load_prepared_forcefield,
    save_prepared_forcefield,
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--output", type=Path)
    p.add_argument("--xml", type=Path)
    p.add_argument("--frc", type=Path)
    p.add_argument("--child", action="store_true")
    args = p.parse_args()
    sources = PreparedForceFieldSources(opls_xml=args.xml, pcff_frc=args.frc)
    loaded = load_prepared_forcefield(args.bundle, sources=sources)
    system, prepared = loaded.system, loaded.prepared
    if args.child:
        result = create_evaluator(system, prepared).evaluate_fresh()
        print(
            json.dumps(
                {
                    "pid": os.getpid(),
                    "family": prepared.metadata["family"],
                    "prepared_identity": prepared.identity,
                    "energy_kj_mol": result.potential_energy,
                    "components_kj_mol": dict(result.energy_components),
                    "model_fingerprint": result.model_fingerprint,
                    "production_validated": False,
                    "simulation_readiness": "not_established",
                },
                indent=2,
            )
        )
        return
    if args.output is None:
        p.error("--output required for save/relocate")
    args.output.mkdir(parents=True, exist_ok=False)
    artifacts = None
    if prepared.metadata["family"] in ("gaff", "gaff2"):
        # Explicitly copy the validated original prmtop, not its historical path.
        artifacts = AmberBundleArtifacts(args.bundle / "result.prmtop")
    save_prepared_forcefield(
        system, prepared, args.output / "before", sources=sources, artifacts=artifacts
    )
    relocated = args.output / "relocated"
    shutil.move(args.output / "before", relocated)
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--bundle",
        str(relocated.resolve()),
        "--child",
    ]
    for key in ("xml", "frc"):
        value = getattr(args, key)
        if value is not None:
            command.extend(["--" + key, str(value.resolve())])
    print(
        "Parent PID:",
        os.getpid(),
        "; child reconstructs the relocated bundle.",
        flush=True,
    )
    subprocess.run(command, cwd=args.output, check=True)


if __name__ == "__main__":
    main()
