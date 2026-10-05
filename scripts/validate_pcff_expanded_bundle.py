"""Declared fresh expanded-PCFF facade and separate-process bundle check.

No QM, dynamics or changes to source records. Always use a new output directory.
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from textwrap import dedent

import numpy as np

from island.chemistry import from_smiles
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldSources,
    create_evaluator,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.workflows.storage import json_bytes, publish


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = Path(args.output).resolve()
    root.mkdir(exist_ok=False)
    source = Path(args.source).resolve()
    publish(
        root / "declaration.json",
        json_bytes(
            {
                "smiles": "C1CCCCC1",
                "seed": 2026,
                "profile": "island_pcff_source_graph_v1",
                "lj": [0, 0, 1],
                "coulomb": [0, 0, 1],
                "atol": 1e-10,
                "rtol": 1e-12,
                "no_QM": True,
            }
        ),
    )
    m = from_smiles("C1CCCCC1", random_seed=2026)
    p = prepare_forcefield(
        m,
        ForceFieldRequest(
            "pcff",
            PCFFOptions(
                source,
                (0, 0, 1),
                (0, 0, 1),
                typing_profile="island_pcff_source_graph_v1",
            ),
        ),
    )
    e = create_evaluator(m, p)
    v = e.evaluate()
    bundle = save_prepared_forcefield(
        m, p, root / "bundle", sources=PreparedForceFieldSources(pcff_frc=source)
    )
    shutil.move(bundle, root / "relocated")
    code = """
    import json, sys
    from island.forcefields import load_prepared_forcefield,PreparedForceFieldSources,create_evaluator
    loaded=load_prepared_forcefield(sys.argv[1],sources=PreparedForceFieldSources(pcff_frc=sys.argv[2]))
    e=create_evaluator(loaded.system,loaded.prepared)
    a=e.evaluate_fresh()
    with e.open_session() as s: b=s.evaluate()
    json.dump({'identity':loaded.prepared.identity,'native':loaded.prepared.native_result.identity,'parameter':a.parameter_fingerprint,'model':a.model_fingerprint,'energy':a.potential_energy,'components':dict(a.energy_components),'forces':dict(a.forces),'session_energy':b.potential_energy,'session_forces':dict(b.forces)},open(sys.argv[3],'w'))
    """
    cmd = [
        sys.executable,
        "-c",
        dedent(code),
        str(root / "relocated"),
        str(source),
        str(root / "child.json"),
    ]
    publish(root / "command.json", json_bytes(cmd))
    subprocess.run(cmd, check=True, timeout=300)
    r = json.loads((root / "child.json").read_text())
    assert r["identity"] == p.identity and r["native"] == p.native_result.identity
    assert (
        r["model"] == v.model_fingerprint and r["parameter"] == v.parameter_fingerprint
    )
    np.testing.assert_allclose(
        [r["energy"], r["session_energy"]], v.potential_energy, atol=1e-10, rtol=1e-12
    )
    for i, f in v.forces.items():
        np.testing.assert_allclose(r["forces"][str(i)], f, atol=1e-10, rtol=1e-12)
        np.testing.assert_allclose(
            r["session_forces"][str(i)], f, atol=1e-10, rtol=1e-12
        )
    assert r["components"] == dict(v.energy_components)
    publish(
        root / "report.json",
        json_bytes(
            {
                "passed": True,
                "identity": p.identity,
                "native": p.native_result.identity,
                "child": r,
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
    )
    print("Fresh facade, relocation, separate-process bundle, fresh/session passed")


if __name__ == "__main__":
    main()
