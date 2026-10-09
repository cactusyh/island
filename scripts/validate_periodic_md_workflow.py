"""Reproduce J23 tests and collect exact workflow evidence, without LAMMPS MD.

Run with --run-tests for fresh execution, or collect existing pytest artifacts.
The receipt records software reproducibility only, never scientific acceptance.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

from island.periodic_md import load_periodic_md_workflow


def collect(root):
    cases = []
    for full in sorted(root.glob("test_four_family_exact*/full")):
        if full.parent.is_symlink():
            continue
        split = load_periodic_md_workflow(full.parent / "split").payload
        expected = load_periodic_md_workflow(full).payload
        actual = load_periodic_md_workflow(full.parent / "moved").payload
        assert actual["manifest"] == expected["manifest"]
        assert actual["frames"] == expected["frames"]
        assert split["frames"] == expected["frames"][: len(split["frames"])]
        m = expected["manifest"]
        frames = expected["frames"]
        cases.append(
            {
                "family": m["config"]["family"],
                "mode": m["config"]["mode"],
                "manifest_identity": m["identity"],
                "backend_identity": m["backend"]["identity"],
                "graph_identity": m["backend"]["final_graph_identity"],
                "packing_identity": m["packing_identity"],
                "source": m["backend"]["config"]["source"],
                "engine": m["engine"],
                "steps": frames[-1]["step"],
                "exact_frames_and_checkpoint_bytes": True,
                "split_boundary_equal": True,
                "relocated_child_process_resume": True,
                "initial_box_vectors_nm": frames[0]["box_vectors_nm"],
                "final_box_vectors_nm": frames[-1]["box_vectors_nm"],
                "final_checkpoint_sha256": frames[-1]["checkpoint_sha256"],
            }
        )
    assert {(c["family"], c["mode"]) for c in cases} == {
        (f, m)
        for f in ("PCFF", "OPLS-AA", "GAFF", "GAFF2")
        for m in ("nvt", "npt", "annealing")
    }
    minima = []
    for path in sorted(root.glob("test_four_family_minimization*/min")):
        if path.parent.is_symlink():
            continue
        run = load_periodic_md_workflow(path).payload
        diagnostic = json.loads((path / "minimization.json").read_text())["payload"]
        assert diagnostic["independent_context_verified"] is True
        minima.append({"family": run["manifest"]["config"]["family"], **diagnostic})
    assert {m["family"] for m in minima} == {"PCFF", "OPLS-AA", "GAFF", "GAFF2"}
    return {
        "schema": "island_phase_4j23_evidence_v1",
        "base": "7b5500dfff913d37b234ba723bca2e3ec36a82c6",
        "branch": "codex/phase-4j23-periodic-md-workflow",
        "openmm_workflow_execution": cases,
        "periodic_minimization": minima,
        "checkpoint_contract": "Exact Reference-platform binary checkpoint bytes and full frame values; same recorded OpenMM/environment required",
        "lammps_export_status": "J22 export unchanged; J23 does not execute LAMMPS",
        "lammps_runtime_status": "blocked: KSPACE package unavailable",
        "lammps_kspace_package": False,
        "lammps_openmm_periodic_comparison_claimed": False,
        "scientific_acceptance": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.run_tests:
        if args.artifacts.exists():
            raise FileExistsError(args.artifacts)
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                "tests/test_phase_4j23_periodic_md.py",
                "--basetemp=" + str(args.artifacts),
            ],
            env={**os.environ, "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"},
            capture_output=True,
            text=True,
            check=False,
        )
        args.artifacts.mkdir(parents=True, exist_ok=True)
        (args.artifacts / "pytest.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise SystemExit(result.returncode)
    receipt = collect(args.artifacts)
    receipt["retained_artifact_root"] = str(args.artifacts)
    receipt["artifact_hashes"] = {
        str(p.relative_to(args.artifacts)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(args.artifacts.glob("test_four_family_*/*/manifest.json"))
        if not p.parent.parent.is_symlink()
    }
    args.output.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(f"12 exact continuation cases and 4 independent minimizations: {args.output}")


if __name__ == "__main__":
    main()
