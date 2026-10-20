"""Real-source network assignment and exact split/continuous workflow acceptance.

Every family starts from the final J20/J21 graph. Failed preparation is reported
as blocked; it never enters dynamics. The separate fixed Foyer environment is
used explicitly for OPLS-AA. No source or force-field implementation is patched.
"""

import argparse
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from island.exceptions import IslandError
from scripts.phase_4j31_network_fixture import (
    FRC,
    ROOT,
    XML,
    build_network_fixture,
    independent_topology,
    prepare_network_backend,
    select_network,
)


def execute(fixture, backend, root):
    from island.analysis import (
        calculate_network_observables,
        load_network_analysis_bundle,
        load_verified_trajectory,
        save_network_analysis_bundle,
    )
    from island.periodic_md import (
        PeriodicMDWorkflowConfig,
        compare_periodic_md_workflows,
        load_periodic_md_workflow,
        save_periodic_md_workflow,
        start_periodic_md_workflow,
    )

    config = PeriodicMDWorkflowConfig(
        backend.config.payload["family"],
        steps=4,
        timestep_fs=0.01,
        minimize=False,
        seed=2026,
        velocity_seed=2027,
    )
    full, split, moved = (root / name for name in ("full", "split", "relocated"))
    start_periodic_md_workflow(backend, config, full)
    start_periodic_md_workflow(backend, config, split, steps=2)
    save_periodic_md_workflow(split, moved)
    identity = load_periodic_md_workflow(moved).identity
    subprocess.run(
        [
            sys.executable,
            "-c",
            "from island.periodic_md import resume_periodic_md_workflow; import sys; resume_periodic_md_workflow(sys.argv[1], expected_identity=sys.argv[2])",
            str(moved),
            identity,
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    comparison = compare_periodic_md_workflows(full, moved)
    if not comparison["passed"]:
        raise ValueError("Network split/uninterrupted checkpoints differ")
    trajectory = load_verified_trajectory(full)
    metadata = replace(
        trajectory.metadata,
        source_metadata={
            **trajectory.metadata.source_metadata,
            "crosslink_plan_identity": fixture.crosslinked.plan.identity,
        },
    )
    trajectory = replace(trajectory, metadata=metadata)
    selection = select_network(trajectory, fixture)
    values = calculate_network_observables(trajectory, selection)
    if any(v["status"] != "ok" for v in values.values()):
        raise ValueError("Executable fixture has blocked network observables")
    bundle = save_network_analysis_bundle(
        trajectory, selection, root / "analysis-original"
    )
    relocated = root / "analysis-relocated"
    bundle.rename(relocated)
    _, _, record = load_network_analysis_bundle(relocated)
    code = "from island.analysis import load_network_analysis_bundle; import sys; print(load_network_analysis_bundle(sys.argv[1],expected_identity=sys.argv[2])[2].to_json(),end='')"
    child = subprocess.run(
        [sys.executable, "-c", code, str(relocated), record.identity],
        check=True,
        capture_output=True,
        text=True,
    )
    if child.stdout != record.to_json():
        raise ValueError("Fresh process network reconstruction differs")
    return {
        "comparison": comparison,
        "analysis_record_identity": record.identity,
        "selection_identity": selection.identity,
        "frame_count": trajectory.frame_count,
        "frame_content_identity": selection.payload["frame_content_identity"],
        "same_process_relocation": "PASS",
        "separate_process_reconstruction": "PASS",
        "network_observables": values,
    }


def independent_assignment_check(backend):
    """Compare regenerated final-graph paths and independently parsed Amber forces."""
    from itertools import combinations, pairwise, permutations

    import numpy as np

    from island.forcefields import create_evaluator

    p = backend.payload
    data = p["numerical_data"]
    adjacency = {s["id"]: set() for s in p["final_graph"]["sites"]}
    for row in p["final_graph"]["bonds"]:
        a, b = row["sites"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    expected_bonds = {tuple(b["sites"]) for b in p["final_graph"]["bonds"]}
    expected_angles = {
        (a, b, c) for b in adjacency for a, c in combinations(sorted(adjacency[b]), 2)
    }
    expected_propers = set()
    # Exhaustive distinct four-site tuples are an independent bounded oracle.
    for path in permutations(adjacency, 4):
        if all(b in adjacency[a] for a, b in pairwise(path)):
            expected_propers.add(min(path, path[::-1]))
    for name, expected in (
        ("bonds", expected_bonds),
        ("angles", expected_angles),
        ("torsions", expected_propers),
    ):
        actual = [
            min(tuple(row["sites"]), tuple(reversed(row["sites"])))
            for row in data[name]
        ]
        if set(actual) != expected or len(actual) != len(expected):
            raise ValueError(f"Independent final-graph {name} coverage mismatch")
    result = {
        "bonds": len(expected_bonds),
        "angles": len(expected_angles),
        "proper_paths": len(expected_propers),
        "source_terms": "all required paths present",
    }
    prepared = backend._prepared
    if prepared._family in {"gaff", "gaff2"}:
        from scripts.validate_singlepoint_references import reference

        native = prepared.native_result
        source = Path(native.record["artifact_dir"]) / "result.prmtop"
        xyz = {i: prepared._system.coordinates.get(i).tolist() for i in adjacency}
        energy, forces, _ = reference(source, dict(native.imported_result.mapping), xyz)
        measured = create_evaluator(prepared._system, prepared).evaluate()
        error = max(
            float(np.max(np.abs(np.asarray(forces[i]) - measured.forces[i])))
            for i in adjacency
        )
        if (
            not np.isclose(energy, measured.potential_energy, rtol=1e-10, atol=1e-6)
            or error > 1e-5
        ):
            raise ValueError("Independent OpenMM AmberPrmtopFile energy/force mismatch")
        result.update(
            independent_reader="OpenMM AmberPrmtopFile versus ParmEd/ISLAND",
            energy_error_kj_mol=abs(energy - measured.potential_energy),
            max_force_error_kj_mol_angstrom=error,
            parmchk2_estimation_lines=native.record["parmchk2_estimation_lines"],
        )
    else:
        result["independent_reader"] = (
            "pinned upstream Foyer/OpenMM complete source force-content verification"
        )
    return result


def worker(family, kind, output, amberhome):
    fixture = build_network_fixture(kind)
    output.mkdir(parents=True, exist_ok=False)
    result = {
        "family": family,
        "fixture": kind,
        "final_graph_identity": fixture.packed.graph.identity,
        "crosslink_plan_identity": fixture.crosslinked.plan.identity,
        "packing_identity": fixture.packed.plan.identity,
        "independent_reference": independent_topology(fixture.packed.graph),
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        backend = prepare_network_backend(
            fixture, family, amberhome=amberhome, work_root=output
        )
    except (IslandError, ValueError, OSError) as error:
        result.update(
            status="blocked",
            reason=f"{type(error).__name__}: {error}",
            dynamics_executed=False,
        )
        (output / "result.json").write_text(
            json.dumps(result, sort_keys=True, indent=2) + "\n"
        )
        return result
    result.update(
        status="complete",
        backend_identity=backend.identity,
        source=backend.config.payload["source"],
        typed_graph_identity=backend.typed_graph.identity,
        charge_identity=backend.graph_charges.identity,
        assignment_identity=backend.parameter_assignment.identity,
        independent_assignment_verification="native source term/charge/inventory validation; OPLS additionally verifies pinned upstream force content",
    )
    # Execution failures are failures, never converted to unsupported-source claims.
    result["independent_assignment"] = independent_assignment_check(backend)
    result["execution"] = execute(fixture, backend, output)
    result["dynamics_executed"] = True
    (output / "result.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n"
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--family", choices=("pcff", "oplsaa", "gaff", "gaff2"))
    parser.add_argument(
        "--kind", choices=("branched", "crosslinked"), default="crosslinked"
    )
    parser.add_argument("--amberhome", type=Path)
    parser.add_argument(
        "--opls-python",
        type=Path,
        default=ROOT / "island-validation/phase4g1-env/bin/python",
    )
    args = parser.parse_args()
    if args.family:
        print(
            json.dumps(
                worker(args.family, args.kind, args.output, args.amberhome),
                sort_keys=True,
            )
        )
        return
    args.output.mkdir(parents=True, exist_ok=False)
    environment = {**os.environ, "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}
    if args.amberhome:
        environment["PATH"] = (
            str(args.amberhome.resolve() / "bin") + os.pathsep + os.environ["PATH"]
        )
    matrix = {}
    for family in ("pcff", "oplsaa", "gaff", "gaff2"):
        for kind in ("branched", "crosslinked"):
            executable = str(args.opls_python) if family == "oplsaa" else sys.executable
            command = [
                executable,
                "-m",
                "scripts.validate_phase_4j31_network",
                "--output",
                str(args.output / f"{family}-{kind}"),
                "--family",
                family,
                "--kind",
                kind,
            ]
            if args.amberhome:
                command.extend(["--amberhome", str(args.amberhome.resolve())])
            process = subprocess.run(
                command, env=environment, capture_output=True, text=True, check=False
            )
            (args.output / f"{family}-{kind}.log").write_text(
                process.stdout + process.stderr
            )
            if process.returncode:
                raise RuntimeError(
                    f"{family}/{kind} acceptance failed: {process.stderr[-3000:]}"
                )
            result = json.loads(
                (args.output / f"{family}-{kind}" / "result.json").read_text()
            )
            matrix[f"{family}/{kind}"] = result
            print(f"{family}/{kind}: {result['status']}", flush=True)
    report = {
        "source_paths": {"pcff": str(FRC), "oplsaa": str(XML)},
        "matrix": matrix,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    (args.output / "matrix.json").write_text(
        json.dumps(report, sort_keys=True, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
