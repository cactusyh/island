"""Real-source network assignment and exact split/continuous workflow acceptance.

Every family starts from the final J20/J21 graph. Failed preparation is reported
as blocked; it never enters dynamics. The separate fixed Foyer environment is
used explicitly for OPLS-AA. No source or force-field implementation is patched.
"""

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import traceback
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

from island.exceptions import PCFFError
from scripts.phase_4j31_contract import (
    CASES,
    CONTRACT,
    OPLS_SHA,
    PCFF_MESSAGE,
    PCFF_SHA,
    SUCCESS_FIELDS,
    case_errors,
    expected_typing,
    matrix_errors,
    typing_evidence,
)
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


class AssetUnavailable(RuntimeError):
    """Required environment or source is absent/invalid, never domain evidence."""


def check_assets(family, amberhome):
    modules = ["rdkit", "openmm", "parmed"]
    if family == "oplsaa":
        modules.append("foyer")
    for module in modules:
        if importlib.util.find_spec(module) is None:
            raise AssetUnavailable(f"Missing required dependency: {module}")
    if family in ("pcff", "oplsaa"):
        source, expected = (FRC, PCFF_SHA) if family == "pcff" else (XML, OPLS_SHA)
        try:
            digest = sha256(source.read_bytes()).hexdigest()
        except OSError as error:
            raise AssetUnavailable(
                f"Required source unavailable: {source}: {error}"
            ) from error
        if digest != expected:
            raise AssetUnavailable(
                f"Pinned {family} source checksum mismatch: {source}"
            )
    else:
        home = (
            Path(amberhome)
            if amberhome
            else Path(os.environ.get("AMBERHOME", ROOT / "miniforge3/envs/md_rg"))
        )
        for name in ("antechamber", "parmchk2", "tleap", "sqm"):
            executable = home / "bin" / name
            if not executable.is_file() or not os.access(executable, os.X_OK):
                raise AssetUnavailable(f"Missing required executable: {executable}")
        source = home / "dat/leap/parm" / (family + ".dat")
        if not source.is_file() or not source.stat().st_size:
            raise AssetUnavailable(f"Missing required parameter source: {source}")


def _expected_pcff_exception(error):
    if type(error) is not PCFFError or str(error) != PCFF_MESSAGE:
        return False
    tb = error.__traceback__
    while tb is not None:
        if (
            tb.tb_frame.f_globals.get("__name__") == "island.forcefields.pcff.automatic"
            and tb.tb_frame.f_code.co_name == "bridge_data"
        ):
            return True
        tb = tb.tb_next
    return False


def _write_result(output, result):
    (output / "result.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n"
    )
    return result


def worker(family, kind, output, amberhome):
    output.mkdir(parents=True, exist_ok=False)
    result = {
        "family": family,
        "fixture": kind,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    stage, typing = "assets", None
    try:
        check_assets(family, amberhome)
        stage = "fixture"
        fixture = build_network_fixture(kind)
        result.update(
            final_graph_identity=fixture.packed.graph.identity,
            crosslink_plan_identity=fixture.crosslinked.plan.identity,
            packing_identity=fixture.packed.plan.identity,
            independent_reference=independent_topology(fixture.packed.graph),
        )
        if family == "pcff":
            stage = "typing_probe"
            typing = typing_evidence(fixture, FRC)
        stage = "preparation"
        backend = prepare_network_backend(
            fixture, family, amberhome=amberhome, work_root=output, frc=FRC, xml=XML
        )
        if family == "pcff":
            raise ValueError(
                "Unexpected PCFF preparation success violates frozen J31 negative control"
            )
        result.update(
            status="complete",
            backend_identity=backend.identity,
            source=backend.config.payload["source"],
            typed_graph_identity=backend.typed_graph.identity,
            charge_identity=backend.graph_charges.identity,
            assignment_identity=backend.parameter_assignment.identity,
            independent_assignment_verification="native source term/charge/inventory validation; OPLS additionally verifies pinned upstream force content",
        )
        stage = "independent_verification"
        result["independent_assignment"] = independent_assignment_check(backend)
        stage = "execution"
        result["execution"] = execute(fixture, backend, output)
        result["dynamics_executed"] = True
        result["diagnostic"] = {"classification": "completed", "stage": "complete"}
        stage = "acceptance_contract"
        errors = case_errors(f"{family}/{kind}", result)
        if errors:
            raise ValueError("; ".join(errors))
    except Exception as error:  # noqa: BLE001 -- preserve diagnostics, then fail acceptance
        log = traceback.format_exc()
        for key in ("stdout", "stderr"):
            value = getattr(error, key, None)
            if value:
                log += str(value)
        (output / "failure.log").write_text(log)
        expected = (
            family == "pcff"
            and stage == "preparation"
            and typing is not None
            and expected_typing(kind, typing)
            and _expected_pcff_exception(error)
        )
        unavailable = isinstance(
            error, (AssetUnavailable, ImportError, FileNotFoundError, PermissionError)
        )
        classification = (
            "expected_typing_rejection"
            if expected
            else "asset_unavailable"
            if unavailable
            else "execution_failure"
            if stage in ("independent_verification", "execution", "acceptance_contract")
            else "unexpected_preparation_failure"
        )
        partial = {k: result.pop(k) for k in SUCCESS_FIELDS if k in result}
        if partial:
            result["partial_evidence"] = partial
        result.update(
            status="blocked" if expected or unavailable else "failed",
            dynamics_executed=result.get(
                "dynamics_executed", None if stage == "execution" else False
            ),
            reason=f"{type(error).__name__}: {error}",
            diagnostic={
                "classification": classification,
                "stage": "typing" if expected else stage,
                "exception_type": type(error).__name__,
                "message": str(error),
            },
        )
        if expected:
            result["diagnostic"].update(
                origin="island.forcefields.pcff.automatic.bridge_data",
                typing_evidence=typing,
            )
    result["acceptance_met"] = not case_errors(f"{family}/{kind}", result)
    return _write_result(output, result)


def run_matrix(output, amberhome, opls_python):
    output.mkdir(parents=True, exist_ok=False)
    environment = {**os.environ, "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}
    if amberhome:
        environment["PATH"] = (
            str(amberhome.resolve() / "bin") + os.pathsep + os.environ["PATH"]
        )
    report = {
        "contract": CONTRACT,
        "matrix": {},
        "case_failures": {},
        "acceptance_met": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }

    def persist():
        (output / "matrix.json").write_text(
            json.dumps(report, sort_keys=True, indent=2) + "\n"
        )

    persist()
    try:
        for name in CASES:
            family, kind = name.split("/")
            destination = output / f"{family}-{kind}"
            executable = str(opls_python) if family == "oplsaa" else sys.executable
            command = [
                executable,
                "-m",
                "scripts.validate_phase_4j31_network",
                "--output",
                str(destination),
                "--family",
                family,
                "--kind",
                kind,
            ]
            if amberhome:
                command.extend(["--amberhome", str(amberhome.resolve())])
            try:
                process = subprocess.run(
                    command,
                    env=environment,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                (output / f"{family}-{kind}.log").write_text(
                    process.stdout + process.stderr
                )
                result = json.loads((destination / "result.json").read_text())
                report["matrix"][name] = result
                errors = case_errors(name, result)
                if process.returncode or errors:
                    report["case_failures"][name] = {
                        "exit_code": process.returncode,
                        "contract_errors": errors,
                    }
            except Exception as error:  # noqa: BLE001 -- preserve diagnostics, then fail acceptance
                with (output / f"{family}-{kind}.log").open("a") as stream:
                    stream.write(traceback.format_exc())
                report["case_failures"][name] = {
                    "classification": "asset_unavailable"
                    if isinstance(error, OSError)
                    else "unexpected_preparation_failure",
                    "exception_type": type(error).__name__,
                    "message": str(error),
                }
            persist()
            print(
                f"{name}: {'unmet' if name in report['case_failures'] else 'accepted'}",
                flush=True,
            )
    except BaseException as error:
        report["interruption"] = {
            "exception_type": type(error).__name__,
            "message": str(error),
            "classification": "interrupted",
        }
        (output / "interruption.log").write_text(traceback.format_exc())
        raise
    finally:
        report["contract_errors"] = matrix_errors(report["matrix"])
        report["acceptance_met"] = (
            not report["case_failures"] and not report["contract_errors"]
        )
        persist()
    return 0 if report["acceptance_met"] else 1


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
        result = worker(args.family, args.kind, args.output, args.amberhome)
        print(json.dumps(result, sort_keys=True))
        return 0 if result["acceptance_met"] else 1
    return run_matrix(args.output, args.amberhome, args.opls_python)


if __name__ == "__main__":
    raise SystemExit(main())
