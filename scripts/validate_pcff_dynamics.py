"""Declared PCFF dynamics/relocated continuation using retained H5/H6 inputs.

Acceptance-only external records, not a new production workflow or restart API.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
from validate_oplsaa_dynamics import (
    KINDS,
    compare_frames,
    continuation_options,
    counted,
    digest,
    load_segment,
    options,
    read_record,
    require,
    save_segment,
)
from validate_pcff_singlepoint import (
    GROUPS,
    aa_reference_coefficients,
    data_sections,
    lammps_input,
    run,
)

from island.charge_references.records import unpack
from island.dynamics import (
    create_dynamics_checkpoint,
    initialize_velocities,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    run_langevin,
    run_nve,
    save_dynamics_checkpoint,
)
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import (
    load_pcff_model,
    load_pcff_parameters,
    load_pcff_source,
)
from island.workflows import storage
from island.workflows.bundle import (
    initialization_from,
    minimum_from,
    record,
    system_data,
    system_from,
)

CASES = ("butane", "ethanol", "PE_DP3", "PEO_DP3")
FILES = {
    "base.json",
    "system.json",
    "minimum.json",
    "initialization.json",
    "parameters.json",
    "model.json",
    "pcff.frc",
}


def verify_files(root, files):
    for name, checksum in files.items():
        require(
            digest(storage.child(root, name)) == checksum,
            "Artifact checksum differs: " + name,
        )


def check_minimum(base, system, minimum, parameters, model):
    base.validate()
    system.validate()
    minimum.validate_integrity()
    require(
        minimum.converged
        and minimum.final_evaluation_verified
        and minimum.final_fmax <= 0.1,
        "Minimum is not independently force converged",
    )
    parameters.validate_integrity(base)
    parameters.validate_integrity(system)
    model.validate_integrity(base)
    model.validate_integrity(system)
    require(
        system_data(minimum.to_system(base)) == system_data(system),
        "Applied minimum/optimized system mismatch",
    )
    policy = unpack(model.json_text)["special_pairs"]
    require(
        policy["lj"] == [0, 0, 1] and policy["coulomb"] == [0, 0, 1],
        "H5 special-pair policy changed",
    )
    evaluator = PCFFSinglePointEvaluator(system, model)
    require(
        evaluator.parameter_fingerprint
        == minimum.final_evaluation.parameter_fingerprint
        == minimum.initial_evaluation.parameter_fingerprint,
        "Minimum parameter identity differs",
    )
    require(
        evaluator.model_fingerprint
        == minimum.final_evaluation.model_fingerprint
        == minimum.initial_evaluation.model_fingerprint,
        "Minimum model identity differs",
    )
    return evaluator


def retained(h5, h6, frc, name):
    source = load_pcff_source(frc)
    base = system_from(read_record(h6 / name / "starting-system.json"))
    system = system_from(read_record(h6 / name / "optimized-system.json"))
    minimum = minimum_from(read_record(h6 / name / "minimum.json"))
    parameters = load_pcff_parameters(h5 / name / "parameters.json", source)
    model = load_pcff_model(h5 / name / "model-0.json", parameters)
    evaluator = check_minimum(base, system, minimum, parameters, model)
    row = next(
        r for r in storage.read_json(h6 / "report.json")["cases"] if r["case"] == name
    )
    require(
        row["passed"]
        and row["model_identity"] == evaluator.model_fingerprint
        and row["parameter_identity"] == evaluator.parameter_fingerprint,
        "H6 acceptance identity differs",
    )
    require(
        row["initial_coordinate_fingerprint"]
        == minimum.initial_evaluation.coordinate_fingerprint
        and row["final_coordinate_fingerprint"]
        == minimum.final_evaluation.coordinate_fingerprint,
        "H6 coordinate identity differs",
    )
    require(
        all(
            row[k] == getattr(minimum, k)
            for k in (
                "initial_energy",
                "final_energy",
                "iterations",
                "evaluations",
                "termination_reason",
            )
        ),
        "H6 summary differs",
    )
    return base, system, parameters, model, minimum, source, evaluator


def declaration(args):
    from island.dynamics._checkpoint_data import environment

    evidence_dir = Path(__file__).resolve().parents[1] / "docs/evidence"
    h5 = storage.read_json(evidence_dir / "phase_4h5.json")
    h6 = storage.read_json(evidence_dir / "phase_4h6.json")
    verify_files(args.retained, h5["files"])
    verify_files(args.minima, h6["files"])
    require(
        digest(args.lammps) == h5["upstream"]["lammps_sha256"],
        "Pinned LAMMPS executable mismatch",
    )
    source = load_pcff_source(args.source)
    cases = {}
    for name in CASES:
        _, s, _p, _m, minimum, _, e = retained(
            args.retained, args.minima, args.source, name
        )
        cases[name] = {
            "sites": s.number_of_sites,
            "parameter_identity": e.parameter_fingerprint,
            "model_identity": e.model_fingerprint,
            "starting_coordinates": minimum.final_evaluation.coordinate_fingerprint,
        }
    return {
        "schema": "island_pcff_dynamics_acceptance_v1",
        "cases": cases,
        "h5_files": h5["files"],
        "h6_files": h6["files"],
        "source": source.identity,
        "lammps_sha256": digest(args.lammps),
        "options": {
            kind: {str(n): asdict(options(kind, n)) for n in (100, 200)}
            for kind in KINDS
        },
        "continuation_options": asdict(continuation_options()),
        "velocity_temperature_kelvin": 300.0,
        "velocity_seed": 78123,
        "split_atol": 1e-10,
        "split_rtol": 1e-12,
        "reference_atol_energy": 1e-5,
        "reference_atol_force": 1e-5,
        "reference_rtol": 2e-10,
        "reference_steps": [0, 100, 200],
        "retained_steps": list(range(0, 201, 20)),
        "lj": [0, 0, 1],
        "coulomb": [0, 0, 1],
        "units": {
            "coordinates": "angstrom",
            "velocities": "angstrom/ps",
            "energy": "kJ/mol",
            "force": "kJ/(mol*angstrom)",
            "time": "ps",
        },
        "reference": "Independent H5 converter inventories + H5/H6 independently reconstructed AA equilibrium roles",
        "environment": environment(),
        "platform": "Reference",
        "worker_timeout_seconds": 300,
        "outer_retries": 0,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def bundle(directory, base, system, parameters, model, minimum, source, initialization):
    directory.mkdir()
    for name, value in {
        "base.json": system_data(base),
        "system.json": system_data(system),
        "minimum.json": record(minimum),
        "initialization.json": record(initialization),
    }.items():
        storage.publish(directory / name, storage.json_bytes(storage.encode(value)))
    storage.publish(directory / "parameters.json", parameters.json_text.encode())
    storage.publish(directory / "model.json", model.json_text.encode())
    storage.publish(directory / "pcff.frc", source.raw)
    storage.publish(
        directory / "inputs.json",
        storage.json_bytes(
            {
                "schema": "pcff_acceptance_inputs_v1",
                "files": {name: digest(directory / name) for name in sorted(FILES)},
                "parameter_identity": minimum.final_evaluation.parameter_fingerprint,
                "model_identity": minimum.final_evaluation.model_fingerprint,
                "source": source.identity,
            }
        ),
    )


def reconstruct(directory):
    manifest = storage.read_json(directory / "inputs.json")
    require(
        manifest["schema"] == "pcff_acceptance_inputs_v1"
        and set(manifest["files"]) == FILES,
        "Unsupported external input inventory",
    )
    verify_files(directory, manifest["files"])
    source = load_pcff_source(directory / "pcff.frc")
    require(source.identity == manifest["source"], "Source identity differs")
    base = system_from(read_record(directory / "base.json"))
    system = system_from(read_record(directory / "system.json"))
    minimum = minimum_from(read_record(directory / "minimum.json"))
    parameters = load_pcff_parameters(directory / "parameters.json", source)
    model = load_pcff_model(directory / "model.json", parameters)
    evaluator = check_minimum(base, system, minimum, parameters, model)
    require(
        evaluator.parameter_fingerprint == manifest["parameter_identity"]
        and evaluator.model_fingerprint == manifest["model_identity"],
        "External model identity differs",
    )
    init = initialization_from(read_record(directory / "initialization.json"))
    from island.dynamics.thermal import mass_inventory
    from island.minimization.models import system_identity

    init.validate_integrity()
    require(
        dict(init.masses) == mass_inventory(system)
        and init.system_fingerprint == system_identity(system),
        "Initialization/system identity differs",
    )
    return system, evaluator, init


def validate_origin(checkpoint, system, initialization):
    """Bind saved setup to original trajectory without redrawing velocities."""
    from island.dynamics._checkpoint_data import inventory

    initial = checkpoint.payload["origin"]["initial_state"]
    require(
        {sid: tuple(v) for sid, v in inventory(initial["coordinates"]).items()}
        == {sid: tuple(system.coordinates.get(sid)) for sid in system.topology.sites},
        "Checkpoint origin/optimized coordinates differ",
    )
    require(
        {sid: tuple(v) for sid, v in inventory(initial["velocities"]).items()}
        == dict(initialization.velocities)
        and initial["velocity_fingerprint"] == initialization.velocity_fingerprint,
        "Checkpoint origin/initialization velocities differ",
    )


def worker(directory):
    request = storage.read_json(directory / "resume-request.json")
    require(
        request["inputs_sha256"] == digest(directory / "inputs.json"),
        "Resume input manifest mismatch",
    )
    require(
        request["checkpoint_sha256"] == digest(directory / "checkpoint.json"),
        "Resume checkpoint file mismatch",
    )
    require(
        request["parent_pid"] != os.getpid(), "Resume did not enter a separate process"
    )
    s, e, initialization = reconstruct(
        directory
    )  # no construction, typing, minimization or new initialization
    cp = load_dynamics_checkpoint(directory / "checkpoint.json")
    validate_origin(cp, s, initialization)
    require(cp.absolute_step == 100, "Expected the declared boundary at step 100")
    require(
        request["options"] == asdict(continuation_options()),
        "Resume budgets differ from declared experiment",
    )
    result, contexts = counted(
        e, lambda session: resume_dynamics(cp, s, session, continuation_options())
    )
    save_segment(result, directory / "continued.json")
    status = {
        "pid": os.getpid(),
        "evaluations": result.evaluations,
        "contexts": contexts,
        "completed": result.completed,
        "loaded_optional_modules": {
            n: n in sys.modules for n in ("foyer", "parmed", "scipy", "rdkit")
        },
    }
    if result.completed:
        save_dynamics_checkpoint(
            create_dynamics_checkpoint(result), directory / "continued-checkpoint.json"
        )
    storage.publish(directory / "worker.json", storage.json_bytes(status))
    return 0 if result.completed else 1


def run_case(root, name, kind, inputs, h5, executable):
    base, s, p, model, m, source, e = inputs
    init = initialize_velocities(s, temperature_kelvin=300.0, seed=78123)
    row = {"case": name, "integrator": kind, "status": "failed"}
    initial = root / f"{name}-{kind}-before-relocation"
    moved = root / f"{name}-{kind}"
    bundle(initial, base, s, p, model, m, source, init)
    method = run_nve if kind == "nve" else run_langevin
    direct, direct_contexts = counted(
        e, lambda session: method(s, session, init.velocities, options(kind, 200))
    )
    storage.publish(
        initial / "direct.json", storage.json_bytes(storage.encode(record(direct)))
    )
    require(
        direct.completed and direct.final_evaluation_verified,
        "Direct dynamics gate failed",
    )
    whole, whole_contexts = counted(
        e,
        lambda session: run_dynamics_segment(
            s, session, init.velocities, options(kind, 200)
        ),
    )
    save_segment(whole, initial / "whole.json")
    require(whole.completed, "Uninterrupted segment failed")
    compare_frames(direct.frames, whole.frames)
    first, first_contexts = counted(
        e,
        lambda session: run_dynamics_segment(
            s, session, init.velocities, options(kind, 100)
        ),
    )
    save_segment(first, initial / "first.json")
    require(first.completed, "First split segment failed")
    cp = create_dynamics_checkpoint(first)
    save_dynamics_checkpoint(cp, initial / "checkpoint.json")
    storage.publish(
        initial / "resume-request.json",
        storage.json_bytes(
            {
                "inputs_sha256": digest(initial / "inputs.json"),
                "checkpoint_sha256": digest(initial / "checkpoint.json"),
                "parent_pid": os.getpid(),
                "options": asdict(continuation_options()),
            }
        ),
    )
    shutil.move(str(initial), str(moved))
    child = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--resume-directory",
            str(moved.resolve()),
        ],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    storage.publish(moved / "worker.log", (child.stdout + child.stderr).encode())
    require(
        child.returncode == 0,
        f"Resume subprocess failed: exit {child.returncode}; see worker.log",
    )
    last = load_segment(moved / "continued.json")
    require(last.completed, "Continued segment failed")
    status = storage.read_json(moved / "worker.json")
    require(
        status["pid"] != os.getpid()
        and not any(status["loaded_optional_modules"].values()),
        "Resume isolation gate failed",
    )
    require(
        first.final_state.coordinates == last.frames[0].coordinates
        and first.final_state.velocities == last.frames[0].velocities
        and first.final_state.step == last.frames[0].step
        and first.final_state.time_ps == last.frames[0].time_ps,
        "Shared synchronized boundary differs",
    )
    compare_frames([first.final_state], [last.frames[0]])
    combined = (
        first.frames + last.frames[1:]
    )  # retain earlier boundary after comparison
    require(
        [f.step for f in combined] == list(range(0, 201, 20)),
        "Unexpected retained schedule",
    )
    errors = compare_frames(whole.frames, combined)
    wp, lp = whole.payload, last.payload
    require(wp["origin"] == lp["origin"], "Original trajectory/reference changed")
    require(wp["rng"] == lp["rng"], "Full PCG64 state differs")
    require(
        lp["lineage"][:-1] == first.payload["lineage"] and len(lp["lineage"]) == 2,
        "Accepted lineage changed",
    )
    require(
        sum(r["end_step"] - r["start_step"] for r in lp["lineage"]) == 200
        and lp["lineage"][-1]["parent_checksum"] == cp.content_checksum,
        "Cumulative steps/parent lineage differ",
    )
    require(
        wp["counters"]["normal_draws"] == lp["counters"]["normal_draws"],
        "RNG draw accounting differs",
    )
    require(
        whole.evaluations == direct.evaluations == 202
        and first.evaluations == last.evaluations == 102
        and lp["counters"]["evaluations"] == 204,
        "Evaluator accounting differs",
    )
    require(
        direct_contexts == 2
        and whole_contexts == first_contexts == status["contexts"] == 3,
        "Unexpected Context construction count",
    )
    if kind == "nve":
        require(
            lp["max_abs_energy_deviation"] <= 1.0
            and wp["max_abs_energy_deviation"] <= 1.0,
            "NVE guard violated",
        )
    raw = h5 / name / "reference.data"
    overrides = aa_reference_coefficients(data_sections(raw))
    comparisons = []
    for frame in (whole.frames[0], whole.frames[5], whole.frames[-1]):
        folder = moved / f"reference-{frame.step}"
        folder.mkdir()
        (folder / "reference.data").write_bytes(raw.read_bytes())
        ids = sorted(frame.coordinates)
        xyz = np.array([frame.coordinates[i] for i in ids])
        (folder / "input.lmp").write_text(
            lammps_input(xyz, [0, 0, 1], [0, 0, 1], overrides)
        )
        storage.publish(folder / "aa-overrides.json", storage.json_bytes(overrides))
        seconds = run([executable, "-in", "input.lmp"], folder, "lammps")
        line = [
            x
            for x in (folder / "lammps.log").read_text().splitlines()
            if x.startswith("ENERGIES ")
        ][-1]
        energies = np.array([float(x) for x in line.split()[1:]]) * 4.184
        dump = (
            (folder / "forces.dump")
            .read_text()
            .split("ITEM: ATOMS id x y z fx fy fz\n")[1]
        )
        forces = (
            np.array(
                [
                    [float(v) for v in line.split()[4:7]]
                    for line in dump.strip().splitlines()
                ]
            )
            * 4.184
        )
        ev = frame.evaluation
        calculated = np.array(
            [ev.potential_energy]
            + [sum(ev.energy_components[k] for k in keys) for keys in GROUPS.values()]
        )
        actual = np.array([ev.forces[i] for i in ids])
        comparison = {
            "step": frame.step,
            "energy_and_component_max_abs_error": float(
                np.max(abs(calculated - energies))
            ),
            "force_max_abs_error": float(np.max(abs(actual - forces))),
            "island_energies": calculated.tolist(),
            "lammps_energies": energies.tolist(),
            "seconds": seconds,
        }
        storage.publish(folder / "comparison.json", storage.json_bytes(comparison))
        require(
            np.allclose(calculated, energies, atol=1e-5, rtol=2e-10),
            "Reference energy/components differ",
        )
        require(
            np.allclose(actual, forces, atol=1e-5, rtol=2e-10),
            "Reference forces differ",
        )
        comparisons.append(comparison)
    reconstruct(moved)  # all external input hashes remain unchanged
    row.update(
        status="passed",
        completed_steps=last.final_state.step,
        time_ps=last.final_state.time_ps,
        evaluator_calls={
            "direct": direct.evaluations,
            "whole": whole.evaluations,
            "first": first.evaluations,
            "resumed": last.evaluations,
            "split_cumulative": lp["counters"]["evaluations"],
        },
        contexts={
            "direct": direct_contexts,
            "whole": whole_contexts,
            "first": first_contexts,
            "resumed": status["contexts"],
            "external_lammps_processes": 3,
        },
        split_max_abs_errors=errors,
        rng_equal=True,
        origin_equal=True,
        maximum_total_energy_deviation={
            "whole": wp["max_abs_energy_deviation"],
            "split": lp["max_abs_energy_deviation"],
            "nve_guard_applied": kind == "nve",
        },
        reference=comparisons,
        reference_inspection="Independent retained H5 converter data with explicit AA overrides",
        worker=status,
        trajectory_identity=wp["origin"]["trajectory_fingerprint"],
        parameter_identity=e.parameter_fingerprint,
        final_coordinate_fingerprint=last.final_state.coordinate_fingerprint,
        final_checkpoint_sha256=digest(moved / "continued-checkpoint.json"),
        normal_draws=lp["counters"]["normal_draws"],
    )
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retained", type=Path)
    parser.add_argument("--minima", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--lammps", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--execute-declared", action="store_true")
    parser.add_argument("--resume-directory", type=Path)
    args = parser.parse_args()
    if args.resume_directory:
        return worker(args.resume_directory)
    if any(
        getattr(args, k) is None
        for k in ("retained", "minima", "source", "lammps", "output")
    ):
        parser.error("--retained --minima --source --lammps --output required")
    plan = declaration(args)
    if not args.execute_declared:
        args.output.mkdir(parents=True, exist_ok=False)
        storage.publish(args.output / "declaration.json", storage.json_bytes(plan))
        return 0
    require(
        storage.read_json(args.output / "declaration.json") == plan,
        "Declaration/input mismatch",
    )
    require(
        {p.name for p in args.output.iterdir()} == {"declaration.json"},
        "Execution directory must contain only declaration",
    )
    outcomes = []
    for name in CASES:
        for kind in KINDS:
            start = time.perf_counter()
            try:
                result = run_case(
                    args.output,
                    name,
                    kind,
                    retained(args.retained, args.minima, args.source, name),
                    args.retained,
                    args.lammps.resolve(),
                )
            except Exception as exc:  # noqa: BLE001 -- retain failures, never retry
                result = {
                    "case": name,
                    "integrator": kind,
                    "status": "failed",
                    "failure": f"{type(exc).__name__}: {exc}",
                }
            result["seconds"] = time.perf_counter() - start
            outcomes.append(result)
            print(name, kind, result["status"], result.get("failure", ""), flush=True)
            storage.publish(
                args.output / "outcomes.json",
                storage.json_bytes({"cases": outcomes, "complete": False}),
                replace=True,
            )
    unchanged = declaration(args) == plan
    passed = (
        unchanged
        and len(outcomes) == 8
        and all(r["status"] == "passed" for r in outcomes)
    )
    storage.publish(
        args.output / "outcomes.json",
        storage.json_bytes(
            {
                "cases": outcomes,
                "complete": passed,
                "inputs_unchanged": unchanged,
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
        replace=True,
    )
    storage.publish(
        args.output / "files.json",
        storage.json_bytes(
            {
                str(p.relative_to(args.output)): digest(p)
                for p in args.output.rglob("*")
                if p.is_file()
            }
        ),
    )
    return int(not passed)


if __name__ == "__main__":
    raise SystemExit(main())
