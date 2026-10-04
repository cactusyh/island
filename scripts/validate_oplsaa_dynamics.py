"""Declare then validate OPLS dynamics and relocated, separate-process continuation.

Data files here are acceptance inputs/results, not a new workflow or checkpoint API.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

import numpy as np

from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import (
    DynamicsOptions,
    DynamicsSegment,
    DynamicsSegmentOptions,
    LangevinOptions,
    create_dynamics_checkpoint,
    initialize_velocities,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    run_langevin,
    run_nve,
    save_dynamics_checkpoint,
)
from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.forcefields.oplsaa import OPLSParameterizationResult, load_oplsaa_source
from island.workflows import storage
from island.workflows.bundle import (
    initialization_from,
    minimum_from,
    record,
    system_data,
    system_from,
)

CASES = ("butane", "ethanol", "pe3", "peo3")
KINDS = ("nve", "baoab")


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(path):
    return storage.checksum(Path(path).read_bytes())


def read_record(path):
    return storage.decode(storage.read_json(path))


def options(kind, steps):
    common = {
        "timestep_fs": 0.1,
        "steps": steps,
        "max_evaluations": steps + 2,
        "max_frames": steps // 20 + 1,
        "recording_interval": 20,
    }
    return (
        DynamicsOptions(**common, max_energy_deviation=1.0)
        if kind == "nve"
        else LangevinOptions(
            **common,
            temperature_kelvin=300.0,
            friction_per_ps=5.0,
            thermostat_seed=99181,
        )
    )


def continuation_options():
    return DynamicsSegmentOptions(100, 102, 6, 20)


def save_segment(segment, path):
    segment.validate_integrity()
    storage.publish(
        path,
        storage.json_bytes(
            {
                "payload_json": segment.payload_json,
                "content_checksum": segment.content_checksum,
            }
        ),
    )


def load_segment(path):
    return DynamicsSegment(**storage.read_json(path))


def retained(g2, g3, xml, name):
    source = load_oplsaa_source(xml)
    report = storage.read_json(g3 / "outcomes.json")
    historical = storage.read_json(
        Path(__file__).resolve().parents[1] / "docs/phase_4g3_acceptance.json"
    )
    require(report == historical, "Retained 4G3 outcomes differ from reviewed evidence")
    row = next(r for r in report["cases"] if r["case"] == name)
    require(row["status"] == "passed", "Required minimum was not accepted")
    require(
        report["declaration"]["source"] == source.identity, "Source differs from 4G3"
    )
    for filename, checksum in report["declaration"]["cases"][name]["files"].items():
        require(
            digest(storage.child(g2, filename)) == checksum,
            "4G2 source file checksum differs",
        )
    base = system_from(read_record(g2 / f"{name}-system.json"))
    parameters = OPLSParameterizationResult(
        (g2 / f"{name}-parameters.json").read_text()
    )
    minimum = minimum_from(read_record(g3 / f"{name}-minimum.json"))
    system = system_from(read_record(g3 / f"{name}-optimized.json"))
    parameters.validate_integrity(base, source)
    parameters.validate_integrity(system, source)
    require(
        parameters.identity
        == row["parameter_identity"]
        == minimum.final_evaluation.parameter_fingerprint,
        "Minimum/parameter mismatch",
    )
    require(
        system_data(minimum.to_system(base)) == system_data(system),
        "Saved optimized system differs from verified minimum",
    )
    require(
        coordinate_hash(minimum.coordinates) == row["final_coordinate_fingerprint"],
        "Minimum final coordinates differ from 4G3 outcome",
    )
    require(
        coordinate_hash(minimum.initial_coordinates)
        == row["initial_coordinate_fingerprint"],
        "Minimum initial coordinates differ",
    )
    require(
        minimum.final_energy == row["final_energy"]
        and minimum.evaluations == row["evaluations"]
        and minimum.iterations == row["iterations"],
        "Minimum summary differs",
    )
    return base, system, parameters, minimum, source


def declaration(g2, g3, xml):
    from island.dynamics._checkpoint_data import environment

    cases = {}
    for name in CASES:
        _, s, p, m, _ = retained(g2, g3, xml, name)
        cases[name] = {
            "g2_files": {
                f"{name}-{suffix}.json": digest(g2 / f"{name}-{suffix}.json")
                for suffix in ("system", "parameters")
            },
            "g3_files": {
                f"{name}-{suffix}.json": digest(g3 / f"{name}-{suffix}.json")
                for suffix in ("minimum", "optimized")
            },
            "parameter_identity": p.identity,
            "starting_coordinates": coordinate_hash(m.coordinates),
            "sites": s.number_of_sites,
        }
    return {
        "schema": "island_opls_dynamics_acceptance_v1",
        "cases": cases,
        "source": load_oplsaa_source(xml).identity,
        "minimum_outcomes_sha256": digest(g3 / "outcomes.json"),
        "options": {
            kind: {str(n): asdict(options(kind, n)) for n in (100, 200)}
            for kind in KINDS
        },
        "continuation_options": asdict(continuation_options()),
        "velocity_temperature_kelvin": 300.0,
        "velocity_seed": 78123,
        "split_atol": 1e-10,
        "split_rtol": 1e-12,
        "units": {
            "coordinates": "angstrom",
            "velocities": "angstrom/ps",
            "energy": "kJ/mol",
            "force": "kJ/(mol*angstrom)",
            "time": "ps",
        },
        "reference_atol_energy": 1e-5,
        "reference_atol_force": 1e-5,
        "reference_rtol": 2e-10,
        "reference_steps": [0, 100, 200],
        "retained_steps": list(range(0, 201, 20)),
        "reference": "4G2 independent Foyer.apply -> ParmEd geometric mixing",
        "environment": environment(),
        "platform": "Reference",
        "worker_timeout_seconds": 300,
        "outer_retries": 0,
        "expected_rejection": "ps3_native_charge",
        "gates": [
            "all eight trajectories complete",
            "split consistency",
            "exact PCG64",
            "independent reference",
            "NVE drift <= 1 kJ/mol",
            "source nonmutation",
            "fresh verification",
            "separate process and relocation",
        ],
    }


def bundle(directory, base, system, parameters, minimum, source, initialization):
    directory.mkdir()
    objects = {
        "base.json": system_data(base),
        "system.json": system_data(system),
        "minimum.json": record(minimum),
        "initialization.json": record(initialization),
    }
    for name, value in objects.items():
        storage.publish(directory / name, storage.json_bytes(storage.encode(value)))
    storage.publish(directory / "parameters.json", parameters.json_text.encode())
    storage.publish(directory / "oplsaa.xml", source.xml)
    files = {
        n: digest(directory / n) for n in (*objects, "parameters.json", "oplsaa.xml")
    }
    storage.publish(
        directory / "inputs.json",
        storage.json_bytes(
            {
                "schema": "opls_acceptance_inputs_v1",
                "files": files,
                "parameter_identity": parameters.identity,
                "source": source.identity,
            }
        ),
    )


def reconstruct(directory):
    manifest = storage.read_json(directory / "inputs.json")
    require(
        manifest["schema"] == "opls_acceptance_inputs_v1",
        "Unsupported acceptance inputs",
    )
    require(
        set(manifest["files"])
        == {
            "base.json",
            "system.json",
            "minimum.json",
            "initialization.json",
            "parameters.json",
            "oplsaa.xml",
        },
        "External file inventory differs",
    )
    for name, checksum in manifest["files"].items():
        require(
            digest(storage.child(directory, name)) == checksum,
            "External artifact checksum differs",
        )
    source = load_oplsaa_source(directory / "oplsaa.xml")
    require(source.identity == manifest["source"], "External source differs")
    system = system_from(read_record(directory / "system.json"))
    base = system_from(read_record(directory / "base.json"))
    minimum = minimum_from(read_record(directory / "minimum.json"))
    require(
        system_data(minimum.to_system(base)) == system_data(system),
        "External starting system/minimum mismatch",
    )
    parameters = OPLSParameterizationResult((directory / "parameters.json").read_text())
    parameters.validate_integrity(system, source)
    require(
        parameters.identity
        == manifest["parameter_identity"]
        == minimum.final_evaluation.parameter_fingerprint,
        "External parameter identity differs",
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
    return system, OPLSSinglePointEvaluator(system, parameters, source), init


def counted(evaluator, call):
    import openmm as mm

    constructor = mm.Context
    count = []

    def context(*args, **kwargs):
        obj = constructor(*args, **kwargs)
        count.append(1)
        return obj

    with patch.object(mm, "Context", context), evaluator.open_session() as session:
        result = call(session)
    return result, len(count)


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
    s, e, _ = reconstruct(
        directory
    )  # no construction, typing, minimization or new initialization
    cp = load_dynamics_checkpoint(directory / "checkpoint.json")
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


def compare_frames(left, right):
    require(len(left) == len(right), "Retained frame count differs")
    maxima = {
        k: 0.0
        for k in (
            "coordinates",
            "velocities",
            "potential_energy",
            "kinetic_energy",
            "total_energy",
            "forces",
            "components",
            "time_ps",
        )
    }
    for a, b in zip(left, right, strict=True):
        require(a.step == b.step, "Absolute step differs")
        require(
            set(a.coordinates) == set(b.coordinates), "Stable-site coverage differs"
        )
        ids = sorted(a.coordinates)
        pairs = {
            "coordinates": (
                [a.coordinates[i] for i in ids],
                [b.coordinates[i] for i in ids],
            ),
            "velocities": (
                [a.velocities[i] for i in ids],
                [b.velocities[i] for i in ids],
            ),
            "forces": (
                [a.evaluation.forces[i] for i in ids],
                [b.evaluation.forces[i] for i in ids],
            ),
            "potential_energy": (
                a.evaluation.potential_energy,
                b.evaluation.potential_energy,
            ),
            "kinetic_energy": (a.kinetic_energy, b.kinetic_energy),
            "total_energy": (a.total_energy, b.total_energy),
            "time_ps": (a.time_ps, b.time_ps),
            "components": (
                list(a.evaluation.energy_components.values()),
                [
                    b.evaluation.energy_components[k]
                    for k in a.evaluation.energy_components
                ],
            ),
        }
        for key, (x, y) in pairs.items():
            require(
                np.allclose(x, y, atol=1e-10, rtol=1e-12),
                f"Split {key} tolerance exceeded at step {a.step}",
            )
            maxima[key] = max(
                maxima[key], float(np.max(np.abs(np.asarray(x) - np.asarray(y))))
            )
        require(
            a.evaluation.model_fingerprint == b.evaluation.model_fingerprint
            and a.evaluation.parameter_fingerprint
            == b.evaluation.parameter_fingerprint,
            "Frame model identity differs",
        )
    return maxima


def run_case(root, name, kind, inputs):
    from validate_oplsaa_energy import reference, reference_value

    base, s, p, m, source = inputs
    init = initialize_velocities(s, temperature_kelvin=300.0, seed=78123)
    e = OPLSSinglePointEvaluator(s, p, source)
    row = {"case": name, "integrator": kind, "status": "failed"}
    initial = root / f"{name}-{kind}-before-relocation"
    moved = root / f"{name}-{kind}"
    bundle(initial, base, s, p, m, source, init)
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
    ref, inspection = reference(s, moved / "oplsaa.xml")
    comparisons = []
    for frame in (whole.frames[0], whole.frames[5], whole.frames[-1]):
        ids = sorted(frame.coordinates)
        energy, forces, components = reference_value(
            ref, np.array([frame.coordinates[i] for i in ids])
        )
        ev = frame.evaluation
        actual = {k: ev.energy_components[k] for k in ("bond", "angle", "rb_proper")}
        actual["nonbonded"] = (
            ev.energy_components["lj"] + ev.energy_components["coulomb"]
        )
        require(
            np.isclose(ev.potential_energy, energy, atol=1e-5, rtol=2e-10),
            "Independent reference energy differs",
        )
        require(
            np.allclose([ev.forces[i] for i in ids], forces, atol=1e-5, rtol=2e-10),
            "Independent reference forces differ",
        )
        require(
            all(
                np.isclose(actual[k], components[k], atol=1e-5, rtol=2e-10)
                for k in components
            ),
            "Reference components differ",
        )
        comparisons.append(
            {
                "step": frame.step,
                "energy_abs_error": abs(ev.potential_energy - energy),
                "force_max_abs_error": float(
                    np.max(np.abs(np.array([ev.forces[i] for i in ids]) - forces))
                ),
                "component_abs_errors": {
                    k: abs(actual[k] - components[k]) for k in components
                },
            }
        )
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
            "oracle": 3,
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
        reference_inspection=inspection,
        worker=status,
        trajectory_identity=wp["origin"]["trajectory_fingerprint"],
        parameter_identity=p.identity,
        final_coordinate_fingerprint=last.final_state.coordinate_fingerprint,
        final_checkpoint_sha256=digest(moved / "continued-checkpoint.json"),
        normal_draws=lp["counters"]["normal_draws"],
    )
    return row


def execute(args, plan):
    from validate_oplsaa_energy import make_system

    from island.forcefields.oplsaa import parameterize_oplsaa

    outcomes = []
    try:
        parameterize_oplsaa(
            make_system("[*:1]CC(c1ccccc1)[*:2]", 3), load_oplsaa_source(args.xml)
        )
    except Exception as error:  # noqa: BLE001 -- retain expected and unexpected failures distinctly
        outcomes.append(
            {
                "case": "ps3",
                "status": "expected_rejection"
                if "Native component neutrality failed" in str(error)
                else "failed",
                "message": str(error),
            }
        )
    else:
        outcomes.append(
            {"case": "ps3", "status": "failed", "message": "Missing expected rejection"}
        )
    for name in CASES:
        for kind in KINDS:
            start = time.perf_counter()
            try:
                result = run_case(
                    args.output,
                    name,
                    kind,
                    retained(args.retained, args.minima, args.xml, name),
                )
            except Exception as error:  # noqa: BLE001 -- preserve every experiment failure, never retry
                result = {
                    "case": name,
                    "integrator": kind,
                    "status": "failed",
                    "failure": f"{type(error).__name__}: {error}",
                }
            result["seconds"] = time.perf_counter() - start
            outcomes.append(result)
            print(name, kind, result["status"], result.get("failure", ""), flush=True)
            storage.publish(
                args.output / "outcomes.json",
                storage.json_bytes(
                    {"declaration": plan, "cases": outcomes, "complete": False}
                ),
                replace=True,
            )
    unchanged = declaration(args.retained, args.minima, args.xml) == plan
    passed = (
        unchanged
        and len(outcomes) == 9
        and all(
            r["status"] == ("expected_rejection" if r["case"] == "ps3" else "passed")
            for r in outcomes
        )
    )
    storage.publish(
        args.output / "outcomes.json",
        storage.json_bytes(
            {
                "declaration": plan,
                "cases": outcomes,
                "complete": passed,
                "inputs_unchanged": unchanged,
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
        replace=True,
    )
    return 0 if passed else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retained", type=Path)
    parser.add_argument("--minima", type=Path)
    parser.add_argument("--xml", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--execute-declared", action="store_true")
    parser.add_argument("--resume-directory", type=Path)
    args = parser.parse_args()
    if args.resume_directory:
        return worker(args.resume_directory)
    if any(getattr(args, k) is None for k in ("retained", "minima", "xml", "output")):
        parser.error("--retained, --minima, --xml and --output required")
    plan = declaration(args.retained, args.minima, args.xml)
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
    return execute(args, plan)


if __name__ == "__main__":
    raise SystemExit(main())
