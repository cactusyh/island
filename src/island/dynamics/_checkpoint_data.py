"""Strict data-only checkpoint codecs and structural validation (no authenticity claim)."""

import hashlib
import io
import json
import platform
import sys
import warnings
from contextlib import redirect_stdout
from copy import deepcopy
from dataclasses import fields
from functools import lru_cache
from math import isclose
from pathlib import Path

import numpy as np

from island.evaluation import EvaluationResult
from island.exceptions import InvalidDynamicsCheckpointError as Invalid
from island.minimization.integrity import number
from island.minimization.models import system_identity

from .integrity import checked_evaluation, same_model, verify_final
from .langevin_models import BAOAB, LangevinFrame, LangevinOptions
from .models import (
    INTEGRATOR,
    DynamicsFrame,
    DynamicsOptions,
    kinetic_energy,
)
from .thermal import RNG_ALGORITHM

SCHEMA = "island_dynamics_checkpoint_v1"
SEGMENT_SCHEMA = "island_dynamics_segment_v1"
UNITS = {
    "coordinates": "angstrom",
    "velocities": "angstrom/ps",
    "mass": "dalton",
    "energy": "kJ/mol",
    "time": "ps",
}


def require(test, message):
    if not test:
        raise Invalid(message)


def keys(value, expected):
    require(
        type(value) is dict and set(value) == set(expected),
        "Unexpected or missing record fields",
    )


def integer(value, minimum=0):
    return type(value) is int and value >= minimum


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def checksum(value):
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()


def _pairs(pairs):
    result = {}
    for k, v in pairs:
        require(k not in result, f"Duplicate JSON key: {k}")
        result[k] = v
    return result


def strict_load(text):
    def reject(value):
        raise Invalid(f"Nonfinite JSON constant: {value}")

    return json.loads(text, object_pairs_hook=_pairs, parse_constant=reject)


def rows(mapping):
    return [
        [
            site,
            list(mapping[site]) if isinstance(mapping[site], tuple) else mapping[site],
        ]
        for site in sorted(mapping)
    ]


def inventory(data):
    require(type(data) is list and bool(data), "Inventory must be a nonempty list")
    result = {}
    for row in data:
        require(
            type(row) is list and len(row) == 2 and type(row[0]) is int,
            "Invalid stable-ID row",
        )
        require(row[0] not in result, "Duplicate stable ID")
        result[row[0]] = row[1]
    require(
        list(result) == sorted(result), "Inventory must be in ascending stable-ID order"
    )
    return result


def evaluation_data(record):
    if record is None:
        return None
    result = {f.name: getattr(record, f.name) for f in fields(record)}
    result["forces"] = rows(record.forces)
    result["energy_components"] = dict(record.energy_components)
    result["settings"] = dict(record.settings)
    return result


def evaluation_from(data):
    if data is None:
        return None
    keys(data, [f.name for f in fields(EvaluationResult)])
    values = dict(data)
    values["forces"] = inventory(values["forces"])
    for f in fields(EvaluationResult):
        if not f.init:
            observed = values.pop(f.name)
            require(
                type(observed) is type(f.default) and observed == f.default,
                "Unsupported evaluation units/readiness",
            )
    return EvaluationResult(**values)


def frame_data(frame):
    return {
        f.name: (
            rows(getattr(frame, f.name))
            if f.name in ("coordinates", "velocities")
            else evaluation_data(frame.evaluation)
            if f.name == "evaluation"
            else getattr(frame, f.name)
        )
        for f in fields(frame)
    }


def frame_from(data, thermal):
    cls = LangevinFrame if thermal else DynamicsFrame
    keys(data, [f.name for f in fields(cls)])
    values = dict(data)
    for name in ("coordinates", "velocities"):
        values[name] = inventory(values[name])
    values["evaluation"] = evaluation_from(values["evaluation"])
    for f in fields(cls):
        if not f.init:
            observed = values.pop(f.name)
            require(
                type(observed) is type(f.default) and observed == f.default,
                "Unsupported frame units",
            )
    result = cls(**values)
    result.validate_integrity()
    return result


def physical(options):
    result = {
        "integrator": BAOAB if type(options) is LangevinOptions else INTEGRATOR,
        "timestep_fs": options.timestep_fs,
    }
    names = (
        ("temperature_kelvin", "friction_per_ps", "thermostat_seed")
        if type(options) is LangevinOptions
        else ("max_energy_deviation",)
    )
    result.update({name: getattr(options, name) for name in names})
    return result


def validate_physical(data):
    require(type(data) is dict, "Physical settings must be an object")
    thermal = data.get("integrator") == BAOAB
    expected = {"integrator", "timestep_fs"} | (
        {"temperature_kelvin", "friction_per_ps", "thermostat_seed"}
        if thermal
        else {"max_energy_deviation"}
    )
    keys(data, expected)
    require(data["integrator"] in (BAOAB, INTEGRATOR), "Unsupported integrator")
    (LangevinOptions if thermal else DynamicsOptions)(
        **{k: v for k, v in data.items() if k != "integrator"},
        steps=1,
        max_evaluations=2,
        max_frames=1,
    )
    return thermal


def validate_rng(state):
    keys(state, {"bit_generator", "state", "has_uint32", "uinteger"})
    require(state["bit_generator"] == "PCG64", "Only PCG64 state is supported")
    keys(state["state"], {"state", "inc"})
    for name in ("state", "inc"):
        value = state["state"][name]
        require(integer(value) and value < 2**128, "Invalid 128-bit PCG64 state")
    require(state["state"]["inc"] % 2 == 1, "PCG64 increment must be odd")
    require(
        type(state["has_uint32"]) is int and state["has_uint32"] in (0, 1),
        "Invalid PCG64 cache flag",
    )
    require(
        integer(state["uinteger"]) and state["uinteger"] < 2**32,
        "Invalid PCG64 cached integer",
    )
    bitgen = np.random.PCG64(0)
    bitgen.state = deepcopy(state)
    require(bitgen.state == state, "PCG64 state did not round trip exactly")


@lru_cache(maxsize=1)
def _environment_json():
    stream = io.StringIO()
    with warnings.catch_warnings(), redirect_stdout(stream):
        warnings.filterwarnings("ignore", message="Install `pyyaml` for better output")
        np.show_config()
    cpu = platform.processor() or platform.machine()
    try:
        labels = {
            "vendor_id",
            "cpu family",
            "model",
            "model name",
            "stepping",
            "microcode",
            "flags",
        }
        cpu = canonical(
            sorted(
                {
                    line.strip()
                    for line in Path("/proc/cpuinfo").read_text().splitlines()
                    if line.split(":", 1)[0].strip() in labels
                }
            )
        )
    except OSError:
        pass
    return canonical(
        {
            "host": platform.node() or "unknown",
            "cpu_sha256": hashlib.sha256(cpu.encode()).hexdigest(),
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "numpy_build_sha256": hashlib.sha256(
                stream.getvalue().encode()
            ).hexdigest(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor() or platform.machine(),
            "byteorder": sys.byteorder,
        }
    )


def environment():
    return strict_load(_environment_json())


def compatibility(system):
    # Only the established coordinate bookkeeping namespace is normalized.
    # Chemical, repeat, preparation, and parameter metadata remain signed here.
    copied = deepcopy(system)
    for k in (
        "coordinate_source",
        "coordinate_generation",
        "coordinate_history",
        "local_minimization",
        "dynamics",
    ):
        copied.metadata.pop(k, None)
    for k, v in (
        ("production_validated", False),
        ("simulation_readiness", "not_established"),
    ):
        if copied.metadata.get(k) == v:
            copied.metadata.pop(k, None)
    polymer = copied.metadata.get("polymer")
    if isinstance(polymer, dict):
        polymer.pop("coordinates", None)
        polymer.pop("coordinate_generation", None)
    return system_identity(copied)


def validate_payload(p, *, segment=False):
    base = {
        "schema",
        "units",
        "physical",
        "environment",
        "system_fingerprint",
        "origin",
        "lineage",
        "state",
        "masses",
        "rng",
        "counters",
        "max_abs_energy_deviation",
        "final_evaluation",
        "production_validated",
        "simulation_readiness",
    }
    if segment:
        base |= {"frames", "diagnostic"}
    keys(p, base)
    require(
        p["schema"] == (SEGMENT_SCHEMA if segment else SCHEMA),
        "Unsupported checkpoint schema",
    )
    require(p["units"] == UNITS, "Unsupported checkpoint units")
    require(
        p["production_validated"] is False
        and p["simulation_readiness"] == "not_established",
        "Unsupported readiness",
    )
    thermal = validate_physical(p["physical"])
    dt = p["physical"]["timestep_fs"] * 0.001
    env = p["environment"]
    keys(env, environment())
    require(
        all(type(v) is str and v for v in env.values()), "Invalid environment identity"
    )
    require(
        type(p["system_fingerprint"]) is str and len(p["system_fingerprint"]) == 64,
        "Invalid system identity",
    )
    masses = inventory(p["masses"])
    require(all(number(m) and m > 0 for m in masses.values()), "Invalid masses")
    frame = frame_from(p["state"], thermal)
    require(set(frame.coordinates) == set(masses), "Mass/state coverage mismatch")
    keys(
        p["origin"], {"trajectory_fingerprint", "initial_total_energy", "initial_state"}
    )
    require(
        type(p["origin"]["trajectory_fingerprint"]) is str
        and len(p["origin"]["trajectory_fingerprint"]) == 64,
        "Invalid trajectory identity",
    )
    initial = frame_from(p["origin"]["initial_state"], thermal)
    require(
        initial.step == 0
        and initial.time_ps == 0
        and set(initial.coordinates) == set(masses),
        "Invalid original state",
    )
    require(
        isclose(
            initial.kinetic_energy,
            kinetic_energy(masses, initial.velocities),
            rel_tol=1e-12,
            abs_tol=1e-10,
        ),
        "Invalid original kinetic energy",
    )
    require(
        p["origin"]["initial_total_energy"] == initial.total_energy,
        "Original energy reference changed",
    )
    require(
        p["origin"]["trajectory_fingerprint"] == trajectory_identity(p),
        "Original trajectory identity changed",
    )
    require(
        initial.energy_deviation == (0.0 if initial.evaluation is not None else None),
        "Invalid initial energy deviation",
    )
    if initial.evaluation is not None and frame.evaluation is not None:
        same_model(initial.evaluation, frame.evaluation)
    elif frame.evaluation is not None:
        require(False, "Missing initial potential identity")
    energy = p["origin"]["initial_total_energy"]
    require(energy is None or number(energy), "Invalid original total energy")
    require(
        type(p["lineage"]) is list and bool(p["lineage"]), "Missing segment lineage"
    )
    total_calls = 0
    end = 0
    for index, row in enumerate(p["lineage"]):
        keys(
            row,
            {
                "index",
                "start_step",
                "end_step",
                "evaluations",
                "parent_checksum",
                "options",
                "termination_reason",
            },
        )
        require(
            type(row["index"]) is int and row["index"] == index, "Invalid segment index"
        )
        require(
            type(row["start_step"]) is int
            and row["start_step"] == end
            and integer(row["end_step"])
            and row["end_step"] >= end,
            "Noncontiguous lineage",
        )
        end = row["end_step"]
        total_calls += row["evaluations"] if integer(row["evaluations"], 1) else 0
        keys(
            row["options"],
            {"steps", "max_evaluations", "max_frames", "recording_interval"},
        )
        require(
            all(integer(v, 1) for v in row["options"].values()),
            "Invalid segment budgets",
        )
        count = end - row["start_step"]
        require(
            count <= row["options"]["steps"]
            and integer(row["evaluations"], 1)
            and row["evaluations"] <= row["options"]["max_evaluations"],
            "Invalid lineage counts",
        )
        require(
            row["parent_checksum"] is None
            if index == 0
            else type(row["parent_checksum"]) is str
            and len(row["parent_checksum"]) == 64,
            "Invalid parent checksum",
        )
        if index < len(p["lineage"]) - 1 or not segment:
            require(
                row["termination_reason"]
                in ("completed", "maximum_evaluations", "maximum_frames")
                and row["evaluations"] == count + 2,
                "Ineligible checkpoint lineage",
            )
            if row["termination_reason"] == "completed":
                require(
                    count == row["options"]["steps"], "Incomplete completed segment"
                )
            if row["termination_reason"] == "maximum_evaluations":
                require(
                    row["evaluations"] == row["options"]["max_evaluations"]
                    and count < row["options"]["steps"],
                    "Unexhausted evaluation budget",
                )
            if row["termination_reason"] == "maximum_frames":
                require(
                    1
                    + (count + 1) // row["options"]["recording_interval"]
                    + int((count + 1) % row["options"]["recording_interval"] != 0)
                    > row["options"]["max_frames"]
                    and count < row["options"]["steps"],
                    "Unexhausted frame budget",
                )
    require(
        frame.step == end and frame.time_ps == end * dt, "Absolute time/step mismatch"
    )
    keys(p["counters"], {"evaluations", "random_steps", "normal_draws"})
    require(
        all(integer(v) for v in p["counters"].values())
        and p["counters"]["evaluations"] == total_calls,
        "Cumulative call count mismatch",
    )
    if thermal:
        keys(p["rng"], {"algorithm", "numpy_version", "state"})
        require(
            p["rng"]["algorithm"] == RNG_ALGORITHM
            and p["rng"]["numpy_version"] == env["numpy"],
            "RNG environment mismatch",
        )
        validate_rng(p["rng"]["state"])
        require(
            p["counters"]["normal_draws"]
            == 3 * len(masses) * p["counters"]["random_steps"],
            "Normal draw count mismatch",
        )
        require(
            p["counters"]["random_steps"]
            == end
            + (int(p["diagnostic"]["failure_stage"] == "trial") if segment else 0),
            "Random step count mismatch",
        )
    else:
        require(
            p["rng"] is None
            and p["counters"]["normal_draws"] == p["counters"]["random_steps"] == 0,
            "Unexpected NVE RNG state",
        )
    frames = [frame] if not segment else [frame_from(f, thermal) for f in p["frames"]]
    require(bool(frames), "Missing segment frames")
    for f in frames:
        require(
            set(f.coordinates) == set(masses) and f.time_ps == f.step * dt,
            "Invalid frame inventory/time",
        )
        require(
            isclose(
                f.kinetic_energy,
                kinetic_energy(masses, f.velocities),
                rel_tol=1e-12,
                abs_tol=1e-10,
            ),
            "Invalid kinetic energy",
        )
        if f.evaluation is not None:
            require(
                energy is not None
                and isclose(
                    f.energy_deviation,
                    f.total_energy - energy,
                    rel_tol=1e-12,
                    abs_tol=1e-10,
                ),
                "Original energy reference mismatch",
            )
            require(
                number(p["max_abs_energy_deviation"], nonnegative=True)
                and abs(f.energy_deviation) <= p["max_abs_energy_deviation"],
                "Invalid cumulative deviation",
            )
            if not thermal:
                require(
                    p["max_abs_energy_deviation"]
                    <= p["physical"]["max_energy_deviation"],
                    "NVE original energy guard violated",
                )
            if frame.evaluation is not None:
                same_model(f.evaluation, frame.evaluation)
    final = evaluation_from(p["final_evaluation"])
    if final is not None:
        checked_evaluation(final, frame.coordinates)
        require(frame.evaluation is not None, "No accepted evaluation")
        verify_final(frame.evaluation, final)
    if not segment:
        require(
            final is not None and frame.evaluation is not None,
            "Checkpoint requires verified accepted boundary",
        )
    else:
        validate_segment(p, frames, frame, final)


def validate_segment(p, frames, frame, final):
    row = p["lineage"][-1]
    opts = row["options"]
    start = row["start_step"]
    count = frame.step - start
    d = p["diagnostic"]
    keys(
        d,
        {
            "completed",
            "startup_verified",
            "final_verified",
            "final_attempted",
            "failed_trial_evaluations",
            "attempted_step",
            "failure_stage",
            "message",
            "final_error",
        },
    )
    require(
        all(
            type(d[k]) is bool
            for k in (
                "completed",
                "startup_verified",
                "final_verified",
                "final_attempted",
            )
        ),
        "Invalid verification flags",
    )
    require(
        type(d["message"]) is str
        and bool(d["message"])
        and (d["final_error"] is None or type(d["final_error"]) is str),
        "Invalid diagnostics",
    )
    require(
        type(d["failed_trial_evaluations"]) is int
        and d["failed_trial_evaluations"] in (0, 1),
        "Invalid failed calls",
    )
    require(
        row["evaluations"]
        == 1 + count + d["failed_trial_evaluations"] + int(d["final_attempted"]),
        "Local call accounting mismatch",
    )
    expected = sorted(
        set(range(start, frame.step + 1, opts["recording_interval"])) | {frame.step}
    )
    require(
        [f.step for f in frames] == expected
        and len(frames) <= opts["max_frames"]
        and frame_data(frames[-1]) == p["state"],
        "Invalid segment recording schedule",
    )
    reason = row["termination_reason"]
    stages = {
        "completed": None,
        "maximum_evaluations": "budget",
        "maximum_frames": "budget",
        "evaluation_failed": "trial",
        "invalid_geometry": "trial",
        "stereochemistry_changed": "trial",
        "energy_guard_exceeded": "trial",
        "startup_verification_failed": "startup",
        "final_evaluation_failed": "final",
    }
    require(
        reason in stages and d["failure_stage"] == stages[reason],
        "Invalid failure stage/reason",
    )
    require(
        d["completed"] == (reason == "completed")
        and d["final_verified"] == (final is not None),
        "Inconsistent completion/verification",
    )
    require(
        not d["final_verified"] or d["final_attempted"] and d["final_error"] is None,
        "Invalid final verification",
    )
    require(
        (frames[0].evaluation is not None)
        if d["startup_verified"]
        else reason == "startup_verification_failed",
        "Invalid startup verification",
    )
    if reason == "completed":
        require(
            count == opts["steps"]
            and d["final_verified"]
            and d["startup_verified"]
            and d["attempted_step"] is None,
            "Incomplete successful segment",
        )
    else:
        attempt = (
            start
            if reason == "startup_verification_failed"
            else frame.step
            if reason == "final_evaluation_failed"
            else frame.step + 1
        )
        require(
            type(d["attempted_step"]) is int and d["attempted_step"] == attempt,
            "Invalid attempted step",
        )
    if reason == "startup_verification_failed":
        require(
            count == 0
            and row["evaluations"] == 1
            and final is None
            and not d["final_attempted"],
            "Invalid startup failure",
        )
    if reason in ("maximum_evaluations", "maximum_frames"):
        require(
            count < opts["steps"] and not d["failed_trial_evaluations"],
            "Invalid budget stop",
        )
        if reason == "maximum_evaluations":
            require(
                row["evaluations"] == opts["max_evaluations"],
                "Unexhausted evaluator budget",
            )
        else:
            require(
                1
                + (count + 1) // opts["recording_interval"]
                + int((count + 1) % opts["recording_interval"] != 0)
                > opts["max_frames"],
                "Unexhausted frame budget",
            )
    if d["startup_verified"] and not d["final_attempted"]:
        require(
            row["evaluations"] == opts["max_evaluations"] == 1,
            "Omitted reserved final check",
        )
    if d["failure_stage"] != "trial":
        require(d["failed_trial_evaluations"] == 0, "Misassigned failed call")
    if d["startup_verified"] and not d["final_verified"]:
        require(
            type(d["final_error"]) is str and bool(d["final_error"]),
            "Missing final verification diagnostic",
        )
    if d["failure_stage"] not in ("startup", "final") and reason != "completed":
        require(count < opts["steps"], "Failure after all requested trials")
    if reason == "final_evaluation_failed":
        require(
            count == opts["steps"] and not d["final_verified"], "Invalid final failure"
        )


def trajectory_identity(p):
    return checksum(
        {
            key: p[key]
            for key in ("physical", "system_fingerprint", "environment", "masses")
        }
        | {"initial_state": p["origin"]["initial_state"]}
    )
