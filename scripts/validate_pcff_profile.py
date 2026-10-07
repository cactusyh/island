"""J8 declared aromatic vertical slice. Native science engines are unchanged."""

import argparse
import builtins
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np
from pcff_j5_reference import raw_inventory, read_source, resolve_raw
from validate_oplsaa_dynamics import compare_frames
from validate_pcff_expanded import independent_charges
from validate_pcff_singlepoint import (
    GROUPS,
    aa_reference_coefficients,
    check_inventory,
    data_sections,
    lammps_input,
    reference_inputs,
    run,
)
from validate_prepared_workflow import count, hashes

from island.charge_references.records import pack, unpack
from island.dynamics import run_dynamics_segment
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedBundleError,
    PreparedForceFieldSources,
    create_evaluator,
    load_prepared_forcefield,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff import PCFFOperationalSelection
from island.minimization import MinimizationOptions
from island.minimization.models import force_metrics
from island.workflows import (
    PreparedWorkflowConfig,
    prepared_workflow_status,
    read_prepared_workflow_frames,
    resume_prepared_workflow,
    start_prepared_bundle_workflow,
    storage,
)
from island.workflows.bundle import system_from
from island.workflows.chain import _load_segment
from island.workflows.prepared import _load_setup


def require(ok, message):
    if not ok:
        raise ValueError(message)


def group_values(evaluation):
    groups = {**GROUPS, "improper": ["angle-angle", "wilson_out_of_plane"]}
    return np.array(
        [evaluation.potential_energy]
        + [sum(evaluation.energy_components[t] for t in fs) for fs in groups.values()]
    )


def selection_check(system, model, source):
    labels = {
        i: "cp" if a.element == "C" else "hc" for i, a in system.topology.sites.items()
    }
    raw = read_source(source.read_bytes())
    assignment = unpack(model.assignment.json_text)

    def key(f, sites):
        s = tuple(sites)
        if f == "angle-angle":
            s = min(s, (s[3], s[1], s[2], s[0]))
        elif f == "wilson_out_of_plane":
            s = (min(s[0], s[2], s[3]), s[1], *sorted([s[0], s[2], s[3]])[1:])
        else:
            s = min(s, s[::-1])
        return f, s

    expected = Counter(key(f, s) for f, s in raw_inventory(system))
    actual = Counter(key(a["family"], a["sites"]) for a in assignment["assignments"])
    require(expected == actual, "Independent authoritative inventory mismatch")
    rows = []
    for a in assignment["assignments"]:
        q = resolve_raw(raw, a["family"], [labels[i] for i in a["sites"]])
        require(
            q["status"] == a["status"] == "assigned",
            "Incomplete independently selected source record",
        )
        require(
            q["selected_ids"] == sorted({r["record_id"] for r in a["selected"]}),
            "Independent selected-row mismatch",
        )
        np.testing.assert_allclose(
            q["values"], a["normalized_values"], atol=1e-12, rtol=0
        )
        rows.append(
            {
                "family": a["family"],
                "sites": a["sites"],
                "source_rows": q["selected_ids"],
                "independent_values": q["values"],
            }
        )
    charges, contributions = independent_charges(source.read_bytes(), system, labels)
    require(abs(sum(charges.values())) <= 1e-12, "Independent native charge residual")
    require(
        charges
        == {n["site"]: n["charge"] for n in unpack(model.json_text)["nonbonded"]},
        "Independent native per-site charge mismatch",
    )
    return (
        labels,
        charges,
        {
            "inventory": dict(Counter(f for f, s in expected.elements())),
            "rows": rows,
            "charge_contributions": contributions,
        },
    )


def reference(exe, work, xyz, ids, aa):
    work.mkdir()
    shutil.copyfile(work.parent / "reference.data", work / "reference.data")
    (work / "input.lmp").write_text(lammps_input(xyz, [0, 0, 1], [0, 0, 1], aa))
    run([str(exe), "-in", "input.lmp"], work, "lammps")
    line = [
        x
        for x in (work / "lammps.log").read_text().splitlines()
        if x.startswith("ENERGIES ")
    ][-1]
    energy = np.array(list(map(float, line.split()[1:]))) * 4.184
    text = (
        (work / "forces.dump").read_text().split("ITEM: ATOMS id x y z fx fy fz\n")[1]
    )
    forces = (
        np.array(
            [list(map(float, x.split()[4:7])) for x in text.splitlines() if x.strip()]
        )
        * 4.184
    )
    return energy, forces


def numerical(evaluator, system, xyz, work, name, exe, aa, tol):
    ids = sorted(system.topology.sites)
    coordinates = dict(zip(ids, xyz, strict=True))
    result = evaluator.evaluate_fresh(coordinates)
    e, f = reference(exe, work / name, xyz, ids, aa)
    actual = group_values(result)
    force = np.array([result.forces[i] for i in ids])
    np.testing.assert_allclose(
        actual, e, atol=tol["energy_atol_kj_mol"], rtol=tol["reference_rtol"]
    )
    np.testing.assert_allclose(
        force, f, atol=tol["force_atol_kj_mol_angstrom"], rtol=tol["reference_rtol"]
    )
    return {
        "configuration": name,
        "energy_component_max_abs": float(np.max(np.abs(actual - e))),
        "force_max_abs": float(np.max(np.abs(force - f))),
        "energy": result.potential_energy,
        "components": dict(result.energy_components),
    }


def tamper_gate(bundle, sources, output):
    before = hashes(bundle)
    checks = []
    output.mkdir()
    for kind in (
        "profile" if (bundle / "pcff-profile.json").exists() else "policy",
        "model_source",
    ):
        copy = output / kind
        shutil.copytree(bundle, copy)
        filename = "pcff-profile.json" if kind == "profile" else "model.json"
        data = unpack((copy / filename).read_text())
        if kind == "profile":
            data["source_sha256"] = (
                "3ad5a1be7334c646ed6cb813b769d0e89a1aa03fe695013e940626b7df922693"
            )
        elif kind == "policy":
            data["compatibility_profile"]["resolution_policy"]["name"] = "invented"
        else:
            data["source"]["sha256"] = (
                "3ad5a1be7334c646ed6cb813b769d0e89a1aa03fe695013e940626b7df922693"
            )
        (copy / filename).write_text(pack(data))
        env = storage.read_json(copy / "manifest.json")
        logical = "pcff_profile" if kind == "profile" else "model"
        env["payload"]["files"][logical]["sha256"] = storage.checksum(
            (copy / filename).read_bytes()
        )
        if kind == "profile":
            env["payload"]["prepared"]["pcff_operational_profile"] = data
        env["sha256"] = storage.checksum(storage.json_bytes(env["payload"]))
        (copy / "manifest.json").write_bytes(storage.json_bytes(env))
        try:
            load_prepared_forcefield(copy, sources=sources)
        except PreparedBundleError as error:
            require(
                "checksum" not in str(error).lower(), "Tamper test hit stale checksum"
            )
            checks.append({"case": kind, "semantic_rejection": str(error)})
        else:
            raise ValueError("Rechecksummed operational contradiction accepted")
    require(hashes(bundle) == before, "Tamper controls changed original bundle")
    return {"passed": len(checks) == 2, "checks": checks, "original_unchanged": True}


def child(args):
    sources = PreparedForceFieldSources(pcff_frc=args.source)
    importer = builtins.__import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if level == 0 and name.split(".")[0] in {
            "openmm",
            "rdkit",
            "scipy",
            "foyer",
            "parmed",
        }:
            raise AssertionError("Scientific import during reconstruction")
        return importer(name, globals, locals, fromlist, level)

    builtins.__import__ = blocked
    try:
        loaded = load_prepared_forcefield(
            args.output / "bundle-relocated", sources=sources
        )
        declaration = storage.read_json(args.output / "declaration.json")
        if declaration.get("operational_profile"):
            require(
                loaded.prepared.operational_profile is not None,
                "Lost operational authorization",
            )
        else:
            require(
                loaded.prepared.native_result.payload["compatibility_profile"][
                    "resolution_policy"
                ]["name"]
                == declaration["resolution_policy"],
                "Lost native resolver policy",
            )
    finally:
        builtins.__import__ = importer
    # Resume runs existing engines, while typing/preparation/minimization/initialization are forbidden.
    from unittest.mock import patch

    import island.workflows.prepared as workflow

    def forbidden(*a, **kw):
        raise AssertionError("Scientific setup during resume")

    with (
        patch.object(workflow, "minimize_geometry", forbidden),
        patch.object(workflow, "initialize_velocities", forbidden),
    ):
        status, contexts = count(
            lambda: resume_prepared_workflow(
                args.output / "workflow-relocated", sources=sources
            )
        )
    storage.publish(
        args.output / "child.json",
        storage.json_bytes(
            {
                "pid": os.getpid(),
                "status": status["status"],
                "contexts": contexts,
                "profile": loaded.prepared.operational_profile.identity
                if loaded.prepared.operational_profile
                else None,
                "prepared": loaded.prepared.identity,
                "native": loaded.prepared.native_result.identity,
                "offline_reconstruction": True,
            }
        ),
    )
    return 0 if status["status"] == "completed" else 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--declaration",
        type=Path,
        default=Path("docs/evidence/phase_4j8_declaration.json"),
    )
    p.add_argument("--child", action="store_true")
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    if a.child:
        return child(a)
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    outcome = {
        "vertical_slice_passed": False,
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        d = storage.read_json(a.declaration)
        storage.publish(root / "declaration.json", storage.json_bytes(d))
        require(
            storage.checksum(a.source.read_bytes()) == d["source"]["sha256"],
            "Wrong source hash",
        )
        for k in ["lammps_path", "converter_path"]:
            require(
                storage.checksum(Path(d["reference"][k]).read_bytes())
                == d["reference"][k + "_sha256"],
                "Wrong executable identity",
            )
        original = Path(d["case"]["system_path"])
        require(
            storage.checksum(original.read_bytes()) == d["case"]["system_sha256"],
            "Historical input changed",
        )
        system = system_from(storage.decode(storage.read_json(original)))
        before = system.to_dict()
        sel = (
            PCFFOperationalSelection(d["operational_profile"], d["source"]["sha256"])
            if d.get("operational_profile")
            else None
        )
        prepared = prepare_forcefield(
            system,
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    a.source,
                    (0, 0, 1),
                    (0, 0, 1),
                    typing_profile=d["typing_profile"],
                    source_profile=sel,
                    resolution_policy=d.get("resolution_policy"),
                ),
            ),
        )
        sources = PreparedForceFieldSources(pcff_frc=a.source)
        save_prepared_forcefield(system, prepared, root / "bundle", sources=sources)
        outcome["tamper"] = tamper_gate(
            root / "bundle", sources, root / "tamper-copies"
        )
        model = prepared.native_result
        ids = sorted(system.topology.sites)
        xyz = np.array([system.coordinates.get(i) for i in ids])
        check_selection = selection_check
        if d.get("reference_label_policy") == "declared_saturated_cho_s_v1":
            from pcff_j13_reference import selection_check as check_selection
        if d.get("reference_label_policy") == "declared_ref1_aliphatic_amines_v1":
            from pcff_j14_reference import selection_check as check_selection
        labels, charges, selected = check_selection(system, model, a.source)
        storage.publish(
            root / "independent-selection.json", storage.json_bytes(selected)
        )
        ref = root / "reference"
        ref.mkdir()
        _rid, neighbors = reference_inputs(system, labels, charges, xyz, ref)
        run(
            [
                str(Path(d["reference"]["converter_path"]).resolve()),
                "reference",
                "-class",
                "II",
                "-frc",
                str(a.source.resolve()),
                "-p",
                "3",
                "-nocenter",
            ],
            ref,
            "converter",
        )
        sections = data_sections(ref / "reference.data")
        check_inventory(sections, ids, neighbors, charges)
        if d.get("reference_label_policy") in (
            "declared_saturated_cho_s_v1",
            "declared_ref1_aliphatic_amines_v1",
        ):
            from pcff_j13_reference import converter_zero_check

            outcome["compiled_converter_bb13"] = converter_zero_check(sections)
        aa = aa_reference_coefficients(sections)
        if d.get("reference_aa_policy") == "raw_source_center_preserving_j14_v1":
            from pcff_j14_reference import aa_raw_reference_coefficients

            aa, evidence = aa_raw_reference_coefficients(sections, labels, a.source)
            storage.publish(
                ref / "aa-raw-source-correction.json", storage.json_bytes(evidence)
            )
        storage.publish(ref / "aa-overrides.json", storage.json_bytes(aa))
        evaluator = create_evaluator(system, prepared)
        pert = xyz + 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3)
        outcome["numerical"] = [
            numerical(
                evaluator,
                system,
                x,
                ref,
                n,
                Path(d["reference"]["lammps_path"]).resolve(),
                aa,
                d["tolerances"],
            )
            for n, x in [("initial", xyz), ("asymmetric", pert)]
        ]
        actual = evaluator.evaluate_fresh(dict(zip(ids, pert, strict=True)))
        force = np.array([actual.forces[i] for i in ids])
        fd = []
        with evaluator.open_session() as session:
            for h in d["tolerances"]["finite_difference_angstrom"]:
                approx = np.zeros_like(pert)
                for i in range(len(ids)):
                    for k in range(3):
                        plus = pert.copy()
                        minus = pert.copy()
                        plus[i, k] += h
                        minus[i, k] -= h
                        ep = session.evaluate(
                            dict(zip(ids, plus, strict=True))
                        ).potential_energy
                        em = session.evaluate(
                            dict(zip(ids, minus, strict=True))
                        ).potential_energy
                        approx[i, k] = -(ep - em) / (2 * h)
                fd.append(
                    {
                        "step_angstrom": h,
                        "maximum_force_error": float(np.max(np.abs(approx - force))),
                    }
                )
        require(
            fd[-1]["maximum_force_error"]
            <= d["tolerances"]["force_atol_kj_mol_angstrom"],
            "Finite-difference gate failed",
        )
        outcome["finite_difference"] = fd
        config = PreparedWorkflowConfig(
            minimization=MinimizationOptions(
                force_tolerance=0.1, max_iterations=5000, max_evaluations=10000
            ),
            total_steps=4,
            segment_steps=2,
            recording_interval=1,
            max_evaluations_per_segment=4,
            max_frames_per_segment=3,
        )
        status, contexts = count(
            lambda: start_prepared_bundle_workflow(
                root / "bundle", root / "workflow", config, sources=sources
            )
        )
        require(
            status["status"] == "paused" and status["accepted_step"] == 2,
            "Minimum/first segment gate failed",
        )
        optimized, restored, initialization, minimum = _load_setup(
            root / "workflow", status, sources, with_minimum=True
        )
        require(
            minimum.converged
            and minimum.final_evaluation_verified
            and force_metrics(minimum.final_evaluation)[0] <= 0.1,
            "Verified force convergence required",
        )
        outcome["minimum"] = {
            "iterations": minimum.iterations,
            "evaluations": minimum.evaluations,
            "initial_energy": minimum.initial_evaluation.potential_energy,
            "final_energy": minimum.final_evaluation.potential_energy,
            "fmax": force_metrics(minimum.final_evaluation)[0],
            "rms": force_metrics(minimum.final_evaluation)[1],
            "reason": minimum.termination_reason,
            "fresh_verified": minimum.final_evaluation_verified,
        }
        opt_eval = create_evaluator(optimized, restored)

        def uninterrupted():
            with opt_eval.open_session() as session:
                return run_dynamics_segment(
                    optimized,
                    session,
                    initialization.velocities,
                    replace(config.langevin(4), max_evaluations=6, max_frames=5),
                )

        whole, whole_contexts = count(uninterrupted)
        require(whole.completed, "Whole dynamics failed")
        storage.publish(
            root / "whole.json",
            storage.json_bytes(
                {"payload": whole.payload, "sha256": whole.content_checksum}
            ),
        )
        shutil.move(root / "bundle", root / "bundle-relocated")
        shutil.move(root / "workflow", root / "workflow-relocated")
        cmd = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--child",
            "--source",
            str(a.source.resolve()),
            "--output",
            str(root),
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=d["reference"]["child_seconds"],
        )
        storage.publish(root / "child.stdout", result.stdout.encode())
        storage.publish(root / "child.stderr", result.stderr.encode())
        require(result.returncode == 0, "Separate-process continuation failed")
        status = prepared_workflow_status(root / "workflow-relocated", sources=sources)
        frames = read_prepared_workflow_frames(
            root / "workflow-relocated", sources=sources
        )
        final = _load_segment(root / "workflow-relocated", status["segments"][-1])
        outcome["differences"] = compare_frames(whole.frames, frames)
        require(whole.payload["rng"] == final.payload["rng"], "Complete RNG mismatch")
        require(whole.payload["origin"] == final.payload["origin"], "Origin mismatch")
        require(
            resume_prepared_workflow(root / "workflow-relocated", sources=sources)
            == status,
            "Completed no-op resume mismatch",
        )
        outcome["numerical"] += [
            numerical(
                opt_eval,
                optimized,
                np.array([f.coordinates[i] for i in ids]),
                ref,
                "step" + str(f.step),
                Path(d["reference"]["lammps_path"]).resolve(),
                aa,
                d["tolerances"],
            )
            for f in frames
            if f.step in (0, 2, 4)
        ]
        require(
            system.to_dict() == before,
            "Caller system mutation",
        )
        require(
            storage.checksum(original.read_bytes()) == d["case"]["system_sha256"],
            "Historical input mutation",
        )
        outcome.update(
            vertical_slice_passed=True,
            profile=prepared.operational_profile.identity
            if prepared.operational_profile
            else None,
            prepared=prepared.identity,
            native=model.identity,
            start_contexts=contexts,
            whole_contexts=whole_contexts,
            whole_evaluations=whole.evaluations,
            split_evaluations=final.payload["counters"]["evaluations"],
            rng_equal=True,
            child=storage.read_json(root / "child.json"),
            child_command=cmd,
            frames=len(frames),
        )
    except Exception as e:  # noqa: BLE001 -- preserve all failed acceptance stages
        outcome["error"] = type(e).__name__ + ": " + str(e)
    finally:
        outcome["artifact_hashes"] = hashes(root)
        storage.publish(root / "report.json", storage.json_bytes(outcome))
    print(
        json.dumps(
            {k: v for k, v in outcome.items() if k != "artifact_hashes"}, indent=2
        )
    )
    return 0 if outcome["vertical_slice_passed"] and not a.require_full_source else 1


if __name__ == "__main__":
    raise SystemExit(main())
