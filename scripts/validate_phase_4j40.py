"""Finite J40 benchmark; convergence and cross-engine failures remain nonzero."""

import argparse
import json
import shutil
import subprocess
import sys
import traceback
from hashlib import sha256
from pathlib import Path

from island.analysis.equilibration import diagnose
from island.analysis.equilibration_run import run_heated_workflow
from island.analysis.npt import FLAGS, NPTConfiguration, _signed, _write_json
from island.analysis.npt_continuation import (
    compare_continuations,
    relocate_continuation,
)
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges
from scripts.phase_4j32_assignment import independent_check, prepare
from scripts.phase_4j38_c3_reference import audit
from scripts.phase_4j38_contract import require
from scripts.phase_4j38_fixture import build
from scripts.validate_phase_4j38 import ANALYSIS
from scripts.validate_phase_4j38_c3 import command

BASE = "2610606f8ee15b6ea64df9790ae8ef126c922b7b"
CASES = {
    "gaff2-pe": ("pe", "RESP/Psi4"),
    "opls-resp-aromatic": ("aromatic", "RESP/Psi4"),
    "opls-peo": ("peo", "native_oplsaa"),
    "opls-network": ("network", "native_oplsaa"),
}


def engine_audit(executable, output):
    executable = executable.resolve()
    report = {
        "schema": "island.j40-engine.v1",
        "executable": str(executable),
        "status": "asset_unavailable",
        "kspace": False,
    }
    if executable.is_file():
        report["sha256"] = sha256(executable.read_bytes()).hexdigest()
        proc = subprocess.run(
            [str(executable), "-h"], capture_output=True, text=True, check=False
        )
        (output / "lammps-help.txt").write_text(proc.stdout + proc.stderr)
        report.update(
            help_sha256=sha256((proc.stdout + proc.stderr).encode()).hexdigest(),
            help_exit_code=proc.returncode,
            kspace="KSPACE" in proc.stdout and "ewald" in proc.stdout,
            version_lines=proc.stdout.splitlines()[:4],
        )
        report["status"] = (
            "available" if report["kspace"] and proc.returncode == 0 else "blocked"
        )
        source = Path("../island-validation/phase4h4-upstream/lammps")
        commit = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "36f7c81ced"], text=True
        ).strip()
        files = (
            "src/KSPACE/ewald.cpp",
            "src/KSPACE/pppm.cpp",
            "src/KSPACE/pair_lj_cut_coul_long.cpp",
            "src/version.h",
        )
        report["declared_source_commit"] = commit
        report["source_sha256"] = {
            name: sha256(
                subprocess.check_output(
                    ["git", "-C", str(source), "show", commit + ":" + name]
                )
            ).hexdigest()
            for name in files
        }
        report["build_provenance"] = (
            "Retained J22 KSPACE executable; Git version from actual help, source objects independently hashed; not rebuilt in J40"
        )
    return _signed(report)


def backend_for(key, fixture, output, raw):
    kind, origin = CASES[key]
    if origin == "native_oplsaa":
        backend = prepare(fixture, "oplsaa", output)
        _write_json(output / "independent-assignment.json", independent_check(backend))
        return backend
    require(
        json.loads((raw / "anchors.json").read_text())
        == json.loads(Path("scripts/phase_4j38_c3_anchors.json").read_text()),
        "Unreviewed raw anchors",
    )
    runtime = json.loads(Path("docs/evidence/phase_4j36/environment.json").read_text())
    review = audit(
        raw,
        kind,
        json.loads((raw / (kind + "-charges.json")).read_text()),
        runtime,
        json.loads((raw / "anchors.json").read_text()),
    )
    _write_json(output / "raw-review.json", review)
    charges = FiniteRESPCharges.load(
        raw / (kind + "-charges.json"),
        fixture.packed.graph,
        expected_identity=review["charge_identity"],
    )
    charges.save(output / "charges.json")
    if key == "gaff2-pe":
        from scripts.validate_phase_4j36_c2 import prepare as gaff

        amber = json.loads(Path("docs/evidence/phase_4j34/runtime.json").read_text())
        backend = gaff(fixture, output, amber, charges)
        _write_json(output / "independent-assignment.json", independent_check(backend))
        return backend
    from island.forcefields.oplsaa.resp import prepare as resp
    from island.forcefields.oplsaa.source import load_oplsaa_source
    from scripts.phase_4j32_real_fixture import XML

    return resp(
        fixture.packed.system, fixture.packed.graph, load_oplsaa_source(XML), charges
    )


def worker(key, output, raw, engine):
    output.mkdir(parents=True, exist_ok=False)
    row = dict(
        schema="island.j40-cell.v1",
        case=key,
        charge_origin=CASES[key][1],
        software_execution="not_executed",
        equilibration_status="not_executed",
        scientific_validation_status="not_established",
        acceptance_met=False,
        **FLAGS,
    )
    stage = "graph_and_charge_source"
    try:
        fixture, evidence = build(CASES[key][0])
        _write_json(output / "fixture.json", evidence)
        backend = backend_for(key, fixture, output, raw)
        _write_json(output / "backend.json", backend.payload)
        replicas = []
        for replica in range(2):
            root = output / ("replica-" + str(replica))
            root.mkdir()
            shutil.copyfile(output / "fixture.json", root / "fixture.json")
            config = NPTConfiguration(
                family=backend.payload["config"]["family"],
                source_identity=backend.payload["config"]["source"]["sha256"],
                final_graph_identity=fixture.packed.graph.identity,
                packing_identity=fixture.packed.plan.identity,
                box_lengths_nm=tuple(x / 10 for x in fixture.packed.system.box.lengths),
                seed=4040 + replica * 10,
                velocity_seed=4041 + replica * 10,
                barostat_seed=4042 + replica * 10,
                nvt_steps=448,
                npt_steps=1024,
                checkpoint_interval=4,
                pressure_interval=25,
                cutoff_nm=0.5,
            )
            stage = "replica_" + str(replica) + "_uninterrupted"
            workflow = run_heated_workflow(
                backend, config, root / "uninterrupted", continuation=True
            )
            stage = "split_relocation_child"
            run_heated_workflow(backend, config, root / "split", stop_npt_step=512)
            restart = relocate_continuation(root / "split", root / "relocated")
            argv = [
                sys.executable,
                "-m",
                "scripts.validate_phase_4j33_c1",
                "--resume",
                str(root / "relocated"),
                "--destination",
                str(root / "resumed"),
                "--restart",
                restart["identity"],
                "--parent",
                restart["parent_workflow_identity"],
                "--config",
                config.identity,
            ]
            require(
                command(argv, root / "resume") == 0, "Fresh-process continuation failed"
            )
            comparison = compare_continuations(
                root / "uninterrupted", root / "resumed", root / "relocated"
            )
            _write_json(root / "comparison.json", comparison)
            require(
                comparison["passed"],
                "Uninterrupted/split histories or checkpoint/RNG differ",
            )
            stage = "mdanalysis"
            argv = [
                str(ANALYSIS),
                "-m",
                "scripts.phase_4j38_observables",
                "--case",
                str(root),
                "--output",
                str(root / "observables.json"),
            ]
            require(command(argv, root / "observables") == 0, "MDAnalysis failed")
            replicas.append(workflow)
        stage = "equilibration_diagnostics"
        diagnostics = diagnose(replicas)
        _write_json(output / "equilibration.json", diagnostics)
        row.update(
            software_execution="passed",
            equilibration_status=diagnostics["equilibration_status"],
            graph_identity=backend.payload["final_graph_identity"],
            packing_identity=fixture.packed.plan.identity,
            backend_identity=backend.payload["identity"],
            source_identity=backend.payload["config"]["source"]["sha256"],
            charge_identity=backend.payload["charge_identity"],
            parameter_completeness="passed",
            diagnostic_identity=diagnostics["identity"],
            replica_count=len(replicas),
        )
        stage = "cross_engine"
        if engine["status"] == "available":
            import openmm as mm
            from openmm import unit

            from scripts.phase_4j40_cross_engine import evaluate

            reference = evaluate(
                backend, mm, unit, output / "cross-engine", Path(engine["executable"])
            )
            reference.update(
                executable_identity=engine["identity"],
                electrostatics="ewald",
                comparison_geometry="original_exact_packed_graph",
                precision=backend.payload["config"]["pme_tolerance"],
            )
        else:
            reference = {
                "passed": False,
                "status": "blocked",
                "reason": "Declared KSPACE executable unavailable",
                "executable_identity": engine["identity"],
            }
        _write_json(output / "cross-engine.json", reference)
        row["cross_engine_status"] = (
            "passed"
            if reference["passed"]
            else "failed"
            if engine["status"] == "available"
            else "blocked"
        )
        row["acceptance_met"] = diagnostics["acceptance_met"] and reference["passed"]
    except Exception as error:  # noqa: BLE001 -- preserve actual failure and partial artifacts
        (output / "failure.log").write_text(traceback.format_exc())
        row["diagnostic"] = {
            "stage": stage,
            "exception": type(error).__name__,
            "message": str(error),
        }
        row["software_execution"] = "failed"
    _write_json(output / "result.json", _signed(row))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    parser.add_argument("--worker", choices=CASES)
    args = parser.parse_args()
    if args.worker:
        engine = engine_audit(args.lammps, args.output.parent)
        return int(
            not worker(args.worker, args.output, args.raw, engine)["acceptance_met"]
        )
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(args.raw, args.output / "raw")
    engine = engine_audit(args.lammps, args.output)
    _write_json(args.output / "engine.json", engine)
    controls = args.output / "controls"
    controls.mkdir()
    from scripts.validate_phase_4j32 import worker as pcff_worker

    pcff = pcff_worker(controls / "pcff", "pcff", "network")
    require(
        pcff["acceptance_met"] and pcff["dynamics_executed"] is False,
        "PCFF rejection contract changed",
    )
    from scripts.phase_4j32_real_fixture import XML
    from scripts.phase_4j38_c1_audit import aromatic
    from scripts.phase_4j38_c2 import aromatic_decision

    _, evidence = build("aromatic")
    _write_json(
        controls / "native-aromatic.json",
        aromatic_decision(aromatic(evidence, XML), evidence),
    )
    limitation = Path("docs/evidence/phase_4j34_c2/receipt.json")
    _write_json(
        controls / "native-am1bcc.json",
        {
            "status": "blocked",
            "dynamics_executed": False,
            "diagnostic": "Pinned AmberTools 24.8 has no supported high-precision native AM1-BCC export; retained neutral fixture residual -0.003 e",
            "historical_evidence": str(limitation),
            "sha256": sha256(limitation.read_bytes()).hexdigest(),
            "reexecuted_quantum_generation": False,
        },
    )
    rows = {}
    from island.forcefields.ambertools.runtime import declared_environment

    amber = json.loads(Path("docs/evidence/phase_4j34/runtime.json").read_text())
    for key in CASES:
        argv = [
            sys.executable,
            "-m",
            "scripts.validate_phase_4j40",
            "--output",
            str(args.output / key),
            "--raw",
            str(args.output / "raw"),
            "--lammps",
            str(args.lammps.resolve()),
            "--worker",
            key,
        ]
        command(
            argv,
            args.output / (key + "-command"),
            env=declared_environment(amber) if key == "gaff2-pe" else None,
        )
        rows[key] = json.loads((args.output / key / "result.json").read_text())
    _write_json(
        args.output / "matrix.json",
        _signed(
            dict(
                schema="island.j40-matrix.v1",
                base_commit=BASE,
                cases=rows,
                engine_identity=engine["identity"],
                unsupported_origins={
                    "native_aromatic_oplsaa": "source_incompatible",
                    "native_am1bcc": "blocked_precision",
                    "unsupported_pcff": "typing_rejection",
                },
                acceptance_met=all(r["acceptance_met"] for r in rows.values()),
                **FLAGS,
            )
        ),
    )
    return int(not all(r["acceptance_met"] for r in rows.values()))


if __name__ == "__main__":
    raise SystemExit(main())
