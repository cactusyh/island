"""Strict J32 native acceptance CLI. A diagnostic failure always exits nonzero."""

import argparse
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

from island.analysis.network import _signed
from island.periodic_md import (
    PeriodicMDWorkflowConfig,
    compare_periodic_md_workflows,
    load_periodic_md_workflow,
    save_periodic_md_workflow,
    start_periodic_md_workflow,
)
from scripts.phase_4j31_contract import PCFF_SHA, typing_evidence
from scripts.phase_4j32_assignment import checksum, independent_check, prepare
from scripts.phase_4j32_contract import CASES, validate_case, validate_matrix
from scripts.phase_4j32_real_fixture import (
    AMBER,
    FOYER,
    FRC,
    ROOT,
    build_fixture,
    require,
)
from scripts.validate_phase_4j31_network import _expected_pcff_exception

FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
ANALYSIS_PYTHON = ROOT / "island-validation/phase4j32-analysis-env/bin/python"


def write(path, data):
    path.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")


def run_logged(command, log):
    process = subprocess.run(
        list(map(str, command)), capture_output=True, text=True, check=False
    )
    log.write_text(process.stdout + process.stderr)
    require(process.returncode == 0, f"Subprocess failed ({process.returncode}): {log}")
    return process.stdout


def pcff_contract(fixture, evidence, error):
    ids = sorted(fixture.packed.system.topology.sites)
    return (
        checksum(FRC) == PCFF_SHA
        and evidence["source_verified"] is True
        and evidence["source_sha256"] == PCFF_SHA
        and evidence["final_graph_identity"] == fixture.packed.graph.identity
        and evidence["profile"] == "island_pcff_acyclic_cho_v1"
        and evidence["complete"] is False
        and evidence["diagnostics"] == [{"sites": ids, "reason": "cyclic_component"}]
        and _expected_pcff_exception(error)
    )


def worker(output, family, kind):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    stage = "fixture"
    result = {
        "family": family,
        "fixture": kind,
        "status": "failed",
        "acceptance_met": False,
        **FLAGS,
    }
    try:
        fixture = build_fixture(kind)
        write(output / "fixture.json", fixture.evidence())
        result.update(
            final_graph_identity=fixture.packed.graph.identity,
            packing_identity=fixture.packed.plan.identity,
            setup=fixture.evidence()["setup"],
        )
        stage = "preparation"
        if family == "pcff":
            shutil.copyfile(FRC, output / "pcff.frc")
            evidence = typing_evidence(fixture, FRC)
            try:
                prepare(fixture, family, output)
            except Exception as error:  # noqa: BLE001 -- persist diagnostics and fail acceptance
                require(
                    pcff_contract(fixture, evidence, error),
                    "Unexpected PCFF failure; frozen typing outcome not met: "
                    + str(error),
                )
                (output / "expected-typing-rejection.log").write_text(
                    traceback.format_exc()
                )
                result.update(
                    status="blocked",
                    acceptance_met=True,
                    dynamics_executed=False,
                    diagnostic={
                        "classification": "expected_typing_rejection",
                        "stage": "typing",
                        "exception_type": type(error).__name__,
                        "message": str(error),
                        "origin": "island.forcefields.pcff.automatic.bridge_data",
                        "typing_evidence": evidence,
                    },
                )
                validate_case(f"{family}/{kind}", _signed(result))
                write(output / "result.json", _signed(result))
                return result
            raise ValueError("Unexpected PCFF preparation success")
        backend = prepare(fixture, family, output)
        stage = "independent_verification"
        assignment = independent_check(backend)
        assignment["environment"] = {
            name: importlib.metadata.version(name)
            for name in ("numpy", "rdkit", "openmm", "parmed")
        }
        assignment["python_executable_sha256"] = checksum(
            Path(sys.executable).resolve()
        )
        if family == "gaff2":
            assignment["source_executables"] = {
                str(p.relative_to(AMBER)): checksum(p)
                for p in [
                    *(
                        AMBER / "bin" / n
                        for n in ("antechamber", "parmchk2", "tleap", "sqm")
                    ),
                    AMBER / "dat/leap/parm/gaff2.dat",
                    AMBER / "dat/antechamber/BCCPARM.DAT",
                    AMBER / "dat/antechamber/ATOMTYPE_GFF.DEF",
                    AMBER / "dat/antechamber/ATOMTYPE_BCC.DEF",
                    AMBER / "bin/am1bcc",
                    AMBER / "bin/atomtype",
                    AMBER / "bin/wrapped_progs/antechamber",
                ]
            }
            for relative, digest in assignment["source_executables"].items():
                if relative.startswith("dat/"):
                    archived = output / "source-assets" / "amber-parameters" / relative
                    archived.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(AMBER / relative, archived)
                    require(
                        checksum(archived) == digest, "Source changed while archiving"
                    )
            assignment["amber_preparation_record"] = dict(
                backend._prepared.native_result.record
            )
        write(output / "assignment.json", assignment)
        config = PeriodicMDWorkflowConfig(
            backend.config.payload["family"],
            mode="nvt",
            steps=1000,
            timestep_fs=0.5,
            minimize=True,
            minimization_tolerance_kj_mol_nm=100.0,
            minimization_max_iterations=10000,
            seed=2032,
            velocity_seed=2033,
        )
        stage = "execution"
        start_periodic_md_workflow(backend, config, output / "workflow")
        result["dynamics_executed"] = True
        start_periodic_md_workflow(backend, config, output / "split", steps=500)
        save_periodic_md_workflow(output / "split", output / "continuation")
        identity = load_periodic_md_workflow(output / "continuation").identity
        run_logged(
            [
                sys.executable,
                "-c",
                "from island.periodic_md import resume_periodic_md_workflow;import sys;resume_periodic_md_workflow(sys.argv[1],expected_identity=sys.argv[2])",
                output / "continuation",
                identity,
            ],
            output / "continuation.log",
        )
        comparison = compare_periodic_md_workflows(
            output / "workflow", output / "continuation"
        )
        require(
            comparison["passed"],
            "Split/uninterrupted numerical or RNG comparison failed",
        )
        workflow = load_periodic_md_workflow(output / "workflow").payload
        execution = {
            "family": config.family,
            "comparison": comparison,
            "config": config.payload,
            "engine": workflow["manifest"]["engine"],
            "initial_box_vectors_nm": workflow["manifest"]["initial_box_vectors_nm"],
            "final_box_vectors_nm": workflow["frames"][-1]["box_vectors_nm"],
            "final_checkpoint_sha256": workflow["frames"][-1]["checkpoint_sha256"],
            "rng_evidence": "full binary OpenMM checkpoint bytes at every step and exact split comparison",
            "npt": {
                "status": "unexecuted",
                "reason": "Small finite oligomer acceptance run; density relaxation is not interpreted from this short sample",
            },
            "equilibrium_status": "not_established",
            "thermodynamic_origin": "J23 frame kinetic/potential energies and box; no inherited measurements",
            "short_nvt_duration_ps": 0.5,
            **FLAGS,
        }
        write(output / "execution.json", execution)
        # Seal the scientific bundle before publication; logs and results stay
        # outside it, so diagnostics do not invalidate its scientific identity.
        bundle = output / "bundle"
        bundle.mkdir()
        for name in (
            "fixture.json",
            "assignment.json",
            "execution.json",
            "workflow",
            "continuation",
            "source-assets",
        ):
            shutil.move(str(output / name), bundle / name)
        stage = "analysis"
        raw = run_logged(
            [
                ANALYSIS_PYTHON,
                "-c",
                'from island.analysis.condensed import save_bundle_record;import sys;print(save_bundle_record(sys.argv[1])["identity"])',
                bundle,
            ],
            output / "analysis.log",
        )
        record_identity = raw.strip()
        stage = "relocation"
        with tempfile.TemporaryDirectory(prefix="island-j32-relocated-") as temporary:
            relocated = Path(temporary) / "bundle"
            shutil.copytree(bundle, relocated)
            code = 'from island.analysis.condensed import load_bundle;import sys;print(load_bundle(sys.argv[1],expected_identity=sys.argv[2])["identity"])'
            # One process loads both locations, then another fresh interpreter
            # reconstructs solely from the moved bundle and packaged code.
            same = 'from island.analysis.condensed import load_bundle;import sys;a=load_bundle(sys.argv[1],expected_identity=sys.argv[3]);b=load_bundle(sys.argv[2],expected_identity=sys.argv[3]);assert a==b;print(a["identity"])'
            require(
                run_logged(
                    [ANALYSIS_PYTHON, "-c", same, bundle, relocated, record_identity],
                    output / "same-process.log",
                ).strip()
                == record_identity,
                "Relocation identity mismatch",
            )
            require(
                run_logged(
                    [ANALYSIS_PYTHON, "-c", code, relocated, record_identity],
                    output / "fresh-process.log",
                ).strip()
                == record_identity,
                "Fresh reconstruction identity mismatch",
            )
        result.update(
            status="complete",
            acceptance_met=True,
            analysis_identity=record_identity,
            same_process_relocation="PASS",
            fresh_process_reconstruction="PASS",
            comparison=comparison,
            diagnostic={"classification": "completed", "stage": "complete"},
        )
        validate_case(f"{family}/{kind}", _signed(result))
    except Exception as error:  # noqa: BLE001 -- persist diagnostics and fail acceptance
        (output / "failure.log").write_text(traceback.format_exc())
        result.update(
            status="failed",
            acceptance_met=False,
            reason=f"{type(error).__name__}: {error}",
            diagnostic={"classification": "unexpected_failure", "stage": stage},
        )
    write(output / "result.json", _signed(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--family", choices=("oplsaa", "gaff2", "pcff"))
    parser.add_argument("--kind", choices=("linear", "network"), default="network")
    args = parser.parse_args()
    os.environ["PATH"] = str(AMBER / "bin") + os.pathsep + os.environ.get("PATH", "")
    if args.family:
        require(
            f"{args.family}/{args.kind}" in CASES, "Case not in acceptance contract"
        )
        return int(not worker(args.output, args.family, args.kind)["acceptance_met"])
    args.output.mkdir(parents=True, exist_ok=False)
    report = {"cases": {}, "acceptance_met": False, **FLAGS}
    try:
        for name in CASES:
            family, kind = name.split("/")
            directory = args.output / name.replace("/", "-")
            executable = FOYER if family == "oplsaa" else sys.executable
            process = subprocess.run(
                [
                    str(executable),
                    "-m",
                    "scripts.validate_phase_4j32",
                    "--output",
                    str(directory),
                    "--family",
                    family,
                    "--kind",
                    kind,
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            (args.output / (name.replace("/", "-") + ".log")).write_text(
                process.stdout + process.stderr
            )
            row = (
                json.loads((directory / "result.json").read_text())
                if (directory / "result.json").exists()
                else {
                    "status": "failed",
                    "acceptance_met": False,
                    "reason": "No diagnostic record",
                }
            )
            report["cases"][name] = {"exit_code": process.returncode, "result": row}
            write(args.output / "matrix.json", _signed(report))
            print(name, row["status"], flush=True)
    finally:
        report["acceptance_met"] = set(report["cases"]) == set(CASES) and all(
            r["exit_code"] == 0 and r["result"]["acceptance_met"] is True
            for r in report["cases"].values()
        )
        write(args.output / "matrix.json", _signed(report))
    if report["acceptance_met"]:
        validate_matrix(_signed(report))
    return int(not report["acceptance_met"])


if __name__ == "__main__":
    raise SystemExit(main())
