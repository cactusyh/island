"""J38 real-polymer benchmark. Incomplete software/scientific gates exit nonzero."""

import argparse
import json
import sys
import traceback
from pathlib import Path

from island.analysis.npt import (
    FLAGS,
    NPTConfiguration,
    _digest,
    _signed,
    _write_json,
    run_npt_workflow,
)
from island.analysis.npt_continuation import (
    compare_continuations,
    load_continuation,
    relocate_continuation,
)
from scripts.phase_4j32_assignment import independent_check, prepare
from scripts.phase_4j32_real_fixture import ROOT
from scripts.phase_4j38_contract import (
    BASE,
    CASES,
    base_row,
    gate,
    require,
    stable,
    validate_matrix,
    validate_row,
)
from scripts.phase_4j38_fixture import build
from scripts.phase_4j38_physics import check
from scripts.validate_phase_4j37 import command

ANALYSIS = ROOT / "island-validation/phase4j32-analysis-env/bin/python"


def input_provenance(fixture, evidence):
    return {
        "psmiles": evidence["psmiles"],
        "degree_of_polymerization": evidence["dp"],
        "end_group_policy": evidence["end_group_policy"],
        "stereochemistry": evidence["stereochemistry"],
        "construction_seeds": evidence["construction_seeds"],
        "final_graph_identity": fixture.packed.graph.identity,
        "packing_identity": fixture.packed.plan.identity,
        "box": {
            "lengths_angstrom": list(fixture.packed.system.box.lengths),
            "periodic": list(fixture.packed.system.box.periodic),
        },
        "repeat_unit_provenance": fixture.packed.graph.payload["repeat_provenance"],
        "crosslink_plan_identity": fixture.crosslinked.plan.identity
        if fixture.crosslinked
        else None,
    }


def finish(row, output):
    row["execution_met"] = row["gates"]["software_execution"]["status"] == "passed"
    row["software_acceptance_met"] = all(
        row["gates"][g]["status"] == "passed"
        for g in (
            "software_execution",
            "source_completeness",
            "numerical_reference",
            "finite_time_stability",
        )
    )
    row = _signed(row)
    validate_row(row["case"], row)
    _write_json(output / "result.json", row)
    return row


def postprocess(output, row, fixture, evidence):
    workflow, restart = load_continuation(output / "uninterrupted")
    comparison = compare_continuations(
        output / "uninterrupted", output / "resumed", output / "relocated"
    )
    require(comparison["passed"], "Split/uninterrupted comparison failed")
    _write_json(output / "comparison-j38.json", comparison)
    argv = [
        str(ANALYSIS),
        "-m",
        "scripts.phase_4j38_observables",
        "--case",
        str(output),
        "--output",
        str(output / "observables.json"),
    ]
    require(
        command(argv, output / "observables") == 0,
        "MDAnalysis reconstruction failed; see observables.log",
    )
    # Actual MD frames exported as portable NPZ, with frame-order/box metadata.
    import numpy as np

    frames = workflow["frames"]
    np.savez_compressed(
        output / "trajectory.npz",
        positions_nm=np.array([f["positions_nm"] for f in frames]),
        box_vectors_nm=np.array([f["box_vectors_nm"] for f in frames]),
        atom_ids=np.array(workflow["manifest"]["atom_ids"]),
    )
    observables = json.loads((output / "observables.json").read_text())
    backend = json.loads((output / "uninterrupted/backend.json").read_text())
    require(
        backend["final_graph_identity"] == fixture.packed.graph.identity
        and backend["config"]["family"]
        == ("OPLS-AA" if row["family"] == "oplsaa" else "GAFF2"),
        "Wrong graph/family",
    )
    row.update(
        final_graph_identity=fixture.packed.graph.identity,
        packing_identity=fixture.packed.plan.identity,
        source_identity=workflow["manifest"]["source_identity"],
        assignment_identity=backend["assignment_identity"],
        typed_graph_identity=backend["typed_graph_identity"],
        charge_identity=backend["charge_identity"],
        workflow_identity=workflow["manifest"]["identity"],
        trajectory_identity=_digest(frames),
        checkpoint_identity=restart["checkpoint_identity"],
        analysis_identity=observables["identity"],
        comparison_identity=comparison["identity"],
        frame_count=len(frames),
        density_blocks=workflow["analysis"]["density"]["blocks"],
        split_equal=True,
        checkpoint_byte_equal=comparison["checkpoint_hashes"]
        == comparison["resumed_checkpoint_hashes"],
        rng_equal=comparison["rng_state_hashes"]
        == comparison["resumed_rng_state_hashes"],
        relocation=True,
        fresh_process=True,
        box_volume_changed=len({f["volume_nm3"] for f in frames if f["stage"] == "npt"})
        > 1,
        setup=evidence["setup"],
        charge_origin="RESP/Psi4" if row["family"] == "gaff2" else "OPLS-AA source",
    )
    require(row["box_volume_changed"], "NPT volume relaxation not demonstrated")
    row["gates"]["software_execution"] = gate(
        "passed",
        "Actual minimize/NVT/NPT, exact split comparison, relocated fresh-process continuation and MDAnalysis completed",
    )
    row["gates"]["finite_time_stability"] = gate(
        "passed" if stable(frames, workflow["analysis"]) else "failed",
        "Declared finite values, 0<T<2000 K, density max/min<2 and four computable blocks; no equilibrium inference",
    )
    return row


def worker(output, key):
    output.mkdir(parents=True, exist_ok=False)
    row = base_row(key)
    stage = "fixture"
    try:
        fixture, evidence = build(row["benchmark"])
        _write_json(output / "fixture.json", evidence)
        row["input_provenance"] = input_provenance(fixture, evidence)

        if row["family"] != "oplsaa":
            row["gates"]["source_completeness"] = gate(
                "blocked",
                "No independently reviewed graph-specific RESP/Psi4 raw ESP/charge artifacts for this new polymer; retained PEO charges cannot be transferred",
            )
            row["diagnostic"] = {
                "stage": "charge_provenance",
                "classification": "unverified_graph_specific_charge_source",
                "final_graph_identity": fixture.packed.graph.identity,
            }
            return finish(row, output)
        stage = "native_assignment"
        backend = prepare(fixture, "oplsaa", output)
        _write_json(output / "independent-assignment.json", independent_check(backend))
        row["gates"]["source_completeness"] = gate(
            "passed",
            "Exact final-graph reassignment; pinned source and independent upstream force inventory/path verification",
        )
        config = NPTConfiguration(
            family="OPLS-AA",
            source_identity=backend.payload["config"]["source"]["sha256"],
            final_graph_identity=fixture.packed.graph.identity,
            packing_identity=fixture.packed.plan.identity,
            box_lengths_nm=tuple(x / 10 for x in fixture.packed.system.box.lengths),
            nvt_steps=64,
            npt_steps=256,
        )
        stage = "uninterrupted"
        workflow = run_npt_workflow(
            backend, config, output / "uninterrupted", continuation=True
        )
        stage = "independent_reference"
        _write_json(output / "physics.json", check(backend, workflow["frames"]))
        row["gates"]["numerical_reference"] = gate(
            "passed",
            "Pinned Foyer.apply/ParmEd independent parameter construction; periodic Ewald energies/forces/components first and last actual frame; same OpenMM Reference engine",
        )
        stage = "split"
        run_npt_workflow(
            backend, config, output / "split", stop_npt_step=config.npt_steps // 2
        )
        stage = "relocation"
        restart = relocate_continuation(output / "split", output / "relocated")
        stage = "fresh_child_resume"
        argv = [
            sys.executable,
            "-m",
            "scripts.validate_phase_4j33_c1",
            "--resume",
            str(output / "relocated"),
            "--destination",
            str(output / "resumed"),
            "--restart",
            restart["identity"],
            "--parent",
            restart["parent_workflow_identity"],
            "--config",
            config.identity,
        ]
        require(command(argv, output / "resume") == 0, "Fresh-process resume failed")
        stage = "analysis"
        row = postprocess(output, row, fixture, evidence)
    except Exception as error:  # noqa: BLE001 -- preserve failed science/software gates
        (output / "failure.log").write_text(traceback.format_exc())
        row["diagnostic"] = {
            "stage": stage,
            "exception": type(error).__name__,
            "message": str(error),
        }
        row["gates"]["software_execution"] = gate("failed", str(error))
        if stage == "native_assignment":
            row["gates"]["source_completeness"] = gate("failed", str(error))
    return finish(row, output)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    matrix = dict(
        schema="island.j38-matrix.v1",
        base_commit=BASE,
        cases={},
        acceptance_met=False,
        software_acceptance_met=False,
        **FLAGS,
    )
    # J37 replay includes independent raw provenance review and external native RESP refit.
    retained = output / "verified-resp"
    code = command(
        [
            sys.executable,
            "-m",
            "scripts.validate_phase_4j37",
            "--output",
            str(retained),
        ],
        output / "verified-resp",
    )
    for key in CASES:
        family, benchmark = key.split("/")
        case = output / key.replace("/", "-")
        if family != "gaff2" or benchmark not in {"peo", "network"}:
            row = worker(case, key)
        else:
            case.mkdir()
            row = base_row(key)
            try:
                require(
                    code == 0,
                    "J37 raw provenance/refit and fresh GAFF2 workflow failed",
                )
                kind = "linear" if benchmark == "peo" else "network"
                import shutil

                shutil.copytree(retained / kind, case, dirs_exist_ok=True)
                fixture, evidence = build(benchmark)
                _write_json(case / "fixture.json", evidence)
                row["input_provenance"] = input_provenance(fixture, evidence)
                row["gates"]["source_completeness"] = gate(
                    "passed",
                    "Fresh GAFF2 preparation gated by J37 independent raw RESP/Psi4 review and native refit; graph-specific unchanged charges",
                )
                physics = json.loads((retained / (kind + "-physics.json")).read_text())
                require(physics["passed"], "Periodic Amber reference failed")
                _write_json(case / "physics.json", physics)
                row["gates"]["numerical_reference"] = gate(
                    "passed",
                    "Independent AmberPrmtopFile periodic Ewald construction; same OpenMM Reference engine",
                )
                row = postprocess(case, row, fixture, evidence)
            except Exception as error:  # noqa: BLE001
                (case / "failure.log").write_text(traceback.format_exc())
                row["diagnostic"] = {"stage": "j37_revalidation", "message": str(error)}
            row = finish(row, case)
        matrix["cases"][key] = row
        _write_json(output / "partial-matrix.json", _signed(matrix))
    # Retain actual negative controls separately from positive benchmark inventory.
    controls = output / "controls"
    controls.mkdir()
    from scripts.validate_phase_4j32 import worker as pcff_worker

    matrix["pcff_control"] = pcff_worker(controls / "pcff", "pcff", "network")
    # The native-charge limitation remains an independently anchored historical diagnostic.
    matrix["native_am1bcc"] = {
        "status": "blocked",
        "acceptance_met": False,
        "diagnostic": "Native AmberTools 24.8 precision limitation retained; no native AM1-BCC success",
        "evidence": "docs/evidence/phase_4j34_c2/receipt.json",
    }
    matrix["software_acceptance_met"] = all(
        r["software_acceptance_met"] for r in matrix["cases"].values()
    )
    matrix = _signed(matrix)
    validate_matrix(matrix)
    _write_json(output / "matrix.json", matrix)
    return matrix


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--case", choices=CASES)
    args = p.parse_args()
    result = worker(args.output, args.case) if args.case else run(args.output)
    return int(not result["acceptance_met"])


if __name__ == "__main__":
    raise SystemExit(main())
