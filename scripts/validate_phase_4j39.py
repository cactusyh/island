"""J39 explicit OPLS-compatible RESP model; native aromatic OPLS stays rejected."""

import argparse
import json
import shutil
import sys
import tarfile
import traceback
from pathlib import Path

from island.analysis.npt import (
    FLAGS,
    NPTConfiguration,
    _digest,
    _signed,
    run_npt_workflow,
)
from island.analysis.npt_continuation import (
    compare_continuations,
    relocate_continuation,
)
from island.forcefields.ambertools.resp_policy import FiniteRESPCharges
from island.forcefields.oplsaa.resp import CHARGES, prepare
from island.forcefields.oplsaa.source import load_oplsaa_source
from scripts.phase_4j32_real_fixture import XML
from scripts.phase_4j38_c1_audit import aromatic
from scripts.phase_4j38_c2 import aromatic_decision
from scripts.phase_4j38_c3_charges import write
from scripts.phase_4j38_c3_reference import audit, sha
from scripts.phase_4j38_contract import require, stable
from scripts.phase_4j38_fixture import build
from scripts.phase_4j39_physics import check
from scripts.validate_phase_4j38 import ANALYSIS
from scripts.validate_phase_4j38_c3 import command

BASE = "78a46d46bd66921f527fc08082903d151f0e4f5e"


def matrix(result, decision):
    # Original cells remain separate, inherited evidence, never relabeled.
    archive_path = Path("docs/evidence/phase_4j38_c3/artifacts.tar.gz")
    require(
        sha(archive_path)
        == "1545fa0df9f15e997ebf10f369a3405d555d5891933195a7d7415766b7b8d1bc",
        "Changed historical C3 matrix archive",
    )
    with tarfile.open(archive_path) as archive:
        previous = json.load(archive.extractfile("run/matrix.json"))
    executed = result.get("software_acceptance_met") is True
    return _signed(
        dict(
            schema="island.j39-matrix.v1",
            benchmark_cell="oplsaa/aromatic",
            historical_c3_matrix_identity=previous["identity"],
            retained_cases=previous["cases"],
            retained_cases_executed_in_j39=False,
            native_source_decision_identity=decision["identity"],
            paths={
                "native_oplsaa": {
                    "status": "source_incompatible",
                    "charge_origin": "native_nonbonded",
                    "dynamics_executed": False,
                },
                "opls_compatible_provided": {
                    "status": "asset_unavailable",
                    "charge_origin": "provided",
                    "dynamics_executed": False,
                    "diagnostic": "No separately reviewed non-RESP provided vector declared",
                },
                "opls_compatible_resp_psi4": {
                    "status": "passed" if executed else result["status"],
                    "charge_origin": "RESP/Psi4",
                    "dynamics_executed": executed,
                    "result_identity": result["identity"],
                },
            },
            software_acceptance_met=executed,
            native_opls_acceptance_met=False,
            acceptance_met=False,
            equilibration_status="insufficient_sampling"
            if executed
            else "not_executed",
            scientific_validation_status="not_established",
            **FLAGS,
        )
    )


def run(output, reviewed):
    output.mkdir(parents=True, exist_ok=False)
    fixture, evidence = build("aromatic")
    evidence = json.loads(json.dumps(evidence))
    write(output / "fixture.json", evidence)
    raw_audit = aromatic(evidence, XML)
    decision = aromatic_decision(raw_audit, evidence)
    write(output / "native-source-audit.json", raw_audit)
    write(output / "native-source-decision.json", decision)
    result = dict(
        schema="island.j39-result.v1",
        label="OPLS-compatible parameters with RESP/Psi4 charges",
        charge_origin="RESP/Psi4",
        native_opls_charges=False,
        native_am1bcc_success=False,
        software_acceptance_met=False,
        dynamics_executed=False,
        status="blocked",
        **FLAGS,
    )
    stage = "raw_source_audit"
    try:
        require(not (output / "raw").exists(), "Source destination exists")
        shutil.copytree(reviewed, output / "raw")
        raw = output / "raw"
        require(
            json.loads((raw / "anchors.json").read_text())
            == json.loads(Path("scripts/phase_4j38_c3_anchors.json").read_text()),
            "Changed reviewed raw source anchors",
        )
        runtime = json.loads(
            Path("docs/evidence/phase_4j36/environment.json").read_text()
        )
        charges = json.loads((raw / "aromatic-charges.json").read_text())
        review = audit(
            raw,
            "aromatic",
            charges,
            runtime,
            json.loads((raw / "anchors.json").read_text()),
        )
        require(
            review == json.loads((raw / "aromatic-review.json").read_text()),
            "Changed raw source review",
        )
        write(output / "raw-review.json", review)
        provided = FiniteRESPCharges.load(
            raw / "aromatic-charges.json",
            fixture.packed.graph,
            expected_identity=CHARGES,
        )
        stage = "parameter_preparation"
        backend = prepare(
            fixture.packed.system,
            fixture.packed.graph,
            load_oplsaa_source(XML),
            provided,
        )
        backend.validate_integrity()
        write(output / "backend.json", backend.payload)
        result.update(
            charge_vector_identity=provided.identity,
            raw_review_identity=review["identity"],
            final_graph_identity=fixture.packed.graph.identity,
            packing_identity=fixture.packed.plan.identity,
            backend_identity=backend.payload["identity"],
        )
        q = provided.validate(fixture.packed.graph, expected_reference_identity=CHARGES)
        require(
            all(
                q[r["id"]].hex() == r["charge"].hex()
                for r in backend.payload["numerical_data"]["sites"]
            ),
            "Charge values changed",
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
        result["dynamics_executed"] = True
        stage = "independent_energy_force_reference"
        physics = check(backend, workflow["frames"])
        write(output / "physics.json", physics)
        stage = "split"
        run_npt_workflow(backend, config, output / "split", stop_npt_step=128)
        stage = "relocation"
        restart = relocate_continuation(output / "split", output / "relocated")
        stage = "fresh_process_continuation"
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
        require(
            command(argv, output / "resume") == 0,
            "Fresh-process continuation failed: see resume.log",
        )
        comparison = compare_continuations(
            output / "uninterrupted", output / "resumed", output / "relocated"
        )
        write(output / "comparison.json", comparison)
        require(comparison["passed"], "Split comparison failed")
        stage = "mdanalysis"
        require(
            command(
                [
                    str(ANALYSIS),
                    "-m",
                    "scripts.phase_4j38_observables",
                    "--case",
                    str(output),
                    "--output",
                    str(output / "observables.json"),
                ],
                output / "observables",
            )
            == 0,
            "MDAnalysis child failed",
        )
        import numpy as np

        frames = workflow["frames"]
        np.savez_compressed(
            output / "trajectory.npz",
            positions_nm=np.array([f["positions_nm"] for f in frames]),
            velocities_nm_ps=np.array([f["velocities_nm_ps"] for f in frames]),
            box_vectors_nm=np.array([f["box_vectors_nm"] for f in frames]),
            atom_ids=workflow["manifest"]["atom_ids"],
        )
        require(
            stable(frames, workflow["analysis"]), "Finite-run stability criteria failed"
        )
        require(
            workflow["analysis"]["density"]["blocks"] >= 4,
            "Missing four density blocks",
        )
        require(
            len({f["volume_nm3"] for f in frames if f["stage"] == "npt"}) > 1,
            "NPT volume did not change",
        )
        result.update(
            status="passed",
            software_acceptance_met=True,
            frame_count=len(frames),
            workflow_identity=workflow["manifest"]["identity"],
            trajectory_identity=_digest(frames),
            analysis_identity=json.loads((output / "observables.json").read_text())[
                "identity"
            ],
            comparison_identity=comparison["identity"],
            checkpoint_byte_equal=comparison["checkpoint_hashes"]
            == comparison["resumed_checkpoint_hashes"],
            rng_equal=comparison["rng_state_hashes"]
            == comparison["resumed_rng_state_hashes"],
            split_boundary=128,
            comparison_tolerances="exact",
            equilibration_status="insufficient_sampling",
            measured_pressure_status="unavailable",
        )
    except Exception as error:  # noqa: BLE001 -- preserve failures, no source substitution
        (output / "failure.log").write_text(traceback.format_exc())
        result.update(
            status="asset_unavailable"
            if isinstance(error, FileNotFoundError)
            else "blocked",
            diagnostic={
                "stage": stage,
                "exception": type(error).__name__,
                "message": str(error),
            },
        )
    result = _signed(result)
    write(output / "result.json", result)
    write(output / "matrix.json", matrix(result, decision))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reviewed", type=Path, required=True)
    parser.add_argument("--software-only", action="store_true")
    args = parser.parse_args()
    result = run(args.output, args.reviewed)
    return 0 if args.software_only and result["software_acceptance_met"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
