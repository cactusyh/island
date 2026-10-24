"""Execute the C1 long protocol; scientific failures never suppress MD attempts."""

import argparse
import sys
import traceback
from pathlib import Path

from island.analysis.equilibration_protocol import CorrectiveProtocol
from island.analysis.equilibration_run import run_heated_workflow
from island.analysis.npt import FLAGS, NPTConfiguration, _signed, _write_json
from island.analysis.npt_continuation import (
    compare_continuations,
    relocate_continuation,
)
from scripts.validate_phase_4j38_c3 import command
from scripts.validate_phase_4j40 import CASES, backend_for

BASE = "c7b0060594db92b7fa1c741690d9bb8c1c825323"


def configuration(backend, fixture, replica):
    p = CorrectiveProtocol()
    return NPTConfiguration(
        family=backend.payload["config"]["family"],
        source_identity=backend.payload["config"]["source"]["sha256"],
        final_graph_identity=fixture.packed.graph.identity,
        packing_identity=fixture.packed.plan.identity,
        box_lengths_nm=tuple(x / 10 for x in fixture.packed.system.box.lengths),
        seed=4040 + replica * 10,
        velocity_seed=4041 + replica * 10,
        barostat_seed=4042 + replica * 10,
        nvt_steps=p.nvt_steps,
        npt_steps=p.npt_window_steps * p.maximum_windows,
        checkpoint_interval=p.frame_interval,
        pressure_interval=p.pressure_interval,
        cutoff_nm=0.5,
    )


def execute(case, replica, output, raw):
    """Independent preparation and integration per replica; no reference trajectory input.

    Execute the full predeclared maximum budget, rather than stop at a scientific
    preflight. Integration failure terminates that trajectory; other replicas
    must still be attempted. Pressure remains unavailable, never target-derived.
    """
    from scripts.phase_4j38_fixture import build

    output.mkdir(parents=True, exist_ok=False)
    row = dict(
        schema="island.j40-c2-replica.v1",
        case=case,
        replica=replica,
        protocol=CorrectiveProtocol().payload,
        charge_origin=CASES[case][1],
        scientific_acceptance=False,
        pressure_status="pressure_unavailable",
        continuation_status="not_executed",
        **FLAGS,
    )
    stage = "graph_charge_parameter_preparation"
    try:
        fixture, evidence = build(CASES[case][0])
        _write_json(output / "fixture.json", evidence)
        backend = backend_for(case, fixture, output, raw)
        _write_json(output / "backend.json", backend.payload)
        config = configuration(backend, fixture, replica)
        _write_json(output / "config.json", config.payload)
        row["backend_identity"] = backend.payload["identity"]
        row["config_identity"] = config.identity
        stage = "uninterrupted_dynamics"
        # Even a failed trajectory is retained by the production executor.
        failures = []
        for label, boundary in (("uninterrupted", None), ("split", 100)):
            try:
                run_heated_workflow(
                    backend,
                    config,
                    output / label,
                    stop_npt_step=boundary,
                    continuation=True,
                    heating=CorrectiveProtocol().heating,
                )
            except Exception as error:  # noqa: BLE001 -- preserve independent failures
                failures.append(
                    {
                        "stage": label,
                        "exception": type(error).__name__,
                        "diagnostic": str(error),
                    }
                )
                (output / (label + "-failure.log")).write_text(traceback.format_exc())
        if (output / "split/restart.json").exists():
            stage = "relocation_fresh_process_continuation"
            restart = relocate_continuation(output / "split", output / "relocated")
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
            code = command(argv, output / "resume")
            row["fresh_process_exit_code"] = code
            if code == 0 and (output / "uninterrupted/restart.json").exists():
                comparison = compare_continuations(
                    output / "uninterrupted", output / "resumed", output / "relocated"
                )
                _write_json(output / "comparison.json", comparison)
                row["continuation_status"] = (
                    "passed" if comparison["passed"] else "failed"
                )
            else:
                row["continuation_status"] = "blocked_integration_failure"
        row["integration_failures"] = failures
        row["software_execution"] = "failed" if failures else "completed"
    except Exception as error:  # noqa: BLE001 -- precise failure evidence
        row["software_execution"] = "failed"
        row["failure"] = {
            "stage": stage,
            "exception": type(error).__name__,
            "diagnostic": str(error),
        }
        (output / "failure.log").write_text(traceback.format_exc())
    _write_json(output / "result.json", _signed(row))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, choices=CASES)
    parser.add_argument("--replica", required=True, type=int, choices=range(3))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--raw", required=True, type=Path)
    args = parser.parse_args()
    result = execute(args.case, args.replica, args.output, args.raw)
    return int(not result["scientific_acceptance"])


if __name__ == "__main__":
    raise SystemExit(main())
