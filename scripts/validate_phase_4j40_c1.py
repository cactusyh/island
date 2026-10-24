"""Reproduce J40 failures, audit exact models, and fail closed before new MD."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

from island.analysis.equilibration_protocol import CorrectiveProtocol, assess
from island.analysis.npt import FLAGS, _signed, _write_json
from island.analysis.npt_continuation import load_continuation
from scripts.phase_4j38_contract import require
from scripts.phase_4j38_fixture import build
from scripts.phase_4j40_c1_finite_difference import audit
from scripts.phase_4j40_c1_parity import audit as parity
from scripts.phase_4j40_cross_engine import evaluate
from scripts.validate_phase_4j40 import CASES
from scripts.validate_phase_4j40_receipt import sha, validate_bundle

BASE = "a7dfb9a5aaf8ef37d94604e7771cb2771572d82c"


def execute(source, output, executable):
    import openmm as mm
    from openmm import unit

    require(not output.exists(), "Fresh empty output directory required")
    validate_bundle(
        source
    )  # strict raw-source, graph, config, trajectory and RNG validation
    output.mkdir(parents=True)
    before = {}
    for key in CASES:
        before[key] = {
            name: json.loads((source / key / (name + ".json")).read_text())
            for name in ("equilibration", "cross-engine")
        }
    _write_json(
        output / "before.json",
        _signed(dict(schema="island.j40-c1-before.v1", cases=before, **FLAGS)),
    )
    engine = json.loads((source / "engine.json").read_text())
    require(sha(executable) == engine["sha256"], "Changed KSPACE executable")
    rows = {}
    for key, (kind, _) in CASES.items():
        destination = output / key
        destination.mkdir()
        fixture, _ = build(kind)
        payload = json.loads((source / key / "backend.json").read_text())
        require(
            payload["final_graph_identity"] == fixture.packed.graph.identity,
            "Stale final graph",
        )
        backend = SimpleNamespace(payload=payload, system=fixture.packed.system)
        # Source was validated offline above; this adapter is only for single points.
        original = json.dumps(payload, sort_keys=True)
        backend.validate_integrity = lambda p=payload, s=original: require(
            json.dumps(p, sort_keys=True) == s, "Changed numerical payload"
        )
        reference = evaluate(
            backend, mm, unit, destination / "cross-engine", executable
        )
        require(
            sha(destination / "cross-engine/system.data")
            == sha(source / key / "cross-engine/system.data")
            and sha(destination / "cross-engine/in.periodic")
            == sha(source / key / "cross-engine/in.periodic"),
            "Different graph, masses, charges, parameters, coordinates, box or engine settings",
        )
        _write_json(destination / "cross-engine.json", reference)
        _write_json(
            destination / "input-parity.json",
            parity(
                payload,
                {
                    i: fixture.packed.system.coordinates.get(i)
                    for i in sorted(fixture.packed.system.topology.sites)
                },
                destination / "cross-engine",
            ),
        )
        if key == "opls-network":
            ids = sorted(fixture.packed.system.topology.sites)
            finite = audit(
                payload,
                [fixture.packed.system.coordinates.get(i) / 10 for i in ids],
                ids,
            )
            _write_json(destination / "packed-finite-difference.json", finite)
            minimized = json.loads(
                (source / key / "replica-0/uninterrupted/frames.json").read_text()
            )["frames"][0]
            require(minimized["stage"] == "minimize", "Wrong geometry stage")
            minimized_fd = audit(payload, minimized["positions_nm"], ids)
            _write_json(destination / "minimized-finite-difference.json", minimized_fd)
        else:
            finite = {
                "passed": reference["checks"]["finite_difference"],
                "scope": "historical_first_and_last_atom_all_axes",
                "not_a_full_component_audit": True,
            }
        workflows = [
            load_continuation(source / key / f"replica-{i}/uninterrupted")[0]
            for i in range(2)
        ]
        result = assess(
            workflows,
            reference=reference,
            finite_difference=finite,
            pressure={
                "status": "pressure_unavailable",
                "reason": "OpenMM state exposes no complete instantaneous pressure; no independently validated virial estimator",
            },
        )
        result = _signed(
            dict(
                schema="island.j40-c1-cell.v1",
                case=key,
                assessment=result,
                new_dynamics="not_executed",
                reason="Prerequisites unresolved; longer dynamics cannot repair reference or pressure measurement",
                charge_origin=CASES[key][1],
                graph_identity=fixture.packed.graph.identity,
                historical_backend_sha256=sha(source / key / "backend.json"),
                **FLAGS,
            )
        )
        _write_json(destination / "result.json", result)
        rows[key] = result
    # Full nested offline validation in a different path and a fresh interpreter.
    # Historical files are copied byte-for-byte only for validation, never claimed as new MD.
    relocated = output / "historical-relocated"
    shutil.copytree(source, relocated)
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_phase_4j40_receipt",
            str(relocated),
            "--bundle",
            "--diagnostic-only",
        ],
        capture_output=True,
        check=False,
        text=True,
    )
    _write_json(
        output / "child.json",
        {
            "argv": proc.args,
            "exit_code": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        },
    )
    require(proc.returncode == 0, "Relocated separate-process reconstruction failed")
    shutil.rmtree(relocated)  # hash-bound historical source remains archived in J40
    record = _signed(
        dict(
            schema="island.j40-c1-matrix.v1",
            base_commit=BASE,
            protocol=CorrectiveProtocol().payload,
            historical_matrix_sha256=sha(source / "matrix.json"),
            engine_sha256=sha(executable),
            cases=rows,
            scientific_acceptance=False,
            new_dynamics="not_executed",
            relocation_reconstruction="passed",
            **FLAGS,
        )
    )
    _write_json(output / "matrix.json", record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    args = parser.parse_args()
    execute(args.source, args.output, args.lammps.resolve())
    return 1  # faithfully blocked scientific acceptance


if __name__ == "__main__":
    raise SystemExit(main())
