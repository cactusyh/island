"""Strict offline J40 bundle and publication receipt validation.

Diagnostic validation returns zero for a faithfully blocked run; acceptance
returns one until all declared gates pass. Invalid evidence returns two.
"""

import argparse
import json
import subprocess
import tarfile
import tempfile
from hashlib import sha256
from pathlib import Path

from island.analysis.equilibration import EquilibrationPolicy, validate_analysis
from island.analysis.equilibration_run import HeatingSchedule
from island.analysis.npt import FLAGS, _body, _digest
from island.analysis.npt_continuation import compare_continuations, load_continuation
from scripts.phase_4j38_contract import require
from scripts.validate_phase_4j40 import BASE, CASES

SOURCES = (
    "src/island/analysis/equilibration.py",
    "src/island/analysis/equilibration_run.py",
    "scripts/phase_4j40_cross_engine.py",
    "scripts/validate_phase_4j40.py",
    "scripts/validate_phase_4j40_receipt.py",
    "tests/test_phase_4j40.py",
)


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(value, schema):
    require(
        value.get("schema") == schema
        and value.get("identity") == _digest(_body(value)),
        "Changed " + schema,
    )
    for k, v in FLAGS.items():
        require(type(value.get(k)) is type(v) and value[k] == v, "Changed readiness")


def validate_bundle(root):
    root = Path(root)
    matrix = read(root / "matrix.json")
    signed(matrix, "island.j40-matrix.v1")
    require(
        matrix["base_commit"] == BASE and set(matrix["cases"]) == set(CASES),
        "Changed benchmark matrix",
    )
    engine = read(root / "engine.json")
    require(
        engine["identity"] == _digest(_body(engine))
        and engine["schema"] == "island.j40-engine.v1",
        "Engine identity changed",
    )
    require(engine["identity"] == matrix["engine_identity"], "Stale engine identity")
    if engine["status"] == "available":
        require(
            engine["kspace"] and sha(root / "lammps-help.txt") == engine["help_sha256"],
            "Missing KSPACE evidence",
        )
    from scripts.phase_4j38_c3_reference import audit

    runtime = read("docs/evidence/phase_4j36/environment.json")
    anchors = read(root / "raw/anchors.json")
    require(
        anchors == read("scripts/phase_4j38_c3_anchors.json"),
        "Changed raw source anchors",
    )
    for kind in ("pe", "aromatic"):
        key = "gaff2-pe" if kind == "pe" else "opls-resp-aromatic"
        review = audit(
            root / "raw",
            kind,
            read(root / ("raw/" + kind + "-charges.json")),
            runtime,
            anchors,
        )
        require(
            review == read(root / key / "raw-review.json"),
            "Stale independent source review",
        )
    pcff = read(root / "controls/pcff/result.json")
    require(
        pcff["acceptance_met"]
        and pcff["dynamics_executed"] is False
        and pcff["diagnostic"]["classification"] == "expected_typing_rejection",
        "PCFF unexpectedly accepted",
    )
    native = read(root / "controls/native-am1bcc.json")
    require(
        native["status"] == "blocked"
        and native["dynamics_executed"] is False
        and sha(native["historical_evidence"]) == native["sha256"],
        "AM1-BCC limitation relabeled",
    )
    for key, row in matrix["cases"].items():
        case = root / key
        require(row["case"] == key, "Changed benchmark cell label")
        from scripts.phase_4j38_fixture import build

        fixture, evidence = build(CASES[key][0])
        require(
            read(case / "fixture.json") == json.loads(json.dumps(evidence)),
            "Changed benchmark fixture",
        )
        require(
            row["graph_identity"] == fixture.packed.graph.identity
            and row["packing_identity"] == fixture.packed.plan.identity,
            "Changed exact benchmark graph",
        )
        signed(row, "island.j40-cell.v1")
        require(
            row == read(case / "result.json") and row["charge_origin"] == CASES[key][1],
            "Changed charge origin or cell result",
        )
        require(
            row["software_execution"] == "passed",
            "Incomplete software workflow: " + key,
        )
        workflows = []
        for i in range(2):
            replica = case / ("replica-" + str(i))
            w, _ = load_continuation(replica / "uninterrupted")
            comparison = compare_continuations(
                replica / "uninterrupted", replica / "resumed", replica / "relocated"
            )
            require(
                comparison["passed"]
                and comparison == read(replica / "comparison.json"),
                "Changed continuation comparison",
            )
            require(
                read(replica / "split/restart.json")
                == read(replica / "relocated/restart.json"),
                "Relocation identity mismatch",
            )
            cfg = w["manifest"]["config"]
            require(
                (cfg["seed"], cfg["velocity_seed"], cfg["barostat_seed"])
                == (4040 + i * 10, 4041 + i * 10, 4042 + i * 10),
                "Changed replica seeds",
            )
            require(
                cfg["nvt_steps"] == 448
                and cfg["npt_steps"] == 1024
                and cfg["checkpoint_interval"] == 4
                and cfg["pressure_interval"] == 25,
                "Changed declared protocol",
            )
            schedule = HeatingSchedule(
                ((64, 100.0), (128, 200.0), (192, 300.0)), 448, 300.0
            )
            require(
                w["manifest"]["heating"] == schedule.payload, "Changed heating schedule"
            )
            for f in w["frames"]:
                if f["stage"] in ("nvt", "minimize"):
                    require(
                        f["target_temperature_kelvin"]
                        == schedule.temperature(f["step"]),
                        "Heating targets changed",
                    )
            require(
                w["manifest"]["backend_identity"] == row["backend_identity"]
                and w["manifest"]["final_graph_identity"] == row["graph_identity"],
                "Stale cell binding",
            )
            obs = read(replica / "observables.json")
            require(
                obs["identity"] == _digest(_body(obs)), "Observable identity changed"
            )
            require(
                obs["frame_content_identity"] == _digest(w["frames"])
                and obs["workflow_identity"] == w["manifest"]["identity"],
                "Stale observable trajectory binding",
            )
            workflows.append(w)
        diagnostic = read(case / "equilibration.json")
        require(
            diagnostic["policy"] == EquilibrationPolicy().payload,
            "Undeclared statistical thresholds",
        )
        validate_analysis(diagnostic, workflows)
        require(
            row["diagnostic_identity"] == diagnostic["identity"]
            and row["equilibration_status"] == diagnostic["equilibration_status"],
            "False equilibration claim",
        )
        reference = read(case / "cross-engine.json")
        require(
            reference["executable_identity"] == engine["identity"],
            "Changed reference executable",
        )
        if reference.get("passed"):
            require(
                reference["lammps_exit_code"] == 0
                and all(reference["checks"].values()),
                "False cross-engine success",
            )
            require(
                reference["energy_error_kj_mol"] <= 1e-5
                and reference["force_max_error_kj_mol_angstrom"] <= 1e-5
                and max(reference["component_errors_kj_mol"].values()) <= 1e-5,
                "Reference mismatch labeled successful",
            )
        require(
            row["acceptance_met"]
            == (diagnostic["acceptance_met"] and reference["passed"]),
            "False cell acceptance",
        )
    require(
        matrix["acceptance_met"]
        == all(r["acceptance_met"] for r in matrix["cases"].values()),
        "False matrix acceptance",
    )
    return matrix


def validate_receipt(path):
    path = Path(path)
    root = path.parent
    receipt = read(path)
    signed(receipt, "island.j40-receipt.v1")
    require(receipt["base_commit"] == BASE, "Wrong base")
    manifest = read(root / "manifest.json")
    require(
        manifest["schema"] == "island.j40-manifest.v1"
        and sha(root / "manifest.json") == receipt["manifest_sha256"],
        "Changed manifest",
    )
    for name, value in manifest["files"].items():
        require(
            not Path(name).is_absolute()
            and ".." not in Path(name).parts
            and sha(root / name) == value,
            "Changed artifact: " + name,
        )
    for name in SOURCES:
        source = subprocess.check_output(
            ["git", "show", receipt["implementation_commit"] + ":" + name]
        )
        require(source == Path(name).read_bytes(), "Untested implementation: " + name)
    for name, digest in read(root / "historical-hashes.json").items():
        require(sha(name) == digest, "Historical evidence changed: " + name)
    for gate in (
        "focused",
        "j22-j39",
        "full_supported",
        "ruff",
        "compileall",
        "diff_check",
        "pip_supported",
        "pip_foyer",
        "pip_analysis",
        "pip_amber",
        "pip_quantum",
    ):
        command = read(root / "validation" / (gate + "-command.json"))
        require(
            command["exit_code"] == 0
            and command["implementation_commit"] == receipt["implementation_commit"],
            "Unpassed validation: " + gate,
        )
    with tempfile.TemporaryDirectory(prefix="j40-receipt-") as temporary:
        with tarfile.open(root / "artifacts.tar.gz") as archive:
            archive.extractall(temporary, filter="data")
        matrix = validate_bundle(Path(temporary) / "run")
    require(
        receipt["matrix_identity"] == matrix["identity"]
        and receipt["acceptance_met"] == matrix["acceptance_met"],
        "Receipt/matrix mismatch",
    )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--bundle", action="store_true")
    parser.add_argument("--diagnostic-only", action="store_true")
    args = parser.parse_args()
    try:
        result = (
            validate_bundle(args.path) if args.bundle else validate_receipt(args.path)
        )
        print(
            json.dumps(
                {"valid": True, "acceptance_met": result["acceptance_met"], **FLAGS}
            )
        )
        return 0 if args.diagnostic_only or result["acceptance_met"] else 1
    except Exception as error:  # noqa: BLE001 -- explicit invalid evidence exit
        print(json.dumps({"valid": False, "diagnostic": str(error)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
