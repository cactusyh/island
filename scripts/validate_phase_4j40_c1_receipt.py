"""Offline corrective evidence integrity (0), blocked acceptance (1), invalid (2)."""

import argparse
import json
import tarfile
import tempfile
from pathlib import Path

from island.analysis.equilibration_protocol import (
    CorrectiveProtocol,
    validate_assessment,
)
from island.analysis.npt import FLAGS, _body, _digest
from island.analysis.npt_continuation import load_continuation
from island.exceptions import IslandError
from scripts.phase_4j38_contract import require
from scripts.validate_phase_4j40 import CASES
from scripts.validate_phase_4j40_c1 import BASE
from scripts.validate_phase_4j40_receipt import sha

SOURCES = (
    "src/island/analysis/equilibration_protocol.py",
    "scripts/phase_4j40_c1_finite_difference.py",
    "scripts/phase_4j40_c1_parity.py",
    "scripts/validate_phase_4j40_c1.py",
    "scripts/validate_phase_4j40_c1_receipt.py",
    "tests/test_phase_4j40_c1.py",
)


def read(path):
    return json.loads(Path(path).read_text())


def signed(record):
    require(record.get("identity") == _digest(_body(record)), "Changed record identity")
    for k, v in FLAGS.items():
        require(type(record.get(k)) is type(v) and record[k] == v, "Changed readiness")


def validate_run(root, source):
    matrix = read(root / "matrix.json")
    signed(matrix)
    require(
        matrix["base_commit"] == BASE
        and matrix["protocol"] == CorrectiveProtocol().payload,
        "Changed base or protocol",
    )
    require(
        matrix["historical_matrix_sha256"] == sha(source / "matrix.json"),
        "Stale historical matrix",
    )
    require(set(matrix["cases"]) == set(CASES), "Changed acceptance matrix")
    require(
        matrix["new_dynamics"] == "not_executed"
        and not matrix["scientific_acceptance"],
        "Unexecuted corrective protocol relabeled successful",
    )
    require(read(root / "child.json")["exit_code"] == 0, "Relocation child failed")
    for key in CASES:
        cell = read(root / key / "result.json")
        signed(cell)
        require(
            cell == matrix["cases"][key] and cell["charge_origin"] == CASES[key][1],
            "Changed cell or charge origin",
        )
        require(
            cell["historical_backend_sha256"] == sha(source / key / "backend.json"),
            "Changed backend",
        )
        reference = read(root / key / "cross-engine.json")
        for name in ("system.data", "in.periodic"):
            require(
                sha(root / key / "cross-engine" / name)
                == sha(source / key / "cross-engine" / name),
                "Different graph, mass, charge, geometry or force-field settings",
            )
        finite = (
            read(root / key / "packed-finite-difference.json")
            if key == "opls-network"
            else {
                "passed": reference["checks"]["finite_difference"],
                "scope": "historical_first_and_last_atom_all_axes",
                "not_a_full_component_audit": True,
            }
        )
        if key == "opls-network":
            for name in ("packed", "minimized"):
                record = read(root / key / (name + "-finite-difference.json"))
                signed(record)
                require(record["bound"] == 1e-5, "Relaxed finite-difference bound")
                for component in record["components"].values():
                    require(
                        component["passed"]
                        == bool(component["converged_intervals_nm"]),
                        "False finite-difference success",
                    )
        workflows = [
            load_continuation(source / key / f"replica-{i}/uninterrupted")[0]
            for i in range(2)
        ]
        validate_assessment(
            cell["assessment"],
            workflows,
            reference=reference,
            finite_difference=finite,
            pressure=cell["assessment"]["pressure"],
        )
        require(
            not cell["assessment"]["scientific_acceptance"],
            "Blocked prerequisites relabeled success",
        )
    return matrix


def validate_receipt(path):
    receipt = read(path)
    signed(receipt)
    require(
        receipt["base_commit"] == BASE and not receipt["scientific_acceptance"],
        "Changed acceptance",
    )
    directory = path.parent
    for name, digest in receipt["files"].items():
        require(sha(directory / name) == digest, "Changed evidence file: " + name)
    for name, digest in receipt["historical_hashes"].items():
        require(sha(name) == digest, "Changed historical file: " + name)
    for name, digest in receipt["source_hashes"].items():
        require(sha(name) == digest, "Changed corrective source: " + name)
    require(set(receipt["source_hashes"]) == set(SOURCES), "Missing corrective sources")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        with tarfile.open("docs/evidence/phase_4j40/artifacts.tar.xz") as archive:
            archive.extractall(tmp / "historical", filter="data")
        with tarfile.open(directory / "artifacts.tar.xz") as archive:
            archive.extractall(tmp / "corrective", filter="data")
        validate_run(tmp / "corrective/run", tmp / "historical/run")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--diagnostic-only", action="store_true")
    args = parser.parse_args()
    try:
        if args.source:
            validate_run(args.path, args.source)
        else:
            validate_receipt(args.path)
    except (
        IslandError,
        ValueError,
        OSError,
        KeyError,
        TypeError,
        tarfile.TarError,
    ) as exc:
        print(type(exc).__name__ + ": " + str(exc))
        return 2
    print(
        "Integrity valid; scientific acceptance blocked; no corrective dynamics claimed"
    )
    return 0 if args.diagnostic_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
