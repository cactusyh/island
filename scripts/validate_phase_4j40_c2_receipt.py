"""C2 integrity is distinct from (currently blocked) scientific acceptance."""

import argparse
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

from island.analysis.equilibration import require
from island.analysis.equilibration_protocol import CorrectiveProtocol
from island.analysis.npt import FLAGS, _body, _digest
from island.analysis.npt_continuation import load_continuation
from island.exceptions import IslandError
from scripts.phase_4j40_c2_analysis import (
    observations,
    read,
    replica_analysis,
    trajectory,
)
from scripts.validate_phase_4j40 import CASES
from scripts.validate_phase_4j40_c2 import BASE, configuration
from scripts.validate_phase_4j40_receipt import sha

SOURCES = (
    "scripts/validate_phase_4j40_c2.py",
    "scripts/phase_4j40_c2_analysis.py",
    "scripts/phase_4j40_c2_resume.py",
    "scripts/phase_4j40_c2_finalize.py",
    "scripts/phase_4j40_c2_finite_difference.py",
    "scripts/validate_phase_4j40_c2_receipt.py",
    "tests/test_phase_4j40_c2.py",
)


def signed(record):
    require(record.get("identity") == _digest(_body(record)), "Changed record identity")
    for key, value in FLAGS.items():
        require(
            type(record.get(key)) is type(value) and record[key] == value,
            "Changed readiness",
        )


def history():
    """Git symlink blobs are link text, not the bytes of their target archive."""
    rows = {}
    entries = subprocess.check_output(
        ["git", "ls-tree", "-rz", BASE, "docs/evidence"]
    ).split(b"\0")
    for entry in entries:
        if not entry:
            continue
        meta, name = entry.split(b"\t")
        mode, _, blob = meta.decode().split()
        path = Path(name.decode())
        original = subprocess.check_output(["git", "cat-file", "blob", blob])
        actual = (
            str(path.readlink()).encode() if mode == "120000" else path.read_bytes()
        )
        require(actual == original, "Changed historical evidence: " + str(path))
        rows[str(path)] = {"git_blob": blob, "mode": mode, "resolved_sha256": sha(path)}
    return rows


def validate_run(root):
    from scripts.phase_4j38_fixture import build

    matrix = read(root / "matrix.json")
    signed(matrix)
    require(
        matrix["schema"] == "island.j40-c2-matrix.v1" and matrix["base_commit"] == BASE,
        "Changed matrix/base",
    )
    require(matrix["protocol"] == CorrectiveProtocol().payload, "Changed long protocol")
    require(set(matrix["cases"]) == set(CASES), "Missing benchmark cell")
    require(
        matrix["scientific_acceptance"] is False, "Blocked scientific gates relabeled"
    )
    for key, rows in matrix["cases"].items():
        require(len(rows) == 3, "Three attempted replicas required")
        fixture, evidence = build(CASES[key][0])
        roots = []
        for i, row in enumerate(rows):
            signed(row)
            root_replica = root / key / f"replica-{i}"
            require(row == read(root_replica / "result.json"), "Changed replica result")
            require(
                row["case"] == key
                and row["replica"] == i
                and row["charge_origin"] == CASES[key][1],
                "Changed cell/charge origin",
            )
            require(
                read(root_replica / "fixture.json") == json.loads(json.dumps(evidence)),
                "Changed exact fixture",
            )
            backend, config, frames, complete = trajectory(root_replica)
            require(
                backend["final_graph_identity"] == fixture.packed.graph.identity,
                "Stale final graph",
            )
            from types import SimpleNamespace

            require(
                config
                == configuration(SimpleNamespace(payload=backend), fixture, i).payload,
                "Changed protocol/seeds/configuration",
            )
            require(
                read(root_replica / "observations.json") == observations(root_replica),
                "Changed actual observations/statistics",
            )
            c = read(root_replica / "component-energies.json")
            signed(c)
            require(
                c["trajectory_identity"] == _digest(frames)
                and c["backend_identity"] == backend["identity"],
                "Stale components",
            )
            require(
                len(c["rows"]) == len(frames) and c["bound_kj_mol"] == 1e-5,
                "Changed component gate",
            )
            for frame, component in zip(frames, c["rows"], strict=True):
                error = abs(
                    sum(component["components_kj_mol"].values())
                    - frame["potential_energy_kj_mol"]
                )
                require(
                    component["step"] == frame["step"]
                    and component["sum_error_kj_mol"] == error
                    and component["passed"] == (error <= 1e-5),
                    "Changed component energy comparison",
                )
            require(
                c["passed"] == all(r["passed"] for r in c["rows"]),
                "False component success",
            )
            attempt = read(root_replica / "resumed-attempt/attempt.json")
            signed(attempt)
            _, restart = load_continuation(root_replica / "relocated")
            require(
                attempt["parent_restart"] == restart
                and attempt["config"] == config
                and attempt["backend_identity"] == backend["identity"],
                "Stale continuation",
            )
            require(
                attempt["frames"] == frames,
                "Split/uninterrupted frame/checkpoint/RNG mismatch",
            )
            require(
                attempt["status"] == ("completed" if complete else "failed"),
                "False continuation status",
            )
            comparison = read(root_replica / "partial-comparison.json")
            signed(comparison)
            require(
                comparison["frames_identity"] == _digest(frames)
                and comparison["byte_exact_retained_history"] is True
                and comparison["completed_workflow"] == complete,
                "False partial/completed equality",
            )
            roots.append(root_replica)
        require(
            matrix["replica_diagnostics"][key] == replica_analysis(roots),
            "Stale replica statistics",
        )
    for partition in ("ordinary", "crosslink"):
        fd = read(root / (partition + "-finite-difference.json"))
        signed(fd)
        require(
            fd["partition"] == partition and fd["bound"] == 1e-5,
            "Changed FD partition/bound",
        )
        require(
            fd["passed"] == all(c["passed"] for c in fd["components"].values()),
            "False FD success",
        )
    return matrix


def validate_receipt(path):
    receipt = read(path)
    signed(receipt)
    require(receipt["schema"] == "island.j40-c2-receipt.v1", "Wrong receipt schema")
    require(receipt["base_commit"] == BASE, "Changed receipt base")
    require(receipt["historical_files"] == history(), "Changed historical hashes")
    for source in SOURCES:
        raw = subprocess.check_output(
            ["git", "show", receipt["tested_implementation_commit"] + ":" + source]
        )
        from hashlib import sha256

        require(
            sha256(raw).hexdigest()
            == receipt["implementation_sources"][source]
            == sha(source),
            "Changed tested source",
        )
    for name, digest in receipt["files"].items():
        require(sha(path.parent / name) == digest, "Changed evidence file " + name)
    with tempfile.TemporaryDirectory(prefix="j40-c2-offline-") as directory:
        with tarfile.open(path.parent / "artifacts.tar.xz") as archive:
            archive.extractall(directory, filter="data")
        root = Path(directory)
        require(
            {str(p.relative_to(root)): sha(p) for p in root.rglob("*") if p.is_file()}
            == receipt["artifacts"],
            "Artifact inventory mismatch",
        )
        validate_run(root / "run")
        # Declared immutable J40 source asset supplies independently audited raw
        # RESP artifacts and native assignment records, never C2 trajectories.
        source_dir = root / "historical-source"
        source_dir.mkdir()
        with tarfile.open("docs/evidence/phase_4j40/artifacts.tar.xz") as archive:
            archive.extractall(source_dir, filter="data")
        from scripts.validate_phase_4j40_c1_receipt import (
            validate_run as validate_reference,
        )

        validate_reference(root / "run/reference", source_dir / "run")
        for case in CASES:
            original = read(source_dir / "run" / case / "backend.json")
            for i in range(3):
                fresh = read(root / "run" / case / f"replica-{i}" / "backend.json")
                require(
                    fresh["numerical_data"] == original["numerical_data"],
                    "Changed charge vector, masses, terms or parameters",
                )
                require(
                    {
                        k: v
                        for k, v in fresh["config"].items()
                        if k not in {"source", "identity"}
                    }
                    == {
                        k: v
                        for k, v in original["config"].items()
                        if k not in {"source", "identity"}
                    },
                    "Changed periodic model/settings",
                )
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("path", type=Path)
    p.add_argument("--bundle", action="store_true")
    p.add_argument("--diagnostic-only", action="store_true")
    a = p.parse_args()
    try:
        record = validate_run(a.path) if a.bundle else validate_receipt(a.path)
        print(
            json.dumps(
                {
                    "integrity": "passed",
                    "scientific_acceptance": record["scientific_acceptance"],
                }
            )
        )
        return 0 if a.diagnostic_only else int(not record["scientific_acceptance"])
    except (
        IslandError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AssertionError,
    ) as error:
        print(json.dumps({"integrity": "failed", "diagnostic": str(error)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
