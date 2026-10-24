"""Reconstruct C12 comparison completions; scientific acceptance stays nonzero."""

import argparse
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from scripts.phase_4j40_c11_numeric import (
    KJ_MOL_NM3_TO_BAR,
    parse_dump,
    require,
    sha,
    signed,
)
from scripts.phase_4j40_c12 import comparison, completion, corrected
from scripts.validate_phase_4j40_c11_receipt import validate_study as validate_c11

BASE = "ac0c38185e7b316b68109cd8be24e00c897f3fa8"
SOURCES = [
    "scripts/phase_4j40_c12.py",
    "scripts/phase_4j40_c12_engine.py",
    "scripts/phase_4j40_c12_study.py",
    "scripts/validate_phase_4j40_c12_receipt.py",
    "tests/test_phase_4j40_c12.py",
]


def read(path):
    return json.loads(Path(path).read_text())


def historical_hashes():
    result = {}
    for row in subprocess.check_output(
        ["git", "ls-tree", "-rz", BASE, "docs/evidence"]
    ).split(b"\0"):
        if not row:
            continue
        meta, path = row.split(b"\t")
        mode, _, blob = meta.decode().split()
        path = Path(path.decode())
        actual = (
            str(path.readlink()).encode() if mode == "120000" else path.read_bytes()
        )
        expected = subprocess.check_output(["git", "cat-file", "blob", blob])
        require(actual == expected, "Changed historical evidence: " + str(path))
        result[str(path)] = {"git_blob": blob, "mode": mode, "sha256": sha(path)}
    return result


def check_native(root, row):
    for engine in ("openmm", "lammps"):
        native = read(root / engine / "point.json")
        require(native == row[engine], "Changed native engine record")
    native = row["lammps"]
    dump = parse_dump(root / "lammps/forces.dump")
    require(
        dump["forces_kj_mol_nm"] == native["forces_kj_mol_nm"],
        "Changed native force vector",
    )
    require(dump["charges_e"] == native["charges_e"], "Changed native charge vector")
    require(dump["atom_ids"] == native["atom_ids"], "Changed native atom mapping")
    require(
        sha(root / "openmm/system.xml") == row["openmm"]["system_xml_sha256"],
        "Changed OpenMM model",
    )
    for name, digest in native["raw_files"].items():
        require(sha(root / "lammps" / name) == digest, "Changed raw LAMMPS output")


def validate_study(root):
    root = Path(root)
    result = read(root / "study.json")
    require(
        result == signed(result) and result["schema"] == "island.j40-c12-study.v1",
        "Changed study",
    )
    require(
        result["scientific_acceptance"] is False
        and result["pressure_status"] == "pressure_unavailable"
        and result["ewald_status"] == "ewald_unvalidated"
        and result["long_npt"] == "not_executed"
        and result["cutoff_recovery"] == "not_executed"
        and result["production_validated"] is False
        and result["simulation_readiness"] == "not_established",
        "False readiness claim",
    )
    validate_c11(root / "c11-rerun")
    native = read(root / "c11-rerun/study.json")
    require(native["identity"] == result["c11_original_identity"], "Changed C11 matrix")
    b = read(root / "c11-rerun/inputs/backend.json")
    completed = read(root / "c11-completed.json")
    require(len(completed) == 38, "Incomplete C11 matrix")
    for item, row in zip(completed, native["rows"], strict=True):
        require(item["name"] == row["name"], "Changed C11 ordering")
        f = read(root / "c11-rerun" / row["name"] / "frame.json")
        c = completion(b, f, row["lammps"]["setting"], row["lammps"])
        require(
            item["completion"] == c
            and item["comparison"] == comparison(row["openmm"], row["lammps"], c),
            "Changed analytical completion or comparison",
        )
    rows = read(root / "additional-matrix.json")
    require(
        len(rows) == 42 and result["additional_rows"] == 42,
        "Incomplete actual-system matrix",
    )
    expected = {
        f"{name}-{method}-scale-{scale}-translated-{translated}"
        for name in (
            "point-charges",
            "counterions",
            "pairwise",
            "fragment-complete-C8H18",
            "opls-peo",
            "opls-resp-aromatic",
            "opls-network",
        )
        for method in ("pme", "ewald")
        for scale, translated in ((1.0, False), (1.5, False), (1.5, True))
    }
    require(
        {r["name"] for r in rows} == expected, "Missing system/box/translation control"
    )
    for row in rows:
        path = root / row["name"]
        check_native(path, row)
        backend = read(root / "additional-inputs" / (row["case"] + "-backend.json"))
        f = read(path / "frame.json")
        c = completion(backend, f, row["lammps"]["setting"], row["lammps"])
        require(
            c == row["completion"]
            and row["comparison"] == comparison(row["openmm"], row["lammps"], c),
            "Changed family-aware completion",
        )
    base = read(root / "c11-rerun/box-1.5/frame.json")
    pressure = []
    for h in (1e-4, 1e-5, 1e-6):
        pairs = []
        for sign in (-1, 1):
            path = root / f"strained-pressure-{h}-{sign}"
            row = read(path / "comparison.json")
            check_native(path, row)
            f = read(path / "frame.json")
            c = completion(b, f, row["lammps"]["setting"], row["lammps"])
            require(
                c == row["completion"]
                and corrected(row["lammps"], c) == row["completed"],
                "Changed pressure completion",
            )
            pairs.append(row)
        derivatives = {
            engine: -(
                pairs[1][engine]["energy_kj_mol"] - pairs[0][engine]["energy_kj_mol"]
            )
            / (2 * h * np.prod(base["box_lengths_nm"]))
            * KJ_MOL_NM3_TO_BAR
            for engine in ("openmm", "lammps", "completed")
        }
        pressure.append(
            {
                "relative_volume_step": h,
                "configurational_volume_derivative_bar": derivatives,
            }
        )
    require(
        pressure == read(root / "strained-pressure.json"), "Changed pressure derivative"
    )
    return result


def validate(path):
    path = Path(path)
    r = read(path)
    require(
        r == signed(r)
        and r["schema"] == "island.j40-c12-receipt.v1"
        and r["base_commit"] == BASE,
        "Changed receipt",
    )
    require(
        r["historical_files"] == historical_hashes(), "Changed historical inventory"
    )
    require(set(r["implementation_sources"]) == set(SOURCES), "Missing tested source")
    gates = {
        "focused",
        "regressions",
        "full",
        "ruff",
        "compileall",
        "diff_check",
        "pip_supported",
        "pip_foyer",
        "pip_analysis",
        "pip_amber",
        "pip_quantum",
    }
    require(
        set(r["validation"]) == gates
        and all(v["exit_code"] == 0 for v in r["validation"].values()),
        "Missing/failed verification gate",
    )
    for name, digest in r["implementation_sources"].items():
        require(sha(name) == digest, "Changed implementation")
        import hashlib

        tested = subprocess.check_output(
            ["git", "show", r["tested_implementation_commit"] + ":" + name]
        )
        require(hashlib.sha256(tested).hexdigest() == digest, "Untested implementation")
    for name, digest in r["files"].items():
        require(sha(path.parent / name) == digest, "Changed evidence: " + name)
    require(
        set(r["files"])
        == {
            "artifacts.tar.xz",
            "artifact-inventory.json",
            "manifest.json",
            "summary.json",
            "validation.tar.xz",
        },
        "Missing evidence assets",
    )
    manifest = read(path.parent / "manifest.json")
    require(
        manifest == signed(manifest) and manifest["identity"] == r["manifest_identity"],
        "Changed manifest",
    )
    require(
        manifest["schema"] == "island.j40-c12-manifest.v1"
        and manifest["tested_implementation_commit"]
        == r["tested_implementation_commit"]
        and manifest["artifacts_archive_sha256"] == r["files"]["artifacts.tar.xz"],
        "Changed manifest binding",
    )
    with tempfile.TemporaryDirectory(prefix="j40-c12-offline-") as temporary:
        root = Path(temporary)
        with tarfile.open(path.parent / "artifacts.tar.xz") as archive:
            archive.extractall(root, filter="data")
        require(
            {str(p.relative_to(root)): sha(p) for p in root.rglob("*") if p.is_file()}
            == r["artifact_inventory"],
            "Changed artifact inventory",
        )
        result = validate_study(root / "study")
        require(result["identity"] == r["study_identity"], "Changed study identity")
    require(
        r["scientific_acceptance"] is False
        and r["production_validated"] is False
        and r["simulation_readiness"] == "not_established",
        "False scientific claim",
    )
    return r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("receipt", type=Path)
    p.add_argument("--diagnostic-only", action="store_true")
    a = p.parse_args()
    validate(a.receipt)
    print(json.dumps({"integrity": "passed", "scientific_acceptance": False}))
    return 0 if a.diagnostic_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
