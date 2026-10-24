"""Offline raw-output and historical-byte validator; acceptance stays nonzero."""

import argparse
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

from scripts.phase_4j40_c11_numeric import compare, parse_dump, require, sha, signed
from scripts.validate_phase_4j40_c11 import characterize

BASE = "d8a4e878d660f5ac1b1965ac077c24ad052b23e6"
SOURCES = (
    "scripts/phase_4j40_c11_numeric.py",
    "scripts/phase_4j40_c11_engine.py",
    "scripts/phase_4j40_c11_study.py",
    "scripts/validate_phase_4j40_c11.py",
    "scripts/validate_phase_4j40_c11_receipt.py",
    "tests/test_phase_4j40_c11.py",
)


def historical_hashes():
    rows = {}
    for item in subprocess.check_output(
        ["git", "ls-tree", "-rz", BASE, "docs/evidence"]
    ).split(b"\0"):
        if not item:
            continue
        meta, name = item.split(b"\t")
        mode, _, blob = meta.decode().split()
        path = Path(name.decode())
        actual = (
            str(path.readlink()).encode() if mode == "120000" else path.read_bytes()
        )
        expected = subprocess.check_output(["git", "cat-file", "blob", blob])
        require(actual == expected, "Changed historical bytes: " + str(path))
        rows[str(path)] = {"git_blob": blob, "mode": mode, "resolved_sha256": sha(path)}
    return rows


def validate_study(root):
    study = json.loads((root / "study.json").read_text())
    require(study == signed(study), "Changed study identity")
    from scripts.phase_4j40_c11_study import settings

    plan = json.loads((root / "plan.json").read_text())
    require(plan["settings"] == settings(), "Changed bounded sweep plan")
    require(
        len(study["rows"]) == 38 and len(study["retained_frames"]) == 8,
        "Missing actual evaluation",
    )
    require(
        study["scientific_acceptance"] is False and not study["long_npt_executed"],
        "False readiness or dynamics claim",
    )
    for name, digest in study["source_files"].items():
        require(
            sha(root / "inputs" / name) == digest, "Changed source artifact: " + name
        )
    for row in study["rows"]:
        path = root / row["name"]
        frame = json.loads((path / "frame.json").read_text())
        require(frame == signed(frame), "Changed physical frame identity")
        require(
            row["openmm"]["frame_identity"]
            == row["lammps"]["frame_identity"]
            == frame["identity"],
            "Different comparison frames",
        )
        for engine in ("openmm", "lammps"):
            native = json.loads((path / engine / "point.json").read_text())
            require(native == row[engine], "Changed engine point")
        require(
            row["comparison"] == compare(row["openmm"], row["lammps"]),
            "Changed errors/bounds",
        )
        dump = parse_dump(path / "lammps/forces.dump")
        require(
            dump["forces_kj_mol_nm"] == row["lammps"]["forces_kj_mol_nm"],
            "Changed native force vector",
        )
        require(
            sha(path / "openmm/system.xml") == row["openmm"]["system_xml_sha256"],
            "Changed actual system",
        )
        for name, digest in row["lammps"]["raw_files"].items():
            require(sha(path / "lammps" / name) == digest, "Changed KSPACE raw output")
    return characterize(root)


def validate(path):
    path = Path(path)
    r = json.loads(path.read_text())
    require(r == signed(r), "Changed receipt identity")
    require(
        r["schema"] == "island.j40-c11-receipt.v2" and r["base_commit"] == BASE,
        "Wrong C11 schema/base",
    )
    require(
        r["scientific_acceptance"] is False
        and r["ewald_status"] == "ewald_unvalidated"
        and r["pressure_status"] == "pressure_unavailable",
        "Blocked result relabeled",
    )
    require(
        r["production_validated"] is False
        and r["simulation_readiness"] == "not_established",
        "Changed readiness",
    )
    require(r["historical_files"] == historical_hashes(), "Changed historical evidence")
    for name, digest in r["files"].items():
        require(sha(path.parent / name) == digest, "Changed evidence asset: " + name)
    for name in SOURCES:
        require(
            sha(name) == r["implementation_sources"][name], "Changed implementation"
        )
        tested = subprocess.check_output(
            ["git", "show", r["tested_implementation_commit"] + ":" + name]
        )
        import hashlib

        require(
            hashlib.sha256(tested).hexdigest() == sha(name), "Untested implementation"
        )
    with tempfile.TemporaryDirectory(prefix="j40-c11-offline-") as tmp:
        root = Path(tmp)
        with tarfile.open(path.parent / "artifacts.tar.xz") as archive:
            archive.extractall(root, filter="data")
        actual = {
            str(p.relative_to(root)): sha(p) for p in root.rglob("*") if p.is_file()
        }
        require(actual == r["artifact_inventory"], "Changed archive file inventory")
        c = validate_study(root / "study")
        require(
            c == json.loads((path.parent / "characterization.json").read_text()),
            "Changed statistics or characterization",
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
