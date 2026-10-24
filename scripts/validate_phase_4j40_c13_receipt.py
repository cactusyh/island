"""Offline byte-integrity and raw numerical reconstruction for J40-C13."""

import argparse
import json
import math
import tarfile
from pathlib import Path

from scripts.phase_4j40_c11_numeric import require, sha, signed
from scripts.phase_4j40_c13 import characterize, validate_summary


def validate(directory):
    receipt = json.loads((directory / "receipt.json").read_text())
    require(signed(receipt) == receipt, "Mutated receipt")
    manifest = json.loads((directory / "manifest.json").read_text())
    require(signed(manifest) == manifest, "Mutated manifest")
    require(
        receipt["manifest_sha256"] == sha(directory / "manifest.json"),
        "Manifest changed",
    )
    for name, digest in manifest["files"].items():
        require(sha(directory / name) == digest, "Changed evidence: " + name)
    summary = json.loads((directory / "summary.json").read_text())
    validate_summary(summary)
    require(receipt["summary_identity"] == summary["identity"], "Changed summary")
    with tarfile.open(directory / "artifacts.tar.xz") as archive:

        def read(name):
            return json.load(archive.extractfile("study/" + name))

        for row in summary["rows"]:
            prefix = row["name"]
            if prefix.startswith("c13-") or (prefix + "/comparison.json") in {
                m.name.removeprefix("study/") for m in archive.getmembers()
            }:
                raw = read(prefix + "/comparison.json")
                frame = read(prefix + "/frame.json")
            else:
                from scripts.phase_4j40_c12 import completion

                raw = next(
                    r
                    for r in read("c11-rerun/study.json")["rows"]
                    if r["name"] == prefix
                )
                frame = read("c11-rerun/" + prefix + "/frame.json")
                backend = read("additional-inputs/gaff2-pe-backend.json")
                raw = dict(
                    raw,
                    frame_identity=frame["identity"],
                    completion=completion(
                        backend, frame, raw["lammps"]["setting"], raw["lammps"]
                    ),
                )
            require(
                characterize(raw, frame) == row,
                "Numerical reconstruction mismatch: " + prefix,
            )
        from scripts.phase_4j40_c11_numeric import KJ_MOL_NM3_TO_BAR
        from scripts.phase_4j40_c12 import corrected

        for row in summary["pressure"]:
            case, h = row["case"], row["relative_volume_step"]
            frame = read(f"{case}-ewald-scale-1.0-translated-False/frame.json")
            energies = {k: [] for k in ("openmm", "completed")}
            for sign in (-1, 1):
                prefix = f"c13-pressure-{case}-{h}-{sign}/"
                a = read(prefix + "openmm/point.json")
                b = read(prefix + "lammps/point.json")
                c = read(prefix + "completion.json")
                energies["openmm"].append(a["energy_kj_mol"])
                energies["completed"].append(corrected(b, c)["energy_kj_mol"])
            for engine, pair in energies.items():
                pressure = (
                    -(pair[1] - pair[0])
                    / (2 * h * math.prod(frame["box_lengths_nm"]))
                    * KJ_MOL_NM3_TO_BAR
                )
                require(
                    math.isclose(
                        pressure,
                        row["finite_volume_pressure_bar"][engine],
                        rel_tol=1e-14,
                        abs_tol=1e-8,
                    ),
                    "Pressure reconstruction mismatch",
                )
            require(
                row["status"] == "pressure_unavailable",
                "Unsupported pressure qualification",
            )
    repo = Path(__file__).resolve().parents[1]
    for name, digest in manifest["historical_files"].items():
        require(sha(repo / name) == digest, "Historical bytes changed: " + name)
    for name, digest in manifest["implementation_files"].items():
        require(sha(repo / name) == digest, "Implementation changed: " + name)
    require(receipt["scientific_acceptance"] is False, "Unsupported acceptance")
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Unsupported readiness",
    )
    return {
        "integrity_exit_code": 0,
        "scientific_acceptance_exit_code": 1,
        "rows_reconstructed": len(summary["rows"]),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    args = p.parse_args()
    print(json.dumps(validate(args.directory), indent=2))


if __name__ == "__main__":
    main()
