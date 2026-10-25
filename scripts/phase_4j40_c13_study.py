"""Extend a freshly executed corrected C12 matrix with bounded PME sweeps."""

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.phase_4j40_c11_numeric import signed, write
from scripts.phase_4j40_c12 import comparison, completion, corrected
from scripts.phase_4j40_c12_engine import lammps_point, openmm_point
from scripts.phase_4j40_c12_study import frame_for
from scripts.phase_4j40_c13 import POLICY, characterize, validate_summary


def read(p):
    return json.loads(p.read_text())


def execute(root, executable):
    rows = read(root / "additional-matrix.json")
    # Fixed independent sweep, declared before results: no parameter selection.
    b = read(root / "additional-inputs/gaff2-pe-backend.json")
    template = (root / "additional-inputs/gaff2-pe.data").read_text()
    base = read(root / "gaff2-pe-pme-scale-1.0-translated-False/frame.json")
    sweeps = []
    plan = [("baseline", 32, 5, 1e-7, 8.0, 0.5)]
    plan += [("mesh", m, 5, 1e-7, 8.0, 0.5) for m in (24, 48)]
    plan += [("order", 32, o, 1e-7, 8.0, 0.5) for o in (4, 6)]
    plan += [("tolerance", 32, 5, t, 8.0, 0.5) for t in (1e-5, 1e-9)]
    plan += [("alpha", 32, 5, 1e-7, a, 0.5) for a in (6.0, 10.0)]
    plan += [("cutoff", 32, 5, 1e-7, 8.0, c) for c in (0.4, 0.48)]
    write(root / "c13-plan.json", plan)
    for index, (axis, mesh, order, tolerance, alpha, cutoff) in enumerate(plan):
        for translation in (False, True):
            xyz = np.array(base["positions_nm"])
            box = np.array(base["box_lengths_nm"])
            if translation:
                xyz = np.mod(xyz + [0.17, 0.23, 0.11], box)
            f = frame_for(b, xyz, box)
            name = f"c13-sweep-{index}-{axis}-translated-{translation}"
            p = root / name
            p.mkdir()
            write(p / "frame.json", f)
            setting = {
                "method": "pme",
                "cutoff_nm": cutoff,
                "tolerance": tolerance,
                "alpha_nm": alpha,
                "mesh": [mesh] * 3,
                "order": order,
            }
            om = openmm_point(b, xyz, box, setting, f["identity"], p / "openmm")
            ls = dict(setting, method="pppm")
            lm = lammps_point(
                template, xyz, box, ls, f["identity"], p / "lammps", executable
            )
            c = completion(b, f, ls, lm)
            row = {
                "name": name,
                "frame_identity": f["identity"],
                "openmm": om,
                "lammps": lm,
                "completion": c,
                "comparison": comparison(om, lm, c),
            }
            write(p / "comparison.json", row)
            rows.append(row)
            sweeps.append(row)
    write(root / "c13-mesh-sweep.json", sweeps)
    results = [characterize(r, read(root / r["name"] / "frame.json")) for r in rows]
    original = read(root / "c11-rerun/study.json")
    # Include all original cutoffs, reciprocal tolerances, replica frames and strains.
    for row in original["rows"]:
        f = read(root / "c11-rerun" / row["name"] / "frame.json")
        c = completion(b, f, row["lammps"]["setting"], row["lammps"])
        results.append(
            characterize(dict(row, frame_identity=f["identity"], completion=c), f)
        )
    # Fresh two-sided, same-frame volume derivatives for every actual polymer.
    pressures = []
    from scripts.phase_4j40_c11_numeric import KJ_MOL_NM3_TO_BAR

    for case in ("gaff2-pe", "opls-peo", "opls-resp-aromatic", "opls-network"):
        backend = read(root / "additional-inputs" / (case + "-backend.json"))
        template = (root / "additional-inputs" / (case + ".data")).read_text()
        base = read(root / f"{case}-ewald-scale-1.0-translated-False/frame.json")
        origin = next(r for r in rows if r.get("case") == case) if case else None
        for h in (1e-4, 1e-5, 1e-6):
            pair = []
            for sign in (-1, 1):
                scale = (1 + sign * h) ** (1 / 3)
                xyz = np.array(base["positions_nm"]) * scale
                box = np.array(base["box_lengths_nm"]) * scale
                f = frame_for(backend, xyz, box)
                p = root / f"c13-pressure-{case}-{h}-{sign}"
                p.mkdir()
                write(p / "frame.json", f)
                setting = origin["lammps"]["setting"]
                om = openmm_point(
                    backend, xyz, box, setting, f["identity"], p / "openmm"
                )
                lm = lammps_point(
                    template, xyz, box, setting, f["identity"], p / "lammps", executable
                )
                c = completion(backend, f, setting, lm)
                write(p / "completion.json", c)
                pair.append({"openmm": om, "completed": corrected(lm, c)})
            derivatives = {
                engine: -(
                    pair[1][engine]["energy_kj_mol"] - pair[0][engine]["energy_kj_mol"]
                )
                / (2 * h * np.prod(base["box_lengths_nm"]))
                * KJ_MOL_NM3_TO_BAR
                for engine in ("openmm", "completed")
            }
            pressures.append(
                {
                    "case": case,
                    "relative_volume_step": h,
                    "finite_volume_pressure_bar": derivatives,
                    "cross_engine_difference_bar": derivatives["openmm"]
                    - derivatives["completed"],
                    "status": "pressure_unavailable",
                    "reason": "Finite-volume agreement alone does not validate independent OpenMM force/virial tensor",
                }
            )
    write(root / "c13-pressure.json", pressures)
    result = signed(
        {
            "schema": "island.j40-c13-study.v1",
            "policy": POLICY,
            "rows": results,
            "pressure": pressures,
            "convergence_limitations": [
                "Fixed-parameter one-variable sweeps do not establish convergence over independent equilibrated replicas",
                "OpenMM interpolation order fixed at 5; LAMMPS order varies",
                "Single/mixed precision and alternate FFT builds unavailable",
                "No independently measured full OpenMM reciprocal virial tensor",
                "Repeated translations are invariance controls, not independent thermodynamic replicas",
            ],
            "ewald_status": "ewald_unvalidated",
            "pressure_status": "pressure_unavailable",
            "cutoff_recovery": "not_executed",
            "long_npt": "not_executed",
            "scientific_acceptance": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    validate_summary(result)
    write(root / "c13-summary.json", result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--lammps", type=Path, required=True)
    a = p.parse_args()
    execute(a.run, a.lammps.resolve())
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
