"""Independent raw-output residual and pressure audit (no comparison imports)."""

import argparse
import json
import math
from pathlib import Path


def audit(root):
    summary = json.loads((root / "c13-summary.json").read_text())
    reports = []
    for row in summary["rows"]:
        p = root / row["name"]
        if not p.exists():
            p = root / "c11-rerun" / row["name"]
        a = json.loads((p / "openmm/point.json").read_text())
        b = json.loads((p / "lammps/point.json").read_text())
        cp = p / "comparison.json"
        if cp.exists() and "completion" in json.loads(cp.read_text()):
            c = json.loads(cp.read_text())["completion"]
        else:
            c = next(
                x["completion"]
                for x in json.loads((root / "c11-completed.json").read_text())
                if x["name"] == row["name"]
            )
        force = []
        for x, y, z in zip(
            a["forces_kj_mol_nm"],
            b["forces_kj_mol_nm"],
            c["forces_kj_mol_nm"],
            strict=True,
        ):
            force.extend(u - (v + w) for u, v, w in zip(x, y, z, strict=True))
        rms = math.sqrt(math.fsum(d * d for d in force) / len(force))
        reported = row["metrics"]["forces"]["rms_component_error"]
        assert math.isclose(rms, reported, rel_tol=2e-15, abs_tol=1e-15)
        energy = a["energy_kj_mol"] - (
            b["energy_kj_mol"] + c["energy_lj_kj_mol"] + c["energy_coulomb_kj_mol"]
        )
        # Alternative addition order may differ by float64 roundoff.
        bound = 8 * math.ulp(max(abs(a["energy_kj_mol"]), abs(b["energy_kj_mol"]), 1))
        assert abs(energy - row["metrics"]["signed_energy_error_kj_mol"]) <= bound
        reports.append(
            {
                "name": row["name"],
                "rms_force_error": rms,
                "energy_error_kj_mol": energy,
                "reconstruction_roundoff_bound": bound,
            }
        )
    pressure = []
    for row in summary["pressure"]:
        case, h = row["case"], row["relative_volume_step"]
        base = json.loads(
            (root / f"{case}-ewald-scale-1.0-translated-False/frame.json").read_text()
        )
        volume = math.prod(base["box_lengths_nm"])
        values = {k: [] for k in ("openmm", "completed")}
        for sign in (-1, 1):
            p = root / f"c13-pressure-{case}-{h}-{sign}"
            a = json.loads((p / "openmm/point.json").read_text())
            b = json.loads((p / "lammps/point.json").read_text())
            c = json.loads((p / "completion.json").read_text())
            values["openmm"].append(a["energy_kj_mol"])
            values["completed"].append(
                b["energy_kj_mol"]
                + (c["energy_lj_kj_mol"] + c["energy_coulomb_kj_mol"])
            )
        # Exact SI constants yield the same conversion independently.
        conversion = 1000 / 6.02214076e23 / 1e-27 / 1e5
        derivatives = {
            k: -(v[1] - v[0]) / (2 * h * volume) * conversion for k, v in values.items()
        }
        for k, derivative in derivatives.items():
            assert math.isclose(
                derivative,
                row["finite_volume_pressure_bar"][k],
                rel_tol=1e-12,
                abs_tol=1e-5,
            )
        pressure.append(
            {"case": case, "h": h, "independent_finite_volume_bar": derivatives}
        )
    convergence = []
    for translated in (False, True):
        for engine in ("openmm", "lammps"):
            points = []
            for index in (1, 0, 2):  # 24, 32, 48 cubed, all other settings fixed
                axis = "baseline" if index == 0 else "mesh"
                p = (
                    root
                    / f"c13-sweep-{index}-{axis}-translated-{translated}"
                    / engine
                    / "point.json"
                )
                points.append(json.loads(p.read_text()))
            ref = points[-1]
            for mesh, point in zip((24, 32, 48), points, strict=True):
                differences = [
                    x - y
                    for a, b in zip(
                        point["forces_kj_mol_nm"], ref["forces_kj_mol_nm"], strict=True
                    )
                    for x, y in zip(a, b, strict=True)
                ]
                convergence.append(
                    {
                        "engine": engine,
                        "translated": translated,
                        "mesh": mesh,
                        "reference_mesh": 48,
                        "energy_change_to_finest_kj_mol": point["energy_kj_mol"]
                        - ref["energy_kj_mol"],
                        "rms_force_change_to_finest": math.sqrt(
                            math.fsum(x * x for x in differences) / len(differences)
                        ),
                        "interpretation": "Finest-mesh self-difference is zero by construction; not a reference-accuracy or acceptance assertion",
                    }
                )
    return {
        "mesh_convergence_by_engine": convergence,
        "schema": "island.j40-c13-independent-audit.v1",
        "all_passed": True,
        "residuals": reports,
        "pressure": pressure,
        "scope": "Independent arithmetic reconstruction of raw numerical residuals and finite-volume derivatives; does not qualify full reciprocal virial pressure",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    args.output.write_text(
        json.dumps(audit(args.run), indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
