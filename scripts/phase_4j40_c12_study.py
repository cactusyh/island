"""Execute C11 again and a family-aware C12 KSPACE/PME control matrix."""

import argparse
import copy
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from scripts.phase_4j40_c11_numeric import sha, signed, write
from scripts.phase_4j40_c11_study import execute as c11_execute
from scripts.phase_4j40_c12 import comparison, completion, corrected
from scripts.phase_4j40_c12_engine import lammps_point, openmm_point


def read(path):
    return json.loads(Path(path).read_text())


def frame_for(b, xyz, box):
    return signed(
        {
            **{
                k: b[k]
                for k in (
                    "final_graph_identity",
                    "assignment_identity",
                    "charge_identity",
                    "typed_graph_identity",
                    "numerical_data",
                )
            },
            "atom_ids": [s["id"] for s in b["numerical_data"]["sites"]],
            "positions_nm": np.asarray(xyz).tolist(),
            "box_lengths_nm": np.asarray(box).tolist(),
        }
    )


def template_for(b):
    from island._periodic_lammps import render

    coords = {int(k): np.asarray(v) for k, v in b["system"]["coordinates"].items()}
    return render(
        SimpleNamespace(payload=b, system=SimpleNamespace(coordinates=coords))
    )[0]


def controls(parent):
    for name, charges, positions in (
        ("point-charges", [1.0, -1.0], [[0.2, 0.3, 0.4], [0.65, 0.4, 0.5]]),
        (
            "counterions",
            [1.0, -0.5, -0.5],
            [[0.2, 0.3, 0.4], [0.7, 0.4, 0.5], [0.8, 0.9, 0.7]],
        ),
        ("pairwise", [0.2, -0.2], [[0.4, 0.4, 0.4], [0.72, 0.4, 0.4]]),
    ):
        b = copy.deepcopy(parent)
        b["numerical_data"] = {
            "family": "GAFF2",
            "sites": [
                {
                    "id": i + 1,
                    "atom_type": "analytic",
                    "mass": 12.0,
                    "sigma": 2.0,
                    "epsilon": 0.0,
                    "charge": q,
                }
                for i, q in enumerate(charges)
            ],
            **{k: [] for k in ("bonds", "angles", "torsions", "impropers", "pairs")},
        }
        for k in (
            "final_graph_identity",
            "assignment_identity",
            "charge_identity",
            "typed_graph_identity",
        ):
            b[k] = signed({"analytic_control": name, "charges": charges, "role": k})[
                "identity"
            ]
        b["system"]["coordinates"] = {
            str(i + 1): (np.array(x) * 10).tolist() for i, x in enumerate(positions)
        }
        b["config"]["box_lengths"] = [20.0, 20.0, 20.0]
        b["final_graph"] = {
            "molecule_membership": {str(i + 1): i + 1 for i in range(len(charges))}
        }
        yield (
            name,
            b,
            template_for(b),
            "analytic_control; no polymer charge-origin claim",
        )
    b = copy.deepcopy(parent)
    ids = [s["id"] for s in b["numerical_data"]["sites"] if s["id"] <= 26]
    b["numerical_data"]["sites"] = [
        s for s in b["numerical_data"]["sites"] if s["id"] in ids
    ]
    for key in ("bonds", "angles", "torsions", "impropers", "pairs"):
        b["numerical_data"][key] = [
            r for r in b["numerical_data"][key] if all(i in ids for i in r["sites"])
        ]
    b["system"]["coordinates"] = {
        k: v for k, v in b["system"]["coordinates"].items() if int(k) in ids
    }
    b["final_graph"] = {"molecule_membership": {str(i): 1 for i in ids}}
    for k in (
        "final_graph_identity",
        "assignment_identity",
        "charge_identity",
        "typed_graph_identity",
    ):
        b[k] = signed(
            {"parent_identity": b[k], "exact_complete_component_atom_ids": ids}
        )["identity"]
    yield (
        "fragment-complete-C8H18",
        b,
        template_for(b),
        "retained identical RESP/Psi4 component; no new charge generation",
    )


def execute(args):
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    write(
        out / "invocation.json",
        {
            "argv": sys.argv,
            "implementation_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            "lammps": str(args.lammps.resolve()),
            "lammps_sha256": sha(args.lammps),
            "source_hashes": {
                str(p): sha(p) for p in Path("scripts").glob("*phase_4j40_c12*.py")
            },
        },
    )
    original = c11_execute(
        args.source / "gaff2-pe/backend.json",
        args.reference / "gaff2-pe",
        args.lammps.resolve(),
        out / "c11-rerun",
        args.source / "gaff2-pe",
    )
    backend = read(out / "c11-rerun/inputs/backend.json")
    completed = []
    for row in original["rows"]:
        f = read(out / "c11-rerun" / row["name"] / "frame.json")
        c = completion(backend, f, row["lammps"]["setting"], row["lammps"])
        completed.append(
            {
                "name": row["name"],
                "completion": c,
                "comparison": comparison(row["openmm"], row["lammps"], c),
            }
        )
    write(out / "c11-completed.json", completed)
    inputs = out / "additional-inputs"
    inputs.mkdir()
    cases = list(controls(backend))
    for key in ("opls-peo", "opls-resp-aromatic", "opls-network"):
        b = read(args.source / key / "backend.json")
        template = (args.reference / key / "cross-engine/system.data").read_text()
        cases.append((key, b, template, b["config"]["provenance"]))
    rows = []
    for name, b, template, origin in cases:
        write(inputs / (name + "-backend.json"), b)
        (inputs / (name + ".data")).write_text(template)
        xyz = (
            np.array(
                [
                    b["system"]["coordinates"][str(s["id"])]
                    for s in b["numerical_data"]["sites"]
                ]
            )
            / 10
        )
        box = np.array(b["config"]["box_lengths"]) / 10
        for method in ("ewald", "pme"):
            for scale, translation in ((1.0, False), (1.5, False), (1.5, True)):
                label = f"{name}-{method}-scale-{scale}-translated-{translation}"
                x, lengths = xyz * scale, box * scale
                if translation:
                    x = np.mod(x + np.array([0.17, 0.23, 0.11]), lengths)
                f = frame_for(b, x, lengths)
                path = out / label
                path.mkdir()
                write(path / "frame.json", f)
                setting = {
                    "method": method,
                    "cutoff_nm": b["config"]["nonbonded_cutoff"] / 10,
                    "tolerance": 1e-7 if method == "pme" else 1e-9,
                    "mixing": "arithmetic"
                    if b["numerical_data"]["family"] == "GAFF2"
                    else "geometric",
                    "dihedral_style": "fourier"
                    if b["numerical_data"]["family"] == "GAFF2"
                    else "multi/harmonic",
                    "lj14": b["config"]["one_four_scaling"][0],
                    "coul14": b["config"]["one_four_scaling"][1],
                }
                om = openmm_point(
                    b, x, lengths, setting, f["identity"], path / "openmm"
                )
                setting = dict(setting, method="pppm" if method == "pme" else "ewald")
                lm = lammps_point(
                    template,
                    x,
                    lengths,
                    setting,
                    f["identity"],
                    path / "lammps",
                    args.lammps.resolve(),
                )
                c = completion(b, f, setting, lm)
                row = {
                    "name": label,
                    "case": name,
                    "charge_origin": origin,
                    "frame_identity": f["identity"],
                    "openmm": om,
                    "lammps": lm,
                    "completion": c,
                    "comparison": comparison(om, lm, c),
                }
                write(path / "comparison.json", row)
                rows.append(row)
    # Explicit completion pressure is checked against a new energy-volume derivative
    # at the first failing geometry, with native calculations on both sides.
    base = read(out / "c11-rerun/box-1.5/frame.json")
    pressure = []
    for h in (1e-4, 1e-5, 1e-6):
        pairs = []
        for sign in (-1, 1):
            factor = (1 + sign * h) ** (1 / 3)
            x = np.array(base["positions_nm"]) * factor
            lengths = np.array(base["box_lengths_nm"]) * factor
            f = frame_for(backend, x, lengths)
            setting = {"method": "ewald", "tolerance": 1e-9, "cutoff_nm": 0.5}
            path = out / f"strained-pressure-{h}-{sign}"
            path.mkdir()
            write(path / "frame.json", f)
            om = openmm_point(
                backend, x, lengths, setting, f["identity"], path / "openmm"
            )
            lm = lammps_point(
                (out / "c11-rerun/inputs/system.data").read_text(),
                x,
                lengths,
                setting,
                f["identity"],
                path / "lammps",
                args.lammps.resolve(),
            )
            c = completion(backend, f, setting, lm)
            row = {
                "openmm": om,
                "lammps": lm,
                "completion": c,
                "completed": corrected(lm, c),
            }
            write(path / "comparison.json", row)
            pairs.append(row)
        volume = np.prod(base["box_lengths_nm"])
        from scripts.phase_4j40_c11_numeric import KJ_MOL_NM3_TO_BAR

        derivatives = {
            engine: -(
                pairs[1][engine]["energy_kj_mol"] - pairs[0][engine]["energy_kj_mol"]
            )
            / (2 * h * volume)
            * KJ_MOL_NM3_TO_BAR
            for engine in ("openmm", "lammps", "completed")
        }
        pressure.append(
            {
                "relative_volume_step": h,
                "configurational_volume_derivative_bar": derivatives,
            }
        )
    write(out / "additional-matrix.json", rows)
    write(out / "strained-pressure.json", pressure)
    result = signed(
        {
            "schema": "island.j40-c12-study.v1",
            "c11_rows": len(completed),
            "additional_rows": len(rows),
            "strained_pressure_pairs": 6,
            "scientific_acceptance": False,
            "ewald_status": "ewald_unvalidated",
            "pressure_status": "pressure_unavailable",
            "cutoff_recovery": "not_executed",
            "long_npt": "not_executed",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "comparison_semantics": "native outputs plus explicitly labelled single-frame analytical completion; production exporter remains conservative",
            "c11_original_identity": original["identity"],
        }
    )
    write(out / "study.json", result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--lammps", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    execute(args)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
