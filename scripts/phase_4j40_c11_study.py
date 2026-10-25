"""Execute a bounded diagnostic sweep, never a production parameter change."""

import argparse
import copy
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from scripts.phase_4j40_c11_engine import finite_strain, lammps_point, openmm_point
from scripts.phase_4j40_c11_numeric import compare, require, sha, signed, write


def settings():
    # This finite matrix is fixed before inspecting the sweep results.
    rows = [
        {"name": "baseline", "method": "ewald", "tolerance": 1e-5, "cutoff_nm": 0.5}
    ]
    rows += [
        {
            "name": f"tolerance-{tol}",
            "method": "ewald",
            "tolerance": tol,
            "cutoff_nm": 0.5,
        }
        for tol in (1e-4, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
    ]
    rows += [
        {"name": f"cutoff-{rc}", "method": "ewald", "tolerance": 1e-9, "cutoff_nm": rc}
        for rc in (0.4, 0.45)
    ]
    rows += [
        {"name": f"pme-pppm-{tol}", "method": "pme", "tolerance": tol, "cutoff_nm": 0.5}
        for tol in (1e-4, 1e-5, 1e-6, 1e-7)
    ]
    rows += [
        {
            "name": f"mesh-{n}",
            "method": "pme",
            "tolerance": 1e-7,
            "cutoff_nm": 0.5,
            "alpha_nm": 7.3,
            "mesh": [n] * 3,
            "order": 5,
        }
        for n in (16, 24, 32, 48, 64)
    ]
    rows += [
        {
            "name": f"order-{n}",
            "method": "pme",
            "tolerance": 1e-7,
            "cutoff_nm": 0.5,
            "alpha_nm": 7.3,
            "mesh": [48] * 3,
            "order": n,
        }
        for n in (3, 7)
    ]
    rows += [
        {
            "name": f"alpha-{a}",
            "method": "ewald",
            "tolerance": 1e-9,
            "cutoff_nm": 0.5,
            "lammps_alpha_nm": a,
        }
        for a in (6.0, 7.3, 9.0)
    ]
    rows += [
        {
            "name": f"kmax-{n}",
            "method": "ewald",
            "tolerance": 1e-9,
            "cutoff_nm": 0.5,
            "lammps_alpha_nm": 7.3,
            "kmax": n,
        }
        for n in (4, 6, 8, 10)
    ]
    rows += [
        {
            "name": f"box-{s}",
            "method": "ewald",
            "tolerance": 1e-9,
            "cutoff_nm": 0.5,
            "box_scale": s,
        }
        for s in (1.25, 1.5)
    ]
    rows += [
        {
            "name": "translation",
            "method": "ewald",
            "tolerance": 1e-9,
            "cutoff_nm": 0.5,
            "translation_nm": [0.17, 0.23, 0.11],
        }
    ]
    return rows


def execute(backend_file, reference, executable, output, frames_root):
    import openmm as mm

    output.mkdir(parents=True, exist_ok=False)
    write(
        output / "invocation.json",
        {
            "argv": sys.argv,
            "implementation_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            "source_hashes": {
                str(p): sha(p) for p in Path("scripts").glob("*phase_4j40_c11*.py")
            },
        },
    )
    inputs = output / "inputs"
    inputs.mkdir()
    shutil.copyfile(backend_file, inputs / "backend.json")
    for n in ("system.data", "in.periodic", "forces.dump", "execution.log"):
        shutil.copyfile(reference / "cross-engine" / n, inputs / n)
    shutil.copyfile(
        reference / "cross-engine.json", inputs / "historical-comparison.json"
    )
    backend = json.loads(backend_file.read_text())
    historical = json.loads((inputs / "historical-comparison.json").read_text())
    require(
        backend["identity"] == historical["backend_identity"],
        "Wrong historical backend",
    )
    template = (inputs / "system.data").read_text()
    sites = sorted(backend["numerical_data"]["sites"], key=lambda s: s["id"])
    xyz = (
        np.asarray([backend["system"]["coordinates"][str(s["id"])] for s in sites]) / 10
    )
    lengths = np.asarray(backend["config"]["box_lengths"]) / 10
    identity_data = {
        k: backend[k]
        for k in (
            "final_graph_identity",
            "assignment_identity",
            "charge_identity",
            "typed_graph_identity",
            "numerical_data",
        )
    }

    def frame_identity(x, box):
        return signed(
            dict(
                identity_data,
                atom_ids=[s["id"] for s in sites],
                positions_nm=x.tolist(),
                box_lengths_nm=box.tolist(),
            )
        )

    plan = settings()
    write(
        output / "plan.json",
        signed(
            {
                "settings": plan,
                "scientific_acceptance": False,
                "policy": "fixed bounded diagnostic sweep; no optimization or production setting selection",
            }
        ),
    )
    rows = []

    def point_pair(name, x, box, setting):
        root = output / name
        root.mkdir()
        binding = frame_identity(x, box)
        write(root / "frame.json", binding)
        oms = copy.deepcopy(setting)
        lms = copy.deepcopy(setting)
        lms["method"] = "pppm" if setting["method"] == "pme" else "ewald"
        if "lammps_alpha_nm" in lms:
            lms["alpha_nm"] = lms.pop("lammps_alpha_nm")
        # OpenMM's interpolation order is fixed by its implementation; an order
        # sweep changes only LAMMPS and records that difference explicitly.
        a = openmm_point(backend, x, box, oms, binding["identity"], root / "openmm")
        b = lammps_point(
            template, x, box, lms, binding["identity"], root / "lammps", executable
        )
        require(
            b["charges_e"] == [s["charge"] for s in sites], "Changed LAMMPS charges"
        )
        row = {
            "name": name,
            "setting": setting,
            "openmm": a,
            "lammps": b,
            "comparison": compare(a, b),
        }
        write(root / "comparison.json", row)
        rows.append(row)
        return row

    for setting in plan:
        scale = setting.get("box_scale", 1.0)
        x, box = xyz * scale, lengths * scale
        if "translation_nm" in setting:
            x = np.mod(x + np.asarray(setting["translation_nm"]), box)
        point_pair(setting["name"], x, box, setting)
    # Same graph/charges/assignment on independent retained replica frames.
    retained = []
    for replica in (0, 1):
        source = frames_root / f"replica-{replica}" / "uninterrupted/frames.json"
        shutil.copyfile(source, inputs / f"replica-{replica}-frames.json")
        frames = json.loads(source.read_text())
        if isinstance(frames, dict):
            frames = frames["frames"]
        for index in (0, len(frames) - 1):
            f = frames[index]
            box = np.asarray(f["box_vectors_nm"])
            require(
                np.array_equal(box, np.diag(np.diag(box))), "Unexpected triclinic frame"
            )
            x, box = np.asarray(f["positions_nm"]), np.diag(box)
            for tol in (1e-5, 1e-9):
                name = f"replica-{replica}-frame-{index}-tol-{tol}"
                row = point_pair(
                    name,
                    x,
                    box,
                    {"method": "ewald", "tolerance": tol, "cutoff_nm": 0.5},
                )
                retained.append(
                    {
                        "name": name,
                        "source": str(source),
                        "step": f["step"],
                        "source_frame_identity": signed(f)["identity"],
                        "comparison": row["comparison"],
                    }
                )
    # Isotropic derivatives use the actual polymer model, no point-charge proxy.
    pressure = {}
    for method in ("ewald", "pme"):
        cfg = {
            "method": method,
            "tolerance": 1e-9 if method == "ewald" else 1e-7,
            "cutoff_nm": 0.5,
        }
        pressure[method] = {}
        for engine in ("openmm", "lammps"):
            if engine == "lammps":
                cfg = dict(cfg, method="ewald" if method == "ewald" else "pppm")
            pressure[method][engine] = finite_strain(
                backend,
                xyz,
                lengths,
                cfg,
                frame_identity(xyz, lengths)["identity"],
                output / (method + "-strain-" + engine),
                engine,
                template,
                executable,
            )
    help_run = subprocess.run(
        [str(executable), "-h"], text=True, capture_output=True, check=False
    )
    (output / "lammps-help.txt").write_text(help_run.stdout + help_run.stderr)
    library = Path(mm.__file__).parent / "_openmm.cpython-311-x86_64-linux-gnu.so"
    inventory = {
        "python": sys.version,
        "python_executable": sys.executable,
        "openmm": mm.version.version,
        "openmm_module": str(Path(mm.__file__).resolve()),
        "openmm_library_hash_if_found": sha(library) if library.exists() else None,
        "platforms": [
            mm.Platform.getPlatform(i).getName()
            for i in range(mm.Platform.getNumPlatforms())
        ],
        "system": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "lammps_executable": str(executable.resolve()),
        "lammps_sha256": sha(executable),
        "environment": {
            k: os.environ.get(k)
            for k in ("OMP_NUM_THREADS", "OPENMM_CPU_THREADS", "PYTHONPATH")
        },
        "precision": "OpenMM Reference double; LAMMPS DOUBLE DOUBLE FFT from raw help; GPU single/mixed not available",
        "unsupported": {
            "openmm_interpolation_order": "fixed at 5 by pinned implementation",
            "fft_precision_sweep": "pinned LAMMPS FFT_DOUBLE build; no alternate precision binary declared",
            "single_gpu": "no CUDA/OpenCL platform available",
        },
    }
    write(output / "environment.json", inventory)
    result = signed(
        {
            "schema": "island.j40-c11-study.v2",
            "baseline": rows[0],
            "rows": rows,
            "retained_frames": retained,
            "pressure_strain": pressure,
            "source_files": {p.name: sha(p) for p in inputs.iterdir()},
            "scientific_acceptance": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "long_npt_executed": False,
        }
    )
    write(output / "study.json", result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--backend", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--lammps", type=Path, required=True)
    p.add_argument("--frames-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    execute(a.backend, a.reference, a.lammps.resolve(), a.output, a.frames_root)
    return 1  # Diagnostic convergence cannot promote scientific acceptance.


if __name__ == "__main__":
    raise SystemExit(main())
