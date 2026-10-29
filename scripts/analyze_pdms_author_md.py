"""Offline analysis of the fixed author's local PDMS diagnostic MD reproduction.

This never constructs an ISLAND preparation or grants an execution profile.
"""

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

PINS = {
    "opls_system.data": "8b03e6e8f3c51525f011d9d64e02d69d81faab78cd17ba66676b208d1c4b9c24",
    "opls_system.in": "dfd22dac270ad11bf50001573e7d581ea924e287036e34aa85cded67913d6cdd",
    "opls_system.in.init": "4ae3754ef4f428f60eff9f1a33f80a040ec17b45e41ac1e48bba1be3d3ee2592",
    "opls_system.in.settings": "d22a19757398bdf9fb04536f4acf74a7ad3dabce686aaba65d58675bba0016e5",
}


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_atoms(path):
    rows, active = {}, False
    for raw in Path(path).read_text().splitlines():
        tokens = raw.split("#", 1)[0].split()
        if tokens == ["Atoms"]:
            active = True
        elif active and tokens:
            if not tokens[0].isdigit():
                break
            rows[int(tokens[0])] = tokens
    return rows


def final_dump(path):
    lines = Path(path).read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    columns = lines[index].split()[2:]
    rows = [
        dict(zip(columns, line.split(), strict=True)) for line in lines[index + 1 :]
    ]
    return {int(row["id"]): row for row in rows}


def thermodynamics(path):
    result, columns = [], None
    for line in Path(path).read_text().splitlines():
        tokens = line.split()
        if tokens[:3] == ["Step", "Time", "Atoms"]:
            columns = tokens
        elif columns:
            try:
                row = dict(zip(columns, map(float, tokens), strict=True))
            except ValueError:
                columns = None
                continue
            if not all(math.isfinite(v) for v in row.values()):
                raise ValueError("nonfinite_thermodynamic_sample")
            result.append(row)
    return result


def analyze(directory):
    root = Path(directory)
    for name, pin in PINS.items():
        if file_sha256(root / name) != pin:
            raise ValueError("source_mutated:" + name)
    source = source_atoms(root / "opls_system.data")
    final = final_dump(root / "final-state.lammpstrj")
    resumed = final_dump(root / "relocated-resume/resumed-final-state.lammpstrj")
    if len(source) != 32448 or set(source) != set(final) or set(final) != set(resumed):
        raise ValueError("atom_count_or_ids_changed")
    for i, row in source.items():
        for state in (final[i], resumed[i]):
            if not all(math.isfinite(float(v)) for v in state.values()):
                raise ValueError("nonfinite_final_atom_state")
            if (
                int(state["mol"]) != int(row[1])
                or int(state["type"]) != int(row[2])
                or float(state["q"]).hex() != float(row[3]).hex()
            ):
                raise ValueError("molecule_type_or_charge_changed")
    frames = []
    with (root / "trajectory.lammpstrj").open() as stream:
        while line := stream.readline():
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError("trajectory_header_invalid")
            step = int(next(stream))
            if next(stream).strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("trajectory_count_header_invalid")
            count = int(next(stream))
            if count != len(source) or not next(stream).startswith("ITEM: BOX BOUNDS"):
                raise ValueError("trajectory_atom_loss_or_box_header_invalid")
            for _ in range(3):
                next(stream)
            columns = next(stream).split()[2:]
            ids = set()
            for _ in range(count):
                row = dict(zip(columns, next(stream).split(), strict=True))
                i = int(row["id"])
                if i not in source or i in ids:
                    raise ValueError("trajectory_atom_identity_changed")
                ids.add(i)
                if (
                    int(row["mol"]) != int(source[i][1])
                    or int(row["type"]) != int(source[i][2])
                    or float(row["q"]).hex() != float(source[i][3]).hex()
                    or not all(math.isfinite(float(v)) for v in row.values())
                ):
                    raise ValueError(
                        "trajectory_identity_charge_or_finite_state_changed"
                    )
            frames.append(step)
    if frames != list(range(0, 1001, 100)):
        raise ValueError("trajectory_incomplete")
    samples = thermodynamics(root / "md-nvt.log")
    if (
        len(samples) != 21
        or samples[0]["Step"] != 0
        or samples[-1]["Step"] != 1000
        or samples[-1]["Time"] != 500
        or any(s["Atoms"] != 32448 for s in samples)
        or "ERROR" in (root / "md-nvt.log").read_text()
    ):
        raise ValueError("MD_not_completed_or_atom_loss")
    differences = {
        key: max(abs(float(final[i][key]) - float(resumed[i][key])) for i in final)
        for key in ("x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz")
    }
    dump_lines = (root / "final-state.lammpstrj").read_text().splitlines()
    bounds = next(
        i for i, line in enumerate(dump_lines) if line.startswith("ITEM: BOX BOUNDS")
    )
    lengths = [
        float(row.split()[1]) - float(row.split()[0])
        for row in dump_lines[bounds + 1 : bounds + 4]
    ]
    unwrapped_differences, minimum_image_differences = {}, {}
    for axis, length in zip(("x", "y", "z"), lengths, strict=True):
        unwrapped_differences[axis] = max(
            abs(
                float(final[i][axis])
                + int(final[i]["i" + axis]) * length
                - float(resumed[i][axis])
                - int(resumed[i]["i" + axis]) * length
            )
            for i in final
        )
        deltas = [float(final[i][axis]) - float(resumed[i][axis]) for i in final]
        minimum_image_differences[axis] = max(
            abs(delta - length * round(delta / length)) for delta in deltas
        )
    fd = json.loads((root / "verified-finite-difference-report.json").read_text())
    report = {
        "schema": "island_author_PDMS_local_diagnostic_MD_v1",
        "diagnostic_MD_completed": True,
        "scientific_status": "scientifically_unqualified",
        "public_ISLAND_PDMS_execution_qualified": False,
        "execution_profiles_created": [],
        "acceptance_gate": 1,
        "model": "Author DFT/NPA-charge OPLS-compatible PDMS input; original switched-LJ/PPPM settings",
        "paper_doi": "10.1021/acs.jpcb.4c08471",
        "author_commit": "0053fda8449797e7e2b0790426cfc52dcd023e3e",
        "source_sha256": PINS,
        "atoms": len(source),
        "chains": len({row[1] for row in source.values()}),
        "chain_sizes": dict(
            Counter(Counter(row[1] for row in source.values()).values())
        ),
        "ensemble": "NVT",
        "temperature_target_K": 300,
        "thermostat": "Nose-Hoover; damping 100 fs; source velocities retained",
        "timestep_fs": 0.5,
        "steps": 1000,
        "elapsed_simulation_ps": 0.5,
        "trajectory_frames_verified": len(frames),
        "pppm_accuracy": 0.0001,
        "pppm_order": 5,
        "coulomb_table_bits": 12,
        "temperature_K": {
            "initial": samples[0]["Temp"],
            "final": samples[-1]["Temp"],
            "minimum": min(s["Temp"] for s in samples),
            "maximum": max(s["Temp"] for s in samples),
            "mean_of_21_samples": sum(s["Temp"] for s in samples) / len(samples),
        },
        "potential_energy_kcal_mol": {
            "initial": samples[0]["PotEng"],
            "final": samples[-1]["PotEng"],
        },
        "fixed_box_density_g_cm3": samples[0]["Density"],
        "density_is_an_imposed_NVT_input_not_an_equilibrium_prediction": True,
        "identity_checks": {
            "atom_ids_unchanged": True,
            "molecule_membership_unchanged": True,
            "type_ids_unchanged": True,
            "charges_bitwise_equal_to_source": True,
            "normalization_or_redistribution": False,
        },
        "sampled_strict_FD_control": {
            "settings": fd["settings"],
            "comparisons": fd["comparisons"],
            "maximum_error_kj_mol_angstrom": fd["maximum_error_kj_mol_angstrom"],
            "threshold_kj_mol_angstrom": fd["threshold_kj_mol_angstrom"],
            "passed": fd["passed"],
            "scope": "Separate tightened-solver control; MD used original author solver settings.",
        },
        "fresh_process_relocated_continuation": {
            "restart_step": 500,
            "final_step": 1000,
            "max_absolute_differences": differences,
            "max_unwrapped_coordinate_differences_angstrom": unwrapped_differences,
            "max_minimum_image_coordinate_differences_angstrom": minimum_image_differences,
            "raw_coordinates_may_differ_by_periodic_images": True,
            "all_dump_fields_equal": final == resumed,
            "checkpoint_bytes_equal": (
                (root / "final.restart").read_bytes()
                == (root / "relocated-resume/resumed-final.restart").read_bytes()
            ),
            "RNG": "No stochastic RNG used: source velocities + deterministic Nose-Hoover.",
        },
        "production_validated": False,
        "simulation_readiness": "not_established",
        "NPT_run": False,
        "minimization_run": False,
        "initial_state": "Unmodified author equilibrated data and velocities; no new packing or graph typing.",
        "scientific_limitations": [
            "Strict sampled finite-difference control did not pass.",
            "No independent actual-PDMS cross-engine numerical qualification.",
            "No equilibrium or experimental agreement established by 0.5 ps.",
            "Native forcefields and existing Ewald/pressure/cutoff/long-NPT gates unchanged.",
        ],
    }
    return report, samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()
    report, samples = analyze(args.run)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    with args.csv.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
        writer.writeheader()
        writer.writerows(samples)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
