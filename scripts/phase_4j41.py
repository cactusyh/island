"""Execute bounded PCFF Si/P source, converter and finite native acceptance."""

import argparse
import json
import os
import shutil
import subprocess
import sys
from hashlib import sha256
from pathlib import Path

import numpy as np
from validate_oplsaa_dynamics import continuation_options, load_segment
from validate_pcff_dynamics import reconstruct, run_case
from validate_pcff_singlepoint import (
    GROUPS,
    aa_reference_coefficients,
    check_inventory,
    data_sections,
    lammps_input,
    reference_inputs,
    run,
)

from island.chemistry import from_smiles
from island.dynamics import (
    create_dynamics_checkpoint,
    load_dynamics_checkpoint,
    resume_dynamics,
    save_dynamics_checkpoint,
)
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import (
    assess_pcff_sip,
    load_pcff_source,
    prepare_pcff_sip,
    save_pcff_model,
    save_pcff_parameters,
)
from island.minimization import MinimizationOptions, minimize_geometry
from island.workflows import storage
from island.workflows.bundle import system_data
from scripts.pcff_j41_reference import audit_queries, inventory, labels_for, require

CASES = [
    ("silane", "[SiH4]", "silane_source"),
    ("methylsilane", "C[SiH3]", "organosilicon_source"),
    ("siloxane", "C[Si](C)(C)O[Si](C)(C)C", "siloxane_source"),
    ("phosphate_ester", "COP(=O)(OC)OC", "phosphate_ester_audit"),
    ("phosphonate", "CP(=O)(OC)OC", "phosphonate_audit"),
]


def write(path, value):
    path.write_bytes(storage.json_bytes(storage.encode(value)))


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def components(result):
    return np.array(
        [result.potential_energy]
        + [sum(result.energy_components[k] for k in keys) for keys in GROUPS.values()]
    )


def finite_differences(evaluator, system):
    """Central energy derivatives for every coordinate and each compiled component."""
    import openmm as mm
    from openmm import unit

    ids = sorted(system.topology.sites)
    xyz = np.array([system.coordinates.get(i) for i in ids])
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(
        mm.XmlSerializer.deserialize(evaluator._xml),
        integrator,
        mm.Platform.getPlatformByName("Reference"),
    )
    records = []
    for group, family in enumerate(evaluator._components):
        context.setPositions(xyz * 0.1)
        force = (
            context.getState(getForces=True, groups=1 << group)
            .getForces(asNumpy=True)
            .value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            * 0.1
        )
        errors = []
        for h in (0.001, 0.0001, 1e-05):
            numerical = np.zeros_like(xyz)
            for i in range(len(ids)):
                for axis in range(3):
                    energies = []
                    for sign in (-1, 1):
                        trial = xyz.copy()
                        trial[i, axis] += sign * h
                        context.setPositions(trial * 0.1)
                        energies.append(
                            context.getState(getEnergy=True, groups=1 << group)
                            .getPotentialEnergy()
                            .value_in_unit(unit.kilojoule_per_mole)
                        )
                    numerical[i, axis] = -(energies[1] - energies[0]) / (2 * h)
            errors.append(float(np.max(abs(numerical - force))))
        records.append(
            {
                "family": family,
                "steps_angstrom": [0.001, 0.0001, 1e-05],
                "max_force_errors": errors,
                "passed": errors[-1] <= 1e-05,
            }
        )
    require(all(r["passed"] for r in records), records)
    return records


def compare_frame(system, evaluator, folder):
    ids = sorted(system.topology.sites)
    xyz = np.array([system.coordinates.get(i) for i in ids])
    result = evaluator.evaluate(dict(zip(ids, xyz)))
    line = next(
        x
        for x in (folder / "lammps.log").read_text().splitlines()
        if x.startswith("ENERGIES ")
    )
    energies = np.array(list(map(float, line.split()[1:]))) * 4.184
    text = (
        (folder / "forces.dump").read_text().split("ITEM: ATOMS id x y z fx fy fz\n")[1]
    )
    forces = (
        np.array(
            [list(map(float, line.split()[4:7])) for line in text.strip().splitlines()]
        )
        * 4.184
    )
    native = np.array([result.forces[i] for i in ids])
    record = {
        "native_energies_kj_mol": components(result).tolist(),
        "reference_energies_kj_mol": energies.tolist(),
        "native_forces_kj_mol_angstrom": native.tolist(),
        "reference_forces_kj_mol_angstrom": forces.tolist(),
        "energy_max_abs": float(max(abs(components(result) - energies))),
        "force_max_abs": float(np.max(abs(native - forces))),
        "energy_force_atol": 1e-05,
        "rtol": 2e-10,
    }
    require(np.allclose(components(result), energies, atol=1e-05, rtol=2e-10), record)
    require(np.allclose(native, forces, atol=1e-05, rtol=2e-10), record)
    return record


def coefficient_audit(model, sections):
    """Silane has one type per interaction; compare all executable coefficients."""
    from math import pi

    result = []
    specification = model.payload
    for family, section in (
        ("quartic_bond", "Bond Coeffs"),
        ("quartic_angle", "Angle Coeffs"),
        ("bond-bond", "BondBond Coeffs"),
        ("bond-angle", "BondAngle Coeffs"),
        ("angle-angle", "AngleAngle Coeffs"),
    ):
        require(
            len(sections[section]) == 1,
            "Expected independently declared single silane interaction type",
        )
        converted = list(map(float, sections[section][0][1:]))
        if family == "quartic_bond":
            native = [
                t["coefficients"]
                for t in specification["terms"]
                if t["family"] == family
            ]
            expected = [converted[0], *[v * 4.184 for v in converted[1:]]]
        elif family == "quartic_angle":
            native = [
                t["coefficients"]
                for t in specification["terms"]
                if t["family"] == family
            ]
            expected = [converted[0] * pi / 180, *[v * 4.184 for v in converted[1:]]]
        elif family in ("bond-bond", "bond-angle"):
            native = [
                t["coefficients"] + t["equilibria"][:2]
                for t in specification["terms"]
                if t["family"] == family
            ]
            count = 1 if family == "bond-bond" else 2
            expected = [*[v * 4.184 for v in converted[:count]], *converted[count:]]
        else:
            require(
                converted[:3] == [converted[0]] * 3
                and converted[3:] == [converted[3]] * 3,
                "Silane AA symmetry mismatch",
            )
            native = [
                t["coefficients"] + t["equilibria"]
                for t in specification["terms"]
                if t["family"] == family
            ]
            expected = [
                converted[0] * 4.184,
                converted[3] * pi / 180,
                converted[3] * pi / 180,
            ]
        require(
            bool(native)
            and all(
                np.allclose(row, expected, atol=1e-10, rtol=1e-13) for row in native
            ),
            "Converter coefficient mismatch: " + family,
        )
        result.append(
            {
                "family": family,
                "converter_coefficients_real_units": converted,
                "native_rows_checked": len(native),
                "normalized_expected": expected,
                "passed": True,
            }
        )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--distribution", type=Path, required=True)
    p.add_argument("--converter", type=Path, required=True)
    p.add_argument("--lammps", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    source_path, converter, lammps = (
        (a.distribution / "pcff.frc").resolve(),
        a.converter.resolve(),
        a.lammps.resolve(),
    )
    source = load_pcff_source(source_path)
    for filename in ("pcff.frc", "pcff.rlb", "pcff_templates.dat"):
        shutil.copyfile(a.distribution / filename, a.output / filename)
    write(a.output / "source-inventory.json", inventory(source.raw))
    report = {
        "schema": "island.j41-execution.v1",
        "argv": sys.argv,
        "pid": os.getpid(),
        "source_sha256": digest(source_path),
        "distribution_hashes": {
            n: digest(a.distribution / n)
            for n in ("pcff.frc", "pcff.rlb", "pcff_templates.dat")
        },
        "converter_sha256": digest(converter),
        "lammps_sha256": digest(lammps),
        "converter_authority": "supplied types and charges only; no chemical perception or charge assignment; no -ignore",
        "cases": [],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    for name, smiles, suffix in CASES:
        folder = a.output / name
        folder.mkdir()
        system = from_smiles(smiles, random_seed=2026)
        ids = sorted(system.topology.sites)
        xyz = np.array([system.coordinates.get(i) for i in ids])
        xyz += 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3)
        for i, position in zip(ids, xyz):
            system.coordinates.set(i, position)
        write(folder / "system.json", system_data(system))
        assessment = assess_pcff_sip(
            system, source, profile="island_pcff_" + suffix + "_v1"
        )
        (folder / "assessment.json").write_text(assessment.json_text)
        row = {
            "name": name,
            "smiles": smiles,
            "seed": 2026,
            "assessment_identity": assessment.identity,
            "graph_identity": assessment.payload["graph_identity"],
            "source_complete": assessment.payload["complete"],
            "diagnostics": assessment.payload["diagnostics"],
            "charge_origin": assessment.payload["charge_origin"],
            "dynamics_executed": False,
        }
        if "P" not in {s.element for s in system.topology.sites.values()}:
            audit = audit_queries(system, assessment, source.raw)
            write(folder / "independent-audit.json", audit)
            labels = labels_for(system)
            charges = {
                int(i): float(q) for i, q in audit["charges"]["partial_charges"].items()
            }
            reference_inputs(system, labels, charges, xyz, folder)
            command = [
                str(converter),
                "reference",
                "-class",
                "II",
                "-frc",
                str(source_path),
                "-p",
                "3",
                "-nocenter",
            ]
            result = subprocess.run(
                command,
                cwd=folder,
                capture_output=True,
                text=True,
                check=False,
                timeout=180,
            )
            (folder / "converter.log").write_text(result.stdout + result.stderr)
            row["converter"] = {"argv": command, "exit_code": result.returncode}
            if result.returncode == 0:
                sections = data_sections(folder / "reference.data")
                sections.setdefault("Dihedrals", [])
                neighbors = {i: [] for i in ids}
                for b in system.topology.bonds.values():
                    neighbors[b.site1].append(b.site2)
                    neighbors[b.site2].append(b.site1)
                row["converter"]["inventory"] = check_inventory(
                    sections, ids, neighbors, charges
                )
                row["converter"]["zero_rows_are_not_missing_source_authority"] = True
        else:
            row["converter"] = {
                "status": "not_executed",
                "reason": "No independently justified P/oxygen atom typing; converter consumes rather than validates supplied chemical types",
            }
        if assessment.payload["complete"]:
            require(
                row["converter"]["exit_code"] == 0, "Independent J41 control failed"
            )
            model = prepare_pcff_sip(system, assessment)
            save_pcff_parameters(model.assignment, folder / "parameters.json")
            save_pcff_model(model, folder / "model.json")
            row["coefficient_audit"] = coefficient_audit(model, sections)
            row["parameter_identity"] = model.assignment.identity
            row["model_identity"] = model.identity
            evaluator = PCFFSinglePointEvaluator(system, model)
            overrides = aa_reference_coefficients(sections)
            write(folder / "aa-overrides.json", overrides)
            (folder / "input.lmp").write_text(
                lammps_input(xyz, [0, 0, 1], [0, 0, 1], overrides)
            )
            run([lammps, "-in", "input.lmp"], folder, "lammps")
            row["singlepoint"] = compare_frame(system, evaluator, folder)
            row["finite_differences"] = finite_differences(evaluator, system)
            minimum = minimize_geometry(
                system,
                evaluator,
                MinimizationOptions(
                    force_tolerance=0.1, max_iterations=5000, max_evaluations=10000
                ),
            )
            require(
                minimum.converged
                and minimum.final_evaluation_verified
                and (minimum.final_fmax <= 0.1),
                "Independent J41 control failed",
            )
            optimized = minimum.to_system(system)
            evaluator = PCFFSinglePointEvaluator(optimized, model)
            row["minimization"] = {
                "converged": minimum.converged,
                "final_fmax": minimum.final_fmax,
            }
            inputs = (
                system,
                optimized,
                model.assignment,
                model,
                minimum,
                source,
                evaluator,
            )
            row["dynamics"] = []
            for kind in ("nve", "baoab"):
                dynamics = run_case(a.output, name, kind, inputs, a.output, lammps)
                moved = a.output / f"{name}-{kind}"
                s, e, _ = reconstruct(moved)
                checkpoint = load_dynamics_checkpoint(moved / "checkpoint.json")
                continued = resume_dynamics(checkpoint, s, e, continuation_options())
                save_dynamics_checkpoint(
                    create_dynamics_checkpoint(continued),
                    moved / "in-process-checkpoint.json",
                )
                require(
                    (moved / "in-process-checkpoint.json").read_bytes()
                    == (moved / "continued-checkpoint.json").read_bytes(),
                    "Independent J41 control failed",
                )
                whole = load_segment(moved / "whole.json")
                require(
                    whole.final_state.coordinates == continued.final_state.coordinates,
                    "Independent J41 control failed",
                )
                require(
                    whole.final_state.velocities == continued.final_state.velocities,
                    "Independent J41 control failed",
                )
                require(
                    whole.payload["rng"] == continued.payload["rng"],
                    "Independent J41 control failed",
                )
                dynamics["fresh_process_checkpoint_byte_equal"] = True
                dynamics["uninterrupted_state_and_rng_exact"] = True
                dynamics["checkpoint_comparison_scope"] = (
                    "identical segmented lineage in-process versus relocated fresh process; whole-run checkpoint retains distinct lineage/counters"
                )
                row["dynamics"].append(dynamics)
            row["dynamics_executed"] = True
        else:
            row["workflow_status"] = "blocked_before_model_publication"
        report["cases"].append(row)
        write(a.output / "report.json", report)
        print(name, "source_complete=", row["source_complete"], flush=True)
    report["acceptance_gate"] = int(
        not all(r["source_complete"] for r in report["cases"])
    )
    write(a.output / "report.json", report)
    return report["acceptance_gate"]


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--reconstruct":
        from scripts.validate_phase_4j41_receipt import validate_study

        root = Path(sys.argv[2])
        result = validate_study(root)
        for kind in ("nve", "baoab"):
            reconstruct(root / ("silane-" + kind))
        print(
            json.dumps(
                {**result, "fresh_pid": os.getpid(), "native_reconstruction": True}
            )
        )
    else:
        raise SystemExit(main())
