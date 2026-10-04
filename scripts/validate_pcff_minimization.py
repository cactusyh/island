"""Bounded PCFF minimization of checksummed H5 inputs; no regeneration/retries."""

import argparse
import json
import time
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path
from unittest.mock import patch

import numpy as np
from validate_pcff_singlepoint import (
    GROUPS,
    aa_reference_coefficients,
    data_sections,
    digest,
    lammps_input,
    run,
)

from island.charge_references.records import unpack
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import (
    load_pcff_model,
    load_pcff_parameters,
    load_pcff_source,
)
from island.forcefields.pcff.automatic import chemical_graph, graph_system
from island.minimization import MinimizationOptions, minimize_geometry
from island.workflows import storage
from island.workflows.bundle import record, system_data

CASES = ("butane", "ethanol", "PE_DP3", "PEO_DP3")
OPTIONS = MinimizationOptions(
    force_tolerance=0.1, max_iterations=5000, max_evaluations=10000
)


def inputs(root, name, source):
    directory = root / name
    assignment = load_pcff_parameters(directory / "parameters.json", source)
    system = graph_system(
        unpack(assignment.json_text)["charge_record"]["automatic_typing"]["graph"]
    )
    for sid, xyz in json.loads(
        (directory / "p0-f0/coordinates.json").read_text()
    ).items():
        system.coordinates.set(int(sid), xyz)
    model = load_pcff_model(directory / "model-0.json", assignment, system=system)
    policy = unpack(model.json_text)["special_pairs"]
    if policy["lj"] != [0, 0, 1] or policy["coulomb"] != [0, 0, 1]:
        raise ValueError("Wrong H5 pair policy")
    return system, model


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--evidence", default="docs/evidence/phase_4h5.json")
    parser.add_argument("--lammps", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    root = Path(args.inputs).resolve()
    out = Path(args.output).resolve()
    evidence = json.loads(Path(args.evidence).read_text())
    # Check ALL published H5 files before consuming any scientific record.
    for name, expected in evidence["files"].items():
        if digest(root / name) != expected:
            raise ValueError("H5 artifact hash mismatch: " + name)
    executable = Path(args.lammps).resolve()
    if digest(executable) != evidence["upstream"]["lammps_sha256"]:
        raise ValueError("Require the verified H5 LAMMPS executable")
    source = load_pcff_source(args.source)
    out.mkdir(parents=True, exist_ok=False)
    plan = {
        "schema": "island_pcff_minimization_acceptance_v1",
        "cases": CASES,
        "inputs": evidence["files"],
        "h5_evidence_sha256": digest(args.evidence),
        "source": source.identity,
        "lammps_sha256": digest(executable),
        "options": asdict(OPTIONS),
        "lj": [0, 0, 1],
        "coulomb": [0, 0, 1],
        "energy_atol_kj_mol": 1e-5,
        "force_atol_kj_mol_angstrom": 1e-5,
        "rtol": 2e-10,
        "reference": "H5 independent converter inventories, with independently reconstructed ABC/ABD/CBD angle-angle roles",
        "versions": {k: version(k) for k in ("numpy", "openmm", "scipy")},
        "retries": 0,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    storage.publish(out / "declaration.json", storage.json_bytes(plan))
    import openmm as mm

    outcomes = []
    h5rows = {r["case"]: r for r in evidence["report"]["results"]}
    for name in CASES:
        folder = out / name
        folder.mkdir()
        row = {"case": name, "passed": False}
        try:
            system, model = inputs(root, name, source)
            before = storage.json_bytes(storage.encode(system_data(system)))
            storage.publish(folder / "starting-system.json", before)
            start = time.perf_counter()
            evaluator = PCFFSinglePointEvaluator(system, model)
            row["binding_compilation_seconds"] = time.perf_counter() - start
            constructions = []
            constructor = mm.Context

            def context(*a, constructions=constructions, constructor=constructor, **kw):
                constructions.append(time.perf_counter())
                return constructor(*a, **kw)

            with patch.object(mm, "Context", context):
                start = time.perf_counter()
                session = evaluator.open_session()
                row["session_construction_seconds"] = time.perf_counter() - start
                with session:
                    start = time.perf_counter()
                    minimum = minimize_geometry(system, session, OPTIONS)
                    row["optimization_including_final_verification_seconds"] = (
                        time.perf_counter() - start
                    )
            minimum.validate_integrity()
            storage.publish(
                folder / "minimum.json",
                storage.json_bytes(storage.encode(record(minimum))),
            )
            row.update(
                {
                    k: getattr(minimum, k)
                    for k in (
                        "initial_energy",
                        "final_energy",
                        "initial_fmax",
                        "final_fmax",
                        "initial_rms_force",
                        "final_rms_force",
                        "iterations",
                        "evaluations",
                        "termination_reason",
                        "converged",
                        "final_evaluation_verified",
                    )
                }
            )
            row["context_constructions"] = len(constructions)
            row["extra_diagnostic_evaluations"] = 0
            row["parameter_identity"] = evaluator.parameter_fingerprint
            row["model_identity"] = evaluator.model_fingerprint
            row["initial_coordinate_fingerprint"] = (
                minimum.initial_evaluation.coordinate_fingerprint
            )
            row["final_coordinate_fingerprint"] = (
                minimum.final_evaluation.coordinate_fingerprint
            )
            expected = h5rows[name]["comparisons"][0]
            if (
                row["initial_coordinate_fingerprint"]
                != expected["coordinate_fingerprint"]
                or row["parameter_identity"] != expected["parameter_fingerprint"]
            ):
                raise ValueError("H5 starting geometry/model identity changed")
            ids = sorted(system.topology.sites)
            xyz = np.array([minimum.coordinates[s] for s in ids])
            raw = root / name / "reference.data"
            (folder / "reference.data").write_bytes(raw.read_bytes())
            overrides = aa_reference_coefficients(data_sections(raw))
            storage.publish(
                folder / "aa_reference_commands.json", storage.json_bytes(overrides)
            )
            (folder / "input.lmp").write_text(
                lammps_input(xyz, [0, 0, 1], [0, 0, 1], overrides)
            )
            row["lammps_seconds"] = run(
                [executable, "-in", "input.lmp"], folder, "lammps"
            )
            line = [
                s
                for s in (folder / "lammps.log").read_text().splitlines()
                if s.startswith("ENERGIES ")
            ][-1]
            energies = np.array(list(map(float, line.split()[1:]))) * 4.184
            dump = (
                (folder / "forces.dump")
                .read_text()
                .split("ITEM: ATOMS id x y z fx fy fz\n")[1]
            )
            forces = (
                np.array(
                    [
                        [float(v) for v in s.split()[4:7]]
                        for s in dump.strip().splitlines()
                    ]
                )
                * 4.184
            )
            final = minimum.final_evaluation
            actual = np.array([final.forces[s] for s in ids])
            calculated = np.array(
                [final.potential_energy]
                + [
                    sum(final.energy_components[k] for k in group)
                    for group in GROUPS.values()
                ]
            )
            row.update(
                energy_max_abs_error_kj_mol=float(np.max(abs(calculated - energies))),
                force_max_abs_error_kj_mol_angstrom=float(np.max(abs(actual - forces))),
                lammps_fmax=float(np.linalg.norm(forces, axis=1).max()),
                lammps_rms_force=float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
                island_energies=calculated.tolist(),
                lammps_energies=energies.tolist(),
            )
            if not np.allclose(
                calculated, energies, atol=1e-5, rtol=2e-10
            ) or not np.allclose(actual, forces, atol=1e-5, rtol=2e-10):
                raise ValueError("Independent LAMMPS gate failed")
            if before != storage.json_bytes(storage.encode(system_data(system))):
                raise ValueError("Input mutation")
            if not (
                minimum.converged
                and minimum.final_evaluation_verified
                and minimum.final_fmax <= 0.1
                and row["lammps_fmax"] <= 0.1
                and minimum.final_energy
                <= minimum.initial_energy + OPTIONS.energy_increase_tolerance
                and len(constructions) == 2
            ):
                raise ValueError(
                    "Minimization/resource gate unmet: " + minimum.termination_reason
                )
            copied = minimum.to_system(system)
            copied.validate()
            model.validate_integrity(copied)
            if chemical_graph(copied) != chemical_graph(system):
                raise ValueError("Applied chemical graph changed")
            storage.publish(
                folder / "optimized-system.json",
                storage.json_bytes(storage.encode(system_data(copied))),
            )
            row["passed"] = True
        except Exception as exc:  # noqa: BLE001 -- retain failed bounded experiments
            row["failure"] = f"{type(exc).__name__}: {exc}"
        row["files"] = {p.name: digest(p) for p in folder.iterdir() if p.is_file()}
        storage.publish(folder / "outcome.json", storage.json_bytes(row))
        outcomes.append(row)
        print(
            name,
            row["passed"],
            row.get("final_fmax"),
            row.get("failure", ""),
            flush=True,
        )
    for name, expected in evidence["files"].items():
        if digest(root / name) != expected:
            raise ValueError("H5 artifact changed during execution")
    passed = all(r["passed"] for r in outcomes)
    storage.publish(
        out / "report.json",
        storage.json_bytes(
            {
                "cases": outcomes,
                "passed": passed,
                "declaration_sha256": digest(out / "declaration.json"),
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
    )
    return int(not passed)


if __name__ == "__main__":
    raise SystemExit(main())
